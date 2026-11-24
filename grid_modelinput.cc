// The read of the input ejecta model:
// - model.txt gives the densities and the mass fractions of the nuclides;
// - abundances.txt gives the mass fractions of the elements;
// - the optional Ye.txt gives the electron fractions.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <format>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "grid_internal.h"
#include "input.h"
#include "inputfilestream.h"
#include "mpi_logging.h"

namespace grid {

namespace {

// The format version of model.txt. A file of version 2 starts with the line "format_version 2" and has a header of
// keys with the cell edges. A file of version 1 has no version line
int model_format_version{1};
// The columns X_<element>_otherstable of a model.txt of format version 2. Such a column holds the mass fraction of
// the stable isotopes of an element that have no nuclide column. The first list gives the atomic number of each
// column. The node leaders hold the values of each model cell and column until read_elem_abundances() uses them.
std::vector<int> model_otherstable_atomic_numbers;
std::vector<float> model_otherstable_massfracs;
constexpr std::string_view otherstable_suffix = "_otherstable";

// The count of the cell lines that the optional key cell_line_count of a model.txt of format version 2 gives. A
// cell with no line is empty, thus only this count finds a file that lost lines at its end.
std::optional<ptrdiff_t> model_cell_line_count;

double min_den{-1.};  // minimum model density

int first_input_cellid{-1};  // auto-determine first cell index in model.txt (usually 1 or 0)

void set_initelectronfrac(const int modelgridindex, const float electronfrac) {
  if (std::isnan(electronfrac) || electronfrac < 0. || electronfrac > 1.001) {
    fatal_crash("input Ye {:g} for cell {} is outside the physical range [0, 1]", electronfrac, modelgridindex);
  }
  modelgrid_input[modelgridindex].initelectronfrac = electronfrac;
}

void read_possible_yefile() {
  if (model_format_version == 2 && inputfile_exists("Ye.txt")) {
    fatal_crash(
        "Ye.txt is present, but model.txt has format version 2. The column Ye of such a model.txt gives the electron "
        "fraction. Remove Ye.txt");
  }
  if (!inputfile_exists("Ye.txt")) {
    printlnlog("Ye.txt is not present, so the model keeps the electron fractions of model.txt");
    return;
  }

  // the electron fractions are written to node-shared memory, so only the node leaders read the file
  // (synchronised by the barrier below)
  if (globals::rank_in_node == 0) {
    auto filein = istream_required("Ye.txt");
    int nlines_in = 0;
    if (!(filein >> nlines_in) || nlines_in < 0) {
      fatal_crash("Ye.txt must start with the count of its cell lines, but the count is unreadable or negative");
    }

    const int last_mgi = get_npts_model() - 1;
    int cells_set = 0;
    bool has_first_or_last_cell = false;
    for (int n = 0; n < nlines_in; n++) {
      int cellnumberin = -1;
      float initelecfrac = 0.;
      if (!(filein >> cellnumberin >> initelecfrac)) {
        fatal_crash("Ye.txt ends or has an unreadable value after {} of its {} cell lines", n, nlines_in);
      }
      // Ye.txt uses the same cell ids as model.txt. read_ejecta_model() detects the id of the first
      // cell (0 or 1) and stores it in first_input_cellid before this function runs.
      const int mgi = cellnumberin - first_input_cellid;
      if (mgi < 0 || mgi > last_mgi) {
        // An out-of-range id is a sign of a cell id mismatch with model.txt. A silently shifted
        // electron fraction would give wrong grey opacities in every cell.
        fatal_crash("Ye.txt: cell {} is outside the model cell range [0..{}]. The cell ids of model.txt start at {}",
                    mgi, last_mgi, first_input_cellid);
      }
      has_first_or_last_cell = has_first_or_last_cell || mgi == 0 || mgi == last_mgi;
      set_initelectronfrac(mgi, initelecfrac);
      cells_set++;
    }
    if (cells_set > 0 && !has_first_or_last_cell) {
      // A file with either boundary id would hit the out-of-range stop above under the other id
      // convention. A sparse file with neither boundary id gives no such confirmation, e.g. a
      // 1-based file for a 0-based model would read with a one-cell shift.
      printlnlog(
          "[warning] Ye.txt holds neither the first cell 0 nor the last cell {}, so the file cannot confirm "
          "that its ids follow model.txt",
          last_mgi);
    }
    printlnlog("Ye.txt: set the initial electron fraction for {} of {} model cells (ids follow model.txt, first id {})",
               cells_set, get_npts_model(), first_input_cellid);
  }
  MPI_Barrier_allranks();
}

// Negative input mass fractions (within roundoff of zero) are counted as they are clamped during
// the model read and summarised once at the end of read_ejecta_model() rather than being logged
// per cell and nuclide, which could produce millions of lines for hydro-derived models.
int n_negative_input_massfrac = 0;
double most_negative_input_massfrac = 0.;
int first_negative_massfrac_mgi = -1;
int first_negative_massfrac_z = -1;  // -1 means the Fe-group mass fraction column
int first_negative_massfrac_a = -1;

// names of model.txt columns that were not recognised, reported after the model read
std::set<std::string> ignored_model_columns;

void count_negative_input_massfrac(const int modelgridindex, const int z, const int a, const double massfrac) {
  n_negative_input_massfrac++;
  most_negative_input_massfrac = std::min(most_negative_input_massfrac, massfrac);
  if (first_negative_massfrac_mgi < 0) {
    first_negative_massfrac_mgi = modelgridindex;
    first_negative_massfrac_z = z;
    first_negative_massfrac_a = a;
  }
}

void report_negative_input_massfracs() {
  if (n_negative_input_massfrac > 0) {
    printlog(
        "[warning] {} negative mass fraction(s) in the input model were clamped to zero (most negative {:g}; first "
        "in cell {} for ",
        n_negative_input_massfrac, most_negative_input_massfrac, first_negative_massfrac_mgi);
    if (first_negative_massfrac_z < 0) {
      printlnlog("X_Fegroup)");
    } else {
      printlnlog("Z={} A={})", first_negative_massfrac_z, first_negative_massfrac_a);
    }
  }
}

void set_ffegrp(const int modelgridindex, float x) {
  if (!(x >= 0.)) {
    if (!(x > -1e-6)) {
      fatal_crash("Fe-group mass fraction {:g} is negative in cell {}", x, modelgridindex);
    }
    count_negative_input_massfrac(modelgridindex, -1, -1, x);
    x = 0.;
  }

  if (x > 1.001) {
    fatal_crash("model.txt gives the Fe-group mass fraction {:g} in cell {}, but a mass fraction must be 1 or less", x,
                modelgridindex);
  }
  modelgrid_input[modelgridindex].ffegrp = x;
}

void set_initenergyq(const int modelgridindex, const float initenergyq) {
  modelgrid_input[modelgridindex].initenergyq = initenergyq;
}

[[nodiscard]] constexpr auto is_otherstable_column(const std::string_view colname) -> bool {
  return colname.starts_with("X_") && colname.ends_with(otherstable_suffix);
}

void parse_model_headerline(const std::string& line, std::vector<int>& zlist, std::vector<int>& alist,
                            std::vector<std::string>& colnames) {
  // custom header line
  std::istringstream iss(line);
  std::string token;

  int columnindex = -1;

  while (iss >> token) {
    columnindex++;

    // The header of a file of format version 2 gives the cell edges, thus such a file has no position column, and
    // rho comes directly after the cell id.
    const bool model_has_edge_lists = (model_format_version == 2);

    // -1 for a column that this model type cannot have
    int required_columnindex = columnindex;
    if (token == "#inputcellid") {
      required_columnindex = 0;
    } else if (token == "velocity_outer" || token == "vel_r_max_kmps") {
      required_columnindex = model_has_edge_lists ? -1 : 1;
    } else if (token.starts_with("pos_") && !model_has_edge_lists) {
      continue;
    } else if (token.starts_with("pos_") || token.starts_with("vel_")) {
      required_columnindex = -1;
    } else if (token == "logrho") {
      // 1D models have log10(rho [g/cm3])
      required_columnindex = (get_modelgridtype() == GridType::SPHERICAL1D && !model_has_edge_lists) ? 2 : -1;
    } else if (token == "rho") {
      // 2D and 3D models have rho [g/cm3] after the cell id and the position columns
      if (model_has_edge_lists) {
        required_columnindex = 1;
      } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
        required_columnindex = -1;
      } else {
        required_columnindex = get_ndim(get_modelgridtype()) + 1;
      }
    } else if (is_otherstable_column(token)) {
      if (!model_has_edge_lists) {
        fatal_crash("the model.txt header line has the column {}, but only format version 2 has such a column", token);
      }
      const int atomic_number =
          decay::get_atomic_number_of_elsymbol(token.substr(2, token.size() - 2 - otherstable_suffix.size()));
      if (atomic_number < 1 || atomic_number > max_atomic_number) {
        fatal_crash("the model.txt header line has the column {}, which does not name an element", token);
      }
      if (std::ranges::find(model_otherstable_atomic_numbers, atomic_number) !=
          model_otherstable_atomic_numbers.end()) {
        fatal_crash("the model.txt header line has two columns {} for the element with Z={}", token, atomic_number);
      }
      model_otherstable_atomic_numbers.push_back(atomic_number);
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    } else if (token.starts_with("X_") && token != "X_Fegroup") {
      const int z = decay::get_nucstring_z(token.substr(2));  // + 2 skips the 'X_'
      const int a = decay::get_nucstring_a(token.substr(2));
      // two names can give one nuclide, e.g. X_Ni56 and X_ni56, and a second column then replaces the mass fraction
      for (auto i = 0Z; i < std::ssize(colnames); i++) {
        if (zlist[i] == z && alist[i] == a) {
          fatal_crash("the model.txt header line has the columns {} and {}, which are the same nuclide (Z={} A={})",
                      colnames[i], token, z, a);
        }
      }
      colnames.push_back(token);
      zlist.push_back(z);
      alist.push_back(a);
    } else {
      if (std::ranges::find(colnames, token) != colnames.end()) {
        fatal_crash("the model.txt header line has the column {} two times", token);
      }
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    }

    if (required_columnindex < 0) {
      fatal_crash(
          "the model.txt header line has the column {}, but a {}D {} model of format version {} cannot have this "
          "column",
          token, get_ndim(get_modelgridtype()), get_grid_type_name(get_modelgridtype()), model_format_version);
    }
    if (columnindex != required_columnindex) {
      fatal_crash("the model.txt header line has the column {} at index {}, but the required index is {}", token,
                  columnindex, required_columnindex);
    }
  }
}

auto get_token_count(std::string const& line) -> int {
  std::string token;
  int abundcolcount = 0;
  auto ssline = std::istringstream{line};
  while (ssline >> token) {
    abundcolcount++;
  }
  return abundcolcount;
}

// The header probes of model.txt read ahead into the data lines. The cell reader takes those lines
// first, so no reader seeks in the file, which a compressed file does not support.
struct ModelFileReader {
  InputFileStream file;
  std::deque<std::string> lines_read_ahead;

  [[nodiscard]] auto getline(std::string& line) -> bool {
    if (!lines_read_ahead.empty()) {
      line = std::move(lines_read_ahead.front());
      lines_read_ahead.pop_front();
      return true;
    }
    return static_cast<bool>(std::getline(file, line));
  }
};

// The limits of a value of a model.txt column, which the reader of a cell with density also applies when it stores the
// value. A mass fraction can be below zero by the rounding of the file. A NaN is outside each limit.
[[nodiscard]] auto is_valid_model_column_value(const std::string& colname, const bool is_tracked_nuclide,
                                               const double value) -> bool {
  constexpr double massfrac_min = -1e-6;
  if (is_tracked_nuclide) {
    return value > massfrac_min && value <= 1.;
  }
  if (colname == "X_Fegroup") {
    return value > massfrac_min && value <= 1.001;
  }
  if (is_otherstable_column(colname)) {
    return value > massfrac_min;
  }
  if (colname == "cellYe" || colname == "Ye") {
    return value >= 0. && value <= 1.001;
  }
  if (colname == "q") {
    return value >= 0.;
  }
  return true;
}

void read_model_radioabundances(ModelFileReader& fmodel, std::string_view& remainder, const int mgi,
                                const bool keepcell, const std::vector<std::string>& colnames,
                                const std::vector<int>& nucindexlist, const bool one_line_per_cell) {
  if (!one_line_per_cell) {
    static std::string line;
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends before the second line of cell {}", mgi);
    }
    remainder = std::string_view{line};
  }

  // A line of an empty cell of format version 2 must also have one valid value for each column. Thus a line that lost
  // its end or that holds a value outside the limits stops the run. The reader stores no value of an empty cell.
  if (!keepcell && model_format_version != 2) {
    return;
  }

  const auto num_otherstable_columns = std::ssize(model_otherstable_atomic_numbers);
  ptrdiff_t otherstable_column = 0;
  for (auto i = 0Z; i < std::ssize(colnames); i++) {
    double valuein = 0.;  // a mass fraction or another column value, e.g. Ye or q
    if (!parse_next_token(remainder, valuein)) {
      fatal_crash("model.txt cell {} has no readable value for the column {}", mgi, colnames[i]);
    }
    if (!keepcell) {
      if (!is_valid_model_column_value(colnames[i], nucindexlist[i] >= 0, valuein)) {
        fatal_crash(
            "model.txt cell {} has the value {:g} in the column {}, which is outside the limits of that column. The "
            "line of a cell with a density of zero must also hold valid values",
            mgi, valuein, colnames[i]);
      }
      continue;
    }

    if (nucindexlist[i] >= 0) {
      set_modelinitnucmassfrac(mgi, nucindexlist[i], static_cast<float>(valuein));
    } else if (is_otherstable_column(colnames[i])) {
      // the columns have the sequence of model_otherstable_atomic_numbers
      if (valuein < 0.) {
        if (valuein <= -1e-6) {
          fatal_crash("model.txt cell {} has the negative value {:g} in the column {}", mgi, valuein, colnames[i]);
        }
        valuein = 0.;
      }
      model_otherstable_massfracs[(static_cast<ptrdiff_t>(mgi) * num_otherstable_columns) + otherstable_column] =
          static_cast<float>(valuein);
      otherstable_column++;
    } else if (colnames[i] == "X_Fegroup") {
      set_ffegrp(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "cellYe" || colnames[i] == "Ye") {
      set_initelectronfrac(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "q") {
      // The q column holds the trapped radiation energy per mass at t_model. It already includes the adiabatic
      // losses before t_model. The energy of a comoving mass element falls as 1/t, so the value is scaled to tmin.
      if (!(valuein >= 0.)) {
        fatal_crash("model.txt cell {} has the value {:g} in the column q, but the value must be zero or positive", mgi,
                    valuein);
      }
      set_initenergyq(mgi, static_cast<float>(valuein * t_model / globals::tmin));
    } else if (colnames[i] == "tracercount") {
      ;
    } else {
      // reported once after the model read (checking only mgi == 0 could miss it entirely if cell 0 is empty)
      ignored_model_columns.insert(colnames[i]);
    }
  }
  double valuein = 0.;
  if (parse_next_token(remainder, valuein)) {
    fatal_crash("model.txt cell {} has more values than the header has columns", mgi);
  }
  // a cell line of format version 2 has only whitespace after its values. Text there shows a damaged line
  remainder = trim_token_whitespace(remainder);
  if (model_format_version == 2 && !remainder.empty()) {
    fatal_crash("model.txt cell {} has the text '{}' after the values of the columns", mgi, remainder);
  }
}

auto read_model_columns(ModelFileReader& fmodel) -> std::tuple<std::vector<std::string>, std::vector<int>, bool> {
  std::vector<int> zlist;
  std::vector<int> alist;
  std::vector<std::string> colnames;

  std::string line;
  if (!fmodel.getline(line)) {
    fatal_crash("model.txt ends before the column header or the first cell");
  }

  std::string headerline;

  const bool header_specified = lineiscommentonly(line);

  if (header_specified) {
    // line is the header
    headerline = line;
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends after the column header, before the first cell");
    }
    if (model_format_version == 2) {
      // a comment line or an empty line can come before each cell line of format version 2, also before the first one
      while (lineiscommentonly(line)) {
        if (!fmodel.getline(line)) {
          fatal_crash("model.txt ends after the key columns, before the first cell");
        }
      }
    }
  } else {
    // line is not a comment, so it must be the first line of data
    // add a default header for unlabelled columns
    switch (get_modelgridtype()) {
      case GridType::SPHERICAL1D:
        headerline = "#inputcellid vel_r_max_kmps logrho";
        break;
      case GridType::CYLINDRICAL2D:
        headerline = "#inputcellid pos_rcyl_mid pos_z_mid rho";
        break;
      case GridType::CARTESIAN3D:
        headerline = "#inputcellid pos_x_min pos_y_min pos_z_min rho";
        break;
      case GridType::SPHERICAL2D:
      case GridType::SPHERICAL3D:
        fatal_crash("model.txt: a 2D spherical or a 3D spherical model needs format version 2");
    }
    headerline += " X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48";
  }

  int colcount = get_token_count(line);
  const bool one_line_per_cell = (colcount >= get_token_count(headerline));

  printlnlog("model.txt has {} line per cell format", one_line_per_cell ? "one" : "two");
  if (!one_line_per_cell && model_format_version == 2) {
    fatal_crash(
        "model.txt: the first cell has {} values, but the key columns has {} names. A cell has one line in format "
        "version 2",
        colcount, get_token_count(headerline));
  }

  std::string secondline;
  if (!one_line_per_cell) {  // add columns from the second line
    if (!fmodel.getline(secondline)) {
      fatal_crash("model.txt ends before the second line of the first cell");
    }
    colcount += get_token_count(secondline);
  }

  // the cell reader takes the lines of the first cell again
  fmodel.lines_read_ahead.push_back(line);
  if (!one_line_per_cell) {
    fmodel.lines_read_ahead.push_back(secondline);
  }

  if (!header_specified && colcount > get_token_count(headerline)) {
    headerline += " X_Ni57 X_Co57";
  }

  if (colcount != get_token_count(headerline)) {
    fatal_crash("model.txt cell 0 has {} values, but the {} '{}' has {} columns", colcount,
                header_specified ? "header line" : "default header for a file without a header line", headerline,
                get_token_count(headerline));
  }

  if (header_specified) {
    printlnlog("model.txt has a header line.");
  } else {
    printlnlog("model.txt has no header line. Using default: {}", headerline);
  }

  parse_model_headerline(headerline, zlist, alist, colnames);

  if (std::ranges::none_of(colnames, [](const std::string& colname) { return colname == "X_Fegroup"; })) {
    printlnlog(
        "model.txt has no column X_Fegroup, thus the Fe-group mass fraction of each cell is 0. This mass fraction "
        "sets the photoelectric opacity and the pair production opacity of the gamma rays, and the grey opacity of "
        "FEGROUP_APPROX");
  }

  // a column that is not a nuclide has Z = -1, and the decay network does not hold the neutron (Z = 0)
  std::vector<int> nuclide_columns_z;
  std::vector<int> nuclide_columns_a;
  for (const auto [z, a] : std::views::zip(zlist, alist)) {
    if (z > 0) {
      nuclide_columns_z.push_back(z);
      nuclide_columns_a.push_back(a);
    }
  }
  decay::init_nuclides(nuclide_columns_z, nuclide_columns_a);

  std::vector<int> nucindexlist(zlist.size());
  for (auto i = 0Z; i < std::ssize(zlist); i++) {
    nucindexlist[i] = (zlist[i] > 0) ? decay::get_nucindex(zlist[i], alist[i]) : -1;
  }

  assert_always(npts_model > 0);

  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells = MPI_shared_array<float>(npts_model * num_nuclides, 0.);
  if (globals::rank_in_node == 0) {
    model_otherstable_massfracs.assign(npts_model * std::ssize(model_otherstable_atomic_numbers), 0.F);
  }
  printlnlog(
      "[info] mem_usage: input abundance data for {} nuclides for {} cells occupies {:.3f} MB (node shared memory)",
      num_nuclides, npts_model, (initnucmassfrac_allcells.size() * sizeof(float)) / 1024. / 1024.);

  return {colnames, nucindexlist, one_line_per_cell};
}

// The largest difference [rad] between the last azimuthal edge and 2 pi that the reader of model.txt replaces by the
// exact angle. Each cell of cos(theta) and each cell of phi must also be wider than this number. The first polar edge
// away from the z axis is then more than 1.4e-4 rad from the axis. Near the axis, a cosine gives the position of its
// cone to approximately 1.1e-16 / theta of the radius, and the position tolerance of a cell boundary is 1e-12 of the
// radius.
constexpr double model_angle_edge_tolerance = 1e-8;

// the largest difference of a cylindrical or a Cartesian cell edge from its place on a grid with cells of equal width,
// as a fraction of vmax. The model grid of these two types has cells of equal width only
constexpr double model_uniform_edge_tolerance = 1e-9;

// Read the header of a model.txt of format version 2, after its first line "format_version 2". Each header line is
// a key and its values, and the keys can come in any sequence. The key columns gives the column names and ends the
// header. The other keys are coordinates, t_model_init_days, and one list of cell edges for each axis of the
// coordinates, e.g. vel_r_edges_kmps. A list has one value more than the count of the cells on its axis. The optional
// key cells_per_axis gives these counts again, in the sequence of the axes, e.g. "5 4 4". A key that
// this version of sn3d does not know gets a warning, thus a later version of sn3d can add a key to format version 2.
// The function gives the column names back to the reader as a header line of format version 1.
auto read_model_keyword_header(ModelFileReader& fmodel) -> GridType {
  using KeyLine = std::pair<std::string, std::string>;  // a key and the text after it, without a comment
  std::vector<KeyLine> keylines;
  const auto has_key = [&keylines](const std::string_view key) {
    return std::ranges::contains(keylines, key, &KeyLine::first);
  };
  std::string line;
  while (!has_key("columns")) {
    if (!std::getline(fmodel.file, line)) {
      fatal_crash("model.txt: the header must end with the key columns, which gives the column names");
    }
    if (lineiscommentonly(line)) {
      continue;
    }
    const auto text = trim_token_whitespace(std::string_view{line}.substr(0, line.find('#')));
    const auto key = text.substr(0, text.find_first_of(token_whitespace));
    if (has_key(key)) {
      fatal_crash("model.txt: the header has the key {} two times", key);
    }
    keylines.emplace_back(std::string{key}, std::string{trim_token_whitespace(text.substr(key.size()))});
  }

  const auto get_text_of_key = [&keylines](const std::string_view key) -> std::string_view {
    const auto keyline = std::ranges::find(keylines, key, &KeyLine::first);
    if (keyline != keylines.end()) {
      return keyline->second;
    }
    fatal_crash("model.txt: the header has no key {}", key);
  };

  const auto get_numbers_of_key = [&](const std::string_view key) {
    auto remainder = get_text_of_key(key);
    std::vector<double> numbers;
    double number{NAN};
    while (parse_next_token(remainder, number)) {
      numbers.push_back(number);
    }
    remainder = trim_token_whitespace(remainder);
    if (!remainder.empty()) {
      fatal_crash("model.txt: the key {} has an unreadable number at '{}'", key, remainder);
    }
    return numbers;
  };

  // An edge list increases. A list of r, rcyl, or phi starts at exactly zero, and a list of cos(theta) starts at
  // exactly -1, or at 0 with a mirror plane. An angular list ends at last_edge_exact: a list of cos(theta) ends at
  // exactly 1, and a list of phi ends at 2 pi within last_edge_tolerance
  const auto require_valid_edges = [](std::vector<double>& edges, const std::string_view key,
                                      const std::optional<double> first_edge_exact,
                                      const std::optional<double> last_edge_exact,
                                      const double last_edge_tolerance = 0.) {
    if (edges.size() < 2) {
      fatal_crash("model.txt: the key {} must have two or more values", key);
    }
    if (first_edge_exact.has_value() && edges.front() != first_edge_exact.value()) {
      fatal_crash("model.txt: the first value of {} is {:.17g}, but it must be {:g}", key, edges.front(),
                  first_edge_exact.value());
    }
    if (last_edge_exact.has_value()) {
      if (!(std::abs(edges.back() - last_edge_exact.value()) <= last_edge_tolerance)) {
        fatal_crash("model.txt: the last value of {} is {:.17g}, but it must be {:.17g}", key, edges.back(),
                    last_edge_exact.value());
      }
      edges.back() = last_edge_exact.value();
    }
    for (auto i = 1Z; i < std::ssize(edges); i++) {
      if (!(edges[i] > edges[i - 1])) {
        fatal_crash("model.txt: the values of {} must increase, but {:.17g} follows {:.17g}", key, edges[i],
                    edges[i - 1]);
      }
      // The position tolerance of a cell boundary is 1e-12 of the radius (see cellbound_tolerance). A packet in an
      // angular cell that is thinner than this tolerance goes through the cell with no path length.
      if (last_edge_exact.has_value() && !(edges[i] - edges[i - 1] > model_angle_edge_tolerance)) {
        fatal_crash(
            "model.txt: the values {:.17g} and {:.17g} of {} are the edges of a cell that is too thin. Each cell of "
            "cos(theta) and of phi [rad] must be wider than {:g}",
            edges[i - 1], edges[i], key, model_angle_edge_tolerance);
      }
    }
  };

  // the cell edges of a velocity axis [cm/s]
  const auto get_velocity_edges = [&](const std::string_view key, const bool first_edge_is_zero) {
    auto edges = get_numbers_of_key(key);
    require_valid_edges(edges, key, first_edge_is_zero ? std::optional<double>{0.} : std::nullopt, std::nullopt);
    for (auto& edge : edges) {
      edge *= 1.e5;  // km/s to cm/s
    }
    // two values in km/s that differ in the last bit can give the same velocity in cm/s
    for (auto i = 1Z; i < std::ssize(edges); i++) {
      if (!(edges[i] > edges[i - 1])) {
        fatal_crash(
            "model.txt: the values of {} at the indices {} and {} give the velocities {:.17g} and {:.17g} [cm/s]. The "
            "values must also increase in cm/s",
            key, i - 1, i, edges[i - 1], edges[i]);
      }
    }
    return edges;
  };

  // sn3d holds a cylindrical grid and a Cartesian grid with cells of equal width from first_edge to vmax [cm/s]
  const auto require_uniform_edges = [](const std::vector<double>& edges, const std::string_view key,
                                        const double first_edge, const double vmax) {
    const auto num_cells = std::ssize(edges) - 1;
    for (auto i = 0Z; i <= num_cells; i++) {
      const double uniform_edge = first_edge + ((vmax - first_edge) * static_cast<double>(i) / num_cells);
      if (!(std::abs(edges[i] - uniform_edge) <= model_uniform_edge_tolerance * vmax)) {
        fatal_crash(
            "model.txt: the values of {} must have equal steps from {:g} to {:g} [km/s], but the value at index {} is "
            "{:.17g}. sn3d reads a cylindrical grid or a Cartesian grid only with cells of equal width and the same "
            "vmax on each axis",
            key, first_edge / 1.e5, vmax / 1.e5, i, edges[i] / 1.e5);
      }
    }
  };

  const auto coordinates = get_text_of_key("coordinates");
  // the name of each coordinate type, with the keys of its edge lists in the sequence of its axes
  using CoordinatesSpec = std::tuple<std::string_view, GridType, std::array<std::string_view, 3>>;
  constexpr std::array<CoordinatesSpec, 5> coordinates_specs{
      CoordinatesSpec{"spherical1d", GridType::SPHERICAL1D, {"vel_r_edges_kmps", "", ""}},
      CoordinatesSpec{"cylindrical2d", GridType::CYLINDRICAL2D, {"vel_rcyl_edges_kmps", "vel_z_edges_kmps", ""}},
      CoordinatesSpec{
          "cartesian3d", GridType::CARTESIAN3D, {"vel_x_edges_kmps", "vel_y_edges_kmps", "vel_z_edges_kmps"}},
      CoordinatesSpec{"spherical2d", GridType::SPHERICAL2D, {"vel_r_edges_kmps", "costheta_edges", ""}},
      CoordinatesSpec{
          "spherical3d", GridType::SPHERICAL3D, {"vel_r_edges_kmps", "costheta_edges", "phi_edges_radians"}},
  };
  std::optional<GridType> named_model_type;
  std::vector<std::string_view> edge_keys_of_model;
  for (const auto& [name, gridtype, edge_keys_of_type] : coordinates_specs) {
    if (name == coordinates) {
      named_model_type = gridtype;
      std::ranges::copy_if(edge_keys_of_type, std::back_inserter(edge_keys_of_model),
                           [](const std::string_view key) { return !key.empty(); });
    }
  }
  if (!named_model_type.has_value()) {
    fatal_crash(
        "model.txt: the coordinates '{}' are not known. The choices are spherical1d, cylindrical2d, cartesian3d, "
        "spherical2d, and spherical3d",
        coordinates);
  }
  const auto detected_model_type = named_model_type.value();

  model_mirror_z = false;
  if (has_key("mirror_plane")) {
    if (get_text_of_key("mirror_plane") != "z=0") {
      fatal_crash("model.txt: the key mirror_plane must have the value z=0, but it has '{}'",
                  get_text_of_key("mirror_plane"));
    }
    if (detected_model_type != GridType::CYLINDRICAL2D && detected_model_type != GridType::SPHERICAL2D &&
        detected_model_type != GridType::SPHERICAL3D) {
      fatal_crash(
          "model.txt: a {} model takes no key mirror_plane. The key applies to cylindrical2d, spherical2d, and "
          "spherical3d models",
          coordinates);
    }
    model_mirror_z = true;
    printlnlog(
        "model.txt has mirror symmetry in the plane z = 0. Each model cell fills its cell with z >= 0 and the mirror "
        "cell with z <= 0");
  }
  printlnlog("Detected {} model of format version 2", coordinates);

  const auto t_model_days = get_numbers_of_key("t_model_init_days");
  if (t_model_days.size() != 1 || !std::isfinite(t_model_days[0]) || t_model_days[0] <= 0.) {
    fatal_crash("model.txt: the key t_model_init_days must have one value above zero, but it has '{}'",
                get_text_of_key("t_model_init_days"));
  }
  t_model = t_model_days[0] * DAY;
  if (t_model > globals::tmin) {
    fatal_crash("The model snapshot time {} d is after tmin {} d", t_model_days[0], globals::tmin / DAY);
  }

  constexpr std::array<std::string_view, 7> edge_keys{
      "vel_r_edges_kmps", "costheta_edges",   "phi_edges_radians", "vel_rcyl_edges_kmps",
      "vel_x_edges_kmps", "vel_y_edges_kmps", "vel_z_edges_kmps",
  };
  for (const auto& keyline : keylines) {
    const std::string_view key = keyline.first;
    if (key == "coordinates" || key == "cells_per_axis" || key == "mirror_plane" || key == "t_model_init_days" ||
        key == "cell_line_count" || key == "columns") {
      continue;
    }
    if (std::ranges::find(edge_keys, key) == edge_keys.end()) {
      printlnlog("[warning] model.txt: this version of sn3d does not know the key {}, and it ignores the key", key);
    } else if (std::ranges::find(edge_keys_of_model, key) == edge_keys_of_model.end()) {
      fatal_crash("model.txt: a {} model takes no key {}", coordinates, key);
    }
  }

  switch (detected_model_type) {
    case GridType::SPHERICAL1D: {
      model_vel_r_edges = get_velocity_edges("vel_r_edges_kmps", true);
      globals::vmax = model_vel_r_edges.back();
      // a version 1 file gives the outer velocity of each shell in a column
      vout_model.assign(model_vel_r_edges.begin() + 1, model_vel_r_edges.end());
      ncoord_model = {static_cast<int>(std::ssize(vout_model)), 0, 0};
      npts_model = std::ssize(vout_model);
      break;
    }

    case GridType::CYLINDRICAL2D: {
      const auto vel_rcyl_edges = get_velocity_edges("vel_rcyl_edges_kmps", true);
      const auto vel_z_edges = get_velocity_edges("vel_z_edges_kmps", false);
      globals::vmax = vel_rcyl_edges.back();
      require_uniform_edges(vel_rcyl_edges, "vel_rcyl_edges_kmps", 0., globals::vmax);
      // with a mirror plane, the z axis of the file starts at the plane z = 0
      require_uniform_edges(vel_z_edges, "vel_z_edges_kmps", model_mirror_z ? 0. : -globals::vmax, globals::vmax);
      ncoord_model = {
          static_cast<int>(std::ssize(vel_rcyl_edges)) - 1,
          static_cast<int>(std::ssize(vel_z_edges)) - 1,
          0,
      };
      npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1];
      break;
    }

    case GridType::CARTESIAN3D: {
      for (int axis = 0; axis < 3; axis++) {
        const auto vel_edges = get_velocity_edges(edge_keys_of_model[axis], false);
        if (axis == 0) {
          globals::vmax = vel_edges.back();
        }
        require_uniform_edges(vel_edges, edge_keys_of_model[axis], -globals::vmax, globals::vmax);
        ncoord_model[axis] = static_cast<int>(std::ssize(vel_edges)) - 1;
      }
      if (ncoord_model[1] != ncoord_model[0] || ncoord_model[2] != ncoord_model[0]) {
        fatal_crash(
            "model.txt: the lists vel_x_edges_kmps, vel_y_edges_kmps, and vel_z_edges_kmps give {} x {} x {} cells, "
            "but sn3d reads a Cartesian grid with the same count of cells on each axis",
            ncoord_model[0], ncoord_model[1], ncoord_model[2]);
      }
      npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * ncoord_model[2];
      break;
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      model_vel_r_edges = get_velocity_edges("vel_r_edges_kmps", true);
      globals::vmax = model_vel_r_edges.back();

      // the polar axis of a model with a mirror plane starts at the plane z = 0
      model_costheta_edges = get_numbers_of_key("costheta_edges");
      require_valid_edges(model_costheta_edges, "costheta_edges", model_mirror_z ? 0. : -1., 1.);

      if (detected_model_type == GridType::SPHERICAL3D) {
        model_phi_edges = get_numbers_of_key("phi_edges_radians");
        require_valid_edges(model_phi_edges, "phi_edges_radians", 0., 2 * PI, model_angle_edge_tolerance);
      } else {
        model_phi_edges = {0., 2 * PI};
      }
      set_spherical_model_cell_counts(detected_model_type);
      break;
    }
  }
  if (npts_model > std::numeric_limits<int>::max()) {
    fatal_crash("model.txt: the cell edges give {} cells, which is more than the largest cell index", npts_model);
  }

  // An edge list that lost a value or got one more value gives a different number of cells on its axis. The key is
  // optional.
  if (has_key("cells_per_axis")) {
    const auto cells_per_axis = get_numbers_of_key("cells_per_axis");
    const auto num_axes = std::ssize(edge_keys_of_model);
    std::string edge_cell_counts;
    bool counts_agree = (std::ssize(cells_per_axis) == num_axes);
    for (int axis = 0; axis < num_axes; axis++) {
      edge_cell_counts += std::format("{}{}", (axis > 0) ? " " : "", ncoord_model[axis]);
      counts_agree = counts_agree && (cells_per_axis[axis] == ncoord_model[axis]);
    }
    if (!counts_agree) {
      fatal_crash("model.txt: the key cells_per_axis gives the cell counts '{}', but the cell edges give '{}'",
                  get_text_of_key("cells_per_axis"), edge_cell_counts);
    }
  }

  model_cell_line_count.reset();
  if (has_key("cell_line_count")) {
    const auto counts = get_numbers_of_key("cell_line_count");
    if (counts.size() != 1 || counts[0] != std::floor(counts[0]) || counts[0] < 1 ||
        counts[0] > static_cast<double>(npts_model)) {
      fatal_crash(
          "model.txt: the key cell_line_count must have one integer from 1 to {}, which is the count of the cells "
          "that the cell edges give, but its text is '{}'",
          npts_model, get_text_of_key("cell_line_count"));
    }
    model_cell_line_count = static_cast<ptrdiff_t>(counts[0]);
  }

  const auto columns = get_text_of_key("columns");
  std::string first_column;
  std::string second_column;
  std::istringstream{std::string{columns}} >> first_column >> second_column;
  if (first_column != "inputcellid" || second_column != "rho") {
    fatal_crash("model.txt: the first two values of the key columns must be inputcellid and rho, but they are '{}'",
                columns);
  }
  // read_model_columns() takes the column names as a header line of format version 1
  fmodel.lines_read_ahead.push_back(std::format("#{}", columns));

  return detected_model_type;
}

}  // anonymous namespace

void set_modelinitnucmassfrac(const int modelgridindex, const int nucindex, float abund) {
  // set the mass fraction of a nuclide in a model grid cell at t=t_model by nuclide index
  // initnucmassfrac array is in node shared memory
  assert_always(nucindex >= 0);
  if (!(abund >= 0.)) {
    if (!(abund > -1e-6)) {
      fatal_crash("mass fraction {:g} for Z={} A={} is negative in cell {}", abund, decay::get_nuc_z(nucindex),
                  decay::get_nuc_a(nucindex), modelgridindex);
    }
    count_negative_input_massfrac(modelgridindex, decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex), abund);
    abund = 0.;
  }

  if (abund > 1.) {
    fatal_crash("the mass fraction {:.9g} of Z={} A={} in cell {} is above 1", abund, decay::get_nuc_z(nucindex),
                decay::get_nuc_a(nucindex), modelgridindex);
  }
  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells[(modelgridindex * num_nuclides) + nucindex] = abund;
}

// Read the element mass fractions of each model cell with a density above zero. A cell with no mass in the
// elements of compositiondata.txt gets zero density, so that the mapping makes it an empty cell. Return the mass
// fractions in model cell order (an empty vector on the ranks other than the node leader) and the number of
// cells that became empty.
auto read_elem_abundances() -> std::tuple<std::vector<float>, int> {
  // a model.txt of format version 2 gives the element mass fractions, and such a model has no abundances.txt
  const bool elements_from_model = (model_format_version == 2);
  const std::string_view massfrac_source = elements_from_model ? "model.txt" : "abundances.txt";
  if (elements_from_model && inputfile_exists("abundances.txt")) {
    fatal_crash(
        "abundances.txt is present, but model.txt has format version 2. Such a model.txt gives the mass fraction of "
        "an element as the sum of its nuclide columns and its column X_<element>_otherstable. Remove abundances.txt");
  }
  if (elements_from_model) {
    printlnlog("calculating the element mass fractions from the columns of model.txt...");
  } else {
    printlnlog("reading abundances.txt...");
  }
  const int nelements = get_nelements();

  // abundances.txt holds one line for each model cell, also for an empty cell. Each line gives the cell number and
  // then the mass fraction of each element from Z=1 upwards. The line can end after the last element that the
  // model needs.

  // the density is in node-shared memory, so only the node leader of each node parses the file and
  // writes the values (synchronised by the barrier below). The other ranks would just discard everything they read
  std::vector<float> elem_massfracs_of_mgi;
  int ncells_emptied = 0;
  if (globals::rank_in_node == 0) {
    elem_massfracs_of_mgi.assign(static_cast<size_t>(get_npts_model()) * nelements, 0.F);
    std::set<int> atomic_numbers_without_atomic_data;  // elements with mass in the cells that became empty
    // Every log line gets a timestamp and a flush, so a large 3D model must not warn per cell.
    constexpr int max_warnings = 10;

    // Check the element mass fractions of one cell with density, and store those of the elements with atomic data.
    // elem_massfracs_in holds the mass fraction of each element from Z=1 upwards.
    const auto store_cell_elem_massfracs = [&](const int mgi,
                                               const std::array<float, max_atomic_number>& elem_massfracs_in) {
      // a sum above one is not physical. The limit permits the rounding of the values in the file.
      const double elem_massfrac_sum = std::accumulate(elem_massfracs_in.begin(), elem_massfracs_in.end(), 0.);
      if (elem_massfrac_sum > 1. + 1e-3) {
        fatal_crash("{}: cell {} has element mass fractions that sum to {:g}, which is more than 1", massfrac_source,
                    mgi, elem_massfrac_sum);
      }

      bool has_mass_in_composition_elements = false;
      const auto tracked_isotope_massfrac_sum_by_atomic_number = get_tracked_isotope_massfrac_sum_by_atomic_number(mgi);

      for (int element = 0; element < nelements; element++) {
        // set the mass fraction of each element that the atomic data includes
        const int atomic_number = get_atomicnumber(element);
        const float elemmassfrac = elem_massfracs_in[atomic_number - 1];
        assert_always(elemmassfrac >= 0.);

        // read_ejecta_model() has already set the nuclide mass fractions. Check here, while
        // the mass fractions are still exactly as the input files gave them, that abundances.txt gives the
        // element at least as much mass as model.txt gives its tracked isotopes. Any later mapping rescale
        // changes the nuclide mass fractions but not the elemental ones, so this cannot be tested afterwards.
        // a small negative remainder is allowed for roundoff error
        const double tracked_isotope_massfrac_sum = tracked_isotope_massfrac_sum_by_atomic_number[atomic_number];
        if ((elemmassfrac - tracked_isotope_massfrac_sum) < -1e-2) {
          fatal_crash(
              "{} cell {} gives the mass fraction {:g} for Z={}, but the nuclide columns of model.txt give this "
              "element the larger mass fraction {:g}",
              massfrac_source, mgi, elemmassfrac, atomic_number, tracked_isotope_massfrac_sum);
        }

        elem_massfracs_of_mgi[(static_cast<size_t>(mgi) * nelements) + element] = elemmassfrac;
        has_mass_in_composition_elements =
            has_mass_in_composition_elements ||
            (get_nions(element) > 0 && (elemmassfrac > 0.F || tracked_isotope_massfrac_sum > 0.));
      }

      if (!has_mass_in_composition_elements) {
        // e.g. a mapped model gives zero abundances to a cell whose mass comes only from particles without
        // nucleosynthesis data. Matter without the elements of compositiondata.txt has no opacity and no cooling
        // processes.
        ncells_emptied++;
        // every element with mass in this cell has no atomic data, so list them for the user
        std::string elements_without_atomic_data;
        for (int elem_z_index = 0; elem_z_index < std::ssize(elem_massfracs_in); elem_z_index++) {
          if (elem_massfracs_in[elem_z_index] > 0.F) {
            atomic_numbers_without_atomic_data.insert(elem_z_index + 1);
            elements_without_atomic_data +=
                std::format(" Z={} ({:.2e})", elem_z_index + 1, elem_massfracs_in[elem_z_index]);
          }
        }
        if (ncells_emptied <= max_warnings) {
          printlnlog(
              "[warning] read_elem_abundances: cell {} has density {:g} [g/cm3] at tmin and no mass in the "
              "elements of compositiondata.txt. The cell becomes empty. Elements with mass and no atomic data:{}",
              mgi, get_rho_tmin(mgi), elements_without_atomic_data.empty() ? " none" : elements_without_atomic_data);
        }
        set_rho_tmin(mgi, 0.);
      }
    };

    if (elements_from_model) {
      // The mass fraction of an element is the sum of its nuclide columns of model.txt and of its column
      // X_<element>_otherstable, which holds the stable isotopes that have no column of their own.
      const auto num_otherstable_columns = std::ssize(model_otherstable_atomic_numbers);
      for (int mgi = 0; mgi < get_npts_model(); mgi++) {
        if (get_rho_tmin(mgi) <= 0.) {
          continue;
        }
        const auto tracked_isotope_massfrac_sum_by_atomic_number =
            get_tracked_isotope_massfrac_sum_by_atomic_number(mgi);
        std::array<double, max_atomic_number> elem_massfracs_double{};
        for (int atomic_number = 1; atomic_number <= max_atomic_number; atomic_number++) {
          elem_massfracs_double[atomic_number - 1] = tracked_isotope_massfrac_sum_by_atomic_number[atomic_number];
        }
        for (auto column = 0Z; column < num_otherstable_columns; column++) {
          elem_massfracs_double[model_otherstable_atomic_numbers[column] - 1] +=
              model_otherstable_massfracs[(static_cast<ptrdiff_t>(mgi) * num_otherstable_columns) + column];
        }
        std::array<float, max_atomic_number> elem_massfracs_in{};
        for (int elem_z_index = 0; elem_z_index < max_atomic_number; elem_z_index++) {
          elem_massfracs_in[elem_z_index] = static_cast<float>(elem_massfracs_double[elem_z_index]);
        }
        store_cell_elem_massfracs(mgi, elem_massfracs_in);
      }
      model_otherstable_massfracs.clear();
      model_otherstable_massfracs.shrink_to_fit();
    } else {
      auto abundance_file = istream_required("abundances.txt");
      std::string line;
      for (int mgi = 0; mgi < get_npts_model(); mgi++) {
        if (!get_noncommentline(abundance_file, line)) {
          fatal_crash("abundances.txt ends before cell {} of the {} cells of model.txt", mgi, get_npts_model());
        }
        auto remainder = std::string_view{line};

        int cellnumberinput = -1;
        if (!parse_next_token(remainder, cellnumberinput) || cellnumberinput != mgi + first_input_cellid) {
          fatal_crash("abundances.txt cell {} must start with the cell id {} as in model.txt, but the line is '{}'",
                      mgi, mgi + first_input_cellid, line);
        }

        // the abundances.txt file specifies the elemental mass fractions for each model cell.
        // The values are true mass fractions and get no normalisation. Their sum can be below one.
        // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
        std::array<float, max_atomic_number> elem_massfracs_in{};
        double abund_in = 0.;
        for (int elem_z_index = 0; elem_z_index < std::ssize(elem_massfracs_in); elem_z_index++) {
          if (!parse_next_token(remainder, abund_in)) {
            // a row with no values is a cell with no elements
            break;
          }

          if (abund_in < std::numeric_limits<float>::min()) {
            if (!(abund_in > -1e-6)) {
              fatal_crash("abundances.txt cell {} gives the negative mass fraction {:g} for Z={}", mgi, abund_in,
                          elem_z_index + 1);
            }
            abund_in = 0.;
          }
          elem_massfracs_in[elem_z_index] = static_cast<float>(abund_in);
        }

        // parse_next_token() gives false at the end of the row and also for a token that is not a
        // number, e.g. "nan". The loop above stops on both, so a bad token would silently zero every
        // later element. Only whitespace may remain here.
        remainder = trim_token_whitespace(remainder);
        if (!remainder.empty()) {
          fatal_crash("read_elem_abundances: cell {} has an unreadable token at '{}'", mgi, remainder);
        }

        if (get_rho_tmin(mgi) > 0.) {
          store_cell_elem_massfracs(mgi, elem_massfracs_in);
        }
      }
    }

    if (!atomic_numbers_without_atomic_data.empty()) {
      std::string atomic_numbers_text;
      for (const int atomic_number : atomic_numbers_without_atomic_data) {
        atomic_numbers_text += std::format(" {}", atomic_number);
      }
      printlnlog(
          "[warning] read_elem_abundances: add atomic data for these elements to keep the mass of the cells that "
          "became empty. Atomic numbers:{}",
          atomic_numbers_text);
    }
  }

  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier_allranks();
  MPI_Bcast_safe(ncells_emptied, 0, globals::mpi_comm_node);
  printlnlog("finished setting the element mass fractions from {}", massfrac_source);
  return {std::move(elem_massfracs_of_mgi), ncells_emptied};
}

// Set the cell counts of a 2D spherical or a 3D spherical model from its cell edges. The third count of a 2D model
// is zero, as for a 2D cylindrical model.
void set_spherical_model_cell_counts(const GridType spherical_model_type) {
  ncoord_model = {
      static_cast<int>(std::ssize(model_vel_r_edges)) - 1,
      static_cast<int>(std::ssize(model_costheta_edges)) - 1,
      (spherical_model_type == GridType::SPHERICAL3D) ? static_cast<int>(std::ssize(model_phi_edges)) - 1 : 0,
  };
  npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * std::max(1, ncoord_model[2]);
}

// read the input ejecta model from model.txt, auto-detecting its format version and its dimensionality (1D spherical,
// 2D cylindrical, 3D Cartesian, 2D spherical, or 3D spherical) and reading the per-cell densities,
// abundances, and any optional extra columns into the model grid
void read_ejecta_model() {
  auto fmodel = ModelFileReader{.file = istream_required("model.txt"), .lines_read_ahead = {}};
  std::string line;
  std::optional<GridType> detected_model_type{};

  // two integers on the first line of the model file
  int npts_0 = 0;  // total model points for 1D/3D, and number of points in r for 2D
  int npts_1 = 0;  // number of points in z for 2D
  // The first line of a file of format version 2 is "format_version 2". Such a file has a header of keys in place
  // of the header of numbers. A file of version 1 has no version line.
  if (!std::getline(fmodel.file, line)) {
    fatal_crash("model.txt is empty");
  }
  // the version line can end with a comment
  const auto firstline = trim_token_whitespace(std::string_view{line}.substr(0, line.find('#')));
  std::string firstword{firstline.substr(0, firstline.find_first_of(token_whitespace))};
  model_format_version = (firstword == "format_version") ? 2 : 1;
  if (model_format_version == 2) {
    // the version is the only text after the key
    const auto versiontext = trim_token_whitespace(firstline.substr(firstword.size()));
    if (versiontext != "2") {
      fatal_crash(
          "model.txt has the format version '{}'. This version of sn3d reads format version 2, and format version 1, "
          "which has no format_version line",
          versiontext);
    }
    detected_model_type = read_model_keyword_header(fmodel);
  } else {
    if (lineiscommentonly(line)) {
      if (!get_noncommentline(fmodel.file, line)) {
        fatal_crash("model.txt holds only comment lines. Its first line with a value must give the cell count");
      }
      firstword.clear();
      std::istringstream{line} >> firstword;
    }
    if (firstword == "format_version") {
      fatal_crash("model.txt: the line '{}' must be the first line of the file", line);
    }
    if (firstword == "coordinates") {
      fatal_crash("model.txt: the key coordinates needs the line 'format_version 2' as the first line of the file");
    }
    auto ssline = std::istringstream{line};
    ssline >> npts_0;
    if (npts_0 <= 0) {
      fatal_crash("model.txt: could not read a positive cell count from the first line '{}'", line);
    }
    if (ssline >> npts_1) {
      // second number on the line for 2D means the line was n_r n_z
      detected_model_type = GridType::CYLINDRICAL2D;
      printlnlog("Detected 2D model");
      npts_model = static_cast<ptrdiff_t>(npts_0) * npts_1;
    } else {
      npts_model = npts_0;
    }

    // Now read the time (in days) at which the model is specified.
    double t_model_days{NAN};
    if (!get_noncommentline(fmodel.file, line)) {
      fatal_crash("model.txt ends after the cell count, before the line of the snapshot time");
    }
    std::istringstream{line} >> t_model_days;
    // a failed extraction stores zero, which would zero all densities via the (t_model / tmin)^3 scaling
    if (!std::isfinite(t_model_days) || t_model_days <= 0.) {
      fatal_crash("model.txt: could not read a positive snapshot time in days from line '{}'", line);
    }
    t_model = t_model_days * DAY;
    if (t_model > globals::tmin) {
      fatal_crash("The model snapshot time {} d is after tmin {} d", t_model_days, globals::tmin / DAY);
    }

    // if the next line is a single float, it is the vmax (so 2D or 3D)
    // otherwise, it is the first line of the model or a header comment (so 1D)
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends after the snapshot time, before the vmax line or the first cell");
    }
    if (!line.starts_with('#')) {
      double num_after_vmax{NAN};
      auto sslinevmax = std::istringstream{line};
      if ((sslinevmax >> globals::vmax) && !(sslinevmax >> num_after_vmax)) {
        // single value on the line is a vmax, so 2D or 3D
        // if it's not already know to be 2D (based on n_r n_z line), then it's 3D
        if (detected_model_type != GridType::CYLINDRICAL2D) {
          detected_model_type = GridType::CARTESIAN3D;
          printlnlog("Detected 3D model");
        }
      }
    }
    if (!detected_model_type.has_value()) {
      detected_model_type = GridType::SPHERICAL1D;
      printlnlog("Detected 1D model");
      fmodel.lines_read_ahead.push_back(line);
    }
  }

  model_type = detected_model_type;

  assert_always(modelgrid_input.empty());
  modelgrid_input = MPI_shared_array<ModelGridCellInput>(npts_model, ModelGridCellInput{});
  modelgrid_numpropcells.resize(npts_model, 0);
  nonemptymgi_of_mgi.resize(npts_model, -1);

  // geometry-specific setup of the model coordinate counts
  if (model_format_version == 2) {
    // read_model_keyword_header() set the cell counts from the cell edges
    printlnlog("model.txt gives {} cells on the first axis{}{}", ncoord_model[0],
               (ncoord_model[1] > 0) ? std::format(", {} on the second axis", ncoord_model[1]) : std::string{},
               (ncoord_model[2] > 0) ? std::format(", and {} on the third axis", ncoord_model[2]) : std::string{});
    if (get_modelgridtype() == GridType::CARTESIAN3D) {
      // for a 3D input model, the propagation cells will match the input cells exactly
      ncoordgrid = ncoord_model;
      ngrid = npts_model;
    }
  } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = 0;
    ncoord_model[2] = 0;
    vout_model.resize(get_npts_model(), NAN);
  } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = npts_1;
    ncoord_model[2] = 0;
  } else {
    ncoord_model[0] = static_cast<int>(round(std::cbrt(npts_model)));
    ncoord_model[1] = ncoord_model[0];
    ncoord_model[2] = ncoord_model[0];
    if (static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * ncoord_model[2] != npts_model) {
      fatal_crash("model.txt: a 3D model needs a cube number of cells, but it has {} cells", npts_model);
    }
    // for a 3D input model, the propagation cells will match the input cells exactly
    ncoordgrid = ncoord_model;
    ngrid = npts_model;
  }

  const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

  // only global rank 0 parses the cell data. The node-shared arrays it fills are broadcast to the leaders of the
  // other nodes below, and the rank-local results are broadcast to all ranks, so that only one rank reads the
  // (potentially very large) file
  if (globals::my_rank == 0) {
    // 3D only: whether the input cell positions match the expected values in x-y-z or z-y-x
    // column order (set false when a mismatch is detected)
    bool posmatch_xyz = true;
    bool posmatch_zyx = true;

    if (model_format_version == 2) {
      first_input_cellid = 0;
    }
    ptrdiff_t cell_lines_read = 0;
    int mgi = 0;
    while (mgi < get_npts_model() && fmodel.getline(line)) {
      auto remainder = std::string_view{line};
      int cellnumberin = 0;
      double rho_tmodel{NAN};  // the cell density [g/cm3] at the model snapshot time t_model

      // parse the geometry-specific part of the cell line to get the density (and for 1D, the
      // outer velocity), and check that any cell midpoint coordinates match the expected values
      if (model_format_version == 2) {
        if (lineiscommentonly(line)) {
          continue;
        }
        // columns: cell number, density. The header gives the cell edges. The cell number is the index of the cell
        // in the grid of the edges, from zero, with the index of the first axis as the fastest index.
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, rho_tmodel))) {
          fatal_crash("model.txt: a cell line must start with 2 numbers (inputcellid rho), but the line is: {}", line);
        }
        if (cellnumberin < mgi || cellnumberin >= get_npts_model()) {
          fatal_crash(
              "model.txt has the cell id {} after {} cells. The cell ids of format version 2 must increase, and the "
              "cell edges of the header give the ids 0 to {}",
              cellnumberin, mgi, get_npts_model() - 1);
        }
        // a cell with no line is empty
        for (; mgi < cellnumberin; mgi++) {
          set_rho_tmin(mgi, 0.F);
        }
        if (get_modelgridtype() == GridType::CARTESIAN3D) {
          // the reader of a version 1 file of a 3D model takes the density as a float, thus the same file in
          // version 2 gives the same density
          const auto rho_model_in = static_cast<float>(rho_tmodel);
          rho_tmodel = rho_model_in;
          if (min_den < 0. || min_den > rho_tmodel) {
            min_den = rho_tmodel;
          }
        }
      } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
        // columns: cell number, outer velocity [km/s], log10 of the density
        double vout_kmps{NAN};
        double log_rho{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, vout_kmps) &&
              parse_next_token(remainder, log_rho))) {
          fatal_crash(
              "model.txt cell {}: expected at least 3 values (inputcellid vel_r_max_kmps log10rho) but could "
              "not parse line: {}",
              mgi, line);
        }
        vout_model[mgi] = vout_kmps * 1.e5;
        // the velocity grid is binary searched by int_index_lowerbound(), so it has to increase
        // outwards. Without this check a non-monotonic model.txt gives silently wrong cell lookups.
        if (mgi > 0 && !(vout_model[mgi] > vout_model[mgi - 1])) {
          fatal_crash("model.txt cell {} has vel_r_max_kmps {:g}, but it must be larger than {:g} of the cell before",
                      mgi, vout_kmps, vout_model[mgi - 1] / 1.e5);
        }
        if (mgi == 0 && !(vout_model[0] > 0.)) {
          fatal_crash("model.txt cell 0 has vel_r_max_kmps {:g}, but it must be above 0", vout_kmps);
        }
        rho_tmodel = (log_rho > -90) ? pow(10., log_rho) : 0.;
      } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
        // columns: cell number, r midpoint, z midpoint, density
        float cell_r_in{NAN};
        float cell_z_in{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cell_r_in) &&
              parse_next_token(remainder, cell_z_in) && parse_next_token(remainder, rho_tmodel))) {
          fatal_crash(
              "model.txt cell {} must start with 4 numbers (inputcellid pos_rcyl_mid pos_z_mid rho), but the line is: "
              "{}",
              mgi, line);
        }

        const int n_rcyl = (mgi % ncoord_model[0]);
        const double pos_r_cyl_mid = (n_rcyl + 0.5) * globals::vmax * t_model / ncoord_model[0];
        if (!(fabs((cell_r_in / pos_r_cyl_mid) - 1) < 1e-3)) {
          fatal_crash("model.txt cell {} has pos_rcyl_mid {:g}, but the expected value is {:g}", mgi, cell_r_in,
                      pos_r_cyl_mid);
        }
        const int n_z = (mgi / ncoord_model[0]);
        if (((2 * n_z) + 1) == ncoord_model[1]) {
          // an odd number of z rows puts the middle row at the origin, where a relative comparison
          // has nothing to divide by. Compare against the cell size instead. The row is identified
          // from the integer index rather than by testing the coordinate against zero, because the
          // expression below is not required to evaluate to exactly zero under -ffast-math.
          if (!(fabs(cell_z_in) < (1e-3 * globals::vmax * t_model / ncoord_model[1]))) {
            fatal_crash("model.txt cell {} has pos_z_mid {:g}, but the expected value is 0", mgi, cell_z_in);
          }
        } else {
          const double pos_z_mid = globals::vmax * t_model * (-1 + (2 * (n_z + 0.5) / ncoord_model[1]));
          if (!(fabs((cell_z_in / pos_z_mid) - 1) < 1e-3)) {
            fatal_crash("model.txt cell {} has pos_z_mid {:g}, but the expected value is {:g}", mgi, cell_z_in,
                        pos_z_mid);
          }
        }
      } else {
        // columns: cell number, x minimum, y minimum, z minimum, density. The position columns hold the corner
        // of the cell with the lowest coordinates. The check below permits a difference of half a cell width, so
        // it cannot tell a file of cell midpoints from a file of cell corners.
        std::array<float, 3> cellpos_in{};
        float rho_model_in{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cellpos_in[0]) &&
              parse_next_token(remainder, cellpos_in[1]) && parse_next_token(remainder, cellpos_in[2]) &&
              parse_next_token(remainder, rho_model_in))) {
          fatal_crash(
              "model.txt cell {} must start with 5 numbers (inputcellid pos_x_min pos_y_min pos_z_min rho), but the "
              "line is: {}",
              mgi, line);
        }
        rho_tmodel = rho_model_in;

        // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
        // however, the cellindex always should increment X first, then Y, then Z
        const double xmax_tmodel = globals::vmax * t_model;
        for (int axis = 0; axis < 3; axis++) {
          const double cellwidth = 2 * xmax_tmodel / ncoordgrid[axis];
          const double cellpos_expected = -xmax_tmodel + (cellwidth * get_cellcoordindex(mgi, axis));
          if (fabs(cellpos_expected - cellpos_in[axis]) > 0.5 * cellwidth) {
            posmatch_xyz = false;
          }
          if (fabs(cellpos_expected - cellpos_in[2 - axis]) > 0.5 * cellwidth) {
            posmatch_zyx = false;
          }
        }

        if (min_den < 0. || min_den > rho_tmodel) {
          min_den = rho_tmodel;
        }
      }

      if (mgi == 0 && model_format_version == 1) {
        first_input_cellid = cellnumberin;
        printlnlog("model.txt: autodetected cell numbering starting from {}", first_input_cellid);
        if (first_input_cellid != 0 && first_input_cellid != 1) {
          fatal_crash("model.txt: the first cell id is {}, but it must be 0 or 1", first_input_cellid);
        }
      }
      if (cellnumberin != mgi + first_input_cellid) {
        fatal_crash("model.txt cell {} has the cell id {}, but the expected id is {}", mgi, cellnumberin,
                    mgi + first_input_cellid);
      }

      if (rho_tmodel < 0) {
        fatal_crash("model.txt cell {} has negative density {:g} [g/cm3] at t_model", mgi, rho_tmodel);
      }

      // The model holds the density at tmin as a float. A larger value becomes infinity in that type, and the mass of
      // the cell is then infinite. The density of a 3D model at t_model is already a float here.
      const double rho_tmin = rho_tmodel * pow3(t_model / globals::tmin);
      if (!(rho_tmin <= std::numeric_limits<float>::max())) {
        fatal_crash(
            "model.txt cell {} has a density that a float cannot hold: {:g} [g/cm3] at t_model, and {:g} [g/cm3] at "
            "tmin",
            mgi, rho_tmodel, rho_tmin);
      }

      const bool keepcell = (rho_tmodel > 0);
      set_rho_tmin(mgi, static_cast<float>(rho_tmin));

      read_model_radioabundances(fmodel, remainder, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      cell_lines_read++;
      mgi++;
    }

    if (model_format_version == 2) {
      // a cell with no line is empty, also after the last line of the file
      for (; mgi < get_npts_model(); mgi++) {
        set_rho_tmin(mgi, 0.F);
      }
    }
    if (mgi != get_npts_model()) {
      fatal_crash("model.txt: found only {} cells instead of {} expected.", mgi, get_npts_model());
    }

    if (model_format_version == 2) {
      // the header of this format has no cell count, thus a cell after the last one is a sign of a wrong edge list
      while (fmodel.getline(line)) {
        if (!lineiscommentonly(line)) {
          fatal_crash("model.txt has more than {} cells, but the cell edges of its header give {} cells",
                      get_npts_model(), get_npts_model());
        }
      }
      if (model_cell_line_count.has_value() && cell_lines_read != model_cell_line_count.value()) {
        fatal_crash(
            "model.txt: the key cell_line_count gives {} cell lines, but the file has {}. The file is possibly not "
            "complete",
            model_cell_line_count.value(), cell_lines_read);
      }
    }

    if (get_modelgridtype() == GridType::SPHERICAL1D) {
      // for 1D models, vmax is the outer velocity of the last cell
      globals::vmax = vout_model[get_npts_model() - 1];
    } else if (get_modelgridtype() == GridType::CARTESIAN3D && model_format_version == 2) {
      printlnlog("minimum model density {:g} [g/cm3]", min_den);
    } else if (get_modelgridtype() == GridType::CARTESIAN3D) {
      // posmatch_xyz and posmatch_zyx should be mutually exclusive: if both column orders matched, an
      // infinity probably occurred in the calculated positions
      if (posmatch_xyz) {
        printlnlog(
            "Cell positions in model.txt are consistent with calculated values when x-y-z column order is used.");
      } else if (posmatch_zyx) {
        printlnlog(
            "Cell positions in model.txt are consistent with calculated values when z-y-x column order is used.");
      } else {
        printlnlog(
            "[warning] Cell positions in model.txt are not consistent with calculated values in either x-y-z or z-y-x "
            "order.");
      }
      printlnlog("minimum model density {:g} [g/cm3]", min_den);
    }
  }

  // share the rank-local results of the cell parsing with the other ranks
  MPI_Bcast_safe(first_input_cellid, 0, MPI_COMM_WORLD);
  MPI_Bcast_safe(globals::vmax, 0, MPI_COMM_WORLD);
  if (get_modelgridtype() == GridType::SPHERICAL1D) {
    MPI_Bcast_safe(vout_model, 0, MPI_COMM_WORLD);
  }

  // send the node-shared model data filled by rank 0 to the node leaders of the other nodes
  if (globals::rank_in_node == 0) {
    MPI_Bcast_safe(modelgrid_input, 0, globals::mpi_comm_internode);
    MPI_Bcast_safe(initnucmassfrac_allcells, 0, globals::mpi_comm_internode);
    if (!model_otherstable_massfracs.empty()) {
      MPI_Bcast_safe(model_otherstable_massfracs, 0, globals::mpi_comm_internode);
    }
  }
  MPI_Barrier_allranks();

  assert_always(get_npts_model() ==
                std::max(1, ncoord_model[0]) * std::max(1, ncoord_model[1]) * std::max(1, ncoord_model[2]));
  printlnlog("npts_model: {}", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printlnlog("vmax {:g} [cm/s] ({:.2f}c)", globals::vmax, globals::vmax / CLIGHT);
  if (!(globals::vmax < CLIGHT)) {
    fatal_crash("model.txt gives vmax = {:g} cm/s, but vmax must be below the speed of light", globals::vmax);
  }
  printlnlog("tmin {:g} [s] = {:.2f} [d]", globals::tmin, globals::tmin / DAY);
  printlnlog("rmax {:g} [cm] (at t=tmin)", globals::rmax);
  if (get_modelgridtype() == GridType::SPHERICAL1D) {
    std::vector<double> vel_r_edges{0.};
    vel_r_edges.insert(vel_r_edges.end(), vout_model.begin(), vout_model.end());
    require_increasing_shell_radii_at_tmin(vel_r_edges);
  } else if (has_spherical_angle_axes(get_modelgridtype())) {
    require_increasing_shell_radii_at_tmin(model_vel_r_edges);
  }

  calc_modelinit_totmassnuclides();

  printlnlog("Total input model mass: {:9.3e} [Msun]", mtot_input / MSUN);
  printlnlog("Nuclide masses at t=t_model_init [Msun]:");
  printlnlog("  56Ni: {:9.3e}  56Co: {:9.3e}  52Fe: {:9.3e}  48Cr: {:9.3e}", get_totmassnuclide_tmodel(28, 56) / MSUN,
             get_totmassnuclide_tmodel(27, 56) / MSUN, get_totmassnuclide_tmodel(26, 52) / MSUN,
             get_totmassnuclide_tmodel(24, 48) / MSUN);
  printlnlog("  Fe-group: {:9.3e}  57Ni: {:9.3e}  57Co: {:9.3e}", mfegroup / MSUN,
             get_totmassnuclide_tmodel(28, 57) / MSUN, get_totmassnuclide_tmodel(27, 57) / MSUN);

  report_negative_input_massfracs();

  if (!ignored_model_columns.empty()) {
    printlog(
        "[warning] ignored model.txt columns not recognised as a known nuclide, X_Fegroup, Ye, q, or tracercount:");
    for (const auto& colname : ignored_model_columns) {
      printlog(" '{}'", colname);
    }
    printlnlog("");
  }

  read_possible_yefile();
}

}  // namespace grid
