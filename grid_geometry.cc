// The geometry of the propagation grid: the setup of each grid type, the coordinates and the volumes of the
// propagation cells, and the distance of a packet to the boundary of its cell.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <iterator>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <string_view>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "grid_internal.h"
#include "mpi_logging.h"
#include "random.h"
#include "sn3d.h"
#include "vectors.h"

namespace grid {

namespace {

// Initial co-ordinates of inner most corner of cell. The polar axis of a spherical grid holds the lower cos(theta) of
// each cell, and the azimuthal axis holds the lower angle [rad]. These values do not change with time.
std::array<std::vector<double>, 3> coord_pos_min_tmin{};

// The upper bound of the last cell on each axis. A length axis ends at globals::rmax. The polar axis of a spherical
// grid ends at cos(theta) = 1, and the azimuthal angle ends at 2 pi.
std::array<double, 3> coord_pos_max_of_grid_tmin{};

// The cosine and the sine of each cell edge of the polar angle and of the azimuthal angle of a spherical propagation
// grid. The polar edges have the sequence of the cells, thus their cosine increases. boundary_distance() runs for
// every packet step, thus setup_grid_spherical() calculates them one time.
std::vector<double> theta_edge_cos{};
std::vector<double> theta_edge_sin{};
std::vector<double> phi_edge_cos{};
std::vector<double> phi_edge_sin{};

// Get the number of axes that are lengths. The cell boundaries on such an axis move with the homologous expansion.
// These axes come first. Each later axis of a spherical grid is an angle, and its cell boundaries do not move.
constexpr auto get_num_length_axes(const GridType gridtype) -> int {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return 1;
    case GridType::CYLINDRICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
      return 3;
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

// maximum out-of-cell position error (in cm) attributable to floating-point rounding of the
// packet position updates. The tolerance scales with the coordinate magnitude because
// rounding errors are proportional to ulp(pos), e.g. about 8 cm at 4e16 cm.
[[nodiscard]] constexpr auto cellbound_tolerance(const double boundarypos) -> double {
  return std::max(10., std::abs(boundarypos) * 1e-12);
}

// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r direction)
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_coordcellindexstride(const int axis) -> int {
  int stride = 1;
  for (int a = 0; a < axis; ++a) {
    stride *= ncoordgrid[a];
  }
  return stride;
}

// Get the index of the cell that holds the value, from the lower edge of each cell of one axis. A cell holds its
// lower edge. A value below the first edge is in the first cell, and the last cell holds each value above its lower
// edge.
[[nodiscard]] DEVICE_FUNC constexpr auto get_cellindex_from_lower_edges(const std::span<const double> lower_edges,
                                                                        const double value) -> int {
  const auto num_cells = static_cast<int>(std::ssize(lower_edges));
  const auto index_above = static_cast<int>(std::ranges::upper_bound(lower_edges, value) - lower_edges.begin());
  return std::clamp(index_above - 1, 0, num_cells - 1);
}

// Get cos(theta) of the polar angle theta from the +z axis and the azimuthal angle phi = atan2(y, x) in [0, 2 pi) of
// a position or a direction. A vector on the z axis has phi = 0. A vector of zero length has cos(theta) = -1, which
// is in the first polar cell.
[[nodiscard]] DEVICE_FUNC inline auto get_costheta_phi_from_xyz(const Vec3d& vec_xyz) -> std::array<double, 2> {
  const double length = vec_len(vec_xyz);
  const double costheta = (length > 0.) ? std::clamp(vec_xyz[2] / length, -1., 1.) : -1.;
  double phi = (vec_xyz[0] == 0. && vec_xyz[1] == 0.) ? 0. : std::atan2(vec_xyz[1], vec_xyz[0]);
  if (phi < 0.) {
    phi += 2 * PI;
  }
  // the sum of 2 pi and a small negative angle can round to 2 pi, which is the angle 0
  if (phi >= 2 * PI) {
    phi = 0.;
  }
  return {costheta, phi};
}

// Stop the run when the grid corners expand faster than light and nothing removes them. The
// propagation cannot treat a superluminal boundary, so the setup must either cut the corners with
// FORCE_SPHERICAL_ESCAPE_SURFACE or use a smaller vmax.
void require_subluminal_corners(const double vmax_corner) {
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  if (!FORCE_SPHERICAL_ESCAPE_SURFACE && vmax_corner >= CLIGHT) {
    fatal_crash(
        "the grid corners expand faster than light. Enable FORCE_SPHERICAL_ESCAPE_SURFACE in "
        "artisoptions.h to remove the superluminal corners, or reduce vmax.");
  }
}

// Get the discrete index of the coordinate value (where pos must be position in grid coordinate system, not
// necessarily xyz)
auto get_poscoordpointnum(const double pos, const double time, const int axis) -> int {
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D: {
      const auto idx = static_cast<int>(((pos / time) + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
      assert_always(idx >= 0);
      assert_always(idx < ncoordgrid[axis]);
      return idx;
    }

    case GridType::CYLINDRICAL2D: {
      if (axis == 0) {
        const auto n_rcyl = static_cast<int>(pos / time / globals::vmax * ncoordgrid[axis]);
        assert_always(n_rcyl >= 0);
        assert_always(n_rcyl < ncoordgrid[axis]);
        return n_rcyl;
      }
      if (axis == 1) {
        const auto n_z = static_cast<int>(((pos / time) + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
        assert_always(n_z >= 0);
        assert_always(n_z < ncoordgrid[axis]);
        return n_z;
      }
      break;
    }

    case GridType::SPHERICAL1D: {
      // radial spacing is non-uniform, so we have to do a search
      const auto trat = time / globals::tmin;
      for (int n_r = 0; n_r < ncoordgrid[0]; n_r++) {
        if ((pos < get_cellcoordmax(n_r, 0) * trat) && (pos >= get_cellcoordmin(n_r, 0) * trat)) {
          return n_r;
        }
      }
      break;
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      // get_cellindex_from_pos() searches the cell edges of these grids itself
      break;
  }

  fatal_crash("Position {:g} on axis {} at time {:g} has no cell index", pos, axis, time);
}

// Convert a position vector from Cartesian xyz to the grid coordinate system
[[nodiscard]] constexpr auto get_gridcoords_from_xyz(const Vec3d& pos_xyz, const GridType gridtype) -> Vec3d {
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return pos_xyz;
    case GridType::CYLINDRICAL2D:
      return {std::sqrt(pow2(pos_xyz[0]) + pow2(pos_xyz[1])), pos_xyz[2], NAN};
    // only the radius of a spherical grid is a length. boundary_distance() treats the angles of a 2D spherical and a
    // 3D spherical grid separately
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return {vec_len(pos_xyz), NAN, NAN};
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

// get the velocity in the grid coordinate system from the xyz position and direction
[[nodiscard]] constexpr auto get_gridcoords_vel_from_xyz_pos_dir(const Vec3d& pos_xyz, const Vec3d& dir_xyz,
                                                                 const std::span<const double> pktposgridcoord,
                                                                 const GridType gridtype) -> Vec3d {
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return {dir_xyz[0] * CLIGHT_PROP, dir_xyz[1] * CLIGHT_PROP, dir_xyz[2] * CLIGHT_PROP};

    case GridType::CYLINDRICAL2D: {
      const double v_rcyl = pktposgridcoord[0] > 0. ? ((pos_xyz[0] * dir_xyz[0]) + (pos_xyz[1] * dir_xyz[1])) /
                                                          pktposgridcoord[0] * CLIGHT_PROP
                                                    : std::hypot(dir_xyz[0], dir_xyz[1]) * CLIGHT_PROP;
      const double v_z = dir_xyz[2] * CLIGHT_PROP;
      return {v_rcyl, v_z, NAN};
    }

    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      const double v_radial =
          (pktposgridcoord[0] > 0. ? dot(pos_xyz, dir_xyz) / pktposgridcoord[0] : vec_len(dir_xyz)) * CLIGHT_PROP;
      return {v_radial, NAN, NAN};
    }
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

[[nodiscard]] constexpr auto distance_cartesian_boundary(const double pktposgridcoord, const double pktvelgridcoord,
                                                         const double cellboundarypos, const double tstart) -> double {
  // numerically stable formulation: compute time difference directly as
  // (pos - boundary_at_tstart) / (boundary_velocity - packet_velocity)
  // to avoid catastrophic cancellation when t_crossing ≈ tstart
  return CLIGHT_PROP * (pktposgridcoord - (cellboundarypos / globals::tmin * tstart)) /
         ((cellboundarypos / globals::tmin) - pktvelgridcoord);
}

// Check if the packet is at or (within rounding-error tolerance) past a cell boundary that it
// is moving towards relative to the boundary's homologous motion. In that case the
// intersection calculation finds no forward crossing (a slightly negative distance, or no
// intersection at all for the curved boundaries), and the crossing must instead be treated as
// immediate (zero distance). If it were silently dropped, the crossing could never be
// detected again (the packet is outrunning the boundary), and the out-of-cell position error
// would grow without bound along the rest of the path until the outside-boundary check fails.
template <BoundaryType boundarytype>
[[nodiscard]] constexpr auto is_boundary_overshoot_within_tolerance(const double pktposgridcoord,
                                                                    const double pktvelgridcoord,
                                                                    const double boundarypos_tmin, const double tstart)
    -> bool {
  const double boundaryvel = boundarypos_tmin / globals::tmin;
  const double boundarypos = boundaryvel * tstart;
  const double overshoot =
      boundarytype == BoundaryType::UPPER ? (pktposgridcoord - boundarypos) : (boundarypos - pktposgridcoord);
  const bool movingtowards =
      boundarytype == BoundaryType::UPPER ? (pktvelgridcoord > boundaryvel) : (pktvelgridcoord < boundaryvel);
  return movingtowards && (overshoot >= 0.) && (overshoot <= cellbound_tolerance(boundarypos));
}

// The distance [cm] to an angular face or to a point before a close approach to the centre or the z axis.
// The cell indices apply after the step. A negative distance means that the ray needs no angular step.
struct AngularBoundaryCrossing {
  double distance;
  int index_costheta;
  int index_phi;
};

// Find the first crossing of a polar or an azimuthal cell boundary on the path from pos along the unit vector dir,
// for a packet in the cell with the given polar and azimuthal index. The cones of constant cos(theta) and the
// half-planes of constant phi are stationary. The time sets the shortest step that can advance the packet clock.
[[nodiscard]] DEVICE_FUNC auto get_first_angular_boundary_crossing(const Vec3d& pos, const Vec3d& dir,
                                                                   const int index_costheta, const int index_phi,
                                                                   const double time) -> AngularBoundaryCrossing {
  const int num_costheta = ncoordgrid[1];
  const int num_phi = ncoordgrid[2];
  AngularBoundaryCrossing first_crossing{.distance = -1., .index_costheta = index_costheta, .index_phi = index_phi};
  const auto take_crossing_if_nearer = [&first_crossing](const double distance, const int index_costheta_after,
                                                         const int index_phi_after) {
    if (distance >= 0. && (first_crossing.distance < 0. || distance < first_crossing.distance)) {
      first_crossing = {.distance = distance, .index_costheta = index_costheta_after, .index_phi = index_phi_after};
    }
  };

  // The angular momentum pos x dir of a unit direction has the length of the smallest distance of the path from the
  // centre. Its z component divided by the length of dir in the xy plane is the smallest distance from the z axis.
  if (num_costheta == 1 && num_phi == 1) {
    return first_crossing;
  }
  const auto angmom = ray_angular_momentum(pos, dir);
  const double dir_xy_squared = pow2(dir[0]) + pow2(dir[1]);
  if (angmom[0] == 0. && angmom[1] == 0. && angmom[2] == 0.) {
    // The path goes through the centre, e.g. the radial path of a new packet. Its cos(theta) and phi are constant on
    // each side of the centre, and after the centre they are the values of dir.
    const auto [costheta_of_dir, phi_of_dir] = get_costheta_phi_from_xyz(dir);
    const int index_costheta_after = get_cellindex_from_lower_edges(coord_pos_min_tmin[1], costheta_of_dir);
    const int index_phi_after = get_cellindex_from_lower_edges(coord_pos_min_tmin[2], phi_of_dir);
    if (index_costheta_after != index_costheta || index_phi_after != index_phi) {
      // At the centre itself, the direction selects the cell for the next positive path length.
      take_crossing_if_nearer(std::max(-dot(pos, dir) / dot(dir, dir), 0.), index_costheta_after, index_phi_after);
    }
    return first_crossing;
  }

  // A long step can lose the small transverse offset when the ray approaches the centre or the z axis.
  // Stop before the angular faces when the distance to the axis is below the square root of machine epsilon
  // times the distance along the ray. The next call then resolves the faces from a nearby Cartesian position.
  const auto approach_if_unresolved = [&](const double distance_to_closest, const double impact_parameter) {
    const double relative_distance = std::sqrt(std::numeric_limits<double>::epsilon());
    if (distance_to_closest > 0. && impact_parameter < relative_distance * distance_to_closest) {
      // A shorter step must advance the packet time. At smaller lengths, stop at closest approach.
      const double distance = distance_to_closest > std::numeric_limits<double>::epsilon() * CLIGHT_PROP * time
                                  ? distance_to_closest * (1. - relative_distance)
                                  : distance_to_closest;
      take_crossing_if_nearer(distance, index_costheta, index_phi);
    }
  };
  const double dir_squared = dot(dir, dir);
  approach_if_unresolved(-dot(pos, dir) / dir_squared, vec_len(angmom) / dir_squared);

  // the tolerance of a boundary that the packet is at, or just went through (see
  // is_boundary_overshoot_within_tolerance). An angle tolerance has no use near the centre, thus the test uses
  // the distance [cm] from the boundary
  const double tolerance = cellbound_tolerance(vec_len(pos));

  // The sign of r^3 d(cos(theta))/ds gives the direction of a polar crossing.
  const double rcyl = std::sqrt(pow2(pos[0]) + pow2(pos[1]));
  const double polar_rate = difference_of_products(pos[1], angmom[0], pos[0], angmom[1]);

  // The polar axis is not a cell boundary, thus the first cell has no lower boundary and the last cell has no upper
  // boundary. cos(theta) can have an extremum on the path, thus both boundaries get a test for each direction.
  const auto take_costheta_boundary_if_nearer = [&]<BoundaryType boundarytype>(const int index_costheta_after) {
    const bool is_upper = (boundarytype == BoundaryType::UPPER);
    const int edgeindex = is_upper ? index_costheta + 1 : index_costheta;
    // r * sin(theta_edge - theta): the distance of the packet from the cone, positive on the side of larger
    // cos(theta)
    const double distance_from_cone = (theta_edge_sin[edgeindex] * pos[2]) - (theta_edge_cos[edgeindex] * rcyl);
    // At an extremum of cos(theta), a non-radial ray moves towards the equator on its next positive step.
    const double outgoing_polar_sense = polar_rate != 0. ? polar_rate : -pos[2];
    const bool moves_through_boundary = is_upper ? outgoing_polar_sense > 0. : outgoing_polar_sense < 0.;
    bool crossed_face = false;
    if (moves_through_boundary && std::abs(distance_from_cone) <= tolerance) {
      double signed_side = distance_from_cone;
      const double mu = theta_edge_cos[edgeindex];
      if (mu == 0.) {
        signed_side = pos[2];
      } else if (pos[2] * mu > 0.) {
        const double cone_product = polar_cone_product(pos, pos, mu);
        signed_side = mu > 0. ? cone_product : -cone_product;
      }
      crossed_face = is_upper ? signed_side >= 0. : signed_side <= 0.;
    }
    const double distance =
        crossed_face ? 0. : costheta_cone_intersection<boundarytype>(pos, dir, theta_edge_cos[edgeindex]);
    take_crossing_if_nearer(distance, index_costheta_after, index_phi);
  };
  if (index_costheta < num_costheta - 1) {
    take_costheta_boundary_if_nearer.template operator()<BoundaryType::UPPER>(index_costheta + 1);
  }
  if (index_costheta > 0) {
    take_costheta_boundary_if_nearer.template operator()<BoundaryType::LOWER>(index_costheta - 1);
  }

  // A direction along the z axis keeps phi constant. A 2D spherical grid has one cell of phi.
  if (num_phi == 1 || dir_xy_squared <= 0.) {
    return first_crossing;
  }

  if (angmom[2] == 0.) {
    // The path goes through the z axis, where phi changes by pi. After the axis, phi is the angle of dir.
    const int index_phi_after =
        get_cellindex_from_lower_edges(coord_pos_min_tmin[2], get_costheta_phi_from_xyz(dir)[1]);
    if (index_phi_after != index_phi) {
      const double distance_to_axis = -((pos[0] * dir[0]) + (pos[1] * dir[1])) / dir_xy_squared;
      take_crossing_if_nearer(std::max(distance_to_axis, 0.), index_costheta, index_phi_after);
    }
    return first_crossing;
  }

  approach_if_unresolved(-((pos[0] * dir[0]) + (pos[1] * dir[1])) / dir_xy_squared,
                         std::abs(angmom[2]) / dir_xy_squared);

  // phi is periodic: the half-plane phi = 2 pi of the last cell is the half-plane phi = 0 of the first cell
  const auto take_phi_boundary_if_nearer = [&]<BoundaryType boundarytype>(const int index_phi_after) {
    const bool is_upper = (boundarytype == BoundaryType::UPPER);
    const int edgeindex = is_upper ? index_phi + 1 : index_phi;
    const double cos_edge = phi_edge_cos[edgeindex];
    const double sin_edge = phi_edge_sin[edgeindex];
    // the distance of the packet from the full plane of the boundary, positive on the side of larger phi
    const double distance_from_plane = difference_of_products(cos_edge, pos[1], sin_edge, pos[0]);
    const double distance_from_plane_rate = difference_of_products(cos_edge, dir[1], sin_edge, dir[0]);
    const bool moves_through_boundary = is_upper ? (angmom[2] > 0. && distance_from_plane_rate > 0.)
                                                 : (angmom[2] < 0. && distance_from_plane_rate < 0.);
    // the other half of the plane has phi = phi_edge + pi
    const bool is_near_boundary_half = ((cos_edge * pos[0]) + (sin_edge * pos[1])) > 0.;
    const double overshoot = is_upper ? distance_from_plane : -distance_from_plane;
    const double distance =
        (moves_through_boundary && is_near_boundary_half && overshoot >= 0. && overshoot <= tolerance)
            ? 0.
            : phi_halfplane_intersection<boundarytype>(pos, dir, cos_edge, sin_edge);
    take_crossing_if_nearer(distance, index_costheta, index_phi_after);
  };
  take_phi_boundary_if_nearer.template operator()<BoundaryType::UPPER>((index_phi + 1) % num_phi);
  take_phi_boundary_if_nearer.template operator()<BoundaryType::LOWER>((index_phi + num_phi - 1) % num_phi);

  return first_crossing;
}

}  // anonymous namespace

[[gnu::pure]] DEVICE_FUNC auto get_propgridtype() -> GridType {
  if constexpr (GRID_TYPE_OVERRIDE.has_value()) {
    return GRID_TYPE_OVERRIDE.value();
  }
  return get_modelgridtype();
}

// convert a cell index into the integer coordinate index (x, y, z, or r) in [0, ncoordgrid[axis])
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordindex(const int cellindex, const int axis) -> int {
  return (cellindex / get_coordcellindexstride(axis)) % ncoordgrid[axis];
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(const int cellindex, const int axis) -> double {
  return coord_pos_min_tmin[axis][get_cellcoordindex(cellindex, axis)];
}

// get the maximum position value of a coordinate axis at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the maximum x position in xyz coords, or the maximum radius in spherical 1D. An angular axis of a spherical
// grid gives the upper angle [rad] of the cell.
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmax(const int cellindex, const int axis) -> double {
  const auto idx = get_cellcoordindex(cellindex, axis);
  const auto idxlast = ncoordgrid[axis] - 1;
  return idx < idxlast ? coord_pos_min_tmin[axis][idx + 1] : coord_pos_max_of_grid_tmin[axis];
}

// return the inner radius of a propagation cell at time tmin
auto get_cell_r_inner(const int cellindex, const GridType prop_gridtype) -> double {
  if (prop_gridtype == GridType::SPHERICAL1D || has_spherical_angle_axes(prop_gridtype)) {
    return get_cellcoordmin(cellindex, 0);
  }

  // minimum |coordinate| within the cell along one axis: zero if the cell straddles the coordinate plane
  const auto axis_mindist = [](const double coordmin, const double coordmax) {
    return (coordmin <= 0. && coordmax >= 0.) ? 0. : std::min(std::abs(coordmin), std::abs(coordmax));
  };

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const auto rcyl_inner = get_cellcoordmin(cellindex, 0);
    const auto z_inner = axis_mindist(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1));
    return std::sqrt(pow2(rcyl_inner) + pow2(z_inner));
  }

  if (prop_gridtype == GridType::CARTESIAN3D) {
    const auto x_inner = axis_mindist(get_cellcoordmin(cellindex, 0), get_cellcoordmax(cellindex, 0));
    const auto y_inner = axis_mindist(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1));
    const auto z_inner = axis_mindist(get_cellcoordmin(cellindex, 2), get_cellcoordmax(cellindex, 2));
    return std::sqrt(pow2(x_inner) + pow2(y_inner) + pow2(z_inner));
  }

  fatal_crash("Unknown propagation grid type {}", static_cast<int>(prop_gridtype));
}

// radius of the midpoint of a propagation cell at time tmin. The midpoint has the middle value of each grid
// coordinate: r for a spherical grid, (rcyl, z) for a cylindrical grid, and (x, y, z) for a Cartesian grid.
auto get_cellradialposmid(const int cellindex) -> double {
  double radius_squared = 0.;
  for (int axis = 0; axis < get_num_length_axes(get_propgridtype()); axis++) {
    radius_squared += pow2((get_cellcoordmin(cellindex, axis) + get_cellcoordmax(cellindex, axis)) / 2);
  }
  return std::sqrt(radius_squared);
}

// volume averaged mean of the squared radius of a propagation cell. The kinetic energy sum needs the
// mean of r^2, which is larger than the square of the radius that get_cellradialposmid() gives.
auto get_cellradialposmeansquared(const int cellindex) -> double {
  const auto prop_gridtype = get_propgridtype();
  if (prop_gridtype == GridType::SPHERICAL1D || has_spherical_angle_axes(prop_gridtype)) {
    // each cell of a spherical grid is a part of a shell with a constant solid angle, thus the angles cancel
    const double r_inner = get_cellcoordmin(cellindex, 0);
    const double r_outer = get_cellcoordmax(cellindex, 0);
    return 3. / 5. * (pow5(r_outer) - pow5(r_inner)) / (pow3(r_outer) - pow3(r_inner));
  }

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const double rcyl_inner = get_cellcoordmin(cellindex, 0);
    const double rcyl_outer = get_cellcoordmax(cellindex, 0);
    const double z_low = get_cellcoordmin(cellindex, 1);
    const double z_high = get_cellcoordmax(cellindex, 1);
    return ((pow2(rcyl_outer) + pow2(rcyl_inner)) / 2.) + ((pow2(z_low) + (z_low * z_high) + pow2(z_high)) / 3.);
  }

  // cubic cell: the mean of x^2 + y^2 + z^2 separates into a sum of the per-axis means
  double meansq = 0.;
  for (int axis = 0; axis < 3; axis++) {
    const double low = get_cellcoordmin(cellindex, axis);
    const double high = get_cellcoordmax(cellindex, axis);
    meansq += (pow2(low) + (low * high) + pow2(high)) / 3.;
  }
  return meansq;
}

// Get the index of the cell of a 2D spherical or a 3D spherical model that holds a position at tmin, or -1 for a
// position outside the outermost shell. The index is i_r + n_r * (i_costheta + n_costheta * i_phi).
[[nodiscard]] auto get_spherical_model_cellindex_tmin(const Vec3d& pos_tmin) -> int {
  const double vel_r = vec_len(pos_tmin) / globals::tmin;
  if (vel_r > model_vel_r_edges.back()) {
    return -1;
  }
  const auto [costheta_pos, phi] = get_costheta_phi_from_xyz(pos_tmin);
  // The polar axis of a mirrored model starts at the plane z = 0, and a position with z < 0 takes its mirror cell.
  // The centre is in the first polar cell.
  const double costheta = (model_mirror_z && vel_r > 0.) ? std::abs(costheta_pos) : costheta_pos;
  // the last value of each list is the upper edge of the last cell
  const auto lower_edges = [](const std::vector<double>& edges) {
    return std::span<const double>{edges}.first(edges.size() - 1);
  };
  const int index_r = get_cellindex_from_lower_edges(lower_edges(model_vel_r_edges), vel_r);
  const int index_costheta = get_cellindex_from_lower_edges(lower_edges(model_costheta_edges), costheta);
  const int index_phi = get_cellindex_from_lower_edges(lower_edges(model_phi_edges), phi);
  return index_r + (ncoord_model[0] * (index_costheta + (ncoord_model[1] * index_phi)));
}

void setup_grid_cartesian_3d() {
  // vmax is per coordinate, but the simulation volume corners will have a higher expansion velocity than the sides
  require_subluminal_corners(sqrt(3 * pow2(globals::vmax)));

  // Set grid size for uniform xyz grid
  // a 3D ejecta model gives the propagation grid its own cell counts
  if (get_modelgridtype() != GridType::CARTESIAN3D) {
    ncoordgrid = {CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, CUBOID_NCOORDGRID_Z};
  }

  // artis assumes in some places that the cells are cubes, not cuboids
  if (ncoordgrid[0] <= 0 || ncoordgrid[0] != ncoordgrid[1] || ncoordgrid[0] != ncoordgrid[2]) {
    fatal_crash(
        "the 3D propagation grid has {} x {} x {} cells, but the three counts must be equal and positive. The counts "
        "come from {}",
        ncoordgrid[0], ncoordgrid[1], ncoordgrid[2],
        get_modelgridtype() == GridType::CARTESIAN3D ? "model.txt"
                                                     : "CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, and "
                                                       "CUBOID_NCOORDGRID_Z in artisoptions.h");
  }

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];
  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  reserve_resize(coord_pos_min_tmin[2], ncoordgrid[2]);

  for (int axis = 0; axis < 3; axis++) {
    for (int i = 0; i < ncoordgrid[axis]; i++) {
      coord_pos_min_tmin[axis][i] = -globals::rmax + (2 * i * globals::rmax / ncoordgrid[axis]);
    }
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
}

// The readers of model.txt check the radial edges as velocities. Two edges that differ in the last bits can give
// radii at tmin that are equal or nearly equal. A packet at a boundary is in the position tolerance of that boundary
// (see cellbound_tolerance). In a shell that is not wider than this tolerance, the packet is at once at the next
// boundary, thus the shell gets no path length and no estimators. A model with such a shell is not correct on each
// propagation grid, thus each spherical model gets this check after the read, also on the 3D Cartesian grid.
void require_increasing_shell_radii_at_tmin(const std::span<const double> vel_r_edges) {
  for (ptrdiff_t n_r = 0; n_r < std::ssize(vel_r_edges) - 1; n_r++) {
    const double radius_inner_tmin = vel_r_edges[n_r] * globals::tmin;
    const double radius_outer_tmin = vel_r_edges[n_r + 1] * globals::tmin;
    const double tolerance = cellbound_tolerance(radius_outer_tmin);
    if (!(radius_outer_tmin - radius_inner_tmin > tolerance)) {
      fatal_crash(
          "model.txt: the radial shell {} of the model has the inner radius {:.17g} [cm] and the outer radius {:.17g} "
          "[cm] at tmin. Each shell must be wider than the position tolerance of a cell boundary, which is {:g} [cm] "
          "at this radius",
          n_r, radius_inner_tmin, radius_outer_tmin, tolerance);
    }
  }
}

void setup_grid_spherical_1d() {
  assert_always(get_modelgridtype() == GridType::SPHERICAL1D);

  ncoordgrid = {get_npts_model(), 1, 1};

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];

  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = cellindex;  // interchangeable in this mode
    const double v_inner = mgi > 0 ? vout_model[mgi - 1] : 0.;
    coord_pos_min_tmin[0][cellindex] = v_inner * globals::tmin;
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
}

void setup_grid_cylindrical_2d() {
  require_subluminal_corners(sqrt(2 * pow2(globals::vmax)));

  assert_always(get_modelgridtype() == GridType::CYLINDRICAL2D);

  ncoordgrid = ncoord_model;
  if (model_mirror_z) {
    // the propagation grid has the z cells of both halves
    ncoordgrid[1] = 2 * ncoord_model[1];
  }

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1];
  assert_always(ngrid == get_npts_model() * (model_mirror_z ? 2 : 1));

  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  for (int n_rcyl = 0; n_rcyl < ncoordgrid[0]; n_rcyl++) {
    coord_pos_min_tmin[0][n_rcyl] = n_rcyl * globals::rmax / ncoord_model[0];
  }

  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  for (int n_z = 0; n_z < ncoordgrid[1]; n_z++) {
    coord_pos_min_tmin[1][n_z] = globals::rmax * (-1 + (n_z * 2. / ncoordgrid[1]));
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
}

// Set up a 2D spherical or a 3D spherical propagation grid with the cell edges of the model. A 2D spherical grid is
// a 3D spherical grid with one phi cell from 0 to 2 pi, thus one geometry code serves both.
void setup_grid_spherical() {
  assert_always(has_spherical_angle_axes(get_modelgridtype()));
  assert_always(get_propgridtype() == get_modelgridtype());

  // The propagation grid of a mirrored model has the polar cells of both halves. The edges of the half with z < 0
  // are the exact opposites of the model edges, thus the two cells of a mirror pair have the same volume.
  std::vector<double> prop_costheta_edges;
  if (model_mirror_z) {
    for (auto i = std::ssize(model_costheta_edges) - 1; i >= 1; i--) {
      prop_costheta_edges.push_back(-model_costheta_edges[i]);
    }
  }
  prop_costheta_edges.insert(prop_costheta_edges.end(), model_costheta_edges.begin(), model_costheta_edges.end());
  const std::array<std::span<const double>, 3> model_edges{model_vel_r_edges, prop_costheta_edges, model_phi_edges};
  for (int axis = 0; axis < 3; axis++) {
    ncoordgrid[axis] = static_cast<int>(std::ssize(model_edges[axis])) - 1;
    // the last edge of each axis is the upper bound of the grid
    coord_pos_min_tmin[axis].assign(model_edges[axis].begin(), model_edges[axis].end() - 1);
  }
  for (auto& radius_inner_tmin : coord_pos_min_tmin[0]) {
    radius_inner_tmin *= globals::tmin;
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, 1., 2 * PI};

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];
  assert_always(ngrid == get_npts_model() * (model_mirror_z ? 2 : 1));

  // The outer radius of the first shell must exceed the minimum radius that the model format permits.
  const double first_shell_outer_radius = (ncoordgrid[0] > 1) ? coord_pos_min_tmin[0][1] : globals::rmax;
  if (!(first_shell_outer_radius > 1e4 * cellbound_tolerance(0.))) {
    fatal_crash(
        "model.txt: the first radial shell of the spherical grid has the outer radius {:g} [cm] at tmin, but it must "
        "be above {:g} [cm]",
        first_shell_outer_radius, 1e4 * cellbound_tolerance(0.));
  }

  theta_edge_cos = prop_costheta_edges;
  theta_edge_sin.resize(theta_edge_cos.size());
  for (auto i = 0Z; i < std::ssize(theta_edge_cos); i++) {
    // this form of sqrt(1 - cos^2) keeps its precision for an edge near the z axis
    theta_edge_sin[i] = std::sqrt((1. - theta_edge_cos[i]) * (1. + theta_edge_cos[i]));
  }

  phi_edge_cos.resize(model_phi_edges.size());
  phi_edge_sin.resize(model_phi_edges.size());
  for (auto i = 0Z; i < std::ssize(model_phi_edges); i++) {
    phi_edge_cos[i] = std::cos(model_phi_edges[i]);
    phi_edge_sin[i] = std::sin(model_phi_edges[i]);
  }
  // the first and the last azimuthal edge are the half-plane phi = 0, which has exact numbers
  phi_edge_cos.front() = 1.;
  phi_edge_sin.front() = 0.;
  phi_edge_cos.back() = 1.;
  phi_edge_sin.back() = 0.;
}

// for a uniform grid get the extent along the x,y,z coordinate (x_2 - x_1, etc.) at time tmin
// for spherical grid get the radial extent (r_outer - r_inner) at time tmin
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto propcell_width_tmin(const int cellindex, const int axis) -> double {
  return get_cellcoordmax(cellindex, axis) - get_cellcoordmin(cellindex, axis);
}

// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(const int cellindex) -> double {
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D: {
      // volume is constant for uniform grid, so any cell index will do
      const double deltax = get_cellcoordmax(0, 0) - get_cellcoordmin(0, 0);
      // deltax = deltay = deltaz for uniform grid
      return pow3(deltax);
    }

    case GridType::CYLINDRICAL2D: {
      // use cellindex = 0 because spacing is uniform
      const auto delta_z = get_cellcoordmax(0, 1) - get_cellcoordmin(0, 1);
      // delta_rcyl is uniform, but the area of a cylindrical shell is proportional to r^2
      const auto delta_rcylsquared = pow2(get_cellcoordmax(cellindex, 0)) - pow2(get_cellcoordmin(cellindex, 0));
      return delta_z * PI * delta_rcylsquared;
    }

    case GridType::SPHERICAL1D: {
      // the cellindex matters here because the radial spacing is non-uniform
      return 4. / 3. * PI * (pow3(get_cellcoordmax(cellindex, 0)) - pow3(get_cellcoordmin(cellindex, 0)));
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      // the integral of r^2 sin(theta) over the cell. A 2D spherical grid has one phi cell from 0 to 2 pi.
      const double radialpart = (pow3(get_cellcoordmax(cellindex, 0)) - pow3(get_cellcoordmin(cellindex, 0))) / 3.;
      return radialpart * propcell_width_tmin(cellindex, 1) * propcell_width_tmin(cellindex, 2);
    }
  }
  fatal_crash("Unknown propagation grid type {}", static_cast<int>(get_propgridtype()));
}

[[nodiscard]] auto get_propcell_random_xyz_position_tmin(int cellindex, rngstate_type& rngstate) -> Vec3d {
  switch (get_propgridtype()) {
    case GridType::SPHERICAL1D: {
      const double r_inner = get_cellcoordmin(cellindex, 0);
      const double r_outer = get_cellcoordmax(cellindex, 0);
      // use equal volume probability distribution to select radius. A small draw of a thin shell can round to the
      // outer radius, which is the exclusive boundary of the next shell, thus the radius stays below it
      const double radius = std::min(std::cbrt(std::lerp(pow3(r_outer), pow3(r_inner), rng_uniform_pos(rngstate))),
                                     std::nextafter(r_outer, r_inner));

      return vec_scale(get_rand_isotropic_unitvec(rngstate), radius);
    }

    case GridType::CYLINDRICAL2D: {
      const double rcyl_inner = get_cellcoordmin(cellindex, 0);
      const double rcyl_outer = get_cellcoordmax(cellindex, 0);
      // use equal area probability distribution to select radius
      const double rcyl_rand = std::sqrt(std::lerp(pow2(rcyl_outer), pow2(rcyl_inner), rng_uniform_pos(rngstate)));
      const double theta_rand = rng_uniform(rngstate) * 2 * PI;
      const double z_rand =
          std::lerp(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1), rng_uniform_pos(rngstate));
      return {std::cos(theta_rand) * rcyl_rand, std::sin(theta_rand) * rcyl_rand, z_rand};
    }

    case GridType::CARTESIAN3D: {
      Vec3d pos;
      for (int axis = 0; axis < 3; axis++) {
        pos[axis] =
            get_cellcoordmin(cellindex, axis) + (rng_uniform_pos(rngstate) * propcell_width_tmin(cellindex, axis));
      }
      return pos;
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      // The volume element is r^2 dr dcos(theta) dphi, thus r^3, cos(theta), and phi have uniform distributions. A
      // draw near the upper edge of a thin cell can round to that edge, which is the exclusive boundary of the next
      // cell, thus each coordinate stays inside the cell
      const double r_inner = get_cellcoordmin(cellindex, 0);
      const double r_outer = get_cellcoordmax(cellindex, 0);
      const double radius = std::min(std::cbrt(std::lerp(pow3(r_outer), pow3(r_inner), rng_uniform_pos(rngstate))),
                                     std::nextafter(r_outer, r_inner));
      const double costheta_lower = get_cellcoordmin(cellindex, 1);
      const double costheta_upper = get_cellcoordmax(cellindex, 1);
      const double costheta = std::min(std::lerp(costheta_upper, costheta_lower, rng_uniform_pos(rngstate)),
                                       std::nextafter(costheta_upper, costheta_lower));
      const double sintheta = std::sqrt((1. - costheta) * (1. + costheta));
      const double phi_lower = get_cellcoordmin(cellindex, 2);
      const double phi_upper = get_cellcoordmax(cellindex, 2);
      const double phi =
          std::min(std::lerp(phi_lower, phi_upper, rng_uniform_pos(rngstate)), std::nextafter(phi_upper, phi_lower));
      return {radius * sintheta * std::cos(phi), radius * sintheta * std::sin(phi), radius * costheta};
    }
  }
  fatal_crash("Unknown propagation grid type {}", static_cast<int>(get_propgridtype()));
}

DEVICE_FUNC auto get_modelgridtype() -> GridType {
  assert_testmodeonly(model_type.has_value());
  return model_type.value();
}

void setup_spherical_grid_from_edges(const GridType spherical_grid_type, const std::span<const double> vel_r_edges,
                                     const std::span<const double> costheta_edges,
                                     const std::span<const double> phi_edges) {
  assert_always(has_spherical_angle_axes(spherical_grid_type));
  model_type = spherical_grid_type;
  model_vel_r_edges.assign(vel_r_edges.begin(), vel_r_edges.end());
  model_costheta_edges.assign(costheta_edges.begin(), costheta_edges.end());
  if (spherical_grid_type == GridType::SPHERICAL3D) {
    model_phi_edges.assign(phi_edges.begin(), phi_edges.end());
  } else {
    model_phi_edges = {0., 2 * PI};
  }
  set_spherical_model_cell_counts(spherical_grid_type);
  globals::vmax = model_vel_r_edges.back();
  globals::rmax = globals::vmax * globals::tmin;
  require_increasing_shell_radii_at_tmin(model_vel_r_edges);
  setup_grid_spherical();
}

// identify the cell index from an (x,y,z) position and a time.
[[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, const double time) -> int {
  const auto prop_gridtype = get_propgridtype();
  if (has_spherical_angle_axes(prop_gridtype)) {
    const double radius = vec_len(pos);
    if (radius > (globals::vmax * time)) {
      // outside grid
      return -99;
    }
    const double trat = time / globals::tmin;
    const auto [costheta, phi] = get_costheta_phi_from_xyz(pos);
    const int index_r = get_cellindex_from_lower_edges(coord_pos_min_tmin[0], radius / trat);
    const int index_costheta = get_cellindex_from_lower_edges(coord_pos_min_tmin[1], costheta);
    const int index_phi = get_cellindex_from_lower_edges(coord_pos_min_tmin[2], phi);
    const int cellindex = index_r + (ncoordgrid[0] * (index_costheta + (ncoordgrid[1] * index_phi)));

    // the division by trat can round the radius to the other side of a shell boundary
    const double radius_inner = get_cellcoordmin(cellindex, 0) * trat;
    const double radius_outer = get_cellcoordmax(cellindex, 0) * trat;
    assert_always(radius >= (radius_inner - cellbound_tolerance(radius_inner)));
    assert_always(radius <= (radius_outer + cellbound_tolerance(radius_outer)));
    return cellindex;
  }

  const auto ndim = get_ndim(prop_gridtype);
  const auto posgridcoords = get_gridcoords_from_xyz(pos, prop_gridtype);
  int cellindex = 0;
  for (int d = 0; d < ndim; d++) {
    if (std::abs(posgridcoords[d]) > (globals::vmax * time)) {
      // outside grid
      return -99;
    }
    cellindex += get_coordcellindexstride(d) * get_poscoordpointnum(posgridcoords[d], time, d);
  }
  assert_always(cellindex >= 0);
  assert_always(cellindex < ngrid);

  // do a check that the position is within the cell
  const double trat = time / globals::tmin;
  for (int n = 0; n < ndim; n++) {
    const double cellposmin = get_cellcoordmin(cellindex, n) * trat;
    const double cellposmax = get_cellcoordmax(cellindex, n) * trat;
    assert_always(posgridcoords[n] >= (cellposmin - cellbound_tolerance(cellposmin)));
    assert_always(posgridcoords[n] <= (cellposmax + cellbound_tolerance(cellposmax)));
  }
  return cellindex;
}

// Re-establish the invariant that pos lies inside propagation cell cellindex at the given
// time by clamping each grid coordinate into the cell bounds. Called after a boundary
// crossing, this snaps the position exactly onto the crossed cell face(s), so that
// floating-point rounding errors from the position updates (about one ulp per move, i.e.
// up to several cm at the outer grid) cannot accumulate over many crossings.
DEVICE_FUNC void snap_pos_to_cell(Vec3d& pos, const double time, const int cellindex) {
  if (get_propgridtype() != GridType::CARTESIAN3D) {
    return;
  }
  for (int d = 0; d < 3; d++) {
    const double cellposmin = get_cellcoordmin(cellindex, d) / globals::tmin * time;
    const double cellposmax = get_cellcoordmax(cellindex, d) / globals::tmin * time;
    const double newpos_d = std::clamp(pos[d], cellposmin, cellposmax);
    // corrections should only ever be at the floating-point rounding error level
    assert_testmodeonly(std::abs(newpos_d - pos[d]) <= cellbound_tolerance(pos[d]));
    pos[d] = newpos_d;
  }
}

// compute distance to a cell boundary.
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, const double tstart,
                                                 const int cellindex) -> std::tuple<double, int> {
  const auto prop_gridtype = get_propgridtype();
  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    // allocate_nonemptymodelcells() precomputed the radius comparison per cell, because this
    // function runs for every packet step
    if (propcell_outside_escape_surface[cellindex] != 0) {
      return {0., -99};
    }
  }

  // the following vector are in grid coordinates, so either x,y,z (3D) or r (1D), or r_xy, z (2D)

  const auto pktposgridcoord = get_gridcoords_from_xyz(pos, prop_gridtype);

  // dir * CLIGHT_PROP converted from xyz to grid coordinates
  const auto pktvelgridcoord = get_gridcoords_vel_from_xyz_pos_dir(pos, dir, pktposgridcoord, prop_gridtype);

  std::array<int, 3> cellcoordidx{};
  std::array<double, 3> cellcoordmin{};
  std::array<double, 3> cellcoordmax{};
  for (int d = 0; d < get_ndim(prop_gridtype); d++) {
    cellcoordidx[d] = get_cellcoordindex(cellindex, d);
    cellcoordmin[d] = get_cellcoordmin(cellindex, d);
    cellcoordmax[d] = get_cellcoordmax(cellindex, d);
  }

  {
    // the boundaries of an angular axis do not move with the flow, thus the check is for the length axes
    for (int d = 0; d < get_num_length_axes(prop_gridtype); d++) {
      // pos_component_vel_relative_to_flow is constant along a ray with a given direction in Cartesian coordinates, but
      // for non-Cartesian coordinates, we still need to check at the current position whether the packet is
      // moving in the positive or negative direction in each grid coordinate direction relative to the homologous grid
      // flow, otherwise we might never enter the cell that we're supposed to be in
      const bool pos_component_vel_relative_to_flow = (pktvelgridcoord[d] * tstart) > pktposgridcoord[d];

      if constexpr (TESTMODE) {
        bool isoutside_error = false;
        double delta = 0.;
        if (pos_component_vel_relative_to_flow) {
          // check if packet pos is above cell max while moving in the positive direction relative to the grid flow
          const double boundaryposmax = cellcoordmax[d] / globals::tmin * tstart;
          delta = pktposgridcoord[d] - boundaryposmax;
          isoutside_error = pktposgridcoord[d] > (boundaryposmax + cellbound_tolerance(boundaryposmax));
        } else {
          // check if packet pos is below cell min while moving in the negative direction relative to the grid flow
          const double boundaryposmin = cellcoordmin[d] / globals::tmin * tstart;
          delta = pktposgridcoord[d] - boundaryposmin;
          isoutside_error = pktposgridcoord[d] < (boundaryposmin - cellbound_tolerance(boundaryposmin));
        }

        if (isoutside_error) {
          fatal_crash(
              "timestep {}: packet outside coord {} {}{} boundary of cell {} by delta {:g}. vel {:g} initpos {:g} "
              "cellcoordmin {:g} cellcoordmax {:g} dir [{:g}, {:g}, {:g}] tmin {:g} s tstart {:g} s",
              globals::timestep, d, pos_component_vel_relative_to_flow ? '+' : '-', get_coordlabel(prop_gridtype, d),
              cellindex, delta, pktvelgridcoord[d], pktposgridcoord[d], cellcoordmin[d] / globals::tmin * tstart,
              cellcoordmax[d] / globals::tmin * tstart, dir[0], dir[1], dir[2], globals::tmin, tstart);
        }
      }
    }
  }

  double distance = std::numeric_limits<double>::max();
  int next_cellindex{-1};

  // A Cartesian axis d has two possible boundary crossings. For x, the packet trajectory is
  // x = x0 + (dir.x) * c * (t - tstart)
  // the boundaries follow
  // x+/- = x+/-(tmin) * (t/tmin)
  // so the crossing occurs when
  // t - tstart = (x0 - x+/-(tmin)/tmin * tstart) / (x+/-(tmin)/tmin - (dir.x)*c)
  // distance = c * (t - tstart)
  const auto take_cartesian_boundary_if_nearer = [&] [[gnu::always_inline]] (const int d) {
    if (pktvelgridcoord[d] > (cellcoordmax[d] / globals::tmin)) {
      const double d_coordmaxboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                      cellcoordmax[d], tstart)
              ? 0.
              : distance_cartesian_boundary(pktposgridcoord[d], pktvelgridcoord[d], cellcoordmax[d], tstart);

      if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
        distance = d_coordmaxboundary;
        next_cellindex = (cellcoordidx[d] == (ncoordgrid[d] - 1)) ? -99 : cellindex + get_coordcellindexstride(d);
      }
    } else if (pktvelgridcoord[d] < (cellcoordmin[d] / globals::tmin)) {
      const double d_coordminboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                      cellcoordmin[d], tstart)
              ? 0.
              : distance_cartesian_boundary(pktposgridcoord[d], pktvelgridcoord[d], cellcoordmin[d], tstart);

      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        next_cellindex = (cellcoordidx[d] == 0) ? -99 : cellindex - get_coordcellindexstride(d);
      }
    }
  };

  // The radius from the origin is coordinate 0 of each spherical grid. Its two boundaries are spheres that expand.
  const auto take_radial_shell_boundary_if_nearer = [&] [[gnu::always_inline]] () {
    const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_coordmaxboundary =
        is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                    cellcoordmax[0], tstart)
            ? 0.
            : expanding_shell_intersection<BoundaryType::UPPER>(pos, dir, speed, r_outer, tstart);

    // upper d coordinate of the current cell
    if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
      distance = d_coordmaxboundary;
      next_cellindex = (cellcoordidx[0] == (ncoordgrid[0] - 1)) ? -99 : cellindex + get_coordcellindexstride(0);
    }

    const double r_inner = cellcoordmin[0] * tstart / globals::tmin;
    if (r_inner > 0.) {
      const double d_coordminboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                      cellcoordmin[0], tstart)
              ? 0.
              : expanding_shell_intersection<BoundaryType::LOWER>(pos, dir, speed, r_inner, tstart);
      // lower d coordinate of the current cell
      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
      }
    }
  };

  if (prop_gridtype == GridType::SPHERICAL1D) {
    // the only coordinate is the radius from the origin
    take_radial_shell_boundary_if_nearer();
  } else if (has_spherical_angle_axes(prop_gridtype)) {
    // coordinate 0 is the radius, coordinate 1 is cos(theta) of the polar angle, and coordinate 2 is the azimuthal
    // angle phi
    take_radial_shell_boundary_if_nearer();

    const auto angular_crossing =
        get_first_angular_boundary_crossing(pos, dir, cellcoordidx[1], cellcoordidx[2], tstart);
    if ((angular_crossing.distance >= 0.) && (angular_crossing.distance < distance)) {
      distance = angular_crossing.distance;
      next_cellindex =
          cellcoordidx[0] +
          (ncoordgrid[0] * (angular_crossing.index_costheta + (ncoordgrid[1] * angular_crossing.index_phi)));
    }
  } else if (prop_gridtype == GridType::CYLINDRICAL2D) {
    // coordinate 0 is cylindrical radius (distance from z=0 in x-y plane), coord 1 is z

    const std::array<double, 2> posnoz = {pos[0], pos[1]};

    // r_cyl component of direction vector
    const double dirxylen = std::sqrt(pow2(dir[0]) + pow2(dir[1]));
    // r_cyl component of velocity
    const double xyspeed = dirxylen * CLIGHT_PROP;

    if (dirxylen > 0.) {
      // make a normalised 2D direction vector in the xy plane
      const std::array<double, 2> dirnoz = {dir[0] / dirxylen, dir[1] / dirxylen};

      const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
      const double d_rcyl_coordmaxboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                      cellcoordmax[0], tstart)
              ? 0.
              : expanding_shell_intersection<BoundaryType::UPPER>(posnoz, dirnoz, xyspeed, r_outer, tstart);
      if (d_rcyl_coordmaxboundary >= 0.) {
        // how far did the packet travel in the z direction during this time?
        const double d_z_coordmaxboundary = d_rcyl_coordmaxboundary / xyspeed * dir[2] * CLIGHT_PROP;
        // distance from two perpendicular components to the r_cyl upper boundary
        const double d_coordmaxboundary_rcyl = std::sqrt(pow2(d_rcyl_coordmaxboundary) + pow2(d_z_coordmaxboundary));
        if ((d_coordmaxboundary_rcyl >= 0.) && (d_coordmaxboundary_rcyl < distance)) {
          distance = d_coordmaxboundary_rcyl;
          next_cellindex = (cellcoordidx[0] == (ncoordgrid[0] - 1)) ? -99 : cellindex + get_coordcellindexstride(0);
        }
      }

      const double r_inner = cellcoordmin[0] * tstart / globals::tmin;
      // don't try to calculate the intersection if the inner radius is zero
      if (r_inner > 0) {
        // calculate the distance in the xy plane to the inner boundary
        const double d_rcyl_coordminboundary =
            is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                        cellcoordmin[0], tstart)
                ? 0.
                : expanding_shell_intersection<BoundaryType::LOWER>(posnoz, dirnoz, xyspeed, r_inner, tstart);
        if (d_rcyl_coordminboundary >= 0.) {
          const double d_z_coordminboundary = d_rcyl_coordminboundary / xyspeed * dir[2] * CLIGHT_PROP;
          // distance from two perpendicular components to the r_cyl lower boundary
          const double d_coordminboundary_rcyl = std::sqrt(pow2(d_rcyl_coordminboundary) + pow2(d_z_coordminboundary));
          if ((d_coordminboundary_rcyl >= 0.) && (d_coordminboundary_rcyl < distance)) {
            distance = d_coordminboundary_rcyl;
            next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
          }
        }
      }
    } else {
      // The packet is moving exactly along the z axis (reachable when the comoving emission direction
      // is sampled at exactly costheta = +/-1 and the local velocity is parallel to z, or for a virtual
      // packet to an observer at costheta = +/-1), so the
      // general-direction code above would divide by zero. The packet's cylindrical radius stays
      // constant while the grid expands: the receding outer r_cyl boundary can never be crossed, and
      // the expanding inner boundary catches up to the packet at t_cross = rcyl / (rcyl_inner_tmin / tmin)
      const double rcyl_inner_tmin = cellcoordmin[0];
      if (rcyl_inner_tmin > 0.) {
        const double d_coordminboundary_rcyl =
            is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                        cellcoordmin[0], tstart)
                ? 0.
                : ((pktposgridcoord[0] * globals::tmin / rcyl_inner_tmin) - tstart) * CLIGHT_PROP;
        if ((d_coordminboundary_rcyl >= 0.) && (d_coordminboundary_rcyl < distance)) {
          distance = d_coordminboundary_rcyl;
          next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
        }
      }
    }

    // the z axis is Cartesian
    take_cartesian_boundary_if_nearer(1);
  } else if (prop_gridtype == GridType::CARTESIAN3D) {
    for (int d = 0; d < 3; d++) {
      take_cartesian_boundary_if_nearer(d);
    }
  } else {
    assert_always(false);
  }

  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    if (next_cellindex == -1) {
      // No receding cell boundary lies ahead of the packet, so it escapes through the slower
      // spherical escape surface. The comoving coordinates of a packet converge to a point with
      // speed CLIGHT. Only a cell whose fastest corner is above CLIGHT can hold this state. A
      // slower cell here is a sign of a defective boundary intersection, and the run then stops.
      double cornerspeed_squared = 0.;
      for (int d = 0; d < get_num_length_axes(prop_gridtype); d++) {
        cornerspeed_squared += pow2(std::max(std::abs(cellcoordmin[d]), std::abs(cellcoordmax[d])) / globals::tmin);
      }
      assert_always(cornerspeed_squared > pow2(CLIGHT_PROP));

      // Truncation must act on whole cells to keep the volumes and the packet statistics valid.
      // The trapped packet never leaves this cell. The arrival time is constant along a free ray,
      // so a mid-cell escape is valid only for a cell with no matter. A matter cell can trap
      // packets only when the cell diagonal spans from below vmax to above CLIGHT, so a finer
      // grid removes this stop.
      if (get_propcell_modelgridindex(cellindex) >= 0) {
        fatal_crash(
            "a packet cannot reach any boundary of matter cell {}, because every boundary recedes faster than light. "
            "The cell reaches from inside the escape surface to beyond the light speed, so the grid is too coarse. Use "
            "more propagation cells per axis.",
            cellindex);
      }

      const double r_surface = globals::rmax * tstart / globals::tmin;
      if (dot(pos, pos) >= pow2(r_surface)) {
        // A cell that straddles the sphere keeps positions outside the surface. The intersection
        // then lies behind the packet, so it escapes at its current position.
        distance = 0.;
      } else {
        const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised
        // A ray that grazes the surface within rounding has no valid forward intersection and
        // gives -1. The packet is then at the surface itself, so the distance clamps to zero.
        distance = std::max(expanding_shell_intersection<BoundaryType::UPPER>(pos, dir, speed, r_surface, tstart), 0.);
      }
      next_cellindex = -99;
    }
  }

  assert_always((next_cellindex == -99) || ((next_cellindex >= 0) && (next_cellindex < ngrid)));

  const double maxboundarydist =
      (prop_gridtype == GridType::CARTESIAN3D)
          ? std::numbers::sqrt3 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin
          : 2 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin;

  assert_always(distance >= 0.);
  assert_always(distance <= maxboundarydist);

  if (distance > globals::max_path_step) {
    return {globals::max_path_step, cellindex};
  }

  return {distance, next_cellindex};
}

}  // namespace grid
