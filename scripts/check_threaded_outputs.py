#!/usr/bin/env python3
# /// script
# requires-python = ">=3.14"
# ///
"""Check physical bounds in the outputs of the threaded test models."""

import math
import sys
from collections import defaultdict
from pathlib import Path

from compression import zstd


def require(condition: bool, message: str) -> None:
    """Stop the check when a required condition fails."""
    if not condition:
        raise ValueError(message)


def read_rows(path: Path, *, header: bool = False) -> list[list[float]]:
    """Read numeric rows from a plain file or a zstd file."""
    opener = zstd.open if path.suffix == ".zst" else Path.open
    with opener(path, mode="rt", encoding="utf-8") as file:
        text = file.read()
    rows = [[float(value) for value in line.split()] for line in text.splitlines()[int(header) :] if line.strip()]
    require(all(math.isfinite(value) for row in rows for value in row), f"{path}: a value is not finite")
    return rows


def check_light_curve(modelpath: Path) -> None:
    """Require finite, non-negative luminosities and at least one positive luminosity."""
    paths = list(modelpath.glob("light_curve.out*"))
    require(len(paths) == 1, "The model must have one light curve")
    rows = read_rows(paths[0])
    finish = int((modelpath / "input-newrun.txt").read_text(encoding="utf-8").splitlines()[2].split()[1])
    require(len(rows) == finish and all(len(row) == 3 for row in rows), "The light curve is incomplete")
    require(all(row[0] > 0 and min(row[1:]) >= 0 for row in rows), "The light curve has a negative value")
    require(any(row[1] > 0 for row in rows), "The light curve has no positive luminosity")


def check_nlte(modelpath: Path) -> None:
    """Check NLTE populations and the sum of the level fractions in each ion."""
    ion_fractions: defaultdict[tuple[float, ...], float] = defaultdict(float)
    has_positive_nlte_population = False
    for path in modelpath.glob("job_from_ts0000/nlte_*.out*"):
        for row in read_rows(path, header=True):
            require(len(row) == 8, f"{path}: an NLTE row has the wrong number of columns")
            require(min(row[5:]) >= 0 and row[7] <= 1.00001, f"{path}: an NLTE population is outside its bounds")
            ion_fractions[tuple(row[:4])] += row[7]
            has_positive_nlte_population |= row[6] > 0
    require(bool(ion_fractions), "The run has no NLTE populations")
    require(has_positive_nlte_population, "The run has no positive NLTE population")
    # The ion population includes autoionising levels that these rows omit. The fraction sum can be less than one.
    # The output has six significant digits. Permit a small error in the sum of the rounded fractions.
    require(all(total <= 1.0005 for total in ion_fractions.values()), "The NLTE level fractions exceed one")


def check_virtual_spectra(modelpath: Path) -> None:
    """Require a positive Stokes I flux and physical Stokes Q and U fluxes."""
    paths = list(modelpath.glob("vspecpol/vspecpol_*.out*"))
    require(len(paths) == 2, "The run must have one virtual spectrum for each of the two ranks")
    positive_flux = False
    for path in paths:
        for row in read_rows(path):
            require(len(row) > 1 and (len(row) - 1) % 3 == 0, f"{path}: a virtual spectrum row is incomplete")
            if row[0] == 0:
                continue
            timebins = (len(row) - 1) // 3
            for timestep in range(timebins):
                flux_i, flux_q, flux_u = (row[1 + timestep + stokes * timebins] for stokes in range(3))
                require(flux_i >= 0, f"{path}: Stokes I is negative")
                require(math.hypot(flux_q, flux_u) <= flux_i * 1.00002, f"{path}: polarisation exceeds Stokes I")
                positive_flux |= flux_i > 0
    require(positive_flux, "The virtual spectra have no positive Stokes I flux")


def main() -> None:
    """Select the physical checks for one test model."""
    modelpath = Path(sys.argv[1])
    check_light_curve(modelpath)
    if modelpath.name.startswith("nebular_"):
        check_nlte(modelpath)
    if modelpath.name.startswith("classicmode_"):
        check_virtual_spectra(modelpath)
    print(f"The numerical checks passed for {modelpath.name}")


if __name__ == "__main__":
    main()
