"""Mutation sentinel test for the uncertainty pipeline (Phase 3.0).

Phase 3 of the uncertainty refactor introduces dataclasses that hold
references to per-scanline numpy arrays pulled from the input
:class:`xr.Dataset` via ``.values``. ``ds[var].values`` returns a *view*
into the backing array, not a copy. Several places in the current code
write into those views in-place (e.g.
``prt[gd] = 0; prt[ifix] = np.interp(...)`` in ``get_vars``). If that
mutation propagates back into ``ds``, downstream callers (or repeat calls
within the same pipeline run) will see different inputs than they
originally had, which is both a latent bug and a landmine for any
refactor that changes call ordering.

This module pins two contracts:

1. End-to-end: :func:`pygac.uncertainty.uncertainty` does not mutate any
   variable in the input dataset (or the mask).
2. Per-helper: :func:`pygac.uncertainty.ir.get_vars` does not mutate
   ``ds["mean_prt_counts"]`` even when the bad-PRT interpolation branch
   fires. This is the specific landmine that motivated the audit -- on
   our committed fixtures the branch happens to never fire, so the
   end-to-end test alone misses it; the targeted helper test triggers
   it by lowering the threshold.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pygac.uncertainty import uncertainty
from pygac.uncertainty.ir import get_uncert_parameter_thresholds, get_vars

FIXTURES_DIR = Path(__file__).parent / "data" / "uncertainty_regression"


def _discover_fixtures() -> list[str]:
    if not FIXTURES_DIR.is_dir():
        return []
    return sorted(p.stem.removesuffix(".input") for p in FIXTURES_DIR.glob("*.input.nc"))


@pytest.mark.parametrize("fixture_name", _discover_fixtures())
def test_uncertainty_does_not_mutate_input_dataset(fixture_name: str) -> None:
    """uncertainty(ds, mask) must not modify any variable in ds in place."""
    input_path = FIXTURES_DIR / f"{fixture_name}.input.nc"

    with xr.open_dataset(input_path) as ds_in:
        ds = ds_in.load()

    snapshot = {name: ds[name].values.copy() for name in ds.data_vars}
    mask = ds["scan_line_mask"].values.astype(bool)
    mask_snapshot = mask.copy()

    uncertainty(ds, mask)

    mutated = [name for name, expected in snapshot.items()
               if not _arrays_byte_equal(ds[name].values, expected)]

    assert not mutated, (
        f"uncertainty() mutated input variables in place: {mutated}. "
        "Phase 3 refactor relies on inputs being read-only."
    )
    assert _arrays_byte_equal(mask, mask_snapshot), (
        "uncertainty() mutated the mask array in place."
    )


@pytest.mark.parametrize("fixture_name", _discover_fixtures())
def test_get_vars_does_not_mutate_prt_when_interpolation_fires(fixture_name: str) -> None:
    """get_vars must not mutate ds["mean_prt_counts"] via .values views.

    The committed fixtures contain already-interpolated PRT data, so the
    ``prt[ifix] = np.interp(...)`` branch in get_vars never fires on them
    naturally. We force the branch by injecting four consecutive
    below-threshold values so every PRT group (iprt 1..4) gets at least
    one bad value to interpolate over.
    """
    input_path = FIXTURES_DIR / f"{fixture_name}.input.nc"
    with xr.open_dataset(input_path) as ds_in:
        ds = ds_in.load()

    window, _, _, prt_threshold, ict_threshold, space_threshold = \
        get_uncert_parameter_thresholds()

    # Inject NaN values to trigger the prt[gd]=0 branch in get_vars (line ~738).
    # We could inject low values to trigger the iprt-aware interpolation branch,
    # but get_prt_nos excludes any value below threshold from PRT numbering
    # (assigns iprt=0), so the interpolation branch never sees them.
    new_prt = ds["mean_prt_counts"].values.copy()
    inject_at = slice(100, 104)
    new_prt[inject_at] = np.nan
    ds = ds.assign(mean_prt_counts=(ds["mean_prt_counts"].dims, new_prt))

    from pygac.calibration.noaa import Calibrator
    cal = Calibrator(ds.attrs["spacecraft_name"])
    mask = ds["scan_line_mask"].values.astype(bool)

    snapshot = ds["mean_prt_counts"].values.copy()

    get_vars(
        ds, 0, None, window,
        prt_threshold, ict_threshold, space_threshold,
        True, cal, mask, out_prt=True,
    )

    assert _arrays_byte_equal(ds["mean_prt_counts"].values, snapshot), (
        "get_vars() mutated ds['mean_prt_counts'] in place when the "
        "interpolation branch fired."
    )


def _arrays_byte_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """Strict equality: same shape, same dtype, NaN-aware element-wise equality."""
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if np.issubdtype(a.dtype, np.floating):
        both_nan = np.isnan(a) & np.isnan(b)
        return bool(np.all(both_nan | (a == b)))
    return bool(np.array_equal(a, b))
