"""Regression tests for the uncertainty pipeline.

Quick tier
----------

Re-runs :func:`pygac.uncertainty.uncertainty` against the small
committed fixtures under ``pygac/tests/data/uncertainty_regression/`` and
asserts the output is *exactly* equal to the captured baseline.

Full tier
---------

Optional. Walks the per-orbit baselines under
``$PYGAC_REGRESSION_BASELINE_DIR`` produced by
``python -m pygac.tests.regression.run_baseline`` and re-runs the whole
pipeline against the original L1B inputs. Marked with
``pytest.mark.regression_full`` and skipped unless
``--run-regression-full`` is passed (see ``conftest.py``).

Known-broken inputs
-------------------

The current uncertainty merge code crashes on KLM files where channel 3a is
present alongside the three IR channels (it tries to write a 3-IR-channel
array into a 2-slot tail). The affected files are documented in
:data:`KNOWN_BROKEN_INPUTS` and covered by an xfail test so that the day
the refactor fixes them, the xfail flips to a hard pass and they enter the
baseline.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pygac.tests.regression.harness import (
    UNCERTAINTY_VARS,
    compute_uncertainty_dataset,
    load_manifest,
)

FIXTURES_DIR = Path(__file__).parent / "data" / "uncertainty_regression"

#: Filenames for which the current uncertainty pipeline raises before
#: producing a baseline. Tracked here so the refactor can flip them to
#: passing inputs as soon as the underlying bug is fixed.
KNOWN_BROKEN_INPUTS: tuple[str, ...] = (
    "ESR.LHRR.NK.D00117.S0546.E0558.B1014848.BL",
    "ESR.LHRR.NK.D03338.S1759.E1812.B2889797.BN",
    "ESR.LHRR.NN.D10329.S1036.E1046.B2841818.MT",
)


def _discover_fixtures() -> list[str]:
    if not FIXTURES_DIR.is_dir():
        return []
    return sorted(p.stem.removesuffix(".input") for p in FIXTURES_DIR.glob("*.input.nc"))


@pytest.mark.parametrize("fixture_name", _discover_fixtures())
def test_quick_uncertainty_regression(fixture_name: str) -> None:
    """Run uncertainty() against a committed input slice and compare with the saved expected output."""
    from pygac.uncertainty import uncertainty

    input_path = FIXTURES_DIR / f"{fixture_name}.input.nc"
    expected_path = FIXTURES_DIR / f"{fixture_name}.expected.nc"

    with xr.open_dataset(input_path) as ds_in:
        ds = ds_in.load()
    with xr.open_dataset(expected_path) as ds_exp:
        expected = ds_exp.load()

    mask = ds["scan_line_mask"].values.astype(bool)
    actual = uncertainty(ds, mask)

    _assert_uncertainty_equal(actual, expected, tolerance="exact")


@pytest.mark.regression_full
def test_full_uncertainty_regression() -> None:
    """Regenerate uncertainty for every baseline in $PYGAC_REGRESSION_BASELINE_DIR and compare."""
    baseline_dir = Path(
        os.environ.get(
            "PYGAC_REGRESSION_BASELINE_DIR",
            ".pytest_cache/pygac-uncertainty-baselines",
        )
    )
    manifest_path = baseline_dir / "manifest.json"
    if not manifest_path.is_file():
        pytest.skip(f"No baseline manifest at {manifest_path}; run run_baseline first")
    entries = load_manifest(manifest_path)
    if not entries:
        pytest.skip("Empty baseline manifest")

    failures: list[str] = []
    for entry in entries:
        with xr.open_dataset(entry.baseline_path) as baseline_ds:
            baseline = baseline_ds.load()
        actual, *_ = compute_uncertainty_dataset(entry.input_path)
        try:
            _assert_uncertainty_equal(actual, baseline, tolerance="exact")
        except AssertionError as err:
            failures.append(f"{Path(entry.input_path).name}: {err}")
    if failures:
        pytest.fail("Regression mismatches:\n" + "\n".join(failures))


@pytest.mark.regression_full
@pytest.mark.parametrize("filename", KNOWN_BROKEN_INPUTS)
@pytest.mark.xfail(strict=True, reason="merge code can't handle 3a + 3 IR channels; tracked for refactor")
def test_known_broken_klm_inputs(filename: str) -> None:
    """Document the 3a-KLM merge bug; flips to passing once the refactor fixes it."""
    input_dir = Path(os.environ.get("PYGAC_REGRESSION_INPUT_DIR", ""))
    if not (input_dir / filename).is_file():
        pytest.skip(f"{filename} not available in {input_dir}")
    compute_uncertainty_dataset(input_dir / filename)


def _assert_uncertainty_equal(actual: xr.Dataset, expected: xr.Dataset, *, tolerance: str) -> None:
    """Compare two uncertainty datasets with the given tolerance policy."""
    actual = _normalise(actual)
    expected = _normalise(expected)

    missing = set(UNCERTAINTY_VARS) - set(actual.data_vars)
    assert not missing, f"actual is missing variables: {missing}"

    for var in UNCERTAINTY_VARS:
        assert actual[var].dims == expected[var].dims, (
            f"{var}: dims {actual[var].dims} != {expected[var].dims}"
        )
        assert actual[var].shape == expected[var].shape, (
            f"{var}: shape {actual[var].shape} != {expected[var].shape}"
        )
        assert actual[var].dtype == expected[var].dtype, (
            f"{var}: dtype {actual[var].dtype} != {expected[var].dtype}"
        )

    if tolerance == "exact":
        for var in UNCERTAINTY_VARS:
            np.testing.assert_array_equal(
                actual[var].values, expected[var].values,
                err_msg=f"{var} mismatch",
            )
    else:  # pragma: no cover - placeholder for arithmetic-changing phases
        raise NotImplementedError(tolerance)


def _normalise(ds: xr.Dataset) -> xr.Dataset:
    """Map the legacy variable names used by the baseline to the canonical ones."""
    rename = {
        "random": "random_uncertainty",
        "systematic": "systematic_uncertainty",
        "chan_covar_ratio": "channel_covariance_ratio",
        "uncert_flags": "uncertainty_flags",
    }
    present = {k: v for k, v in rename.items() if k in ds.variables}
    return ds.rename(present) if present else ds
