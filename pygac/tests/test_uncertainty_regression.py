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

Previously, the uncertainty merge code crashed on KLM files where channel 3a
is present alongside the three IR channels (it tried to write a 3-IR-channel
array into a 2-slot tail). This was fixed in combine.py by deriving
``nb_refl_channels`` from ``ds.sizes["channel_name"] - irdata["random"].shape[-1]``
instead of the buggy ``"3a" in ds["channels"]`` check (which tested data values,
not coordinate labels). ``KNOWN_BROKEN_INPUTS`` is now empty.
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
KNOWN_BROKEN_INPUTS: tuple[str, ...] = ()


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


class TestCombineReflChannelCount:
    """combine.py reflective-channel counting — must handle 5-channel (POD) and 6-channel (KLM) datasets."""

    def _make_ds(self, channel_names):
        """Build a minimal synthetic dataset with the given channel_name coordinate."""
        n_scans, n_pixels = 3, 5
        n_ch = len(channel_names)
        import xarray as xr

        return xr.Dataset(
            {"channels": xr.DataArray(
                np.zeros((n_scans, n_pixels, n_ch), dtype=np.float32),
                dims=["scan_line_index", "columns", "channel_name"],
                coords={"channel_name": channel_names},
            )},
        )

    def _make_uncertainty_output(self, n_scans, n_pixels, n_ch):
        """Return a synthetic (random, systematic, chan_covar_ratio, uncert_flags) dataset."""
        import xarray as xr

        n_ir = 3
        random = xr.DataArray(
            np.zeros((n_scans, n_pixels, n_ch), dtype=np.float32),
            dims=["scan_line_index", "columns", "channel_name"],
        )
        systematic = xr.DataArray(
            np.zeros((n_scans, n_pixels, n_ch), dtype=np.float32),
            dims=["scan_line_index", "columns", "channel_name"],
        )
        chan_covar_ratio = xr.DataArray(
            np.zeros((n_scans, n_pixels, n_ir), dtype=np.float32),
            dims=["scan_line_index", "columns", "ir_channel_name"],
        )
        solar_fov_contam = xr.DataArray(
            np.zeros((n_scans, n_pixels), dtype=np.int8),
            dims=["scan_line_index", "columns"],
        )
        uncert_flags = xr.DataArray(
            np.zeros(n_scans, dtype=np.uint8),
            dims=["scan_line_index"],
        )
        return xr.Dataset(dict(
            random=random,
            systematic=systematic,
            chan_covar_ratio=chan_covar_ratio,
            solar_fov_contam=solar_fov_contam,
            uncert_flags=uncert_flags,
        ))

    def test_five_channel_pod_dataset_does_not_crash(self, monkeypatch):
        """5-channel POD dataset (ch1,2,3b,4,5): 2 reflective + 3 IR → shape (N,P,5)."""
        from pygac.uncertainty import combine

        channel_names = ["1", "2", "3b", "4", "5"]
        ds = self._make_ds(channel_names)
        n_scans, n_pixels = ds.sizes["scan_line_index"], ds.sizes["columns"]
        ir_out = self._make_uncertainty_output(n_scans, n_pixels, 3)
        vis_out = self._make_uncertainty_output(n_scans, n_pixels, 2)

        monkeypatch.setattr(combine, "ir_uncertainty", lambda ds, mask: ir_out)
        monkeypatch.setattr(combine, "vis_uncertainty", lambda ds, mask: vis_out)

        result = combine.uncertainty(ds, mask=None)
        assert result["random"].shape == (n_scans, n_pixels, len(channel_names))

    def test_six_channel_klm_dataset_does_not_crash(self, monkeypatch):
        """6-channel KLM dataset (ch1,2,3a,3b,4,5): 3 reflective + 3 IR → shape (N,P,6)."""
        from pygac.uncertainty import combine

        channel_names = ["1", "2", "3a", "3b", "4", "5"]
        ds = self._make_ds(channel_names)
        n_scans, n_pixels = ds.sizes["scan_line_index"], ds.sizes["columns"]
        ir_out = self._make_uncertainty_output(n_scans, n_pixels, 3)
        vis_out = self._make_uncertainty_output(n_scans, n_pixels, 3)

        monkeypatch.setattr(combine, "ir_uncertainty", lambda ds, mask: ir_out)
        monkeypatch.setattr(combine, "vis_uncertainty", lambda ds, mask: vis_out)

        result = combine.uncertainty(ds, mask=None)
        assert result["random"].shape == (n_scans, n_pixels, len(channel_names))


class TestCombineFlagAssembly:
    """combine.py flag loop: IR per-scanline bit OR-ed with VIS per-pixel bit."""

    def _make_inputs(self, n=8, p=5):
        ir_flags = np.zeros(n, dtype=np.uint8)
        ir_flags[2] = 1   # bad space view
        solar = np.zeros((n, p), dtype=np.int8)
        solar[4, 3] = 1   # solar contamination of FOV
        return ir_flags, solar

    def test_ir_flag_broadcast_to_all_pixels(self):
        from pygac.uncertainty.combine import _assemble_flags
        ir_flags, solar = self._make_inputs()
        result = _assemble_flags(ir_flags, solar)
        assert np.all(result[2, :] & 1)

    def test_solar_bit_set_on_contaminated_pixel(self):
        from pygac.uncertainty.combine import _assemble_flags
        ir_flags, solar = self._make_inputs()
        result = _assemble_flags(ir_flags, solar)
        assert result[4, 3] & 8
        assert not (result[4, 0] & 8)

    def test_clean_scanline_has_zero_flags(self):
        from pygac.uncertainty.combine import _assemble_flags
        ir_flags, solar = self._make_inputs()
        result = _assemble_flags(ir_flags, solar)
        assert result[0, 0] == 0
