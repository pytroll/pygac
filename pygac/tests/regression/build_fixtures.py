"""Build the small committed quick-tier fixtures.

Given a single LHRR file we

1. read & calibrate it through the standard ``Reader`` path,
2. slice the calibrated dataset to a small scan-line / column window,
3. run the *current* uncertainty pipeline on that slice,
4. persist both the sliced calibrated dataset and the expected uncertainty
   output to ``pygac/tests/data/uncertainty_regression/``.

The quick-tier regression test (see
``pygac/tests/test_uncertainty_regression.py``) replays step 3 against the
saved input slice and asserts equality with the saved expected output.

This is *not* a CLI users run by hand on every commit -- it is the
mechanism used to refresh the committed fixtures when (and only when) we
deliberately change behaviour.

Usage::

    PYGAC_CONFIG_FILE=etc/pygac.cfg.template \\
    .venv/bin/python -m pygac.tests.regression.build_fixtures \\
        --input /home/a001673/Downloads/NSS.LHRR.NJ.D00320.S1516.E1527.B3029898.GC \\
        --name noaa14_lac_pod \\
        --scan-slice 600 900 \\
        --column-slice 0 256
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import xarray as xr

from pygac.tests.regression.harness import (
    DEFAULT_TLE_DIR,
    DEFAULT_TLE_NAME,
    compute_calibrated_dataset,
)

LOG = logging.getLogger("pygac.regression.fixtures")

FIXTURES_DIR = Path(__file__).resolve().parents[1] / "data" / "uncertainty_regression"
ENCODING_FLOAT = {"zlib": True, "complevel": 4, "dtype": "float32"}
ENCODING_INT = {"zlib": True, "complevel": 4}


def _slice_dataset(
    ds: xr.Dataset,
    scan_slice: tuple[int, int],
    column_slice: tuple[int, int],
) -> xr.Dataset:
    return ds.isel(
        scan_line_index=slice(*scan_slice),
        columns=slice(*column_slice),
    )


def build_fixture(
    *,
    input_path: Path,
    name: str,
    scan_slice: tuple[int, int],
    column_slice: tuple[int, int],
    tle_dir: str = DEFAULT_TLE_DIR,
    tle_name: str = DEFAULT_TLE_NAME,
    fixtures_dir: Path = FIXTURES_DIR,
) -> tuple[Path, Path]:
    """Write the (input, expected) fixture pair for ``name``."""
    from pygac.uncertainty import uncertainty

    fixtures_dir.mkdir(parents=True, exist_ok=True)

    calibrated_ds, mask, reader_class_name = compute_calibrated_dataset(
        input_path, tle_dir=tle_dir, tle_name=tle_name
    )

    sliced_ds = _slice_dataset(calibrated_ds, scan_slice, column_slice)
    sliced_mask = mask[slice(*scan_slice)]
    sliced_ds = sliced_ds.assign_coords(
        scan_line_mask=("scan_line_index", sliced_mask.astype("uint8"))
    )
    sliced_ds.attrs["reader_class"] = reader_class_name
    sliced_ds.attrs["source_file"] = str(input_path)
    sliced_ds.attrs["scan_slice"] = list(scan_slice)
    sliced_ds.attrs["column_slice"] = list(column_slice)

    expected = uncertainty(sliced_ds, sliced_mask)

    for ds_to_save in (sliced_ds, expected):
        for key, val in list(ds_to_save.attrs.items()):
            if val is None:
                ds_to_save.attrs[key] = "None"

    input_path_out = fixtures_dir / f"{name}.input.nc"
    expected_path_out = fixtures_dir / f"{name}.expected.nc"

    sliced_ds.to_netcdf(input_path_out)
    expected.to_netcdf(expected_path_out)

    LOG.info(
        "Wrote fixture %s (input %.2f MB, expected %.2f MB)",
        name,
        input_path_out.stat().st_size / 1e6,
        expected_path_out.stat().st_size / 1e6,
    )
    return input_path_out, expected_path_out


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--name", required=True)
    p.add_argument("--scan-slice", nargs=2, type=int, required=True, metavar=("START", "STOP"))
    p.add_argument("--column-slice", nargs=2, type=int, required=True, metavar=("START", "STOP"))
    p.add_argument("--tle-dir", default=DEFAULT_TLE_DIR)
    p.add_argument("--tle-name", default=DEFAULT_TLE_NAME)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args(argv)
    build_fixture(
        input_path=args.input,
        name=args.name,
        scan_slice=tuple(args.scan_slice),
        column_slice=tuple(args.column_slice),
        tle_dir=args.tle_dir,
        tle_name=args.tle_name,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
