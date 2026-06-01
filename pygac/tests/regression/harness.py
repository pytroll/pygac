"""Helpers to drive the current uncertainty pipeline and persist baselines.

These helpers deliberately call into ``pygac.calibration.uncertainty`` -- the
same code path used by ``pygac.reader.Reader`` when
``compute_uncertainties=True`` -- so that the regression artefacts produced
here are bit-comparable with what ``pygac-fdr-run --with-uncertainties``
produces.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from pygac.runner import get_reader_class

UNCERTAINTY_VARS: tuple[str, ...] = (
    "random_uncertainty",
    "systematic_uncertainty",
    "channel_covariance_ratio",
    "uncertainty_flags",
)

DEFAULT_TLE_DIR = "gapfilled_tles"
DEFAULT_TLE_NAME = "TLE_%(satname)s.txt"


@dataclass(frozen=True)
class BaselineEntry:
    """Metadata for a single regression baseline."""

    input_path: str
    input_sha256: str
    spacecraft: str
    reader_class: str
    baseline_path: str
    sizes: dict[str, int]


def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compute_calibrated_dataset(
    filename: str | Path,
    *,
    tle_dir: str = DEFAULT_TLE_DIR,
    tle_name: str = DEFAULT_TLE_NAME,
) -> tuple[xr.Dataset, np.ndarray, str]:
    """Read ``filename`` and return ``(calibrated_ds, mask, reader_class_name)``.

    Uncertainties are *not* added to the returned dataset; the caller is
    expected to invoke :func:`pygac.calibration.uncertainty.uncertainty`
    explicitly so that any failure surfaces (the ``Reader`` swallows
    uncertainty exceptions).
    """
    reader_cls = get_reader_class(str(filename))
    reader = reader_cls(
        tle_dir=tle_dir,
        tle_name=tle_name,
        compute_uncertainties=False,
    )
    reader.read(str(filename))
    calibrated_ds = reader.get_calibrated_dataset()
    return calibrated_ds, reader.mask, reader_cls.__name__


def compute_uncertainty_dataset(
    filename: str | Path,
    *,
    tle_dir: str = DEFAULT_TLE_DIR,
    tle_name: str = DEFAULT_TLE_NAME,
) -> tuple[xr.Dataset, str, str]:
    """Run the current uncertainty pipeline against ``filename``.

    Returns ``(uncertainty_ds, spacecraft_name, reader_class_name)``.

    Raises ``RuntimeError`` if any step fails -- we never want a silent
    fallback in the baseline.
    """
    from pygac.calibration.uncertainty import uncertainty

    calibrated_ds, mask, reader_class_name = compute_calibrated_dataset(
        filename, tle_dir=tle_dir, tle_name=tle_name
    )
    spacecraft = calibrated_ds.attrs.get("spacecraft_name", "")
    ucs = uncertainty(calibrated_ds, mask)
    ucs = ucs.rename(
        random="random_uncertainty",
        systematic="systematic_uncertainty",
        chan_covar_ratio="channel_covariance_ratio",
        uncert_flags="uncertainty_flags",
    )
    ucs.attrs["spacecraft_name"] = spacecraft
    ucs.attrs["source_file"] = str(filename)
    return ucs, spacecraft, reader_class_name


def write_baseline(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    tle_dir: str = DEFAULT_TLE_DIR,
    tle_name: str = DEFAULT_TLE_NAME,
) -> BaselineEntry:
    """Generate and persist a baseline for ``input_path``."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ucs, spacecraft, reader_class_name = compute_uncertainty_dataset(
        input_path, tle_dir=tle_dir, tle_name=tle_name
    )
    baseline_path = output_dir / (Path(input_path).name + ".uncertainty.nc")
    encoding = {var: {"zlib": True, "complevel": 4} for var in UNCERTAINTY_VARS}
    ucs.to_netcdf(baseline_path, encoding=encoding)
    return BaselineEntry(
        input_path=str(input_path),
        input_sha256=_sha256(input_path),
        spacecraft=spacecraft,
        reader_class=reader_class_name,
        baseline_path=str(baseline_path),
        sizes={str(k): int(v) for k, v in ucs.sizes.items()},
    )


def write_manifest(entries: list[BaselineEntry], manifest_path: str | Path) -> None:
    """Persist a JSON manifest summarising all baselines generated."""
    payload = {
        "tle_dir": DEFAULT_TLE_DIR,
        "tle_name": DEFAULT_TLE_NAME,
        "entries": [entry.__dict__ for entry in entries],
    }
    Path(manifest_path).write_text(json.dumps(payload, indent=2, sort_keys=True))


def load_manifest(manifest_path: str | Path) -> list[BaselineEntry]:
    """Read back a manifest written by :func:`write_manifest`."""
    payload = json.loads(Path(manifest_path).read_text())
    return [BaselineEntry(**entry) for entry in payload["entries"]]
