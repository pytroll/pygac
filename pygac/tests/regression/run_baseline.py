"""Generate full-tier regression baselines for every ``*LHRR*`` file.

Usage::

    PYGAC_CONFIG_FILE=etc/pygac.cfg.template \\
    PYGAC_REGRESSION_INPUT_DIR=/home/a001673/Downloads \\
    .venv/bin/python -m pygac.tests.regression.run_baseline

The output directory defaults to ``.pytest_cache/pygac-uncertainty-baselines/``
and is never committed; override with ``PYGAC_REGRESSION_BASELINE_DIR``.
"""

from __future__ import annotations

import logging
import os
import sys
import time
import traceback
from pathlib import Path

from pygac.tests.regression.harness import (
    BaselineEntry,
    write_baseline,
    write_manifest,
)

LOG = logging.getLogger("pygac.regression")


def _input_files(input_dir: Path) -> list[Path]:
    return sorted(p for p in input_dir.iterdir() if "LHRR" in p.name and p.is_file())


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    input_dir = Path(os.environ.get("PYGAC_REGRESSION_INPUT_DIR", ""))
    if not input_dir.is_dir():
        LOG.error("Set PYGAC_REGRESSION_INPUT_DIR to a directory containing *LHRR* files")
        return 2

    baseline_dir = Path(
        os.environ.get(
            "PYGAC_REGRESSION_BASELINE_DIR",
            ".pytest_cache/pygac-uncertainty-baselines",
        )
    )
    baseline_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = baseline_dir / "manifest.json"

    files = _input_files(input_dir)
    if not files:
        LOG.error("No *LHRR* files found in %s", input_dir)
        return 2

    entries: list[BaselineEntry] = []
    failures: list[tuple[str, str]] = []
    for path in files:
        LOG.info("Generating baseline for %s", path.name)
        t0 = time.time()
        try:
            entry = write_baseline(path, baseline_dir)
        except Exception as err:  # noqa: BLE001
            LOG.error("FAILED %s: %s", path.name, err)
            failures.append((path.name, traceback.format_exc()))
            continue
        entries.append(entry)
        LOG.info("  -> %s (%.1fs, sizes=%s)", entry.baseline_path, time.time() - t0, entry.sizes)

    write_manifest(entries, manifest_path)
    LOG.info("Manifest written to %s (%d entries)", manifest_path, len(entries))

    if failures:
        LOG.error("%d file(s) failed:", len(failures))
        for name, tb in failures:
            LOG.error("--- %s ---\n%s", name, tb)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
