"""Shared pytest configuration for the pygac test suite."""

from __future__ import annotations

import os


def pytest_addoption(parser):
    parser.addoption(
        "--run-regression-full",
        action="store_true",
        default=False,
        help=(
            "Run the full-tier uncertainty regression suite against the LHRR "
            "baselines in $PYGAC_REGRESSION_BASELINE_DIR. The quick tier always "
            "runs; this flag opts in to the (multi-GB, multi-minute) end-to-end "
            "comparison."
        ),
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-regression-full"):
        return
    if os.environ.get("PYGAC_REGRESSION_FULL"):
        return
    import pytest

    skip_full = pytest.mark.skip(
        reason="full-tier uncertainty regression: pass --run-regression-full to enable"
    )
    for item in items:
        if "regression_full" in item.keywords:
            item.add_marker(skip_full)
