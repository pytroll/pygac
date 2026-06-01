"""Deprecated -- moved to :mod:`pygac.uncertainty.combine`.

Imported here for backwards compatibility. Will be removed in a future
release. Update your code to import from ``pygac.uncertainty`` instead.
"""
import warnings

from pygac.uncertainty.combine import *  # noqa: F401,F403
from pygac.uncertainty.combine import uncertainty  # noqa: F401

warnings.warn(
    "pygac.calibration.uncertainty has moved to pygac.uncertainty.combine; "
    "import from pygac.uncertainty instead.",
    DeprecationWarning,
    stacklevel=2,
)
