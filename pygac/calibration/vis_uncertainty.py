"""Deprecated -- moved to :mod:`pygac.uncertainty.vis`.

Imported here for backwards compatibility. Will be removed in a future
release. Update your code to import from ``pygac.uncertainty.vis`` instead.
"""
import warnings

from pygac.uncertainty.vis import *  # noqa: F401,F403
from pygac.uncertainty.vis import vis_uncertainty  # noqa: F401

warnings.warn(
    "pygac.calibration.vis_uncertainty has moved to pygac.uncertainty.vis; "
    "import from pygac.uncertainty.vis instead.",
    DeprecationWarning,
    stacklevel=2,
)
