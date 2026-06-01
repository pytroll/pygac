"""Deprecated -- moved to :mod:`pygac.uncertainty.ir`.

Imported here for backwards compatibility. Will be removed in a future
release. Update your code to import from ``pygac.uncertainty.ir`` instead.
"""
import warnings

from pygac.uncertainty.ir import *  # noqa: F401,F403
from pygac.uncertainty.ir import (  # noqa: F401
    allan_deviation,
    get_bad_space_counts,
    get_uncert_parameter_thresholds,
    ir_uncertainty,
)

warnings.warn(
    "pygac.calibration.ir_uncertainty has moved to pygac.uncertainty.ir; "
    "import from pygac.uncertainty.ir instead.",
    DeprecationWarning,
    stacklevel=2,
)
