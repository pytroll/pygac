"""Per-pixel calibration uncertainty estimation for AVHRR GAC/LAC data.

This package decomposes into three modules:

* :mod:`pygac.uncertainty.ir`      -- thermal-channel uncertainty (3.7/11/12 µm)
* :mod:`pygac.uncertainty.vis`     -- reflective-channel uncertainty (1/2/3a)
* :mod:`pygac.uncertainty.combine` -- merges both into the final dataset

The merged dataset is what :class:`pygac.reader.Reader` attaches to its
calibrated output when ``compute_uncertainties=True``.
"""

from pygac.uncertainty.combine import uncertainty
from pygac.uncertainty.ir import ir_uncertainty
from pygac.uncertainty.vis import vis_uncertainty

__all__ = ["uncertainty", "ir_uncertainty", "vis_uncertainty"]
