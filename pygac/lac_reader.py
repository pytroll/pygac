"""The LAC reader."""

import logging

import numpy as np

try:
    from pyorbital.geoloc_instrument_definitions import avhrr_from_times
except ImportError:
    # In pyorbital 1.9.2 and earlier avhrr_gac returned LAC geometry
    from pyorbital.geoloc_instrument_definitions import avhrr_gac as avhrr_from_times

from pygac.reader import Reader, ReaderError

LOG = logging.getLogger(__name__)


class LACReader(Reader):
    """Reader for LAC data."""

    # Scanning frequency (scanlines per millisecond)
    scan_freq = 6.0 / 1000.0
    # Max scanlines
    max_scanlines = 65535
    lonlat_sample_points = np.arange(24, 2048, 40)

    def __init__(self, *args, **kwargs):
        """Init the LAC reader."""
        super(LACReader, self).__init__(*args, **kwargs)
        self.scan_width = 2048
        self.geoloc_definition = avhrr_from_times

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader."""
        # call super to enter the Method Resolution Order (MRO)
        super(LACReader, cls)._validate_header(header)
        LOG.debug("validate header")
        data_set_name = header["data_set_name"].decode()
        # split header into parts
        creation_site, transfer_mode, platform_id = (
            data_set_name.split(".")[:3])
        if transfer_mode not in ["LHRR", "HRPT", "FRAC"]:
            raise ReaderError('Improper transfer mode "%s"!' % transfer_mode)
