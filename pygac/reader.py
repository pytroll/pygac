#!/usr/bin/env python

# Copyright (c) 2014, 2019 Pygac Developers

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Carlos Horn <carlos.horn@external.eumetsat.int>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Generic reader for GAC and LAC data.

Can't be used as is, has to be subclassed to add specific read functions.
"""

import datetime
import logging
import os
import re
import types
import warnings
from abc import ABC, abstractmethod
from contextlib import suppress
from functools import cached_property
from importlib.metadata import entry_points

import geotiepoints as gtp
import numpy as np
import pyorbital
import xarray as xr
from packaging.version import Version
from pyorbital import astronomy
from pyorbital.geoloc import compute_pixels, get_lonlatalt
from pyorbital.orbital import Orbital

from pygac import gac_io
from pygac.utils import calculate_sun_earth_distance_correction, centered_modulus, get_absolute_azimuth_angle_diff

LOG = logging.getLogger(__name__)

# rpy values from
# here:http://yyy.rsmas.miami.edu/groups/rrsl/pathfinder/Processing/proc_app_a.html
rpy_coeffs = {
    "noaa7": {
        "roll": 0.000,
        "pitch": 0.000,
        "yaw": 0.000,
    },
    "noaa9": {
        "roll": 0.000,
        "pitch": 0.0025,
        "yaw": 0.000,
    },
    "noaa10": {
        "roll": 0.000,
        "pitch": 0.000,
        "yaw": 0.000,
    },
    "noaa11": {
        "roll": -0.0019,
        "pitch": -0.0037,
        "yaw": 0.000,
    },
    "noaa12": {
        "roll": 0.000,
        "pitch": 0.000,
        "yaw": 0.000,
    },
    "noaa14": {
        "roll": 0.000,
        "pitch": 0.000,
        "yaw": 0.000,
    },
}


class ReaderError(ValueError):
    """Raised in Reader.read if the given file does not correspond to it."""

    pass


class NoTLEData(IndexError):
    """Raised if no TLE data available within time range."""


class DecodingError(ValueError):
    """Raised when decoding of some value fails."""


class TimestampMismatch(IndexError):
    """Raise when matching timestamps doesn't work."""


class Reader(ABC):
    """Reader for GAC and LAC format, POD and KLM platforms."""

    # data set header format, see _validate_header for more details
    data_set_pattern = re.compile(r"\w{3}\.\w{4}\.\w{2}.D\d{5}\.S\d{4}\.E\d{4}\.B\d{7}\.\w{2}")

    def __init__(
        self,
        interpolate_coords=True,
        adjust_clock_drift=True,
        tle_dir=None,
        tle_name=None,
        tle_thresh=7,
        creation_site=None,
        custom_calibration=None,
        calibration_file=None,
        header_date="auto",
        calibration_method="noaa",
        calibration_parameters=None,
        correct_scanlines=True,
        reference_image=None,
        dem=None,
        compute_lonlats_from_tles: bool = False,
        compute_uncertainties: bool = False,
    ):
        """Init the reader.

        Args:
            interpolate_coords: Interpolate coordinates from tiepoint grid
                to all pixels.
            adjust_clock_drift: Adjust the geolocation to compensate for the
                clock error (POD satellites only).
            tle_dir: Directory holding TLE files
            tle_name: Filename pattern of TLE files.
            tle_thresh: Maximum number of days between observation and nearest
                TLE
            creation_site: The three-letter identifier of the creation site (eg 'NSS')
            custom_calibration: dictionary with a subset of user defined satellite specific
                                calibration coefficients
            calibration_file: path to json file containing default calibrations
            header_date: the date to use for pod header choice. Defaults to "auto".
            correct_scanlines: Remove corrrupt scanline numbers. Defaults to True
            reference_image: the reference image to use for georeferencing
            dem: the digital elevation model to use for orthocorrection
            compute_lonlats_from_tles: Do not use the longitudes and latitudes provided in the file, rather compute them
                                       from the TLE.
            compute_uncertainties: Whether to add uncertainty estimates in the calibrated_dataset.

        """
        self.meta_data = {}
        self.interpolate_coords = interpolate_coords
        self.adjust_clock_drift = adjust_clock_drift
        self.tle_dir = tle_dir
        self.tle_name = tle_name
        self.tle_thresh = tle_thresh
        self.creation_site = (creation_site or "NSS").encode("utf-8")
        self.header_date = header_date
        self.correct_scanlines = correct_scanlines
        self.head = None
        self.scans = None
        self.spacecraft_name = None
        self.spacecraft_id = None
        self._times_as_np_datetime64 = None
        self.lats = None
        self.lons = None
        self.tle_lines = None
        self.filename = None
        self._mask = None
        self._rpy = None
        self.reference_image = reference_image
        self.dem = dem
        self.compute_lonlats_from_tles: bool = compute_lonlats_from_tles
        self.compute_uncertainties: bool = compute_uncertainties

        self.clock_drift_correction_applied = False

        if calibration_method.lower() not in ["noaa"]:
            raise ValueError(f"Unknown calibration method {calibration_method}.")
        self.calibration_method = calibration_method
        self.calibration_parameters = calibration_parameters or dict()
        if custom_calibration is not None or calibration_file is not None:
            warnings.warn(
                "`custom_calibration` and `calibration_file` parameters will be deprecated soon. "
                "Please use `calibration_parameters` dictionary instead.",
                PendingDeprecationWarning,
            )
            if not self.calibration_parameters:
                self.calibration_parameters = dict(custom_coeffs=custom_calibration, coeffs_file=calibration_file)

    @property
    def times(self):
        """Get the UTCs as datetime.datetime."""
        return self.to_datetime(self._times_as_np_datetime64)

    @property
    def filename(self):
        """Get the property 'filename'."""
        return self._filename

    @filename.setter
    def filename(self, filepath):
        """Set the property 'filename'."""
        if filepath is None:
            self._filename = None
        else:
            filepath = os.fspath(filepath)
            match = self.data_set_pattern.search(filepath)
            if match:
                self._filename = match.group()
            else:
                self._filename = os.path.basename(filepath)

    @abstractmethod
    def read(self, filename, fileobj=None):  # pragma: no cover
        """Read the GAC/LAC data.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)
        """
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def read_header(cls, filename, fileobj=None):  # pragma: no cover
        """Read the file header.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)

        Returns:
            archive_header (struct): archive header
            header (struct): file header

        Note:
            This is a classmethod to avoid throwaway instances while
            checking if the reader corresponds to the input file.
        """
        raise NotImplementedError

    @classmethod
    def _correct_data_set_name(cls, header, filename):
        """Replace invalid data_set_name from header with filename.

        Args:
            header (struct): file header
            filename (str): path to file
        """
        filename = str(filename)
        data_set_name = header["data_set_name"]
        try:
            header["data_set_name"] = cls._decode_data_set_name(data_set_name)
        except DecodingError:
            LOG.debug(f"The data_set_name in header {header['data_set_name']} does not match. Use filename instead.")
            match = cls.data_set_pattern.search(filename)
            if match:
                data_set_name = match.group()
                LOG.debug(f"Set data_set_name, to filename {data_set_name}")
                header["data_set_name"] = data_set_name.encode()
            else:
                LOG.debug(f"header['data_set_name']={header['data_set_name']}; filename='{filename}'")
                raise ReaderError("Cannot determine data_set_name!")
        return header

    @classmethod
    def _decode_data_set_name(cls, data_set_name):
        for encoding in "utf-8", "cp500":
            try:
                data_set_name = cls._decode_data_set_name_for_encoding(data_set_name, encoding)
            except DecodingError as err:
                LOG.debug(str(err))
            else:
                return data_set_name
        else:
            raise DecodingError("Could not reliably decode the dataset name.")

    @classmethod
    def _decode_data_set_name_for_encoding(cls, data_set_name, encoding):
        data_set_name = data_set_name.decode(encoding, errors="ignore")
        if not cls.data_set_pattern.match(data_set_name):
            raise DecodingError(
                f"The data_set_name in header {data_set_name} does not seem correct using encoding {encoding}."
            )
        else:
            data_set_name = data_set_name.encode()
        return data_set_name

    @classmethod
    def _validate_header(cls, header):
        """Check if the header belongs to this reader.

        Note:
            according to https://www1.ncdc.noaa.gov/pub/data/satellite/
            publications/podguides/TIROS-N%20thru%20N-14/pdf/NCDCPOD2.pdf
            and https://www1.ncdc.noaa.gov/pub/data/satellite/
            publications/podguides/N-15%20thru%20N-19/pdf/
            2.5%20Section%208.0%20NOAA%20Level%201B%20Database.pdf
            the data set name splits into
            PROCESSING-CENTER.DATA-TYPE.SPACECRAFT-UNIQUE-ID.
            YEAR-DAY.START-TIME.STOP-TIME.PROCESSING-BLOCK-ID.SOURCE
            This should be sufficient information to determine the
            reader.
            Global Area Coverage (GAC):
                DATA-TYPE = GHRR
            Local Area Coverage (LAC):
                DATA-TYPE = LHRR
            Polar Orbiter Data (POD):
                SPACECRAFT-UNIQUE-ID in [TN, NA, NB, NC, ND, NE, NF,
                                         NG, NH, NI, NJ]
            NOAA-K, -L, -M system, but also newer satellites (KLM):
                SPACECRAFT-UNIQUE-ID in [NK, NL, NM, NN, NP, M2, M1]
        """
        # This method does not need to be implemented in all subclasses.
        # It is intended for cooperative multiple inheritance, i.e.
        # each child class which implements this method, should call the
        # super method to enter into the method resolution order.
        # See https://docs.python.org/3/library/functions.html#super
        # second use case "diamond diagrams".
        # Check if the data set name matches the pattern
        LOG.debug("validate header")
        data_set_name = header["data_set_name"].decode(errors="ignore")
        if not cls.data_set_pattern.match(data_set_name):
            raise ReaderError("Data set name %s does not match!" % header["data_set_name"])

    def _read_scanlines(self, buffer, count):
        """Read the scanlines from the given buffer.

        Args:
            buffer (bytes, bytearray): buffer to read from
            count (int): number of expected scanlines
        """
        # Calculate the actual number of complete scanlines. The integer divisoin
        # may strip a potentially incomplete line at the end of the file.
        line_count = len(buffer) // self.scanline_type.itemsize
        if line_count != count:
            LOG.warning("Expected %d scan lines, but found %d!" % (count, line_count))
            warnings.warn("Unexpected number of scanlines!", category=RuntimeWarning, stacklevel=2)
        self.scans = np.frombuffer(buffer, dtype=self.scanline_type, count=line_count)

    @classmethod
    def can_read(cls, filename, fileobj=None):
        """Read the GAC/LAC data.

        Args:
            filename (str): Path to GAC/LAC file
            fileobj: An open file object to read from. (optional)

        Retruns:
            result (bool): True if the reader can read the input
        """
        if fileobj:
            pos = fileobj.tell()
        try:
            archive_header, header = cls.read_header(filename, fileobj=fileobj)
            result = True
        except (ReaderError, ValueError) as exception:
            LOG.debug("%s failed to read the file! %s" % (cls.__name__, repr(exception)))
            result = False
        finally:
            if fileobj:
                fileobj.seek(pos)
        return result

    @classmethod
    def fromfile(cls, filename, fileobj=None):
        """Create Reader from file, alternative constructor.

        Args:
            filename (str): Path to GAC/LAC file

        Kwargs:
            fileobj (file object): Open file object to read from

        Note:
            The fileobj is useful when dealing with tar archives,
            where the filename is given by the tarinfo object, but
            the extracted file object's property 'name' is set to
            the filename of the archive.
        """
        instance = cls()
        instance.read(filename, fileobj=fileobj)
        return instance

    def _get_calibrated_channels_uniform_shape(self):
        """Prepare the channels as input for gac_io.save_gac."""
        channels = self.get_calibrated_channels()
        assert channels.shape[-1] == 6
        return channels

    def save(
        self,
        start_line,
        end_line,
        output_file_prefix="PyGAC",
        output_dir="./",
        avhrr_dir=None,
        qual_dir=None,
        sunsatangles_dir=None,
    ):
        """Convert the Reader instance content into hdf5 files."""
        avhrr_dir = avhrr_dir or output_dir
        qual_dir = qual_dir or output_dir
        sunsatangles_dir = sunsatangles_dir or output_dir
        self.get_lonlat()
        channels = self._get_calibrated_channels_uniform_shape()
        sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = self.get_angles()
        qual_flags = self.get_qual_flags()
        if np.all(self.mask):
            print("ERROR: All data is masked out. Stop processing")
            raise ValueError("All data is masked out.")
        gac_io.save_gac(
            self.spacecraft_name,
            self._times_as_np_datetime64,
            self.lats,
            self.lons,
            channels[:, :, 0],
            channels[:, :, 1],
            channels[:, :, 2],
            channels[:, :, 3],
            channels[:, :, 4],
            channels[:, :, 5],
            sun_zen,
            sat_zen,
            sun_azi,
            sat_azi,
            rel_azi,
            qual_flags,
            start_line,
            end_line,
            self.filename,
            self.meta_data,
            output_file_prefix,
            avhrr_dir,
            qual_dir,
            sunsatangles_dir,
        )

    @abstractmethod
    def get_header_timestamp(self):  # pragma: no cover
        """Read start timestamp from the header.

        Returns:
            datetime.datetime: Start timestamp

        """
        raise NotImplementedError

    def get_counts(self):
        """Get the counts.

        Returns:
            np.array: The counts, with channel 3a and 3b split if necessary.

        """
        packed_data = self.scans["sensor_data"]
        counts = np.zeros((len(self.scans), self.scan_width * 5))
        counts_nb = (self.scan_width * 5) // 3
        remainder = (self.scan_width * 5) % 3
        if remainder == 0:
            nb1 = nb2 = nb3 = counts_nb
        elif remainder == 1:
            nb1 = counts_nb + 1
            nb2 = nb3 = counts_nb
        elif remainder == 2:
            nb1 = nb2 = counts_nb + 1
            nb3 = counts_nb
        counts[:, 0::3] = ((packed_data >> 20) & 1023)[:, :nb1]
        counts[:, 1::3] = ((packed_data >> 10) & 1023)[:, :nb2]
        counts[:, 2::3] = (packed_data & 1023)[:, :nb3]
        counts = counts.reshape((-1, self.scan_width, 5))
        try:
            switch = self.get_ch3_switch()
        except AttributeError:
            return counts
        else:
            channels = np.zeros((len(self.scans), self.scan_width, 6), dtype=counts.dtype)
            channels[:, :, :2] = counts[:, :, :2]
            channels[:, :, -2:] = counts[:, :, -2:]
            channels[:, :, 2][switch == 1] = counts[:, :, 2][switch == 1]
            channels[:, :, 3][switch == 0] = counts[:, :, 2][switch == 0]
            return channels

    @abstractmethod
    def _get_times_from_file(self):  # pragma: no cover
        """Specify how to read scanline timestamps from data.

        Returns:
            int: year
            int: day of year
            int: milliseconds since 00:00

        """
        raise NotImplementedError

    def get_times(self):
        """Read scanline timestamps and try to correct invalid values.

        Returns:
            UTC timestamps

        """
        if self._times_as_np_datetime64 is None:
            # Read timestamps
            year, jday, msec = self._get_times_from_file()
            # Correct invalid values
            year, jday, msec = self.correct_times_median(year=year, jday=jday, msec=msec)
            self._times_as_np_datetime64 = self.to_datetime64(year=year, jday=jday, msec=msec)
            try:
                self._times_as_np_datetime64 = self.correct_times_thresh()
            except TimestampMismatch as err:
                LOG.error(str(err))

        return self._times_as_np_datetime64

    @staticmethod
    def to_datetime64(year, jday, msec):
        """Convert timestamps to numpy.datetime64.

        Args:
            year: Year
            jday: Day of the year (1-based)
            msec: Milliseconds since 00:00

        Returns:
            numpy.datetime64: Converted timestamps

        """
        return (
            year.astype(str).astype("datetime64[Y]")
            + (jday - 1).astype("timedelta64[D]")
            + msec.astype("timedelta64[ms]")
        )

    @staticmethod
    def to_datetime(datetime64):
        """Convert numpy.datetime64 to datetime.datetime.

        Args:
            datetime64 (numpy.datetime64): Numpy timestamp to be converted.

        Returns:
            datetime.datetime: Converted timestamp

        """
        return datetime64.astype(datetime.datetime)

    def lineno2msec(self, scan_line_number):
        """Compute ideal scanline timestamp based on the scanline number.

        Assumes a constant scanning frequency.

        Args:
            scan_line_number: Specifies the scanline number (1-based)

        Returns:
            Corresponding timestamps in milliseconds since 1970-01-01 00:00,
            i.e. the first scanline has timestamp 0.

        """
        return (scan_line_number - 1) / self.scan_freq

    def get_sun_earth_distance_correction(self):
        """Get the julian day and the sun-earth distance correction."""
        self.get_times()
        # Shouldn't we use the start time from the header if available?
        start_time = self._times_as_np_datetime64[0]
        jday = (start_time - start_time.astype("datetime64[Y]")).astype("timedelta64[D]").astype(int) + 1
        return calculate_sun_earth_distance_correction(jday)

    def update_meta_data(self):
        """Add some meta data to the meta_data dicitonary."""
        meta_data = self.meta_data
        self._update_meta_data_object(meta_data)
        if "gac_header" not in meta_data:
            meta_data["gac_header"] = self.head.copy()

    def _update_meta_data_object(self, meta_data):
        if "sun_earth_distance_correction_factor" not in meta_data:
            meta_data["sun_earth_distance_correction_factor"] = self.get_sun_earth_distance_correction()
        if "midnight_scanline" not in meta_data:
            meta_data["midnight_scanline"] = self.get_midnight_scanline()
        if "missing_scanlines" not in meta_data:
            meta_data["missing_scanlines"] = self.get_miss_lines()

    def read_as_dataset(self, file_to_read):
        self.read(file_to_read)
        return self.create_counts_dataset()

    def create_counts_dataset(self):
        """Create output xarray dataset containing counts and information relevant to calibration.

        Contents of the dataset:
            channels: Earth scene counts
            prt_counts: Counts from PRTs on ICT
            ict_counts: Counts from observed ICT
            space_counts: Counts from space observation
            bad_space_scans: Scanlines with suspect space view information
            noise: Noise estimates in counts
            ict_noise: ICT noise estimates in counts
            Longitude/latitude: pixel position
        """
        head = dict(zip(self.head.dtype.names, self.head.item(), strict=False))
        scans = self.scans

        counts = self.get_counts()
        line_numbers = scans["scan_line_number"]
        scan_size = counts.shape[1]
        columns = np.arange(scan_size)
        longitudes, latitudes = self._get_lonlat_dataarrays(line_numbers, columns)

        times = self.get_times()

        if counts.shape[-1] == 5:
            channel_names = ["1", "2", "3", "4", "5"]
            ir_channel_names = ["3", "4", "5"]
            vis_channel_names = ["1", "2","3"]         # added
        else:
            channel_names = ["1", "2", "3a", "3b", "4", "5"]
            ir_channel_names = ["3b", "4", "5"]
            vis_channel_names = ["1", "2", "3a"]

        channels = xr.DataArray(
            counts,
            dims=["scan_line_index", "columns", "channel_name"],
            coords=dict(
                scan_line_index=line_numbers,
                channel_name=channel_names,
                columns=columns,
                times=("scan_line_index", times),
            ),
        )

        #
        # Added total (10 per scan line) arrays
        # J.Mittaz University of Reading
        #
        prt, ict, space, total_ict, total_space \
            = self._get_telemetry_dataarrays(line_numbers, ir_channel_names)

        #
        # Added vis_space and total_vis_space to outputs
        # N.Yaghnam, NPL
        #
        vis_space, total_vis_space \
            = self._get_vis_telemetry_dataarrays(line_numbers, vis_channel_names)
        #
        # Get angles - need sun_zen for solar contamination in
        # calibration routines - Added by J.Mittaz / UoR
        #
        sat_azi, sat_zen, sun_azi, sun_zen, rel_azi = self.get_angles()
        sun_zen = xr.DataArray(sun_zen,
                                 dims=["scan_line_index", "columns"],
                                   coords=dict(scan_line_index=line_numbers,columns=columns))


        if self.interpolate_coords:
            channels = channels.assign_coords(
                longitude=(("scan_line_index", "columns"), longitudes.reindex_like(channels).data),
                latitude=(("scan_line_index", "columns"), latitudes.reindex_like(channels).data),
            )

        #
        # Added space_views and noise entries
        #
        # Added vis_space and total_vis_space - N.Yagnam / NPL
        #
        ds = xr.Dataset(dict(channels=channels, prt_counts=prt,
                             ict_counts=ict, space_counts=space,
                             total_ict_counts=total_ict,
                             total_space_counts=total_space,
                             vis_space_counts=vis_space,
                             total_vis_space_counts=total_vis_space,
                             longitude=longitudes, latitude=latitudes,
                             sun_zen=sun_zen),
                             attrs=head)


        ds.attrs["spacecraft_name"] = self.spacecraft_name
        ds.attrs["max_scan_angle"] = 55.25 if self.spacecraft_name == "noaa16" else 55.37
        self._update_meta_data_object(ds.attrs)
        return ds

    def _get_lonlat_dataarrays(self, line_numbers, columns):
        lons, lats = self.get_lonlat()
        lon_lat_columns = columns if self.interpolate_coords is True else self.lonlat_sample_points

        if self.interpolate_coords:
            longitudes = xr.DataArray(
                lons,
                dims=["scan_line_index", "columns"],
                coords={"scan_line_index": line_numbers, "columns": lon_lat_columns},
            )
            latitudes = xr.DataArray(
                lats,
                dims=["scan_line_index", "columns"],
                coords={"scan_line_index": line_numbers, "columns": lon_lat_columns},
            )

        else:
            longitudes = xr.DataArray(
                lons,
                dims=["scan_line_index", "subsampled_columns"],
                coords={"scan_line_index": line_numbers, "subsampled_columns": lon_lat_columns},
            )
            latitudes = xr.DataArray(
                lats,
                dims=["scan_line_index", "subsampled_columns"],
                coords={"scan_line_index": line_numbers, "subsampled_columns": lon_lat_columns},
            )

        return longitudes, latitudes

    def _get_telemetry_dataarrays(self, line_numbers, ir_channel_names):
        """Get data from lower telemetry including bad_scans and noise added
        by J.Mittaz UoR"""
        prt, ict, space, total_space, total_ict = self.get_telemetry()

        prt = xr.DataArray(prt, dims=["scan_line_index"], coords=dict(scan_line_index=line_numbers))
        ict = xr.DataArray(ict, dims=["scan_line_index", "ir_channel_name"],
                           coords=dict(ir_channel_name=ir_channel_names, scan_line_index=line_numbers))
        space = xr.DataArray(space, dims=["scan_line_index", "ir_channel_name"],
                             coords=dict(ir_channel_name=ir_channel_names, scan_line_index=line_numbers))
        #
        # New entries for calibration uncertainty work
        #
        pixel_index = np.arange(10,dtype=np.int8)
        total_ict = xr.DataArray(total_ict,
                                 dims=["scan_line_index", "pixel_index", "ir_channel_name"],
                                 coords=dict(ir_channel_name=ir_channel_names, scan_line_index=line_numbers))
        total_space = xr.DataArray(total_space,
                                   dims=["scan_line_index", "pixel_index", "ir_channel_name"],
                                   coords=dict(ir_channel_name=ir_channel_names, scan_line_index=line_numbers,
                                               pixel_index=pixel_index))

        return prt, ict, space, total_ict, total_space

    def _get_vis_telemetry_dataarrays(self, line_numbers, vis_channel_names):
        """Get data from lower telemetry for the visible channels. Added by N.Yaghnam, NPL"""
        vis_space, total_vis_space = self.get_vis_telemetry()

        pixel_index = np.arange(10, dtype=np.int8)
        vis_space = xr.DataArray(vis_space, dims=["scan_line_index", "vis_channel_name"],
                             coords=dict(vis_channel_name=vis_channel_names, scan_line_index=line_numbers))
        total_vis_space = xr.DataArray(total_vis_space,
                                   dims=["scan_line_index", "pixel_index", "vis_channel_name"],
                                   coords=dict(vis_channel_name=vis_channel_names, scan_line_index=line_numbers,
                                               pixel_index=pixel_index))

        return vis_space, total_vis_space

    def get_calibrated_channels(self):
        """Calibrate and return the channels."""
        calibrated_ds = self.calibrated_dataset

        channels = calibrated_ds["channels"].data
        with suppress(KeyError):
            self.meta_data["calib_coeffs_version"] = calibrated_ds.attrs["calib_coeffs_version"]

        return channels

    @cached_property
    def calibrated_dataset(self):
        """Return the calibrated dataset."""
        return self.get_calibrated_dataset()

    def get_calibrated_dataset(self):
        """Create and calibrate the dataset for the pass."""
        ds = self.create_counts_dataset()
        #
        # Make sure earth counts are kept for uncertainty calculation
        #
        counts = xr.DataArray(name="counts",
                              data=np.copy(ds["channels"].data),
                              dims=ds["channels"].dims,
                              coords=ds["channels"].coords)

        # calibration = {"1": "mitram", "2": "mitram", "4": {"method":"noaa", "coeff_file": "myfile.json"}}

        calibration_entrypoints = entry_points(group="pygac.calibration")
        calibration_function = calibration_entrypoints[self.calibration_method].load()
        calibrated_ds = calibration_function(ds, **self.calibration_parameters)
        calibrated_ds["counts"] = counts

        # Mask out corrupt values
        mask = xr.DataArray(self.mask == False, dims=["scan_line_index"])  # noqa
        calibrated_ds = calibrated_ds.where(mask)

        # Apply KLM/POD specific postprocessing
        self.postproc(calibrated_ds)

        calibrated_ds.attrs["tle"] = self.get_tle_lines()

        # Mask pixels affected by scan motor issue
        if self.is_tsm_affected():
            LOG.info("Correcting for temporary scan motor issue")
            self.mask_tsm_pixels(calibrated_ds)
        if self.reference_image:
            try:
                self._georeference_data(calibrated_ds)
                calibrated_ds.attrs["georeferenced"] = True
            except:  # noqa
                LOG.exception("Could not georeference!")
                calibrated_ds.attrs["georeferenced"] = False
        if self.compute_uncertainties:
            try:
                from pygac.calibration.uncertainty import uncertainty
                ucs = uncertainty(calibrated_ds, self.mask)

                calibrated_ds["random_uncertainty"] = ucs["random"]
                calibrated_ds["systematic_uncertainty"] = ucs["systematic"]
                calibrated_ds["channel_covariance_ratio"] = ucs["chan_covar_ratio"]
                calibrated_ds["uncertainty_flags"] = ucs["uncert_flags"]

                calibrated_ds.attrs["uncertainties_computed"] = True
            except:  # noqa
                LOG.exception("Could not compute uncertainties!")
                calibrated_ds.attrs["uncertainties_computed"] = False
        return calibrated_ds

    def _georeference_data(self, calibrated_ds):
        if not self.adjust_clock_drift:
            self._correct_time_offset(calibrated_ds)

        from georeferencer.georeferencer import get_swath_displacement

        _, sat_zen, _, sun_zen, _ = self.get_angles()
        time_diff, (roll, pitch, yaw), (odistances, mdistances) = get_swath_displacement(
            calibrated_ds, sun_zen, sat_zen, self.reference_image, self.dem
        )

        if (mdist := np.median(mdistances)) > 5000:
            raise RuntimeError("Displacement minimization did not produce convincing improvements")
        calibrated_ds.attrs["median_gcp_distance"] = mdist

        self._rpy = roll, pitch, yaw
        time_diff = np.timedelta64(int(time_diff * 1e9), "ns")
        lons, lats = self._compute_lonlats(time_offset=time_diff)
        self._times_as_np_datetime64 += time_diff
        calibrated_ds["longitude"].data = lons
        calibrated_ds["latitude"].data = lats
        if self.dem:
            from georeferencer.georeferencer import orthocorrection

            calibrated_ds = orthocorrection(calibrated_ds, sat_zen, self.dem)
        calibrated_ds["times"].data = self._times_as_np_datetime64

    def _correct_time_offset(self, calibrated_ds) -> None:
        thinned_lons, thinned_lats = self._get_lonlat_from_file()
        gcps = np.array(
            (
                (0, self.lonlat_sample_points[0]),
                (0, self.lonlat_sample_points[-1]),
                (len(self.scans) - 1, self.lonlat_sample_points[0]),
                (len(self.scans) - 1, self.lonlat_sample_points[-1]),
            )
        )
        ref_lons = (thinned_lons[0, 0], thinned_lons[0, -1], thinned_lons[-1, 0], thinned_lons[-1, -1])
        ref_lats = (thinned_lats[0, 0], thinned_lats[0, -1], thinned_lats[-1, 0], thinned_lats[-1, -1])
        from pyorbital.geoloc_avhrr import estimate_time_offset

        time_diff, _ = estimate_time_offset(
            gcps,
            ref_lons,
            ref_lats,
            calibrated_ds["times"][0].values,
            calibrated_ds.attrs["tle"],
            calibrated_ds.attrs["max_scan_angle"],
        )
        time_diff = np.timedelta64(int(time_diff * 1e9), "ns")
        lons, lats = self._compute_lonlats(time_offset=time_diff)
        self._times_as_np_datetime64 += time_diff

        calibrated_ds["longitude"].data = lons
        calibrated_ds["latitude"].data = lats
        calibrated_ds["times"].data = self._times_as_np_datetime64

    @abstractmethod
    def get_telemetry(self):  # pragma: no cover
        """KLM/POD specific readout of telemetry."""
        raise NotImplementedError

    def lonlat_interpolator(self, lons, lats, cols_subset=None, cols_full=None):
        """Interpolate from lat-lon tie-points to pixel locations

        Args:
            lons: Longitude tie-points
            lats: Latitude tie-points

        Returns:
            pixel_longitudes, pixel_latitudes
        """
        if cols_subset is None and cols_full is None:
            cols_subset = self.lonlat_sample_points
            cols_full = np.arange(self.scan_width)
        rows = np.arange(len(lats))

        along_track_order = 1
        cross_track_order = 3

        satint = gtp.SatelliteInterpolator(
            (lons, lats), (rows, cols_subset), (rows, cols_full), along_track_order, cross_track_order
        )

        return satint.interpolate()

    def get_lonlat(self):
        """Compute lat/lon coordinates.

        TODO: Switch to faster interpolator?
        """
        if self.lons is None or self.lats is None:
            return self._compute_lonlats(mask_scanlines=not self.reference_image)

        return self.lons, self.lats

    def _compute_lonlats(self, time_offset=None, mask_scanlines=True):
        if not self.compute_lonlats_from_tles:
            self.lons, self.lats = self._get_lonlat_from_file()
            # Adjust clock drift
            if self.adjust_clock_drift and not self.clock_drift_correction_applied:
                self._adjust_clock_drift()
                self.clock_drift_correction_applied = True
        else:
            self.get_times()
            if self.adjust_clock_drift and not self.clock_drift_correction_applied:
                try:
                    offsets = self.compute_clock_offsets()
                    self._times_as_np_datetime64 -= (offsets * 1000).astype("timedelta64[ms]")
                    self.clock_drift_correction_applied = True
                    LOG.debug("Applied clock drift correction")
                except AttributeError:
                    LOG.debug("No clock drift correction applied")
            if time_offset:
                new_times = self._times_as_np_datetime64 + time_offset
                self.lons, self.lats = self._compute_lonlat_from_tles(new_times.astype("datetime64[ms]"))
            else:
                self.lons, self.lats = self._compute_lonlat_from_tles(self._times_as_np_datetime64)
        self.update_meta_data()

        # Interpolate from every eighth pixel to all pixels.
        if self.interpolate_coords:
            self.lons, self.lats = self.lonlat_interpolator(self.lons, self.lats)

        # Mask out corrupt scanlines
        if mask_scanlines:
            self.lons[self.mask] = np.nan
            self.lats[self.mask] = np.nan

        # Mask values outside the valid range
        self.lats[np.fabs(self.lats) > 90.0] = np.nan
        self.lons[np.fabs(self.lons) > 180.0] = np.nan

        return self.lons, self.lats

    @abstractmethod
    def _get_lonlat_from_file(self):  # pragma: no cover
        """KLM/POD specific readout of lat/lon coordinates."""
        raise NotImplementedError

    @property
    def mask(self):
        """Mask for corrupt scanlines."""
        if self._mask is None:
            self._mask = self._get_corrupt_mask()
        return self._mask

    @property
    @abstractmethod
    def QFlag(self):  # pragma: no cover
        """KLM/POD specific quality indicators."""
        raise NotImplementedError

    @property
    @abstractmethod
    def _quality_indicators_key(self):  # pragma: no cover
        raise NotImplementedError

    def _get_corrupt_mask(self, flags=None):
        """Readout of corrupt scanline mask.

        Args:
            flags (QFlag.flag): An "ORed" bitmask that defines corrupt values.
                                Defauts to (QFlag.FATAL_FLAG | QFlag.CALIBRATION
                                | QFlag.NO_EARTH_LOCATION)

        Note:
            The Quality flags mapping (QFlag) is KLM/POD specific.
        """
        QFlag = self.QFlag
        if flags is None:
            flags = QFlag.FATAL_FLAG | QFlag.CALIBRATION | QFlag.NO_EARTH_LOCATION
        return (self.scans[self._quality_indicators_key] & int(flags)).astype(bool)

    def get_qual_flags(self):
        """Read quality flags."""
        number_of_scans = self.scans["telemetry"].shape[0]
        qual_flags = np.zeros((int(number_of_scans), 7))
        qual_flags[:, 0] = self.scans["scan_line_number"]
        qual_flags[:, 1] = self._get_corrupt_mask(flags=self.QFlag.FATAL_FLAG)
        qual_flags[:, 2] = self._get_corrupt_mask(flags=self.QFlag.CALIBRATION)
        qual_flags[:, 3] = self._get_corrupt_mask(flags=self.QFlag.NO_EARTH_LOCATION)
        qual_flags[:, 4] = self._get_corrupt_mask(flags=self.QFlag.CH_3_CONTAMINATION)
        qual_flags[:, 5] = self._get_corrupt_mask(flags=self.QFlag.CH_4_CONTAMINATION)
        qual_flags[:, 6] = self._get_corrupt_mask(flags=self.QFlag.CH_5_CONTAMINATION)
        return qual_flags

    @abstractmethod
    def postproc(self, ds):  # pragma: no cover
        """Apply KLM/POD specific postprocessing."""
        raise NotImplementedError

    @abstractmethod
    def _adjust_clock_drift(self):  # pragma: no cover
        """Adjust clock drift."""
        raise NotImplementedError

    @staticmethod
    def tle2datetime64(times):
        """Convert TLE timestamps to numpy.datetime64.

        Args:
           times (float): TLE timestamps as %y%j.1234, e.g. 18001.25
        """
        # Convert %y%j.12345 to %Y%j.12345 (valid for 1950-2049)
        times = np.where(times > 50000, times + 1900000, times + 2000000)

        # Convert float to datetime64
        doys = (times % 1000).astype("int") - 1
        years = (times // 1000).astype("int")
        msecs = np.rint(24 * 3600 * 1000 * (times % 1))
        times64 = (years - 1970).astype("datetime64[Y]").astype("datetime64[ms]")
        times64 += doys.astype("timedelta64[D]")
        times64 += msecs.astype("timedelta64[ms]")

        return times64

    def get_tle_file(self):
        """Find TLE file for the current satellite."""
        tle_dir, tle_name = self.tle_dir, self.tle_name
        if tle_dir is None:
            raise RuntimeError("TLE directory not specified!")
        if tle_name is None:
            raise RuntimeError("TLE name not specified!")
        values = {
            "satname": self.spacecraft_name,
        }
        tle_filename = os.path.join(tle_dir, tle_name % values)
        LOG.info("TLE filename = " + str(tle_filename))

        return tle_filename

    def read_tle_file(self, tle_filename):
        """Read TLE file."""
        with open(tle_filename, "r") as fp_:
            return [line for line in fp_.readlines() if line.strip()]

    def get_tle_lines(self):
        """Find closest two line elements (TLEs) for the current orbit.

        Raises:
            IndexError, if the closest TLE is more than :meth:`pygac.GACReader.tle_thresh` days apart

        """
        if self.tle_lines is not None:
            return self.tle_lines
        self.get_times()
        tle_data = self.read_tle_file(self.get_tle_file())
        sdate = self._times_as_np_datetime64[0]
        dates = self.tle2datetime64(np.array([float(line[18:32]) for line in tle_data[::2]]))

        # Find index "iindex" such that dates[iindex-1] < sdate <= dates[iindex]
        # Notes:
        #     1. If sdate < dates[0] then iindex = 0
        #     2. If sdate > dates[-1] then iindex = len(dates), beyond the right boundary!
        iindex = np.searchsorted(dates, sdate)

        if iindex in (0, len(dates)):
            if iindex == len(dates):
                # Reset index if beyond the right boundary (see note 2. above)
                iindex -= 1
        elif abs(sdate - dates[iindex - 1]) < abs(sdate - dates[iindex]):
            # Choose the closest of the two surrounding dates
            iindex -= 1

        # Make sure the TLE we found is within the threshold
        delta_days = abs(sdate - dates[iindex]) / np.timedelta64(1, "D")
        if delta_days > self.tle_thresh:
            raise NoTLEData(
                "Can't find tle data for %s within +/- %d days around %s"
                % (self.spacecraft_name, self.tle_thresh, sdate)
            )

        if delta_days > 3:
            LOG.warning("Found TLE data for %s that is %f days apart", sdate, delta_days)
        else:
            LOG.debug("Found TLE data for %s that is %f days apart", sdate, delta_days)

        # Select TLE data
        tle1 = tle_data[iindex * 2]
        tle2 = tle_data[iindex * 2 + 1]
        self.tle_lines = tle1, tle2
        return tle1, tle2

    def get_sat_angles(self):
        """Get satellite angles.

        Returns:
            Azimuth, elevation (degrees)
        """
        try:
            return self._get_sat_angles_with_tle()
        except NoTLEData:
            LOG.warning("No TLE data available. Falling back to approximate calculation of satellite angles.")
            return self._get_sat_angles_without_tle()

    def _get_sat_angles_with_tle(self):
        tle1, tle2 = self.get_tle_lines()
        orb = Orbital(self.spacecrafts_orbital[self.spacecraft_id], line1=tle1, line2=tle2)
        sat_azi, sat_elev = orb.get_observer_look(self._times_as_np_datetime64[:, np.newaxis], self.lons, self.lats, 0)
        return sat_azi, sat_elev

    def _get_sat_angles_without_tle(self):
        """Get satellite angles using lat/lon from data to approximate satellite postition instead of TLE."""
        from pyorbital.orbital import get_observer_look as get_observer_look_no_tle

        LOG.warning("Approximating satellite height to 850km (TIROS-N OSCAR)!")
        sat_alt = 850.0  # km  TIROS-N OSCAR
        mid_column = int(0.5 * self.lons.shape[1])
        sat_azi, sat_elev = get_observer_look_no_tle(
            self.lons[:, mid_column][:, np.newaxis],
            self.lats[:, mid_column][:, np.newaxis],  # approximate satellite position
            sat_alt,  # approximate satellite altitude
            self._times_as_np_datetime64[:, np.newaxis],
            self.lons,
            self.lats,
            0,
        )
        # Sometimes (pyorbital <= 1.6.1) the get_observer_look_not_tle returns nodata instead of 90.
        # Problem solved with https://github.com/pytroll/pyorbital/pull/77
        if Version(pyorbital.__version__) <= Version("1.6.1"):
            sat_elev[:, mid_column] = 90
        return sat_azi, sat_elev

    def get_angles(self):
        """Get azimuth and zenith angles.

        Azimuth angle definition is the same as in pyorbital, but with
        different units (degrees not radians for sun azimuth angles)
        and different ranges.

        Returns:
            sat_azi: satellite azimuth angle degree clockwise from north in
            range ]-180, 180]

            sat_zenith: satellite zenith angles in degrees in range [0,90]

            sun_azi: sun azimuth angle degree clockwise from north in range
            ]-180, 180]

            sun_zenith: sun zenith angles in degrees in range [0,90]

            rel_azi: absolute azimuth angle difference in degrees between sun
            and sensor in range [0, 180]

        """
        self.get_times()
        self.get_lonlat()
        times = self._times_as_np_datetime64
        sat_azi, sat_elev = self.get_sat_angles()

        sat_zenith = 90 - sat_elev
        sun_zenith = astronomy.sun_zenith_angle(times[:, np.newaxis], self.lons, self.lats)

        alt, sun_azi = astronomy.get_alt_az(times[:, np.newaxis], self.lons, self.lats)
        del alt
        sun_azi = np.rad2deg(sun_azi)
        rel_azi = get_absolute_azimuth_angle_diff(sun_azi, sat_azi)

        # Scale angles range to half open interval ]-180, 180]
        sat_azi = centered_modulus(sat_azi, 360.0)
        sun_azi = centered_modulus(sun_azi, 360.0)

        # Mask corrupt scanlines
        for arr in (sat_azi, sat_zenith, sun_azi, sun_zenith, rel_azi):
            arr[self.mask] = np.nan

        return sat_azi, sat_zenith, sun_azi, sun_zenith, rel_azi

    def correct_times_median(self, year, jday, msec):
        """Replace invalid timestamps with statistical estimates (using median).

        Assumes that the majority of timestamps is ok.

        Args:
            year: Year
            jday: Day of the year
            msec: Milliseconds since 00:00

        Returns:
            Corrected year
            Corrected day of the year
            Corrected milliseconds

        """
        # Estimate ideal timestamps based on the scanline number. Still without
        # offset, e.g. the first scanline has timestamp 1970-01-01 00:00
        msec_lineno = self.lineno2msec(self.scans["scan_line_number"])

        jday = np.where(np.logical_or(jday < 1, jday > 366), np.median(jday), jday)
        if_wrong_jday = np.ediff1d(jday, to_begin=jday.dtype.type(0))
        jday = np.where(if_wrong_jday < 0, max(jday), jday)

        if_wrong_msec = np.where(msec < 1)
        if_wrong_msec = if_wrong_msec[0]
        if len(if_wrong_msec) > 0:
            if if_wrong_msec[0] != 0:
                msec = msec[0] + msec_lineno
            else:
                msec0 = np.median(msec - msec_lineno)
                msec = msec0 + msec_lineno

        if_wrong_msec = np.ediff1d(msec, to_begin=msec.dtype.type(0))
        msec = np.where(
            np.logical_and(np.logical_or(if_wrong_msec < -1000, if_wrong_msec > 1000), if_wrong_jday != 1),
            msec[0] + msec_lineno,
            msec,
        )

        # checking if year value is out of valid range
        if_wrong_year = np.where(np.logical_or(year < 1978, year > datetime.datetime.now().year))
        if_wrong_year = if_wrong_year[0]
        if len(if_wrong_year) > 0:
            # if the first scanline has valid time stamp
            if if_wrong_year[0] != 0:
                year = year[0]
                jday = jday[0]
                msec = msec[0] + msec_lineno
            # Otherwise use median time stamp
            else:
                year = np.median(year)
                jday = np.median(jday)
                msec0 = np.median(msec - msec_lineno)
                msec = msec0 + msec_lineno

        return year.astype(int), jday.astype(int), msec

    def correct_scan_line_numbers(self):
        """Remove scanlines with corrupted scanline numbers.

        This includes:
            - Scanline numbers outside the valid range
            - Scanline numbers deviating more than a certain threshold from the
              ideal case (1,2,3,...N)

        Example files having corrupt scanline numbers:
            - NSS.GHRR.NJ.D96144.S2000.E2148.B0720102.GC
            - NSS.GHRR.NJ.D96064.S0043.E0236.B0606162.WI
            - NSS.GHRR.NJ.D99286.S1818.E2001.B2466869.WI

        Returns:
            Intermediate and final results (for plotting purpose)

        """
        along_track = np.arange(1, len(self.scans["scan_line_number"]) + 1)
        results = {"along_track": along_track, "n_orig": self.scans["scan_line_number"].copy()}

        # Remove scanlines whose scanline number is outside the valid range
        within_range = np.logical_and(
            self.scans["scan_line_number"] < self.max_scanlines, self.scans["scan_line_number"] >= 0
        )
        self.scans = self.scans[within_range]

        # Remove scanlines deviating more than a certain threshold from the
        # ideal case (1,2,3,...N).
        ideal = np.arange(1, len(self.scans["scan_line_number"]) + 1)

        # ... Estimate possible offset (in case some scanlines are missing in
        # the beginning of the scan)
        offsets = self.scans["scan_line_number"] - ideal
        med_offset = np.median(offsets)

        # ... Compute difference to ideal case (1,2,3,...N) + offset
        diffs = np.abs(self.scans["scan_line_number"] - (ideal + med_offset))

        # ... Remove those scanlines whose difference is larger than a certain
        # threshold. For the threshold computation we only regard nonzero
        # differences.
        nz_diffs = diffs[diffs > 0]
        if len(nz_diffs) < 50:
            # Not enough differences for reliable statistics. Use fixed
            # threshold.
            thresh = 500
        else:
            mean_nz_diffs = np.mean(nz_diffs)
            std_nz_diffs = np.std(nz_diffs)
            med_nz_diffs = np.median(nz_diffs)
            mad_nz_diffs = np.median(np.abs(nz_diffs - med_nz_diffs))
            if mean_nz_diffs / float(med_nz_diffs) < 3:
                # Relatively small variation, keep (almost) everything
                thresh = mean_nz_diffs + 3 * std_nz_diffs
            else:
                # Large variation, filter more aggressively. Use median and
                # median absolute deviation (MAD) as they are less sensitive to
                # outliers. However, allow differences < 500 scanlines as they
                # occur quite often.
                thresh = max(500, med_nz_diffs + 3 * mad_nz_diffs)
        self.scans = self.scans[diffs <= thresh]

        LOG.debug("Removed %s scanline(s) with corrupt scanline numbers", str(len(along_track) - len(self.scans)))

        results.update(
            {
                "n_corr": self.scans["scan_line_number"],
                "within_range": within_range,
                "diffs": diffs,
                "thresh": thresh,
                "nz_diffs": nz_diffs,
            }
        )
        return results

    def correct_times_thresh(
        self, max_diff_from_t0_head=6 * 60 * 1000, min_frac_near_t0_head=0.01, max_diff_from_ideal_t=10 * 1000
    ):
        """Correct corrupted timestamps using a threshold approach.

        The threshold approach is based on the scanline number and the header
        timestamp. It also works if the majority of scanlines has corrupted
        timestamps.

        The header timestamp is used as a guideline to estimate the offset
        between timestamps computed from the scanline number and the actual
        scanline timestamps in the data. If header timestamp and scanline
        timestamps do not match, no correction is applied.

        Once the offset has been estimated, one can calculate the ideal
        timestamps based on the scanline number. Timestamps deviating more than
        a certain threshold from the ideal timestamps are replaced by
        the ideal timestamps.

        Example files having corrupt timestamps:
            - NSS.GHRR.NA.D81193.S2329.E0116.B1061214.WI
            - NSS.GHRR.NL.D01035.S2342.E0135.B0192627.WI

        Args:
            max_diff_from_t0_head (int): Threshold for offset estimation: A
                scanline timestamp matches the header timestamp t0_head if it is
                within the interval

                    [t0_head - max_diff_from_t0_head,
                    t0_head + max_diff_from_t0_head]

                around the header timestamp.
            min_frac_near_t0_head (float): Specifies the minimum fraction of
                scanline timestamps matching the header timestamp required for
                applying the correction.
            max_diff_from_ideal_t (float): Threshold for timestamp correction:
                If a scanline timestamp deviates more than max_diff_from_ideal_t
                from the ideal timestamp, it is regarded as corrupt and will be
                replaced with the ideal timestamp.

        Returns:
            Intermediate and final results (for plotting purpose)

        """
        results = {}
        dt64_msec = ">M8[ms]"

        # Check whether scanline number increases monotonically
        nums = self.scans["scan_line_number"]
        results.update({"t": self._times_as_np_datetime64.copy(), "n": nums})
        if np.any(np.diff(nums) < 0):
            LOG.error("Cannot perform timestamp correction. Scanline number does not increase monotonically.")
            results["fail_reason"] = "Scanline number jumps backwards"
            return results

        # Convert time to milliseconds since 1970-01-01
        t = self._times_as_np_datetime64.astype("i8")
        try:
            t0_head = np.array([self.get_header_timestamp().isoformat()], dtype="datetime64[ms]").astype("i8")[0]
        except ValueError as err:
            LOG.error("Cannot perform timestamp correction: %s", err)
            return

        # Compute ideal timestamps based on the scanline number. Still
        # without offset, i.e. scanline 0 has timestamp 1970-01-01 00:00
        tn = self.lineno2msec(nums)

        # Try to determine the timestamp t0 of the first scanline. Since none
        # of the actual timestamps is trustworthy, use the header timestamp
        # as a guideline. However, the header timestamp may also be corrupted,
        # so we only apply corrections if there is a minimum fraction of
        # scanlines whose timestamps match the header timestamp.
        #
        # 1) Compute offsets between actual timestamps and idealized timestamps
        offsets = t - tn

        # 2) If the offsets of a certain minimum fraction of scanlines are
        #    within a certain interval around the header timestamp, estimate
        #    t0 by calculating the median offset among these timestamps. If not,
        #    we do not have reliable information and cannot proceed.
        near_t0_head = np.where(np.fabs(offsets - t0_head) <= max_diff_from_t0_head)[0]
        results.update({"offsets": offsets, "t0_head": t0_head, "max_diff_from_t0_head": max_diff_from_t0_head})
        if near_t0_head.size / float(nums.size) >= min_frac_near_t0_head:
            t0 = np.median(offsets[near_t0_head])
        else:
            results["fail_reason"] = "Timestamp mismatch"
            raise TimestampMismatch("Timestamp mismatch. Cannot perform timestamp correction.")

        # Add estimated offset to the ideal timestamps
        tn += t0

        # Replace timestamps deviating more than a certain threshold from the
        # ideal timestamp with the ideal timestamp.
        corrupt_lines = np.where(np.fabs(t - tn) > max_diff_from_ideal_t)
        self._times_as_np_datetime64[corrupt_lines] = tn[corrupt_lines].astype(dt64_msec)
        LOG.debug("Corrected %s timestamp(s)", str(len(corrupt_lines[0])))

        return self._times_as_np_datetime64

    @property
    @abstractmethod
    def tsm_affected_intervals(self):  # pragma: no cover
        """Specify time intervals being affected by the scan motor problem.

        Returns:
            dict: Affected time intervals. A dictionary containing a list of
                (start, end) tuples for each affected platform. Both start and
                end must be datetime.datetime objects.

        """
        raise NotImplementedError

    def is_tsm_affected(self):
        """Determine whether this orbit is affected by the scan motor problem.

        Returns:
            bool: True if the orbit is affected, False otherwise.

        """
        self.get_times()
        ts, te = self.to_datetime(self._times_as_np_datetime64[[0, -1]])
        try:
            for interval in self.tsm_affected_intervals[self.spacecraft_id]:
                if ts >= interval[0] and te <= interval[1]:
                    # Found a matching interval
                    return True

            # No matching interval, orbit is not affected
            return False
        except KeyError:
            # Platform is not affected at all
            return False

    def get_midnight_scanline(self):
        """Find the scanline where the UTC date increases by one day.

        Returns:
            int: The midnight scanline if it exists and is unique.
                 None, else.

        """
        self.get_times()
        d0 = np.datetime64(datetime.date(1970, 1, 1), "D")
        days = (self._times_as_np_datetime64.astype("datetime64[D]") - d0).astype(int)
        incr = np.where(np.diff(days) == 1)[0]
        if len(incr) != 1:
            if len(incr) > 1:
                LOG.warning("Unable to determine midnight scanline: UTC date increases more than once. ")
            return None
        else:
            return incr[0]

    def get_miss_lines(self):
        """Find missing scanlines.

        I.e. scanlines which were dropped for some reason or were never recorded.

        Returns:
            Indices of missing scanlines

        """
        # Compare scanline number against the ideal case (1, 2, 3, ...) and
        # find the missing line numbers.
        ideal = set(range(1, self.scans["scan_line_number"][-1] + 1))
        missing = sorted(ideal.difference(set(self.scans["scan_line_number"])))
        return np.array(missing, dtype=int)

    def mask_tsm_pixels(self, ds):
        """Mask pixels affected by the scan motor issue."""
        idx = self.get_tsm_pixels(ds["channels"].values)
        # This is because fancy/pixel indexing doesn't seem to work with xarray's `where` or `loc`
        ds["channels"].data[idx] = np.nan

    @abstractmethod
    def get_tsm_pixels(self, channels):  # pragma: no cover
        """Determine pixels affected by the scan motor issue.

        Channel selection is POD/KLM specific.
        """
        raise NotImplementedError

    def get_attitude_coeffs(self):
        """Return the roll, pitch, yaw values."""
        if self._rpy is None:
            if "constant_yaw_attitude_error" in self.head.dtype.fields:
                rpy = np.deg2rad(
                    [
                        self.head["constant_roll_attitude_error"] / 1e3,
                        self.head["constant_pitch_attitude_error"] / 1e3,
                        self.head["constant_yaw_attitude_error"] / 1e3,
                    ]
                )
            else:
                try:
                    # This needs to be checked thoroughly first
                    # rpy_spacecraft = rpy_coeffs[self.spacecraft_name]
                    # rpy = np.array([rpy_spacecraft['roll'],
                    #                 rpy_spacecraft['pitch'],
                    #                 rpy_spacecraft['yaw']])
                    # LOG.debug("Using static attitude correction")
                    raise KeyError
                except KeyError:
                    LOG.debug("Not applying attitude correction")
                    rpy = np.zeros(3)
            LOG.info("Using rpy: %s", str(rpy))
            self._rpy = rpy
        return self._rpy

    def _compute_lonlat_from_tles(self, utcs):
        """Compute lon lat values using pyorbital."""
        tic = datetime.datetime.now()
        sgeom = self.geoloc_definition(utcs.astype(datetime.datetime), self.lonlat_sample_points)
        t0 = utcs[0].astype(datetime.datetime)
        s_times = sgeom.times(t0)
        tle1, tle2 = self.get_tle_lines()

        rpy = self.get_attitude_coeffs()
        LOG.debug(f"Computing lon/lats with attitude {rpy}")
        pixels_pos = compute_pixels((tle1, tle2), sgeom, s_times, rpy)
        pos_time = get_lonlatalt(pixels_pos, s_times)

        lons, lats = pos_time[:2]

        pixels_per_line = len(self.lonlat_sample_points)
        lons = lons.reshape(-1, pixels_per_line)
        lats = lats.reshape(-1, pixels_per_line)
        # todo: adjust time using build in lon lats
        toc = datetime.datetime.now()
        LOG.warning("Computation of geolocation: %s", str(toc - tic))

        return lons, lats


def inherit_doc(cls):
    """Make a class method inherit its docstring from the parent class.

    Copied from http://stackoverflow.com/a/8101598/5703449 .
    """
    for name, func in vars(cls).items():
        if isinstance(func, types.FunctionType) and not func.__doc__:
            for parent in cls.__bases__:
                parfunc = getattr(parent, name, None)
                if parfunc and getattr(parfunc, "__doc__", None):
                    func.__doc__ = parfunc.__doc__
                    break
    return cls
