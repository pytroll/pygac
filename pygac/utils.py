#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>
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

import gzip
import logging
import numpy as np
import sys
from contextlib import contextmanager

LOG = logging.getLogger(__name__)


def is_file_object(filename):
    """Check if the input is a file object.

    Args:
        filename - object to check

    Note:
        This method only check if the object implements the
        interface of a file object to allow duck types like
        gzip.GzipFile instances.
    """
    has_close = hasattr(filename, 'close')
    has_read = hasattr(filename, 'read')
    if hasattr(filename, 'seekable'):
        is_seekable = filename.seekable()
    else:
        is_seekable = False
    return has_close and has_read and is_seekable


@contextmanager
def _file_opener(file):
    """Open a file depending on the input.

    Args:
        file - path to file or file object
    """
    # open file if necessary
    if is_file_object(file):
        open_file = file
        close = False
    else:
        open_file = open(file, mode='rb')
        close = True
    # check if it is a gzip file
    try:
        file_object = gzip.open(open_file)
        file_object.read(1)
    except OSError:
        file_object = open_file
    finally:
        file_object.seek(0)
    # provide file_object with the context
    try:
        yield file_object
    finally:
        if close:
            file_object.close()


@contextmanager
def _file_opener_py2(file):
    """Open a file depending on the input.

    Args:
        file - path to file
    """
    close = True
    # check if it is a gzip file
    try:
        file_object = gzip.open(file)
        file_object.read(1)
    # Note: in python 2, this is an IOError, but we keep the
    #       OSError for testing.
    except (OSError, IOError):
        file_object = open(file, mode='rb')
    except TypeError:
        # In python 2 gzip.open cannot handle file objects
        LOG.debug("Gzip cannot open file objects in python2!")
        if is_file_object(file):
            file_object = file
            close = False
    finally:
        file_object.seek(0)
    # provide file_object with the context
    try:
        yield file_object
    finally:
        if close:
            file_object.close()


if sys.version_info.major < 3:
    file_opener = _file_opener_py2
else:
    file_opener = _file_opener


def get_absolute_azimuth_angle_diff(sat_azi, sun_azi):
    """Calculates absolute azimuth difference angle. """
    rel_azi = abs(sat_azi - sun_azi)
    rel_azi = rel_azi % 360
    # Not using np.where to avoid copying array
    rel_azi[rel_azi > 180] = 360.0 - rel_azi[rel_azi > 180]
    return rel_azi


def centered_modulus(array, divisor):
    """Transform array to half open range ]-divisor/2, divisor/2]."""
    arr = array % divisor
    arr[arr > divisor / 2] -= divisor
    return arr


def calculate_sun_earth_distance_correction(jday):
    """Calculate the sun earth distance correction.

    In 2008 3-4 different equations of ESD were considered.
    This one was chosen as it at the time gave reflectances most closely
    matching the PATMOS-x data provided then by Andy Heidinger.

    Formula might need to be reconsidered if jday is updated to a float.

    """
    # Earth-Sun distance correction factor
    corr = 1.0 - 0.0334 * np.cos(2.0 * np.pi * (jday - 2) / 365.25)
    return corr


def check_user_scanlines(start_line, end_line, first_valid_lat=None,
                         last_valid_lat=None, along_track=None):
    """Check user-defined scanlines.

    Can be used by both pygac and satpy.

    Args:
        start_line: User-defined start line (afer stripping, if enabled)
        end_line: User-defined end line (afer stripping, if enabled)
        first_valid_lat: First scanline with valid latitudes
        last_valid_lat: Last scanline with valid latitudes
        along_track: Number of scanlines (only needed if stripping
            is disabled)
    """
    if first_valid_lat is not None and last_valid_lat is not None:
        num_valid_lines = last_valid_lat - first_valid_lat + 1
    else:
        if along_track is None:
            raise ValueError('Need along_track')
        num_valid_lines = along_track

    start_line = int(start_line)
    end_line = int(end_line)
    if end_line == 0:
        # If the user specifies 0 as the last scanline, process all
        # scanlines with valid coordinates
        end_line = num_valid_lines - 1
    elif end_line >= num_valid_lines:
        end_line = num_valid_lines - 1
        LOG.warning('Given end line exceeds scanline range, resetting '
                    'to {}'.format(end_line))
    if start_line > num_valid_lines:
        raise ValueError('Given start line {} exceeds scanline range {}'
                         .format(start_line, num_valid_lines))
    return start_line, end_line


def strip_invalid_lat(lats):
    """Strip invalid latitudes at the end and beginning of the orbit."""
    no_wrong_lat = np.where(np.logical_not(np.isnan(lats)))
    return min(no_wrong_lat[0]), max(no_wrong_lat[0])


def slice_channel(ch, start_line, end_line, first_valid_lat=None, last_valid_lat=None):
    """Slice channel data using user-defined start/end line.

    If valid_lat_start/end are given, strip scanlines with invalid
    coordinates at the beginning and end of the orbit.

    Can be used by both pygac and satpy.

    Args:
        ch: Channel data
        start_line: User-defined start line (afer stripping, if enabled)
        end_line: User-defined end line (after stripping, if enabled)
        first_valid_lat: First scanline with valid latitudes
        last_valid_lat: Last scanline with valid latitudes.
    """
    if first_valid_lat is not None and last_valid_lat is not None:
        # Strip invalid coordinates
        ch = _slice(ch, start_line=first_valid_lat, end_line=last_valid_lat)

        # Reset user-defined end line, if it has been removed
        end_line = min(end_line, ch.shape[0] - 1)
        start_line = min(start_line, ch.shape[0] - 1)

    # Slice data using user-defined start/end lines
    ch_slc = _slice(ch, start_line=start_line, end_line=end_line)

    return ch_slc


def _slice(ch, start_line, end_line, update=None):
    """Slice the given channel.

    Args:
        start_line: New start line
        end_line: New end line
        update: List of scanlines to be updated to the new range
    """
    # Slice data using new start/end lines
    if len(ch.shape) == 1:
        ch_slc = ch[start_line:end_line + 1].copy()
    else:
        ch_slc = ch[start_line:end_line + 1, :].copy()

    if update:
        updated = [_update_scanline(l, start_line, end_line)
                   if l is not None else None
                   for l in update]
        return ch_slc, updated

    return ch_slc


def _update_scanline(scanline, new_start_line, new_end_line):
    """Update the given scanline to the new range.

    Set scanline to None if it lies outside the new range.
    """
    scanline -= new_start_line
    num_lines = new_end_line - new_start_line + 1
    if scanline < 0 or scanline >= num_lines:
        scanline = None
    return scanline


def plot_correct_times_thresh(res, filename=None):
    """Visualize results of GACReader.correct_times_thresh."""
    import matplotlib.pyplot as plt

    t = res['t']
    tcorr = res.get('tcorr')
    n = res['n']
    offsets = res.get('offsets')
    t0_head = res.get('t0_head')
    max_diff_from_t0_head = res.get('max_diff_from_t0_head')
    fail_reason = res.get('fail_reason', 'Failed for unknown reason')

    # Setup figure
    along_track = np.arange(t.size)
    _, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True,
                                      figsize=(8, 10))

    # Plot original vs corrected timestamps
    ax0.plot(along_track, t, "b-", label="original")
    if tcorr is not None:
        ax0.plot(along_track, tcorr, color="red", linestyle="--",
                 label="corrected")
    else:
        ax0.set_title(fail_reason)

    ax0.set_ylabel("Time")
    ax0.set_ylim(t.min() - np.timedelta64(30, "m"),
                 t.max() + np.timedelta64(30, "m"))
    ax0.legend(loc="best")

    # Plot offset (original time - ideal time)
    if offsets is not None:
        ax1.plot(along_track, offsets)
        ax1.fill_between(
            along_track,
            t0_head - np.ones(along_track.size) * max_diff_from_t0_head,
            t0_head + np.ones(along_track.size) * max_diff_from_t0_head,
            facecolor="g", alpha=0.33)
        ax1.axhline(y=t0_head, color="g", linestyle="--",
                    label="Header timestamp")
        ax1.set_ylim(t0_head - 5 * max_diff_from_t0_head,
                     t0_head + 5 * max_diff_from_t0_head)
        ax1.set_ylabel("Offset t-tn [ms]")
        ax1.legend(loc="best")

    # Plot scanline number
    ax2.plot(along_track, n)
    ax2.set_ylabel("Scanline number")
    ax2.set_xlabel("Along Track")

    if filename:
        plt.savefig(filename, bbox_inches="tight", dpi=100)
    else:
        plt.show()


def plot_correct_scanline_numbers(res, filename=None):
    """Visualize results of GACReader.correct_scanline_numbers."""
    import matplotlib.pyplot as plt

    along_track = res['along_track']
    n_orig = res['n_orig']
    n_corr = res['n_corr']
    within_range = res['within_range']
    thresh = res['thresh']
    diffs = res['diffs']
    nz_diffs = res['nz_diffs']

    # Setup figure
    _, (ax0, ax1) = plt.subplots(nrows=2)

    # Plot original vs corrected scanline numbers
    ax0.plot(along_track, n_orig, "b-", label="original")
    along_track_corr = along_track.copy()
    along_track_corr = along_track_corr[within_range]
    along_track_corr = along_track_corr[diffs <= thresh]
    ax0.plot(along_track_corr, n_corr, "r--", label="corrected")
    ax0.set_ylabel("Scanline Number")
    ax0.set_xlabel("Along Track")
    ax0.legend(loc="best")

    # Plot difference from ideal
    ax1.plot(np.arange(len(nz_diffs)), nz_diffs)
    ax1.axhline(thresh, color="r", label="thresh={0:.2f}"
                .format(thresh))
    ax1.set_xlabel("Index")
    ax1.set_ylabel("nonzero |n - n'|")
    ax1.legend()

    plt.tight_layout()

    if filename:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()
