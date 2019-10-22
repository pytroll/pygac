#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author(s):

#   Stephan Finkensieper <stephan.finkensieper@dwd.de>

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

import logging
import numpy as np


LOG = logging.getLogger(__name__)


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


def slice_channel(ch, start_line, end_line, first_valid_lat=None,
                  last_valid_lat=None, midnight_scanline=None,
                  miss_lines=None, qual_flags=None):
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
        midnight_scanline: If given, update midnight scanline to the new
            scanline range.
        miss_lines: If given, update list of missing lines with the ones
            that have been stripped due to invalid coordinates
        qual_flags: Quality flags, needed to updated missing lines.
    """
    if first_valid_lat is not None and last_valid_lat is not None:
        # Strip invalid coordinates and update midnight scanline as well as
        # user-defined start/end lines
        ch, updated = _slice(ch,
                             start_line=first_valid_lat,
                             end_line=last_valid_lat,
                             update=[midnight_scanline])
        midnight_scanline = updated[0]

        # Reset user-defined end line, if it has been removed
        end_line = min(end_line, ch.shape[0] - 1)
        start_line = min(start_line, ch.shape[0] - 1)

        # Update missing scanlines
        if miss_lines is not None:
            miss_lines = _update_missing_scanlines(
                miss_lines=miss_lines,
                qual_flags=qual_flags,
                start_line=first_valid_lat,
                end_line=last_valid_lat)

    # Slice data using user-defined start/end lines
    ch_slc, updated = _slice(ch, start_line=start_line, end_line=end_line,
                             update=[midnight_scanline])
    midnight_scanline = updated[0]

    return ch_slc, miss_lines, midnight_scanline


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


def _update_missing_scanlines(miss_lines, qual_flags, start_line, end_line):
    """Add scanlines excluded by slicing to the list of missing scanlines.

    Args:
        miss_lines: List of missing scanlines
        qual_flags: Quality flags
        start_line: New start line of the slice
        end_line: New end line of the slice
    """
    return np.sort(np.unique(
        qual_flags[0:start_line, 0].tolist() +
        miss_lines.tolist() +
        qual_flags[end_line + 1:, 0].tolist()
    ))


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
