import logging
import numpy as np


LOG = logging.getLogger(__name__)


def slice_channel(ch, start_line, end_line, valid_lat_start=None,
                  valid_lat_end=None, midnight_scanline=None,
                  miss_lines=None, qual_flags=None):
    """Slice channel data using user-defined start/end line.

    If valid_lat_start/end are given, strip scanlines with invalid
    coordinates at the beginning and end of the orbit.

    This method is doing too much at once, but it ensures that the
    same slicing method is being used by save_gac and the satpy reader.

    Args:
        ch: Channel data
        start_line: User-defined start line
        end_line: User-defined end line
        valid_lat_start: First scanline with valid latitudes
        valid_lat_end: Last scanline with valid latitudes.
        midnight_scanline: If given, update midnight scanline to the new
            scanline range.
        miss_lines: If given, update list of missing lines with the ones
            that have been stripped due to invalid coordinates
        qual_flags: Quality flags, needed to updated missing lines.
    """
    # Strip invalid coordinates
    if valid_lat_start is None or valid_lat_end is None:
        valid_lat_start, valid_lat_end = 0, lats.shape[0]

    # Update start/end lines
    start_line, end_line = update_start_end_line(user_start=start_line,
                                                 user_end=end_line,
                                                 valid_lat_start=valid_lat_start,
                                                 valid_lat_end=valid_lat_end)

    # Slice data using new start/end lines
    if len(ch.shape) == 1:
        ch_slc = ch[start_line:end_line + 1].copy()
    else:
        ch_slc = ch[start_line:end_line + 1, :].copy()

    if miss_lines is None and midnight_scanline is None:
        return ch_slc
    else:
        # Add stripped lines to list of missing scanlines
        if miss_lines is not None:
            if qual_flags is None:
                raise ValueError('Need qual_flags, too')
            miss_lines = update_missing_scanlines(miss_lines=miss_lines,
                                                  qual_flags=qual_flags,
                                                  valid_lat_start=valid_lat_start,
                                                  valid_lat_end=valid_lat_end)

        # Update midnight scanline
        if midnight_scanline is not None:
            midnight_scanline = update_scanline(midnight_scanline,
                                                new_start_line=start_line,
                                                new_end_line=end_line)

        return ch_slc, miss_lines, midnight_scanline


def strip_invalid_lat(lats):
    """Strip invalid latitudes at the end and beginning of the orbit."""
    no_wrong_lat = np.where(np.logical_not(np.isnan(lats)))
    return min(no_wrong_lat[0]), max(no_wrong_lat[0])


def update_start_end_line(user_start, user_end, valid_lat_start, valid_lat_end):
    """Update user start/end lines after stripping invalid latitudes.

    Returns:
        Updated start_line, updated end_line
    """
    new_start_line = max(user_start, valid_lat_start)
    new_end_line = min(user_end, valid_lat_end)
    return new_start_line, new_end_line


def update_scanline(scanline, new_start_line, new_end_line):
    """Update the given scanline to the new range.

    Set scanline to None if it lies outside the new range.
    """
    scanline -= new_start_line
    num_lines = new_end_line - new_start_line + 1
    if scanline < 0 or scanline >= num_lines:
        scanline = None
    return scanline


def update_missing_scanlines(miss_lines, qual_flags, valid_lat_start, valid_lat_end):
    return np.sort(np.array(
        qual_flags[0:valid_lat_start, 0].tolist() +
        miss_lines.tolist() +
        qual_flags[valid_lat_end+1:, 0].tolist()
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
