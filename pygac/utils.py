import numpy as np


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
