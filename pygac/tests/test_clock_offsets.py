"""Tests for the tables of measured spacecraft clock error."""

import numpy as np
import pytest

from pygac.clock_offsets_converter import MEASUREMENTS, clock_table_covers, get_offsets, txt


def test_a_time_inside_the_measurements_is_covered():
    """NOAA-9's table was measured across 1993, so a pass from then is not guesswork."""
    assert clock_table_covers("noaa9", np.datetime64("1993-06-01T12:00:00"))


def test_a_platform_with_no_table_is_covered_nowhere():
    """noaa6 has no measured clock error at all, so nothing about it is known."""
    assert not clock_table_covers("noaa6", np.datetime64("1981-05-15T09:00:00"))


def test_a_moment_past_the_end_of_a_table_is_not_covered():
    """noaa9's measurements stop in August 1995; after that nothing was recorded."""
    assert not clock_table_covers("noaa9", np.datetime64("1998-06-01T00:00:00"))


def test_every_table_runs_forwards_in_time():
    """np.interp needs its points in order; given them out of order it returns nonsense."""
    for spacecraft_name in sorted(txt):
        measured, _ = get_offsets(spacecraft_name)
        backwards = [(a, b) for a, b in zip(measured, measured[1:]) if b < a]
        assert not backwards, f"{spacecraft_name} steps backwards at {backwards[:1]}"


def test_a_table_that_has_lost_a_measurement_is_refused(monkeypatch):
    """These tables are edited by hand, and a dropped line changes navigation silently."""
    monkeypatch.setitem(txt, "noaa9", "\n".join(txt["noaa9"].strip().split("\n")[:-1]))
    with pytest.raises(ValueError, match="measurement"):
        get_offsets("noaa9")


def test_a_moment_before_a_table_begins_is_not_covered():
    """noaa9's measurements start in January 1986; before that nothing was recorded."""
    assert not clock_table_covers("noaa9", np.datetime64("1980-01-01T00:00:00"))


def test_every_table_records_how_many_measurements_it_holds():
    """A table nobody counted can lose a line without anything failing."""
    uncounted = sorted(set(txt) - set(MEASUREMENTS))
    assert not uncounted, f"no measurement count recorded for {uncounted}"


def test_noaa10_clock_error_is_known_across_its_avhrr_years():
    """NOAA-10 drifted some six milliseconds a day between corrections.

    Left unmeasured, a pass halfway through one of those intervals is navigated
    with most of a second of clock error in it, which is several kilometres along
    the track.
    """
    for when in ("1989-09-01T00:00:00", "1990-06-01T12:00:00", "1991-09-01T00:00:00"):
        assert clock_table_covers("noaa10", np.datetime64(when))
