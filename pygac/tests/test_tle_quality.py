#!/usr/bin/env python

# Copyright (c) 2025 Pytroll Developers

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

"""Unit tests for finding and repairing unusable element sets in a two-line element file."""

from pygac.tle_quality import consistent_elements, corrected_elements

#: Five consecutive NOAA-11 element sets from a quiet week of January 1993.
NOAA11_QUIET_WEEK = [
    ("1 19531U 88089  A 93001.17990934  .00000197 +00000-0 +11722-3 0  9990",
     "2 19531 099.1102 333.8698 0011469 185.1327 174.9728 14.12802666220098"),
    ("1 19531U 88089  A 93002.10057474  .00000265 +00000-0 +15314-3 0  9994",
     "2 19531 099.1106 334.8066 0011460 182.4064 177.7071 14.12803963220229"),
    ("1 19531U 88089  A 93002.87959859  .00000342 +00000-0 +19522-3 0  9997",
     "2 19531 099.1113 335.5993 0011480 180.2707 179.8480 14.12805221220334"),
    ("1 19531U 88089  A 93003.30452059  .00000339 +00000-0 +20366-3 0  9990",
     "2 19531 099.1113 336.0313 0011525 178.8010 181.3184 14.12805643220399"),
    ("1 19531U 88089  A 93003.65862264  .00000346 +00000-0 +20672-3 0  9999",
     "2 19531 099.1111 336.3911 0011556 177.9525 182.1696 14.12805903220441"),
]


def test_mutually_consistent_elements_are_all_kept():
    """Elements that agree with each other must survive untouched."""
    assert consistent_elements(NOAA11_QUIET_WEEK) == NOAA11_QUIET_WEEK


#: Two real NOAA-10 element sets from September 1992, followed by a catalogue
#: entry filed under the same object that cannot describe it: 2.13 revolutions a
#: day at 64.8 degrees inclination is a twelve-hour orbit, not a polar one.
NOAA10_GOOD_PAIR = [
    ("1 16969U 86073  A 92260.04271183  .00000103 +00000-0 +48210-4 0  9993",
     "2 16969 098.5333 277.0138 0014111 039.3874 320.8307 14.24700897311598"),
    ("1 16969U 86073  A 92260.53432119  .00000104 +00000-0 +48210-4 0  9996",
     "2 16969 098.5326 277.4920 0013959 038.4102 321.8080 14.24701023311665"),
]
NOAA10_IMPOSSIBLE_ORBIT = (
    "1 16969U 86073  A 92260.87302325  .00000018 +00000-0 +99999-4 0  9992",
    "2 16969 064.8348 132.0720 0003144 216.8779 143.2746 02.13102711046693",
)


def test_a_physically_impossible_orbit_is_rejected():
    """An element set that no AVHRR platform could fly must not survive."""
    assert consistent_elements(NOAA10_GOOD_PAIR + [NOAA10_IMPOSSIBLE_ORBIT]) == NOAA10_GOOD_PAIR


#: Seven real NOAA-10 element sets spanning August 1989. Six of them agree on
#: where the platform was to within 0.6 km. The first, listed before them here,
#: places it 179 km along track from that agreement, which at orbital speed is
#: 24 seconds of flight: exactly the TAI minus UTC offset in force in 1989.
NOAA10_AUGUST_1989 = [
    ("1 16969U 86073  A 89216.24584550  .00000500 +00000-0 +22873-3 0  9994",
     "2 16969 098.6373 245.8217 0014449 108.1012 252.1736 14.23106590150770"),
    ("1 16969U 86073  A 89217.37050559  .00000494 +00000-0 +22638-3 0  9997",
     "2 16969 098.6373 246.9253 0014488 104.9789 255.2993 14.23107603150933"),
    ("1 16969U 86073  A 89217.58143139  .00000492 +00000-0 +22537-3 0  9993",
     "2 16969 098.6373 247.1319 0014492 104.4648 255.8133 14.23107775150968"),
    ("1 16969U 86073  A 89218.07359174  .00000489 +00000-0 +22440-3 0  9995",
     "2 16969 098.6373 247.6148 0014507 103.1080 257.1719 14.23108213151038"),
    ("1 16969U 86073  A 89219.05791088  .00000486 +00000-0 +22316-3 0  9997",
     "2 16969 098.6373 248.5802 0014513 100.3928 259.8887 14.23109137151178"),
    ("1 16969U 86073  A 89221.23747386  .00000490 +00000-0 +22498-3 0  9998",
     "2 16969 098.6373 250.7183 0014578 094.4440 265.8403 14.23111455151489"),
    ("1 16969U 86073  A 89223.41702944  .00000506 +00000-0 +23129-3 0  9991",
     "2 16969 098.6372 252.8564 0014593 088.6193 271.6659 14.23114085151795"),
]


def test_an_element_that_contradicts_the_majority_is_rejected():
    """An element placing the platform far from where its neighbours agree it was must not survive."""
    assert consistent_elements(NOAA10_AUGUST_1989) == NOAA10_AUGUST_1989[1:]


def test_an_element_on_the_wrong_time_scale_is_corrected_rather_than_discarded():
    """A contaminated epoch is a repairable defect: the element must be kept, put right."""
    repaired = corrected_elements(NOAA10_AUGUST_1989)
    assert len(repaired) == len(NOAA10_AUGUST_1989)
    assert consistent_elements(repaired) == repaired


#: Five real NOAA-6 element sets from December 1986. The first, listed before the
#: others, lags them by 171 km. TAI ran 23 seconds ahead of UTC in 1986, and a
#: shift of that size brings it to within 0.25 km of the other four; a shift of
#: the 24 seconds that applied from 1988 leaves it 7.69 km out.
NOAA6_DECEMBER_1986 = [
    ("1 11416U 79057A   86339.86646901  .00000112 +00000-0 +51913-4 0  9997",
     "2 11416 098.4976 349.9025 0011124 298.0771 061.9284 14.24968040386244"),
    ("1 11416U 79057A   86340.91945807  .00000117 +00000-0 +54073-4 0  9997",
     "2 11416 098.4974 350.9217 0011084 294.6222 065.3803 14.24968456386397"),
    ("1 11416U 79057A   86340.98967501  .00000119 +00000-0 +54557-4 0  9998",
     "2 11416 098.4974 350.9897 0011089 294.3982 065.6044 14.24968497386403"),
    ("1 11416U 79057A   86341.27054280  .00000118 +00000-0 +54323-4 0  9992",
     "2 11416 098.4974 351.2616 0011076 293.5439 066.4578 14.24968556386441"),
    ("1 11416U 79057A   86341.55141068  .00000112 +00000-0 +51756-4 0  9995",
     "2 11416 098.4973 351.5331 0011076 292.5518 067.4490 14.24968546386481"),
]


def _epoch_of(element):
    """Return the epoch recorded on *element*'s first line, in days."""
    return float(element[0][18:32])


def test_the_correction_uses_the_leap_second_count_of_the_elements_own_date():
    """TAI ran 23 seconds ahead of UTC in 1986, not the 24 that applied from 1988."""
    repaired = corrected_elements(NOAA6_DECEMBER_1986)
    moved_seconds = (_epoch_of(NOAA6_DECEMBER_1986[0]) - _epoch_of(repaired[0])) * 86400
    assert abs(moved_seconds - 23.0) < 0.05


#: Four real NOAA-10 element sets from January 1987. Here the contaminated group
#: is the MAJORITY: the three listed last lag the first by 23 seconds, the TAI
#: minus UTC offset of 1987. Shifting those three back by 23 s brings them to
#: within 0.14 km of the first, which is the element that was right all along.
#: Repairing whatever disagrees with the median would move the wrong one, and
#: put it 342 km out instead of 171 km.
NOAA10_JANUARY_1987 = [
    ("1 16969U 86073  A 87012.25905882  .00000049 +00000-0 +25947-4 0  9996",
     "2 16969 098.7386 044.3114 0013174 301.8062 058.1830 14.22483979016570"),
    ("1 16969U 86073  A 87012.61102108  .00000046 +00000-0 +25146-4 0  9994",
     "2 16969 098.7385 044.6596 0013164 300.7238 059.2637 14.22483997016621"),
    ("1 16969U 86073  A 87014.43984069  .00000041 +00000-0 +22870-4 0  9996",
     "2 16969 098.7383 046.4721 0013111 295.2685 064.7134 14.22484091016882"),
    ("1 16969U 86073  A 87016.62035598  .00000041 +00000-0 +22787-4 0  9990",
     "2 16969 098.7381 048.6330 0013011 288.6482 071.3278 14.22484336017190"),
]


def test_the_lagging_elements_are_repaired_even_when_they_are_the_majority():
    """It is lagging by the leap-second offset that marks an element, not being outnumbered."""
    repaired = corrected_elements(NOAA10_JANUARY_1987)
    assert repaired[0] == NOAA10_JANUARY_1987[0]
    assert consistent_elements(repaired) == repaired


#: Six real NOAA-9 element sets from August 1987, when the platform was losing
#: height fast: their drag terms are a hundred times those of the other fixtures,
#: and the first four trail the last by about 12 seconds. TAI led UTC by 23
#: seconds in 1987, so a 12 second lag is some other trouble entirely. Shifting
#: these by 23 seconds would overshoot by 11, and push them further out than
#: they already are.
NOAA9_AUGUST_1987 = [
    ("1 15427U 84123  A 87239.62344547  .00059071 +00000-0 +32512-1 0  9997",
     "2 15427 099.0572 205.8903 0016246 009.8487 350.2982 14.11556334139392"),
    ("1 15427U 84123  A 87239.97785637  .00059115 +00000-0 +32512-1 0  9993",
     "2 15427 099.0576 206.2481 0016028 009.7625 350.3842 14.11597945139448"),
    ("1 15427U 84123  A 87240.19049738  .00059143 +00000-0 +32512-1 0  9995",
     "2 15427 099.0581 206.4646 0016210 008.7944 351.3500 14.11622952139477"),
    ("1 15427U 84123  A 87240.82839625  .00059058 +00000-0 +32512-1 0  9992",
     "2 15427 099.0569 207.1071 0015549 008.2933 351.3293 14.11549639139569"),
    ("1 15427U 84123  A 87241.96257553  .00059597 +00000-0 +32512-1 0  9990",
     "2 15427 099.0572 208.2541 0015682 001.8363 357.8676 14.11527534139725"),
    ("1 15427U 84123  A 87242.17531257  .00059039 +00000-0 +32512-1 0  9991",
     "2 15427 099.0568 208.4733 0015118 000.6933 359.4294 14.11535815139757"),
]


def test_a_lag_that_is_not_a_leap_second_is_left_alone():
    """Only a lag of exactly the leap-second offset marks an epoch on the wrong time scale."""
    assert corrected_elements(NOAA9_AUGUST_1987) == NOAA9_AUGUST_1987
