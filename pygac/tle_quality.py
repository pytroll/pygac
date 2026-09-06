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

"""Finding and repairing unusable element sets in a two-line element file."""

import numpy as np
from pyorbital.orbital import Orbital

#: An AVHRR platform circles the Earth about fourteen times a day. Anything
#: appreciably slower is in a different orbit, and so describes a different object.
SLOWEST_AVHRR_MEAN_MOTION = 12.0

#: How far an element may place the platform from where its neighbours agree it
#: was before it is taken to describe some other flight. Elements of one flight
#: agree to a couple of kilometres; an epoch kept on the wrong time scale puts
#: the platform some 160 km along track, so the two are never close to this line.
FURTHEST_AGREEING_POSITION_KM = 50.0


def _mean_motion(element):
    """Return the revolutions per day *element* records for its object."""
    second_line = element[1]
    return float(second_line[52:63])


def _orbit_of(element):
    """Return the orbit *element* describes."""
    return Orbital("", line1=element[0], line2=element[1])


def _midmost_epoch(orbits):
    """Return the instant lying midway through the epochs of *orbits*."""
    epochs = np.array([orbit.tle.epoch for orbit in orbits], dtype="datetime64[s]")
    return np.datetime64(int(np.median(epochs.astype("int64"))), "s").astype(object)


def _propagated_together(orbits):
    """Return where and how fast *orbits* are, all at the middle of their epochs."""
    when = _midmost_epoch(orbits)
    states = [orbit.get_position(when, normalize=False) for orbit in orbits]
    return (np.array([position for position, _ in states]),
            np.array([velocity for _, velocity in states]))


def consistent_elements(elements):
    """Return those *elements* that agree with one another about the flight they describe."""
    flyable = [element for element in elements
               if _mean_motion(element) >= SLOWEST_AVHRR_MEAN_MOTION]
    positions, _ = _propagated_together([_orbit_of(element) for element in flyable])
    agreed = np.median(positions, axis=0)
    return [element for element, position in zip(flyable, positions)
            if np.linalg.norm(position - agreed) <= FURTHEST_AGREEING_POSITION_KM]


#: How far a lag may sit from the leap-second offset and still be taken for one.
#: Repairs measured across the record land within 1.4 s of the exact value, while
#: other troubles miss it by ten seconds or more.
LEAP_SECOND_TOLERANCE_S = 2.0

#: Seconds by which TAI ran ahead of UTC, by the date each value took effect.
LEAP_SECONDS = ((np.datetime64("1985-07-01"), 23),
                (np.datetime64("1988-01-01"), 24))


def _with_epoch_shifted(element, seconds):
    """Return *element* with its epoch moved by *seconds* and its checksum made good."""
    line = element[0]
    moved = float(line[18:32]) + seconds / 86400.0
    body = f"{line[:18]}{moved:14.8f}{line[32:68]}"
    digits = sum(int(c) if c.isdigit() else 1 if c == "-" else 0 for c in body)
    return (body + str(digits % 10), element[1])


def _tai_minus_utc(orbit):
    """Return the seconds by which TAI led UTC when *orbit* was catalogued."""
    when = orbit.tle.epoch
    return [seconds for start, seconds in LEAP_SECONDS if when >= start][-1]


def corrected_elements(elements):
    """Return *elements* with those kept on the wrong time scale put right."""
    orbits = [_orbit_of(element) for element in elements]
    positions, velocities = _propagated_together(orbits)
    heading = np.median(velocities, axis=0)
    speed = np.linalg.norm(heading)
    along_track = positions @ (heading / speed)
    lags = (along_track.max() - along_track) / speed
    offsets = [_tai_minus_utc(orbit) for orbit in orbits]
    return [_with_epoch_shifted(element, -offset)
            if abs(lag - offset) <= LEAP_SECOND_TOLERANCE_S else element
            for element, lag, offset in zip(elements, lags, offsets)]
