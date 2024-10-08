#!/usr/bin/env python

# Author(s):

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

"""Convert PATMOS-x calibration file tarballs to PyGAC calibration json format.

The official tarballs are available on the PATMOS-x webpage "https://cimss.ssec.wisc.edu/patmosx/avhrr_cal.html".
"""

import argparse
import datetime as dt
import json
import logging
import pathlib
import re
import tarfile

from scipy.optimize import bisect


class PatmosxReader:
    """Read PATMOS-x coefficient files tarballs."""
    # regular expression with named capturing groups to read an entire patmosx file
    regex = re.compile(
        r"\s*(?P<sat_name>\w+)[^\n]*\n"
        r"\s*(?P<solar_3b>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ew_3b>[eE0-9\.-]+)[^\n]*\n"
        r'\s*(?P<ch3b_Ns>[eE0-9\.-]+),?\s*(?P<ch3b_b0>[eE0-9\.-]+),?\s*(?P<ch3b_b1>[eE0-9\.-]+),?\s*(?P<ch3b_b2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch4_Ns>[eE0-9\.-]+),?\s*(?P<ch4_b0>[eE0-9\.-]+),?\s*(?P<ch4_b1>[eE0-9\.-]+),?\s*(?P<ch4_b2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch5_Ns>[eE0-9\.-]+),?\s*(?P<ch5_b0>[eE0-9\.-]+),?\s*(?P<ch5_b1>[eE0-9\.-]+),?\s*(?P<ch5_b2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r"\s*(?P<nu_3b>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a1_3b>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a2_3b>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<nu_4>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a1_4>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a2_4>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<nu_5>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a1_5>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<a2_5>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch1_dark_count>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch2_dark_count>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch3a_dark_count>[eE0-9\.-]+)[^\n]*\n"
        r"(?:[a-z]+[^\n]*\n)?"
        r'\s*(?P<ch1_low_gain_S0>[eE0-9\.-]+)\s*(?P<ch1_low_gain_S1>[eE0-9\.-]+)\s*(?P<ch1_low_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch1_high_gain_S0>[eE0-9\.-]+)\s*(?P<ch1_high_gain_S1>[eE0-9\.-]+)\s*(?P<ch1_high_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch2_low_gain_S0>[eE0-9\.-]+)\s*(?P<ch2_low_gain_S1>[eE0-9\.-]+)\s*(?P<ch2_low_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch2_high_gain_S0>[eE0-9\.-]+)\s*(?P<ch2_high_gain_S1>[eE0-9\.-]+)\s*(?P<ch2_high_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch3a_low_gain_S0>[eE0-9\.-]+)\s*(?P<ch3a_low_gain_S1>[eE0-9\.-]+)\s*(?P<ch3a_low_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<ch3a_high_gain_S0>[eE0-9\.-]+)\s*(?P<ch3a_high_gain_S1>[eE0-9\.-]+)\s*(?P<ch3a_high_gain_S2>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r"\s*(?P<date_of_launch>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch1_gain_switches_count>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch2_gain_switches_count>[eE0-9\.-]+)[^\n]*\n"
        r"\s*(?P<ch3a_gain_switches_count>[eE0-9\.-]+)[^\n]*\n"
        r'\s*(?P<PRT1_0>[eE0-9\.-]+)\,*\s*(?P<PRT1_1>[eE0-9\.-]+)\,*\s*(?P<PRT1_2>[eE0-9\.-]+)\,*\s*(?P<PRT1_3>[eE0-9\.-]+)\,*\s*(?P<PRT1_4>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<PRT2_0>[eE0-9\.-]+)\,*\s*(?P<PRT2_1>[eE0-9\.-]+)\,*\s*(?P<PRT2_2>[eE0-9\.-]+)\,*\s*(?P<PRT2_3>[eE0-9\.-]+)\,*\s*(?P<PRT2_4>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<PRT3_0>[eE0-9\.-]+)\,*\s*(?P<PRT3_1>[eE0-9\.-]+)\,*\s*(?P<PRT3_2>[eE0-9\.-]+)\,*\s*(?P<PRT3_3>[eE0-9\.-]+)\,*\s*(?P<PRT3_4>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<PRT4_0>[eE0-9\.-]+)\,*\s*(?P<PRT4_1>[eE0-9\.-]+)\,*\s*(?P<PRT4_2>[eE0-9\.-]+)\,*\s*(?P<PRT4_3>[eE0-9\.-]+)\,*\s*(?P<PRT4_4>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'\s*(?P<PRT_weight_1>[eE0-9\.-]+)\s*(?P<PRT_weight_2>[eE0-9\.-]+)\s*(?P<PRT_weight_3>[eE0-9\.-]+)\s*(?P<PRT_weight_4>[eE0-9\.-]+)[^\n]*\n'  # noqa
        r'(?:\s*(?P<sst_mask_day_0>[eE0-9\.-]+),\s*(?P<sst_mask_day_1>[eE0-9\.-]+),\s*(?P<sst_mask_day_2>[eE0-9\.-]+),\s*(?P<sst_mask_day_3>[eE0-9\.-]+)[^\n]*\n)?'  # noqa
        r'(?:\s*(?P<sst_mask_night_0>[eE0-9\.-]+),\s*(?P<sst_mask_night_1>[eE0-9\.-]+),\s*(?P<sst_mask_night_2>[eE0-9\.-]+),\s*(?P<sst_mask_night_3>[eE0-9\.-]+)[^\n]*\n)?'  # noqa
        r"(?:\![^v][^\n]*\n)*"
        r"(?:\!(?P<version>v\w+))?"
    )

    def __init__(self, tarball):
        self.tarball = tarball
        self.coeffs = []
        with tarfile.open(tarball) as archive:
            for tarinfo in archive:
                if tarinfo.isfile():
                    # open satellite specific coefficient file
                    filename = tarinfo.name
                    fileobj = archive.extractfile(filename)
                    content = fileobj.read().decode()
                    match = self.regex.match(content)
                    sat_coeffs = {key: self.parse_types(value) for key, value in match.groupdict().items()}
                    self.coeffs.append(sat_coeffs)

    def __iter__(self):
        yield from self.coeffs  # noqa

    @staticmethod
    def parse_types(value):
        """parse the types of coefficients"""
        try:
            try:
                _value = int(value)
            except ValueError:
                _value = float(value)
        except ValueError:
            _value = value
        except TypeError:
            _value = None
        return _value


class Translator:
    """Translate PATMOS-x coefficients to PyGAC format."""
    sat_names = {"m01": "metopb", "m02": "metopa", "n05": "tirosn", "m03": "metopc"}
    sat_names.update({"n{0:02d}".format(i): "noaa{0}".format(i) for i in range(6,20)})
    description = {
        "visible": {
            "channels": ["1", "2", "3a"],
            "coefficients": {
                "dark_count": "instrument counts under dark conditions []",
                "gain_switch": "dual-gain switch count, set to 'null' for single-gain instruments []",
                "s0": "single-gain calibration slope at launch date [%]",
                "s1": "linear single-gain calibration slope parameter [% years^{-1}]",
                "s2": "quadratic single-gain calibration slope parameter [% years^{-2}]",
                "date_of_launch": "timestamp of launch date [UTC]"
            },
            "method": 'Heidinger, A.K., W.C. Straka III, C.C. Molling, J.T. Sullivan, and X. Wu, 2010: Deriving an inter-sensor consistent calibration for the AVHRR solar reflectance data record. International Journal of Remote Sensing, 31:6493-6517'  # noqa
        },
        "thermal": {
            "channels": ["3b", "4", "5"],
            "coefficients": {
                "centroid_wavenumber": "centroid wavenumber [cm^{-1}]",
                "b0": "constant non-linear radiance correction coefficient [mW m^{-2} sr cm^{-1}]",
                "b1": "linear non-linear radiance correction coefficient []",
                "b2": "quadratic non-linear radiance correction coefficient [(mW^{-1} m^2 sr^{-1} cm)]",
                "space_radiance": "radiance of space [mW m^{-2} sr cm^{-1}]",
                "to_eff_blackbody_intercept": "thermal channel temperature to effective blackbody temperature intercept [K]",  # noqa
                "to_eff_blackbody_slope": "thermal channel temperature to effective blackbody temperature slope []",
                "d0": "constant thermometer counts to temperature conversion coefficient [K]",
                "d1": "linear thermometer counts to temperature conversion coefficient [K]",
                "d2": "quadratic thermometer counts to temperature conversion coefficient [K]",
                "d3": "cubic thermometer counts to temperature conversion coefficient [K]",
                "d4": "quartic thermometer counts to temperature conversion coefficient [K]"
            },
            "method": "Goodrum, G., Kidwell, K.B. and W. Winston, 2000: NOAA KLM User's Guide. U.S. Department of Commerce, National Oceanic and Atmospheric Administration, National Environmental Satellite, Data and Information Service; Walton, C. C., J. T. Sullivan, C. R. N. Rao, and M. P. Weinreb, 1998: Corrections for detector nonlinearities and calibration inconsistencies of the infrared channels of the Advanced Very High Resolution Radiometer. J. Geophys. Res., 103, 3323-3337; Trishchenko, A.P., 2002: Removing Unwanted Fluctuations in the AVHRR Thermal Calibration Data Using Robust Techniques. Journal of Atmospheric and Oceanic Technology, 19:1939-1954"  # noqa
        }
    }

    def __init__(self, patmosx_coeffs):
        self.coeffs = {"description": self.description}
        for patmosx_sat_coeffs in patmosx_coeffs:
            sat_name = self.sat_names[patmosx_sat_coeffs["sat_name"]]
            pygac_sat_coeffs = self.convert(patmosx_sat_coeffs)
            self.coeffs[sat_name] = pygac_sat_coeffs

    @classmethod
    def convert(cls, patmosx_sat_coeffs):
        pygac_sat_coeffs = {}
        # visible calibration
        for ch in ("1", "2", "3a"):
            s0l = patmosx_sat_coeffs["ch{0}_low_gain_S0".format(ch)]
            s0h = patmosx_sat_coeffs["ch{0}_high_gain_S0".format(ch)]
            if s0l == s0h:
                gain_switch = None
                s0 = s0l
            else:
                gain_switch = patmosx_sat_coeffs["ch{0}_gain_switches_count".format(ch)]
                s0 = cls.find_s0(s0l, s0h, ch)
            pygac_sat_coeffs["channel_{0}".format(ch)] = {
                "dark_count": float(patmosx_sat_coeffs["ch{0}_dark_count".format(ch)]),
                "gain_switch": gain_switch,
                "s0": s0,
                "s1": patmosx_sat_coeffs["ch{0}_high_gain_S1".format(ch)],
                "s2": patmosx_sat_coeffs["ch{0}_high_gain_S2".format(ch)]
            }
        date_of_launch = cls.float2date(patmosx_sat_coeffs["date_of_launch"])
        pygac_sat_coeffs["date_of_launch"] = date_of_launch.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
        # thermal channels
        for ch in ("3b", "4", "5"):
            pygac_sat_coeffs["channel_{0}".format(ch)] = {
                "b0": patmosx_sat_coeffs["ch{0}_b0".format(ch)],
                "b1": patmosx_sat_coeffs["ch{0}_b1".format(ch)],
                "b2": patmosx_sat_coeffs["ch{0}_b2".format(ch)],
                "centroid_wavenumber": patmosx_sat_coeffs["nu_{0}".format(ch)],
                "space_radiance": patmosx_sat_coeffs["ch{0}_Ns".format(ch)],
                "to_eff_blackbody_intercept": (-patmosx_sat_coeffs["a1_{0}".format(ch)]
                                                / patmosx_sat_coeffs["a2_{0}".format(ch)]),
                "to_eff_blackbody_slope": 1/patmosx_sat_coeffs["a2_{0}".format(ch)]
            }
        for t in range(1, 5):
            pygac_sat_coeffs["thermometer_{0}".format(t)] = {
                "d{0}".format(d): float(patmosx_sat_coeffs["PRT{0}_{1}".format(t, d)])
                for d in range(5)
            }
        return pygac_sat_coeffs

    @staticmethod
    def find_s0(s0_low, s0_high, ch):
        """Find the single-gain calibration slope at launch date.

        Arguments
            s0_low - low gain calibration slope at launch date
            s0_high - high gain calibration slope at launch date
            ch - channel name ("1", "2", "3a")

        Note:
            In case of a single-gain instrument, s0_low is equal to s0_high.
        """
        if s0_low == s0_high:
            # single gain case
            return s0_low
        if ch == "3a":
            g_low, g_high = 0.25, 1.75
        else:
            g_low, g_high = 0.5, 1.5

        # Note: the PATMOS-x coefficients are rounded to three decimals.
        def diff(s0): return s0_low - round(s0*g_low, 3) + s0_high - round(s0*g_high, 3)

        s0_l = s0_low / g_low
        s0_h = s0_high / g_high
        if diff(s0_l) == 0:
            s0 = s0_l
        elif diff(s0_h) == 0:
            s0 = s0_h
        else:
            s0 = bisect(diff, s0_l, s0_h)
        return s0

    @staticmethod
    def float2date(date_float):
        """Convert date float to date.

        Argument
            date_float (float) - date as year

        Return
            date (datetime.datetime) - date

        Note
            This is the reverse function of date2float.
        """
        year = int(date_float)
        days_in_year = (dt.datetime(year+1, 1, 1) - dt.datetime(year, 1, 1)).days
        seconds = date_float*days_in_year*24*3600 - year*days_in_year*24*3600
        diff = dt.timedelta(seconds=seconds)
        date = dt.datetime(year, 1, 1) + diff
        return date

    def save(self, filepath):
        """Save coefficients as PyGAC json file."""
        with open(filepath, mode="w") as json_file:
            json.dump(self.coeffs, json_file, indent=4, sort_keys=True)


def main():
    """The main function."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tarball", type=str, help="path to PATMOS-x coefficients tarball")
    parser.add_argument("-o", "--output", type=str, metavar="JSON",
                        help='path to PyGAC json file, defaults to tarball path with suffix ".json"')
    parser.add_argument("-v", "--verbose", action="store_true", help="explain what is being done")
    args = parser.parse_args()
    if args.verbose:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING
    logging.basicConfig(level=loglevel, format="[%(asctime)s] %(message)s")
    tarball = pathlib.Path(args.tarball)
    logging.info('Read PATMOS-x tarball "%s".', tarball)
    patmosx_coeffs = PatmosxReader(tarball)
    logging.info("Translate PATMOS-x coefficients to PyGAC format.")
    pygac_coeffs = Translator(patmosx_coeffs)
    output = args.output or tarball.with_suffix(".json")
    logging.info('Write PyGAC calibration json file "%s".', output)
    pygac_coeffs.save(output)
    logging.info("Done!")

if __name__ == "__main__":
    main()
