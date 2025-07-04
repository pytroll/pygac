import numpy as np
from datetime import datetime, timedelta
from pygac import get_reader_class
import argparse
from pygac.calibration.vis_uncertainty import get_random, vis_uncertainty

# patmos calibration coefficients for noaa19
ch1_S0_low = 0.054
ch1_S0_high = 0.0163
ch1_S1_low = 0.286
ch1_S1_high = 0.286
ch1_S2_low = 0.012
ch1_S2_high = 0.012
ch2_S0_low = 0.061
ch2_S0_high = 0.183
ch2_S1_low = 0.478
ch2_S1_high = 0.478
ch2_S2_low = 0.052
ch2_S2_high = 0.052
ch3_S0_low = 0.0272
ch3_S0_high = 0.188
ch3_S1_low = 0.0
ch3_S1_high = 0.0
ch3_S2_low = 0.0
ch3_S2_high = 0.0

date_of_launch = "2009-02-05T00:57:36.000000Z"
min_year = 2009.167
max_year = 2023.583

def get_patmos_coeffs(time):
    #
    # Get time in year fraction and then get delta_time
    #
    year = time.year
    start_year = datetime(year, 1, 1)
    end_year = datetime(year, 12, 31)
    delta_y = (time - start_year).total_seconds() * 1 / \
              (end_year - start_year).total_seconds()
    time_frac = year + delta_y
    time_diff = time_frac - date_of_launch
    time_diff_m = time_diff * 365
    print(time, time_frac, time_diff)


    slope_low_1 = ch1_S0_low * (100. + ch1_S1_low * time_diff + \
                                   ch1_S2_low * time_diff * time_diff) / \
                100.
    slope_high_1 = ch1_S0_high * (100. + ch1_S1_high * time_diff + \
                                     ch1_S2_high * time_diff * time_diff) / \
                 100.
    slope_1 = (slope_high_1 + slope_low_1) / 2

    slope_low_2 = ch2_S0_low * (100. + ch2_S1_low * time_diff + \
                                   ch2_S2_low * time_diff * time_diff) / \
                100.
    slope_high_2 = ch2_S0_high * (100. + ch1_S1_high * time_diff + \
                                     ch2_S2_high * time_diff * time_diff) / \
                 100.
    slope_2 = (slope_high_2 + slope_low_2) / 2

    slope_low_3 = ch3_S0_low * (100. + ch3_S1_low * time_diff + \
                                   ch3_S2_low * time_diff * time_diff) / \
                100.
    slope_high_3 = ch3_S0_high * (100. + ch3_S1_high * time_diff + \
                                     ch3_S2_high * time_diff * time_diff) / \
                 100.
    slope_3 = (slope_high_3 + slope_low_3) / 2

    return slope_1, slope_2, slope_3


# retrieve gain slopes from vis_uncertainty.py and plot against patmos gains
def plot_gain(ds, mask, plot=False):

    vis_uncert = vis_uncertainty(ds, mask, plot=False)
    gain_ch1 = vis_uncert.gain_1
    gain_ch2 = vis_uncert.gain_2
    gain_ch3 = vis_uncert.gain_3

    outrad = None

    ndata = (int(max_year) + 1 - int(date_of_launch)) * 4
    yeararr = int(date_of_launch) + np.arange(ndata) * 0.25
    gd = (yeararr > date_of_launch) & (yeararr <= max_year)
    yeararr = yeararr[gd]
    year_arr = np.array([date_of_launch])
    year_arr = np.append(year_arr, yeararr)

    c1 = []
    c2 = []
    c3 = []
    time_arr = []
    for yr in year_arr:
        iyear = int(yr)
        remainder = yr - iyear
        boy = datetime(iyear, 1, 1)
        eoy = datetime(iyear + 1, 1, 1)
        seconds = remainder * (eoy - boy).total_seconds()
        dtime = boy + timedelta(seconds=seconds)
        time_arr.append(dtime)

        c1.append(get_patmos_coeffs(dtime)[0])
        c2.append(get_patmos_coeffs(dtime)[1])
        c3.append(get_patmos_coeffs(dtime)[2])

    outrad = np.zeros((2, 3, len(c1)))
    outrad[0, 0, :] = np.array(c1)
    outrad[0, 1, :] = np.array(c2)
    outrad[0, 2, :] = np.array(c3)

    # define pygac gains here

    outrad[1, 0, :] = np.array(gain_ch1)
    outrad[1, 1, :] = np.array(gain_ch2)
    outrad[1, 2, :] = np.array(gain_ch3)


    if plot:
        import matplotlib.pyplot as plt
        cols = ['blue', 'orange', 'red', 'green', 'magenta', 'pink']
        labels = ['patmos','pygac']
        time_arr = np.array(time_arr)

        for j in range(outrad.shape[0]):
            plt.plot(time_arr, (outrad[j, 1, :]), label=labels[j])
        plt.title('NOAA-19')
        ax = plt.gca()
        plt.ylabel('Calibration Slope')
        plt.legend()
        plt.show()

    return gain_ch1, gain_ch2, gain_ch3

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()

    reader_cls = get_reader_class(args.filename)

    reader = reader_cls(tle_dir="/gws/nopw/j04/npl_eo/users/nyaghnam/pygac/TLE",
                        tle_name="TLE_%(satname).txt",
                        calibration_method="noaa",
                        adjust_clock_drift=False)
    reader.read(args.filename)
    ds = reader.get_calibrated_dataset()
    mask = reader.mask

    args = parser.parse_args()
    cal_slopes = plot_gain(ds,mask,plot=True)
    print(cal_slopes)






