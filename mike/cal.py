from astropy.io import fits
from scipy.special import erfc
import numpy as np
from file_exception import MyException
import warnings
from mike import Mike
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

class Cal:
    def __init__(self, file):
        self.cal = None
        self.no_cal = None
        self.file = file

        if self.file.validated_header == False:
            raise MyException('FITS header has not been validated!')
        # if self.file.validated_data == False:
        #     raise MyException('FITS data has not been validated!')


    def split_calibration(self):
        """Split data into calibration and non-calibration sets based on CAL (-1, 1 = cal, 0 = no cal)."""

        non_cal = file.data[file.data["CALSTATE"] == 0]
        cal = file.data[(file.data["CALSTATE"] == 1) | (file.data["CALSTATE"] == -1)]

        return non_cal, cal
    

    def split_polarity(self, data):
        """Split data into X and Y polarization channels based on PLNUM (0 = X, 1 = Y)."""

        x_pol = data[data["PLNUM"] == 0]
        y_pol = data[data["PLNUM"] == 1]

        return x_pol, y_pol
    

    def split_slp(self, data):
        """Dynamically split data into channels based on IFNUM (e.g., pol0, pol1, ...)."""

        split = {}
        ifnums = np.unique(data["IFNUM"])
        print(ifnums)
        for ifnum in ifnums:
            split[f"pol{ifnum}"] = data[data["IFNUM"] == ifnum]

        return split
    

    def integrate(self, data, axis):
        """Integrate across time (axis = 0) or integrate across frequency (axis = 1)."""

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means
    

    def chauvenet_rejection(self):
        _, cal = self.split_calibration()

        off = file.data[file.data["CALSTATE"] == 1]
        on = file.data[file.data["CALSTATE"] == -1]
        print(off)

        x_pol, y_pol = self.split_polarity(cal)
        xxs = self.split_slp(x_pol)

        for i, key in enumerate(xxs):
            xx = xxs[key]
            xx_freq = self.integrate(xx, axis=1)

            times_xx = [datetime.fromisoformat(t) for t in xx["DATE-OBS"]]
            t0 = times_xx[0]
            time_xx = [(t - t0).total_seconds() for t in times_xx]

            xx_freq = np.array(xx_freq)
            time_xx = np.array(time_xx)

            N = len(xx_freq)
            mean = np.mean(xx_freq)
            std = np.std(xx_freq)
            criterion = 1.0 / (2 * N)

            deviations = np.abs(xx_freq - mean) / std

            prob = erfc(deviations / np.sqrt(2))

            mask = prob >= criterion
            filtered_freq = xx_freq[mask]
            filtered_time = time_xx[mask]

            plt.figure()
            plt.scatter(time_xx, xx_freq, label="Original", alpha=0.4)
            plt.scatter(filtered_time, filtered_freq, color='red', label="Accepted")
            plt.legend()
            plt.xlabel("Time [s]")
            plt.ylabel("Integrated Signal")
            plt.title(f"Chauvenet Rejection on {key}")
            plt.savefig(f"chauvenet_filtered_{key}.png")





file = Mike("C:/Users/starb/Downloads/0136645.fits")
file.validate_primary_header()
file2 = Cal(file)
file2.chauvenet_rejection()