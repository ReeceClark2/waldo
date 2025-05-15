from astropy.io import fits
from scipy.special import erfc
from scipy.ndimage import label
from scipy.stats import linregress
import numpy as np
from file_exception import MyException
from datetime import datetime
import rcr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from file_init import Mike
from val import Val


class Cal:
    def __init__(self, file):
        '''Initialize file and confirm header and data have been validated.'''

        self.file = file

        if self.file.validated_header == False:
            raise MyException('FITS header has not been validated!')
        # if self.file.validated_data == False:
        #     raise MyException('FITS data has not been validated!')


    def split_calibration(self):
        """
        Split data into alternating calibration and non-calibration sets.
        Each contiguous calibration (CALSTATE == ±1) or science (CALSTATE == 0) block becomes a separate element in the returned lists.
        Returns:
            non_cal_list: List of non-calibration DataFrames
            cal_list: List of calibration DataFrames
        """
        calstate = self.file.data["CALSTATE"]
        cal_mask = (calstate == 1) | (calstate == -1)

        labels, num_labels = label(cal_mask)

        chunks = []
        cal_list = []
        non_cal_list = []

        start = 0
        for i in range(1, len(labels)):
            if labels[i] != labels[i - 1]:
                end = i
                chunk = file.data[start:end]
                if labels[start] == 0:
                    non_cal_list.append(chunk)
                else:
                    cal_list.append(chunk)
                start = i

        # Handle final segment
        chunk = file.data[start:]
        if labels[start] == 0:
            non_cal_list.append(chunk)
        else:
            cal_list.append(chunk)

        print(len(cal_list[0]), len(cal_list[1]))
        return non_cal_list, cal_list
    

    def split_polarity(self, data):
        """Split data into X and Y polarization channels based on PLNUM (0 = X, 1 = Y)."""

        x_pol = data[data["PLNUM"] == 0]
        y_pol = data[data["PLNUM"] == 1]

        return x_pol, y_pol
    

    def split_slp(self, data):
        """Dynamically split data into channels based on IFNUM (e.g., pol0, pol1, ...)."""

        split = {}
        ifnums = np.unique(data["IFNUM"])

        for ifnum in ifnums:
            split[f"pol{ifnum}"] = data[data["IFNUM"] == ifnum]

        return split
    

    def integrate(self, data, axis):
        """Integrate across time (axis = 0) or integrate across frequency (axis = 1)."""

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means
    

    def linear(self, x, params): # model function
        return params[0] + x * params[1]


    def d_linear_1(self, x, params): # first model parameter derivative
        return 1


    def d_linear_2(self, x, params): # second model parameter derivative
        return x


    def rcr(self, x, y):
        """Perform Robust Chauvenet Rejection on a given dataset."""
        x = x[1:]
        y = y[1:]

        result = linregress(x, y)
        m = result.slope
        b = result.intercept

        guess = [m, b]
        model = rcr.FunctionalForm(self.linear,
            x,
            y,
            [self.d_linear_1, self.d_linear_2],
            guess
        )

        r = rcr.RCR(rcr.LS_MODE_68) 
        r.setParametricModel(model)
        r.performBulkRejection(y)

        best_fit_parameters = model.result.parameters

        plt.scatter(x, y, color='black')
        x = np.linspace(x[0], x[-1], 1000)
        y = m * x + b
        plt.plot(x, y, label='Pre RCR', color='red')

        print('Best fit parameters:', best_fit_parameters)
        b = best_fit_parameters[0]
        m = best_fit_parameters[1]
        y = m * x + b
        plt.plot(x, y, label='Post RCR', color='green')

        plt.legend()
        plt.savefig('RCR Fit')

        return


    def sdfits_to_array(self):
        """Convert SDFITS data to arrays to be used by RCR."""

        data, cal = self.split_calibration()

        x_pol, y_pol = self.split_polarity(cal[0])
        xxs = self.split_slp(x_pol)
        yys = self.split_slp(y_pol)

        xxs_cal = []
        yys_cal = []

        for i, key in enumerate(xxs):
            xx = xxs[key]
            xx_freq = self.integrate(xx, axis=1)

            times_xx = [datetime.fromisoformat(t) for t in xx["DATE-OBS"]]
            t0 = times_xx[0]
            time_xx = [(t - t0).total_seconds() for t in times_xx]

            xx_freq = np.array(xx_freq)
            time_xx = np.array(time_xx)

            xxs_cal.append([xx_freq, time_xx])

        for i, key in enumerate(yys):
            yy = yys[key]
            yy_freq = self.integrate(yy, axis=1)

            times_yy = [datetime.fromisoformat(t) for t in yy["DATE-OBS"]]
            t0 = times_yy[0]
            time_yy = [(t - t0).total_seconds() for t in times_yy]

            yy_freq = np.array(yy_freq)
            time_yy = np.array(time_yy)

            yys_cal.append([yy_freq, time_yy])

        return xxs_cal, yys_cal
    

    def clean_data(self):
        """Clean data of outliers via RCR post extracting values from SDFITS file. Remove outliers from original data."""
        
        xxs_cal, yys_cal = self.sdfits_to_array()

        new_xxs = []
        new_yys = []
        
        plt.scatter(xxs_cal[0][0][1:],xxs_cal[0][1][1:])
        plt.savefig('test cal')

        values = self.rcr(xxs_cal[0][0], xxs_cal[0][1])
        new_xxs.append(values[0])

        print(len(xxs_cal[0]), len(new_xxs[0]))
        

if __name__ == "__main__":
    '''Test function to implement calibration.'''

    file = Mike("C:/Users/starb/Downloads/0136645.fits")
    v = Val(file)
    v.validate_primary_header()
    # v.validate_data()

    c = Cal(file)
    c.clean_data()
