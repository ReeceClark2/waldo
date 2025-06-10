from contourpy import SerialContourGenerator
import numpy as np
from scipy.stats import linregress
from datetime import datetime
import rcr

from file_exception import MyException
from file_init import Mike
from val import Val
from sort import Sort
from cal import Cal

import matplotlib
import re
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Spec:
    def __init__(self, file):
        '''
        Initialize file.
        '''

        self.file = file


    def average(self, data, axis):
        '''
        Integrate across time (axis = 0) or integrate across frequency (axis = 1).

        params: data: astropy FITS table

        returns: array of summed intensities for each frequency
        '''

        intensities = data['DATA']
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means


    def make_spec(self):
        '''
        Generate a spectrum plot for each channel in the data.

        params: file: Mike class file

        returns: populates the file's spectrum field with frequency and intensity data
        '''
        #go through each channel in the data
        for i in range(len(self.file.data)):
            #get feed and polarization numbers
            feednum = np.unique(self.file.data[i]["IFNUM"])[0]
            polnum = np.unique(self.file.data[i]["PLNUM"])[0]   
            # print(f"Processing data for IFNUM: {feednum}, PLNUM: {polnum}")

            # only get actual data, not the calibration data
            data = self.file.data[i][self.file.data_indices[i][0]:self.file.data_indices[i][1]]

            # get the summed intensities for each frequency
            result = self.average(data, 0)
            #reverse the list becaise the data is in reverse order
            result = result[::-1]

            #get the frequency range for the channel based on the feed
            frequencies = np.linspace(self.file.freqs[feednum][0], self.file.freqs[feednum][1], len(result))

            # If the feednum is 0, then it's XX polarization, if 1 then YY, otherwise
            if polnum == 0:
                pol = "XX"
            elif polnum == 1:
                pol = "YY"
            #placeholder for other polarizations
            else:
                pol = f"Pol{feednum+1}"

            # Add the frequency and result to the file's spectrum field
            self.file.continuum.append(np.array([np.array(frequencies), np.array(result)]))


            # Plot the spectrum
            plt.plot(frequencies, result, label=f"{self.file.labels[i]}")
            # Check if there are any more entries with the same feednum left
            feeds_left = any(np.unique(self.file.data[j]["IFNUM"]) == feednum for j in range(i + 1, len(self.file.data)))
            if not feeds_left:
                plt.title(f"Feed {feednum+1} Spectrum")
                plt.xlabel("Frequency (MHz)")
                plt.ylabel("Intensity")
                plt.grid()
                plt.legend()
                plt.savefig(f"Spectrum_feed_Num_{feednum+1}.png")
                plt.close()
                print(f"Saved spectrum plot for Feed {feednum+1} with polarization {pol}.")

if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''

    file = Mike("TrackingHighRes/0136376.fits")
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.get_startend_freqs()
    s.get_startstop_channels()

    c = Cal(file)
    c.compute_gain_deltas()

    spec = Spec(file)
    spec.make_spec()

    # print(file.continuum[0])
