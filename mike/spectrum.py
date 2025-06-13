import numpy as np
from datetime import datetime
from file_init import Mike
from val import Val
from sort import Sort
from cal import Cal
import matplotlib
import os
from glob import glob
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from child_init import Sully


class Spectrum:
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


    def make_spec(self, indices=None, feeds=None):
        '''
        Generate a spectrum plot for each channel in the data.

        params: file: Mike class file

        returns: populates the file's spectrum field with frequency and intensity data
        '''
        #go through each channel in the data
        for ind, i in enumerate(self.file.data):
            #get feed and polarization numbers
            feednum = np.unique(self.file.data[ind]["IFNUM"])[0]
            polnum = np.unique(self.file.data[ind]["PLNUM"])[0]
            if feeds is not None and feednum not in feeds:
                continue   

            # only get actual data, not the calibration data
            data = self.file.data[ind][self.file.data_indicies[ind][0]:self.file.data_indicies[ind][1]]

            # get the summed intensities for each frequency
            result = self.average(data, 0)
            #reverse the list becaise the data is in reverse order
            result = result[::-1]

            #get the frequency range for the channel based on the feed
            # if the spectrum field has not been fully populated yet create the frequency range from the original start and stop frequencies
            if len(self.file.spectrum) == ind:
                frequencies = np.linspace(self.file.freqs[ind][0], self.file.freqs[ind][1], len(result))
            else:
                frequencies = self.file.spectrum[ind][0]

            # Add the frequency and result to the file's spectrum field
            self.file.spectrum.append([np.array(frequencies), np.array(result)])

        return


if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''
    keep_indices = [[1300, 1400], [1406, 1410], [1412, 1420]]  # Specify the indices you want to keep
    feed= [1]  # Specify the feeds you want to keep
    fits_path = "TrackingHighRes/0136484.fits"
    file = Mike(fits_path)
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.get_startend_freqs()
    s.get_startstop_channels()

    spec = Spectrum(file)
    spec.make_spec()

    if keep_indices != []:
        child = Sully(file)
        sortchild = Sort(child)
        specchild = Spectrum(child)

        sortchild.user_cuts(keep_indices, "spectrum", "cut", feed)
        specchild.make_spec(keep_indices, feed)
        print(file.data[0]["DATA"].shape)
        print(child.data[0]["DATA"].shape)

    c = Cal(file)
    c.compute_gain_deltas()


 
