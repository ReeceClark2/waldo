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
        '''

        intensities = data
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means

    

    def make_spec(self):
        '''        
        Generate a spectrum plot for each channel in the data.

        params: file: Mike class file

        returns: saves a spectrum plot for each channel in the data
        '''
        # Generate a spectrum plot for each channel in the data.
        center = []
        for i in range(len(self.file.data)):
            #get feed and polarization numbers
            feednum = np.unique(self.file.data[i]["IFNUM"])
            polnum = np.unique(self.file.data[i]["PLNUM"])
            print(f"Processing data for IFNUM: {feednum}, PLNUM: {polnum}")

            # only get actual data, not the calibration data
            data = self.file.data[i][self.file.data_indices[i][0]:self.file.data_indices[i][1]]

            #print(repr(self.file.header))
            #scour the header for the bandwidth, center frequencies, and start/stop channels
            for key, value in self.file.header.items():
                #bandwidth
                if key == ("OBSBW"):
                    band = value
                #lowres center frequency
                if key == ("OBSFREQ"):
                    center = [value]
                #center frequencies and start/stop channels
                if key == ("HISTORY"):
                    #if HIRES BANDS exist replace center with the HIRES center frequencies
                    if value.startswith("HIRES bands"):
                        # print(f"Found HISTORY key: {key} with value: {value}")
                        # Extract all integers from the string
                        value = value.replace(",", " ").strip()
                        value = value.split(" ")
                        center = []
                        #split the value into individual words and numbers
                        for k in value:
                            k = str(k).strip()
                            try:
                                #if it's a float add it to the center list
                                k = float(k)
                                center.append(k)
                            except ValueError:
                                #otherwise skip it
                                continue
                    # start/stop channels
                    elif value.startswith("START,STOP"):
                        # Extract all integers from the string
                        value = value.replace(",", " ").strip()
                        value = value.split(" ")
                        #split the value into individual words and numbers
                        channels = []
                        for k in value:
                            k = str(k).strip()
                            try:
                                #if it's an integer add it to the channels list
                                k = int(k)
                                channels.append(k)
                            except ValueError:
                                #if it can't become a integer, skip it
                                continue
                        # Remove all string type characters from the channels list
                        channel1 = int(channels[0])
                        channel2 = int(channels[1])
                        #only keep the correct frequencies
                        freqs = np.array([arr[channel1:channel2] for arr in data["DATA"]])
            #sum the intensities across frequencies
            result = self.average(freqs, 0)
            #reverse the list becaise the data is in reverse order
            result = result[::-1]

            # Generate An x axis as long as the intensities per frequency array 
            Xs = np.arange(len(result))
            # Scale the X values to the bandwidth
            Xs = Xs*band/len(Xs)
            # Add the center frequency minus half the bandwidth to the X values
            Xs = Xs + (center[int(feednum)]-band/2)  # Add the center frequency to the X values

            # If the feednum is 0, then it's XX polarization, if 1 then YY, otherwise
            if polnum == 0:
                pol = "XX"
            elif polnum == 1:
                pol = "YY"
            #placeholder for other polarizations
            else:
                pol = f"Pol{feednum+1}"

            self.file.continuum.append(np.array([np.array(Xs), np.array(result)]))


            # Plot the spectrum
            plt.plot(Xs, result, label=f"{pol} {feednum+1}")
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

    c = Cal(file)
    c.gain_calibration()

    spec = Spec(file)
    spec.make_spec()

    # print(file.continuum[0])
