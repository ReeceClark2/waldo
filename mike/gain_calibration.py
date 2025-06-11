from calendar import day_abbr
from astropy.io import fits
from astropy.table import Table
from datetime import datetime
import numpy as np
from file_exception import MyException
import warnings
from file_init import Mike
from astropy.table import Column
import os
import time
from cal import Cal
from val import Val
from sort import Sort
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from child_init import Sully



class gain_calibration:
    """
    Class to handle continuum data from FITS files and perform gain calibration.
    Attributes:
        file: An object representing the FITS file, expected to have a 'file_path' attribute and other calibration-related attributes.
        dataH: The header of the FITS file's primary data extension.
        data: The data table extracted from the FITS file.
    Methods:
        __init__(file):
            Initializes the gain_calibration object by loading the FITS file header and data table.
        calib_Heights(file):
            Calibrates the continuum data using gain start and end values.
            For each data entry, applies gain calibration based on the presence of gain start and/or gain end values.
            Updates the data in-place and sets the 'Gain_Calibrated' flag accordingly.
    """
    """Class to handle continuum data from FITS files."""
    def __init__(self, file):
        """
        Initialize the Gain_Calibration object with a FITS file.
            """
        self.file = file
        with fits.open(self.file.file_path) as hdul:
            self.dataH = hdul[1].header
            self.data = Table(hdul[1].data)


    def calib_Heights(self, file, feeds=None):
        """
        calibrate the continuum data.

        param file: Mike class file

        returns: adds calibrated height data to the file's continuum
        """
        # print(file.data["DATA"])
        f = Cal(file)
        calibrations = 0
        datas = len(file.data)
        #go through each data channel and calibrate the heights
        for ind, i in enumerate(file.data):
            feednum = np.unique(self.file.data[ind]["IFNUM"])[0]
            pol = np.unique(self.file.data[ind]["PLNUM"])[0]
            if feeds is not None and feednum not in feeds:
                continue  

            calib_height_data = []
            #first check if the gain start and end values are present
            if file.gain_start[ind] is not None and file.gain_end[ind] is not None:
                calibrations += 1
                #get all the gain values
                delta1 = file.gain_start[ind][0]
                delta2 = file.gain_end[ind][0]
                time1 = file.gain_start[ind][1]
                time2 = file.gain_end[ind][1]

                #get an array of the continuum data
                data = file.data[ind]
                data = f.sdfits_to_array(data)

                #for the time array in the data find the calibrated height for each intensity
                calib_height = []
                for idx, i in enumerate(data[0]):
                    delta = delta1 + (delta2 - delta1) * (data[0][idx] - time1) / (time2 - time1)
                    #add each calibration height to the calib_height list
                    data[1][idx] = (data[1][idx]/delta)
                calib_height_data = (data)

            # If gain_start is None and gain_end is not None, use gain_end for calibration
            elif file.gain_start[ind] is None and file.gain_end[ind] is not None:
                calibrations += 1
                delta = file.gain_end[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(data)

                temp = data[1] / delta
                calib_height_data = ([data[0], temp])

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif file.gain_start[ind] is not None and file.gain_end[ind] is None:
                calibrations += 1
                delta = file.gain_start[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(data)

                data[1] = data[1] / delta
                calib_height_data = (data)

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif file.gain_start[ind] is not None and file.gain_end[ind] is None:
                calibrations += 1
                delta = file.gain_start[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(data)

                data[1] = data[1] / delta
                calib_height_data.append(data)
                # If both gain_start and gain_end are None, just copy the original data
                #But gotta turn it into a list of time, intensity lists
                
                self.file.continuum.append(calib_height_data)
                continue

            # add the calibrated height data to continuum
            self.file.continuum.append(calib_height_data)
            #check if all calibrations are done
            if calibrations == datas:
                self.file.Gain_Calibrated = True
            # Plot the spectrum
            plt.plot(calib_height_data[0], calib_height_data[1], label=f"{pol} {feednum+1}")
            # Check if there are any more entries with the same feednum left
            feeds_left = any(np.unique(self.file.data[j]["IFNUM"]) == feednum for j in range(ind + 1, len(self.file.data)))
            if not feeds_left:
                plt.title(f"Feed {feednum+1} Continuum Spectrum")
                plt.xlabel("Times (sec)")
                plt.ylabel("Intensity")
                plt.grid()
                plt.legend()
                plt.savefig(f"Continuum_{feednum+1}.png")
                plt.close()
                print(f"Saved continuum plot for Feed {feednum+1}")

if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    keep_times = [[8,12], [1404, 1410], [1412, 1420]]  # Specify the indices you want to keep
    feed= [0]  # Specify the feeds you want to keep
    file = Mike("C:\\Users\\anshm\\Downloads\\0126929.fits")
    Data = gain_calibration(file)

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.get_startend_freqs()
    s.get_startstop_channels()
    

    c = Cal(file)
    c.gain_calibration()

    g = gain_calibration(file)
    g.calib_Heights(file)

    if keep_times != []:
        child = Sully(file)
        sortchild = Sort(child)
        contChild = gain_calibration(child)

        sortchild.user_cuts(keep_times, "continuum", "cut", feed)
        g.calib_Heights(file, feed)
        # s.back_to_original_data(org_data)
        print(file.data[0]["DATA"].shape)
        print(child.data[0]["DATA"].shape)
