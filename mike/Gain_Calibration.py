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


    def calib_Heights(self, file):
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
                data = f.sdfits_to_array(file, data)

                #for the time array in the data find the calibrated height for each intensity
                for idx, i in enumerate(data[0]):
                    delta = delta1 + (delta2 - delta1) * (data[0][idx] - time1) / (time2 - time1)
                    calib_height = data[1][idx]/delta
                    calib_height_data.append([data[0][idx], calib_height])
                #print("HEre", calib_height_data)

            # If gain_start is None and gain_end is not None, use gain_end for calibration
            elif file.gain_start[ind] is None and file.gain_end[ind] is not None:
                calibrations += 1
                delta = file.gain_end[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(file, data)


                for idx, i in enumerate(data[0]):
                    calib_height = data[1][idx]/delta
                    calib_height_data.append([data[0][idx], calib_height])

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif file.gain_start[ind] is not None and file.gain_end[ind] is None:
                calibrations += 1
                delta = file.gain_start[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(file, data)


                for idx, i in enumerate(data[0]):
                    calib_height = data[1][idx]/delta
                    calib_height_data.append([data[0][idx], calib_height])

            # If both gain_start and gain_end are None, don't calibrate
            elif file.gain_start[ind] is None and file.gain_end[ind] is None:
                self.file.Gain_Calibrated = False

                data = file.data[ind]
                data = f.sdfits_to_array(file, data)

                # If both gain_start and gain_end are None, just copy the original data
                #But gotta turn it into a list of time, intensity lists
                calib_height_data = []
                for t, intensity in zip(data[0], data[1]):
                    calib_height_data.append([t, intensity])
                
                self.file.continum.append(calib_height_data)
                continue

            # add the calibrated height data to continum
            self.file.continum.append(calib_height_data)
            #check if all calibrations are done
            if calibrations == datas:
                self.file.Gain_Calibrated = True

if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    file = Mike("TrackingHighRes/0136376.fits")
    Data = gain_calibration(file)

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()

    c = Cal(file)
    c.gain_calibration(file)
    Data.calib_Heights(file)

    # print(file.continum)