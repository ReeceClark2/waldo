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



class Gain_Calibration:
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

        f = Cal(file)
        calibrations = 0
        datas = len(file.data)

        # Go through each data channel and calibrate the heights
        for ind1, i in enumerate(file.data):
            calib_height_data = []
            # First check if the gain start and end values are present
            if file.gain_start[ind1] is not None and file.gain_end[ind1] is not None:
                calibrations += 1
                # Get all the gain values
                delta1 = file.gain_start[ind1][0]
                delta2 = file.gain_end[ind1][0]
                time1 = file.gain_start[ind1][1]
                time2 = file.gain_end[ind1][1]

                # Get an array of the continuum data
                data = file.data[ind1]
                data = f.sdfits_to_array(data)

                # For the time array in the data find the calibrated height for each intensity
                for ind2, j in enumerate(data[0]):
                    delta = delta1 + (delta2 - delta1) * (data[0][ind2] - time1) / (time2 - time1)
                    # Add each calibration height to the calib_height list
                    data[1][ind2] = (data[1][ind2] / delta)

                calib_height_data = data

            # If gain_start is None and gain_end is not None, use gain_end for calibration
            elif file.gain_start[ind1] is None and file.gain_end[ind1] is not None:
                calibrations += 1
                delta = file.gain_end[ind1][0]

                data = file.data[ind1]
                data = f.sdfits_to_array(data)

                data[1] = data[1] / delta
                calib_height_data = (data)

            # If gain_start is present but gain_end is not, use gain_start for calibration
            elif file.gain_start[ind1] is not None and file.gain_end[ind1] is None:
                calibrations += 1
                delta = file.gain_start[ind1][0]

                data = file.data[ind1]
                data = f.sdfits_to_array(data)

                data[1] = data[1] / delta
                calib_height_data = (data)

            # If gain_start is not present and gain_end is not pass default data for continuum
            elif file.gain_start[ind1] is None and file.gain_end[ind1] is None:
                data = file.data[ind1]
                data = f.sdfits_to_array(data)
                
                self.file.continuum.append(data)
                continue

            # Add the calibrated height data to continuum
            self.file.continuum.append(calib_height_data)

            # Check if all calibrations are done
            if calibrations == datas:
                self.file.Gain_Calibrated = True


if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    file = Mike("TrackingHighRes/0136376.fits")
    Data = Gain_Calibration(file)

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()

    c = Cal(file)
    c.compute_gain_deltas()
    Data.calib_Heights(file)

    print(file.continuum[1][0][1])
