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
        self.file.Gain_Calibrated = False
        """calibrate the continuum data."""
        # print(file.data["DATA"])
        f = Cal(file)
        #go through each data channel and calibrate the heights
        for ind, i in enumerate(file.data):
            
            calib_height_data = []
            #first check if the gain start and end values are present
            if file.gain_start[ind] is not None and file.gain_end[ind] is not None:
                #get all the gain values
                delta1 = file.gain_start[ind][0]
                delta2 = file.gain_end[ind][0]
                time1 = file.gain_start[ind][1]
                time2 = file.gain_end[ind][1]

                #get an array of the continuum data
                data = file.data[ind]
                data = f.sdfits_to_array(file, data)
                
                #for the time array in the data find the calibrated height for each intensity
                for i in range(len(data[0])):
                    delta = delta1 + (delta2 - delta1) * (data[0][i] - time1) / (time2 - time1)
                    calib_height = data[1][i]/delta
                    calib_height_data.append(calib_height)
                #print("HEre", calib_height_data)

            elif file.gain_start[ind] is None and file.gain_end[ind] is not None:
                delta = file.gain_end[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(file, data)


                for i in range(len(data[0])):
                    calib_height = data[1][i]/delta
                    calib_height_data.append(calib_height)

            elif file.gain_start[ind] is not None and file.gain_end[ind] is None:
                delta = file.gain_start[ind][0]

                data = file.data[ind]
                data = f.sdfits_to_array(file, data)


                for i in range(len(data[0])):
                    calib_height = data[1][i]/delta
                    calib_height_data.append(calib_height)

            elif file.gain_start[ind] is None and file.gain_end[ind] is None:
                self.file.Gain_Calibrated = True
                continue

            file.data[ind] = calib_height_data
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
    s.clean_sections()

    c = Cal(file)
    c.gain_calibration(file)
    Data.calib_Heights(file)

    # print(file.data[0])


    # Further processing can be added here