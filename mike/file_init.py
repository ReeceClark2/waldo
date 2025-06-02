from astropy.io import fits
from file_exception import MyException
import numpy as np


class Mike:
    def __init__(self, file_path):
        '''
        Initialize binary fits cube by reading its header and data from a given file.
        
        File fields:
        missing values: empty primary header values
        file_path: path to file in directory

        validated_header: check whether header has been validated
        validated_data: check whether data has been validated

        data_indices: list of arrays of shape (1,2)
            1st value is the first data point in the corresponding indexed data array
            2nd value is the second data point in the corresponding indexed data array
        gain_start: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        gain_end: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        '''

        try:
            with fits.open(file_path) as hdul:
                self.header = hdul[0].header
                self.data = hdul[1].data
        except Exception as e:
            raise MyException(f"Error reading FITS file: {e}")
        
        self.missing_values = []
        self.file_path = file_path

        self.validated_header = False
        self.validated_data = False
        self.Gain_Calibrated = False

        self.data_indices = []
        self.gain_start = []
        self.gain_end = []

        self.continum = []
        self.spectrum = [[]]


if __name__ == "__main__":
    '''
    Test function to test uploaded file.
    '''

    file = Mike("C:/Users/starb/Downloads/0136940.fits")

    np.set_printoptions(threshold=100000)

    print(repr(file.header))
    print(repr(file.data[0]))