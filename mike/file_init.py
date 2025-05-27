from astropy.io import fits
import numpy as np
from file_exception import MyException
import warnings

class Mike:
    def __init__(self, file_path):
        '''Initialize binary fits cube by reading its header and data from a given file.'''

        with fits.open(file_path) as hdul:
            self.header = hdul[0].header
            self.data = hdul[1].data
        self.missing_values = []
        self.file_path = file_path

        self.validated_header = False
        self.validated_data = False

        self.data_indices = []
        self.gain_start = []
        self.gain_end = []



if __name__ == "__main__":
    '''Test function to test uploaded file.'''

    file = Mike("C:/Users/starb/Downloads/0115701.fits")

    print(file.header)
    np.set_printoptions(threshold=100000)
    print(len(file.data['SWPVALID']))
    # print(file.data['CALSTATE'])