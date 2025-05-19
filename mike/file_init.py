from astropy.io import fits
import numpy as np
from file_exception import MyException
import warnings

class Mike:
    def __init__(self, file_path):
        '''Initialize binary fits cube by reading its header and data from a given file.'''

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


if __name__ == "__main__":
    file = Mike("0136645V2.fits")
    
    print(file.header)
    print(file.data[0])