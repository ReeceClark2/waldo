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


    def validate_primary_header(self):
        '''Validate the primary header of the given FITS file by:
            1) checking if it complies to 2880 byte header standard
            2) checking primary header cards
            3) checking for duplicate cards
            4) checking that each card is 80 character long
            5) recording missing values in the header'''
        
        self.validate_header_size()
        self.validate_primary_header_cards()
        self.validate_header_cards()

        self.validated_header = True

        return
        

    def validate_header_size(self):
        '''Validate that the header size is a multiple of 2880 bytes.'''

        header_size = len(self.header.tostring())
        if 0 != header_size % 2880:
            raise MyException("File header does not conform to 2880 byte standard!")
        
        return


    def validate_primary_header_cards(self):
        '''Validate primary header cards 'SIMPLE, BITPIX, NAXIS, and END' exist and are valid.'''

        if self.header.cards[0].keyword != 'SIMPLE':
            raise MyException("The first keyword in the header must be SIMPLE.")

        required_keys = ['SIMPLE', 'BITPIX', 'NAXIS']
        for key in required_keys:
            if key not in self.header:
                raise MyException(f"Required keyword '{key}' missing in header.")
            
        if self.header.get('BITPIX') not in [8, 16, 32, 64, -32, -64]:
            raise MyException("BITPIX has an invalid value.")
        
        naxis = self.header.get('NAXIS')
        for i in range(1, naxis + 1):
            axis_key = f"NAXIS{i}"
            if axis_key not in self.header:
                raise MyException(f"Missing {axis_key} keyword.")
            if self.header[axis_key] < 0:
                raise MyException(f"{axis_key} must be non-negative.")
            
        with open(self.file_path, 'rb') as f:
            block_size = 2880
            header_bytes = b""
            while True:
                chunk = f.read(block_size)
                if not chunk:
                    break
                header_bytes += chunk
                if b'END' in chunk:
                    break

            if b'END' not in header_bytes:
                raise MyException("FITS header does not contain the END keyword.")
            
        return


    def validate_header_cards(self):
        '''Validate all other header cards to ensure 80 byte length and there are no duplicates.'''

        seen = set()
        for card in self.header.cards:
            key, value, comment = card
            if key in seen and key != "COMMENT" and key != "HISTORY":
                raise MyException(f"Duplicate keyword found: {key}")
            seen.add(key)

        types = []
        for card in self.header.cards:
            key, value, comment = card
            types.append(type(value))

            card_str = card.__str__()
            card_length = len(card_str)

            if card_length != 80:
                raise MyException(f"Card '{key}' is {card_length} characters long: {card_str}")
            
            if value is None or (isinstance(value, str) and not value.strip()):
                self.missing_values.append(key)

        if len(self.missing_values) > 0:
            warnings.warn(f"Values do not exist for {self.missing_values}", stacklevel=2)

        return


    def validate_data(self):
        return
    
file = Mike("C:/Users/starb/Downloads/0136645.fits")
file.validate_primary_header()
file.validate_data()