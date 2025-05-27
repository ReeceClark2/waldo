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

class Val:
    def __init__(self, file):
        '''Initialize binary fits cube by reading its header and data from a given file.'''

        self.file = file

        pass


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

        self.file.validated_header = True

        return
        

    def validate_header_size(self):
        '''Validate that the header size is a multiple of 2880 bytes.'''

        header_size = len(self.file.header.tostring())
        if 0 != header_size % 2880:
            raise MyException("File header does not conform to 2880 byte standard!")
        
        return


    def validate_primary_header_cards(self):
        '''Validate primary header cards 'SIMPLE, BITPIX, NAXIS, and END' exist and are valid.'''

        if self.file.header.cards[0].keyword != 'SIMPLE':
            raise MyException("The first keyword in the header must be SIMPLE.")

        required_keys = ['SIMPLE', 'BITPIX', 'NAXIS']
        for key in required_keys:
            if key not in self.file.header:
                raise MyException(f"Required keyword '{key}' missing in header.")
            
        if self.file.header.get('BITPIX') not in [8, 16, 32, 64, -32, -64]:
            raise MyException("BITPIX has an invalid value.")
        
        naxis = self.file.header.get('NAXIS')
        for i in range(1, naxis + 1):
            axis_key = f"NAXIS{i}"
            if axis_key not in self.file.header:
                raise MyException(f"Missing {axis_key} keyword.")
            if self.file.header[axis_key] < 0:
                raise MyException(f"{axis_key} must be non-negative.")
            
        with open(self.file.file_path, 'rb') as f:
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
        for card in self.file.header.cards:
            key, value, comment = card
            if key in seen and key != "COMMENT" and key != "HISTORY":
                raise MyException(f"Duplicate keyword found: {key}")
            seen.add(key)

        types = []
        for card in self.file.header.cards:
            key, value, comment = card
            types.append(type(value))

            card_str = card.__str__()
            card_length = len(card_str)

            if card_length != 80:
                raise MyException(f"Card '{key}' is {card_length} characters long: {card_str}")
            
            if value is None or (isinstance(value, str) and not value.strip()):
                self.file.missing_values.append(key)

        if len(self.file.missing_values) > 0:
            warnings.warn(f"Values do not exist for {self.file.missing_values}", stacklevel=2)

        return


    def validate_data(self):
        try:
            with fits.open(self.file.file_path) as hdul:
                self.dataH = hdul[1].header
                self.data = Table(hdul[1].data)
            #print(repr(self.dataH))

        except Exception as e:
            raise MyException(f"Error reading data from FITS file: {e}")


        
        # Check if the data has bad values (NaN or negative)
        Data = np.ma.array(self.data["DATA"])
        # Check each array in Data for NaN or 0s and create a masked array
        mask = np.isnan(Data) | (Data < 0)
        if np.any(mask):
            warnings.warn("🚫 Data contains negative or NaN values.", stacklevel=2)
        self.data["DATA"] = np.ma.masked_array(Data, mask=mask)

        # Check if the data types are correct and converts the time columns to be datetime or floats
        self.validate_types()

        self.file.validated_data = True

        return


    def validate_types(self):
        '''Checks that all values in a given column are the correct type.'''

        #check everry column
        for column in self.data.colnames:
            # check the types
            column_data = self.data[column]
                #print(f"🚫 Column '{column}' contains values that do not match its data type. Expected data type: {column_dtype}.")
                #fix the type mismatch
            
            self.match_types(column)
            self.convert_to_datetime(column)
            #else:
                #print(f"✅ Column '{column}' data types are valid.")
            #convert the time columns to datetime
            
            #check if the numeric columns have weird numbers
            self.check_numbers(column)


    def match_types(self, column):
        '''Make sure all the column types match the data types of the values.'''

        column_data = self.data[column]
        #find how many unique types there are in the column
        unique_types = set(type(value) for value in self.data[column])
        
        #if there are more than one type of value in a column, do smth about it
        if len(unique_types) > 1:
            #try to convert the values to the most common type
            common_type = max(unique_types, key=lambda t: sum(isinstance(value, t) for value in self.data[column]))
            for i, value in enumerate(self.data[column]):
                if not isinstance(value, common_type):
                    try:
                        self.data[column][i] = common_type(value)
                    except Exception:
                        # Mask the messed up values if conversion fails
                        if not hasattr(self.data[column], 'mask'):
                            self.data[column] = np.ma.array(self.data[column])
                        self.data[column].mask[i] = True
            print(f"🚫 Column '{column}' contains mixed data types: {unique_types}.")
        else:
            # otherwise convert the column to the common type
            common_type = unique_types.pop()
            try:

                if common_type is str:
                    # automatic conversion to str doesn't work some of the time so we need to try UTF-8 first
                    if column_data.dtype.char == 'S':  # Check if it's a byte string
                        self.data[column] = np.char.decode(self.data[column], encoding='utf-8', errors='replace')
                    else:
                        self.data[column] = np.array(self.data[column], dtype=str)

                # make sure the arrays in the column are all floats
                elif common_type is np.ndarray:
                    if not all(isinstance(float(item), float) for subarray in column_data for item in subarray):
                        #if not try to make them floats, otherwise raise an error
                        for subarray in column_data:
                            for idx, item in enumerate(subarray):
                                try:
                                    float(item)
                                except Exception:
                                    raise TypeError(f"🚫 Column '{column}' contains arrays with non-float elements: {item}")
                    #else:
                        #print(f"✅ All elements in arrays of column '{column}' are floats.")

                # otherwise just convert the column to the common type
                else:
                    self.data[column] = self.data[column].astype(common_type)

                #print(f"🔄 Column '{column}' data type changed to {common_type}.")

            except Exception as e:
                #if that all fails, raise an error
                MyException(f"⚠️ Failed to change column '{column}' data type to {common_type}: {e}")

        return
    

    def convert_to_datetime(self, column):
        '''Convert the dates and times to datetime objects from the date time library'''
        
        #check if the column name contains any of the keywords
        if any(keyword in column.upper() for keyword in ["DATE", "TIME", "DURATION", "EXPOSURE", "MJD", "UTC", "UTSECS"]) or column.upper() == "LST":
                try:
                    # Attempt to convert to datetime object
                    self.data[column] = np.vectorize(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S.%f"))(self.data[column])
                    self.data[column] = self.data[column].astype(type(self.data[column][0]))
                except (ValueError, TypeError):
                    try:
                        # If conversion fails, convert to seconds (float)
                      self.data[column] = Column(self.data[column].astype(float), dtype='f8')

                    except ValueError:
                        try:
                            # If conversion fails, make it a string
                            self.data[column] = Column(self.data[column].astype(str), dtype='U')
                        except Exception:
                            raise MyException(f"⚠️ Failed to convert value '{self.data[column]}' in column '{column}' to datetime or seconds.")
                

        return

    def check_numbers(self, column):
        '''Check if certain columns have negative values and remove them'''

        if self.data[column].dtype.kind in {'U', 'S', 'O'}:
            # Vectorized replacement of 'nan' or 'NaN' (case-insensitive) with np.nan in string columns
            def replace_nan(val):
                if isinstance(val, str) and val.strip().lower() == 'nan':
                    return np.nan
                return val
            self.data[column] = np.vectorize(replace_nan)(self.data[column])

        if column.upper() in ["DURATION", "EXPOSURE", "TSYS", "TCAL", "LST", "ELEVATION", "TAMBIENT", "PRESSURE", "HUMIDITY", "RESTFREQ", "FREQRES", "TRGTLONG", "MJD", "UTSECS" ]:
            if np.any(self.data[column] < 0):
                num_negatives = np.sum(self.data[column] < 0)
                print(f"Found {num_negatives} negative values in column '{column}'. out of {len(self.data[column])} total values.")
                self.data = self.data[self.data[column] >= 0]
        
        return
    

if __name__ == "__main__":
    '''Test function to implement validation.'''
    start = time.time()
    # #multiple file testing
    folder_path = "TrackingLowRes"  # Replace with your folder path
    average_time = 0
    for fileN in os.listdir(folder_path):
        start = time.time()
        if fileN.lower().endswith(".fits"):
            file_path = os.path.join(folder_path, fileN)
            try:
                file = Mike(file_path)
                v = Val(file)
                v.validate_primary_header()
                v.validate_data()
                print(f"{fileN}: Headers ☺ - {file.validated_header}")
                print(f"{fileN}: Data ☺ - {file.validated_data}")
            except Exception as e:
                print(f"{fileN}: Validation failed - {e}")            
        end = time.time()
        print(f"Validation took {end - start:.2f} seconds.")
        average_time += (end - start) / len(os.listdir(folder_path))
    print(f"Average validation time: {average_time:.2f} seconds.")
    #single file testing
    fileN = "TrackingGl.fits" 
    file = Mike(fileN)
    #print(file.data.header)
    v = Val(file)
    v.validate_primary_header()
    v.validate_data()

    print(str(file.validated_header) +" " + fileN + " Headers ☺")
    print(str(file.validated_data) + " " + fileN + " Data ☺")
    end = time.time()
    print(f"Validation took {end - start:.2f} seconds.")
