from astropy.io import fits
from file_exception import MyException
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Mike:
    def __init__(self, file_path):
        '''
        Initialize binary fits cube by reading its header and data from a given file.
        
        File fields:
        missing values: empty primary header values
        file_path: path to file in directory

        validated_header: check whether header has been validated
        validated_data: check whether data has been validated

        data_indices: list of arrays of shape (1,2) (track)
            1st value is the first data point in the corresponding indexed data array
            2nd value is the first post calibration point in the corresponding indexed data array
        data_indices: list of arrays of shape (1,2,3,4) (on/off)
            1st value is the first on data point in the corresponding indexed data array
            2nd value is the first transitioning point in the corresponding indexed data array
            3rd value is the first off data point in the corresponding indexed data array
            4th value is the first post calibration point in the corresponding indexed data array
        gain_start: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        gain_end: list of arrays of shape (1,2)
            1st value is the gain delta for the given index
            2nd value is the time for the gain delta
        continuum:
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

        self.labels = []
        self.data_indices = []
        self.gain_start = []
        self.gain_end = []

        self.continuum = []
        self.spectrum = [[]]


if __name__ == "__main__":
    '''
    Test function to test uploaded file.
    '''

    file = Mike("C:/Users/starb/Downloads/0136870.fits")

    np.set_printoptions(threshold=100000)

    print(repr(file.header))
    print(repr(file.data['OBSID']))

    # file = Mike("C:/Users/starb/Downloads/0136869.fits")

    # np.set_printoptions(threshold=100000)

    # print(repr(file.header))

    fig, ax1 = plt.subplots()

    # Plot ELEVATIO on the primary y-axis
    ax1.plot(file.data["CALSTATE"], color='tab:blue', label='CALSTATE')
    ax1.set_ylabel("CALSTATE", color='tab:blue')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.set_xlabel("Index")
    ax1.set_ylim(-1.1, 1.1)  # Set SWPVALID y-axis from 0 to 1

    # Create a second y-axis sharing the same x-axis
    ax2 = ax1.twinx()
    ax2.plot(file.data["SWPVALID"], color='tab:red', label='SWPVALID', alpha=0.5)
    ax2.set_ylabel("SWPVALID", color='tab:red')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    ax2.set_ylim(-0.1, 1.1)  # Set SWPVALID y-axis from 0 to 1

    # Optional: Add a legend
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)

    # Save the figure
    plt.title("CALSTATE and SWPVALID")
    plt.savefig("CALSTATE.png", dpi=300, bbox_inches='tight')
    plt.close()