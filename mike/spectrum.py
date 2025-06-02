import numpy as np
from scipy.stats import linregress
from datetime import datetime
import rcr

from file_exception import MyException
from file_init import Mike
from val import Val
from sort import Sort
from cal import Cal

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Spec:
    def __init__(self, file):
        '''
        Initialize file.
        '''

        self.file = file


    def average(self, data, axis):
        '''
        Integrate across time (axis = 0) or integrate across frequency (axis = 1).
        '''

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means

    

    def make_spec(self):
        result = self.average(self.file.data[0], 0)
        plt.plot(result)
        plt.savefig("Spectrum")



if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''

    file = Mike("C:/Users/starb/Downloads/0136981.fits")
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()

    c = Cal(file)
    c.gain_calibration()

    spec = Spec(file)
    spec.make_spec()

