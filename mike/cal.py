from astropy.io import fits
from scipy.special import erfc
from scipy.ndimage import label
from scipy.stats import linregress
import numpy as np
from file_exception import MyException
from datetime import datetime
import rcr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from file_init import Mike
from val import Val
from sort import Sort


class Cal:
    def __init__(self, file):
        '''
        Initialize file and confirm header and data have been validated.
        '''

        self.file = file

        # if self.file.validated_header == False:
        #     raise MyException('FITS header has not been validated!')
        # if self.file.validated_data == False:
        #     raise MyException('FITS data has not been validated!')
    

    def average(self, data, axis):
        '''
        Average across time (axis = 0) or integrate across frequency (axis = 1).
        '''

        intensities = np.array([row[6] for row in data])
        count = intensities.shape[axis]
        channel_means = np.sum(intensities, axis=axis) / count

        return channel_means
    

    def linear(self, x, params): # model function
        return params[0] + x * params[1]


    def d_linear_1(self, x, params): # first model parameter derivative
        return 1


    def d_linear_2(self, x, params): # second model parameter derivative
        return x


    def rcr(self, array):
        '''
        Perform functional Robust Chauvenet Rejection on a given dataset.
        '''

        x = array[0]
        y = array[1]

        result = linregress(x, y)
        m = result.slope
        b = result.intercept

        guess = [m, b]
        model = rcr.FunctionalForm(self.linear,
            x,
            y,
            [self.d_linear_1, self.d_linear_2],
            guess
        )

        r = rcr.RCR(rcr.LS_MODE_68) 
        r.setParametricModel(model)
        r.performBulkRejection(y)

        best_fit_parameters = model.result.parameters

        return best_fit_parameters


    def sdfits_to_array(self, file, data):
        '''
        Convert sdfits data to more accessible, lighter arrays.
        '''

        freq = self.average(data, axis=1)

        times = [datetime.fromisoformat(t) for t in data["DATE-OBS"]]
        t0 = datetime.fromisoformat(file.header["DATE"])
        time_rel = [(t - t0).total_seconds() for t in times]

        result = ([np.array(time_rel), np.array(freq)])

        return result
    

    def compute_gain_deltas(self, file, ind):
        subset_data = file.data[ind]
        subset_indices = file.data_indices[ind]

        try:
            pre_cal = subset_data[
                (np.arange(len(subset_data)) < subset_indices[0]) &
                (subset_data["SWPVALID"] == 0)
            ]
        except:
            pre_cal = None
            
        try:
            post_cal = subset_data[
                (np.arange(len(subset_data)) >= subset_indices[1]) &
                (subset_data["SWPVALID"] == 0)
            ]
        except:
            post_cal = None


        def get_delta(cal):
            cal_on_array = self.sdfits_to_array(file, cal[cal["CALSTATE"] == 1])
            cal_on_params = self.rcr(cal_on_array)
            cal_off_array = self.sdfits_to_array(file, cal[cal["CALSTATE"] == 0])
            cal_off_params = self.rcr(cal_off_array)

            time = (np.mean(cal_on_array[0]) + np.mean(cal_off_array[0])) / 2
            if time < (cal_on_array[0][0] + cal_off_array[0][-1]) / 2:
                time = (cal_on_array[0][0] + cal_off_array[0][-1]) / 2
            elif time > (cal_on_array[0][-1] + cal_off_array[0][0]) / 2:
                time = (cal_on_array[0][-1] + cal_off_array[0][0]) / 2

            delta = (cal_on_params[1] * time + cal_on_params[0]) - (cal_off_params[1] * time + cal_off_params[0])
            
            return delta, time
        

        if pre_cal is not None:
            delta1, t1 = get_delta(pre_cal)
            file.gain_start.append([delta1, t1])
        else:
            file.gain_start.append(None)

        if post_cal is not None:
            delta2, t2 = get_delta(post_cal)
            file.gain_end.append([delta2, t2])
        else:
            file.gain_end.append(None)

        return


    def gain_calibration(self, file):
        '''
        Carry out gain calibration for a given file.
        '''

        for ind, i in enumerate(file.data):
            self.compute_gain_deltas(file, ind)

        return
        
                
    def source_calibration(self, cal_file):
        '''
        Carry out source calibration with currently selected file and uploaded calibration file.
        '''



        return



if __name__ == "__main__":
    '''
    Test function to implement calibration.
    '''

    file = Mike("C:/Users/starb/Downloads/0115701.fits")
    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()
    
    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.clean_sections()

    cal_file = Mike("C:/Users/starb/Downloads/0115701.fits")
    
    s = Sort(cal_file)
    s.split_slp_feed()
    s.sort_data()
    s.clean_sections()

    c = Cal(file)
    # c.test_plot(file.data)
    c.gain_calibration(file)
    
    print(file.gain_start)
    print(file.gain_end)
