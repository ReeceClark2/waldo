from collections import defaultdict
from scipy.ndimage import label, median_filter
import numpy as np
from file_init import Mike
from file_exception import MyException


class Sort:
    def __init__(self, file):
        '''
        Initialize file and confirm header and data have been validated.
        '''

        self.file = file

        # if self.file.validated_header == False:
        #     raise MyException('FITS header has not been validated!')
        # if self.file.validated_data == False:
        #     raise MyException('FITS data has not been validated!')


    def split_slp_feed(self):
        '''
        Split the data of the provided file by channel and feed.
        '''

        ifnums = np.unique(self.file.data["IFNUM"])
        plnums = np.unique(self.file.data["PLNUM"])
        print(f"IFNUMs: {ifnums}, PLNUMs: {plnums}")

        ind = 0
        data = []
        for i in ifnums:
            for j in plnums:
                subset_mask = (self.file.data["IFNUM"] == i) & (self.file.data["PLNUM"] == j)
                subset_data = self.file.data[subset_mask]
                data.append(subset_data)
                ind += 1
        self.file.data = data

        return 


    def divide_section(self, section):
        '''
        Locate indices within some channel of the data.

        Params:
        section: astropy FITS table

        Returns:
        array: 1st index of data, 1st index of post calibration
        '''

        data_start_ind = None
        post_cal_start_ind = None

        cal_started = False
        pre_cal_complete = False

        for ind, point in enumerate(section):
            if point["CALSTATE"] == 1:
                cal_started = True

            if cal_started and point["CALSTATE"] == 0 and point["SWPVALID"] == 1 and not pre_cal_complete:
                data_start_ind = ind
                pre_cal_complete = True

            if pre_cal_complete and point["SWPVALID"] == 0 and section[ind-1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind
            else:
                post_cal_start_ind = None

            if pre_cal_complete and point["SWPVALID"] == 0 and point["CALSTATE"] == 1:
                break


        if data_start_ind is None:
            for ind, point in enumerate(section):
                if point["SWPVALID"] == 1 and data_start_ind is None:
                    data_start_ind = ind
                    continue  

                if data_start_ind is not None and point["SWPVALID"] == 0:
                    post_cal_start_ind = ind
                    break

        return np.array([data_start_ind, post_cal_start_ind])


    def sort_data(self):
        '''
        Iterate through each channel of data and append to data_indices.
        '''

        for i in self.file.data:
            result = self.divide_section(i)
            self.file.data_indices.append(result)

        return


    def section_debug(self):
        '''
        Iterate through data and check the pre calibration, data, and post calibration of each.
        '''

        for ind, i in enumerate(self.file.data):
            current_indices = self.file.data_indices[ind]
            data_start_ind = current_indices[0]
            post_cal_start_ind = current_indices[1]

            pre_cal = i[:data_start_ind]
            data = i[data_start_ind:post_cal_start_ind]
            post_cal = i[post_cal_start_ind:]

            print(len(pre_cal),len(data),len(post_cal))


if __name__ == "__main__":
    file = Mike("C:/Users/starb/Downloads/0136870.fits")

    np.set_printoptions(threshold=100000)

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.section_debug()
    
    print(file.data[0]['CALSTATE'])
    print(file.data[0]['SWPVALID'])
