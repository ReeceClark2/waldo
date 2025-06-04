from collections import defaultdict
from tabnanny import check
from scipy.ndimage import label, median_filter
import numpy as np
from file_init import Mike
from file_exception import MyException
from astropy.table import Table


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
    
    def get_startend_freqs(self):
        '''
        Get the start and stop frequencies for each channel in the data.
        param file: Mike class file
        returns: populates the file's freqs field with start and stop frequencies
        '''
        freqs = []
        for i in range(len(self.file.data)):
            #scour the header for the bandwidth and center frequencies
            for key, value in self.file.header.items():
                #bandwidth
                if key == ("OBSBW"):
                    band = value
                #lowres center frequency
                elif key == ("OBSFREQ"):
                    center = [value]
                #HIRES bands center frequencies
                elif key == ("HISTORY"):
                    #if HIRES BANDS exist replace center with the HIRES center frequencies
                    if value.startswith("HIRES bands"):
                        # Extract all integers from the string
                        value = value.replace(",", " ").strip()
                        value = value.split(" ")
                        center = []
                        #split the value into individual words and numbers
                        for k in value:
                            k = str(k).strip()
                            try:
                                #if it's a float add it to the center list
                                k = float(k)
                                center.append(k)
                            except ValueError:
                                #otherwise skip it
                                continue
                    # start/stop channels
            for c in center:
                start_f = c - (band / 2)
                stop_f = c + (band / 2)
            freqs.append(np.array([start_f, stop_f]))
        self.file.freqs = freqs
    
    def get_startstop_channels(self):
        '''
        Get the start and stop channels for each channel in the data.
        param file: Mike class file
        returns: cuts out the channels in the data that are not in the start and stop channels
        '''

        for i in range(len(self.file.data)):
            #search through the header 
            for key, value in self.file.header.items():
                if key == ("HISTORY"):
                    if value.startswith("START,STOP"):
                        # Extract all integers from the string
                        value = value.replace(",", " ").strip()
                        value = value.split(" ")
                        #split the value into individual words and numbers
                        channels = []
                        for k in value:
                            k = str(k).strip()
                            try:
                                #if it's an integer add it to the channels list
                                k = int(k)
                                channels.append(k)
                            except ValueError:
                                #if it can't become a integer, skip it
                                continue
                        # Remove all string type characters from the channels list
                        channel1 = int(channels[0])
                        channel2 = int(channels[1])
            # Cut the DATA column to only include the channels in the start and stop channels
            t = Table(self.file.data[i])
            t.replace_column('DATA', np.array([row[channel1:channel2] for row in t['DATA']]))
            self.file.data[i] = t

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
    file = Mike("TrackingHighRes/0136376.fits")

    np.set_printoptions(threshold=100000)

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()
    s.section_debug()
    s.get_startend_freqs()
    s.get_startstop_channels()
    
    print(file.data[0]['CALSTATE'])
    print(file.data[0]['SWPVALID'])
