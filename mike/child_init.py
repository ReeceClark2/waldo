import copy
from file_init import Mike

class Sully(Mike):
    def __init__(self, file):
        super().__init__(file.file_path)
        self.data = copy.deepcopy(file.data)
        # Copy other attributes as needed
        self.data_indices = copy.deepcopy(file.data_indices)
        self.freqs = copy.deepcopy(file.freqs)

        #TO be saved
        self.params = []
        self.continuum = copy.deepcopy(file.continuum)
        self.spectrum = copy.deepcopy(file.spectrum)
        # self.parent = file.file_path
