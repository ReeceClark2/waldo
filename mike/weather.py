import numpy as np
from scipy.stats import linregress
from astropy import units as u
from datetime import datetime
import rcr
import itur

from file_exception import MyException
from file_init import Mike
from val import Val
from sort import Sort
from cal import Cal

import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


class Weather:
    def __init__(self, file):
        '''
        Initialize binary fits cube by reading its header and data from a given file.
        '''

        self.file = file


    def water_vapor_density(self, P, T_K, RH):
        """
        Estimate surface water vapor density (rho_wv) in g/m³.
        Inputs:
            P: Total pressure (not used here, can be removed)
            T_K: Temperature in Kelvin
            RH: Relative humidity in %
        """
        T_C = T_K - 273.15
        e_s = 6.112 * np.exp((17.67 * T_C) / (T_C + 243.5))  # hPa
        e = RH / 100 * e_s  # hPa
        rho_wv = 216.7 * e / T_K  # g/m³

        return rho_wv
    

    def estimate_iwv(self, T_K, RH, H=2000):
        """
        Estimate integrated water vapor (IWV) in mm (equivalent to kg/m²).
        Inputs:
            T_K: Surface temperature in Kelvin
            RH: Relative humidity in %
            H: Water vapor scale height in meters (default 2000 m)
        """
        T_C = T_K - 273.15
        e_s = 6.112 * np.exp((17.67 * T_C) / (T_C + 243.5))  # hPa
        e = RH / 100 * e_s  # hPa
        rho_wv = 216.7 * e / T_K  # g/m³
        iwv = (rho_wv * H) / 1000  # mm or kg/m²

        return iwv
    

    def compute_gaseous_attenuation(self, f, P, T, rh, site_elevation, theta_rad):
        # Compute water vapor density (assumed to be in g/m^3)
        rho = self.water_vapor_density(P, T, rh)

        # Path lengths adjusted for elevation (assumed in km)
        H_o = (8 - site_elevation / 1000)  
        H_w = max((2 - site_elevation / 1000), 0.1) 

        # Slant path factors
        airmass = 1 / np.sin(theta_rad)
        L_o = H_o * airmass  # Path length for oxygen
        L_w = H_w * airmass  # Path length for water vapor

        # Method 2: Summing individual oxygen and water vapor attenuation
        gamma_o = itur.models.itu676.gamma0_exact(f, P, rho, T)
        gamma_w = itur.models.itu676.gammaw_exact(f, P, rho, T)

        A_o = gamma_o.value * L_o  # Oxygen attenuation in dB
        A_w = gamma_w.value * L_w  # Water vapor attenuation in dB
        A = A_o + A_w  # Total attenuation from individual gases in dB

        return A


    def weather_correction(self):
        for ind1, i in enumerate(self.file.data):
            subset_data = self.file.data[ind1]

            P = subset_data['PRESSURE'] * 1.33322 
            T = subset_data['TAMBIENT'] + 273.15
            rh = subset_data['HUMIDITY']   
            theta_rad = np.radians(subset_data['ELEVATIO'])
            
            f = np.linspace(self.file.freqs[ind1][0], self.file.freqs[ind1][1], len(subset_data['DATA'][0]))
            site_elevation = self.file.header['SITEELEV'] / 1000

            for ind2, j in enumerate(subset_data):
                A_gas = self.compute_gaseous_attenuation(f, P[ind2], T[ind2], rh[ind2], site_elevation, theta_rad[ind2])

                # self.compute_rain_attenuation()
                # self.compute_scintillation_attenuation()
                # self.compute_cloud_attenuation()

                transmission = 1 - 10**(-A_gas / 10)

                subset_data['DATA'][ind2] = subset_data['DATA'][ind2] * transmission
        
        return


            
if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    file = Mike("C:/Users/starb/Downloads/0136375.fits")

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.sort()

    w = Weather(file)
    print(file.data[0]['DATA'][0])
    w.weather_correction()
    print(file.data[0]['DATA'][0])