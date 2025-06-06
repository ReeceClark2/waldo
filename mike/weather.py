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
        
        self.site_elevation = self.file.header['SITEELEV'] / 1000
        
        self.oxygen_data = np.loadtxt('tables/oxygen_coefficients.txt')
        self.oxygen_frequencies = self.oxygen_data[:, 0]
        self.oxygen_coefficients = self.oxygen_data[:, 1:5]

        self.water_vapor_data = np.loadtxt('tables/water_vapor_coefficients.txt')
        self.water_vapor_frequencies = self.water_vapor_data[:, 0]
        self.water_vapor_coefficients = self.water_vapor_data[:, 1:5]

        pass


    def debye_spec(self, f, p, T):
        """
        Compute the imaginary part of the oxygen Debye spectrum N''_D(f).

        Parameters:
            f : frequency in GHz (scalar or array)
            p : total pressure in hPa
            T : temperature in Kelvin

        Returns:
            N_imag : imaginary part of the refractivity (dimensionless)
        """

        theta = 300 / T
        d = 5.6 + 3.75 * theta

        first_term = 6.14e-5 / (d * (1 + (f / d) ** 2))
        second_term = (1.4e-12 * p * theta ** 1.5) / (1 + 1.9e-5 * f ** 1.5)

        N_imag = f * p * theta ** 2 * (first_term + second_term)

        return N_imag


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


    def get_coefficients(self, f, medium):
        if medium == 'O':
            aO = np.interp(f, self.oxygen_frequencies, self.oxygen_coefficients[:, 0])
            bO = np.interp(f, self.oxygen_frequencies, self.oxygen_coefficients[:, 1])
            cO = np.interp(f, self.oxygen_frequencies, self.oxygen_coefficients[:, 2])
            dO = np.interp(f, self.oxygen_frequencies, self.oxygen_coefficients[:, 3])
            return aO, bO, cO, dO
        elif medium == 'V':
            aV = np.interp(f, self.water_vapor_frequencies, self.water_vapor_coefficients[:, 0])
            bV = np.interp(f, self.water_vapor_frequencies, self.water_vapor_coefficients[:, 1])
            cV = np.interp(f, self.water_vapor_frequencies, self.water_vapor_coefficients[:, 2])
            dV = np.interp(f, self.water_vapor_frequencies, self.water_vapor_coefficients[:, 3])
            return aV, bV, cV, dV
        else:
            return


    def oxygen_attenuation(self, f, p, T, rho, theta):
        gamma = 0.1820 * self.debye_spec(f, p, T)

        aO, bO, cO, dO = self.get_coefficients(f, 'O')

        hO = aO + bO * T + cO * p + dO * rho

        AO = (gamma * hO / np.sin(theta)) * ((8 - self.site_elevation) / 8)

        return AO
    

    def water_vapor_attenuation(self, f, p, T, rho, theta, rh):
        iwv = self.estimate_iwv(T, rh)

        aV, bV, cV, dV = self.get_coefficients(f, 'V')

        kV = aV + bV * rho + cV * T + dV * p

        AW = (kV * iwv / np.sin(theta)) * ((2 - self.site_elevation) / 2)

        return AW


    def total_attenuation(self):
        # Frequency range (GHz)
        f = np.linspace(1, 12, 100) * u.GHz

        # Environmental parameters
        p = self.file.data[0][0]['PRESSURE'] * 1.33322 * u.hPa  # Convert mmHg to hPa
        T = 280 * u.K
        rh = 100  # Relative humidity %
        theta_rad = np.radians(40)
        site_elevation = 23 * u.m

        # Compute water vapor density
        rho = self.water_vapor_density(p.value, T.value, rh) * u.g / u.m**3

        # Path lengths adjusted for elevation
        H_o = (8 * u.km - site_elevation).to(u.km)
        H_w = max((2 * u.km - site_elevation).to(u.km).value, 0.1) * u.km  # prevent negative or zero

        # Slant path factors
        airmass = 1 / np.sin(theta_rad)
        L_o = H_o * airmass
        L_w = H_w * airmass

        # Specific attenuation (dB/km)
        gamma_o = itur.models.itu676.gamma0_exact(f, p, rho, T)
        gamma_w = itur.models.itu676.gammaw_exact(f, p, rho, T)

        # Total attenuation (dB)
        A_o = gamma_o * L_o
        A_w = gamma_w * L_w
        A_total = (A_o + A_w).to(u.dB).value

        # Convert to transmission %
        transmission = 10 ** (-A_total / 10) * 100

        # Plot
        plt.plot(f.value, transmission)
        plt.xlabel("Frequency (GHz)")
        plt.ylabel("Transmission (%)")
        plt.title("Transmission vs Frequency")
        plt.grid(True)
        plt.savefig('thing')

            


    def total_attenuation_in_house(self):
        for subset_data in self.file.data:
            for point in subset_data:
                intensities = point['DATA']

                f = 8
                p = self.file.data[0][0]['PRESSURE'] * 1.33322 # Convert to hPa
                T = self.file.data[0][0]['TAMBIENT'] + 273.15 # C
                rh = self.file.data[0][0]['HUMIDITY']

                rho = self.water_vapor_density(p, T, rh)
                theta = np.radians(self.file.data[0][0]['ELEVATIO'])

                AO = self.oxygen_attenuation(f, p, T, rho, theta)
                AW = self.water_vapor_attenuation(f, p, T, rho, theta, rh)
                print("AO:", AO, "  AW:", AW)

                A = AO + AW
                
                trans = 10 ** (-A / 10)

                print(trans)


    def total_attenuation_plotting(self):
        # Normalize over the full range of θ (deg), T (K), RH (%)
        norm = Normalize(vmin=30 + 266 + 0, vmax=90 + 306 + 95)
        cmap = cm.viridis

        angles = np.linspace(30, 90, 2)      # elevation in degrees
        temps = np.linspace(266, 306, 2)     # temperature in K
        rhs = np.linspace(0, 95, 2)          # relative humidity in %

        for theta_deg in angles:
            theta_rad = np.radians(theta_deg)
            for T in temps:
                for rh in rhs:
                    f = np.linspace(1, 8, 100)
                    p = self.file.data[0][0]['PRESSURE'] * 1.33322

                    rho = self.water_vapor_density(p, T, rh)
                    AO = self.oxygen_attenuation(f, p, T, rho, theta_rad)
                    AW = self.water_vapor_attenuation(f, p, T, rho, theta_rad, rh)

                    A = AO + AW
                    # Compute specific attenuations
                    gamma_o = itur.models.itu676.gamma0_exact(f, p, rho, T) * (8 / np.sin(theta_rad)) # Oxygen
                    gamma_w = itur.models.itu676.gammaw_exact(f, p, rho, T) * (2 / np.sin(theta_rad))  # Water vapor

                    # Total specific attenuation
                    gamma_total = (gamma_o + gamma_w).value
                    
                    trans = (10 ** (-A / 10))
                    trans_itur = (10 ** (-gamma_total / 10))

                    # Color is determined by combined T + RH + θ
                    scalar = T + rh + theta_deg
                    color_value = norm(scalar)
                    plt.plot(f, np.array(trans) * 100, color='blue',
                            label=f'In-House')
                    plt.plot(f, np.array(trans_itur) * 100, color='red',
                            label=f'ITUR')

            plt.xlabel('Frequency (GHz)')
            plt.ylabel('Attenuation')
            plt.title(f'Attenuation ITUR vs In-House Models')
            # plt.legend(fontsize='small')
            plt.savefig(f'thing.png', dpi=300)
            plt.close()


if __name__ == "__main__":
    """Test function to implement continuum data handling."""

    file = Mike("C:/Users/starb/Downloads/0136375.fits")

    v = Val(file)
    # v.validate_primary_header()
    # v.validate_data()

    s = Sort(file)
    s.split_slp_feed()
    s.sort_data()

    w = Weather(file)
    w.total_attenuation()