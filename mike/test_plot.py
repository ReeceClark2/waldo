from astropy.io import fits
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

class CalibrationPlotter:
    def __init__(self):
        self.fig, self.axes = plt.subplots(2, 2, figsize=(12, 8))
        self.twin_axes = [ax.twinx() for ax in self.axes.flat]
        self.fig.subplots_adjust(hspace=0.4, wspace=0.3)
    
    def easy_graph(self, uncalibrated_data, calibrated_data, ind):
        """
        Plot calibrated vs uncalibrated data on subplot indexed 1 to 4.
        Each subplot uses twin y-axes.
        """
        if not (1 <= ind <= 4):
            raise ValueError("Index must be between 1 and 4")
        
        # Convert 1-based index to 0-based
        i = ind - 1
        ax1 = self.axes.flat[i]
        ax2 = self.twin_axes[i]

        # Clear in case reused
        ax1.cla()
        ax2.cla()

        # Plot calibrated
        ax1.plot(calibrated_data[0], calibrated_data[1], color='tab:blue', label='Calibrated Intensity')
        ax1.set_ylabel("Calibrated", color='tab:blue')
        ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax1.set_xlabel("Index")
        # ax1.set_ylim(59, 63)

        # Plot uncalibrated
        ax2.plot(uncalibrated_data[0], uncalibrated_data[1], color='tab:red', label='Uncalibrated Intensity', alpha=0.5)
        ax2.set_ylabel("Uncalibrated", color='tab:red')
        ax2.tick_params(axis='y', labelcolor='tab:red')
        # ax2.set_ylim(295, 320)

        # Title per subplot
        ax1.set_title(f"Plot {ind}: Calibrated vs Uncalibrated")

    def save(self, filename="DuelCalibration.png"):
        """Call this after all plots have been filled to save the final figure."""
        self.fig.suptitle("Calibration Comparisons", fontsize=16)
        self.fig.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close(self.fig)


def split_polarity(data):
    """Split data into X and Y polarization channels based on PLNUM (0 = X, 1 = Y)."""
    x_pol = data[data["PLNUM"] == 0]
    y_pol = data[data["PLNUM"] == 1]
    return x_pol, y_pol


def split_calibration(data):
    """Split data into calibration and non-calibration sets based on CAL (-1, 1 = cal, 0 = no cal)."""
    non_cal = data[data["CALSTATE"] == 0]
    cal = data[(data["CALSTATE"] == 1) | (data["CALSTATE"] == -1)]
    return non_cal, cal


def split_slp(data):
    """Dynamically split data into channels based on IFNUM (e.g., pol0, pol1, ...)."""
    split = {}
    ifnums = np.unique(data["IFNUM"])
    print(ifnums)
    for ifnum in ifnums:
        split[f"pol{ifnum}"] = data[data["IFNUM"] == ifnum]
    return split


def integrate(data, axis):
    """Integrate across time (axis = 0) or integrate across frequency (axis = 1)."""
    intensities = np.array([row[6] for row in data])
    count = intensities.shape[axis]
    channel_means = np.sum(intensities, axis=axis) / count
    return channel_means


def graph(xxs: dict, yys: dict):
    """Plot continuum and spectrum data for each polarization IFNUM channel (e.g., pol0, pol1, ...)."""
    n = len(xxs)
    fig, axs = plt.subplots(n, 2, figsize=(10, 4 * n))

    if n == 1:
        axs = np.expand_dims(axs, axis=0) 

    for i, key in enumerate(xxs):
        xx = xxs[key]
        yy = yys.get(key, None) 

        # Integrate across frequency (to get time series)
        xx_freq = integrate(xx, axis=1)
        yy_freq = integrate(yy, axis=1) if yy is not None else None

        # Integrate across time (to get spectrum)
        xx_time = integrate(xx, axis=0)
        yy_time = integrate(yy, axis=0) if yy is not None else None

        # Time axis in seconds
        times_xx = [datetime.fromisoformat(t) for t in xx["DATE-OBS"]]
        t0 = times_xx[0]
        time_xx = [(t - t0).total_seconds() for t in times_xx]

        time_yy = None
        if yy is not None:
            times_yy = [datetime.fromisoformat(t) for t in yy["DATE-OBS"]]
            time_yy = [(t - t0).total_seconds() for t in times_yy]

        # Left: Continuum (intensity over time)
        axs[i][0].plot(time_xx, xx_freq, label=f"{key.upper()} (X)", marker='o', markersize=3)
        if yy_freq is not None:
            axs[i][0].plot(time_yy, yy_freq, label=f"{key.upper()} (Y)", marker='s', markersize=3)
        axs[i][0].set_title(f"Continuum ({key.upper()})")
        axs[i][0].set_xlabel("Time (s)")
        axs[i][0].set_ylabel("Total Intensity")
        axs[i][0].legend()
        axs[i][0].grid(True)

        # Right: Spectrum (intensity over channel)
        axs[i][1].plot(xx_time, label=f"{key.upper()} (X)", marker='o', markersize=3)
        if yy_time is not None:
            axs[i][1].plot(yy_time, label=f"{key.upper()} (Y)", marker='s', markersize=3)
        axs[i][1].set_title(f"Spectrum ({key.upper()})")
        axs[i][1].set_xlabel("Channel")
        axs[i][1].set_ylabel("Integrated Intensity")
        axs[i][1].legend()
        axs[i][1].grid(True)

    plt.tight_layout()
    plt.savefig("3C286_polarization_all_channels.png")

# with fits.open("C:/Users/starb/Downloads/0136645.fits") as hdul:
with fits.open("C:/Users/starb/Downloads/0136376.fits") as hdul:
    header = hdul[0].header
    data = hdul[1].data



non_cal, cal = split_calibration(data[1:])
x_pol, y_pol = split_polarity(data[5:])

xxs = split_slp(x_pol)
yys = split_slp(y_pol)


np.set_printoptions(threshold=np.inf)

graph(xxs, yys)
