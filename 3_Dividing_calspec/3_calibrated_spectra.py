
from astropy.table import Table
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import os


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def tanspec_data(q=0):
    file_path = '/home/varghese/Desktop/DP2_TANSPEC/\
Data/Output_spectra/hip41815/s0.5/tanspec_slopes'
    flux = fits.getdata(os.path.join(file_path,
                                     'op_hip41815_s0.5_grating1_A.ms.wlc.fits')
                        )
    order = flux[q]
    wavelength = np.load(
        os.path.join(file_path,
                     'op_hip41815_s0.5_grating1_A_combarc.OutputWavlFile{}.npy'
                     .format(q))
    )
    return order, wavelength


calspec_data = Table.read(
    '/home/varghese/Desktop/DP2_TANSPEC/Calspec_data/hd205905_mod_004.fits'
)

calspec_wl = np.array(calspec_data['WAVELENGTH'])
calspec_flux = np.array(calspec_data['FLUX'])
smoothed_flux = smooth(calspec_flux, 19)
interpolation = interpolate.interp1d(calspec_wl, smoothed_flux)

fig, axs = plt.subplots(5, 2, figsize=(16, 9))
for order in range(0, 10, 1):
    tanspec_flux, wavelength = tanspec_data(order)
    flux_ratio = np.load('smoothed_ratio_{}.npy'.format(order))
    corrected_flux = flux_ratio * tanspec_flux
    cal_orders = interpolation(wavelength)
    axs[order//2][order % 2].plot(wavelength, corrected_flux)
    # axs[order//2][order % 2].plot(wavelength, cal_orders)
    # axs[order//2][order % 2].plot(wavelength, flux_ratio)
    axs[order//2][order % 2].set_title('order {}'.format(order))
plt.savefig('corrected_spectra.png')

# Code ends here
