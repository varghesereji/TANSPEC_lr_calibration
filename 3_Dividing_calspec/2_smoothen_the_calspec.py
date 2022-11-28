
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
Data/Output_spectra/hd37725/s0.5/tanspec_slopes'
    flux = fits.getdata(os.path.join(file_path,
                                     'op_hd37725_grating1_A.ms.wlc.fits')
                        )
    order = flux[q]
    wavelength = np.load(
        os.path.join(file_path,
                     'op_hd37725_grating1_A_combarc.OutputWavlFile{}.npy'
                     .format(q))
    )
    return order, wavelength


calspec_data = Table.read(
    '/home/varghese/Desktop/DP2_TANSPEC/Calspec_data/hd37725_mod_004.fits'
)

calspec_wl = np.array(calspec_data['WAVELENGTH'])
calspec_flux = np.array(calspec_data['FLUX'])
smoothed_flux = smooth(calspec_flux, 19)
interpolation = interpolate.interp1d(calspec_wl, smoothed_flux)

fig, axs = plt.subplots(5, 2, figsize=(16, 9))
for order in range(0, 10, 1):
    tanspec_flux, wavelength = tanspec_data(order)
    cal_orders = interpolation(wavelength)
    flux_ratio = cal_orders / tanspec_flux
    # axs[order//2][order % 2].plot(wavelength, tanspec_flux)
    # axs[order//2][order % 2].plot(wavelength, cal_orders)
    axs[order//2][order % 2].plot(wavelength, flux_ratio)
    axs[order//2][order % 2].set_title('order {}'.format(order))
    np.save('smoothed_ratio_{}.npy'.format(order), flux_ratio)
plt.savefig('smoothed_flux_ratio.png')

# Code ends here
