
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from astropy.table import Table
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
interpolation = interpolate.interp1d(calspec_wl, calspec_flux)


pdfplots = PdfPages('Images/compare_with_calspec.pdf')
# fig, axs = plt.subplots(1, 2, figsize=(16, 9))
for order in range(0, 10, 1):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle('Order {}'.format(order))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    tanspec_flux, wavelength = tanspec_data(order)
    flux_ratio = np.load('ratio_files/smoothed_ratio_{}.npy'.format(order))
    corrected_flux = flux_ratio * tanspec_flux
    cal_orders = interpolation(wavelength)
    axs[0].plot(wavelength, tanspec_flux, label='tanspec flux')
    axs[2].plot(wavelength, cal_orders, color='green', label='calspec')
    axs[1].plot(wavelength,
                corrected_flux, color='red', label='calibrated tanspec')
    # axs[order//2][order % 2].plot(wavelength, cal_orders)
    # axs[order//2][order % 2].plot(wavelength, flux_ratio)
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('Flux')
    pdfplots.savefig()
pdfplots.close()

# Code ends here
