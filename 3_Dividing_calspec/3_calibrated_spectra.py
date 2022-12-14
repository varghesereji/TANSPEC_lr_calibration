
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
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


pdfplots = PdfPages('Images/Cleaned_sigmal.pdf')
# fig, axs = plt.subplots(1, 2, figsize=(16, 9))
for order in range(0, 10, 1):
    fig, axs = plt.subplots(2, figsize=(16, 9))
    tanspec_flux, wavelength = tanspec_data(order)
    flux_ratio = np.load('ratio_files/smoothed_ratio_{}.npy'.format(order))
    corrected_flux = flux_ratio * tanspec_flux
    axs[0].plot(wavelength, tanspec_flux, label='tanspec flux')
    axs[1].plot(wavelength, corrected_flux, color='red', label='calibrated_flux')
    # axs[order//2][order % 2].plot(wavelength, cal_orders)
    # axs[order//2][order % 2].plot(wavelength, flux_ratio)
    axs[0].legend()
    axs[1].legend()
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('Flux')
    pdfplots.savefig()
pdfplots.close()

# Code ends here
