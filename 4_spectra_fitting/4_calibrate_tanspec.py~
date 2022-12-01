# This code is to Filter out the absorption lines in the spectra

from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy import interpolate
from scipy import ndimage
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os


def smoothening_rank(flux, rank=-1, size=20):
    '''
    Parameters
    --------------------------------
    flux: The input flux that have to smoothen
    rank, size: defined in the function scipy.ndimage.rank_filters
    ================================
    Return
    --------------------------------
    smoothend_flux: The spectra which absorption spectra are removed
    '''
    smoothend_flux = ndimage.rank_filter(flux, rank, size)
    return smoothend_flux


def smoothening_median(flux, size=20):
    '''
    Parameters
    --------------------------------
    flux: The input flux that have to smoothen
    rank, size: defined in the function scipy.ndimage.rank_filters
    ================================
    Return
    --------------------------------
    smoothend_flux: The spectra which absorption spectra are removed
    '''
    smoothend_flux = ndimage.median_filter(flux, size)
    return smoothend_flux


def clipping(signal, high, low):
    return np.clip(signal, low, high)


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
smoothend_flux = smoothening_rank(calspec_flux)
interpolation_smoothed = interpolate.interp1d(calspec_wl, smoothend_flux)
interpolation_nsm = interpolate.interp1d(calspec_wl, calspec_flux)

pdfplots = PdfPages('Images/test_images.pdf')

for order in range(0, 10, 1):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle('Order {}'.format(order))
    gs = fig.add_gridspec(5, hspace=0)
    axs = gs.subplots(sharex=True)
    tanspec_flux, wavelength = tanspec_data(order)
    cal_orders_nsm = interpolation_nsm(wavelength)
    cal_orders_sm = interpolation_smoothed(wavelength)
    tanspec_sm = smoothening_median(tanspec_flux)
    ratio = cal_orders_sm / tanspec_sm
    clipped_ratio = smoothening_median(smoothening_median(ratio))
    axs[0].plot(wavelength, cal_orders_nsm)
    axs[1].plot(wavelength, cal_orders_sm)
    axs[2].plot(wavelength, tanspec_flux)
    axs[3].plot(wavelength, tanspec_sm)
    axs[4].plot(wavelength, clipped_ratio)
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('ratio')
    pdfplots.savefig()
    np.save('ratio_files/smoothed_ratio_{}.npy'.format(order), clipped_ratio)
pdfplots.close()


# End
