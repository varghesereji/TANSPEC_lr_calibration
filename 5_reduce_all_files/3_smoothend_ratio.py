# This code is to Filter out the absorption lines in the spectra

from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy import interpolate
from scipy import ndimage
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import re


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


def smoothening_mean(flux, size=20):
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
    smoothend_flux = ndimage.uniform_filter1d(flux, size)
    return smoothend_flux


def find_file(dirc):
    '''
    Parameters
    --------------------------------
    dirc: path to directory
    ================================
    filename: Name of the file
    --------------------------------
    '''
    files_list = os.listdir(dirc)
    for filename in files_list:
        match = re.search("\.final\.fits", filename)
        if match:
            return filename


def tanspec_data(star='hd37725', slit='s0.5', q=0,):
    '''
    Parameters
    --------------------------------
    star: string. Name of star, lower case.
    slit: string. Name of slit.
    q: order of spectra starting from 0
    ================================
    Return
    --------------------------------
    order: Spectra from reduced tanspec data for given order
    wavelength: Array of wavelength.
    '''
    data_dir_path = '/home/varghese/Desktop/DP2_TANSPEC/\
Data/Output_spectra'
    file_path = os.path.join(data_dir_path, star, slit, 'tanspec_slopes')
    data_file = find_file(file_path)
    flux = fits.getdata(os.path.join(
        file_path, data_file),
        ext=0)
    order = flux[q]
    wl = fits.getdata(
        os.path.join(file_path,
                     data_file), ext=1)
    wavelength = wl[q]
    return order, wavelength


def calspec_data(wavelength, star='hd37725'):
    '''
    Parameters
    --------------------------------
    wavelength: array or float, values of wavelength that you wanted the flux
    star: calspec file of which star
    ================================
    Return
    --------------------------------
    interpolated_flux: Interpolated flux data.
    '''
    path = '/home/varghese/Desktop/DP2_TANSPEC/Calspec_data/'
    file_name = star + '_mod_004.fits'
    data_file = Table.read(os.path.join(path, file_name))
    calspec_wl = np.array(data_file['WAVELENGTH'])
    calspec_flux = np.array(data_file['FLUX'])
    interpolated = interpolate.interp1d(calspec_wl, calspec_flux)
    interpolated_flux = interpolated(wavelength)
    return interpolated_flux


star = 'hd37725'
slit = 's4.0'
calspec = 'hd37725'

pdfplots = PdfPages('Images/Finding_ratio_{0}_{1}_images.pdf'.format(star, slit))

smoothed_ratio = None
for order in range(0, 10, 1):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle('Order {}'.format(order))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
    cal_order = calspec_data(wavelength, calspec)
    ratio = tanspec_flux / cal_order
    clipped_ratio = smoothening_mean(smoothening_median(ratio))
    axs[0].plot(wavelength, cal_order)
    axs[1].plot(wavelength, tanspec_flux)
    axs[2].plot(wavelength, clipped_ratio)
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('ratio')
    pdfplots.savefig()
    if smoothed_ratio is None:
        smoothed_ratio = clipped_ratio
    else:
        smoothed_ratio = np.vstack((smoothed_ratio, clipped_ratio))

np.save('ratio_files/flux_cal_ratio_{0}_{1}.npy'.format(star, slit),
        smoothed_ratio)
pdfplots.close()


# End
