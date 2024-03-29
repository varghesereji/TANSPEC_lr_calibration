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


def telluric_bands():
    '''
    Telluric bands
    ---------------------
    r' : 5550 - 6700
    i' : 6810 - 8050
    Y  : 9700 - 10700
    J  : 11700 - 13300
    H  : 14900 - 17800
    Ks : 19900 - 23100
    H2 : 20970 - 2138
    Brg : 21600 - 21810
    '''
    bands = {
        'r': [5550, 6700],
        'i': [6810, 8050],
        'Y': [9700, 10700],
        'J': [11700, 13300],
        'H': [14900, 17800],
        'Ks': [19900, 23100],
        'H2': [20970, 2138],
        'Brg': [21600, 21810],
        }
    return bands


def filter_wavelength(wl):
    bands = telluric_bands().keys()
    ranges = telluric_bands().values()

    wl_max = max(wl)
    wl_min = min(wl)

    for index, ranges in enumerate(ranges):
        if wl_min > ranges[0] and wl_max < ranges[1]:
            return bands[index]


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


def plotting(star, slit, calspec):
    pdfplots = PdfPages(
        'Images/compare_with_calspec_{}_{}.pdf'.format(star, slit))
    flux_ratio_array = np.load(
        'ratio_files/flux_cal_ratio_hd37725_{}.npy'.format(slit))
    for order in range(0, 10, 1):
        tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
        print(wavelength)
        flux_ratio = flux_ratio_array[order]
        corrected_flux = tanspec_flux / flux_ratio  # flux_ratio / tanspec_flux
        cal_orders = calspec_data(wavelength, star=calspec)

        fig = plt.figure(figsize=(16, 9))
        fig.suptitle('{0} Order {1}'.format(star, order))
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True)
        axs[0].plot(wavelength, tanspec_flux, label='tanspec flux')
        axs[1].plot(wavelength, cal_orders, color='green', label='calspec')
        axs[1].plot(wavelength,
                    corrected_flux, color='red', label='calibrated tanspec')
        # axs[1].set_ylim(min(cal_orders), max(cal_orders))
        # axs[order//2][order % 2].plot(wavelength, cal_orders)
        # axs[order//2][order % 2].plot(wavelength, flux_ratio)
        axs[0].legend()
        axs[1].legend()
        axs[1].set_ylim(min(cal_orders), max(cal_orders))
        # axs[2].legend()
        plt.xlabel('$\lambda \AA$')
        plt.ylabel('Flux')
        pdfplots.savefig()
    pdfplots.close()


star = 'hd37725'   # hd37725, hip41815
slit = 's0.5_by_s4.0'
calspec = 'hd37725'  # hd37725, hd205905

plotting(star, slit, calspec)
# End
