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


HYDROGEN = [
    8440.3,
    8469.6,
    8504.8,
    8547.7,
    8600.8,
    8667.4,
    8752.9,
    8865.2,
    9017.4,
    9231.5,
    9548.6,
    10052.1,
    10941.1,
    12821.6,
    15137.4,
    15196.,
    15264.7,
    15346.,
    15443.1,
    15560.7,
    15705.,
    15884.9,
    16113.7,
    16411.7,
    16811.1,
    17366.9,
    18179.1,
    18756.1,
    19450.9,
    21661.2,
    23297.9,
    23329.6,
    23364.4,
    23402.8,
    23445.3,
    23492.4,
    23544.8,
    23603.5,
    23669.4,
    23743.8,
    23828.2,
    23924.8,
    24035.5,
    24163.9,
    24313.6,
    24490.,
    24699.9,
    24952.5,
    25260.9,
    25643.3,
    26126.5,
    26258.7,
    26751.3,
    27582.7,
]


star = 'hd37725'   # hd37725, hip41815
slit = 's0.5'
calspec = 'hd37725'  # hd37725, hd205905

pdfplots = PdfPages('Images/compare_with_calspec_{}.pdf'.format(slit))

flux_ratio_array = np.load(
    'ratio_files/flux_cal_ratio_hd37725_{}.npy'.format(slit))

for order in range(0, 10, 1):
    tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
    print(wavelength)
    flux_ratio = flux_ratio_array[order]
    corrected_flux = tanspec_flux / flux_ratio  # flux_ratio / tanspec_flux
    smoothend_flux = smoothening_median(corrected_flux)
    cal_orders = calspec_data(wavelength, star=calspec)

    fig = plt.figure(figsize=(16, 4))
    fig.suptitle('{0} Order {1}'.format(star, order))
    plt.plot(wavelength, cal_orders, color='green', label='calspec')
    plt.plot(wavelength,
             corrected_flux, color='red', label='calibrated tanspec')
    for wl in HYDROGEN:
        plt.axvline(x=wl, color='k', linestyle='--')
    plt.legend()
    plt.ylim(min(cal_orders), max(cal_orders))
    plt.xlim(min(wavelength), max(wavelength))

    plt.xlabel('$\lambda \AA$')
    plt.ylabel('Flux')
    pdfplots.savefig()
    np.save('Data_Files/corrected_flux_{}_{}.npy'.format(order, slit),
            corrected_flux)
    np.save('Data_Files/wavelength_{}_{}.npy'.format(order, slit), wavelength)
    np.save('Data_Files/cal_orders_{}_{}.npy'.format(order, slit), cal_orders)
pdfplots.close()

np.save('Data_Files/hydrogen_lines.npy', np.array(HYDROGEN))

# End