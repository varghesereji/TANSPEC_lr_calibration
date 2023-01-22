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
        'H2': [20970, 21380],
        'Brg': [21600, 21810],
        }
    return bands


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


def filter_wavelength(wl):
    bands = telluric_bands().keys()
    ranges = telluric_bands().values()
    print(bands)
    wl_max = max(wl)
    wl_min = min(wl)
    spectra_band = {}
    for index, ranges in enumerate(ranges):
        print(index, wl_min, ranges, wl_max)
        if wl_min < ranges[0] and wl_max > ranges[1]:
            print(list(bands)[index])
            band_mask = (ranges[0] < wl) & (wl < ranges[1])
            spectra_band[list(bands)[index]] = band_mask, ranges
    return spectra_band


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
    pdfplots = PdfPages('./Images/{0}_{1}_filtered.pdf'.format(star, slit))
    flux_ratio_array = np.load(
        'ratio_files/flux_cal_ratio_hd37725_{}.npy'.format(slit))

    for order in range(0, 10, 1):
        tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
        flux_ratio = flux_ratio_array[order]
        corrected_flux = tanspec_flux / flux_ratio  # flux_ratio / tanspec_flux
        cal_orders = calspec_data(wavelength, star=calspec)
        print(order)
        filter_dict = filter_wavelength(wavelength)
        print(filter_dict, 'going to plot')
        if filter_dict != {}:
            band_name = list(filter_dict.keys())
            range_and_mask = list(filter_dict.values())

            for q, band in enumerate(band_name):
                fig = plt.figure(figsize=(16, 9))
                fig.suptitle('{0} Order {1} band {2}'.format(star, order, band))
                #gs = fig.add_gridspec(1, hspace=0)
                #axs = gs.subplots(sharex=True)
                print('bands:', band)
                print('values:', range_and_mask[q])
                mask, ranges = range_and_mask[q]
                print('mask', mask)
                # axs[0].plot(wavelength[mask], tanspec_flux[mask],
                #          label='tanspec flux')
                # axs[1].plot(wavelength[mask], flux_ratio[mask], label='ratio')
                # axs[0].plot([ranges[0], ranges[0]], [-1, 1000], '--',
                #             color='black')
                # axs[0].plot([ranges[1], ranges[1]], [-1, 1000], '--',
                #             color='black')
                # axs[0].text((ranges[0]+ranges[1])/2, 0, band, fontsize=22)
                plt.plot(wavelength[mask], cal_orders[mask], color='green',
                         label='calspec')
                plt.plot(wavelength[mask],
                         corrected_flux[mask], color='red',
                         label='calibrated tanspec')
                for wl in HYDROGEN:
                    plt.axvline(x=wl, color='k', linestyle='--')
                plt.legend()
                plt.ylim(min(cal_orders[mask]), max(cal_orders[mask]))
                plt.xlim(min(wavelength[mask]), max(wavelength[mask]))
                # axs[0].legend()
                # axs[2].legend()
                # axs[3].legend()
            plt.xlabel('$\lambda \AA$')
            plt.ylabel('Flux')
            pdfplots.savefig()

    pdfplots.close()
    return 0


star = 'hd37725'   # hd37725, hip41815
slit = 's0.5_by_s4.0'
calspec = 'hd37725'  # hd37725, hd205905

plotting(star, slit, calspec)


# End
