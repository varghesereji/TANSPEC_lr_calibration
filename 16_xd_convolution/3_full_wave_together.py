from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.signal import convolve
from astropy.io import fits
import os
import re


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


def gaussian(mu, sigma, lmin, lmax, N):
    print(mu, sigma)
    lambdas = np.linspace(lmin, lmax, N)
    normalization = 1 / (sigma * np.sqrt(2*np.pi))
    arguement = (lambdas - mu)**2 / (2 * sigma**2)
    func = np.exp(-arguement)
    return func


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


project_path = '/home/varghese/Desktop/DP2_TANSPEC'
response_path = os.path.join(project_path, 'Data/Instrument_response/smoothed')
mask = np.load(os.path.join(response_path, 'Trim_spectrum.npy'))

star = 'hd37725'
slit = 's4.0'
calspec = 'hd37725'

ratio_file = np.load(os.path.join(response_path,
                                  'Normalized_response_curve.npy'))

flux_wavelength = None

resolutions = [110, 105, 105, 100, 100, 110, 120, 140, 160, 300]
plt.figure()
for order, ratio in enumerate(ratio_file):
    print(order)
    tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
    mask_portion = mask[order]
    tanspec_flux = tanspec_flux[mask_portion] / ratio[mask_portion]
    wavelength = wavelength[mask_portion]
    # axs[0].plot(wavelength, tanspec_flux)
    tanspec_flux = smoothening_median(tanspec_flux, 10)
    N = wavelength.size
    mu = np.mean(wavelength)
    sigma = mu / (resolutions[order] * 2.355)
    lmin = min(wavelength)
    lmax = max(wavelength)
    gaussian_func = gaussian(mu, sigma, lmin, lmax, N)
    convolved_spectra = convolve(tanspec_flux, gaussian_func, mode='same')
    # axs[0].plot(wavelength, tanspec_flux)
    # axs[1].plot(wavelength, gaussian_func)
    plt.plot(wavelength, convolved_spectra)
    np.save('flux_wavelength_{}.npy'.format(order),
            np.vstack((convolved_spectra, wavelength)))
plt.show()


# End
