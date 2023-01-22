from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import re


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


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


star = 'hd37725'
slit = 's0.5'
# slit = 's4.0'
mask_full = np.load('Mask_for_Instrument_response.npy')
response = np.load('Response_function_{}.npy'.format(slit))
trimmer = np.load('Trim_spectrum.npy')
pixels_actual = np.arange(0, 2048, 1)
pixels = pixels_actual/1024 - 1
pdfplots = PdfPages('response_fn_{}.pdf'.format(slit))

plt.figure(figsize=(16, 9))
for order, ratio in enumerate(response):
    trim = trimmer[order]
    tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
    plt.plot(wavelength[trim], ratio[trim])

plt.savefig('trimmed_response_{}.png'.format(slit))

# End
