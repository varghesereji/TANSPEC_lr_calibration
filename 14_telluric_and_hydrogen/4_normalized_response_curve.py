from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from astropy.table import Table
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


def interpolate_pixels(pixels, wavelength, wl):
    interp = interp1d(wavelength, pixels, fill_value='extrapolate')
    wl_position = interp(wl)
    return wl_position


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


def telluric_regions():
    telluric = np.arange(3000, 5550, 1)
    telluric = np.concatenate((telluric,
                               np.arange(6700, 6810, 1)))
    telluric = np.concatenate((telluric,
                               np.arange(8050, 9700, 1)))
    telluric = np.concatenate((telluric,
                               np.arange(10700, 11700, 1)))
    telluric = np.concatenate((telluric,
                               np.arange(13300, 14900, 1)))
    telluric = np.concatenate((telluric,
                               np.arange(17800, 19900, 1)))
    print(telluric)
    return telluric


'''
telluric_mask = (wavelength > 5550) & (wavelength < 6700)
    | (wavelength > 6810) & (wavelength < 8050)
    | (wavelength > 9700) & (wavelength < 10700)
    | (wavelength > 11700) & (wavelength < 13300)
    | (wavelength > 14900) & (wavelength < 17800)
    | (wavelength > 19900) & (wavelength < 23100)
    | (wavelength > 20970) & (wavelength < 2138)
    | (wavelength > 21600) & (wavelength < 21810)
'''

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


project_path = '/home/varghese/Desktop/DP2_TANSPEC'
response_path = os.path.join(project_path, 'Data/Instrument_response/smoothed')
mask = np.load(os.path.join(response_path, 'Trim_spectrum.npy'))

# file = open('atmabs.txt', 'r')
# telluric_lines = []
# for line in file:
#     stripped_line = line.strip()
#     telluric_lines.append(float(stripped_line))

# star = 'hip41815'
# slit = 's1.0'
# calspec = 'hd205905'

star = 'hd37725'
slit = 's4.0'
calspec = 'hd37725'

ratio_file = np.load(os.path.join(response_path,
                                  'Response_function_{}.npy'.format(slit)))
# ratio_file = ratio_file / np.median(ratio_file)
scaling_factor = np.median(ratio_file)
pdfplots = PdfPages('normalized_ratio_{0}_{1}.pdf'.format(star,
                                                          slit)
                    )


file = open('atmabs.txt', 'r')
telluric_lines = []
for line in file:
    stripped_line = line.strip()
    telluric_lines.append(float(stripped_line))


with open('telluric_lines/mktrans_zm_30_10.txt') as f:
    polyshape = []
    for line in f:
        line = line.split()
        if line:
            line = [float(i) for i in line]
            polyshape.append(line)

polyshape = np.array(polyshape)

for order, ratio in enumerate(ratio_file):
    tanspec_flux, wavelength = tanspec_data(star=star, slit=slit, q=order)
    mask_portion = mask[order]
    pixels_actual = np.arange(0, 2048, 1)
    fig, ax1 = plt.subplots(figsize=(16, 9))
    pixels = pixels_actual/1024 - 1
    master_mask = mask[order]
    # telluric_positions = telluric_regions()
    calib_tanspec = tanspec_flux / (ratio)
    color = 'red'
    ax1.set_xlabel('wavelength $\AA$')
    ax1.set_ylabel('TANSPEC Flux', color=color)
    # for windex in telluric_positions:
    #     ax[0].axvline(x=windex, color='black', linestyle='-', alpha=0.1)
    #     ax[1].axvline(x=windex, color='black', linestyle='-', alpha=0.1)
    # for wl in telluric_lines:
    for wl in HYDROGEN:
        ax1.axvline(x=wl, color='blue', linestyle='--')
    ax1.plot(wavelength[mask_portion],
             calib_tanspec[mask_portion] * scaling_factor,
             color='red')
    ax1.tick_params(axis='y', labelcolor=color)
    ax2 = ax1.twinx()
    color = 'green'
    cal_orders = calspec_data(wavelength, star=calspec)
    # scaling_factor = np.median(calib_tanspec) / np.median(cal_orders)
    # cal_orders = cal_orders * scaling_factor
    ax2.plot(wavelength[mask_portion], cal_orders[mask_portion],
             color='green')
    ax1.vlines(x=10000 * polyshape[:, 0], colors='black',
               ymin=0,
               ymax=(1-polyshape[:, 1]) * max(cal_orders) * scaling_factor,
               alpha=0.3)

    ax2.set_ylabel('CALSPEC Flux', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim(min(cal_orders), max(cal_orders))
    ax2.set_xlim(min(wavelength[mask_portion]),
                 max(wavelength[mask_portion]))
    ax1.set_xlim(min(wavelength[mask_portion]),
                 max(wavelength[mask_portion]))
    ax1.set_ylim(min(cal_orders) * scaling_factor,
                 max(cal_orders) * scaling_factor,)
    ax1.vlines(x=telluric_lines, colors='black', linestyle=':', alpha=0.2,
               ymin=min(cal_orders) * scaling_factor,
               ymax=max(cal_orders) * scaling_factor)
    pdfplots.savefig()
    print(min(wavelength[mask_portion]), max(wavelength[mask_portion]))
pdfplots.close()


# End
