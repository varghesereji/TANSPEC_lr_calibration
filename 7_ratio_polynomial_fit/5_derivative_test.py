# This code is to fit the ratio files with six degree polynonial

from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt


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


def mask_dictionary(pixels):
    masks = {
        0: {1: (pixels < 930) & (pixels > 800),
            2: (pixels > 1750) & (pixels < 1800)},
        1: {1: (pixels < 820) & (pixels > 580),
            2: (pixels > 1080) & (pixels < 1400)},
        2: {1: (pixels < 1400) & (pixels > 400)},
        3: {1: (pixels < 1050) & (pixels > 490),
            2: (pixels > 1135) & (pixels < 1540)},
        4: {1: (pixels < 1600) & (pixels > 320)},
        5: {1: (pixels < 750) & (pixels > 350),
            2: (pixels > 1100) & (pixels < 1700)},
        6: {1: (pixels < 820) & (pixels > 200),
            2: (pixels > 1225) & (pixels < 1730)},
        7: {1: (pixels > 0) & (pixels < 600),
            2: (pixels > 760) & (pixels < 860),
            3: (pixels > 1700) & (pixels < 2000)},
        8: {1: (pixels < 1240) & (pixels > 0),
            2: (pixels > 1990) & (pixels < 2040)},
        9: {1: (pixels < 250) & (pixels > 50),
            2: (pixels > 600) & (pixels < 1700)}
        }
    return masks


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


pdfplots = PdfPages('test_derivative.pdf')

ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')

for order, ratio in enumerate(ratio_file):
    pixels_i = np.arange(0, 2048, 1)
    # mask_dict = mask_dictionary(pixels)
    fig, axs = plt.subplots(2, figsize=(16, 9))
    pixels = pixels_i/1024 - 1
    derivative = np.gradient(ratio)
    smoothed_derivative = smooth(derivative, 150)
    median_smooth_derivative = smoothening_median(derivative, 150)
    integral = integrate.cumtrapz(derivative, pixels_i, initial=pixels_i[0])
    smoothed_integral = integrate.cumtrapz(smoothed_derivative, pixels_i,
                                           initial=pixels_i[0])
    med_smoothed_integral = integrate.cumtrapz(median_smooth_derivative,
                                               pixels_i, initial=pixels_i[0])
    correction_factor = np.mean(ratio - integral)
    integral += correction_factor
    axs[0].plot(pixels, ratio, label='ratio')
    axs[1].plot(pixels, derivative, label='derivative')
    axs[1].plot(pixels, smoothed_derivative,
                label='smooth by convolution')
    axs[1].plot(pixels, median_smooth_derivative, label='smooth by median')
    axs[0].plot(pixels, integral, label='Integral')
    axs[0].plot(pixels, smoothed_integral + correction_factor, '--',
                label='conv deriv integral')
    axs[0].plot(pixels, med_smoothed_integral + correction_factor, '--',
                label='med deriv integral')
    axs[0].legend()
    axs[1].legend()
    # print(mask_dict)
    # for segm in mask_dict[order].values():
    #     pixels_segm = pixels[segm]
    #     smoothed_ratio = smooth(ratio[segm], 35)
    #     pixels_interp = np.linspace(pixels_segm[0], pixels_segm[-1],
    #                                 len(smoothed_ratio))
    #     interp_ratio = interp1d(pixels_interp, smoothed_ratio)
    #     smooth_ratio = interp_ratio(pixels[segm])
    #     plt.plot(pixels_segm, ratio[segm], color='red')
    #     plt.plot(pixels_segm, smooth_ratio, '--', color='green')
    pdfplots.savefig()
pdfplots.close()
