import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import fits
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


def scale_pixels(pixels):
    return pixels/1024 - 1


def reverse_scaling(sc_pixels):
    return (sc_pixels+1) * 1024


def add_more_points(ratio, actual_pixels,
                    first_pixel, last_pixel):
    added_pixels = np.linspace(first_pixel, last_pixel, 1000)
    ratio_interp = interpolate.interp1d(
        scale_pixels(actual_pixels[first_pixel-10:last_pixel+10]),
        ratio[first_pixel-10:last_pixel+10])
    ratio_more_points = ratio_interp(scale_pixels(added_pixels))
    return ratio_more_points, added_pixels


order = 1
ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')
ratio = ratio_file[order]
mask_full = np.load('Mask_for_Instrument_response.npy')
pixels_actual = np.arange(0, 2048, 1)
scaled_pixels = scale_pixels(pixels_actual)


mask = (pixels_actual < 500) | (pixels_actual > 1500)
    # | (pixels_actual > 850) & (pixels_actual < 1080)
masked_reg = ~mask

median_ratio = np.median(ratio)
ratio = ratio / median_ratio

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 850, 855)
ratio_to_fit = np.concatenate((ratio[~mask], added_ratio))
pixels_to_fit = np.concatenate((scaled_pixels[~mask],
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 822, 830)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1030, 1055)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1025, 1027)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1385, 1390)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1334, 1338)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1360, 1362)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1330, 1335)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1140, 1145)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 755, 760)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 590, 595)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 720, 725)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 680, 682)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 700, 702)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 810, 812)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 920, 922)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 640, 660)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 940, 950)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1050, 1070)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1210, 1215)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 798, 800)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1480, 1482)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1260, 1262)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1290, 1292)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1460, 1455)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 895, 900)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 565, 570)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))
added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1600, 1602)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))









project_path = '/home/varghese/Desktop/DP2_TANSPEC'
response_path = os.path.join(project_path, 'Data/Instrument_response/smoothed')
mask_trim = np.load(os.path.join(response_path, 'Trim_spectrum.npy'))
mask_portion = mask_trim[order]

back_scaled_pixels = reverse_scaling(pixels_to_fit)
fig, ax = plt.subplots(2, figsize=(16, 9))
diff_mask = np.diff(mask)
positions = list(np.where(diff_mask)[0])
mask_positions = list(np.where(~mask)[0])
for x_pos in mask_positions:
    ax[0].axvline(x=x_pos, color='yellow', linestyle='-', alpha=0.1)
for mark in positions:
    ax[0].text(mark, 0, str(mark+1), fontsize=15)
ax[0].plot(pixels_actual, ratio, '-',)
polyfit = np.polyfit(pixels_to_fit,
                     ratio_to_fit, 7)
fitted = np.poly1d(polyfit)
ax[0].plot(pixels_actual[masked_reg], ratio[masked_reg], color='red')
ax[0].plot(back_scaled_pixels, ratio_to_fit, '.', color='black')
ax[0].plot(pixels_actual, fitted(scaled_pixels),
           '--', color='green')
ax[0].set_xlabel('pixels')
ax[0].set_ylim(-0.1, 14)
ax[0].set_xlim(min(pixels_actual[mask_portion]),
               max(pixels_actual[mask_portion]))
ax[1].set_xlim(min(pixels_actual[mask_portion]),
               max(pixels_actual[mask_portion]))

tanspec_flux, wavelength = tanspec_data(star='hd37725', slit='s4.0', q=order)
calib_tanspec = tanspec_flux /(median_ratio * ratio)  # fitted(scaled_pixels))   # fitted(pixels)
cal_orders = calspec_data(wavelength, star='hd37725')
ax[1].plot(pixels_actual[mask_portion],
           calib_tanspec[mask_portion],
           color='red')
ax[1].plot(pixels_actual[mask_portion], cal_orders[mask_portion],
           color='green')
ax[1].set_ylim(min(cal_orders), max(cal_orders))
plt.show()
# np.save('instrum_resp/Instrument_resp_1.npy', ratio)
print('done', positions)

'''
Page 2 masked_red = (pixels_actual < 500) \
        | (pixels_acual > 850) & (pixels_actual < 1080)

Page 3 masked_reg = (pixels_actual < 370) \
        | (pixels_actual > 1500) \
        | (pixels_actual > 1070) & (pixels_actual < 1110)

Page 4    masked_reg = (pixels_actual < 300) \
        | (pixels_actual > 1045) & (pixels_actual < 1130) \
        | (pixels_actual > 1700)

Page 5     masked_reg = (pixels_actual < 290) \
        | (pixels_actual > 1740) \
        | (pixels_actual > 800) & (pixels_actual < 870)

Page 6 
Page 8    masked_reg = (pixels_actual > 600) & (pixels_actual < 760) \
        | (pixels_actual > 860) & (pixels_actual < 1701)

Page 9 masked_reg = (pixels_actual > 1240) & (pixels_actual < 1990) \
        | (pixels_actual > 920) & (pixels_actual < 1181) \
        | (pixels_actual > 255) & (pixels_actual < 310) \
        | (pixels_actual > 632) & (pixels_actual < 760)

Page 10.     masked_reg = (pixels_actual > 355) & (pixels_actual < 700) \
        | (pixels_actual > -1) & (pixels_actual < 137) \
        | (pixels_actual > 1490) & (pixels_actual < 1530) \
        | (pixels_actual > 1459) & (pixels_actual < 1486) \
        | (pixels_actual > 1388) & (pixels_actual < 1450) \
        | (pixels_actual > 1551) & (pixels_actual < 1635) \
        | (pixels_actual > 1296) & (pixels_actual < 1377) \
        | (pixels_actual > 1213) & (pixels_actual < 1285) \
        | (pixels_actual > 865) & (pixels_actual < 914) \
        | (pixels_actual > 1164) & (pixels_actual < 1200) \
        | (pixels_actual > 1103) & (pixels_actual < 1160) \
        | (pixels_actual > 1712)

'''

# End

