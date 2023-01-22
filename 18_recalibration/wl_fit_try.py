import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import least_squares
import csv


def polynomial(p, a0, a1, a2, a3):
    return a0 * p**3 + a1 * p**2 + a2 * p + a3


wls = []
pixels = []

argon = []
argon_pixels = []
with open('wl_matching/argon.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        argon.append(float(lines[0]))
        argon_pixels.append(float(lines[1]))

neon = []
neon_pixels = []
with open('wl_matching/neon.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        neon.append(float(lines[0]))
        neon_pixels.append(float(lines[1]))

sky = []
sky_pixels = []
with open('wl_matching/sky.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        sky.append(float(lines[0]))
        sky_pixels.append(float(lines[1]))

telluric = []
telluric_pixels = []
with open('wl_matching/telluric.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        telluric.append(float(lines[0]))
        telluric_pixels.append(float(lines[1]))

# tanspec_flux = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/Data/Output_spectra/hd37725/s4.0/tanspec_slopes/op_hd37725_s4.0_grating1_A.ms.wlc.fits', ext=0)
# bkg = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/Data/Output_spectra/hd37725/s4.0/tanspec_slopes/op_hd37725_s4.0_grating1_A.ms.wlc.fits', ext=1)
# tanspec_flux = tanspec_flux - bkg
# wavelength = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/Data/Output_spectra/hd37725/s4.0/tanspec_slopes/op_hd37725_s4.0_grating1_A.ms.wlc.fits', ext=3)

# lr_spectra = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/lr_background_reduction/spectrum_extraction/extracted_spectra/Ext_hd37725_lr_s0.5hg.fits', ext=0)
# Sky = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/lr_background_reduction/spectrum_extraction/extracted_spectra/Ext_hd37725_lr_s0.5hg.fits', ext=1)

pixels = np.array(pixels)
wl = np.array(wls)

full_pixels = np.arange(700, 1600, 1)

fitting = np.polyfit(pixels, wl, 3)
curve = np.poly1d(fitting)
fitted_wls = curve(full_pixels)


mask = (np.gradient(fitted_wls) < 0)
neg_der_points = full_pixels[mask]
der_value = np.ones(np.sum(mask))*1e70
print(neg_der_points.size, der_value.size)
print(wl)

all_pixels = np.arange(0, 2048, 1)
mask = (all_pixels > 699) & (all_pixels < 1600)
wl_array = np.zeros(2048)
wl_array[mask] = fitted_wls
plt.figure(figsize=(16, 9))
# fig, axs = plt.subplots(3, figsize=(16, 9))
# axs[0].plot(np.arange(800, 1600, 1), lr_spectra[0, 800:1600], color='red', label='Star LR spectra extracted')
# axs[0].plot(np.arange(800, 1600, 1), Sky[0, 800:1600], color='black', label='Sky LR spectra extracted')
# axs[1].plot(wavelength[0], tanspec_flux[0], color='red', label='Star XD spectra')
# axs[1].plot(wavelength[0], bkg[0], color='black', label='Sky XD spectra')
# _ = axs[1].plot(wavelength[1:].T, tanspec_flux[1:].T, color='red')
# _ = axs[1].plot(wavelength[1:].T, bkg[1:].T, color='black')
plt.plot(full_pixels, fitted_wls)
plt.plot(argon_pixels, argon, '<', label='argon')
plt.plot(neon_pixels, neon, 'o', label='neon')
plt.plot(sky_pixels, sky, 's', label='sky')
plt.plot(telluric_pixels, telluric, '*', label='telluric')
plt.xlabel('pixels')
plt.ylabel('wavelength $\AA$')
plt.legend()
# axs[1].legend()
# axs[0].legend()
# print(np.max(tanspec_flux))
# axs[0].set_xlabel('pixels')
# axs[0].set_ylabel('flux')
# axs[1].set_xlabel('wavelength $\AA$')
# axs[1].set_ylabel('flux')
# axs[1].set_ylim(-2, np.max(tanspec_flux))
# np.save('wavelength_calibrated.npy', wl_array)
# plt.savefig('/home/varghese/MEGAsync/DP2_report/Images/wl_matching.png')
plt.show()
