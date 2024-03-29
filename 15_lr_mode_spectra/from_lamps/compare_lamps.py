import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from astropy.io import fits
from scipy.optimize import curve_fit
import os


project = '/home/varghese/Desktop/DP2_TANSPEC'
# science_files = 'Data/Output_spectra/hd37725/s4.0/tanspec_slopes'
# path = os.path.join(project, science_files)


# lr mode spectra

# argon
data_path = 'lr_background_reduction/spectrum_extraction/extracted_spectra'
file_path = os.path.join(project, data_path, 'Ext_argon_lr_s0.5hg.fits')
fitted_wavelength = np.load('../wavelength_calibrated.npy')
argon = fits.getdata(file_path, ext=0)[0]
arg_sky = fits.getdata(file_path, ext=1)[0]
argon1 = argon - arg_sky


argon_data = np.load('argon_data.npy')
argon_visible = np.load('argon_full.npy')

# neon
file_path = os.path.join(project, data_path, 'Ext_neon_lr_s0.5hg.fits')
neon = fits.getdata(file_path, ext=0)[0]
neon_sky = fits.getdata(file_path, ext=1)[0]

neon_data = np.load('neon_data.npy')

# params, err = curve_fit(bb_spectra, wl_si, bg1, p0=np.array([1160, 10**(-6)]))
# print(params)
# T, scale = params
# fig, axs = plt.subplots(2, figsize=(16, 9))
# # axs[0].plot(fitted_wavelength, argon)
# axs[0].plot(fitted_wavelength, argon1)
# axs[1].plot(argon_data[0], argon_data[1])
# axs[1].plot(argon_visible[0], argon_visible[1])

fig, axs = plt.subplots(2, figsize=(16, 9))
axs[0].plot(fitted_wavelength, argon)
axs[0].plot(argon_data[0], argon_data[1] * np.max(argon) / np.max(argon_data[1]),
            color='red')
axs[0].plot(argon_visible[0],
            argon_visible[1] * np.max(argon) / np.max(argon_visible[1]),
            color='red')
axs[1].plot(fitted_wavelength, neon)
axs[1].plot(neon_data[0], neon_data[1] * np.max(neon) / np.max(neon_data[1]),
            color='red')
plt.show()


# End
