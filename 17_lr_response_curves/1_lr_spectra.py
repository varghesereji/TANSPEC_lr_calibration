import matplotlib.pyplot as plt
from astropy.io import fits
import os
import numpy as np

project = '/home/varghese/Desktop/DP2_TANSPEC'

# lr mode spectra
data_path = 'lr_background_reduction/spectrum_extraction/extracted_spectra'
file_path = os.path.join(project, data_path, 'Ext_hd37725_lr_s0.5hg.fits')
fitted_wavelength = np.load('wavelength_calibrated.npy')

star = fits.getdata(file_path, ext=0)[0]
sky = fits.getdata(file_path, ext=1)[0]

star_sky_sbtr = star - sky
mask = fitted_wavelength > 0

plt.figure(figsize=(16,9))
plt.plot(fitted_wavelength[mask], star[mask])
plt.plot(fitted_wavelength[mask], sky[mask])
plt.plot(fitted_wavelength[mask], star_sky_sbtr[mask])
plt.show()
