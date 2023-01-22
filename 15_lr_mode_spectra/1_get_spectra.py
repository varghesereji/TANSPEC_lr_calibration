import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os


def return_filename(star, slit, source):
    project_path = '/home/varghese/Desktop/DP2_TANSPEC/'
    file_path = 'lr_background_reduction/spectrum_extraction/extracted_spectra'

    files_dictionary = {
        'hd37725': {
            's0.5': {
                'star': 'Ext_hd37725_lr_s0.5hg.fits',
                'argon': 'Ext_argon_lr_s0.5hg.fits',
                'neon': 'Ext_neon_lr_s0.5hg.fits'
                },
            's4.0': {
                'star': 'Ext_hd37725_lr_s4.0_hg.fits',
                'argon': 'Ext_argon_lr_s4.0_hg.fits',
                'neon': 'Ext_neon_lr_s4.0_hg.fits'
                }
            },
        'hip41815': {
            's0.5': {
                'star': 'Ext_hip41815_lr_s0.5_hg.fits',
                'argon': 'Ext_argon_lr_s0.5_hg.fits',
                'neon': 'Ext_neon_lr_s0.5_hg.fits'
                },
            's1.0': {
                'star': 'Ext_hip41815_lr_s1.0_hg.fits',
                'argon': 'Ext_argon_lr_s1.0_hg.fits',
                'neon': 'Ext_neon_lr_s1.0_hg.fits'
                }
            }

    }
    print(os.listdir(os.path.join(project_path, file_path)))
    filename = files_dictionary[star][slit][source]
    full_path = os.path.join(project_path, file_path, filename)
    return full_path


def tanspec_data(q=0):
    file_path = '/home/varghese/Desktop/DP2_TANSPEC/\
Data/Output_spectra/hd37725/s0.5_by_s4.0/tanspec_slopes'
    flux = fits.getdata(os.path.join(
        file_path, 'op_hd37725_s0.5_by_s4.0_grating1_A.ms.wlc.final.fits'),
        ext=0)
    order = flux[q]
    wl = fits.getdata(
        os.path.join(file_path,
                     'op_hd37725_s0.5_by_s4.0_grating1_A.ms.wlc.final.fits'),
        ext=1)
    wavelength = wl[q]
    return order, wavelength


star = 'hd37725'
slit = 's0.5'
source = 'star'


filename = return_filename(star, slit, source)
source_spectra = fits.getdata(filename, ext=0)
bkg_1 = fits.getdata(filename, ext=1)
bkg_2 = fits.getdata(filename, ext=2)
source_spectra = source_spectra[0] - (bkg_1[0] + bkg_2[0])/2
fig, axs = plt.subplots(2, figsize=(16, 9))
axs[0].plot(source_spectra, color='blue')
for i in range(10):
    xd_flux, wavelength = tanspec_data(i)
    axs[1].plot(wavelength, xd_flux, color='blue')
plt.show()
# End
