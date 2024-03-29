import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from astropy.io import fits
from scipy.optimize import curve_fit
import os


project = '/home/varghese/Desktop/DP2_TANSPEC'
science_files = 'Data/Output_spectra/hd37725/s0.5_by_s4.0/tanspec_slopes'
path = os.path.join(project, science_files)


def bb_spectra(lam, T, scale):
    a = (np.exp(const.h * const.c / (lam * const.k * T)) - 1)
    b = 2 * const.h * (const.c)**2 / (lam**5)
    bb_eq = scale * b * (1/a)
    return bb_eq


science_file = os.path.join(
    path, 'op_hd37725_s0.5_by_s4.0_grating1_A.ms.wlc.fits')
bg1 = fits.getdata(science_file, ext=1)
wl_sc = fits.getdata(science_file, ext=3)

fmax = np.max(bg1)

mask = bg1 == fmax
wl_max = wl_sc[mask]
print(wl_max)

wl_si = wl_sc * 10**(-10)

# lr mode spectra

sky = fits.getdata('/home/varghese/Desktop/DP2_TANSPEC/lr_background_reduction/spectrum_extraction/extracted_spectra/Ext_hd37725_lr_s0.5hg.fits',
                   ext=1)

# params, err = curve_fit(bb_spectra, wl_si, bg1, p0=np.array([1160, 10**(-6)]))
# print(params)
# T, scale = params
fig, ax = plt.subplots(2, figsize=(16,9))

_ = ax[0].plot(wl_sc.T, bg1.T, color='blue')
ax[1].plot(sky[0], color='blue')
plt.show()

# End
