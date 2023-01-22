import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from astropy.io import fits
from scipy.optimize import curve_fit
import os


project = '/home/varghese/Desktop/DP2_TANSPEC'
science_files = 'Data/Output_spectra/hd37725/s4.0/tanspec_slopes'
path = os.path.join(project, science_files)


def bb_spectra(lam, T, scale):
    a = (np.exp(const.h * const.c / (lam * const.k * T)) - 1)
    b = 2 * const.h * (const.c)**2 / (lam**5)
    bb_eq = scale * b * (1/a)
    return bb_eq


science_file = os.path.join(
    path, 'op_hd37725_s4.0_grating1_A.ms.wlc.fits')
bg1 = fits.getdata(science_file, ext=1)[9, 1024:]
wl_sc = fits.getdata(science_file, ext=3)[9, 1024:]

fmax = np.max(bg1)

mask = bg1 == fmax
wl_max = wl_sc[mask]
print(wl_max)

# wl_si = wl_sc * 10**(-10)  

# bb_spec = bb_spectra(wl_si, 1160, 1)

# params, err = curve_fit(bb_spectra, wl_si, bg1, p0=np.array([1160, 10**(-6)]))
# print(params)
# T, scale = params
# recreate_bb = bb_spectra(wl_si, T, scale)
# plt.figure(figsize=(16, 9))
# plt.plot(wl_sc, bg1)
# plt.show()

# End
