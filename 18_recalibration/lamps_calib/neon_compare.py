from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
from astropy.io import fits
import os

project = '/home/varghese/Desktop/DP2_TANSPEC'

# lr mode spectra
data_path = 'lr_background_reduction/spectrum_extraction/extracted_spectra'
file_path = os.path.join(project, data_path, 'Ext_neon_lr_s0.5hg.fits')
fitted_wavelength = np.load('../wavelength_calibrated.npy')

star = fits.getdata(file_path, ext=0)[0]
sky = fits.getdata(file_path, ext=1)[0]

star_sky_sbtr = star - sky
mask = fitted_wavelength > 0

lr_wavelength = np.load('../wavelength_calibrated.npy')

fig, axs = plt.subplots(5, 2)


xd_response_curves = np.load(os.path.join(project,
                                          'Data/Instrument_response/smoothed',
                                          'Normalized_response_curve.npy'))
xd_extracted_spectra = fits.getdata(
    os.path.join(project,
                 'Data/Output_spectra/hd37725/s0.5_by_s4.0/tanspec_slopes',
                 'op_hd37725_s0.5_by_s4.0_grating1_A_arc2.fits'), ext=0)

xd_extracted_wls = fits.getdata(
    os.path.join(project,
                 'Data/Output_spectra/hd37725/s0.5_by_s4.0/tanspec_slopes',
                 'op_hd37725_s0.5_by_s4.0_grating1_A.ms.wlc.final.fits'), ext=1)


def gaussian(mu, sigma, lmin, lmax, N):
    print(mu, sigma, N)
    lambdas = np.linspace(lmin, lmax, N)
#     normalization = 1 / (sigma * np.sqrt(2*np.pi))
    arguement = (lambdas - mu)**2 / (2 * sigma**2)
    print(arguement)
    print((lambdas-mu)**2/sigma**2)
    func = np.exp(-arguement)
    print(func)
    return func


resolutions = [110, 105, 105, 100, 100, 110, 120, 140, 160, 300]

#pdfplots = PdfPages('lr_calib_xd.pdf')
fig, ax = plt.subplots(figsize=(16, 4))
# ax[0].plot(lr_wavelength[mask], star_sky_sbtr[mask], color='blue', label='LR mode spectra')
# ax[0].legend()
for order in range(0, 10, 1):

    xd_masked_wl = np.load(
        '../xd_flux_wavelength/flux_wavelength_{}.npy'.format(order)
    )[1]
    wl_mask = (lr_wavelength > xd_masked_wl[0]) & (lr_wavelength < xd_masked_wl[-1])
    xd_resp = xd_response_curves[order]
    xd_wls = xd_extracted_wls[order]
    xd_flux = xd_extracted_spectra[order]

    lr_array = np.load('../lr_response_wl_{}.npy'.format(order))
    lr_resp = lr_array[0]
    lr_wl = lr_array[1]

    interp_xd_flux = interp1d(xd_wls, xd_flux, fill_value='extrapolate')
    interp_xd_resp = interp1d(xd_wls, xd_resp, fill_value='extrapolate')

    updated_xd_flux = interp_xd_flux(lr_wl)
    updated_xd_resp = interp_xd_resp(lr_wl)

    response_ratio = lr_resp / updated_xd_resp
    multiplied_xd = updated_xd_flux * response_ratio

    mu = np.mean(lr_wl)
    sigma = mu / (resolutions[order] * 2.355)
    lmin = min(lr_wl)
    lmax = max(lr_wl)
    print(order, lr_wl.size)
    gaussian_func = gaussian(mu, sigma, lmin, lmax, lr_wl.size)

    convolved_xd_spectra = convolve(multiplied_xd, gaussian_func, mode='same')
    if order == 0:
        ax.plot(lr_wl, convolved_xd_spectra, color='red', label='XD flux multiplied by ratio of XD and LR response functions and Convolved with gaussian')
    else:
        ax.plot(lr_wl, convolved_xd_spectra, color='red')
   # ax[0].set_ylim(-10, 5000)
ax.set_ylim(-10, 1000)
ax.legend()
ax.set_xlabel('Wavelength $\AA$')
ax.set_ylabel('Convolved XD flux', color='red')
ax.tick_params(axis='y', labelcolor='red')

ax2 = ax.twinx()
ax2.plot(lr_wavelength[mask], star[mask], color='blue', label='LR mode spectra')
ax2.set_ylabel('LR flux', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
plt.show()
# plt.savefig('Final_XD_LR_comparison.png')
#     pdfplots.savefig()
# pdfplots.close()

# End
