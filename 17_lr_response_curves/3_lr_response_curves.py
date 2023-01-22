from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import os


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
    interpolated = interp1d(calspec_wl, calspec_flux)
    interpolated_flux = interpolated(wavelength)
    return interpolated_flux


project = '/home/varghese/Desktop/DP2_TANSPEC'

# lr mode spectra
data_path = 'lr_background_reduction/spectrum_extraction/extracted_spectra'
file_path = os.path.join(project, data_path, 'Ext_hd37725_lr_s0.5hg.fits')
fitted_wavelength = np.load('wavelength_calibrated.npy')

star = fits.getdata(file_path, ext=0)[0]
sky = fits.getdata(file_path, ext=1)[0]

star_sky_sbtr = star - sky
mask = fitted_wavelength > 0

lr_wavelength = np.load('wavelength_calibrated.npy')

fig, axs = plt.subplots(5, 2)

# pdfplots = PdfPages('xd_to_lr_connection.pdf')
fig, ax = plt.subplots(3, figsize=(16, 9))
for order in range(0, 10, 1):

    data = np.load(
        'xd_flux_wavelength/flux_wavelength_{}.npy'.format(order)
    )
    xd_interp = interp1d(data[1], data[0])
    wl_mask = (lr_wavelength > data[1, 0]) & (lr_wavelength < data[1, -1])
    wl_region = fitted_wavelength[wl_mask]
    xd_updated = xd_interp(wl_region)
    lr_region = star_sky_sbtr[wl_mask]
    ax[0].plot(wl_region, lr_region)
    ax[1].plot(data[1], data[0])
    ax[1].plot(wl_region, xd_updated)
    response_curves = lr_region / xd_updated
    ax[1].plot(wl_region, lr_region / response_curves)
    ax[2].plot(wl_region, response_curves)
#    pdfplots.savefig()
    np.save('lr_response_wl_{}.npy'.format(order),
            np.vstack((response_curves, wl_region)))
plt.savefig('lr_response_curves.png')
#pdfplots.close()

# normalaizing
resp_together = None
for order in range(10):
    resp = np.load('lr_response_wl_{}.npy'.format(order))[0]
    print(resp)
    if resp_together is None:
        resp_together = resp
    else:
        resp_together = np.concatenate((resp_together, resp))
median = np.median(resp_together)
for order in range(10):
    resp = np.load('lr_response_wl_{}.npy'.format(order))
    resp_curve = resp[0]
    resp[0] = resp_curve / median
    np.save('lr_response_wl_{}.npy'.format(order), resp)

# End
