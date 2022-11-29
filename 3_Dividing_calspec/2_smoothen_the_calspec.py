
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import os


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def tanspec_data(q=0):
    file_path = '/home/varghese/Desktop/DP2_TANSPEC/\
Data/Output_spectra/hd37725/s0.5/tanspec_slopes'
    flux = fits.getdata(os.path.join(file_path,
                                     'op_hd37725_grating1_A.ms.wlc.fits')
                        )
    order = flux[q]
    wavelength = np.load(
        os.path.join(file_path,
                     'op_hd37725_grating1_A_combarc.OutputWavlFile{}.npy'
                     .format(q))
    )
    return order, wavelength


calspec_data = Table.read(
    '/home/varghese/Desktop/DP2_TANSPEC/Calspec_data/hd37725_mod_004.fits'
)

calspec_wl = np.array(calspec_data['WAVELENGTH'])
calspec_flux = np.array(calspec_data['FLUX'])
smoothed_flux = smooth(calspec_flux, 2048)
interpolation = interpolate.interp1d(calspec_wl, smoothed_flux)

pdfplots = PdfPages('Images/calspec_ratio.pdf')
for order in range(0, 10, 1):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle('Order {}'.format(order))
    gs = fig.add_gridspec(3, hspace=0)
    axs = gs.subplots(sharex=True)
    tanspec_flux, wavelength = tanspec_data(order)
    cal_orders = interpolation(wavelength)
    smoothed_tanspec = smooth(tanspec_flux, 2048)
    flux_ratio = cal_orders / smoothed_tanspec
    axs[0].plot(wavelength, cal_orders, label='calspec')
    axs[1].plot(wavelength, tanspec_flux, color='red', label='tanspec')
    axs[1].plot(wavelength, smoothed_tanspec, color='black', label='smoothend')
    axs[2].plot(wavelength, flux_ratio, color='green', label='ratio')
    np.save('smoothed_ratio_{}.npy'.format(order), flux_ratio)
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('Flux')
    pdfplots.savefig()
pdfplots.close()


# Code ends here
