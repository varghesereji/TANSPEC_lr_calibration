# This code is to Filter out the absorption lines in the spectra

from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy import interpolate
from scipy import ndimage
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os


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
smoothend_flux = ndimage.rank_filter(calspec_flux, rank=-10, size=2000)
interpolation_smoothed = interpolate.interp1d(calspec_wl, smoothend_flux)
interpolation_nsm = interpolate.interp1d(calspec_wl, calspec_flux)

pdfplots = PdfPages('test_images.pdf')
for order in range(0, 10, 1):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle('Order {}'.format(order))
    gs = fig.add_gridspec(2, hspace=0)
    axs = gs.subplots(sharex=True)
    tanspec_flux, wavelength = tanspec_data(order)
    cal_orders_sm = interpolation_smoothed(wavelength)
    cal_orders_nsm = interpolation_nsm(wavelength)
    axs[0].plot(wavelength, cal_orders_nsm)
    axs[1].plot(wavelength, cal_orders_sm)
    plt.xlabel('$\lambda \AA$')
    plt.ylabel('Flux')
    pdfplots.savefig()
pdfplots.close()


# End
