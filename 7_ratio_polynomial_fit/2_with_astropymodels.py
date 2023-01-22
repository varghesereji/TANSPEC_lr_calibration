# This code is to fit the ratio with standard models in astropy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.modeling import models
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_generic_continuum
from astropy import units as u

pdfplots = PdfPages('test_astropy_fitting.pdf')
ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')*u.adu
pixels = np.arange(0, 2048, 1)*u.um
for order, ratio in enumerate(ratio_file):
    spectrum = Spectrum1D(flux=ratio, spectral_axis=pixels)

    # g_init = models.Polynomial1D(6, (0, 2048), (-1, 1))
    g_fit = fit_generic_continuum(spectrum)
    ratio_fit = g_fit(pixels)

    plt.figure(figsize=(16, 9))
    plt.plot(pixels, ratio)
    plt.plot(pixels, ratio_fit, '--')
    pdfplots.savefig()
pdfplots.close()
# End
