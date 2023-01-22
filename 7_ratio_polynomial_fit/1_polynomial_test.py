# This code is to fit the ratio files with six degree polynonial

from scipy.special import legendre
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt


def generate_legendre(n):
    return legendre(n)


def legendre_function(x, n, a):
    coeffs = legendre(n)
    func = None
    for deg, coeff in enumerate(np.flip(coeffs)):
        if func is None:
            func = x**deg * coeff
        else:
            func += x**deg * coeff
    return a*func


def polynomial(x, a0, a1, a2, a3, a4, a5, a6):
    function = legendre_function(x, 0, a0) + legendre_function(x, 1, a1)
    + legendre_function(x, 2, a2) + legendre_function(x, 3, a3)
    + legendre_function(x, 4, a4) + legendre_function(x, 5, a5)
    + legendre_function(x, 6, a6)
    return function


def ordinary_polynomial(x, a0, a1, a2, a3, a4, a5, a6):
    return a6 + a5*x + a4*x**2 + a3*x**3 + a2*x**4 + a1*x**5 + a0*x**6


pdfplots = PdfPages('test_fig_polynomial_fit.pdf')

ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')

for order, ratio in enumerate(ratio_file):
    pixels = np.arange(0, 2048, 1)
    params, sigma = curve_fit(ordinary_polynomial, pixels[600:1400],
                              ratio[600:1400])
    print(params)
    a0, a1, a2, a3, a4, a5, a6 = params
    fig = plt.figure(figsize=(16, 9))
    plt.plot(pixels, ratio)
    plt.plot(pixels[600:1400],
             ordinary_polynomial(pixels[600:1400], a0, a1, a2, a3, a4, a5, a6),
             '--')
    pdfplots.savefig()
pdfplots.close()
