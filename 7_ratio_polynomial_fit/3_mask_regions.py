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


def mask_dictionary(pixels):
    masks = {
        0: {1: (pixels < 930) & (pixels > 800),
            2: (pixels > 1750) & (pixels < 1800)},
        1: {1: (pixels < 820) & (pixels > 580),
            2: (pixels > 1080) & (pixels < 1400)},
        2: {1: (pixels < 1400) & (pixels > 400)},
        3: {1: (pixels < 1050) & (pixels > 490),
            2: (pixels > 1135) & (pixels < 1540)},
        4: {1: (pixels < 1600) & (pixels > 320)},
        5: {1: (pixels < 750) & (pixels > 350),
            2: (pixels > 1100) & (pixels < 1700)},
        6: {1: (pixels < 820) & (pixels > 200),
            2: (pixels > 1225) & (pixels < 1730)},
        7: {1: (pixels > 0) & (pixels < 600),
            2: (pixels > 760) & (pixels < 860),
            3: (pixels > 1700) & (pixels < 2000)},
        8: {1: (pixels < 1240) & (pixels > 0),
            2: (pixels > 1990) & (pixels < 2040)},
        9: {1: (pixels < 250) & (pixels > 50),
            2: (pixels > 600) & (pixels < 1700)}
        }
    return masks


def ordinary_polynomial(x, a0, a1, a2, a3, a4, a5, a6):
    return a6 + a5*x + a4*x**2 + a3*x**3 + a2*x**4 + a1*x**5 + a0*x**6


pdfplots = PdfPages('test_mask_fit.pdf')

ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')

for order, ratio in enumerate(ratio_file):
    pixels = np.arange(0, 2048, 1)
    mask_dict = mask_dictionary(pixels)
    fig = plt.figure(figsize=(16, 9))
    pixels = pixels/1024 - 1
    plt.plot(pixels, ratio)
    # print(mask_dict)
    for segm in mask_dict[order].values():
        params, sigma = curve_fit(ordinary_polynomial, pixels[segm],
                                  ratio[segm])
        a0, a1, a2, a3, a4, a5, a6 = params

        plt.plot(pixels[segm], ratio[segm])
        plt.plot(pixels[segm],
                 ordinary_polynomial(pixels[segm], a0, a1, a2, a3, a4, a5, a6),
                 '--')
    pdfplots.savefig()
pdfplots.close()
