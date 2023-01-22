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
        0: {1: (pixels > -1)},  # {1: (pixels < 930) & (pixels > -1),
        # 2: (pixels > 1750) & (pixels < 2048)},
        1: {1: (pixels < 850) & (pixels > 500),
            2: (pixels > 1080) & (pixels < 2048)},
        2: {1: (pixels < 2048) & (pixels > -1)},
        3: {1: (pixels < 1060) & (pixels > -1),
            2: (pixels > 1130) & (pixels < 2048)},
        4: {1: (pixels < 2048) & (pixels > -1)},
        5: {1: (pixels < 750) & (pixels > -1),
            2: (pixels > 1100) & (pixels < 2048)},
        6: {1: (pixels < 820) & (pixels > 0),
            2: (pixels > 1225) & (pixels < 2048)},
        7: {1: (pixels > 0) & (pixels < 600),
            2: (pixels > 760) & (pixels < 860),
            3: (pixels > 1700) & (pixels < 2048)},
        8: {1: (pixels < 1240) & (pixels > 0),
            2: (pixels > 1990) & (pixels < 2048)},
        9: {1: (pixels < 350) & (pixels > -1),
            2: (pixels > 600) & (pixels < 2048)}
        }
    return masks


def ordinary_polynomial(x, a0, a1, a2, a3, a4, a5, a6):
    return a6 + a5*x + a4*x**2 + a3*x**3 + a2*x**4 + a1*x**5 + a0*x**6


pdfplots = PdfPages('test_polyfit.pdf')

ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')

for order, ratio in enumerate(ratio_file):
    pixels_actual = np.arange(0, 2048, 1)
    mask_dict = mask_dictionary(pixels_actual)
    fig = plt.figure(figsize=(16, 9))
    pixels = pixels_actual/1024 - 1
    plt.plot(pixels_actual, ratio)
    # print(mask_dict)
    master_mask = np.sum(list(mask_dict[order].values()), axis=0,
                         dtype=np.bool)
    print(master_mask)
    params = np.polyfit(pixels[master_mask],
                        ratio[master_mask], 6)
    a0, a1, a2, a3, a4, a5, a6 = params
    plt.plot(pixels_actual[master_mask], ratio[master_mask])
    plt.plot(pixels_actual,
             ordinary_polynomial(pixels,
                                 a0, a1, a2, a3, a4, a5, a6), '--')
    plt.ylim(-1e14, max(ratio)+1000)
    pdfplots.savefig()
pdfplots.close()
