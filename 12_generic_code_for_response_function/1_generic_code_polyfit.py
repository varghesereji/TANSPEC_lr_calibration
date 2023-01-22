from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


star = 'hd37725'
slit = 's0.5_by_s4.0'
# slit = 's4.0'
mask_full = np.load('Mask_for_Instrument_response.npy')
ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_{}.npy'.format(slit))

if slit == 's0.5_by_s4.0':
    slit = 's0.5'

pixels_actual = np.arange(0, 2048, 1)
pixels = pixels_actual/1024 - 1
pdfplots = PdfPages('response_fn_{}.pdf'.format(slit))
degrees = [0, 8, 8, 8, 8, 7, 8, 8, 6, 5]

response = None

for order, ratio in enumerate(ratio_file):
    masked_reg = mask_full[order]
    if order < 2:
        derivative = np.gradient(ratio)
        smoothed_derivative = smooth(derivative, 45)
        integral = integrate.cumtrapz(derivative, np.arange(0, 2048, 1),
                                      initial=0)
        correction_factor = np.mean(ratio - integral)
        fitted = integrate.cumtrapz(
            smoothed_derivative, np.arange(0, 2048, 1),
            initial=0)
    else:
        masked_reg = mask_full[order]
        polyfit = np.polyfit(pixels[masked_reg],
                             ratio[masked_reg],
                             degrees[order])
        fitted = np.poly1d(polyfit)
        fitted = fitted(pixels)
    if response is None:
        response = fitted
    else:
        response = np.vstack((response, fitted))
    plt.figure(figsize=(16, 9))
    plt.plot(pixels_actual, ratio, '-',)
    plt.plot(pixels_actual[masked_reg], ratio[masked_reg], color='red')
    plt.plot(pixels_actual, fitted, '--', color='green')
    print(order, 'done')
    plt.xlabel('pixels')
    plt.ylim(-1e14, max(ratio)*1.2)
    plt.xlim(0, 2048)
    plt.title('Order {0}, page {1}'.format(12 - order, order+1))
    pdfplots.savefig()
pdfplots.close()

np.save('Response_function_{}.npy'.format(slit), response)


# End
