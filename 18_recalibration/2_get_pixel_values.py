import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import least_squares
import csv


def polynomial(p, a0, a1, a2, a3):
    return a0 * p**3 + a1 * p**2 + a2 * p + a3


wls = []
pixels = []

argon = []
argon_pixels = []
with open('wl_matching/argon.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        argon.append(float(lines[0]))
        argon_pixels.append(float(lines[1]))

neon = []
neon_pixels = []
with open('wl_matching/neon.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        neon.append(float(lines[0]))
        neon_pixels.append(float(lines[1]))

sky = []
sky_pixels = []
with open('wl_matching/sky.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        sky.append(float(lines[0]))
        sky_pixels.append(float(lines[1]))

telluric = []
telluric_pixels = []
with open('wl_matching/telluric.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        wls.append(float(lines[0]))
        pixels.append(float(lines[1]))
        telluric.append(float(lines[0]))
        telluric_pixels.append(float(lines[1]))

recalib_star_lr = []
recalib_star_xd = []
with open('recal_star.csv', mode='r') as file:
    csvFile = csv.reader(file)
    for lines in csvFile:
        recalib_star_lr.append(float(lines[0]))
        recalib_star_xd.append(float(lines[1]))
       # wls.append(float(lines[1]))


full_pixels = np.arange(700, 1600, 1)
fitting = np.polyfit(wls, pixels, 3)
curve = np.poly1d(fitting)
recal_pixl = []
for ind, lr_wl in enumerate(recalib_star_lr):
    print(lr_wl)
    print(lr_wl, curve(lr_wl))
    recal_pixl.append(curve(lr_wl))
    pixels.append(curve(lr_wl))
    wls.append(recalib_star_xd[ind])

print(len(wls), len(pixels))
pixels = np.array(pixels)
wl = np.array(wls)

fitting = np.polyfit(pixels, wl, 3)
curve = np.poly1d(fitting)
fitted_wls = curve(full_pixels)


mask = (np.gradient(fitted_wls) < 0)
neg_der_points = full_pixels[mask]
der_value = np.ones(np.sum(mask))*1e70
# print(neg_der_points.size, der_value.size)
# print(wl)

all_pixels = np.arange(0, 2048, 1)
mask = (all_pixels > 699) & (all_pixels < 1600)
wl_array = np.zeros(2048)
wl_array[mask] = fitted_wls
plt.figure(figsize=(16, 9))
plt.plot(full_pixels, fitted_wls)
plt.plot(argon_pixels, argon, '<', label='argon')
plt.plot(neon_pixels, neon, 'o', label='neon')
plt.plot(sky_pixels, sky, 's', label='sky')
plt.plot(recal_pixl, recalib_star_lr, '>', color='black',
         label='recalib points lr')
plt.plot(recal_pixl, recalib_star_xd, '>', color='blue',
         label='recalib points xd')
plt.plot(telluric_pixels, telluric, '*', label='telluric')
plt.xlabel('pixels')
plt.ylabel('wavelength $\AA$')
plt.legend()
plt.savefig('recalibrated_fitting.png')
np.save('recalibrated_wavelength.npy', wl_array)
# plt.show()
