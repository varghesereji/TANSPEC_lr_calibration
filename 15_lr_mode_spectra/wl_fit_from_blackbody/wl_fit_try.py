import numpy as np
import matplotlib.pyplot as plt

full_pixels = np.arange(0, 2048, 1)

fitting = np.polyfit(pixels, wl, 3)
curve = np.poly1d(fitting)
fitted_wls = curve(full_pixels)
print(fitted_wls)

plt.figure()
plt.plot(full_pixels, fitted_wls)
plt.plot(pixels, wl, 'o')
plt.show()
