import numpy as np
import matplotlib.pyplot as plt


def gaussian(mu, sigma):
    print(mu, sigma)
    lammin = mu - 5*sigma
    lammax = mu + 5*sigma
    lambdas = np.linspace(lammin, lammax, 100)
    normalization = 1 / (sigma * np.sqrt(2*np.pi))
    arguement = (lambdas - mu)**2 / (2 * sigma**2)
    func = normalization * np.exp(-arguement)
    # plt.figure()
    # plt.plot(lambdas, func)
    # plt.show()
    return func, lambdas


def gaussian_sets(wl_values):
    R = 230
    gaussian_sets = None
    wls = None
    for wl in wl_values:
        print(wl)
        dlambda = wl / R
        fn, lambdas = gaussian(wl, dlambda)
        if gaussian_sets is None:
            gaussian_sets = fn
            wls = lambdas
        else:
            gaussian_sets = np.vstack((gaussian_sets, fn))
            wls = np.vstack((wls, lambdas))
    return gaussian_sets, wls


wavelengths = np.load('wavelength_calibrated.npy')
m_wl = wavelengths > 0
wavelengths = wavelengths[m_wl]

gauss_sets, wls = gaussian_sets(wavelengths)

plt.figure(figsize=(16, 9))
_ = plt.plot(wls.T, gauss_sets.T)
plt.savefig('gaussians.png')

# End
