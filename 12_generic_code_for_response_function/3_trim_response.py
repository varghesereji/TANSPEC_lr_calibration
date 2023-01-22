import numpy as np

pixels = np.arange(0, 2048, 1)
mask0 = (pixels > 0)
mask1 = (pixels > 499) & (pixels < 2000)
mask2 = (pixels > 369) & (pixels < 1500)
mask3 = (pixels > 279) & (pixels < 1700)
mask4 = (pixels > 289) & (pixels < 1750)
mask5 = (pixels > -1)
mask6 = (pixels > -1)
mask7 = (pixels > -1)
mask8 = (pixels > -1)
mask9 = (pixels > 160) & (pixels < 1712)

Full_mask = mask0
Full_mask = np.vstack((Full_mask, mask1))
Full_mask = np.vstack((Full_mask, mask2))
Full_mask = np.vstack((Full_mask, mask3))
Full_mask = np.vstack((Full_mask, mask4))
Full_mask = np.vstack((Full_mask, mask5))
Full_mask = np.vstack((Full_mask, mask6))
Full_mask = np.vstack((Full_mask, mask7))
Full_mask = np.vstack((Full_mask, mask8))
Full_mask = np.vstack((Full_mask, mask9))
np.save('/home/varghese/Desktop/DP2_TANSPEC/Data/Instrument_response/smoothed/Trim_spectrum.npy', Full_mask)

# end
