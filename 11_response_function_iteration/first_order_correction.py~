import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

order = 9
ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')
ratio = ratio_file[order]
mask_full = np.load('Mask_for_Instrument_response.npy')
pixels_actual = np.arange(0, 2048, 1)
pixels = pixels_actual/1024 - 1


masked_reg = mask_full[order]
mask = ~masked_reg

interpolation_mask = (pixels_actual > 290) & (pixels_actual < 350)

interpolation_region = np.linspace(min(pixels_actual[interpolation_mask]),
                                   max(pixels_actual[interpolation_mask]),
                                   100)


plt.figure(figsize=(16, 9))
diff_mask = np.diff(mask)
positions = list(np.where(diff_mask)[0])
mask_positions = list(np.where(~mask)[0])
for x_pos in mask_positions:
    plt.axvline(x=x_pos, color='yellow', linestyle='-', alpha=0.1)
for mark in positions:
    plt.text(mark, 0, str(mark+1), fontsize=15)
plt.plot(pixels_actual, ratio, '-',)
polyfit = np.polyfit(pixels[masked_reg],
                     ratio[masked_reg], 5)
fitted = np.poly1d(polyfit)
plt.plot(pixels_actual[masked_reg], ratio[masked_reg], color='red')
plt.plot(pixels_actual, fitted(pixels), '--', color='green')
plt.xlabel('pixels')
plt.ylim(-1e14, max(ratio)+1e15)
plt.xlim(0, 2048)
plt.show()
print('done', positions)

'''
Page 2 masked_red = (pixels_actual < 500) \
        | (pixels_acual > 850) & (pixels_actual < 1080)

Page 3 masked_reg = (pixels_actual < 370) \
        | (pixels_actual > 1500) \
        | (pixels_actual > 1070) & (pixels_actual < 1110)

Page 4    masked_reg = (pixels_actual < 300) \
        | (pixels_actual > 1045) & (pixels_actual < 1130) \
        | (pixels_actual > 1700)

Page 5     masked_reg = (pixels_actual < 290) \
        | (pixels_actual > 1740) \
        | (pixels_actual > 800) & (pixels_actual < 870)

Page 6 
Page 8    masked_reg = (pixels_actual > 600) & (pixels_actual < 760) \
        | (pixels_actual > 860) & (pixels_actual < 1701)

Page 9 masked_reg = (pixels_actual > 1240) & (pixels_actual < 1990) \
        | (pixels_actual > 920) & (pixels_actual < 1181) \
        | (pixels_actual > 255) & (pixels_actual < 310) \
        | (pixels_actual > 632) & (pixels_actual < 760)

Page 10.     masked_reg = (pixels_actual > 355) & (pixels_actual < 700) \
        | (pixels_actual > -1) & (pixels_actual < 137) \
        | (pixels_actual > 1490) & (pixels_actual < 1530) \
        | (pixels_actual > 1459) & (pixels_actual < 1486) \
        | (pixels_actual > 1388) & (pixels_actual < 1450) \
        | (pixels_actual > 1551) & (pixels_actual < 1635) \
        | (pixels_actual > 1296) & (pixels_actual < 1377) \
        | (pixels_actual > 1213) & (pixels_actual < 1285) \
        | (pixels_actual > 865) & (pixels_actual < 914) \
        | (pixels_actual > 1164) & (pixels_actual < 1200) \
        | (pixels_actual > 1103) & (pixels_actual < 1160) \
        | (pixels_actual > 1712)

'''

# End

