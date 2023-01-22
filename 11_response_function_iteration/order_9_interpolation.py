import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


def scale_pixels(pixels):
    return pixels/1024 - 1


def reverse_scaling(sc_pixels):
    return (sc_pixels+1) * 1024


def add_more_points(ratio, actual_pixels,
                    first_pixel, last_pixel):
    added_pixels = np.arange(first_pixel, last_pixel, 0.001)
    ratio_interp = interpolate.interp1d(
        scale_pixels(actual_pixels[first_pixel-10:last_pixel+10]),
        ratio[first_pixel-10:last_pixel+10])
    ratio_more_points = ratio_interp(scale_pixels(added_pixels))
    return ratio_more_points, added_pixels


order = 9
ratio_file = np.load('ratio_files/flux_cal_ratio_hd37725_s4.0.npy')
ratio = ratio_file[order]
mask_full = np.load('Mask_for_Instrument_response.npy')
pixels_actual = np.arange(0, 2048, 1)
scaled_pixels = scale_pixels(pixels_actual)
mask = mask_full[order]

masked_reg = (pixels_actual > 355) & (pixels_actual < 700) \
        | (pixels_actual > -1) & (pixels_actual < 150) \
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
mask = ~masked_reg

median_ratio = np.median(ratio)
ratio = ratio / median_ratio


added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 340, 356)
ratio_to_fit = np.concatenate((ratio[mask], added_ratio))
pixels_to_fit = np.concatenate((scaled_pixels[mask],
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 250, 260)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 200, 220)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 700, 720)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 150, 155)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 808, 812)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 740, 750)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1380, 1385)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1160, 1164)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1206, 1210)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1290, 1294)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 915, 925)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1085, 1095)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1455, 1458)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1533, 1558)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))

added_ratio, added_pixels = add_more_points(ratio, pixels_actual, 1204, 1206)
ratio_to_fit = np.concatenate((ratio_to_fit, added_ratio))
pixels_to_fit = np.concatenate((pixels_to_fit,
                                scale_pixels(added_pixels)))


back_scaled_pixels = reverse_scaling(pixels_to_fit)
plt.figure(figsize=(16, 9))
diff_mask = np.diff(mask)
positions = list(np.where(diff_mask)[0])
mask_positions = list(np.where(~mask)[0])
for x_pos in mask_positions:
    plt.axvline(x=x_pos, color='yellow', linestyle='-', alpha=0.1)
for mark in positions:
    plt.text(mark, 0, str(mark+1), fontsize=15)
plt.plot(pixels_actual, ratio, '-')
polyfit = np.polyfit(pixels_to_fit,
                     ratio_to_fit, 8)
fitted = np.poly1d(polyfit)
plt.plot(pixels_actual[mask], ratio[mask], color='red')
plt.plot(back_scaled_pixels, ratio_to_fit, '.', color='black')
plt.plot(pixels_actual, fitted(scaled_pixels), '--', color='green')
plt.xlabel('pixels')
plt.ylim(-0.1, 2)
plt.xlim(0, 2048)
plt.grid()
plt.show()
np.save('instrum_resp/Instrument_resp_9.npy',
        median_ratio * fitted(scaled_pixels))

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

