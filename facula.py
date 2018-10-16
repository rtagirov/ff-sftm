import importlib
import itertools
import auxfunc
import mask
import sys
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool
from tqdm import tqdm

importlib.reload(auxfunc)
importlib.reload(mask)

conv = np.pi / 180.0

# defining all the constants
norm = 90 * 90 * 4 / np.pi**2 * np.pi

# set position of observer and Bsat
x_c = 0
y_c = 0
B_sat = 484.0
B_spot = 1000.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

#magnetograms = sorted(glob.glob('./mag_test/CalcMagnetogram.2000.*'))
magnetograms = sorted(glob.glob('./mag/CalcMagnetogram.2000.*'))

#spot_mask = np.load('spot_mask_4proc.npy').item()
spot_mask = np.load('spot_mask.npy').item()

start = 170049

f = open('ff_fac.out','w')

for k, mag in enumerate(tqdm(magnetograms, ncols = auxfunc.term_width(), desc = 'Masking faculae', position = 0)):

    name = re.findall('2000.(\d+)', mag)

    visibility = []

    data = np.loadtxt(mag)

    date = int(name[0])

    r = np.zeros(11)

    ff = np.zeros((180, 360))

    spot_x = np.concatenate((spot_mask[date]['xp'], spot_mask[date]['xn']))
    spot_y = np.concatenate((spot_mask[date]['yp'], spot_mask[date]['yn']))

    for i in range(180):

        for j in range(360):

            x_min = j
            x_max = j + 1

            y_min = i
            y_max = i + 1

            spot = spot_x[np.where((spot_x >= x_min) & (spot_x < x_max) & (spot_y >= y_min) & (spot_y < y_max))]

            B = abs(data[i][j])

            if np.shape(spot) != (0, ):

                helper = B - B_spot * len(spot) * 0.1 * 0.1

                if helper <= 0:

                    ff[i, j] = 0.0

                if helper > 0:

                    ff[i, j] = (1 - len(spot) * 0.1 * 0.1) * helper / B_sat

            if np.shape(spot) == (0, ) and B < B_sat:

                ff[i, j] = abs(data[i][j]) / B_sat

            if np.shape(spot) == (0, ) and B >= B_sat:

                ff[i, j] = 1.0

            x_rot = (j + 13.28 * (int(date) - start)) % 359

            x_pos = 180.0 - x_rot

            y_pos = 90.0 - i

            delta_lambda = abs(x_pos - x_c)

            distance = np.arccos(np.sin(y_c * conv) * np.sin(y_pos * conv) + np.cos(y_c * conv) * \
                       np.cos(y_pos * conv) * np.cos(delta_lambda * conv)) / conv

            vis = np.cos(distance * conv)

            idx = np.where((vis > mu_low) & (vis <= mu_up))

            r[idx] += ff[i, j] * vis * np.cos(y_pos * conv)

            if distance <= 90.0:

                visibility.append(ff[i, j] * np.cos(distance * conv) * np.cos(y_pos * conv))

    r /= norm

    f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' \
            %(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], date, sum(r), sum(visibility) / norm))

f.close()

