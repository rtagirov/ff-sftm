import importlib
import itertools
import auxfunc
import auxsys
import sys
import glob
import math
import os

import numpy as np

from multiprocessing import Pool
from tqdm import tqdm

importlib.reload(auxfunc)
importlib.reload(auxsys)

def get_args(args):

    D = '26.5'

    B_sat = 484.0

    nproc = 4

    for i, arg in enumerate(args):

        if arg == '--Bsat':

            B_sat = float(args[i + 1])

        if arg == '--D':

            D = args[i + 1]

        if arg == '--np':

            nproc = int(args[i + 1])

    return D, B_sat, nproc

mag = './mag/'

if not os.path.isdir(mag):

    auxsys.abort('The directory with magnetograms is missing.')

D, B_sat, nproc = get_args(sys.argv[1:])

conv = np.pi / 180.0

norm = 90 * 90 * 4 / np.pi**2 * np.pi

x_c = 0.0
y_c = 0.0

B_spot = 1000.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

spot_mask = np.load('./out/spot_mask_D' + D + '.npy').item()

times = np.array(list(spot_mask.keys()))

#step = times[1] - times[0]

start = 170049

def scan_mag(date):

    B0 = np.loadtxt(mag + 'CalcMagnetogram.2000.' + str(date[0]))
    B1 = np.loadtxt(mag + 'CalcMagnetogram.2000.' + str(date[1]))

#    t = np.arange(0.0, 1.0, step)

    t = times[np.where((times >= date[0]) & (times < date[1]))] - date[0]

    v = np.zeros(len(t))

    r = np.zeros((len(t), 11))

    for k in range(len(t)):

        time = date[0] + t[k]

        spot_x = np.concatenate((spot_mask[time]['xp'], spot_mask[time]['xn']))
        spot_y = np.concatenate((spot_mask[time]['yp'], spot_mask[time]['yn']))

        for i, j in itertools.product(range(180), range(360)):

            B = abs((B1[i, j] - B0[i, j]) * t[k] + B0[i, j])

#            spot = spot_x[np.where((spot_x >= j) & (spot_x < j + 1) & (spot_y >= i) & (spot_y < i + 1))]
            n = len(np.where((spot_x >= j) & (spot_x < j + 1) & (spot_y >= i) & (spot_y < i + 1))[0])

            ff = 0.0

            if n != 0:

                helper = B - B_spot * n * 0.1 * 0.1

                if helper > 0:

                    ff = (1 - n * 0.1 * 0.1) * helper / B_sat

            elif n == 0 and B < B_sat:

                ff = B / B_sat

            elif n == 0 and B >= B_sat:

                ff = 1.0

            x_rot = (j + 13.28 * (time - start)) % 359

            x_pos = 180.0 - x_rot

            y_pos = 90.0 - i

            delta_lambda = abs(x_pos - x_c)

            distance = np.arccos(np.sin(y_c * conv) * np.sin(y_pos * conv) +
                                 np.cos(y_c * conv) * np.cos(y_pos * conv) * np.cos(delta_lambda * conv)) / conv

            vis = np.cos(distance * conv)

            l = np.where((vis > mu_low) & (vis <= mu_up))

            r[k, l] += ff * vis * np.cos(y_pos * conv)

            if distance <= 90.0:

                v[k] += ff * np.cos(distance * conv) * np.cos(y_pos * conv)

    return t + date[0], r / norm, v / norm

dates = [[i, i + 1] for i in range(math.floor(min(times)), math.ceil(max(times)) + 1)]

fpath = './out/' + D + '_' + str(int(B_sat))

f = open(fpath, 'w')

os.system('chmod 754 ' + fpath)

fmt = '%9.2f ' + '%10.6f ' * 12 + '%10.6f\n'

with Pool(processes = nproc) as p:

    maximum = len(dates)

    with tqdm(total = maximum, \
              ncols = auxfunc.term_width(), \
              desc = 'D = ' + D + ', Bsat = ' + str(B_sat), \
              position = 0) as pbar:

        results = p.imap(scan_mag, dates)

        for result in results:

            t, r, v = result

            for k in range(len(t)):

                f.write(fmt % (t[k], \
                               r[k, 0], \
                               r[k, 1], \
                               r[k, 2], \
                               r[k, 3], \
                               r[k, 4], \
                               r[k, 5], \
                               r[k, 6], \
                               r[k, 7], \
                               r[k, 8], \
                               r[k, 9], \
                               r[k, 10], \
                               sum(r[k, :]), \
                               v[k]))

            pbar.update()

    p.close()
    p.join()

f.close()

