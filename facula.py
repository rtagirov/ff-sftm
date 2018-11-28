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

    nproc = 4

    C = 'C22'

    D = '26.5'

    B_sat = 484.0

    B_spot = 1000.0

    y_c = 0.0

    for i, arg in enumerate(args):

        if arg == '--np':

            nproc = int(args[i + 1])

        if arg == '--C':

            C = args[i + 1]

        if arg == '--D':

            D = args[i + 1]

        if arg == '--Bsat':

            B_sat = float(args[i + 1])

        if arg == '--Bspot':

            B_spot = float(args[i + 1])

        if arg == '--i':

            y_c = float(args[i + 1])

    return nproc, C, D, B_sat, B_spot, y_c

mag = './inp/mag/'

if not os.path.isdir(mag):

    auxsys.abort('The magnetograms directory is missing.')

if not os.listdir(mag):

    auxsys.abort('The magnetograms directory is empty.')

nproc, C, D, B_sat, B_spot, y_c = get_args(sys.argv[1:])

conv = np.pi / 180.0

norm = 90 * 90 * 4 / np.pi**2 * np.pi

x_c = 0.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

spot_mask = np.load('./out/npy/' + C + '_D' + D + '.npy').item()

times = np.array(list(spot_mask.keys()))

start = 170049

def scan_mag(date):

    B0 = np.loadtxt(mag + 'CalcMagnetogram.2000.' + str(date[0]))
    B1 = np.loadtxt(mag + 'CalcMagnetogram.2000.' + str(date[1]))

    t = times[np.where((times >= date[0]) & (times < date[1]))] - date[0]

    v = np.zeros(len(t))

    r = np.zeros((len(t), 11))

    B_tot = np.zeros(len(t))
    h_tot = np.zeros(len(t))

    h_cnt = np.zeros(len(t), dtype = 'int')

    for k in range(len(t)):

        time = date[0] + t[k]

        spot_x = np.concatenate((spot_mask[time]['xp'], spot_mask[time]['xn']))
        spot_y = np.concatenate((spot_mask[time]['yp'], spot_mask[time]['yn']))

        for i, j in itertools.product(range(180), range(360)):

            B = abs((B1[i, j] - B0[i, j]) * t[k] + B0[i, j])

            B_tot[k] += B

            n = 0

            if len(spot_x) != 0:

                n = len(np.where((spot_x >= j) & (spot_x < j + 1) & (spot_y >= i) & (spot_y < i + 1))[0])

            ff = 0.0

            if n != 0:

                helper = B - B_spot * n * 0.1 * 0.1

                h_tot[k] += B - helper

                if helper > 0 and helper <= B_sat:

#                    ff = (1 - n * 0.1 * 0.1) * helper / B_sat
                    ff = helper / B_sat

                if helper > B_sat:

                    ff = (1 - n * 0.1 * 0.1)

                if helper <= 0.0:

                    h_cnt[k] += 1

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

    return t + date[0], r / norm, v / norm, B_tot, h_tot, h_cnt

dates = [[i, i + 1] for i in range(math.floor(min(times)), math.ceil(max(times)) + 1)]

fpath = './out/' + C + '_' + D + '_' + str(int(B_sat)) + '_' + str(int(B_spot)) + '_' + str(int(90 - y_c))

f = open(fpath, 'w')

os.system('chmod 754 ' + fpath)

fmt = '%9.2f ' + '%10.6f ' * 13 + '%9.2f ' * 2 + '%i\n'

with Pool(processes = nproc) as p:

    maximum = len(dates)

    with tqdm(total = maximum, \
              ncols = auxfunc.term_width(), \
              desc = C + ', D = ' + D + \
              ', Bsat = ' + str(int(B_sat)) + \
              ', Bspot = ' + str(int(B_spot)) + \
              ', i = ' + str(int(90 - y_c)), \
              position = 0) as pbar:

        results = p.imap(scan_mag, dates)

        for result in results:

            t, r, v, B_tot, h_tot, h_cnt = result

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
                               v[k], \
                               B_tot[k], \
                               h_tot[k], \
                               h_cnt[k]))

            pbar.update()

    p.close()
    p.join()

f.close()

