import importlib
import itertools
import auxfunc
import auxsys
import sys
import glob
import re
import os

import numpy as np

from multiprocessing import Pool
from tqdm import tqdm

importlib.reload(auxfunc)
importlib.reload(auxsys)

nproc = 4

if len(sys.argv) == 2:

    nproc = int(sys.argv[1])

conv = np.pi / 180.0

norm = 90 * 90 * 4 / np.pi**2 * np.pi

x_c = 0.0
y_c = 0.0

B_sat = 484.0
B_spot = 1000.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

#spot_mask = np.load('spot_mask_nproc_' + str(nproc) + '.npy').item()
spot_mask = np.load('spot_mask.npy').item()

start = 170049

def scan_mag(mag):

    name = re.findall('2000.(\d+)', mag)

    data = np.loadtxt(mag)

    date = int(name[0])

    visibility = 0.0

    r = np.zeros(11)

    spot_x = np.concatenate((spot_mask[date]['xp'], spot_mask[date]['xn']))
    spot_y = np.concatenate((spot_mask[date]['yp'], spot_mask[date]['yn']))

    for i, j in itertools.product(range(180), range(360)):

        spot = spot_x[np.where((spot_x >= j) & (spot_x < j + 1) & (spot_y >= i) & (spot_y < i + 1))]

        B = abs(data[i][j])

        ff = 0.0

        if np.shape(spot) != (0, ):

            helper = B - B_spot * len(spot) * 0.1 * 0.1

            if helper > 0:

                ff = (1 - len(spot) * 0.1 * 0.1) * helper / B_sat

        if np.shape(spot) == (0, ) and B < B_sat:

            ff = B / B_sat

        if np.shape(spot) == (0, ) and B >= B_sat:

            ff = 1.0

        x_rot = (j + 13.28 * (date - start)) % 359

        x_pos = 180.0 - x_rot

        y_pos = 90.0 - i

        delta_lambda = abs(x_pos - x_c)

        distance = np.arccos(np.sin(y_c * conv) * np.sin(y_pos * conv) + 
                             np.cos(y_c * conv) * np.cos(y_pos * conv) * np.cos(delta_lambda * conv)) / conv

        vis = np.cos(distance * conv)

        idx = np.where((vis > mu_low) & (vis <= mu_up))

        r[idx] += ff * vis * np.cos(y_pos * conv)

        if distance <= 90.0:

            visibility += ff * np.cos(distance * conv) * np.cos(y_pos * conv)
                    
    r /= norm

    visibility /= norm

    return date, r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], sum(r), visibility

mag_dir = './mag_21/'

if not os.path.isdir(mag_dir):

    auxsys.abort('The directory with magnetograms is missing. Abort.')

magnetograms = sorted(glob.glob(mag_dir + 'CalcMagnetogram.2000.*'))

#f = open('ff_fac_nproc_' + str(nproc) + '.out', 'w')
f = open('ff_fac.out', 'w')

fmt = '%i ' + '%10.6f ' * 12 + '%10.6f\n'

with Pool(processes = nproc) as p:

    maximum = len(magnetograms)

    with tqdm(total = maximum, \
              ncols = auxfunc.term_width(), \
              desc = 'Masking faculae, nproc = ' + str(nproc), \
              position = 0) as pbar:

#        results = p.imap(scan_mag, magnetograms, chunksize = 4)
        results = p.imap(scan_mag, magnetograms)

        for _, result in enumerate(results):

            f.write(fmt % result)

            pbar.update()

    p.close()
    p.join()

f.close()

