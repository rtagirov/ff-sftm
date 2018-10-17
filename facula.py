import importlib
import itertools
import auxfunc
import sys
import glob
import re

import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

#from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Pool, freeze_support
#from multiprocessing import freeze_support

from tqdm import tqdm
from functools import partial
from itertools import repeat

importlib.reload(auxfunc)

conv = np.pi / 180.0

# set position of observer and Bsat
x_c = 0
y_c = 0
B_sat = 484.0
B_spot = 1000.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

start = 170049

def scan_mag(px, spot_x, spot_y, data, date):
#def scan_mag(pair):

    i = px[0]

    j = px[1]

#    print(i)

#    j = arg[0][1]

#    spot_x = arg[2]

#    spot_y = arg[3]

#    data = arg[4]
#    B = arg[4]

#    date = arg[5]

#    for i in range(180):

#        for j in range(360):

#    for j in range(360):

    x_min = j
    x_max = j + 1

    y_min = i
    y_max = i + 1

    spot = spot_x[np.where((spot_x >= x_min) & (spot_x < x_max) & (spot_y >= y_min) & (spot_y < y_max))]

    B = abs(data[i][j])

    if np.shape(spot) != (0, ):

        helper = B - B_spot * len(spot) * 0.1 * 0.1

        if helper > 0:

            ff[i, j] = (1 - len(spot) * 0.1 * 0.1) * helper / B_sat

    if np.shape(spot) == (0, ) and B < B_sat:

        ff[i, j] = B / B_sat

    if np.shape(spot) == (0, ) and B >= B_sat:

        ff[i, j] = 1.0

    x_rot = (j + 13.28 * (date - start)) % 359

    x_pos = 180.0 - x_rot

    y_pos = 90.0 - i

    delta_lambda = abs(x_pos - x_c)

    distance = np.arccos(np.sin(y_c * conv) * np.sin(y_pos * conv) + np.cos(y_c * conv) * \
                         np.cos(y_pos * conv) * np.cos(delta_lambda * conv)) / conv

    vis = np.cos(distance * conv)

    idx = np.where((vis > mu_low) & (vis <= mu_up))

#    lock.acquire()

    r[idx] += ff[i, j] * vis * np.cos(y_pos * conv)

#    lock.release()

    if distance <= 90.0:

#        lock.acquire()

        visibility.append(ff[i, j] * np.cos(distance * conv) * np.cos(y_pos * conv))

#        lock.acquire()

#def init():

#    global r

#    global visibility

#    r = np.zeros(11)

#    visibility = []

#    global lock

#    lock = l

nproc = 4

if len(sys.argv) == 2:

    nproc = int(sys.argv[1])

norm = 90 * 90 * 4 / np.pi**2 * np.pi

magnetograms = sorted(glob.glob('./mag/CalcMagnetogram.2000.*'))

spot_mask = np.load('spot_mask.npy').item()

f = open('ff_fac.out','w')

px = list(itertools.product(range(180), range(360)))

for _, mag in enumerate(tqdm(magnetograms, \
                             ncols = auxfunc.term_width(), \
                             desc = 'Masking faculae', \
                             position = 0)):

    name = re.findall('2000.(\d+)', mag)

    visibility = []

    data = np.loadtxt(mag)

    date = int(name[0])

    r = np.zeros(11)

    ff = np.zeros((180, 360))

    spot_x = np.concatenate((spot_mask[date]['xp'], spot_mask[date]['xn']))
    spot_y = np.concatenate((spot_mask[date]['yp'], spot_mask[date]['yn']))

#    i = range(180)

#    j = range(360)

#    args = [[i, j] for i in range(180) for j in range(360)]

#    for arg in args:

#        arg.append(spot_x)
#        arg.append(spot_y)
#        arg.append(data)
#        arg.append(date)

#    print(arg[0])

#    p.terminate()

#    sys.exit()

#    freeze_support()

#    l = mp.Lock()

#    p = Pool(processes = nproc, initializer = init, initargs = (r, visibility,))
#    p = Pool(processes = nproc, initializer = init)
    p = Pool(processes = nproc)

#    with Pool(processes = nproc, initializer = init, initargs = (l,)) as p:

#        i = range(180)

#        j = range(360)

#        p.imap(scan_mag, ([i, spot_x, spot_y, data, date] for i in range(180)))
#        p.map(scan_mag, list(itertools.product(i, j)))
#    p.map(partial(scan_mag, spot_x = spot_x, spot_y = spot_y, data = data, date = date), args)
    p.starmap(scan_mag, zip(px, repeat(spot_x), repeat(spot_y), repeat(data), repeat(date)))
#        p.imap(f, range(180))
#        p.imap(f, ([pair, spot_x, spot_y, data, date] for i, j in itertools.product(range(180), range(360))))
#        p.imap(f, ([i, j, data[i][j], date] for i, j in itertools.product(range(180), range(360))))

    p.close()
    p.join()

    r /= norm

    f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' \
            %(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], \
              date, sum(r), sum(visibility) / norm))

f.close()

