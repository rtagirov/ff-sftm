import importlib
import itertools
import auxfunc
import mask
import sys
import glob
import re

import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from multiprocessing import Manager

from multiprocessing.pool import ThreadPool as Pool

importlib.reload(auxfunc)
importlib.reload(mask)

factor = 10

xvalues = np.around(np.linspace(0.0, 359.0, num = 360 * factor), decimals = 1)
yvalues = np.around(np.linspace(0.0, 180.0, num = 181 * factor), decimals = 1)

x, y = np.meshgrid(xvalues, yvalues)

conv = np.pi / 180.0

data = np.loadtxt('spot_evol.inp')

long_pos = data[:, 2]
long_neg = data[:, 4]

lat_pos = data[:, 1]
lat_neg = data[:, 3]

grid_max = 359

def f(i):

    date = int(data[i, 0])

    r = np.sqrt(data[i, 5])

    phi_pos = 90.0 - lat_pos[i]
    phi_neg = 90.0 - lat_neg[i]

    x_pos_min = long_pos[i] - r / 2.0 * 1.0 / np.cos(phi_pos * conv)
    x_pos_max = long_pos[i] + r / 2.0 * 1.0 / np.cos(phi_pos * conv)

    y_pos_min = lat_pos[i] - r / 2.0
    y_pos_max = lat_pos[i] + r / 2.0

    x_neg_min = long_neg[i] - r / 2.0 * 1.0 / np.cos(phi_neg * conv)
    x_neg_max = long_neg[i] + r / 2.0 * 1.0 / np.cos(phi_neg * conv)

    y_neg_min = lat_neg[i] - r / 2.0
    y_neg_max = lat_neg[i] + r / 2.0

    x_pos, y_pos = mask.spot_patch(x, y, x_pos_min, x_pos_max, y_pos_min, y_pos_max, grid_max)

#    spot[date]['xp'].append(x_pos)
#    spot[date]['yp'].append(y_pos)
#    spot_mask[date]['xp'] += x_pos
#    spot_mask[date]['yp'] += y_pos
    spot_mask[date]['xp'] = np.concatenate((spot_mask[date]['xp'], x_pos))
    spot_mask[date]['yp'] = np.concatenate((spot_mask[date]['yp'], y_pos))
                                                                          
    x_neg, y_neg = mask.spot_patch(x, y, x_neg_min, x_neg_max, y_neg_min, y_neg_max, grid_max)

#    spot[date]['xn'].append(x_neg)
#    spot[date]['yn'].append(y_neg)
#    spot_mask[date]['xn'] += x_neg
#    spot_mask[date]['yn'] += y_neg
    spot_mask[date]['xn'] = np.concatenate((spot_mask[date]['xn'], x_neg))
    spot_mask[date]['yn'] = np.concatenate((spot_mask[date]['yn'], y_neg))

sdate = int(min(data[:, 0]))
edate = int(max(data[:, 0]))

spot_mask = {}

for date in range(sdate, edate + 1):

    spot_mask[date] = {'xp': np.array([]), 'yp': np.array([]), 'xn': np.array([]), 'yn': np.array([])}

#pool = mp.Pool(4)
#pool = Pool(4)
#pool.imap_unordered(f, tqdm(range(len(data[:, 0])), total = len(data[:, 0]), \
#ncols = auxfunc.term_width(), desc = 'Masking spots'))

#pool.imap(f, ([i, spot_mask] for i in tqdm(range(1000), total = 1000, \
#ncols = auxfunc.term_width(), desc = 'Masking spots', position = 0)))

#pool.imap(f, ([i, spot_mask] for i in tqdm(range(len(data[:, 0])), \
#total = len(data[:, 0]), ncols = auxfunc.term_width(), desc = 'Masking spots', position = 0)))

#pool.imap(f, tqdm(range(100), total = 100, \
#pool.imap_unordered(f, tqdm(range(100), total = 100, \
#ncols = auxfunc.term_width(), desc = 'Masking spots'))

nproc = 4

if len(sys.argv) == 2:

    nproc = int(sys.argv[1])

with Pool(processes = nproc) as p:

    maximum = len(data[:, 0])

    with tqdm(total = maximum, ncols = auxfunc.term_width(), desc = 'Masking spots', position = 0) as pbar:

        for i, _ in enumerate(p.imap_unordered(f, range(maximum))):

            pbar.update()

    p.close()
    p.join()

np.save('spot_mask.npy', spot_mask)

