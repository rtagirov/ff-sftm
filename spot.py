import importlib
import auxfunc
import mask
import sys

import numpy as np
import multiprocessing as mp

from tqdm import tqdm

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

def line_contrib(i):

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
    x_neg, y_neg = mask.spot_patch(x, y, x_neg_min, x_neg_max, y_neg_min, y_neg_max, grid_max)

    lock.acquire()

    spot_mask[date]['xp'] = np.concatenate((spot_mask[date]['xp'], x_pos))
    spot_mask[date]['yp'] = np.concatenate((spot_mask[date]['yp'], y_pos))

    spot_mask[date]['xn'] = np.concatenate((spot_mask[date]['xn'], x_neg))
    spot_mask[date]['yn'] = np.concatenate((spot_mask[date]['yn'], y_neg))

    lock.release()

def init(l):

    global lock

    lock = l

sdate = int(min(data[:, 0]))
edate = int(max(data[:, 0]))

spot_mask = {}

for date in range(sdate, edate + 1):

    spot_mask[date] = {'xp': np.array([]), 'yp': np.array([]), 'xn': np.array([]), 'yn': np.array([])}

nproc = 4

if len(sys.argv) == 2:

    nproc = int(sys.argv[1])

l = mp.Lock()

with Pool(processes = nproc, initializer = init, initargs = (l,)) as p:

    maximum = len(data[:, 0])

    with tqdm(total = maximum, ncols = auxfunc.term_width(), desc = 'Masking spots', position = 0) as pbar:

        for i, _ in enumerate(p.imap(line_contrib, range(maximum))):

            pbar.update()

    p.close()
    p.join()

np.save('spot_mask.npy', spot_mask)

