import importlib
import auxfunc
import mask
import sys
import math

import numpy as np

from tqdm import tqdm
from multiprocessing import Pool

importlib.reload(auxfunc)
importlib.reload(mask)

nproc = 4

inp = 'input_C22_26.5'

args = sys.argv[1:]

for i, arg in enumerate(args):

    if arg == '--inp':

        inp = args[i + 1]

    if arg == '--np':

        nproc = int(args[i + 1])

D = inp.split('_')[2]

#sys.exit()

factor = 10

xvalues = np.around(np.linspace(0.0, 359.0, num = 360 * factor), decimals = 1)
yvalues = np.around(np.linspace(0.0, 180.0, num = 181 * factor), decimals = 1)

x, y = np.meshgrid(xvalues, yvalues)

conv = np.pi / 180.0

data = np.loadtxt(inp)

long_pos = data[:, 2]
long_neg = data[:, 4]

lat_pos = data[:, 1]
lat_neg = data[:, 3]

grid_max = 359

def line_contrib(i):

#    date = data[i, 0]

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

#    return date, x_pos, y_pos, x_neg, y_neg
    return x_pos, y_pos, x_neg, y_neg

#sdate = min(data[:, 0])
#edate = max(data[:, 0])

times = data[:, 0]

spot_mask = {}

#for date in range(sdate, edate + 1):
for time in times:

    spot_mask[time] = {'xp': np.array([]), 'yp': np.array([]), 'xn': np.array([]), 'yn': np.array([])}

#sys.exit()

with Pool(processes = nproc) as p:

    maximum = len(times)

#    n_chunks = math.ceil(maximum / nproc)
#    n_chunks = 100
    n_chunks = 1

    with tqdm(total = maximum, \
              ncols = auxfunc.term_width(), \
              desc = 'Masking spots, nproc = ' + str(nproc), \
              position = 0) as pbar:

        results = p.imap(line_contrib, range(maximum), chunksize = n_chunks)

        for i, result in enumerate(results):

#            date, x_pos, y_pos, x_neg, y_neg = result
            x_pos, y_pos, x_neg, y_neg = result

            time = times[i]

            spot_mask[time]['xp'] = np.concatenate((spot_mask[time]['xp'], x_pos))
            spot_mask[time]['yp'] = np.concatenate((spot_mask[time]['yp'], y_pos))

            spot_mask[time]['xn'] = np.concatenate((spot_mask[time]['xn'], x_neg))
            spot_mask[time]['yn'] = np.concatenate((spot_mask[time]['yn'], y_neg))

            pbar.update()

    p.close()
    p.join()

spot_mask = np.load('spot_mask_D' + D + '.npy').item()

