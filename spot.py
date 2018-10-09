import importlib
import auxfunc
import mask
import sys

import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool
from tqdm import tqdm

importlib.reload(auxfunc)
importlib.reload(mask)

factor = 10

xvalues = np.around(np.linspace(0.0, 359.0, num = 360 * factor), decimals = 1)
yvalues = np.around(np.linspace(0.0, 180.0, num = 181 * factor), decimals = 1)

x, y = np.meshgrid(xvalues, yvalues)

#sys.exit()

conv = np.pi / 180.0

data = np.loadtxt('spot_evolution.inp')

long_pos = data[:, 2]
long_neg = data[:, 4]

lat_pos = data[:, 1]
lat_neg = data[:, 3]

grid_max = 359

sdate = 100893
edate = 104479

spot = {}

#unique_dates = np.unique(data[:, 0].astype(int))

for date in range(sdate, edate + 1):

    spot[date] = {'xp': [], 'yp': [], 'xn': [], 'yn': []}

#sys.exit()

#for i in tqdm(range(len(data[:, 0])), desc = 'Masking spots'):
#for i in tqdm(range(len(data[:, 0])), ncols = auxfunc.term_width(), desc = 'Masking spots'):
for i in tqdm(range(len(data[:, 0])), ncols = auxfunc.term_width(), desc = 'Masking spots', leave = True):
#for i, elem in np.ndenumerate(data[:, 0]):

#    if i > 15:

#        sys.exit()

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

    x_pos, y_pos = mask.patch(x, y, x_pos_min, x_pos_max, y_pos_min, y_pos_max, grid_max)
                                                                          
    spot[date]['xp'].append(x_pos)                                        
    spot[date]['yp'].append(y_pos)                                        
                                                                          
    x_neg, y_neg = mask.patch(x, y, x_neg_min, x_neg_max, y_neg_min, y_neg_max, grid_max)

    spot[date]['xn'].append(x_neg)
    spot[date]['yn'].append(y_neg)

