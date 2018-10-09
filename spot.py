import importlib
import auxfunc
import sys

import numpy as np
import matplotlib.pyplot as plt

from multiprocessing import Pool
from tqdm import tqdm

importlib.reload(auxfunc)

#setting up the grid for the calculations
factor = 10.0

xvalues = np.around(np.linspace(0.0, 359.0, num = 360.0 * factor), decimals = 1)
yvalues = np.around(np.linspace(0.0, 180.0, num = 181.0 * factor), decimals = 1)

x, y = np.meshgrid(xvalues, yvalues)

#for doing the sin / cos calculations
conv = np.pi / 180.0

#load in the new input file that also contains the positions on the days of decay
data = np.loadtxt('spot_evolution.inp')

#defining the coordinates
long_pos = data[:, 2]
long_neg = data[:, 4]

lat_pos = data[:, 1]
lat_neg = data[:, 3]

grid_max = 359

spot = {}

unique_dates = np.unique(data[:, 0].astype(int))

for idx, date in np.ndenumerate(unique_dates):

    spot[date] = {'xp': [], 'yp': [], 'xn': [], 'yn': []}

#sys.exit()

#def f(i):
for i in tqdm(range(0, len(data[:, 0])), ncols = auxfunc.term_width(), desc = 'Masking spots'):
#for i in range(0, len(data[:, 0])):

    date = int(data[i, 0])

    r = np.sqrt(data[i, 5])

#   for taking into account that pixels become smaller towards the pole
    phi_pos = 90.0 - lat_pos[i]
    phi_neg = 90.0 - lat_neg[i]

#   constrain negative polarity patch
    x_min_pos = long_pos[i] - r / 2.0 * 1.0 / np.cos(phi_pos * conv)
    x_max_pos = long_pos[i] + r / 2.0 * 1.0 / np.cos(phi_pos * conv)
    y_min_pos = lat_pos[i] - r / 2.0
    y_max_pos = lat_pos[i] + r / 2.0

#   constrain negative polarity patch
    x_min_neg = long_neg[i] - r / 2.0 * 1.0 / np.cos(phi_neg * conv)
    x_max_neg = long_neg[i] + r / 2.0 * 1.0 / np.cos(phi_neg * conv)
    y_min_neg = lat_neg[i] - r / 2.0
    y_max_neg = lat_neg[i] + r / 2.0

#   condition for a spot crossing the left border of the disk
    if x_min_pos < 0 and x_max_pos > 0:

#       move the left part of the spot to the right side of the grid
        x_min_pos1 = grid_max + x_min_pos

#        x_pos_pos  = x[np.where((x >= x_min_pos1) & (x <= x_max_pos))]
#        y_pos_pos  = y[np.where((y >= y_min_pos)  & (y <= y_max_pos))]

#       where function for the left part of the spot
        x_pos_pos  = x[np.where((x >= x_min_pos1) & (y >= y_min_pos) & (y <= y_max_pos))]
        y_pos_pos  = y[np.where((x >= x_min_pos1) & (y >= y_min_pos) & (y <= y_max_pos))]

        spot[date]['xp'].append(x_pos_pos)
        spot[date]['yp'].append(y_pos_pos)

#       where function for the right part of the spot
        x_pos_pos1 = x[np.where((x <= x_max_pos) & (y >= y_min_pos) & (y <= y_max_pos))]
        y_pos_pos1 = y[np.where((x <= x_max_pos) & (y >= y_min_pos) & (y <= y_max_pos))]

        spot[date]['xp'].append(x_pos_pos1)
        spot[date]['yp'].append(y_pos_pos1)

#   outside of the grid on the left
    if x_min_pos < 0. and x_max_pos < 0:

        x_min_pos2 = grid_max + x_min_pos
        x_max_pos2 = grid_max + x_max_pos

        x_pos_pos2 = x[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos2 = y[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos2:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos2:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xp'].append(x_pos_pos2)
        spot[date]['yp'].append(y_pos_pos2)

#   normal condition for the positive patch
#   (change: should be bigger or equal and we also say that x_max_pos <= grid_max)
    if x_min_pos > 0.:

        x_pos_pos3 = x[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos3 = y[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos3:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos3:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xp'].append(x_pos_pos3)
        spot[date]['yp'].append(y_pos_pos3)

#   left boundary for the negative patch is outside the grid
    if x_min_neg < 0. and x_max_neg > 0:

        x_min_neg1 = grid_max + x_min_neg

        x_pos_neg1 = x[np.where((x >= x_min_neg1) & (y >= y_min_neg) & (y <= y_max_neg))]
        y_pos_neg  = y[np.where((x >= x_min_neg1) & (y >= y_min_neg) & (y <= y_max_neg))]

#        for item in x_pos_neg1:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg1)
        spot[date]['yn'].append(y_pos_neg)

        x_pos_neg1 = x[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg1 = y[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]

#        for item in x_pos_neg1:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg1:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg1)
        spot[date]['yn'].append(y_pos_neg1)

#   the whole negative patch is outside the grid
    if x_min_neg < 0. and x_max_neg < 0:

        x_min_neg2 = grid_max + x_min_neg
        x_max_neg2 = grid_max + x_max_neg

        x_pos_neg2 = x[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg2 = y[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]

#        for item in x_pos_neg2:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg2:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg2)
        spot[date]['yn'].append(y_pos_neg2)

#   normal condition for the negative patch
#   (change: should be bigger or equal and we also say that x_max_neg <= grid_max)
    if x_min_neg > 0.:

        x_pos_neg3 = x[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg3 = y[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]

#        for item in x_pos_neg3:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg3:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg3)
        spot[date]['yn'].append(y_pos_neg3)

#   right boundary of the positive patch is outside the grid
    if x_max_pos > grid_max and x_min_pos < grid_max:

        x_max_pos4 = x_max_pos - grid_max

        x_pos_pos = x[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos = y[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time)) 

        spot[date]['xp'].append(x_pos_pos)
        spot[date]['yp'].append(y_pos_pos)

        x_pos_pos4 = x[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos4 = y[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos4:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos4:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xp'].append(x_pos_pos4)
        spot[date]['yp'].append(y_pos_pos4)

    if x_max_pos > grid_max and x_min_pos > grid_max:

        x_min_pos5 = x_min_pos - grid_max
        x_max_pos5 = x_max_pos - grid_max

        x_pos_pos5 = x[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos5 = y[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos5:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos5:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xp'].append(x_pos_pos5)
        spot[date]['yp'].append(y_pos_pos5)

#   redundant condition
    if x_max_pos > grid_max:

        x_pos_pos6 = x[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos6 = y[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]

#        for item in x_pos_pos6:
#            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_pos6:
#            positivy.write("%f  \t %f \t %f\n" % (item,r,time))
    
        spot[date]['xp'].append(x_pos_pos6)
        spot[date]['yp'].append(y_pos_pos6)

    if x_max_neg > grid_max and x_min_neg < grid_max:

        x_max_neg4 = x_max_pos - grid_max

        x_pos_neg  = x[np.where((x >= x_min_neg) & (y >= y_min_neg) & (y <= y_max_neg))]
        y_pos_neg  = y[np.where((x >= x_min_neg) & (y >= y_min_neg) & (y <= y_max_neg))]

#        for item in x_pos_neg:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg)
        spot[date]['yn'].append(y_pos_neg)

        x_pos_neg4 = x[np.where((x <= x_max_neg4) & (y >= y_min_neg) & (y <= y_max_neg))]
        y_pos_neg4 = y[np.where((x <= x_max_neg4) & (y >= y_min_neg) & (y <= y_max_neg))]

#        for item in x_pos_neg4:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg4:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg4)
        spot[date]['yn'].append(y_pos_neg4)

    if x_max_neg > grid_max and x_min_neg > grid_max:

        x_min_neg5 = x_min_neg - grid_max
        x_max_neg5 = x_max_neg - grid_max

        x_pos_neg5 = x[np.where((x >= x_min_neg5) & (x <= x_max_neg5) & (y >= y_min_neg) & (y <= y_max_neg))]
        y_pos_neg5 = y[np.where((x >= x_min_neg5) & (x <= x_max_neg5) & (y >= y_min_neg) & (y <= y_max_neg))]

#        for item in x_pos_neg5:
#            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
#        for item in y_pos_neg5:
#            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

        spot[date]['xn'].append(x_pos_neg5)
        spot[date]['yn'].append(y_pos_neg5)

#pool = Pool(4)
#pool.map(f, range(start, end + 1))

