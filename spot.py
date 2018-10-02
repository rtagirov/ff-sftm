#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:28:59 2018

@author: nemec
"""
import importlib

import numpy as np
from multiprocessing import Pool

#the next two modules are not needed for the masking, but for the next steps
import matplotlib.pyplot as plt
from scipy import integrate

from tqdm import tqdm

import auxfunc

importlib.reload(auxfunc)

#setting up the grid for the calculations
factor = 10.0

xvalues = np.around(np.linspace(0.0, 359.0, num = 360.0 * factor), decimals = 1)
yvalues = np.around(np.linspace(0.0, 180.0, num = 181.0 * factor), decimals = 1)

x, y = np.meshgrid(xvalues, yvalues)

#for doing the sin/cos calculations
conv = np.pi / 180.0

#load in the new input file that also contains the positions on the days of decay
data = np.loadtxt("spot_evolution.inp")

#defining the coordinates
long_pos = data[:, 2]
long_neg = data[:, 4]

lat_pos = data[:, 1]
lat_neg = data[:, 3]

#spot area
area = data[:, 5]

#define which part of the input should then be used
#start = 70724
#end =   74519

grid_max = 359

#def f(i):
for i in tqdm(range(0, len(data[:, 0])), ncols = auxfunc.term_width(), desc = 'Masking spots'):
#for i in range(0, len(data[:, 0])):

    positivx= open("positivx{}.txt".format(i), "w")
    positivy= open("positivy{}.txt".format(i), "w")
    negativx= open("negativx{}.txt".format(i), "w")
    negativy= open("negativy{}.txt".format(i), "w")

    time = data[i, 0]

#    area = spot[i]
    r = np.sqrt(area[i])

    #for taking into account that pixels become smaller towards the pole
#    phi_pos = 90.0 - data[:, 1]
#    phi_neg = 90.0 - data[:, 3]
    phi_pos = 90.0 - lat_pos[i]
    phi_neg = 90.0 - lat_neg[i]

    x_min_pos = long_pos[i] - r / 2.0 * 1.0 / np.cos(phi_pos * conv)
    x_max_pos = long_pos[i] + r / 2.0 * 1.0 / np.cos(phi_pos * conv)
    y_min_pos = lat_pos[i] - r / 2.0
    y_max_pos = lat_pos[i] + r / 2.0

    #define negative polarity patch
    x_min_neg = long_neg[i] - r / 2.0 * 1.0 / np.cos(phi_neg * conv)
    x_max_neg = long_neg[i] + r / 2.0 * 1.0 / np.cos(phi_neg * conv)
    y_min_neg = lat_neg[i] - r / 2.0
    y_max_neg = lat_neg[i] + r / 2.0

    #the if-statements are just for the cases if part of the spots would be masked
    #outside the grid
    #I need it for both the positive and the negative polarity patch

    #For the next steps, it would be good if we would save the masked pixels into arrays
    #and use the arrays further and not write them out.

    if x_min_pos < 0 and x_max_pos > 0:
        x_min_pos1= grid_max+x_min_pos
        x_pos_pos = x[np.where((x >= x_min_pos1) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos = y[np.where((x >= x_min_pos1) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time)) 
        x_pos_pos1 = x[np.where((x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos1 = y[np.where((x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos1:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos1:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_min_pos < 0. and x_max_pos <0:
        x_min_pos2= grid_max+x_min_pos
        x_max_pos2 = grid_max+x_max_pos
        x_pos_pos2 = x[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos2 = y[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos2:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos2:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_min_pos > 0.:
        x_pos_pos3 = x[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos3 = y[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos3:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos3:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_min_neg < 0. and x_max_neg >0:
        x_min_neg1= grid_max+x_min_neg
        x_pos_neg1 = x[np.where((x >= x_min_neg1) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg= y[np.where((x >= x_min_neg1) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg1:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time)) 
        x_pos_neg1 = x[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg1 = y[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg1:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg1:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_min_neg < 0. and x_max_neg <0:
        x_min_neg2= grid_max+x_min_neg
        x_max_neg2 = grid_max +x_max_neg
        x_pos_neg2 = x[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg2 = y[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg2:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg2:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_min_neg > 0.:
        x_pos_neg3 = x[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg3 = y[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg3:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg3:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_max_pos >grid_max and x_min_pos <grid_max:
        x_max_pos4= x_max_pos-grid_max
        x_pos_pos = x[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos = y[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time)) 
        x_pos_pos4 = x[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos4 = y[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos4:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos4:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_max_pos >grid_max and x_min_pos >grid_max:
        x_min_pos5= x_min_pos-grid_max
        x_max_pos5 =x_max_pos-grid_max
        x_pos_pos5 = x[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos5 = y[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos5:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos5:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))

    if x_max_pos > grid_max:
        x_pos_pos6 = x[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        y_pos_pos6 = y[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
        for item in x_pos_pos6:
            positivx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_pos6:
            positivy.write("%f  \t %f \t %f\n" % (item,r,time))
    
    if x_max_neg >grid_max and x_min_neg <grid_max:
        x_max_neg4= x_max_pos-grid_max
        x_pos_neg = x[np.where((x >= x_min_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg = y[np.where((x >= x_min_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time)) 
        x_pos_neg4 = x[np.where((x<=x_max_neg4) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg4 = y[np.where((x<=x_max_neg4) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg4:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg4:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time))
            
    if x_max_neg >grid_max and x_min_neg >grid_max:
        x_min_neg5= x_min_neg-grid_max
        x_max_neg5 =x_max_neg-grid_max
        x_pos_neg5 = x[np.where((x >= x_min_neg5) & (x<=x_max_neg5) & (y>=y_min_neg) & (y <=y_max_neg))]
        y_pos_neg5 = y[np.where((x >= x_min_neg5) & (x<=x_max_neg5) & (y>=y_min_neg) & (y <=y_max_neg))]
        for item in x_pos_neg5:
            negativx.write("%f  \t %f \t %f\n" % (item,r,time))
        for item in y_pos_neg5:
            negativy.write("%f  \t %f \t %f\n" % (item,r,time))

#pool = Pool(4)
#pool.map(f, range(start, end + 1))

