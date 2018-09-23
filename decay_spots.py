#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 19:20:14 2018

@author: nemec
"""

import numpy as np
from multiprocessing import Pool

#calculating the life time of the spots according to the choosen decay rate

def decay(spot_area,time,D):
    t = spot_area/D +time
    return t

#calculate the meridional flow
def merflow(lat):
    if abs(lat-90) <= 75:
        u = 22*np.sin(2.4*(90-lat)*np.pi/180)*7.0922e-3
    if abs(lat-90) > 75:
        u = 0
    return u

#calculate the differential rotation
def diffrot(lat):
    rot = 0.1813 - 2.3*np.sin((90-lat)*np.pi/180)**2.-1.62*np.sin((90-lat)*np.pi/180)**4.
    return rot


#define the decay rate
D = 30.9 #MHS per day


#setting up the grid on which to mask the spots
#in this case pixels have size of 0.1 by 0.1 degree
factor = 10.
xvalues = np.around(np.linspace(0,359,num=360*factor),decimals=1)
yvalues = np.around(np.linspace(0,180,num=180*factor),decimals=1)
x,y = np.meshgrid(xvalues,yvalues)

#for doing the sin/cos calculations
conv = np.pi/180.


# =============================================================================
# if you want to experiment with random parameters for the positions and areas,
#just comment out the line, where I read in the input file and define the coordinates
# and area yourself
# =============================================================================
#reading in the file
data = np.loadtxt("AR-mod.txt")

#defining the coordinates
long_pos = data[:, 2]
long_neg = data[:,4]

#need to redifine grid so that north pole is a + 90 degree and south pole at -90 degree
lat_pos = 90-data[:,1]
lat_neg = 90 -data[:,3]

#define the area at the time of emergence
spot = data[:,5]

#define which part of the input should then be used
start = 70724
end =   74519
grid_max = 359


#only use this for calculations that should not be run parallel!
#positivx= open("positivx3.txt","w")
#positivy= open("positivy3.txt","w")
#negativx= open("negativx3.txt","w")
#negativy= open("negativy3.txt","w")

#for i in range(start,end):

#starting doing the spot masking parallel
def f(i):
    #for i in range(start,end):
    #print(i)
    positivx= open("positivx{}.txt".format(i),"w")
    positivy= open("positivy{}.txt".format(i),"w")
    negativx= open("negativx{}.txt".format(i),"w")
    negativy= open("negativy{}.txt".format(i),"w")
    
    spot_area = spot[i]
    time = data[i,0]
    t = decay(spot_area,time,D)
    phi_pos = 90-lat_pos[i]
    phi_neg = 90-lat_neg[i]
	
    #print(t)
    if np.int(t-time) == 0:
        area = spot[i]/(30.81*np.pi/2.*np.pi/2.)
        r = area**(1./2.)
        
        #define positive polarity patch
        x_min_pos = long_pos[i]-r/2.*1./np.cos(phi_pos*conv)
        x_max_pos = long_pos[i]+r/2.*1./np.cos(phi_pos*conv)
        y_min_pos = lat_pos[i]-r/2.
        y_max_pos = lat_pos[i]+r/2.
       
        #define negative polarity patch
        x_min_neg = long_neg[i]-r/2.*1./np.cos(phi_neg*conv)
        x_max_neg = long_neg[i]+r/2.*1./np.cos(phi_neg*conv)
        y_min_neg = lat_neg[i]-r/2.
        y_max_neg = lat_neg[i]+r/2.
        
        if x_min_pos < 0 and x_max_pos >0:
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
#           
                        
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
#           
                        
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

        

    if np.int(t-time) > 0:
        lat_pos_new = np.zeros(np.int(t-time)+1)
        lat_neg_new = np.zeros(np.int(t-time)+1)
        long_pos_new = np.zeros(np.int(t-time)+1)
        long_neg_new = np.zeros(np.int(t-time)+1)
        lat_pos_new[0]= lat_pos[i]
        lat_neg_new[0]= lat_neg[i]
        long_pos_new[0]= long_pos[i]
        long_neg_new[0]= long_neg[i]
        n = time

        for n in range(np.int(time),np.int(t)+1):
            #update the are according to the decay law
            #in that case it's a linear!
            area =(spot[i]-D*(n-time))/(30.81*np.pi/2.*np.pi/2.) #in degree
            if area <= 0.:
                r = 0.
            else:
                r = area**(1./2.)
            
            if n == time:
        #define positive polarity patch
                 x_min_pos = long_pos[i]-r/2.*1./np.cos(phi_pos*conv)
                 x_max_pos = long_pos[i]+r/2.*1./np.cos(phi_pos*conv)
                 y_min_pos = lat_pos[i]-r/2.
                 y_max_pos = lat_pos[i]+r/2.
       
        #define negative polarity patch
                 x_min_neg = long_neg[i]-r/2.*1./np.cos(phi_neg*conv)
                 x_max_neg = long_neg[i]+r/2.*1./np.cos(phi_neg*conv)
                 y_min_neg = lat_neg[i]-r/2.
                 y_max_neg = lat_neg[i]+r/2.

                 
                 if x_min_pos < 0 and x_max_pos >0:
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
                     x_max_pos2 = grid_max +x_max_pos
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
                     #print(x_min_neg)
                     x_min_neg1= grid_max+x_min_neg
                     #print(x_min_neg1)
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

            if n >time:
                a = np.int(n-time-1)
                b = np.int(n-time)
                
                #I think that this if-statements could be skipped and the positions just be
                #updated, but I need to think about that again....
                if lat_pos_new[a] <= 90:
                    u = merflow(lat_pos_new[a])
                    lat_pos_new[b] = lat_pos_new[a]-u
                    rot= diffrot(lat_pos_new[a])
                    long_pos_new[b] = long_pos_new[a]+rot
                if lat_neg_new[a] <= 90:
                    u = merflow(lat_neg_new[a])
                    lat_neg_new[b] = lat_neg_new[a]-u
                    rot= diffrot(lat_neg_new[a])
                    long_neg_new[b] = long_neg_new[a]+rot
                if lat_pos_new[a] > 90:
                    u = merflow(lat_pos_new[a])
                    lat_pos_new[b] = lat_pos_new[a]-u
                    rot= diffrot(lat_pos_new[a])
                    long_pos_new[b] = long_pos_new[a]+rot
                if lat_neg_new[a] > 90:
                    u = merflow(lat_neg_new[a])
                    lat_neg_new[b] = lat_neg_new[a]-u
                    rot= diffrot(lat_neg_new[a])
                    long_neg_new[b] = long_neg_new[a]+rot
           
                
                    
                x_min_pos = long_pos_new[b]-r/2.*1./np.cos(phi_pos*conv)
                x_max_pos = long_pos_new[b]+r/2.*1./np.cos(phi_pos*conv)
                y_min_pos = lat_pos_new[b]-r/2.
                y_max_pos = lat_pos_new[b]+r/2.
                
                x_min_neg = long_neg_new[b]-r/2.*1./np.cos(phi_neg*conv)
                x_max_neg = long_neg_new[b]+r/2.*1./np.cos(phi_neg*conv)
                y_min_neg = lat_neg_new[b]-r/2.
                y_max_neg = lat_neg_new[b]+r/2.
                    
                
                 
                if x_min_pos < 0 and x_max_pos >0:
                    x_min_pos1= grid_max+x_min_pos
                    x_pos_pos = x[np.where((x >= x_min_pos1) & (y>=y_min_pos) & (y <=y_max_pos))]
                    y_pos_pos = y[np.where((x >= x_min_pos1) & (y>=y_min_pos) & (y <=y_max_pos))]
                    for item in x_pos_pos:
                        positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_pos:
                        positivy.write("%f  \t %f \t %f\n" % (item,r,n)) 
                    x_pos_pos1 = x[np.where((x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                    y_pos_pos1 = y[np.where((x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                    for item in x_pos_pos1:
                        positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_pos1:
                        positivy.write("%f  \t %f \t %f\n" % (item,r,n))
                    
                
                      
                if x_min_pos < 0. and x_max_pos <0:
                    x_min_pos2= grid_max+x_min_pos
                    x_max_pos2 =grid_max+x_max_pos
                    x_pos_pos2 = x[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]
                    y_pos_pos2 = y[np.where((x >= x_min_pos2) & (x<=x_max_pos2) & (y>=y_min_pos) & (y <=y_max_pos))]
                    for item in x_pos_pos2:
                        positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_pos2:
                        positivy.write("%f  \t %f \t %f\n" % (item,r,n))
#           
                        
                if x_min_pos > 0.:
                    x_pos_pos3 = x[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                    y_pos_pos3 = y[np.where((x >= x_min_pos) & (x<=x_max_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                    for item in x_pos_pos3:
                        positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_pos3:
                        positivy.write("%f  \t %f \t %f\n" % (item,r,n))

            
                if x_min_neg < 0. and x_max_neg >0:
                    #print(x_min_neg,n)
                    #print(x_max_neg)
                    x_min_neg1= grid_max+x_min_neg
                    x_pos_neg1 = x[np.where((x >= x_min_neg1) & (y>=y_min_neg) & (y <=y_max_neg))]
                    y_pos_neg= y[np.where((x >= x_min_neg1) & (y>=y_min_neg) & (y <=y_max_neg))]
                    for item in x_pos_neg1:
                        negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_neg:
                        negativy.write("%f  \t %f \t %f\n" % (item,r,n)) 
                    x_pos_neg1 = x[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                    y_pos_neg1 = y[np.where((x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                    for item in x_pos_neg1:
                        negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_neg1:
                        negativy.write("%f  \t %f \t %f\n" % (item,r,n))
                        
                if x_min_neg < 0. and x_max_neg <0:
                   # print(x_min_neg,n)
                    x_min_neg2= grid_max+x_min_neg
                    x_max_neg2 = grid_max +x_max_neg
                    #print(x_min_neg2,x_max_neg2,n)
                    x_pos_neg2 = x[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]
                    y_pos_neg2 = y[np.where((x >= x_min_neg2) & (x<=x_max_neg2) & (y>=y_min_neg) & (y <=y_max_neg))]
                    for item in x_pos_neg2:
                        negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_neg2:
                        negativy.write("%f  \t %f \t %f\n" % (item,r,n))
                        
                if x_min_neg > 0.:
                    x_pos_neg3 = x[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                    y_pos_neg3 = y[np.where((x >= x_min_neg) & (x<=x_max_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                    for item in x_pos_neg3:
                        negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_neg3:
                        negativy.write("%f  \t %f \t %f\n" % (item,r,n))
                
                if x_max_pos >grid_max and x_min_pos <grid_max:
                     x_max_pos4= x_max_pos-grid_max
                     x_pos_pos = x[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                     y_pos_pos = y[np.where((x >= x_min_pos) & (y>=y_min_pos) & (y <=y_max_pos))]
                     for item in x_pos_pos:
                         positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                     for item in y_pos_pos:
                         positivy.write("%f  \t %f \t %f\n" % (item,r,n)) 
                     x_pos_pos4 = x[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]
                     y_pos_pos4 = y[np.where((x<=x_max_pos4) & (y>=y_min_pos) & (y <=y_max_pos))]
                     for item in x_pos_pos4:
                         positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                     for item in y_pos_pos4:
                         positivy.write("%f  \t %f \t %f\n" % (item,r,n))
                    
                
                if x_max_pos >grid_max and x_min_pos >grid_max:
                    x_min_pos5= x_min_pos-grid_max
                    x_max_pos5 =x_max_pos-grid_max
                    x_pos_pos5 = x[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]
                    y_pos_pos5 = y[np.where((x >= x_min_pos5) & (x<=x_max_pos5) & (y>=y_min_pos) & (y <=y_max_pos))]
                    for item in x_pos_pos5:
                        positivx.write("%f  \t %f \t %f\n" % (item,r,n))
                    for item in y_pos_pos5:
                        positivy.write("%f  \t %f \t %f\n" % (item,r,n))

        
                if x_max_neg >grid_max and x_min_neg <grid_max:
                     x_max_neg4= x_max_pos-grid_max
                     x_pos_neg = x[np.where((x >= x_min_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                     y_pos_neg = y[np.where((x >= x_min_neg) & (y>=y_min_neg) & (y <=y_max_neg))]
                     for item in x_pos_neg:
                         negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                     for item in y_pos_neg:
                         negativy.write("%f  \t %f \t %f\n" % (item,r,n)) 
                     x_pos_neg4 = x[np.where((x<=x_max_neg4) & (y>=y_min_neg) & (y <=y_max_neg))]
                     y_pos_neg4 = y[np.where((x<=x_max_neg4) & (y>=y_min_neg) & (y <=y_max_neg))]
                     for item in x_pos_neg4:
                         negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                     for item in y_pos_neg4:
                         negativy.write("%f  \t %f \t %f\n" % (item,r,n))
                    
                
                if x_max_neg >grid_max and x_min_neg >grid_max:
                     x_min_neg5= x_min_neg-grid_max
                     x_max_neg5 =x_max_neg-grid_max
                     x_pos_neg5 = x[np.where((x >= x_min_neg5) & (x<=x_max_neg5) & (y>=y_min_neg) & (y <=y_max_neg))]
                     y_pos_neg5 = y[np.where((x >= x_min_neg5) & (x<=x_max_neg5) & (y>=y_min_neg) & (y <=y_max_neg))]
                     for item in x_pos_neg5:
                          negativx.write("%f  \t %f \t %f\n" % (item,r,n))
                     for item in y_pos_neg5:
                          negativy.write("%f  \t %f \t %f\n" % (item,r,n))
#           

pool = Pool(4)
pool.map(f, range(start, end + 1))
#                
#           
#           
        
        
        
                
                
        
        
        
        
        
        
    


