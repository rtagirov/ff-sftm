#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 11:45:16 2018

@author: nemec
"""

# =============================================================================
# importing all the packages
import numpy as np
import glob
import re
# =============================================================================

#defining all the constants
norm = 90*90*4/np.pi**2*np.pi #we need to normalize the areas on the sphere
area = 0.1**2. # I don't know if I ever use this....
conv = np.pi/180.
# =============================================================================

#set position of observer and Bsat
#x_c is the central longitude of the observer
#y_c is the central latitude of the observer
# 90 is the north pole, -90 the south pole
x_c = 0
y_c = 0
#we might need to experiment with this two values more, but for now, we just take
#those two as fixed
B_sat = 484
B_spot = 1500
# =============================================================================

#set up the grid for the mu calculations
#I think its probably better if Sasha explains you at same point, how SATIRE works
mu_grid = [1.0000, 0.9000, 0.8000, 0.7000, 0.6000, 0.5000, 0.4000, 0.3000 ,0.2000, 0.1000, 0.0500]
mu_low = [0.95,0.85,0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0]
mu_up = [1,0.95,0.85,0.75,0.65,0.55,0.45,0.35,0.25,0.15,0.075]

# =============================================================================
#reading in the files
#the filenames are saved into "files"
files = sorted(glob.glob("CalcMagnetogram.2000.*"))
spotx = np.loadtxt("pos_x.txt")
spoty = np.loadtxt("pos_y.txt")
start = 100749

g = open("fac_old"+str(B_sat)+".txt","w")
for files in files:
    name = re.findall("2000.(\d+)",files)
    visibility = []
    #I'm setting up all kinds of lists that are needed for the next steps
    a = []
    b = []
    c = []
    mu1 = []
    mu2 = []
    mu3 = []
    mu4 = []
    mu5 = []
    mu6 = []
    mu7 = []
    mu8 = []
    mu9 = []	
    mu10 = []
    mu11 = []
    factor = []
    data = np.loadtxt(files) #this array now contains the information of the flux in each pixel
    #of the magnetogram, so the array has 360X181 entries
    time= int(name[0])
    ff = np.zeros((180,360))
    spot_x= spotx[np.where(int(name[0])==spotx[:,2]),0]
    spot_y= spoty[np.where(int(name[0])==spoty[:,2]),0]
    #if np.shape(spot_x) == np.shape(spot_y):
    for i in range(0,180): #for latitudes
        for j in range(0,360): #for longitudes
            #setting up the boundaries of the pixels
            x_min = j
            x_max = j+1
            y_min = i
            y_max = i+1
            #look for spot-minipixels in the 1x1 deg pixels 
            spot = (spot_x[np.where((x_min <= spot_x) & (spot_x < x_max) & (y_min<= spot_y) & (spot_y <  y_max))])
            B = abs(data[i][j])
            #here we have two methods to calculate the faculae filling factors, maybe stick to that one for now
            if np.shape(spot) !=(0,):
               # ff[i,j] =(1 - len(spot)*0.1*0.1)
               helper =B- B_spot*len(spot)*0.1*0.1
               if helper < 0:
                    ff[i,j] = 0
               else:
                    ff[i,j] = helper/B_sat
                
                   # a.append(1-len(spot)*0.1*0.1)
            if np.shape(spot) == (0,) and B < B_sat:
                ff[i,j] =abs(data[i][j])/B_sat
            if np.shape(spot) == (0,) and B >= B_sat:
                ff[i,j]=1
            
            #here I am rotating the grid and calculate the mu-positions
            x_rot = []   
            conv = np.pi/180.        
            x_rot=(j+13.28*(int(time)-start))%359
            x_pos = 180-x_rot
            y_pos = 90-i
            delta_lambda = abs(x_pos-x_c)
            distance = np.arccos((np.sin((y_c)*conv)*np.sin((y_pos)*conv)+np.cos((y_c)*conv) \
                              *np.cos((y_pos)*conv)*np.cos(delta_lambda*conv)))/conv
            vis = np.cos(distance*conv)
                #print(distance)
                #if distance <=90:
                 #   visibility.append(ff[i,j]*np.cos(distance*conv)*np.cos(y_pos*conv))
            #if vis <=1 and vis >= 0 and   np.shape(spot) !=(0,): #I don't remember what these two
            #lines were for....
               # a.append(1-len(spot)*0.1*0.1)
            if vis <=mu_up[0] and vis >mu_low[0]:
                mu1.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[1] and vis >mu_low[1]:
                mu2.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[2] and vis >mu_low[2]:
                mu3.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[3] and vis >mu_low[3]:
                mu4.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[4] and vis >mu_low[4]:
                mu5.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[5] and vis >mu_low[5]:
                mu6.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[6] and vis >mu_low[6]:
                mu7.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[7] and vis >mu_low[7]:
                mu8.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[8] and vis >mu_low[8]:
                mu9.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[9] and vis >mu_low[9]:
                mu10.append(ff[i,j]*vis*np.cos(y_pos*conv))
            if vis <=mu_up[10] and vis >mu_low[10]:
                mu11.append(ff[i,j]*vis*np.cos(y_pos*conv))
            #the following one is just for testing, if the mu-ring distribution gives the same
            #result, as if we would just use the full "visible" disc
           # if distance <=90:
               # visibility.append(ff[i,j]*np.cos(distance*conv)*np.cos(y_pos*conv))
                    
    r1=sum(mu1)/norm
    r2=sum(mu2)/norm
    r3= sum(mu3)/norm
    r4=sum(mu4)/norm
    r5=sum(mu5)/norm
    r6=sum(mu6)/norm
    r7=sum(mu7)/norm
    r8=sum(mu8)/norm
    r9=sum(mu9)/norm
    r10=sum(mu10)/norm
    r11=sum(mu11)/norm
    total = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11
    #print(a)
    g.write("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f  \n" %(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,time,total))
 #   f.write("%f \t \n" %(sum(a)))
#f.close()


g.close()

         


        
