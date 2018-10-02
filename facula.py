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
norm = 90*90*4/np.pi**2*np.pi
area = 0.1**2.
conv = np.pi/180.
# =============================================================================

#set position of observer and Bsat
x_c = 0
y_c = 0
B_sat = 484.
B_spot = 1000
# =============================================================================

#set up the grid for the mu calculations
mu_grid = [1.0000, 0.9000, 0.8000, 0.7000, 0.6000, 0.5000, 0.4000, 0.3000 ,0.2000, 0.1000, 0.0500]
mu_low = [0.95,0.85,0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0]
mu_up = [1,0.95,0.85,0.75,0.65,0.55,0.45,0.35,0.25,0.15,0.075]

# =============================================================================
#reading in the files
files = sorted(glob.glob("CalcMagnetogram.2000.*"))
spotx = np.loadtxt("pos_x.txt")
spoty = np.loadtxt("pos_y.txt")
start = 170049

#g = open("fac_"+str(x_c)+str(y_c)+str(B_sat)+str(B_spot)+".txt","w")
#f = open("facstats_old.txt","w")
g = open("fac_spots_484","w")
for files in files:
    name = re.findall("2000.(\d+)",files)
    visibility = []
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
  #  print(name)
    data = np.loadtxt(files)
    time= int(name[0])
   # print(time)
    #factor.extend(abs(data/B_sat))
    ff = np.zeros((180,360))
    spot_x= spotx[np.where(int(name[0])==spotx[:,2]),0]
    spot_y= spoty[np.where(int(name[0])==spoty[:,2]),0]
        #if np.shape(spot_x) == np.shape(spot_y):
    for i in range(0,180): #for latitudes
        for j in range(0,360): #for longitudes
    #distinguish now between the 3 different kind of pixels": spot+faculae, pure faculae, faculae+QS  
            x_min = j
            x_max = j+1
            y_min = i
            y_max = i+1
            spot = (spot_x[np.where((x_min <= spot_x) & (spot_x < x_max) & (y_min<= spot_y) & (spot_y <  y_max))])#!!!!!!!!
            B = abs(data[i][j])
            if np.shape(spot) !=(0,):
#            #ff[i,j] =(1 - len(spot)*0.1*0.1)
                helper = B - B_spot*len(spot)*0.1*0.1
                if helper <= 0:
                    ff[i,j] = 0
                if helper > 0:
                    ff[i,j] = (1-len(spot)*0.1*0.1)*helper/B_sat
#               # helper =B- B_spot*len(spot)*0.1*0.1
               #if helper < 0:
               #    ff[i,j] = 0
               # else:
               #     ff[i,j] = helper/B_sat
                
                   # a.append(1-len(spot)*0.1*0.1)
            if np.shape(spot) == (0,) and B < B_sat:
                ff[i,j] =abs(data[i][j])/B_sat
            if np.shape(spot) == (0,) and B >= B_sat:
                ff[i,j]=1
#         if B < B_sat:
#              ff[i,j] =abs(data[i][j])/B_sat
#          ssif B >= B_sat:
#              ff[i,j]=1
            
                    
            x_rot = []   
         #   x_c = 0
        #    y_c = 0
            conv = np.pi/180.        
            x_rot=(j+13.28*(int(time)-start))%359
# x_rot = j
            x_pos = 180-x_rot
            y_pos = 90-i
            delta_lambda = abs(x_pos-x_c)
            distance = np.arccos((np.sin((y_c)*conv)*np.sin((y_pos)*conv)+np.cos((y_c)*conv) \
                              *np.cos((y_pos)*conv)*np.cos(delta_lambda*conv)))/conv
            vis = np.cos(distance*conv)
                #print(distance)
                #if distance <=90:
                 #   visibility.append(ff[i,j]*np.cos(distance*conv)*np.cos(y_pos*conv))
#            if vis <=1 and vis >= 0 and   np.shape(spot) !=(0,):
#                a.append(1-len(spot)*0.1*0.1)
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
            if distance <=90:
                visibility.append(ff[i,j]*np.cos(distance*conv)*np.cos(y_pos*conv))
                    
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
    g.write("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n" %(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,time,total,sum(visibility)/norm))
 #   f.write("%f \t \n" %(sum(a)))
#f.close()


g.close()

         


        
