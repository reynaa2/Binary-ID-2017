import numpy as np
import apogee.tools.read as apread
from matplotlib import pyplot as plt
import pandas as pd
import csv
from apogee.tools import bitmask
import math
from astropy.io import fits
import os.path
from pathlib import Path

#Calculate R-Values for given ranges
def calcR(x,pm):
    ccfCenter = max(x)
    primary = np.where(x == ccfCenter)
    peak_loc = primary[0][0]
    if peak_loc < pm: 
        pm = peak_loc
    if peak_loc > 401 - pm:
        pm = 401 - peak_loc
    if peak_loc == 0:
        r = 0
        return r 
    endpoint = peak_loc+pm
    startpoint= peak_loc-pm
    Mirror = (x[peak_loc:endpoint])[::-1]
    sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[startpoint:peak_loc] - Mirror)**2))
    r = np.max(x) / (np.sqrt(2.0) * sigmaA)
    return r
#Calculate the bisector points for a CCF (uses 4 points)
def bisector(xccf,yccf):
    height = max(yccf) - min(yccf)
    slices = height/4.0
    bounds = np.arange(min(yccf),height,slices)
    if len(bounds) != 0:
        z1 = (bounds[0] + bounds[1])/2.0
        z2 = (bounds[1] + bounds[2])/2.0
        z3 = (bounds[2] + bounds[3])/2.0
        z4 = (bounds[3] + bounds[4])/2.0
        y_bisector = np.array([z1,z2,z3,z4])

        x_bisector = []
        x0 = []
        x1 = []
        x2 = []
        x3 = []
        for i in range(len(yccf)):
            if yccf[i] <= bounds[4] and yccf[i] > bounds[3]:
                x0.append(xccf[i])
        x_0 = (np.mean(x0))
        x_bisector.append(x_0)

        i = 0
        for i in range(len(yccf)):
            if yccf[i] <= bounds[3] and yccf[i] >= bounds[2]:
                x1.append(xccf[i])
        x_1=(np.mean(x1))
        x_bisector.append(x_1)

        i = 0
        for i in range(len(yccf)):
            if yccf[i] <= bounds[2] and yccf[i] >= bounds[1]:
                x2.append(xccf[i])
        x_2=(np.mean(x2))
        x_bisector.append(x_2)

        i = 0
        for i in range(len(yccf)):
            if yccf[i] <= bounds[1] and yccf[i] >= bounds[0]:
                x3.append(xccf[i])
        x_3=(np.mean(x3))
        x_bisector.append(x_3)

        bisector_pts = np.vstack([x_bisector,y_bisector])
        #print(bisector_pts)
        return(bisector_pts)
    #else:
        #x_bisector = 0.0
        #y_bisector = 0.0
        #error = np.vstack([x_bisector,y_bisector])
        #return(error)

def xrange(x_bisector):
    #print(x_bisector)
    xr = max(x_bisector) - min(x_bisector)
    return xr

'''
This function was being used for finding the width of a CCF by locating the maxima, find the midpoint (y-direction),
and take the difference in lag space (x-direction) of the indicides that told where the sides of the CCF was.
The later steps still need to be implimented. 
'''
'''
def FindMaxima(xccf,CCF,peak_loc):
    thresh_down = 1
    if peak_loc > 7:
        peak_low = peak_loc-2
    else:
        peak_low = peak_loc+1
    while thresh_down > 0:
        if peak_low > 2 and (CCF[peak_low-1] - CCF[peak_low] < 0) and CCF[peak_low-1] > 0:
            peak_low = peak_loc-1
        else:
            thresh_down = 0
    thresh_up = 1
    if peak_loc < 394:
        peak_high = peak_loc+2
    else:
        peak_high = peak_loc-1
    while thresh_up > 0:
        if peak_high < 399 and (CCF[peak_high+1] - CCF[peak_high] < 0) and (CCF[peak_high+1] > 0):
            peak_high = peak_high+1
        else:
            thresh_up = 0
    n_lag = np.arange(0,402,1)
    walk = np.zeroes(n_lag)
    
    if peak_low < 0:
        peak_low = 0
    if else peak_high > n_lag-1:
        peak_high = n_lag-1
        
    walk[peak_low:peak_high] = 1
    ''';

#Calculate the R-ratios for the likely_binary function
def r_ratio(r51,r151,r101):
        r1_ratio = r151/r101
        r2_ratio = r101/r51
        R1_ratio = math.log10(r1_ratio)
        R2_ratio = math.log10(r2_ratio)
        ratios = [round(R1_ratio,3),round(R2_ratio,3)]
        return ratios

def idSB2s(R1_ratio, R2_ratio,r51,r151,r101,xr): # cuts to identify SB2s from Kevin's IDL Routine
    min_r51 = r51
    min_r101 = r101
    min_r151 = r151
    r1_ratio = R1_ratio
    r2_ratio = R2_ratio
    max_xr = xr
    
    likely_sb2s = np.where((math.log10(r1_ratio) > 0.06 and (math.log10(r1_ratio) < 0.13 and 
                            math.log10(min_r101) < 0.83)) or (math.log10(r2_ratio) > 0.05 and 
                            math.log10(r2_ratio) < 0.02 and math.log10(min_r51) < 0.83) and
                            math.log10(min_r51) > 0.25 and math.log10(min_r101) > 0.22 and
                            math.log10(peak_401) > -0.5 and math.log10(max_xr) < 2.3 and 
                            math.log10(max_xr) > 0.7
                          )
    return likely_sb2s

# Read in DR14 allStar File

allStarDR14 = apread.allStar(rmcommissioning=False,main=False,ak=True,akvers='targ',adddist=False)

locationIDs = allStarDR14['LOCATION_ID']
apogeeIDs = allStarDR14['APOGEE_ID']
apogeeIDs = [s.decode('utf-8') for s in apogeeIDs]

#Begin output for stats of all stars in DR14

with open('DR14_Stats_Catalog.csv','w') as output:
    column = ['Location_ID','Apogee_ID','x_range','R51','R101','R151','R401','R151/R101','R101/R51','Visit']
    writer = csv.DictWriter(output,delimiter='\t',fieldnames=column)
    writer.writeheader()
    for i in range(len(locationIDs)):
        locationID = locationIDs[i]
        apogeeID = apogeeIDs[i]
        print(i,locationID,apogeeID)
        my_file = Path('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits')
        try: 
            path = '/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits'
    
        except:
            path = '/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStarC-r8-'+str(apogeeID)+'.fits'
        #print(path)
        data = fits.open(path)
        point = data[9]
        xccf = point.data[0][32]
        CCF = point.data[0][27]
        HDU0 = fits.getheader(path,0)
        nvisits = HDU0['NVISITS']
        for visit in range(0,nvisits):
            if nvisits != 1:
                ccf = CCF[visit+2]
                nonzeroes = np.count_nonzero(ccf) # This condition is meant to eliminate visits that are empty
                if nonzeroes >= 1:
                    bs_pt = bisector(xccf, ccf)
                    x_range = xrange(bs_pt[0])
                    R151 = calcR(ccf,75)
                    R101 = calcR(ccf,50)
                    R51 = calcR(ccf,25)
                    Ratios = r_ratio(R51,R151,R101)
                    r1 = Ratios[0]
                    r2 = Ratios[1]
                    writer.writerow({'Location_ID':locationID,'Apogee_ID':apogeeID,'x_range':round(x_range,3),
                                         'R51':round(R51,3),'R101':round(R101,3),'R151':round(R151,3),'R151/R101':r1,'R101/R51':r2,'Visit':visit})


# Read in the file and find how many are SB2s. 
stats = pd.read_csv('DR14_Stats_Catalog.csv',delimiter='/t')
FieldID = stats['Location_ID']
TwoMassID = stats['Apogee_ID']
Star_Visit = stats['Visit']
logR1 = stats['R151/R101']
logR2 = stats['R101/R51']
R51s = stats['R51']
R101s = stats['R101']
R151s = stats['R151']
xr_value = stats['x_range']

with open('DR14_SB2_Catalog.csv','w') as files:
    column = ['Location_ID','Apogee_ID']
    writer = csv.DictWriter(files,delimiter='\t',fieldnames=column)
    writer.writeheader()
    i = 0
    r1 = []
    r2 = []
    r51 = []
    r101 = []
    r151 = []
    xranges = []
    SB2 = []
    for i in range(len(FieldID)):
        if Visit[i] == Visit[i]:
            r1.append(logR1[i])
            r2.append(logR2[i])
            r51.append(R51s[i])
            r101.append(R101s[i])
            r151.append(R151s[i])
            xranges.append(xr_value[i])
        likely_binary = likely_SB2s(r1, r2,r51,r151,r101,xranges)
        index_locID = FieldID[likely_binary]
        index_apoID = TwoMassID[likely_binary]
        SB2.append([index_locID,index_apoID])
        writer.writerow({'Location_ID': SB2[0] ,'Apogee_ID':SB2[1]})


