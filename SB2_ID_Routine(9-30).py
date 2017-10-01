import numpy as np
import apogee.tools.read as apread
from matplotlib import pyplot as plt
import pandas as pd
import peakutils
import csv
from apogee.tools import bitmask
import math

allStarDR14 = apread.allStar(rmcommissioning=True,main=False,ak=True,akvers='targ',adddist=False)
locationIDs = allStarDR14['LOCATION_ID']
apogeeIDs = allStarDR14['APOGEE_ID']
apogeeIDs = [s.decode('utf-8') for s in apogeeIDs]

#Calculate R-Values for given ranges
def calcR401(x, pos1=0, pos2=401):
        # Calculates the value of R with the given array x
        # Returns:  The value of R for whole CCF
        # Assupmtion: the center peak lies in CCF lag space 201
        ccfCenter = 201
        pos1+= 1
        Mirror = (x[ccfCenter:pos2])[::-1]
        sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[pos1:ccfCenter] - Mirror)**2))
        r401 = np.max(x) / (np.sqrt(2.0) * sigmaA)
        #print('R401 = '+str(r401))
        return r401

def calcR151(x, pos1=125, pos2=276):
        ccfCenter = 201
        pos1+= 1
        Mirror = (x[ccfCenter:pos2])[::-1]
        sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[pos1:ccfCenter] - Mirror)**2))
        r151 = np.max(x) / (np.sqrt(2.0) * sigmaA)
        #print('R151 = '+str(r151))
        return r151

def calcR101(x, pos1=150, pos2=251):
        ccfCenter = 201
        pos1+= 1
        Mirror = (x[ccfCenter:pos2])[::-1]
        sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[pos1:ccfCenter] - Mirror)**2))
        r101 = np.max(x) / (np.sqrt(2.0) * sigmaA)
        #print('R101 = '+str(r101))
        return r101
    
def calcR51(x, pos1=175, pos2=226):
        ccfCenter = 201
        pos1+= 1
        Mirror = (x[ccfCenter:pos2])[::-1]
        sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[pos1:ccfCenter] - Mirror)**2))
        r51 = np.max(x) / (np.sqrt(2.0) * sigmaA)
        #print('R51 = '+str(r51))
        return r51

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
        return(error)

def xrange(x_bisector):
    print(x_bisector)
    xr = max(x_bisector) - min(x_bisector)
    return xr

def array_builder(locationID,apogeeID,nvisits):
    header = apread.apStar(locationID,apogeeID,ext=0,header=True)
    data = apread.apStar(locationID,apogeeID,ext=9,header=False)
    nvisits=header[1]['NVISITS']
    y = []
    for visit in range(nvisits):
        if nvisits != 1:
            CCF = data['CCF'][0][2+visit]
            a = calcR151(CCF)
            y.append(a)
    #print(y)
    #print(min(y))
    return y

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

def idSB2s(CCF,xr,r51,r151,r101): # cuts to identify SB2s
    min_r51 = min(r51)
    min_r101 = min(r101)
    min_r151 = min(r151)
    peak_401 = max(CCF) - min(CCF)
    max_xr = max(xr)
    
    r1_ratio = min_r151/min_r101
    r2_ratio = min_r101/min_r51
    
    likely_sb2s = np.where((math.log10(r1_ratio) > 0.06 and (math.log10(r1_ratio) < 0.13 and 
                            math.log10(min_r101) < 0.83)) or (math.log10(r2_ratio) > 0.05 and 
                            math.log10(r2_ratio) < 0.02 and math.log10(min_r51) < 0.83) and
                            math.log10(min_r51) > 0.25 and math.log10(min_r101) > 0.22 and
                            math.log10(peak_401) > -0.5 and math.log10(max_xr) < 2.3 and 
                            math.log10(max_xr) > 0.7
                          )
    return likely_sb2s

# Alternative way to find the maxima. This utilizes peakutils python package.
'''
def detect_peaks(CCF):
    ccf_max = max(CCF)
    ccf_min = min(CCF)
    height = ccf_max - ccf_min
    half_height = height/2.0
    #print(half_height)
    
    #Find the indicies of the peaks
    indices = peakutils.indexes(CCF,thres=(ccf_max - ccf_min)/2.0,min_dist=0.1)
    #print(indicies)
    
    #Make sure to find the y-values of the peaks
    y_indicies = CCF[indicies]
    #print(y_indicies)
    
    #Make sure that the peaks are the global maxima
    y_indx = []
    x_indx = []
    for i in range(len(indicies)):
        if y_indicies[i] >= half_height:
            y_indx.append(y_indicies[i])
            x_indx.append(indicies[i])
    return x_indx, y_indx

def index(CCF,indicies,half_height):
    y_indicies=CCF[indicies]
    y_indx = []
    x_indx = []
    for j in range(len(indicies)):
        if y_indicies[j] >= half_height:
            y_indx.append(y_indicies[j])
            x_indx.append(indicies[j])
    return x_indx,  y_indx
'''

# ----- Main Routine that calls in all of DR14 -----

# Read in the training set SB2s identified from IDL. The main focus is to make sure the parameter functions are working correctly for SB2s (R, max x-range,peak finder, etc). 
k = pd.read_csv('KC_Binaries.csv',delimiter='\t')
vid = k['ID']
locID = k['Location_ID']
apoID = k['Apogee_ID']

locID_bins = []
apoID_bins = []
for i in range(len(vid)):
    if vid[i] != 0:
        locID_bins.append(locID[i])
        apoID_bins.append(apoID[i])

        
xccf = np.arange(0,402,1)
with open('Binary_Maxes.csv','w') as output:
    names = ['Location_ID','Apogee_ID','Maxes']
    writer = csv.DictWriter(output,delimiter='\t',fieldnames=names)
    writer.writeheader()
    try:
        for i in range(len(apoID_bins)):
            location_ID = locID_bins[i]
            apogee_ID = apoID_bins[i]
            header = apread.apStar(location_ID,apogee_ID,ext=0,header=True)
            data = apread.apStar(location_ID,apogee_ID,ext=9,header=False)
            nvisits=header[1]['NVISITS']

            for visit in range(nvisits):
                if visit != 1:
                    CCF = data['CCF'][0][2+visit]
                    #half_height = (max(CCF) - min(CCF))/2.0
                    #indicies = peakutils.indexes(CCF,thres=(max(CCF) - min(CCF))/2.0,min_dist=0.5)
                    #y_indicies = CCF[indicies]
                    #indx = index(CCF,indicies,half_height)
                    vs_indx = detect_peaks(CCF)
                    print(indx)
                    print(vs_indx)
                else:
                    CCF = data['CCF'][0]

                #writer.writerow({'Location_ID':location_ID,'Apogee_ID':apogee_ID,'Maxes':x_indx})
    except FileNotFoundError:
            pass

'''
#This is the routine for reading in all of DR14. The issue here is that the apstarC file isn't reading correctly. 

from astropy.io import fits

for i in range(len(locationIDs)):
        apogeeID = apogeeIDs[i]
        locationID = locationIDs[i]
    
        path = '/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits'
        
        if len(path) == 34: #If Statement designed to read in apStarC files
            if locationID != 1:
                cdata = fits.open(path)
                cpt = Cdata[9]
                cCCFs = cpt.data[0][28]
                cHDU0 = fits.getheader(path,0)
                cvisits = cHDU0['NVISITS']

                if cnvisits >= 3:
                    for visit in range(0,cvisits):
                        snr = cHDU0['SNVIS'+str(visit+1)]
                        if cvisits !=1:
                            ccf = cCCFs[visit+2]
                            r401 = calcR401(ccf)
                            r151 = calcR151(ccf)
                            r101 = calcR101(ccf)
                            r51 = calcR51(ccf)
                            bs_pts = bisector(x,ccf,cvisits)
                            mxr = xrange(cvisits, bs_pts[0])
                        else:
                            ccf = cCCFs
        else:
            if locationID !=1:
                Data = fits.open(path)
                stpt = Data[9]
                CCFs = stpt.data[0][28]
                HDU0 = fits.getheader(path,0)
                nvisits = HDU0['NVISITS']

                if nvisits >= 3:
                    for visit in range(0,nvisits):
                        SNR = HDU0['SNRVIS'+str(visit+1)]
                        if nvisits!=1:
                            CCF = CCFs[visit+2]
                            r401 = calcR401(CCF)
                            r151 = calcR151(CCF)
                            r101 = calcR101(CCF)
                            r51 = calcR51(CCF)
                            bs_pts = bisector(xccf,CCF,nvisits)
                            mxr = xrange(nvisits, bs_pts[0])
                        else:
                            CCF = CCFs
'''
