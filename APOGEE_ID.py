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
        r = 1000
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

def xrange(x_bisector):
    #print(x_bisector)
    xr = max(x_bisector) - min(x_bisector)
    xR = abs(xr)
    return xR


#Find the R ratios of the SB2s
def r_ratio(r51,r151,r101):
        #print(r51, r101,r151)
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

#Turn lists into arrays and then log arrays
def arrays(x):
    x = np.array(x)
    new = x.astype(np.float)
    newer = np.log10(x)
    return newer

# Read in visually id SB2s and generate catalog for binaries
bins = pd.read_csv('KC_Binaries.csv',delimiter='\t')
locID = bins['Location_ID']
apoID = bins['Apogee_ID']
ids = bins['ID']

binApoID = []
binLocID = []
for i in range(len(ids)):
    if ids[i] != 0:
        binApoID.append(apoID[i])
        binLocID.append(locID[i])

#Run the routine on the identified binaries from IDl routine (the training set)
bin_SNR = []
R151s = []
R101s = []
R51s = []
bin_xr = []
binR1 = []
binR2 = []
Visit = []
loc = []
apo= []

for j in range(len(binLocID)):
        locationID = binLocID[j]
        apogeeID = binApoID[j]
        my_file = Path('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits')
        try: 
            path = '/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits'
            data = fits.open(path)
            point = data[9]
            xccf = point.data[0][29] # Proper location of x values of the CCF
            CCF = point.data[0][27]
            HDU0 = fits.getheader(path,0)
            nvisits = HDU0['NVISITS']
            for visit in range(0,nvisits):
                    ccf = CCF[visit+2]
                    snr = HDU0['SNRVIS'+str(visit+1)]
                    nonzeroes = np.count_nonzero(ccf) # This condition is meant to eliminate visits that are empty
                    if nonzeroes >= 1:
                        loc.append(locationID)
                        apo.append(apogeeID)
                        bs_pt = bisector(xccf, ccf)
                        x_range = xrange(bs_pt[0])
                        bin_xr.append(x_range)
                        bin_SNR.append(snr) 
                        R151 = calcR(ccf,75)
                        R101 = calcR(ccf,50)
                        R51 = calcR(ccf,25)
                        R151s.append(round(R151,3))
                        R101s.append(round(R101,3))
                        R51s.append(round(R51,3))
                        Ratios = r_ratio(R51,R151,R101)
                        r1 = Ratios[0]
                        r2 = Ratios[1]
                        binR1.append(r1)
                        binR2.append(r2)
                        Visit.append(visit)
                    else:
                        pass

        except FileNotFoundError:
            pass

#Binary arrays 
binR51 = arrays(R51s)
binR101 = arrays(R101s)
binR151 = arrays(R151s)
BinR1 = arrays(binR1)
BinR2 = arrays(binR2)
binxrange = arrays(bin_xr)

#Write out the results to a file via pandas
cols = ['LocationID', 'ApogeeID']
df = pd.DataFrame(columns = cols)
df['LocationID'] = loc
df['ApogeeID'] = apo
df['log(R51)'] = binR51
df['log(R101)'] = binR101
df['log(R151)'] = binR151
df['log(xr)'] = binxrange

df.to_csv('Binary_Stats.csv')
