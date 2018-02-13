import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv
from astropy.io import fits
import os.path
from pathlib import Path

'''
Open DR14_Stats_Catalog.csv. Over all visits for a star we will need to:
1. Find the minimum R values: R 51, R 101, R151
2. Find the maximum x-range values
3. Find the minimum R ratios: R 151/R 101 and R 101/R 51
 '''
  
def matches(xx,y,r51,r101,r151,xr,R1,R2):
    match51 = []
    match101 = []
    match151 = []
    matchxr = []
    matchR1 = []
    matchR2 = []
    matchlocid = []
    matchapoid = []
    #Command for finding unique values [in this case, of APOGEE IDs]
    uniqueID = np.unique(y)

    ## Do operation [of what kind?] on each star (as identified by unique APOGEE IDs)
    for i in range(len(uniqueID)):
        #find indices that identify visits to star w/ uniqueID[i]
        newID = np.where(y == uniqueID[i])
        #find first object that satisfies this criteria?  Or multiple objects?
        newer = newID[0]
        #now save r51 values associated with one or more objects?
        small = r51[newer]
        #find minimum of r51 values selected above
        little = min(small)
        #round the R51 value to four digits?
        x = round(little,4)
        #save the min R51 value, and the APOGEE ID and Location ID associated for this object
        match51.append(x)
        matchapoid.append(y[newer])
        matchlocid.append(xx[newer])

    for j in range(len(uniqueID)):
        new = np.where(y == uniqueID[j])
        element = new[0]
        R101 = r101[element]
        lil = min(R101)
        x2 = round(lil,4)
        match101.append(x2)

    for k in range(len(uniqueID)):
        New = np.where(y == uniqueID[k])
        part = New[0]
        R151 = r151[part]
        bit = min(R151)
        x3 = round(bit,4)
        match151.append(x3)

    for p in range(len(uniqueID)):
        mew = np.where(y == uniqueID[p])
        elmnt = mew[0]
        XR = xr[elmnt]
        bit2 = max(XR)
        x4 = round(bit2,4)
        matchxr.append(x4)

    for o in range(len(uniqueID)):
        nID = np.where(y == uniqueID[o])
        news = nID[0]
        s = R1[news]
        l = min(s)
        x5 = round(l,4)
        matchR1.append(x5)

    for q in range(len(uniqueID)):
        mID = np.where(y == uniqueID[q])
        ns = mID[0]
        pico = R2[ns]
        tiny = min(pico)
        x6 = round(tiny,4)
        matchR2.append(x6)
    return match51, match101, match151, matchxr,matchR1,matchR2, matchlocid,matchapoid

#Turn lists into arrays 
def arrays(x):
    x = np.array(x)
    new = x.astype(np.float)
    return new

# From Cuts.py, the parameters for an SB2 to be identified have been quantified. This function holds those values
def idSB2(R51,R101,R151,xr,ratio1,ratio2):
    likely_sb2s = np.where((0.7 < R101 < 1.2 and 0.8 < ratio1 < 0.95) or 
                             (0.6 < R51 < 0.70 and 0.80 < R151 < 0.95) or 
                             (0.6 < xr < 1.1 and 0.10 < ratio2 < 0.15) or
                             (0.6 < xr < 1.1 and 0.08 < ratio1 < 0.090) or 
                             (0.6))
    return likely_sb2s

## THIS IS WHERE THE SCRIPT ACTUALLY STARTS

#Read in csv Dr14_Stats_Catalog and require targets to pass a SNR > 10
dr14 = pd.read_csv('DR14_Stats_Catalog2.csv')
locID = dr14['LocationID']
apoID = dr14['ApogeeID']
visit = dr14['Visit']
snr = dr14['SNR']
R51 = dr14['log(R51)']
R101 = dr14['log(R101)']
R151 = dr14['log(R151)']
ratio1 = dr14['log(Ratio1)']
ratio2 = dr14['log(Ratio2)']
xranges = dr14['log(xr)']

# Create new lists to be turned into arrays if they pass the SNR > 10 resuirement
locationid = []
apogeeid = []
r51 = []
r101 = []
r151 = []
R1 = []
R2 = []
vis = []
xr = []

#creating arrays that only hold objects that pass SNR > 10 requirement 
for i in range(len(locID)):
    if snr[i] > 10:
        locationid.append(locID[i])
        apogeeid.append(apoID[i])
        r51.append(R51[i])
        r101.append(R101[i])
        r151.append(R151[i])
        R1.append(ratio1[i])
        R2.append(ratio2[i])
        vis.append(visit[i])
        xr.append(xranges[i])

#Turn these lists above to arrays by passing them to the 'arrays' function (which converts to type float?  Maybe not appropriate for LocID and ApoID?)
LocID = arrays(locationid)
ApoID = arrays(apogeeid)
r_51 = arrays(r51)
r_101 = arrays(r101)
r_151 = arrays(r151)
ratio_1 = arrays(R1)
ratio_2 = arrays(R2)
Visit = arrays(vis)
x_range = arrays(xr)

#Find the minimum R values and the maximum x-range values of a given star. 
# Pass each array given above into the 'matches' definition.  
m = matches(LocID,ApoID,r_51,r_101,r_151,x_range)
min51 = arrays(m[0])
min101 =arrays(m[1])
min151 = arrays(m[2])
maxXR = arrays(m[3])
minratio1 = arrays(m[4])
minratio2 = arrays(m[5])
matchedapoid = arrays(m[6])
matchedlocid = arrays(m[7])

''' This is more of a thought but, here is a wider range for identifying sb2s
def idsb2s(r51,r101,r151,xr,ratio1,ratio2):
    likely_sb2s = np.where( (0.4 < r101 < 1.4 and 0.02 < ratio1 < 0.95) or
                            (0.25 < r151 < 1.50 and 0.5 < r51 < 1.50) or
                            (0.25 < r51 < 1.50 and 0.05 < ratio2 < 0.15) or
                            (0.05 < xr < 1.75 and 0.05 < ratio1 < 0.09) or
                            (0.05 < xr < 1.5 and 0.05 < ratio2 < 0.15)
                            )
'''


# From the provided indicies for binaries (generated in the idSB2 function) we can find the assocated location ID, 
# Apogee ID, and other parameters of the star. This will be used to output into a .csv file

likely_sb2 = likely_sb2s(min51,min101,min151,maxXR,minratio1,minratio2)

bin_apoID = []
bin_locID = []
bin_r51 = []
bin_r101 = []
bin_r151 = []
bin_r1 = []
bin_r2 = []
bin_xr = []
 
for i in range(len(likely_sb2)):
     apo_elem = ApoID[i]
     loc_elem = LocID[i]
     r51_elem = min51[i]
     r101_elem = min101[i]
     r151_elem = min151[i]
     rat1_elem = minratio1[i]
     rat2_elem = minratio2[i]
     xr_elem = maxXR[i]
     bin_apoID.append(apo_elem)
     bin_locID.append(loc_elem)
     bin_r51.append(r51_elem)
     bin_r101.append(r101_elem)
     bin_r151.append(r151_elem)
     bin_r1.append(rat1_elem)
     bin_r2.append(rat2_elem)
     bin_xr.append(xr_elem)

#Construct catalog of identified binaries with defined arrays above
cols = ['LocationID', 'ApogeeID']
df = pd.DataFrame(columns=cols)
df['LocationID'] = bin_locID
df['ApogeeID'] = bin_apoID
df['log(R51)'] = bin_r51
df['log(R101'] = bin_r101
df['log(R151)'] = bin_r151
df['log(R151/R101)'] = bin_r1
df['log(R101/R51)'] = bin_r2
df['log(x-range)'] = bin_xr

#Write to a csv file
df.to_csv('DR14_Binary_Candidates.csv', header=True, index=False) 
     
     
