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
  
'''
NOTE: There is an unknown KeyError:2 and I am not sure why... try to investigate into this...
'''

def matches(xx,y,r51,r101,r151,xr,R1,R2,peak_val):
   # Generate empty lists to hold values of matched APOGEE IDs
    match51 = []
    match101 = []
    match151 = []
    matchxr = []
    matchR1 = []
    matchR2 = []
    matchlocid = []
    matchapoid = []
    matchpeak = []
    #Command for finding unique values (in this case, APOGEE IDs)
    uniqueID = np.unique(y)

    # Loop through all the unique apogeeIDs to find the smallest R values, largest x-range, and the apogee and location IDs
    # that correspond to these values
    for i in range(len(uniqueID)):
        # find indicies that identify visits to star w/ uniqueID[i]
        newID = np.where(y == uniqueID[i])
        # find objects that satisfies criteria of having matched apogee IDs
        newer = newID[0]
        # Plug in found indicies (from newer) and call R 51 values and the apogee and location IDs that correspond to those indicies
        small = r51[newer] 
        y2 = y[newer]
        x2 = xx[newer]
        peak = peak_val[newer]
        # Take smallest R 51 value of the list found above
        little = min(small)
        # Report back the apogee and location ID that correspond with the smallest R 51 value
        loc_of_small = np.where(small == little)
        #Extract value from tuple array
        small_loc = loc_of_small[0][0]
        # Round R 51 smallest value to 4 digits
        x = round(little,4)
        # Append the smallest R 51 values to empty list identified above as match51
        match51.append(x)
        # Save associated apogee ID and location ID for this object
        matchapoid.append(y2[small_loc])
        matchlocid.append(x2[small_loc])
        matchpeak.append(peak[small_loc])

    # Process above is repeated for other empty lists
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

    return match51, match101, match151, matchxr,matchR1,matchR2, matchlocid,matchapoid,matchpeak

#Turn lists into arrays 
def arrays(x):
    x = np.array(x)
    new = x.astype(np.float)
    return new

# From Cuts.py, the parameters for an SB2 to be identified have been quantified. This function holds those values
# Also, this function will return the indicies of where the targets meet these requirements.
def idSB2(R51,R101,R151,xr,ratio1,ratio2):
    likely_sb2s = np.where((0.7 < R101 < 1.2 and 0.8 < ratio1 < 0.95) or 
                             (0.6 < R51 < 0.70 and 0.80 < R151 < 0.95) or 
                             (0.6 < xr < 1.1 and 0.10 < ratio2 < 0.15) or
                             (0.6 < xr < 1.1 and 0.08 < ratio1 < 0.090))# or 
                             #(0.6 <)) #FIX ME
                             # Incorporate the peak_loc requirement as well. Kevin has this as peak_401 > -0.5
    return likely_sb2s


####### THIS IS WHERE THE ROUTINE BEGINS! THIS IS WHERE THE TOTAL .CSV FILE OF DR14 AND STATS ARE READ IN ########

#Read in csv Dr14_Stats_Catalog and require targets to pass a SNR > 10
dr14 = pd.read_csv('DR14_Stats_Catalog2a.csv')
locID = np.array(dr14['LocationID'])
apoID = np.array(dr14['ApogeeID'])
visit = dr14['Visit']
snr = dr14['SNR']
R51 = dr14['log(R51)']
R101 = dr14['log(R101)']
R151 = dr14['log(R151)']
ratio1 = dr14['log(Ratio1)']
ratio2 = dr14['log(Ratio2)']
xranges = dr14['log(xr)'] 
peak_location = dr14['Peak_value']

### TEST FOR SMALL SAMPLE OF DR14 ONLY ####

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
peaks = []

# Create new lists to only hold objects if they pass the SNR > 10 resuirement
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
        peaks.append(peak_location[i])

LocID = np.array(locationid)
ApoID = np.array(apogeeid)
peak_locs = np.array(peaks)
R_51 = arrays(r51)
R_101 = arrays(r101)
R_151 = arrays(r151)
ratio_1 = arrays(R1)
ratio_2 = arrays(R2)
x_ranges = arrays(xr)

m = matches(LocID,ApoID,R_51,R_101,R_151,x_ranges,ratio_1,ratio_2,peak_locs)
min51 = arrays(m[0])
min101 =arrays(m[1])
min151 = arrays(m[2])
maxXR = arrays(m[3])
minratio1 = arrays(m[4])
minratio2 = arrays(m[5])
matchedapoid = np.array(m[6])
matchedlocid = np.array(m[7])
matched_peaks = np.array(m[8])

## FOR TESTING, USE BINARY_STATS.CSV ####
def binary_test(x,y):
    data = pd.read_csv('Binary_Stats.csv')
    locID = data['LocationID']
    apoID = data['ApogeeID']
    #visit = data['Visit']
    #snr = data['SNR']
    R51 = data['log(R51)']
    R101 = data['log(R101)']
    R151 = data['log(R151)']
    ratio1 = data['log(Ratio1)']
    ratio2 = data['log(Ratio2)']
    xranges = data['log(xr)'] 

    ### TEST FOR BINARIES ONLY ####
    r_51 = arrays(R51)
    r_101 = arrays(R101)
    r_151 = arrays(R151)
    ratio_1 = arrays(ratio1)
    ratio_2 = arrays(ratio2)
    x_range = arrays(xranges)

    m = matches(locID,apoID,r_51,r_101,r_151,x_range,ratio_1,ratio_2)
    min51 = arrays(m[0])
    min101 =arrays(m[1])
    min151 = arrays(m[2])
    maxXR = arrays(m[3])
    minratio1 = arrays(m[4])
    minratio2 = arrays(m[5])
    matchedapoid = np.array(m[6])
    matchedlocid = np.array(m[7])
    return

#This is more of a thought but, here is a wider range for identifying sb2s
def idsb2s(r51,r101,r151,xr,ratio1,ratio2):
    likely_sb2s = np.where( (0.4 < r101 < 1.4 and 0.02 < ratio1 < 0.95) or
                            (0.25 < r151 < 1.50 and 0.5 < r51 < 1.50) or
                            (0.25 < r51 < 1.50 and 0.05 < ratio2 < 0.15) or
                            (0.05 < xr < 1.75 and 0.05 < ratio1 < 0.09) or
                            (0.05 < xr < 1.5 and 0.05 < ratio2 < 0.15)
                            )


# Transform np.where definition to a series of if statements to avoid ambiguity
def idSB2(locid,apoid,r51,r101,r151,xr,ratio1,ratio2,peak_val):
    candidates_apoid = []
    candidates_locid = []
    
    for i in range(len(r51)):
        if (0.4<r51[i]<1.4) and (0.02<ratio1[i]<0.95):
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])

        if (0.25<r151[i]<1.5) and (0.5<r51[i]<0.15):
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])

        if (0.25<r51[i]<1.5) and (0.05<ratio2[i]<0.09):
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])

        if (0.05<xr[i]<1.75) and (0.05<ratio1[i]<0.09):
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])

        if (0.05<xr[i]<1.4) and (0.05<ratio2[i]<0.95):
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])
        
        if peak_val[i] > -0.5:
            candidates_locid.append(locid[i])
            candidates_apoid.append(apoid[i])

    return candidates_locid,candidates_apoid

likely_binaries = idSB2(matchedlocid,matchedapoid,min51,min101,min151,maxXR,minratio1,minratio2,matched_peaks)      
       


# From the provided indicies for binaries (generated in the idSB2 function) we can find the assocated location ID, 
# Apogee ID, and other parameters of the star. This will be used to output into a .csv file

#likely_sb2 = idsb2s(min51,min101,min151,maxXR,minratio1,minratio2)

'''
new_apoID = []
new_locID = []
new_r51 = []
new_r101 = []
new_r151 = []
new_r1 = []
new_r2 = []
new_xr = []
 
for i in range(len(likely_sb2)):
     apo_elem = matchedapoid[i]
     loc_elem = matchedlocid[i]
     r51_elem = min51[i]
     r101_elem = min101[i]
     r151_elem = min151[i]
     rat1_elem = minratio1[i]
     rat2_elem = minratio2[i]
     xr_elem = maxXR[i]
     new_apoID.append(apo_elem)
     new_locID.append(loc_elem)
     new_r51.append(r51_elem)
     new_r101.append(r101_elem)
     new_r151.append(r151_elem)
     new_r1.append(rat1_elem)
     new_r2.append(rat2_elem)
     new_xr.append(xr_elem)

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
df.to_csv('DR14_Binary_Candidates.csv', header=True, index=False) '''
     
     