import numpy as np
import apogee.tools.read as apread
from matplotlib import pyplot as plt
import pandas as pd
import csv
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
    if peak_loc == 0: #Give false data that is easily distinguished as such
        r = 1000 
        ccfCenter = 1000 #need to add a condition for peak too
        return r, ccfCenter
    endpoint = peak_loc+pm
    startpoint= peak_loc-pm
    Mirror = (x[peak_loc:endpoint])[::-1]
    sigmaA = np.sqrt(1.0 / (2.0 * len(Mirror)) * np.sum((x[startpoint:peak_loc] - Mirror)**2))
    r = np.max(x) / (np.sqrt(2.0) * sigmaA)
    return r, ccfCenter

#Calculate the bisector points of a given CCF and create the line bisector
def bisector(xccf,yccf):
    height = max(yccf) - min(yccf)
    slices = height/4.0
    new_max = max(yccf)+0.05
    bounds = np.linspace(min(yccf),new_max,5)
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

def x_range(x_bisector):
    xr = max(x_bisector) - min(x_bisector)
    xR = abs(xr)
    return xR

#Find the R ratios of the SB2s in log space
def r_ratio(r51,r151,r101):
        r1_ratio = r151/r101
        r2_ratio = r101/r51
        R1_ratio = math.log10(r1_ratio)
        R2_ratio = math.log10(r2_ratio)
        ratios = [round(R1_ratio,4),round(R2_ratio,4)]
        return ratios

#Turn lists into arrays and then log arrays
def arrays(x):
    x = np.array(x)
    new = x.astype(np.float)
    newer = np.log10(x)
    return newer

def routine(locationID,apogeeID,a,b):
    # Assign the interval over which to run the routine over
    # a = starting point in allStar file
    # b = ending point in allStar file
    a = a
    b = b
    #Create arrays for storing values and then to output into a .csv file
    SNR = []
    oldR151 = []
    oldR101 = []
    oldR51 = []
    xr = []
    R1 = []
    R2 = []
    Visit = []
    peak_value = []
    #create an array for storing HJDs
    HJDs = []
    loc = []
    apo= []
    #Run routine on DR14 to find R values, R ratios, x-ranges and HJDs in batches
    for j in range(len(locationIDs[a:b])):
            print(j)
            locationID = locationIDs[j]
            apogeeID = apogeeIDs[j]
            #File path to open .fits 
            my_file = Path('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits')
            try: 
                path = '/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/dr14/apogee/spectro/redux/r8/stars/apo25m/'+str(locationID)+'/'+'apStar-r8-'+str(apogeeID)+'.fits'
                data = fits.open(path)
                point = data[9]
                xccf = point.data[0][29] # Proper location of x values of the CCF
                CCF = point.data[0][27]
                HDU0 = fits.getheader(path,0)
                nvisits = HDU0['NVISITS']
                #edit this to more seamlessly handle single-visit sources (KRC had one single-visit source fail)
                for visit in range(0,nvisits):
                        if nvisits == 1:
                            ccf = CCF
                        else:
                            ccf = CCF[visit+2]
                        snr = HDU0['SNRVIS'+str(visit+1)]
                        #read in HJD identifier to more uniquely identify individual visits
                        HJD = HDU0['HJD'+str(visit+1)]
                        nonzeroes = np.count_nonzero(ccf) # This condition is meant to eliminate visits that are empty
                        if snr > 10:
                            if nonzeroes >= 1:
                                loc.append(locationID)
                                apo.append(apogeeID)
                                bs_pt = bisector(xccf, ccf)
                                x_bs = bs_pt[0]
                                y_bs = bs_pt[1]
                                xranges = x_range(x_bs)
                                xr.append(xranges)
                                SNR.append(snr) 
                                R151 = calcR(ccf,75)
                                partr151 = R151[1]
                                R101 = calcR(ccf,50)
                                R51 = calcR(ccf,25)
                                #Ensure 3 decimal places reported in .csv
                                oldR151.append(round(R151[0],3))
                                oldR101.append(round(R101[0],3))
                                oldR51.append(round(R51[0],3))
                                # Add to list of CCF peak Values
                                peak_value.append(partr151)
                                #Generate the Ratios 
                                Ratios = r_ratio(R51[0],R151[0],R101[0])
                                r1 = Ratios[0]
                                r2 = Ratios[1]
                                R1.append(r1)
                                R2.append(r2)
                                Visit.append(visit)
                                #store HJDs in an array
                                HJDs.append(HJD)

                            else:
                                pass

            except FileNotFoundError:
                pass
    return loc,apo,Visit,HDJs, peak_value ,oldR51, oldR101, oldR151, R1, R2, xr, SNR

## ---- ROUTINE BEGINS HERE ------ ##
'''
The Routine is stored into a function called 'routine' to split the computations done on the apStar files 
into two large chunks of the allStar file (0 to 150,000 and 150,000 to 266,094). This is done to prevent
the chance of a RunTime Error. 
How the 'routine' function works:
1. Send it locations ids and apogee ids. 'a' and 'b' are the start and end values you want the function
to go through in the allStar list.
2. Calls on the other functions defined before it to perform the calculations for R, R ratios, and x-range
3. Returns the arrays for :
  a) location ids, apogee ids, visit numbers, HDJs, peak values, R_51, R_101, R_151
  b) R151/R_101 (R1), R_101/R_51 (R2) [in log space]
  c) x-ranges, and signal-to-noise ratios (snr)
4. If there are sources it cannot find on the server (such as apStarC files), the routine will skip it
'''
# Read in allStar list for DR14 to get .fits of all stars in APOGEE
allStarDR14 = apread.allStar(rmcommissioning=False,main=False,ak=True,akvers='targ',adddist=False)
locationIDs = allStarDR14['LOCATION_ID']
apogeeIDs = allStarDR14['APOGEE_ID']
apogeeIDs = [s.decode('utf-8') for s in apogeeIDs]

#Run a chunk of the allStar file through the routine function (0 to 150,000)
chunk1 = routine(locationIDs, apogeeIDs,0,150000)
#Extract return arrays
c1_locationid = chunk1[0]
c1_apogeeid = chunk1[1]
c1_visit =chunk1[2]
c1_hdjs = chunk1[3]
c1_peakval =chunk1[4]
c1_R51 = chunk1[5]
c1_R101 = chunk1[6]
c1_R151 = chunk1[7]
c1_R1 = chunk1[8]
c1_R2 = chunk1[9]
c1_xr = chunk1[10]
print(c1_xr)
c1_snr = chunk1[11]

#Run another chunk of the allStar file through the routine function (150,000 to 266,094)
chunk2 = routine(locationIDs,apogeeIDs,150000,266094)
#Extract return arrays
c2_locationid = chunk2[0]
c2_apogeeid = chunk2[1]
c2_visit =chunk2[2]
c2_hdjs = chunk2[3]
c2_peakval =chunk2[4]
c2_R51 = chunk2[5]
c2_R101 = chunk2[6]
c2_R151 = chunk2[7]
c2_R1 = chunk2[8]
c2_R2 = chunk2[9]
c2_xr = chunk2[10]
c2_snr = chunk2[11]

#Combine the two chunks together to form a larger array. We do this by appending the second chunk 
# of arrays to the first chunk. This preserves the order of the elements in the array
for i in range(len(c2_locationid)):
    c1_locationid.append(c2_locationid[i])
    c1_apogeeid.append(c2_apogeeid[i])
    c1_visit.append(c2_visit)
    c1_hdjs.append(c2_hdjs)
    c1_peakval.append(c2_peakval)
    c1_R51.append(c2_R51)
    c1_R101.append(c2_R101)
    c1_R151.append(c2_R151)
    c1_R1.append(c2_R1)
    c1_R2.append(c2_R2)
    c1_xr.append(c2_xr)
    c1_snr.append(c2_snr)

#Replace all zero values to avoid non-float values (-inf, nan, etc) with an outlier # 11.
new_xr = [1000 if math.isinf(y) else y for y in c1_xr]

#Have DR14 sent to arrays function that also convert arrays into log space
newR51 = arrays(c1_R51)
newR101 = arrays(c1_R101)
newR151 = arrays(c1_R151)
new_Xrange = arrays(new_xr)
#peak_val = arrays(c1_peakval)

#Convert snr into an array for panda writing
SNRs = np.array(c1_snr)
SNRS = SNRs.astype(np.float)

#To not over-account for log space converstions, just convert the ratios into arrays.
# # The function r_ratios already converted the elements in these arrays into log space 
RatioR1 = np.array(c1_R1)
RatioR2 = np.array(c1_R2)
newR1 = RatioR1.astype(np.float)
newR2 = RatioR2.astype(np.float)
peak_val = c1_peakval.astype(np.float)

#Write the data frame to add columns of arrays to
cols = ['LocationID', 'ApogeeID']
df = pd.DataFrame(columns = cols)
# Assign the arrays to the column lavels of Location ID and Apogee ID
df['LocationID'] = c1_locationid
df['ApogeeID'] = c1_apogeeid
# Add more columns to the .csv file by defining a label df['column name'] and assigning it an array
# dataframe = df -> df['column name'] = array that holds desired data
df['SNR'] = SNRS
df['Visit'] = c1_visit
#add HJD info to output table and assign more labels associated arrays
df['HJD'] = c1_hdjs
df['log(R51)'] = newR51
df['log(R101)'] = newR101
df['log(R151)'] = newR151
df['log(xr)'] = new_Xrange
df['log(Ratio1)'] = newR1
df['log(Ratio2)'] = newR2
df['Peak_value'] = peak_val
#Write dataframe out to a .csv file
df.to_csv('DR14_Catalog_Revised.csv')
