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

#Calculate the bisector points for a CCF (uses 4 points)
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

def xrange(x_bisector):
    xr = max(x_bisector) - min(x_bisector)
    xR = abs(xr)
    return xR


#Find the R ratios of the SB2s
def r_ratio(r51,r151,r101):
        r1_ratio = r151/r101
        r2_ratio = r101/r51
        R1_ratio = math.log10(r1_ratio)
        R2_ratio = math.log10(r2_ratio)
        ratios = [round(R1_ratio,4),round(R2_ratio,4)]
        return ratios

# def idSB2s(R1_ratio, R2_ratio,r51,r151,r101,xr): # cuts to identify SB2s from Kevin's IDL Routine
#     min_r51 = r51
#     min_r101 = r101
#     min_r151 = r151
#     r1_ratio = R1_ratio
#     r2_ratio = R2_ratio
#     max_xr = xr
    
#     likely_sb2s = np.where((math.log10(r1_ratio) > 0.06 and (math.log10(r1_ratio) < 0.13 and 
#                             math.log10(min_r101) < 0.83)) or (math.log10(r2_ratio) > 0.05 and 
#                             math.log10(r2_ratio) < 0.02 and math.log10(min_r51) < 0.83) and
#                             math.log10(min_r51) > 0.25 and math.log10(min_r101) > 0.22 and
#                             math.log10(peak_401) > -0.5 and math.log10(max_xr) < 2.3 and 
#                             math.log10(max_xr) > 0.7
#                           )
#     return likely_sb2s

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

#Placed main routine in a function to compact space.
#def main_routine(apogeeID,locationID):

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

# Read in allStar list for DR14 to get .fits of all stars in APOGEE
allStarDR14 = apread.allStar(rmcommissioning=False,main=False,ak=True,akvers='targ',adddist=False)
locationIDs = allStarDR14['LOCATION_ID']
apogeeIDs = allStarDR14['APOGEE_ID']
apogeeIDs = [s.decode('utf-8') for s in apogeeIDs]

#Run routine on DR14 to find R values, R ratios, x-ranges and HJDs
#for j in range(len(binLocID[1:100])):
for j in range(len(locationIDs)):
        print(j)
        #locationID = binLocID[j]
        #apogeeID = binApoID[j]
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
                    if nonzeroes >= 1:
                        loc.append(locationID)
                        apo.append(apogeeID)
                        bs_pt = bisector(xccf, ccf)
                        #split up 2D array to get x and y coord. for bisector pts
                        x_bs = bs_pt[0]
                        y_bs = bs_pt[1]
                        #x_range = xrange(bs_pt[0])
                        x_range = xrange(x_bs)
                        #Plotting CCFs
                        #plt.plot(xccf,ccf, label='Visit: '+str(visit))
                        #plt.plot(x_bs,y_bs,'o')
                        #plt.legend(loc='upper right')
                        #plt.show()
                        xr.append(x_range)
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


#Find and replace all nan values with 9 which will be prominent in log space
#x_ranges = [9 if math.isnan(x) else x for x in xr]
# Replace all -inf values with an outlier # which we will assign to be 1000
for Xranges, i in enumerate(x_ranges):
    if i == 0:
        x_ranges[Xranges] = 1000
#Have DR14 sent to arrays function that also convert arrays into log space
newR51 = arrays(oldR51)
newR101 = arrays(oldR101)
newR151 = arrays(oldR151)
new_Xrange = arrays(x_ranges)
peak_val = arrays(peak_value)

#Find and replace all nan values with 9 which will be prominent in log space
#x_ranges = [9 if math.isnan(x) else x for x in Xrange]
#new_Xranges = [9 if math.isinf(y) else y for y in x_ranges]

#Convert snr into an array for panda writing
SNRs = np.array(SNR)
SNRS = SNRs.astype(np.float)

#To not over account for log space converstions, just convert the ratios into arrays 
RatioR1 = np.array(R1)
RatioR2 = np.array(R2)
newR1 = RatioR1.astype(np.float)
newR2 = RatioR2.astype(np.float)

#Write out the results to a file via pandas
cols = ['LocationID', 'ApogeeID']
df = pd.DataFrame(columns = cols)
df['LocationID'] = loc
df['ApogeeID'] = apo
df['SNR'] = SNRS
df['Visit'] = Visit

#add HJD info to output table
df['HJD'] = HJDs
df['log(R51)'] = newR51
df['log(R101)'] = newR101
df['log(R151)'] = newR151
df['log(xr)'] = new_Xrange
df['log(Ratio1)'] = newR1
df['log(Ratio2)'] = newR2
df['Peak_value'] = peak_val

#df.to_csv('TrainingSet_Binary_Stats.csv')
df.to_csv('DR14_Catalog_Stats.csv')

     
