import numpy as np
import apogee.tools.read as apread
from matplotlib import pyplot as plt
import pandas as pd
import math
from astropy.io import fits
import os.path
from pathlib import Path

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
        return(bisector_pts)

def XRANGE(x_bisector):
    xr = max(x_bisector) - min(x_bisector)
    xR = abs(xr)
    return xR

# Generate Function to make CCFs with their respective line bisector
def make_ccfs(apogeeid,locationid, xccf,yccf, x_bs,y_bs, xrange_value, y_bisector, n ,ID):
    # Round x-range value to 3 digits after decimal
    xrange_val = round(xrange_value,3)
    ID = round(ID,3)
    # Plotting CCFs
    if n != 0:
        plt.grid(alpha=0.70)
        plt.plot(xccf,yccf, color = 'black',label='Visit: '+str(visit))
        plt.annotate('binary % = '+str(ID),xy=(-200,max(yccf)))
        plt.xlabel('CCF Lag',fontsize=17)
        plt.ylabel('CCF Units',fontsize=17)
        plt.legend(loc='upper left')
        # Save to the folder on Server under Bnary IS
        plt.savefig('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/BinaryID/MachineLearning_IDStars_June2018/'+str(locationid)+'_'+str(apogeeid)+'_'+str(visit)+'.pdf', dpi=650,bbox_to_inches='tight')
        plt.close('all')
    if n == 0:
        plt.plot(xccf,yccf, color = 'black',label='Visit: '+str(visit))
        plt.xlabel('CCF Lag',fontsize=17)
        plt.ylabel('CCF Units',fontsize=17)
        plt.scatter(x_bs, y_bs, color='red', s=25)
        plt.plot(x_bs,y_bs, color='red',label = '$\Delta$ X: '+str(xrange_val))
        # Create horizontal lines where the bisector y lines are to make distinct regions
        plt.axhline(y=y_bisector[0],color='gray', linestyle='--') # plot horizontal line at first y bisector
        plt.axhline(y=y_bisector[1],color='gray', linestyle='--')
        plt.axhline(y=y_bisector[2],color='gray', linestyle='--')
        plt.axhline(y=y_bisector[3],color='gray', linestyle='--')
        plt.legend(loc='upper left')
        # Save to the folder in vc_workspace/JMR Stars CCFs
        plt.savefig('JMR Star CCFs/Training Set/'+str(locationid)+'_'+str(apogeeid)+'_'+str(visit)+'.png', dpi=650,bbox_to_inches='tight')
        plt.close('all')
    return 

# Make CCFs that overplot one another for the same apogeeID
def ccf_overplot(xccf,yccf,locationid,apogeeid,xrange_value,ID):
    if apogeeid == apogeeid:
        plt.grid(alpha=0.70)
        plt.plot(xccf,yccf+visit, color = 'black',label='Visit: '+str(visit))
        # plt.annotate('binary % = '+str(ID),xy=(-200,max(yccf)))
        plt.xlabel('CCF Lag',fontsize=17)
        plt.ylabel('CCF Units',fontsize=17)
        plt.legend(loc='upper right')
        plt.title(apogeeid)
        # Save to the folder on Server under Bnary ID
        plt.savefig('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/BinaryID/Matched_Methods_BinaryCandidates/'+str(locationid)+'_'+str(apogeeid)+'_'+str(visit)+'.pdf', dpi=650,bbox_to_inches='tight')
        # plt.savefig('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/BinaryID/MachineLearning_IDStars_June2018/90-100/'+str(locationid)+'_'+str(apogeeid)+'_'+str(visit)+'.pdf', dpi=650,bbox_to_inches='tight')
        plt.close('all')
        
    return 


## --- MAIN ROUTINE BEGINS HERE --- ##
# Read in the catalog of stars made by JMR visual cuts
# jmr_data = pd.read_csv('JMR_Cuts_RecentRevision.csv')
# locationIDs = np.asarray(jmr_data['LocationID'])
# apogeeIDs = np.asarray(jmr_data['ApogeeID'])
# data = pd.read_csv('TrainingSet_SmallR_LargeXR.csv')
# data = pd.read_csv('KMEVP1-10R51R401NewXR51R401LogAvREDO.csv')
# locationIDs = np.asarray(data['ID'])
# apogeeIDs = np.asarray(data['Star_Name'])
# percentage = np.asarray(data['1']) # Assuming that the index means that the star is binary

# Read in the matched binary candidates from the methods ML, JMR Cuts, and Gaia Cuts
# matches = pd.read_csv('Matched_BinaryCandidates_Updated.csv')
matches = pd.read_csv('BinaryCandidate_Matches.csv') # This is the list that is not constrained to ML targets that have a % > 90
locationIDs = np.asarray(matches['LocationID'])
apogeeIDs = np.asarray(matches['ApogeeID'])
IDs = np.asarray(matches['Bin%'])

#Run routine on DR14 to find R values, R ratios, x-ranges and HJDs
for j in range(len(locationIDs)):
    locationID = locationIDs[j]
    apogeeID = apogeeIDs[j]
    ID = IDs[j]
    print(j)
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
        plt.figure(figsize=(8,8))
        plt.grid(alpha=0.80)
        for visit in range(0,nvisits):
            if nvisits == 1:
                ccf = CCF
            else:
                ccf = CCF[visit+2]
            bs_pt = bisector(xccf, ccf)
            # # split up 2D array to get x and y coord. for bisector pts
            # x_bs = bs_pt[0]
            # y_bs = bs_pt[1]
            # x_range = XRANGE(x_bs)
            plt.plot(ccf+visit,label='Visit: '+str(visit))
            # plt.annotate('binary % = '+str(round(ID,3)),xy=(0,visit+0.07))
            plt.xlabel('CCF Lag',fontsize=17)
            plt.ylabel('CCF Units',fontsize=17)
            plt.legend(loc='upper right')
            plt.title(apogeeID+' Binary%: '+str(round(ID,3)))
            # Save to the folder on Server under Bnary IS
        plt.savefig('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/BinaryID/MatchedCandidates_NotRestricted/'+str(locationID)+'_'+str(apogeeID)+'.pdf', dpi=650,bbox_to_inches='tight')
        # plt.savefig('/Volumes/coveydata/APOGEE_Spectra/APOGEE2_DR14/BinaryID/Matched_Methods_BinaryCandidates/'+str(locationID)+'_'+str(apogeeID)+'.pdf', dpi=650,bbox_to_inches='tight')
        plt.close('all')
            # plot_ccf = ccf_overplot(xccf,ccf,locationID,apogeeID,x_range,ID)
            # plot_ccfs = make_ccfs(apogeeID, locationID, xccf, ccf, x_bs, y_bs, x_range,y_bs,1,ID)
    except FileNotFoundError:
        pass




