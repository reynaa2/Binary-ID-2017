import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv
from astropy.io import fits
import os.path
from pathlib import Path

def matches(x,y,z,switch,w):
    #x = parameter in question that you want to find the minimum value of over all visits for a given apogee ID
    #y = apogeeID
    #z = location id
    # Switch = parameter to switch from min to max; if switch == 1, use max else use min
    # w = peak values
    match_param_array = []
    match_peakvalue = []
    # Return indices of unique sources
    uniqueID,indexes = np.unique(y,return_index=True)
    # Use indexes to find relating apogee and location ids
    match_apogeeid = y[indexes]
    match_locationid = z[indexes]
    for i in range(len(uniqueID)):
        newID = np.where(y == uniqueID[i])
        new_array = newID[0]
        new_param_array = x[new_array]
        new_apo = y[new_array]
        if switch ==0:
            small_val = min(new_param_array)
            loc_small_val = np.where(new_param_array == small_val)
            print('Switch 0: '+str(i))
        else:
            small_val = max(new_param_array)
            loc_small_val = np.where(new_param_array == small_val)
            print('Switch 1: '+str(i))
        small_loc = loc_small_val[0][0]
        small_param = round(small_val,4)
        #array for small parameters
        match_param_array.append(small_param)
        #array for matching apogee ids and location ids
        match_peakvalue.append(w[small_loc])
    return match_param_array,match_apogeeid,match_locationid, match_peakvalue
# Generate 2d histograms using corner   
def corner_hist2d(x,y,a,b,xlabel,ylabel,title):
    import corner
    xlabel = str(xlabel)
    y = str(ylabel)
    title = str(title)
    ## Set a figure size
    plt.figure(figsize=(8,8))
    # Use corner 2d histogram package. Choose 500 bins to emphasize the dense region
    corner.hist2d(x,y,bins=500,plt_contours=True,fill_contours=True,smooth=0.5,plot_datapoints=True)
    # Apply the scatter plot of training set binaries over the contours
    plt.scatter(a,b,marker='^',facecolor='None',edgecolor='red',s=3,label='Training Set Binaries')
    plt.legend(loc='lower right')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #Save the figure before showing it
    plt.savefig(title + '.pdf',dpi=800,bbox_to_inches='tight')
    plt.show()
    return
#Generates 2d histograms with overlapping scatter plot of comparision sources
def twod_hist(x,y,a,b,title,xlabel,ylabel,c,d,e,f):
    # x = x parameter
    # y = y parameter
    # a = comparison x parameter 
    # b = comparision y parameter
    # title, xlabel, and ylabel are strings for the plot labels.
    # c and d = y limits for the plot
    # e and f = x limits for the plot
    title = str(title)
    xlabel = str(xlabel)
    ylabel = str(ylabel)
    c.astype(float)
    d.astype(float)
    e.astype(float)
    f.astype(float)
    plt.figure(figsize=(8,8))
    plt.hist2d(x,y,bins=(50,50),normed=True,cmap=plt.cm.PuBu)
    cb = plt.colorbar()
    cb.set_label('Counts in Bin')
    plt.scatter(a,b,s=8,marker = '^',facecolor='None',edgecolor='black',label='Training Set Binaries')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='lower right')
    plt.ylim(c,d)
    plt.xlim(e,f)
    plt.savefig(title+'.pdf',dpi=800,bbox_to_inches='tight')
    plt.show()
    return

## Read in DR14 catalog generated by APOGEE_ID.py
filepath2 = 'DR14_Stats_Catalog_Revision.csv'
data_dr14 = pd.read_csv(filepath2)
# Read in columns of data we want for the function "matches"
locationid = np.asarray(data_dr14['LocationID'])
apogeeid = np.asarray(data_dr14['ApogeeID'])
visits = np.asarray(data_dr14['Visit'])
snr = np.asarray(data_dr14['SNR'])
R51 = np.asarray(data_dr14['log(R51)'])
R101 = np.asarray(data_dr14['log(R101)'])
R151 = np.asarray(data_dr14['log(R151)'])
xranges = np.asarray(data_dr14['log(xr)'])
ratio1 = np.asarray(data_dr14['log(Ratio1)'])
ratio2 = np.asarray(data_dr14['log(Ratio2)'])
peak_value = np.asarray(data_dr14['Peak_value'])

#Send parameters to function "matches"
print('FINDING MIN R 51')
dr14 = matches(R51,apogeeid,locationid,0,peak_value)
dr14_peakvalues = dr14[3]
dr14_locationid = dr14[2]
dr14_apogeeid = dr14[1]
dr14_R51 = dr14[0]
print('FINDING MIN R 101')
r101 = matches(R101,apogeeid,locationid,0,peak_value)
dr14_R101= r101[0]
print('FINDING MIN R 151')
r151 = matches(R151,apogeeid,locationid,0,peak_value)
dr14_R151 = r151[0]
print('FINDING MIN R151/R101')
ratio_1 = matches(ratio1,apogeeid,locationid,0,peak_value)
dr14_R1 = ratio_1[0]
print('FINDING MIN R101/R51')
ratio_2 = matches(ratio2,apogeeid,locationid,0,peak_value)
dr14_R2 = ratio_2[0]
# print('FINDING MAX X RANGE')
# xr = matches(xranges,apogeeid,locationid,1,peak_value)
# dr14_xr = xr[0]

#Impliment if statement to eliminate targets with log(xranges) > 2.5 due to this being our flag for bad data
dr14_locationid = dr14_locationid[dr14_xr < 2.5]
dr14_apogeeid = dr14_apogeeid[dr14_xr < 2.5]
dr14_R51 = dr14_R51[dr14_xr < 2.5]
dr14_R101 = dr14_R101[dr14_xr < 2.5]
dr14_R151 = dr14_R151[dr14_xr < 2.5]
dr14_R1 = dr14_R1[dr14_xr < 2.5]
dr14_R2 = dr14_R2[dr14_xr < 2.5]
dr14_xr = dr14_xr[dr14_xr < 2.5]
dr14_peakvalues = dr14_peakvalues[dr14_xr < 2.5]

#Construct a .csv file to hold smallest R, R ratios, and max x-ranges
cols = ['LocationID', 'ApogeeID']
df = pd.DataFrame(columns = cols)
df['LocationID'] = dr14_locationid
df['ApogeeID'] = dr14_apogeeid
#add HJD info to output table
df['log(R51)'] = dr14_R51
df['log(R101)'] = dr14_R101
df['log(R151)'] = dr14_R151
#df['log(xr)'] = dr14_xr
df['log(R151/R101)'] = dr14_R1
df['log(R101/R51)'] = dr14_R2
df['log(Peak_value)'] = dr14_peakvalues
df.to_csv('Revised_DR14_SmallR_LargeXR.csv')

### -------- CODE FOR MAKING ACTUAL PLOTS ----- ###

# # Read in the DR14 catalog of minimum R, min R ratios, and max x-range values
# filename = 'DR14_SmallR_LargeXR.csv'
# data = pd.read_csv(filename)
# dr14_locationid = np.asarray(data['LocationID'])
# dr14_apogeeid = np.asarray(data['ApogeeID'])
# dr14_r51 = np.asarray(data['log(R51)'])
# dr14_r101 = np.asarray(data['log(R101)'])
# dr14_r151 = np.asarray(data['log(R151)'])
# dr14_ratio1 = np.asarray(data['log(Ratio1)'])
# dr14_ratio2 = np.asarray(data['log(Ratio2)'])
# dr14_xr = np.asarray(data['log(xr)'])
# dr14_peak_value = np.asarray(data['Peak_value'])

# # Read in training set binaries catalog of min R, min R ratios, and max x-range
# otherfile = 'TrainingSet_SmallR_LargeXR.csv'
# tsb_data = pd.read_csv(otherfile)
# tsb_locationid = np.asarray(tsb_data['LocationID'])
# tsb_apogeeid = np.asarray(tsb_data['ApogeeID'])
# tsb_r51 = np.asarray(tsb_data['log(R51)'])
# tsb_r101 = np.asarray(tsb_data['log(R101)'])
# tsb_r151 = np.asarray(tsb_data['log(R151)'])
# tsb_ratio1 = np.asarray(tsb_data['log($R_{151}/ R_{101}$)'])
# tsb_ratio2 = np.asarray(tsb_data['log($R_{101}/R_{51}$)'])
# tsb_xr = np.asarray(tsb_data['log(xr)'])
# tsb_peakvalue = np.asarray(tsb_data['Peak_value'])


# #Generate a graph for xrange vs R_101
# xr_vs_r101 = twod_hist(dr14_xr,dr14_r101,tsb_xr,tsb_r101,'X ranges vs $R_{101}$','log(max x-range)','log($R_{101}$)')

# #Generate a graph for max x-range vs R_101/R_51 (Ratio 2)
# xr_vs_R2 = twod_hist(dr14_xr,dr14_ratio2,tsb_xr,tsb_ratio2,'X ranges vs $R_{101}/R_{51}$','log(max x-range)','log($R_{101}/R_{51}$)')

# #Generate a graph for max x-range vs R_151/R_101 (Ratio 1)
# xr_vs_R1 = twod_hist(dr14_xr,dr14_ratio1,tsb_xr,tsb_ratio1,'X ranges vs $R_{151}/R_{101}$','log(max x-range)','log($R_{151}/R_{101}$)')

# #Generate a graph for R_101 vs R_151/R_101 (Ratio 1)
# r101_vs_R1 = twod_hist(dr14_r101,dr14_ratio1,tsb_r101,tsb_ratio1,'R101 vs R151_R101','log($R_{101}$)', 'log($R_{151}/R_{101}$)')

# #Generate a graph for R_51 vs R_101/R_51
# r51_vs_R2 = twod_hist(dr14_r51,dr14_ratio2,tsb_r51,tsb_ratio2,'R51 vs R101_R51','log($R_{51}$)','log($R_{101}/R_{51}$)')

# ### ---- TRY MAKING A DENSITY PLOT ----- ###
# #Generate a graph for R_51 vs R_101/R_51 (Ratio 2)
# r51_vs_R2 = corner_hist2d(r14_r51,dr14_ratio2,tsb_r51,tsb_ratio2,'R51 vs R101_R51','log($R_{51}$)','log($R_{101}/R_{51}$)',-3,3,-2,5)

# #Generate a graph for R_101 vs R_151/R_101 (Ratio 1)
# r101_vs_R1 = corner_hist2d(dr14_r101,dr14_ratio1,tsb_r101,tsb_ratio1,'R101 vs R151_R101','log($R_{101}$)', 'log($R_{151}/R_{101}$)',-3,3,-2,5)

# #Generate a graph for xrange vs R_101
# xr_vs_r101 = twod_hist(dr14_xr,dr14_r101,tsb_xr,tsb_r101,'X ranges vs R_101','log(max x-range)','log($R_{101}$)',-3,3,-2,5)

# #Generate a graph for max x-range vs R_101/R_51 (Ratio 2)
# xr_vs_R2 = twod_hist(dr14_xr,dr14_ratio2,tsb_xr,tsb_ratio2,'X ranges vs R101_R51','log(max x-range)','log($R_{101}/R_{51}$)',-3,3,-2,5)

# #Generate a graph for max x-range vs R_151/R_101 (Ratio 1)
# xr_vs_R1 = twod_hist(dr14_xr,dr14_ratio1,tsb_xr,tsb_ratio1,'X ranges vs R151_R101','log(max x-range)','log($R_{151}/R_{101}$)',-3,3,-2,5)