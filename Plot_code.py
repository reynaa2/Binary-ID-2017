import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import corner

def corner_hist2d(x,y,a,b,title,xlabel,ylabel):
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

def twod_hist(x,y,a,b,title,xlabel,ylabel):
    # x = x parameter
    # y = y parameter
    # a = comparison x parameter 
    # b = comparision y parameter
    # title, xlabel, and ylabel are strings for the plot labels.
    title = str(title)
    xlabel = str(xlabel)
    ylabel = str(ylabel)
    plt.figure(figsize=(8,8))
    plt.hist2d(x,y,bins=(30,30),normed=True,cmap=plt.cm.PuBu)
    # plt.hist2d(x,y,bins=(500,500),normed=True,cmap=plt.cm.RdYlBu_r)
    #plt.hist2d(x,y,bins=(50,50),norm=colors.LogNorm(),cmap='Blues')
    cb = plt.colorbar()
    cb.set_label('Counts in Bin')
    plt.scatter(a,b,s=8,marker = '^',facecolor='None',edgecolor='grey',label='Training Set Binaries')
    #plt.hist2d(a,b,bins=(100,100),normed=True,cmap=plt.cm.RdYlBu_r)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.legend(loc='upper left')
    plt.legend(loc='lower right')
    plt.ylim(-0.2,0.2)
    #plt.xlim(-0.5,2.5)
    #plt.savefig(title+'.pdf',dpi=800,bbox_to_inches='tight')
    plt.show()
    return    

filename = 'Updated_DR14_Stats.csv'
x = pd.read_csv(filename)
dr14_r51 = np.asarray(x['log(R51)'])
dr14_r101 = np.asarray(x['log(R101)'])
dr14_r151 = np.asarray(x['log(R151)'])
dr14_xr = np.asarray(x['log(xr)'])
dr14_ratio1 = np.asarray(x['log(Ratio1)'])
dr14_ratio2 = np.asarray(x['log(Ratio2)'])

# Read in training set binaries catalog of min R, min R ratios, and max x-range
#otherfile = 'TrainingSet_SmallR_LargeXR.csv'
otherfile = 'TrainingSet_Binary_Stats.csv'
tsb_data = pd.read_csv(otherfile)
tsb_locationid = np.asarray(tsb_data['LocationID'])
tsb_apogeeid = np.asarray(tsb_data['ApogeeID'])
tsb_r51 = np.asarray(tsb_data['log(R51)'])
tsb_r101 = np.asarray(tsb_data['log(R101)'])
tsb_r151 = np.asarray(tsb_data['log(R151)'])
tsb_ratio1 = np.asarray(tsb_data['log(Ratio1)'])
tsb_ratio2 = np.asarray(tsb_data['log(Ratio2)'])
tsb_xr = np.asarray(tsb_data['log(xr)'])
tsb_peakvalue = np.asarray(tsb_data['Peak_value'])

#Impliment if statement to eliminate targets with log(xranges) > 2.5 due to this being our flag for bad data
dr14_r101 = dr14_r101[dr14_xr < 2.5]
dr14_ratio2 = dr14_ratio2[dr14_xr < 2.5]
dr14_ratio1 = dr14_ratio1[dr14_xr < 2.5]
dr14_xr = dr14_xr[dr14_xr < 2.5]


plt.figure(figsize=(8,8))
# Use corner 2d histogram package. Choose 500 bins to emphasize the dense region
corner.hist2d(dr14_xr,dr14_ratio2,bins=100,plt_contours=True,fill_contours=True,smooth=2.0,plot_datapoints=True)
# Apply the scatter plot of training set binaries over the contours
plt.scatter(tsb_xr,tsb_ratio2,marker='^',facecolor='None',edgecolor='red',s=3,label='Training Set Binaries')
plt.legend(loc='lower right')
plt.xlabel('log(max x-range)')
plt.ylabel('log($R_{101}/R_{51}$)')
#Save the figure before showing it
plt.savefig('Comp_XR_vs_R2.pdf',dpi=800,bbox_to_inches='tight')
plt.show()

#Plot that highlights the area binaries tend to populate
fig, ax,

# #Generate a graph for xrange vs R_101
# xr_vs_r101 = twod_hist(dr14_xr,dr14_r101,tsb_xr,tsb_r101,'X ranges vs $R_{101}$','log(max x-range)','log($R_{101}$)')

# #Generate a graph for max x-range vs R_101/R_51 (Ratio 2)
# xr_vs_R2 = twod_hist(dr14_xr,dr14_ratio2,tsb_xr,tsb_ratio2,'X ranges vs $R_{101}_R_{51}$','log(max x-range)','log($R_{101}/R_{51}$)')

# #Generate a graph for max x-range vs R_151/R_101 (Ratio 1)
# xr_vs_R1 = twod_hist(dr14_xr,dr14_ratio1,tsb_xr,tsb_ratio1,'X ranges vs $R_{151}_R_{101}$','log(max x-range)','log($R_{151}/R_{101}$)')


# #Generate a graph for xrange vs R_101
# xr_vs_r101 = corner_hist2d(dr14_xr,dr14_r101,tsb_xr,tsb_r101,'X ranges vs $R_{101}$','log(max x-range)','log($R_{101}$)')

# #Generate a graph for max x-range vs R_101/R_51 (Ratio 2)
# xr_vs_R2 = corner_hist2d(dr14_xr,dr14_ratio2,tsb_xr,tsb_ratio2,'X ranges vs $R_{101}/R_{51}$','log(max x-range)','log($R_{101}/R_{51}$)')

# #Generate a graph for max x-range vs R_151/R_101 (Ratio 1)
# xr_vs_R1 = corner_hist2d(dr14_xr,dr14_ratio1,tsb_xr,tsb_ratio1,'X ranges vs $R_{151}/R_{101}$','log(max x-range)','log($R_{151}/R_{101}$)')

