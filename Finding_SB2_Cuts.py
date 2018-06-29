import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def twod_hist(x,y,a,b,tsb_x,tsb_y,title,xlabel,ylabel):#,c,d,e,f):
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
    plt.figure(figsize=(8,8))
    plt.hist2d(x,y,bins=(80,80),normed=True,cmap=plt.cm.PuBu)
    cb = plt.colorbar()
    cb.set_label('Counts in Bin')
    plt.scatter(a,b,s=8,marker = '^',facecolor='None',edgecolor='orange',label='Binaries from JMR Cuts')
    plt.scatter(tsb_x, tsb_y,s=4,label='Training Set Binaries')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='lower right')
    # plt.ylim(c,d)
    # plt.xlim(e,f)
    # Use a rectangle to highlight the region of interest for population comparisions
    # fc = facecolor of the rectangle, ec = edge color of the rectangle
    # j and k are the coordinates of the bottom left point of the rectangle. The shape builds from this pt
    # l = width of the rectangle
    # m = height of the rectangle
    #rectangle = plt.Rectangle((j,k),l,m,ec='red',fc='None',linewidth=3)
    #plt.gca().add_patch(rectangle)
    #plt.savefig(title+'.pdf',dpi=800,bbox_to_inches='tight')
    plt.show()
    return

# Function to implement Jacob Skinner's cuts
def js_cuts(apogeeid, locationid,R51,R101,R_ratio1,R_ratio2,max_xrange,peak):
    # Define empty space for storing quantities that pass the SB2 cut-off
    sb2_locationid = []
    sb2_apogeeid = []
    sb2_minR51 = []
    sb2_minR101 = []
    sb2_minR_ratio1 = []
    sb2_minR_ratio2 = []
    sb2_maxXR = []
    sb2_peak = []
    for i in range(len(R51)):
        if (R101[i] < 0.83 and 0.06 < R_ratio1[i] < 0.13) or (R51[i] < 0.83 and 0.05 < R_ratio2[i] < 0.2):
            if peak[i] > -0.5:
                if -0.7 < max_xrange[i] < 2.3:
                    if R51[i] > 0.25:
                        if R101[i] > 0.22: 
                            sb2_locationid.append(locationid[i])
                            sb2_apogeeid.append(apogeeid[i])
                            sb2_minR51.append(R51[i])
                            sb2_minR101.append(R101[i])
                            sb2_minR_ratio1.append(R_ratio1[i])
                            sb2_minR_ratio2.append(R_ratio2[i])
                            sb2_maxXR.append(max_xrange[i])
                            sb2_peak.append(peak[i])
    return sb2_locationid, sb2_apogeeid, sb2_minR51, sb2_minR101, sb2_minR_ratio1, sb2_minR_ratio2, sb2_maxXR, sb2_peak

# Function for by-eye cuts of the non-conservative bounaries: JESSICA'S CUTS!
def visual_cuts(apogeeid, locationid,R51,R101,R151,R_ratio1,R_ratio2,max_xrange):
    # Define empty space for storing quantities that pass the SB2 cut-off
    sb2_locationid = []
    sb2_apogeeid = []
    sb2_minR51 = []
    sb2_minR101 = []
    sb2_minR_ratio1 = []
    sb2_minR_ratio2 = []
    sb2_maxXR = []
    sb2_minR151 = []
    for i in range(len(R51)):
        if (0.29 < R101[i] < 0.95 and -0.02 < R_ratio1[i] < 0.9) or ( 0.25<R51[i]<1.0 and -0.77<R_ratio2[i]<0.13):
            if (0.35<max_xrange[i]<1.75 and -0.12<R_ratio2[i]<0.15) or (0.50<max_xrange[i]<1.75 and -0.2<R_ratio1[i]<0.9):
                if (0.50 < max_xrange[i] < 1.75) and  (-0.2 < R_ratio1[i]<0.9 or 0.25<R101[i]<1.15): # COMBINED OR STATEMENT! CHECK THAT THIS IS WORKING!
                        sb2_locationid.append(locationid[i])
                        sb2_apogeeid.append(apogeeid[i])
                        sb2_minR51.append(R51[i])
                        sb2_minR101.append(R101[i])
                        sb2_minR151.append(R151[i])
                        sb2_minR_ratio1.append(R_ratio1[i])
                        sb2_minR_ratio2.append(R_ratio2[i])
                        sb2_maxXR.append(max_xrange[i])
    return sb2_locationid, sb2_apogeeid, sb2_minR51, sb2_minR101, sb2_minR151, sb2_minR_ratio1, sb2_minR_ratio2, sb2_maxXR

# Function to extract column data from csv file
def read_csvfile(filename,header1,header2,header3,header4,header5,header6,header7,header8,header9):
    data = pd.read_csv(filename)
    colm_1 = np.asarray(data[header1])
    colm_2 = np.asarray(data[header2])
    colm_3 = np.asarray(data[header3])
    colm_4 = np.asarray(data[header4])
    colm_5 = np.asarray(data[header5])
    colm_6 = np.asarray(data[header6])
    colm_7 = np.asarray(data[header7])
    colm_8 = np.asarray(data[header8])
    colm_9 = np.asarray(data[header9])
    return colm_1, colm_2, colm_3, colm_4, colm_5, colm_6, colm_7, colm_8, colm_9

# Find the apogee IDs that are in both catalogs by searching for un-unique apogee ids
def matches(x,y,z,a):#,b,c):
    array_a = []
    array_b = []
    array_c = []
    array_d = []
    #array_e = []
    #array_f = []
    x = np.asarray(x) # Bigger array of strings
    y = np.asarray(y) # Smaller array of strings
    for i in range(len(y)):
        if y[i] in x:
            array_a.append(x[i]) # bigger list apogeeid
            array_b.append(y[i])
            array_c.append(z[i]) # bigger list location id
            array_d.append(a[i])
            #array_e.append(b[i])
            #array_f.append(c[i])
    return array_a, array_c # array_b, array_c, array_d, array_e, array_f

# Report candidate stars in a csv file
def csv_writer(filename,header1,header2,header3,header4,header5,header6,header7,header8,header9,x,y,w,v,u,t,s,b,delimiter_choice):
    cols = [header1, header2]
    # Construct dataframe to store variable data
    dataframe = pd.DataFrame(columns=cols)
    dataframe[header1] = x # Location ID
    dataframe[header2] = y # Apogee ID
    dataframe[header3] = w # min R51 (in log space)
    dataframe[header4] = v # min R101 (in log space)
    dataframe[header5] = u # min R151 (in log space)
    dataframe[header6] = t # min R151/R101 (in log space)
    dataframe[header7] = s # min R101/R51 (in log space)
    dataframe[header8] = b # max xrange (in log space)
    dataframe[header9] = 'JRC' # Identifier for JMR Cut
    dataframe.to_csv(filename,sep=delimiter_choice,index_label=False)
    return

## --- ROUTINE BEGINS HERE ---- ##
# Read in the training set catalog of min R, min R ratios, and max x-range.
# All quantities are in log space
dr14_locationid,dr14_apogeeid, dr14_minr51, dr14_minr101, dr14_minr151, dr14_minratio1, dr14_minratio2, dr14_maxXR, dr14_peak = read_csvfile('DR14_SmallR_LargeXR_Revised.csv','LocationID','ApogeeID','log(R51)','log(R101)',
  'log(R151)','log(R151/R101)','log(R101/R51)','log(xr)','log(Peak_value)')

DR14_apogeeid = np.copy(dr14_locationid)
DR14_locationID = np.copy(dr14_apogeeid)
DR14_minR51 = np.copy(dr14_minr51)
DR14_minR101 = np.copy(dr14_minr101)
DR14_minR151 = np.copy(dr14_minr151)
DR14_minRatio1 = np.copy(dr14_minratio1)
DR14_minRatio2 = np.copy(dr14_minratio2)
DR14_MaxXR = np.copy(dr14_maxXR)

# Read in the DR14 catalog of min R, min R rations, and max x-range
# All quantities are in log space
tsb_locationid,tsb_apogeeid, tsb_minr51, tsb_minr101, tsb_minr151, tsb_minratio1, tsb_minratio2, tsb_maxXR, tsb_peak = read_csvfile('TrainingSet_SmallR_LargeXR.csv','LocationID','ApogeeID','log(R51)','log(R101)',
  'log(R151)','log(R151/R101)','log(R101/R51)','log(xr)','log(Peak_value)')


# Take out the bad data via outlier assignment of 2.5
dr14_locationid = dr14_locationid[dr14_maxXR < 2.5]
dr14_apogeeid = dr14_apogeeid[dr14_maxXR < 2.5]
dr14_minr51 = dr14_minr51[dr14_maxXR < 2.5]
dr14_minr101 = dr14_minr101[dr14_maxXR < 2.5]
dr14_minr151 = dr14_minr151[dr14_maxXR < 2.5]
dr14_minratio1 = dr14_minratio1[dr14_maxXR < 2.5]
dr14_minratio2 = dr14_minratio2[dr14_maxXR < 2.5]
dr14_peak_value = dr14_peak[dr14_maxXR < 2.5]
dr14_maxXR = dr14_maxXR[dr14_maxXR < 2.5]

print(len(tsb_apogeeid)) # Total number training set binaries that are unique
print(len(dr14_apogeeid)) # Total number of stars in DR14 that are non-outlier and are unique

# Find the number of stars in DR14 that pass Jacob Skinner's SB2 cuts
jsCuts_locationid, jsCuts_apogeeid, jsCuts_minR51,jsCuts_minR101, jsCuts_minR_ratio1, jsCuts_minR_ratio2,jsCuts_maxXR, jsCuts_peak = js_cuts(dr14_apogeeid,dr14_locationid,dr14_minr51, dr14_minr101,dr14_minratio1,dr14_minratio2,dr14_maxXR,dr14_peak_value)
print(len(jsCuts_apogeeid)) # 2,462 stars found

# Find the number of stars in DR14 that pass Jessica's visually determined SB2 Cuts
jmr_locationid, jmr_apogeeid, jmr_minR51, jmr_minR101, jmr_minR151,jmr_minR_ratio1, jmr_minR_ratio2, jmr_maxXR = visual_cuts(dr14_apogeeid, dr14_locationid, dr14_minr51, dr14_minr101, dr14_minr151, dr14_minratio1, dr14_minratio2, dr14_maxXR)
print(len(jmr_apogeeid)) # 6,714 stars found

# Find the amount of training set binaries that were classified in the two different cuts (JS v JMR) as a ratio of TSB_selected/ TSB
# Send lists found from both cuts into "matches" Function, start with JMR Cuts
jmrCut_match_apogeeid, jmrCut_match_locationid = matches(jmr_apogeeid,tsb_apogeeid, jmr_locationid,tsb_locationid)#,jmr_minr51,jmr_minr101,jmr_minr151,jmr_maxXR)
print(len(jmrCut_match_apogeeid)) # 1025/1025 training set stars identified


# # Send JS Cuts to the "matches" function
# jsCut_match_apogeeid, jsCut_match_locationid = matches(jsCuts_apogeeid, tsb_apogeeid,jsCuts_locationid, tsb_locationid)
# #print(jsCut_match_locationid)
# print(len(jsCut_match_apogeeid)) # 89/ 1025 training set stars identified

# Stars from JMR cuts
jmr_cut_file = csv_writer('JMR_Cuts_RecentRevision.csv','LocationID', 'ApogeeID','min_log(R51)','min_log(R101)','min_log(R151)','min_log(R151/R101)','min_log(R101/R51)','max_log(xr)','Identifier',jmr_locationid, jmr_apogeeid,jmr_minR51, jmr_minR101,jmr_minR151 ,jmr_minR_ratio1, jmr_minR_ratio2, jmr_maxXR,',')

# Stars from JS cuts
# js_cut_file = csv_writer('JS_Cuts_RecentRevision.csv','LocationID', 'ApogeeID','min_log(R51)','min_log(R101)','min_log(R151/R101)','min_log(R101/R51)','max_log(xr)','Identifier',jsCuts_locationid, jsCuts_apogeeid,jsCuts_minR51, jsCuts_minR101, jsCuts_minR_ratio1, jsCuts_minR_ratio2, jsCuts_maxXR,',')


#### PLOTTING RESULTS ######
# Use corner to make population maps: This one is for R101 vs R151/R101
# import corner
# plt.figure(figsize=(8,8))
# plt.grid(alpha=0.5)
# title = '$R_{101}$ vs $R_{151}/R_{101}$'
# xlabel = 'min log($R_{101}$)'
# ylabel = 'min log($R_{151}/R_{101}$)'
# print('FIRST PLOT')
# corner.hist2d(DR14_minR101, DR14_minRatio1,bins=500,plt_contours=True,fill_contours=True,smooth=1.,plot_datapoints=True)
# plt.scatter(jmr_minR101,jmr_minR_ratio1,color='blue',alpha=0.6,s=6,label='Binaries from JMR Cuts')#,marker = '^',facecolor='None',edgecolor='orange',label='Binaries from JMR Cuts')
# plt.scatter(tsb_minr101, tsb_minratio1,s=6,color='red',label='Training Set Binaries')
# plt.xlabel(xlabel,fontsize=18)
# plt.ylabel(ylabel,fontsize=18)
# plt.title(title,fontsize=18)
# plt.axvline(x=0.29) # plot vertical line at 0.29
# plt.axvline(x=1.15) # plot vertical line at 1.04
# plt.axhline(y=-0.2) # plot horizontal line at -0.02
# plt.axhline(y=1.2) # plot horizontal line at 0.9
# plt.xlim(0,1.5)
# plt.ylim(-0.8,0.1)
# plt.legend(loc='lower right',fontsize=18)
# plt.savefig('JMR_Cut_(R101_v_Ratio1).pdf',dpi=800)
# # plt.show()

# # # print('SECOND PLOT')
# # Use corner to make population maps: This one is for R51 vs R101/R51 (JMR Cut)
# plt.figure(figsize=(8,8))
# plt.grid(alpha=0.5)
# title = '$R_{51}$ vs $R_{101}/R_{51}$ for JMR Cuts'
# xlabel = 'min log($R_{51}$)'
# ylabel = 'min log($R_{101}/R_{51}$)'
# corner.hist2d(DR14_minR51, DR14_minRatio2,bins=500,plt_contours=True,fill_contours=True,smooth=1.,plot_datapoints=True)
# plt.scatter(jmr_minR51,jmr_minR_ratio2,color='blue',alpha=0.6,s=6,label='Binaries from JMR Cuts')#,marker = '^',facecolor='None',edgecolor='orange',label='Binaries from JMR Cuts')
# # plt.scatter(jsCuts_minR51, jsCuts_minR_ratio2, s=6,color='purple',label='Binaries from JS Cuts')
# plt.scatter(tsb_minr51, tsb_minratio2,s=6,color='red',label='Training Set Binaries')
# plt.xlabel(xlabel,fontsize=18)
# plt.ylabel(ylabel,fontsize=18)
# plt.title(title,fontsize=18)
# plt.axvline(x=0.25) # plot vertical line at 0.25
# plt.axvline(x=1.05) # plot vertical line at 1.0
# plt.axhline(y=-0.77) # plot horizontal line at -0.77
# plt.axhline(y=0.15) # plot horizontal line at 0.13
# plt.xlim(0,1.5)
# plt.ylim(-1.4,0.18)
# plt.xlim(0.1,1.45)
# plt.legend(loc='lower right',fontsize=20)
# plt.savefig('JMR_Cut_(R51 vs Ratio2).pdf',dpi=800)
# # plt.show()

# # # print('THIRD PLOT')
# # Use corner to make histogram for XR vs R101 for JMR Cuts
# plt.figure(figsize=(8,8))
# plt.grid(alpha=0.5)
# title = 'x-range vs $R_{101}$'
# xlabel = 'max log($\Delta x$)'
# ylabel = 'min log($R_{101}$)'
# corner.hist2d(dr14_maxXR, dr14_minr101,bins=50,plt_contours=True,fill_contours=True,smooth=1.2,plot_datapoints=True)
# plt.scatter(jmr_maxXR,jmr_minR101,color='blue',alpha=0.6,s=6,label='Binaries from JMR Cuts')#,marker = '^',facecolor='None',edgecolor='orange',label='Binaries from JMR Cuts')
# plt.scatter(tsb_maxXR, tsb_minr101,s=6,color='red',label='Training Set Binaries')
# plt.xlabel(xlabel,fontsize=18)
# plt.ylabel(ylabel,fontsize=18)
# plt.title(title,fontsize=18)
# plt.axvline(x= 0.5) 
# plt.axvline(x= 1.8) 
# plt.axhline(y= 0.25) 
# plt.axhline(y= 1.15) 
# # plt.xlim(-0.2,2.0)
# # plt.ylim(-1.4,0.22)
# plt.legend(loc='upper left',fontsize=20)
# plt.savefig('JMR_Cut_(XR vs R101).pdf',dpi=800)
# # plt.show()

# # print("FOURTH PLOT")
# # Use corner to make histogram of xr vs ratio 1 for JMR Cuts
# plt.figure(figsize=(8,8))
# plt.grid(alpha=0.5)
# title = 'x-range vs $R_{151}/R_{101}$'
# xlabel = 'max log($\Delta x$)'
# ylabel = 'min log($R_{151}/R_{101}$)'
# corner.hist2d(dr14_maxXR, dr14_minratio1,bins=200,plt_contours=True,fill_contours=True,smooth=1.5,plot_datapoints=True)
# plt.scatter(jmr_maxXR, jmr_minR_ratio1, s=6,color='blue',alpha=0.7,label='Binaries from JMR Cuts')
# plt.scatter(tsb_maxXR, tsb_minratio1,s=6,color='red',label='Training Set Binaries')
# plt.xlabel(xlabel,fontsize=18)
# plt.ylabel(ylabel,fontsize=18)
# plt.title(title,fontsize=18)
# plt.axvline(x=0.45) 
# plt.axvline(x=1.75) 
# plt.axhline(y=-1.1) 
# plt.axhline(y=0.15) 
# plt.xlim(-0.5,2.5)
# plt.ylim(-1.0,0.2)
# plt.legend(loc='lower right',fontsize=18)
# plt.savefig('JMR_Cut_(XR vs Ratio1).pdf',dpi=800)
# # plt.show()

# # print("FIFTH PLOT")
# # Use corner to make histogram of xr vs ratio 2 for JMR Cuts
# plt.figure(figsize=(8,8))
# plt.grid(alpha=0.5)
# title = 'x-range vs $R_{101}/R_{51}$'
# xlabel = 'max log($\Delta x$)'
# ylabel = 'min log($R_{101}/R_{51}$)'
# corner.hist2d(dr14_maxXR, dr14_minratio2,bins=200,plt_contours=True,fill_contours=True,smooth=1.5,plot_datapoints=True)
# plt.scatter(jmr_maxXR, jmr_minR_ratio2, s=6,color='blue',alpha=0.7,label='Binaries from JMR Cuts')
# plt.scatter(tsb_maxXR, tsb_minratio2,s=6,color='red',label='Training Set Binaries')
# plt.xlabel(xlabel,fontsize=18)
# plt.ylabel(ylabel,fontsize=18)
# plt.title(title,fontsize=18)
# plt.axvline(x=0.45) 
# plt.axvline(x=1.75) 
# plt.axhline(y=-1.1) 
# plt.axhline(y=0.15) 
# plt.xlim(-0.5,2.5)
# plt.ylim(-1.0,0.2)
# plt.legend(loc='lower right',fontsize=18)
# plt.savefig('JMR_Cut_(XR vs Ratio2).pdf',dpi=800)
# plt.show()