import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

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
    plt.figure(figsize=(8,8))
    plt.hist2d(x,y,bins=(60,60),normed=True,cmap=plt.cm.PuBu)
    cb = plt.colorbar()
    cb.set_label('Counts in Bin')
    plt.scatter(a,b,s=8,marker = '^',facecolor='None',edgecolor='red',label='Training Set Binaries')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc='lower right')
    plt.ylim(c,d)
    plt.xlim(e,f)
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
                #if max_xrange[i] > 2.3:
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

# Function for by-eye cuts of the non-conservative bounaries
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
        if (0.29 < R101[i] < 1.04 and -0.02 < R_ratio1[i] < 0.9) or ( 0.25<R51[i]<1.0 and -0.77<R_ratio2[i]<0.13):
            if (0.35<max_xrange[i]<1.75 and -0.12<R_ratio2[i]<0.15) or (0.50<max_xrange[i]<1.75 and -0.2<R_ratio1[i]<0.09):
                if (0.50 < max_xrange[i] < 1.75 and -0.2 < R_ratio1[i]<0.09):
                    if (0.05<max_xrange[i]<2.00 and 0.25<R101[i]<1.15):
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

# Read in the training set catalog of min R, min R ratios, and max x-range.
# All quantities are in log space
dr14_locationid,dr14_apogeeid, dr14_minr51, dr14_minr101, dr14_minr151, dr14_minratio1, dr14_minratio2, dr14_maxXR, dr14_peak = read_csvfile('DR14_SmallR_LargeXR_Revised.csv','LocationID','ApogeeID','log(R51)','log(R101)',
  'log(R151)','log(R151/R101)','log(R101/R51)','log(xr)','log(Peak_value)')

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

print(len(tsb_apogeeid))
# Find the number of stars in DR14 that pass Jacob Skinner's SB2 cuts
# js_results = js_cuts(dr14_apogeeid,dr14_locationid, dr14_minr51, dr14_minr101,dr14_minratio1,dr14_minratio2,dr14_maxXR,dr14_peak_value)
# jsCuts_locationid = js_results[0]
# jsCuts_apogeeid = js_results[1]
# print(len(jsCuts_apogeeid))

jsCuts_locationid, jsCuts_apogeeid, jsCuts_minR51,
jsCuts_minR101, jsCuts_minR_ratio1, jsCuts_minR_ratio2,
jsCuts_maxXR, jsCuts_peak = js_cuts(dr14_apogeeid,dr14_locationid,
dr14_minr51, dr14_minr101,dr14_minratio1,dr14_minratio2,dr14_maxXR,dr14_peak_value)


# JMR_results = visual_cuts(dr14_apogeeid, dr14_locationid, dr14_minr51, dr14_minr101, dr14_minr151, dr14_minratio1, dr14_minratio2, dr14_maxXR)
# jmr_locationid = JMR_results[0]
# jmr_apogeeid = JMR_results[1]
# print(len(jmr_apogeeid))
jmr_locationid, jmr_apogeeid, jmr_minR51, jmr_minR101, jmr_minR151,
jmr_minR_ratio1, jmr_minR_ratio2, jmr_maxXR = visual_cuts(dr14_apogeeid, dr14_locationid, dr14_minr51, dr14_minr101, dr14_minr151, dr14_minratio1, dr14_minratio2, dr14_maxXR)

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
            array_a.append(x[i])
            array_b.append(y[i])
            array_c.append(z[i])
            array_d.append(a[i])
            #array_e.append(b[i])
            #array_f.append(c[i])
    return array_a, array_b #, array_c, array_d, array_e, array_f

# Find the amount of training set binaries that were classified in the two different cuts (JS v JMR)
# Send lists found from both cuts into "matches" Function
jmrCut_match_apogeeid, jmrCut_match_locationid = matches(jmr_apogeeid,tsb_apogeeid, jmr_locationid,tsb_locationid)#,jmr_minr51,jmr_minr101,jmr_minr151,jmr_maxXR)



# Find the associated R, R ratios, x-ranges, and peaks of the matched apogee IDs
#match_apogeeid, match_locationid, match_r51, match_r101, match_r151, match_xr = matches(dr14_apogeeid,tsb_apoid,dr14_minr51,dr14_minr101,dr14_minr151,dr14_maxXR)

# match_locationID, match_apogeeID, match_minR51, match_minR101, match_minR151, match_minRatio1, match_minRatio2, match_xr, match_peaks = matches(all_minR51, all_minR101,all_minR151,all_minR1,all_minR2,all_maxXR,all_peaks,all_locationIDs,all_apogeeIDs)
# Remove the matched elements from the array that contains all the parameter values

# Find number of training stars that were identified by the implimented cuts as a ratio of TSB_selected/ TSB
# Find the number of training set binaries that were found by the implimented cuts as a ratio of the totoal stars selected: TS_selected/ Total Selected
# Report these stars in a csv file
