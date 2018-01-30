import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import csv
import math

# Find the min R values and the max x-ranges for a given star
def matches(y,r51,r101,r151,xr):
    match51 = []
    match101 = []
    match151 = []
    matchxr = []
    uniqueID = np.unique(y)
    for i in range(len(uniqueID)):
    #for i in range(len(uniqueID[0:4])):
        newID = np.where(y == uniqueID[i])
        #print(newID[0])
        #print(r51[newID[0]])
        newer = newID[0]
        small = r51[newer]
        little = min(small)
        x = round(little,4)
        match51.append(x)

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

    return match51, match101, match151, matchxr

#Calculate the ratio arrays

#Calculate the ratio arrays
def ratios(logr51,logr101,logr151):
    r1 = logr51/logr101
    r2 = logr101/logr151
    ratio = r1,r2 
    #ratio = [round(r1,4),round(r2,4)]
    return ratio


#Function for generating denisty population plots
def histPlot(x,y,z,w,v):
    z = str(name) + str(.pdf)
    v = str(parameter1)
    w = str(parameter2)
    plt.figure(figsize=(5,5))
    corner.hist2d(x,y bins=80,plot_contours=True,fill_contours=True,smooth=1.2,plot_datapoints=True)
    plt.xlabel(v)
    plt.ylabel(w)
    plt.title('IDl SB2'+ name)
    plt.show()
    plt.savefig(name,dpi=600)

#Function for writing csv files
def csv_writer(u,v,w,x,y,z,s,t):
    # u,v,w are strings to make new column headers
    # x,y,z,s,t are the arrays to be assigned to columns
    name = str(u)
    name2 = str(v)
    name3 = str(w)
    
    cols = ['LocationID', 'ApogeeID']
    df = pd.DataFrame(columns=cols)
    df['LocationID'] = x
    df['ApogeeID'] = y
    df[name] = z
    df[name2] = s
    df[name3] = t
    df.to_csv('Catalog_'+str(m)+str(.pdf),sep='\t', index=False)
    return 

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
#Write a file for binary parameters
df.to_csv('Binary_Stats.csv',sep='\t')

# Read in the catalog of binaries
b = pd.read_csv('KC_Binary_Stats_Revised.csv',delimiter='\t')
bApoID = b['ApogeeID']
bLocID = b['LocationID']
bxr = b['xRng']
b51 = b['R51']
b101 = b['R101']
b151 = b['R151']

# Find the minimum Rs and maximum x-range
m = matches(bApoID,b51,b101,b151,bxr)
min51 = arrays(m[0])
min101 =arrays(m[1])
min151 = arrays(m[2])
maxXR = arrays(m[3])

#Isolate the arrays from one another and make them log arrays 
logmin51 = arrays(m[0])
logmin101 =arrays(m[1])
logmin151 = arrays(m[2])
logmaxXR = arrays(m[3])

# Find R ratios
Rat = ratios(logmin51,logmin101,logmin151)
logbinR1 = Rat[0]
logbinR2 = Rat[1]

# Make a refined csv with only R101 and x-ranges
col = ['R101'] 
data =  pd.DataFrame(columns=col)
data['R101'] = np.asarray(logmin101)
data['XR'] = np.asarray(logmaxXR)
data.to_csv('BinTest2.csv',index=False,sep='\t')
#print(type(logbinR1),type(logbinR2))

#Check array lengths (NOTE: 122 of the 1071 binaries are missing! the np.unique didn't identify them for some odd reason?)
#print(len(logbinR1),len(logbinR2))
#print(len(logmin51), len(logmin101), len(logmin151))

#Generate figures for visual determination of binary requirements
import corner
plt.figure(figsize=(5,5))
corner.hist2d(logmin51,logbinR2, bins=80,plot_contours=True,fill_contours=True,smooth=1.2,plot_datapoints=True)
plt.xlabel('Min log(51)')
plt.ylabel(' Log R101/R51')
plt.title('IDl SB2 R51 vs Ratio 2')
plt.show()

