import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

# Function to extract column data from csv file
def read_csvfile(filename,header1,header2,delimiter_choice,n): # ,header3,header4,header5,header6,header7,header8):
    data = pd.read_csv(filename,delimiter=delimiter_choice)
    if n == 1:
        colm_1 = np.asarray(data[header1])
        colm_2 = np.asarray(data[header2])
        return colm_1, colm_2
    if n != 1:
        header3 = '1'
        colm_1 = np.asarray(data[header1])
        colm_2 = np.asarray(data[header2])
        colm_3 = np.asarray(data[header3])
        return colm_1, colm_2, colm_3

def matches(x,y,z,a,b,c,id1,id2,id3):
    # x = biggest list of candidates apogee IDs
    # y = 2nd biggest list of candidates apogee IDs
    # z = smallest list of candidates apogee IDs
    # a = big list location id
    # b = 2ng big list location id
    # c = smallest list location id
    match_locationid = []
    match_apogeeid = []
    match_identifier = []
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    print(len(x), len(y), len(z))
    for i in range(len(z)):
        if z[i] in y: # and z[i] in y:
            match_locationid.append(a[i])
            match_apogeeid.append(x[i])
            match_identifier.append(id1[i])
    return match_locationid, match_apogeeid, match_identifier

#
# # Find the number of apogee IDs that are in Gaia list, JR cuts list, and ML list.
# def matches(x,y,z,a,b,c,ID1,ID2,ID3):
#     # x = biggest list of candidates apogee IDs
#     # y = 2nd biggest list of candidates apogee IDs
#     # z = smallest list of candidates apogee IDs
#     # a = big list location id
#     # b = 2ng big list location id
#     # c = smallest list location id
#     match_locationid = []
#     match_apogeeid = []
#     match_identifier = []
#     x = np.asarray(x)
#     y = np.asarray(y)
#     z = np.asarray(z)
#     print(len(x),len(y),len(z))
#
#     for i in range(len(z)):
#         if z[i] in x:
#             # print('I am in ML + JRC')
#             match_locationid.append(c[i])
#             match_apogeeid.append(z[i])
#             if ID1[i] != ID2[i] or ID1[i] != ID3[i]:
#                 match_identifier.append(ID1[i])
#             else:
#                 match_identifier.append(ID2[i])
#         if z[i] in y:
#             print('I am in ML + JRC + Gaia')
#             match_locationid.append(c[i])
#             match_apogeeid.append(z[i])
#             if ID1[i] != ID2[i] or ID1[i] != ID3[i]:
#                 match_identifier.append(ID1[i])
#             else:
#                 match_identifier.append(ID2[i])
#         elif z[i] in x and z[i] in y:
#             print('In both!!')
#
#     return match_locationid, match_apogeeid, match_identifier

# Report these stars in a csv file
def csv_writer(filename,header1,header2,header3,header4,header5,header6,header7,x,y,w,v,u,t,s,delimiter_choice):
    cols = [header1, header2]
    # Construct dataframe to store variable data
    dataframe = pd.DataFrame(columns=cols)
    dataframe[header1] = x # Location ID
    dataframe[header2] = y # Apogee ID
    dataframe[header3] = w # min R51 (in log space)
    dataframe[header4] = v # min R101 (in log space)
    dataframe[header5] = u # min R151/R101 (in log space)
    dataframe[header6] = t # min R101/R51 (in log space)
    dataframe[header7] = s # max xrange (in log space)
    dataframe.to_csv(filename,sep=delimiter_choice,index_label=False)
    return

## ---- MAING ROUTINE ---- ##
# Read in the training set binraies via read_csvfile Function
tsb_locationid,tsb_apogeeid = read_csvfile('TrainingSet_SmallR_LargeXR.csv','LocationID','ApogeeID',',',1)
# tsb_locationid,tsb_apogeeid, tsb_minr51, tsb_minr101, tsb_minr151, tsb_minratio1, tsb_minratio2, tsb_maxXR = read_csvfile('TrainingSet_Binaries_SmallR_LargeXR.csv',
# 'LocationID','ApogeeID','log(R51)','log(R101)','log(R151)','log(R151/R101)','log(R101/R51)','log(xr)')
# Give training set an identifier tag of TS
ts_id = ['TS' for i in range(len(tsb_apogeeid))]

# Read in JMR binary candidates via read_csvfile Function
jrc_locationid,jrc_apogeeid = read_csvfile('JMR_Cuts_RecentRevision.csv','LocationID','ApogeeID',',',1)
# jrc_locationid,jrc_apogeeid, jrc_minr51, jrc_minr101, jrc_minr151, jrc_minratio1, jrc_minratio2, jrc_maxXR = read_csvfile('JMR_Cuts_RecentRevision.csv','LocationID','ApogeeID') #,'min_log(R51)','min_log(R101)','min_log(R151)','min_log(R151/R101)','min_log(R101/R51)','max_log(xr)')
# Assign JMR Cut List the identifier: V (for visual cuts)
jmr_id = ['V' for j in range(len(jrc_apogeeid))]

# Read in ML binary candidates via pandas python package for 90% and above confidence %
# ml_locationid, ml_apogeeid = read_csvfile('ML_90_above.csv','Location_ID', 'Apogee_ID','\t')
# print(len(ml_apogeeid))
ml_loc, ml_apo,ml_percent = read_csvfile('KMEVP1-10R51R401NewXR51R401LogAvREDO.csv','ID','Star_Name',',',3)
ml_locationid = []
ml_apogeeid = []
for n in range(len(ml_loc)):
    if ml_percent[n] >= 0.90:
        ml_locationid.append(ml_loc[n])
        ml_apogeeid.append(ml_apo[n])
print(len(ml_locationid))
# Assign machine learning to have the identifier: M
ml_id = ['M' for k in range(len(ml_apogeeid))]

# Read in KC binary candidates (Gaia list) via pandas python package
gaia_apogeeid = pd.read_csv('binaries_from_gaia.tbl', header=None, delimiter='\t', usecols=[0])
gaia_locationid = pd.read_csv('binaries_from_gaia.tbl', header=None, delimiter='\t', usecols=[2])
gaia_flag = pd.read_csv('binaries_from_gaia.tbl', header=None, delimiter='\t', usecols=[3])

# Assign the Gaia candidates the identifier: G
gaia_id = ['G' for n in range(len(gaia_apogeeid))]

# Find the number of matched apogee IDs and location IDs in KC list, JR List, and ML list
match_locationid, match_apogeeid, match_ID = matches(ml_apogeeid, gaia_apogeeid, jrc_apogeeid, ml_locationid, gaia_locationid, jrc_locationid, ml_id, jmr_id, gaia_id)
print(len(match_apogeeid))
