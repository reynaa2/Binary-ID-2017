import numpy as np
import pandas as pd 
import apogee.tools.read as apread
from matplotlib import pyplot as plt
import pandas as pd
import math
from astropy.io import fits
import os.path
from pathlib import Path

# Report these stars in a csv file
def csv_writer(filename,header1,header2,header3,header4,header5,x,y,w,v,delimiter_choice):
    cols = [header1, header2]
    # Construct dataframe to store variable data
    dataframe = pd.DataFrame(columns=cols)
    dataframe[header1] = x # Location ID
    dataframe[header2] = y # Apogee ID
    dataframe[header3] = w # Visit 
    dataframe[header4] = v # Binary % 
    dataframe[header5] = 1 # ID begins with all stars on 1 out of 0,1,2,3 scale
    dataframe.to_csv(filename,sep=delimiter_choice,index_label=False)
    return

# Assign function to generate lists of stars in particular bounds of ML confidence %
def ML_percentage_select(apogeeIDs, locationIDs, percentage, a,b, n):
    # a = lower % limit
    # b = upper % limit
    # Set up arrays for varaible storage
    binary_percent = []
    star_name = []
    field = []
    visits = []
    if n ==1:
        #Run routine on DR14 to find R values, R ratios, x-ranges and HJDs
        for j in range(len(locationIDs)):
            if percentage[j] >= a:
                locationID = locationIDs[j]
                apogeeID = apogeeIDs[j]
                ID = percentage[j]
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
                    for visit in range(0,nvisits):
                        if nvisits == 1:
                            ccf = CCF
                        else:
                            ccf = CCF[visit+2]
                        binary_percent.append(ID)
                        star_name.append(apogeeID)
                        field.append(locationID)
                        visits.append(visit)
                except FileNotFoundError:
                    pass

    if n != 1:
        for j in range(len(locationIDs)):
            if a <= percentage[j] < b:
                locationID = locationIDs[j]
                apogeeID = apogeeIDs[j]
                ID = percentage[j]
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
                    for visit in range(0,nvisits):
                        if nvisits == 1:
                            ccf = CCF
                        else:
                            ccf = CCF[visit+2]
                        binary_percent.append(ID)
                        star_name.append(apogeeID)
                        field.append(locationID)
                        visits.append(visit)
                except FileNotFoundError:
                    pass
    return field, star_name, visits, binary_percent

## --- MAIN ROUTINE BEGINS HERE --- ##
# Read in the catalog of stars made by ML idneifications
data = pd.read_csv('KMEVP1-10R51R401NewXR51R401LogAvREDO.csv')
locationIDs = np.asarray(data['ID'])
apogeeIDs = np.asarray(data['Star_Name'])
percentage = np.asarray(data['1']) # Assuming that the index means that the star is binary

# Generate a csv file that will store stars that have a confidence % > 90
above90_loc, above90_apo, above90_visit, above90_binPercent = ML_percentage_select(apogeeIDs, locationIDs, percentage, 0.90,1.00,1)
above_90 = csv_writer('ML_90_above.csv','Location_ID','Apogee_ID','Visit', 'Binary%','ID',above90_loc, above90_apo,above90_visit,above90_binPercent,'\t')

# Generate a csv file that will store stars with confidence % 80 < % < 90
Eighty_89_loc, Eighty_89_apo, Eighty_89_visit, Eighty_89_percent = ML_percentage_select(apogeeIDs, locationIDs, percentage, 0.80, 0.89,2)
bw_80_89= csv_writer('ML_80_89.csv','Location_ID','Apogee_ID','Visit', 'Binary%','ID',Eighty_89_loc, Eighty_89_apo,Eighty_89_visit, Eighty_89_percent,'\t')

# Generate csv file for stars with confidence 70< % <80
Seventy_79_loc, Seventy_79_apo, Seventy_79_visit, Seventy_79_percent = ML_percentage_select(apogeeIDs, locationIDs, percentage, 0.70, 0.79,2)
middle_70_79 = csv_writer('ML_70_79.csv','Location_ID','Apogee_ID','Visit', 'Binary%','ID',Seventy_79_loc, Seventy_79_apo,Seventy_79_visit, Seventy_79_percent,'\t')
