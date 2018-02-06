import pandas as pd

import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import kendalltau
import corner

#Import catalog of traingset binaries identified from the IDL routine
# We want to find appropriate thresholds for stars to be considered SB2s
filepath='Binary_Stats.csv'
openfile = pd.read_csv(filepath)
binlocID = openfile['LocationID']
binapoID = openfile['ApogeeID']

#All values are in log space
data1 = pd.DataFrame(columns=['log(R51)','log(R101)','log(R151)','log(max x-range)','log(R 151/R 101)','log(R 101/R 51)'])

R51 = np.asarray(openfile['log(R51)'])
R101 = np.asarray(openfile['log(R101)'])
R151 = np.asarray(openfile['log(R151)'])
XRange = np.asarray(openfile['log(xr)'])
Ratio1 = np.asarray(openfile['log(Ratio1)'])
Ratio2 = np.asarray(openfile['log(Ratio2)'])

data1['log(R51)'] = R51
data1['log(R101)'] = R101
data1['log(R151)'] = R151
data1['log(max x-range)'] = XRange
data1['log(R151/R101)'] = Ratio1
data1['log(R101/R51)'] = Ratio2

'''
Define function to construct kernel density and histogram plots for the following:
1. R 51 vs R101/R51 (Ratio2)
2. R 101 vs R151/R101 (Ratio1)
3. max x-range vs Ratio1
4. max x-range vs Ratio2
5. R 151 vs R 51'''

def XR_vs_R2(x,y):
    #Using seaborn, create density plots
    sns.set_style("whitegrid")

    g = sns.JointGrid(x='log(max x-range)',y='log(R101/R51)',data=data1,xlim=(-0.4,2.5),ylim=(-0,0.2))
    g = g.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',n_levels=15)
    g = g.plot_joint(plt.scatter,s=1,color='black')
    g = g.plot_marginals(sns.distplot,kde=True)
    plt.savefig('XRvsR2.pdf',dpi=900)
    plt.show()

def R51_vs_R2(x,y):
    #Using seaborn, create density plots, create #1
    sns.set_style("whitegrid")
    g = sns.JointGrid(x='log(R51)',y='log(R101/R51)',data=data1,xlim=(0,2.0),ylim=(-0.6,0.3))
    g = g.plot_joint(sns.kdeplot,shade=False,shade_lowest=False,cmap='RdGy',n_levels=20)
    g = g.plot_joint(plt.scatter,s=1,color='black')
    g = g.plot_marginals(sns.distplot,kde=True)
    plt.savefig('R51vsR2.pdf',dpi=800)
    plt.show()

def R101_vs_R1(x,y):
    #Using seaborn, create density plots, create #1
    sns.set_style("whitegrid")
    g = sns.JointGrid(x='log(R101)',y='log(R151/R101)',data=data1,xlim=(0,1.75),ylim=(0.0,0.1))
    g = g.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',n_levels=20)
    g = g.plot_joint(plt.scatter,s=1,color='black')
    g = g.plot_marginals(sns.distplot,kde=True)
    plt.savefig('R101vsR1_zoom1.pdf',dpi=800)
    plt.show()

def XR_vs_R1(x,y): 
    #Create plot for #3
    sns.set_style("whitegrid")
    g = sns.JointGrid(x='log(max x-range)',y='log(R151/R101)',data=data1,xlim=(-0.2,3.0),ylim=(-0.4,0.15))
    g = g.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',n_levels=20)
    g = g.plot_joint(plt.scatter,s=1,color='black')
    g = g.plot_marginals(sns.distplot,kde=True)
    plt.savefig('XRvsR1.pdf',dpi=800)
    plt.show()

#def R151_vs_R51:
sns.set_style("whitegrid")
g = sns.JointGrid(x='log(R151)',y='log(R51)',data=data1)#,xlim=(0,2.0),ylim=(0,2.0))
g = g.plot_joint(sns.kdeplot,shade=True,shade_lowest=False,cmap='RdGy',n_levels=20)
g = g.plot_joint(plt.scatter,s=1,color='black')
g = g.plot_marginals(sns.distplot,kde=True)
plt.savefig('R151vsR51.pdf',dpi=800)
plt.show()