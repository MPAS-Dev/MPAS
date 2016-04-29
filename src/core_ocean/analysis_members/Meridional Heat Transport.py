
# coding: utf-8

# In[1]:

###CONFIGURE
# determines the layout of the regional subplots. The default is 3.
subPlotsPerRow = 3
# determines how many of the bins outside of the actual latitude range of the region is being displayed. 
# The default is 3.
regionPlotPadding = 3
# The path to the datafile used. The datafile must contain all the variables of the mocStreamfunction output stream.
dataFilePath = '/Users/nilsfeige/mocStreamfunction.0010-01-01_regionWithRegionBounds_8Months.nc'
# The path of the output pdf file.
outputFilePath = 'mocDaily.pdf'
# Determines if, and, if yes, how many, timesteps are averaged for each plot. The resulting number of plots will be
# numberOfTimestepsInTheDatafileMinusTheValueOfTheTimeAverageVariable.
# The default is 1 (no time averaging)
timeAverage = 5
###END CONFIGURE


# In[2]:

get_ipython().magic('matplotlib inline')

import sys, math
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
from array import *
from math import ceil
set_matplotlib_formats('png', 'pdf')
from colormap import *
from matplotlib.backends.backend_pdf import PdfPages
from datetime import *


# In[3]:

ct = datetime.now()
ds = nc.Dataset(dataFilePath)
plt.rcParams['figure.figsize'] = 18, 12


# In[4]:

mocStreamvalLatAndDepth = ds['mocStreamvalLatAndDepth'][:]
binBoundaries = ds['binBoundaryMocStreamfunction'][:]*180/math.pi
xtime = nc.chartostring(ds['xtime'][:])

regionsInGroup = ds['regionsInGroup'][:]
nRegionsInGroup = ds['nRegionsInGroup'][:]
regionNames = nc.chartostring(ds['regionNames'][:])
regionGroupNames = nc.chartostring(ds['regionGroupNames'][:])
regionMocData = ds['mocStreamvalLatAndDepthRegion'][:]
regionBoundaries = ds['minMaxLatRegion'][:]*180/math.pi

additionalGroupName = getattr(ds,'config_AM_mocStreamfunction_additionalRegion')
nRegionGroups = ds.dimensions['nRegionGroups'].size

for i in range(len(regionNames)):
    regionNames[i] = regionNames[i].strip()
    
for i in range(len(regionGroupNames)):
    regionGroupNames[i] = regionGroupNames[i].strip()
    
regionNumber = -1

for i in range(len(regionGroupNames)):
    if (regionGroupNames[i].decode('utf8') == additionalGroupName):
        curRegionGroup = i

numRegionsInCurGroup = nRegionsInGroup[curRegionGroup]

for i in range(len(xtime)):
    xtime[i] = xtime[i].strip()

rBD = ds['refBottomDepth'][:]*-1

regionBinNumber = [[0 for x in range(numRegionsInCurGroup)] for y in range(2)] 

for i in range(len(binBoundaries)-1): # for every bin
    for j in range(numRegionsInCurGroup): # and every group
        if regionBoundaries[j][0] > binBoundaries[i]: # find min
            regionBinNumber[0][j] = max(0, i - 1 - regionPlotPadding)
        if regionBoundaries[j][1] > binBoundaries[i + 1]: # find max
            regionBinNumber[1][j] = min(len(binBoundaries) - 1, i + regionPlotPadding)
            
print(regionBoundaries.shape)
print(regionBinNumber)
print(xtime.shape)
print(mocStreamvalLatAndDepth.shape)
print(binBoundaries.shape)
print(rBD.shape)
print(regionBoundaries.shape)
print(regionMocData.shape)
print(additionalGroupName)
print(regionsInGroup)
print(nRegionsInGroup)
print(regionNames)
print(regionGroupNames)
print(nRegionGroups)
print(curRegionGroup)
print(numRegionsInCurGroup)


# In[5]:

pp = PdfPages(outputFilePath)
for i in range(len(xtime) - timeAverage):
    mainAverage = mocStreamvalLatAndDepth[i, :, :]
    for j in range(1, timeAverage):
        mainAverage = mainAverage + mocStreamvalLatAndDepth[i + j,:,:]
    mainAverage /= timeAverage
    plt.figure(1, tight_layout=True)
    numSubplotRows = math.ceil(numRegionsInCurGroup / subPlotsPerRow) * 2 + 4
    plt.subplot2grid((numSubplotRows, subPlotsPerRow), (0, 0), colspan=subPlotsPerRow, rowspan=4)
    
    contourSet = plt.contour(binBoundaries[:], rBD[:], mainAverage[:,:], linewidths = 0.5,                              colors="black")
    plt.contourf(binBoundaries[:], rBD[:], mainAverage[:,:])
    plt.set_cmap('plasma')
    cb = plt.colorbar()
    plt.clabel(contourSet, colors="black")
    plt.title('Global MOC by Latitude and Depth [time: ' + xtime[i].decode('utf8') + ']')
    plt.xlabel('latitude [deg]')
    plt.ylabel('Depth [m]')
    subplotCounter = 0
    for j in range(numRegionsInCurGroup):
        mainAverage = regionMocData[i, j, :, :]
        for k in range(1, timeAverage):
            mainAverage = mainAverage + regionMocData[i + k, j, :, :]
        mainAverage /= timeAverage
        curRow = math.floor(j / subPlotsPerRow) * 2 + 4
        curColumn = math.floor(j % subPlotsPerRow)
        plt.subplot2grid((numSubplotRows, subPlotsPerRow), (curRow, curColumn), rowspan=2)
        cs = plt.contour(binBoundaries[regionBinNumber[0][j]:regionBinNumber[1][j]], rBD[:],                          mainAverage[:,regionBinNumber[0][j]:regionBinNumber[1][j]],                          levels=contourSet.levels, linewidths = 0.5, colors="black",                          extent=contourSet.extent, extend='both')
        plt.contourf(binBoundaries[regionBinNumber[0][j]:regionBinNumber[1][j]], rBD[:],                      mainAverage[:,regionBinNumber[0][j]:regionBinNumber[1][j]],                      levels=contourSet.levels, extent=contourSet.extent, extend='both')
        plt.set_cmap('plasma')
        #cb = plt.colorbar()
        plt.clabel(cs, colors="black")
        plt.title('Region:"' + regionNames[regionsInGroup[curRegionGroup][j] - 1].decode('utf8') + '"')
        plt.xlabel('latitude [deg]')
        plt.ylabel('Depth [m]')
    pp.savefig()
    plt.clf()
pp.close()


# In[6]:

ct2 = datetime.now()


# In[7]:

print(ct2 -ct)


# In[ ]:



