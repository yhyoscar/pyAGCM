from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import warnings;  warnings.filterwarnings('ignore')

import glob
import os
import sys

casename = 'HS_T42L26'
pdata = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'+casename+'/'
sys.path.append(pdata)
from namelist import *

pfig  = pdata
vars  = ['U','V','T','OMEGA']
clevs = [np.linspace(-20,48,18)] + [np.linspace(-4,4,9)] + [np.linspace(180,320,15)] + [np.linspace(-0.04,0.04,9)]
ncol  = 2; nrow = int(np.ceil(1.0*len(vars)/ncol))
time  = [ [200,0,0,0], [1200,0,0,0] ]

#======================================================================
tstr = 'd'+format(time[0][0],'04')+'-'+format(time[0][1],'02')+'-'+format(time[0][2],'02')+'-'+format(time[0][3],'02')
if len(time)==2: tstr += 'd'+format(time[1][0],'04')+'-'+format(time[1][1],'02')+'-'+format(time[1][2],'02')+'-'+format(time[1][3],'02')

t0 = time[0][0]*86400 + time[0][1]*3600 + time[0][2]*60 + time[0][3] 
t1 = time[-1][0]*86400 + time[-1][1]*3600 + time[-1][2]*60 + time[-1][3] 
pics = [np.zeros([nlev, nlat]) for i in range(len(vars))]
units = ['' for i in range(len(vars))]
names = ['' for i in range(len(vars))]

files = glob.glob(pdata+casename+'_d*.nc')
files.sort(key=os.path.basename)
nfile = 0
for i in range(len(files)):
    fid = netcdf(files[i], 'r')
    t   = fid.variables['time'][:]
    if (t>=t0 and t<=t1):
        for ivar in range(len(vars)):
            print 'from: ', files[i], ', read ',vars[ivar]
            data = fid.variables[vars[ivar]]
            pics[ivar] += np.mean(data[:], axis=(0,3))
            if nfile==0:
                lat  = fid.variables['lat'][:]
                lev  = fid.variables['lev'][:]
                units[ivar] = data.units
                names[ivar] = data.long_name
        nfile += 1
    fid.close()

#=======================================================================
plt.figure(1, figsize=(8*ncol, 4*nrow))
for ivar in range(len(vars)):
    plt.subplot(nrow,ncol,ivar+1)
    if len(clevs[ivar])>1:
        plt.contourf(lat,lev, pics[ivar]/nfile, clevs[ivar], 
            norm = mc.BoundaryNorm(clevs[ivar], 256), cmap = cm.jet)
        plt.colorbar()
        plt.contour(lat,lev, pics[ivar]/nfile, clevs[ivar], colors='k')        
    else:
        plt.contourf(lat,lev,pics[ivar]/nfile)
        plt.colorbar()
        
    plt.title(tstr+': '+names[ivar]+' ('+units[ivar]+')', loc='left')
    plt.gca().invert_yaxis()
    
plt.tight_layout()
ffig = pfig+'fig_levlat_'+tstr+'.png'
print ffig
plt.savefig(ffig)
plt.close(1)
    
    
    
    
    