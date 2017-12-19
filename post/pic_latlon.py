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

casename = 'initERA_e01'
pdata = '/disk2/yhy/ModelingGroup/XYZ_cases/'+casename+'/'
sys.path.append(pdata)
from namelist import nlon,nlat,nlev


pfig  = pdata
vars  = ['U','V','T','OMEGA']
levs  = [500,500,500,500]
clevs = ['','','',''] 
ncol  = 1; nrow = int(np.ceil(1.0*len(vars)/ncol))

time  = [ [0,0,0,0] ]
alim  = [0,360,-90,90]

figw = 9; figh = 5
#======================================================================
tstr = casename+'_'
tstr += 'd'+format(time[0][0],'04')+'-'+format(time[0][1],'02')+'-'+format(time[0][2],'02')+'-'+format(time[0][3],'02')
if len(time)==2: tstr += 'd'+format(time[1][0],'04')+'-'+format(time[1][1],'02')+'-'+format(time[1][2],'02')+'-'+format(time[1][3],'02')

print tstr

t0 = time[0][0]*86400 + time[0][1]*3600 + time[0][2]*60 + time[0][3] 
t1 = time[-1][0]*86400 + time[-1][1]*3600 + time[-1][2]*60 + time[-1][3] 
pics = [np.zeros([nlat, nlon]) for i in range(len(vars))]
units = ['' for i in range(len(vars))]
names = ['' for i in range(len(vars))]

files = glob.glob(pdata+'*_d*.nc')
files.sort(key=os.path.basename)
nfile = 0
for i in range(len(files)):
    fid = netcdf(files[i], 'r')
    t   = fid.variables['time'][:]
    if (t>=t0 and t<=t1):
        for ivar in range(len(vars)):
            print 'from: ', files[i], ', read ',vars[ivar]
            data = fid.variables[vars[ivar]]
            hyam = fid.variables['hyam'][:]
            hybm = fid.variables['hybm'][:]
            ps = np.exp(fid.variables['PS'][:][0])
            if len(data[:][0].shape) == 2:
                pics[ivar] += data[:][0]
            else:
                for ilev in range(nlev-1):
                    flag1 = np.zeros([nlat,nlon]).astype(int)
                    flag2 = np.zeros([nlat,nlon]).astype(int)
                    p1 = hyam[ilev]*1e5+hybm[ilev]*ps
                    p2 = hyam[ilev+1]*1e5+hybm[ilev+1]*ps
                    flag1[p1 <= levs[ivar]*100] = 1
                    flag2[p2 >= levs[ivar]*100] = 1
                    pics[ivar][flag1*flag2==1] += data[:][0,ilev][flag1*flag2==1] + (data[:][0,ilev+1][flag1*flag2==1] - data[:][0,ilev][flag1*flag2==1]) / (p2[flag1*flag2==1]-p1[flag1*flag2==1]) * (levs[ivar]*100-p1[flag1*flag2==1])
                    
##                for ilat in range(nlat):
##                    for ilon in range(nlon):
##                        pics[ivar][ilat,ilon] = np.interp(levs[ivar]*100, 
##                            hyam*1e5+hybm*ps[ilat,ilon], data[:][0,:,ilat,ilon])
            if nfile==0:
                lat  = fid.variables['lat'][:]
                lon  = fid.variables['lon'][:]
                units[ivar] = data.units
                names[ivar] = data.long_name
        nfile += 1
    fid.close()
for ivar in range(len(vars)):
    pics[ivar] /= nfile
#=======================================================================
plt.figure(1, figsize=(figw*ncol, figh*nrow))
for ivar in range(len(vars)):
    plt.subplot(nrow,ncol,ivar+1)
    if len(clevs[ivar]) > 1:
        plt.contourf(lon,lat, pics[ivar], clevs[ivar], 
            norm = mc.BoundaryNorm(clevs[ivar], 256), cmap = cm.jet)
    else:
        plt.contourf(lon,lat, pics[ivar], 21, cmap = cm.jet)        
    plt.colorbar()
    plt.contour(lon,lat, pics[ivar], [0], colors='k')    
    plt.axis(alim)
    if levs[ivar]==0: plt.title(tstr+': '+names[ivar]+' ('+units[ivar]+')', loc='left')
    if levs[ivar]>0:  plt.title(tstr+': '+names[ivar]+' @'+format(levs[ivar])+'hPa ('+units[ivar]+')', loc='left')
    
plt.tight_layout()
ffig = pfig+'fig_latlon_'+tstr+'.png'
print ffig
plt.savefig(ffig)
plt.close(1)



