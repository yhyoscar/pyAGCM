from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import warnings;  warnings.filterwarnings('ignore')

#import glob
#import os

import sys
sys.path.append('/disk2/yhy/Work/tests/pymods')
from mod_std_read import *
from mod_std_fdiag import *
from mod_std_fplot import *
from mod_std_OOP import *


#==================================================================
titles = ['heating at 250hPa','heating at 500hPa','heating at 750hPa']
cases = ['TropicHeat_e03','TropicHeat_e01','TropicHeat_e02']
heats = [[180,0,250], [180,0,500], [180,0,750]]
fstr = 'tropicheat_cross_difflev'

#titles = ['heating at lat15','heating at lat10','heating at lat5', 'heating at EQ']
#cases = ['TropicHeat_a03','TropicHeat_a02','TropicHeat_a01','TropicHeat_e01']
#heats = [[180], [10,180], [5,180], [0,180]]
#fstr = 'tropicheat_cross_difflev'

var  = 'OMEGA'
fac  = 1
clev = np.linspace(-1,1,21) * 0.1

for iday in range(31):
    tstr = 'd'+format(iday,'04')+'_00-00-00'
    ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+fstr+'_'+tstr+'.png'
    print ffig
    plt.figure(1, figsize=(2*8, len(cases)*4))
    for icase in range(len(cases)):
        pdata = '/disk2/yhy/ModelingGroup/XYZ_cases/'+cases[icase]+'/'
        fdata = pdata+cases[icase]+'_'+tstr+'.nc'
        lon = fncread(fdata, 'lon')
        lat = fncread(fdata, 'lat')
        lev = fncread(fdata, 'lev')
        nlat = len(lat);  nlon = len(lon)
        temp = fncread(fdata, var) 
        for ipic in range(2):
            if ipic==0:
                pic  = 0.5*(temp[0,:,nlat/2-1,:] + temp[0,:,nlat/2,:])
                wind = fncread(fdata, 'U')
                picw = 0.5*(wind[0,:,nlat/2-1,:] + wind[0,:,nlat/2,:])
                x = lon + 0
                y = lev + 0
                dpicx = 2
            if ipic==1:
                pic = temp[0,:,:,nlon/2] + 0
                wind = fncread(fdata, 'V')
                picw = wind[0,:,:,nlon/2] + 0
                x = lat + 0
                y = lev + 0
                dpicx = 1

            plt.subplot(len(cases), 2, icase*2+ipic+1)
            plt.contourf(x, y, pic, clev, 
                norm = mc.BoundaryNorm(clev, 256), cmap=cm.jet)
            plt.colorbar()
            plt.plot([heats[icase][ipic]-15, heats[icase][ipic]+15], [heats[icase][2], heats[icase][2]], 'k-', linewidth=2)
            plt.plot([heats[icase][ipic], heats[icase][ipic]], [heats[icase][2]-200, heats[icase][2]+200], 'k-', linewidth=2)            
            plt.title(var+': '+titles[icase], loc='left')
            plt.gca().invert_yaxis()
            
            
            dpicy = 1
            wind  = plt.quiver( x[::dpicx], y[::dpicy], 
                picw[::dpicy,::dpicx], pic[::dpicy,::dpicx]*(-50), 
                color='m',scale=1/0.003)
            #qk = plt.quiverkey(np.sqrt(picu**2+picv**2), 0.92, 1.02, 10, r'10 m/s')        

    plt.tight_layout()
    plt.savefig(ffig)
    plt.close(1)
    
    