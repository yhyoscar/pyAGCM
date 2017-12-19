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
heats = [[0,180],[0,180], [0,180]]
fstr = 'tropicheat_difflev'

#titles = ['heating at lat20', 'heating at lat10', 'heating at EQ']
#cases = ['TropicHeat_a04', 'TropicHeat_a02', 'TropicHeat_e01']
#heats = [[20,180],[10,180], [0,180]]
#fstr = 'tropicheat_difflat'

vars  = ['PS','DIV']
facs  = [0.01, 1e6]
windflag  = [False,True]
ilevs = [999, 12]
clevs = [np.array([-10,-6,-4,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,4,6,10]) + 1000, 
    np.array([-10,-6,-4,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,4,6,10])/2 ]

for iday in range(31):
    tstr = 'd'+format(iday,'04')+'_00-00-00'
    ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+fstr+'_'+tstr+'.png'
    print ffig
    plt.figure(1, figsize=(len(vars)*8, len(cases)*4))
    for icase in range(len(cases)):
        pdata = '/disk2/yhy/ModelingGroup/XYZ_cases/'+cases[icase]+'/'
        fdata = pdata+cases[icase]+'_'+tstr+'.nc'
        lon = fncread(fdata, 'lon')
        lat = fncread(fdata, 'lat')
        lev = fncread(fdata, 'lev')
        for ivar in range(len(vars)):
            if ilevs[ivar]<900:
                pic = fncread(fdata, vars[ivar])[0,ilevs[ivar]]
                vstr = vars[ivar] + '('+format(lev[ilevs[ivar]], '.0f')+')'
            else:
                pic = fncread(fdata, vars[ivar])[0]
                vstr = vars[ivar]

            pic *= facs[ivar]
            clev = clevs[ivar] + 0
            plt.subplot(len(cases), len(vars), icase*len(vars)+ivar+1)
            plt.contourf(lon, lat, pic, clev, 
                norm = mc.BoundaryNorm(clev, 256), cmap=cm.jet)
            plt.colorbar()
            plt.plot([heats[icase][1]-15, heats[icase][1]+15], [heats[icase][0], heats[icase][0]], 'k-', linewidth=2)
            plt.plot([heats[icase][1], heats[icase][1]], [heats[icase][0]-15, heats[icase][0]+15], 'k-', linewidth=2)            
            plt.title(vstr+' '+tstr+': '+titles[icase], loc='left')
            
            if windflag[ivar]:
                picu = fncread(fdata, 'U')[0,ilevs[ivar]]
                picv = fncread(fdata, 'V')[0,ilevs[ivar]]
                xx,yy = np.meshgrid(lon,lat)
                dpiclon = 2
                dpiclat = 2
                wind  = plt.quiver( lon[::dpiclon],lat[::dpiclat], 
                    picu[::dpiclat,::dpiclon], picv[::dpiclat,::dpiclon], 
                    color='m',scale=1/0.006)
                #qk = plt.quiverkey(np.sqrt(picu**2+picv**2), 0.92, 1.02, 10, r'10 m/s')        

    plt.tight_layout()
    plt.savefig(ffig)
    plt.close(1)
    
    