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

from datetime import datetime, timedelta

#==================================================================
titles = ['deep mode','shallow mode']
cases = ['fp_T30L18_s_dp','fp_T30L18_s_sh']
heats = [[0,180],[0,180]]
fstr = 'mar565_fig1'

#titles = ['heating at lat20', 'heating at lat10', 'heating at EQ']
#cases = ['TropicHeat_a04', 'TropicHeat_a02', 'TropicHeat_e01']
#heats = [[20,180],[10,180], [0,180]]
#fstr = 'tropicheat_difflat'

vars  = ['PS','DIV']
facs  = [0.01, 1e6]
windflag  = [False,True]
ilevs = [999, -1]
clevs = [np.array([-12,-8,-6,-4,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,4,6,8, 12]) + 1000, 
    np.array([-10,-6,-4,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,4,6,10]) *0.8 ]
cms = [cm.jet, cm.coolwarm]

for iday in range(299,301):
    time = datetime(1,1,1) + timedelta(days=iday)
    tstr = format(time.year,'04')+format(time.month,'02')+ \
                format(time.day,'02')+'_'+format(time.hour,'02')+ \
                '-'+format(time.minute,'02')+'-'+format(time.second,'02')

    ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+fstr+'_'+tstr+'.png'
    print ffig
    plt.figure(1, figsize=(len(vars)*6, len(cases)*3))
    for icase in range(len(cases)):
        pdata = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'+cases[icase]+'/'
        fdata = pdata+cases[icase]+'_'+tstr+'.nc'
        lon = fncread(fdata, 'lon')
        lat = fncread(fdata, 'lat')
        lev = fncread(fdata, 'lev')
        for ivar in range(len(vars)):
            if ilevs[ivar]<900:
                pic = fncread(fdata, vars[ivar])[0,ilevs[ivar]]
                vstr = vars[ivar] + '('+format(lev[ilevs[ivar]], '.0f')+' hPa)'
            else:
                pic = fncread(fdata, vars[ivar])[0]
                vstr = vars[ivar]

            pic *= facs[ivar]
            clev = clevs[ivar] + 0
            plt.subplot(len(cases), len(vars), icase*len(vars)+ivar+1)
            plt.contourf(lon, lat, pic, clev, 
                norm = mc.BoundaryNorm(clev, 256), cmap=cms[ivar])
            #plt.colorbar()
            plt.plot([heats[icase][1]-15, heats[icase][1]+15], [heats[icase][0], heats[icase][0]], 'm-', linewidth=2)
            plt.plot([heats[icase][1], heats[icase][1]], [heats[icase][0]-15, heats[icase][0]+15], 'm-', linewidth=2)            
            plt.title(vstr+' '+tstr+': '+titles[icase], loc='left')
            
            if windflag[ivar]:
                picu = fncread(fdata, 'U')[0,ilevs[ivar]]
                picv = fncread(fdata, 'V')[0,ilevs[ivar]]
                xx,yy = np.meshgrid(lon,lat)
                dpiclon = 2
                dpiclat = 2
                wind  = plt.quiver( lon[::dpiclon],lat[::dpiclat], 
                    picu[::dpiclat,::dpiclon], picv[::dpiclat,::dpiclon], 
                    color='k',scale=1/0.006)
                #qk = plt.quiverkey(np.sqrt(picu**2+picv**2), 0.92, 1.02, 10, r'10 m/s')        

            if ivar == 0:
                plt.axis([0,360,-90,90])
            if ivar == 1:
                plt.axis([130,280,-30,30])

    plt.tight_layout()
    plt.savefig(ffig)
    plt.close(1)
    
    
