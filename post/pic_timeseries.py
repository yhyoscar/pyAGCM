import pygrib

from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import os


import sys
sys.path.append('/disk2/yhy/Work/tests/pymods')
from mod_std_read import *
from mod_std_fdiag import *
from mod_std_fplot import *
from mod_std_OOP import *

from mod_classes import *

#==============================================================================        

runtime  = 0
nt       = 7*4
#box    = [250,260,-10,-2]
#rstr   = 'SEP_ITCZ'
box  = [ 152, 157, -4, 0 ]
rstr = 'TOGA_COARE'
case = 'initERA_e03'

vars   = ['OMEGA']#, 'Relative humidity']
facs   = [1, 1]
units  = ['Pa/s', '%']
clevs  = [ np.linspace(-0.5,0.5,21)/2, 
    np.linspace(0,100,21) ]

nvar = len(vars)
plt.figure(1, figsize=(10, 3*nvar))
ncol = 4
#fig, axes = plt.subplots(nrows=nvar, ncols=ncol, figsize=(5.0*nt/4/30+4, 3*nvar))


for ivar in range(nvar):
    runtime = 0
    for itime in range(nt):
        model = ModelFile(case, runtime)
        pic = np.zeros([nt, len(model.lev)])
        data = model.freadvar(vars[ivar]) * facs[ivar]
        lon = model.lon
        lat = model.lat
        lev = model.lev
        pic[itime] = fboxmean(data,lon,lat,box,landmask=False)
        
        runtime += 6*3600

    plt.subplot2grid((nvar,ncol), (ivar,1), colspan=ncol-2)
    plt.contourf(range(nt), lev, np.swapaxes(pic,0,1), clevs[ivar], 
        norm = mc.BoundaryNorm(clevs[ivar], 256), cmap = cm.jet)
    plt.colorbar()
    #plt.title(rstr+': '+vars[ivar]+' ('+units[ivar]+') [1992110106, +30days]', loc='left')
    plt.axis([0,nt-1,100,1000])
    plt.gca().invert_yaxis()

    plt.subplot2grid((nvar,ncol), (ivar, 0), colspan=1)
    #cbar_ax = fig.add_axes([0.85, 1-1.0/nvar*(ivar+1), 0.05, 1.0/(nvar+1)])
    plt.axis('off')
    #    plt.colorbar()
    
    plt.subplot2grid((nvar,ncol), (ivar,ncol-1), colspan=1)
    plt.plot(np.mean(pic, axis=0), lev, linewidth=2)
    plt.axis([min(np.mean(pic, axis=0)), max(np.mean(pic, axis=0)), 100,1000])
    plt.gca().invert_yaxis()

plt.tight_layout()
ffig = '../figures/timeseries/'+case+'_'+rstr+'.png'
plt.savefig(ffig)
plt.close(1)


