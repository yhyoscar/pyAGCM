from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import warnings;  warnings.filterwarnings('ignore')

import sys
sys.path.append('/disk2/yhy/Work/tests/pymods')
from mod_std_read import *
from mod_std_fdiag import *
from mod_std_fplot import *
from mod_std_OOP import *

from mod_classes import *

#==================================================================
cases = ['ERA','initERA_e04','initERA_e03','initERA_e05']
runtimes = [0, 0, 0, 0]
vars  = [ ['U component of wind', 'V component of wind', 'Temperature', 'OMEGA'], 
    ['U','V','T','OMEGA'], 
    ['U','V','T','OMEGA'], 
    ['U','V','T','OMEGA']]
plevs = [500,500,500,500]
clevs = [np.linspace(-50,50,21), np.linspace(-50,50,21), np.linspace(220,280,21), 
    np.linspace(-1,1,21)]

nvar  = len(vars[0])
ncase = len(cases)

inittime = 1992110106
era = ERA(inittime)
for itime in range(4*120):
    plt.figure(1, figsize=(8*nvar, 4*ncase))
    for icase in range(ncase):
        for ivar in range(nvar):
            var = vars[icase][ivar]
            if cases[icase]=='ERA':
                pic = era.freadvar(var)[list(era.lev).index(plevs[ivar])]
                lon = era.lon
                lat = era.lat
            else:
                model = ModelFile(cases[icase], runtimes[icase])
                pic = model.freadvar(var, plevs[ivar])
                lon = model.lon
                lat = model.lat
            plt.subplot(ncase, nvar, icase*nvar+ivar+1)
            plt.contourf(lon,lat, pic, clevs[ivar], 
                norm = mc.BoundaryNorm(clevs[ivar], 256), cmap = cm.jet) 
            plt.colorbar()
            fplotcoast(lon)
            plt.title(cases[icase]+': '+var+format(plevs[ivar])+' ('+era.str+')', loc='left')

        runtimes[icase] += 6*3600
    
    plt.tight_layout()
    ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/initERA_latlon_'+era.str+'.png';  print ffig
    plt.savefig(ffig)
    plt.close(1)        

    era.ftimeforward()
    