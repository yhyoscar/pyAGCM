from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
matplotlib.rcParams.update({'font.size':14})
import warnings;  warnings.filterwarnings('ignore')
import os

#import glob
#==================================================================
#for case in ['fp_T21L8_s303','fp_T21L8_s304']:
#for case in ['fp_T30L18_d_dp','fp_T30L18_d_sh']:
#for case in ['fp_T30L18_o_dp','fp_T30L18_o_sh']:
#for case in ['fp_T30L18_p_dp','fp_T30L18_p_sh_ptop850']:
for case in ['fp_T30L18_o_md_p_sh_heatmag10']:

    pcases = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'
    #pcases = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'
    import sys;  sys.path.append(pcases+case+'/')
    from namelist import *
    import commands

    casename = case
    #vars = ['T','U', 'PS', 'OMEGA']
    #piclevs = [500,300, -1, 500]
    vars = ['T','U', 'PS', 'OMEGA', 'PREC']
    piclevs = [500,300, -1, 500, -1]
    facs = [1,1, 0.01, 1, 1]
    units = ['K', 'm/s', 'hPa', 'Pa/s', 'mm/day']
    cms = [cm.RdYlBu_r]*4 + [cm.hot_r]
    #clevs = [np.linspace(251,258,41)] + [ np.linspace(-18,18,41)] + \
    #    [np.linspace(995,1001, 41)*100]  + [np.linspace(-0.15,0.03,41)] + [np.linspace(0,2.5,21)]

    clevs = [np.linspace(251.5,257.5,25)] + [ np.linspace(-16,16,33)] + \
        [np.linspace(990,1001, 23)]  + [np.linspace(-0.2,0.04,25)] + [np.linspace(0,2.7,28)]

    picdata = {var:[] for var in vars}; reallev={}
    pictime = []
    if timeend*foutfreq < 0: ntime = (-timeend*tunit//timestep)//foutfreq
    else: ntime = timeend//foutfreq
    for itime in range(ntime+1):
    #for itime in range(18):
        if foutfreq<0: dsec = -foutfreq*tunit*itime
        else: dsec = foutfreq*timestep*itime
        time = inittime + timedelta(seconds=dsec)
        tstr = format(time.year,'04')+format(time.month,'02')+ \
                format(time.day,'02')+'_'+format(time.hour,'02')+ \
                '-'+format(time.minute,'02')+'-'+format(time.second,'02')
        fn = pcases+casename+'/'+casename+'_'+tstr+'.nc'
        status, fn   = commands.getstatusoutput('ls '+fn)        
        if status == 0:
            pictime += [time]
            print fn
            fid = netcdf(fn, 'r')
            lev = fid.variables['lev'][:]
            lon = fid.variables['lon'][:]
            for ivar in range(len(vars)):
                if piclevs[ivar]<0:
                    data = np.mean(fid.variables[vars[ivar]][:,(nlat/2-1):(nlat/2+1),:], axis=1)
                    reallev[vars[ivar]] = -1
                else:
                    ilev = np.argmin(np.abs(lev - piclevs[ivar]))
                    reallev[vars[ivar]] = lev[ilev]
                    data = np.mean(fid.variables[vars[ivar]][:,ilev,(nlat/2-1):(nlat/2+1),:], axis=1)
                if len(pictime) == 1:
                    picdata[vars[ivar]] = data + 0
                else:
                    picdata[vars[ivar]] = np.append(picdata[vars[ivar]], data, axis=0)

    #=====================================================================
    plt.figure(1, figsize=(5*len(vars),8))
    for i in range(len(vars)):
        plt.subplot(1,len(vars),i+1)
        plt.contourf(lon, range(len(pictime)), picdata[vars[i]], 21, cmap=cm.jet)
        #plt.contourf(lon, range(len(pictime)), picdata[vars[i]] * facs[i], clevs[i], 
        #        norm = mc.BoundaryNorm(clevs[i], 256), cmap=cms[i])
        plt.colorbar(orientation='vertical')
        plt.plot([180,180], [0,300], 'k--', lw=2)
        plt.contour(lon, range(len(pictime)), picdata[vars[i]], [0], colors='m')
        if reallev[vars[i]] > 0:
            plt.title(vars[i]+' ('+units[i]+') '+' @ '+format(reallev[vars[i]],'.0f') +' hPa')
        else:
            plt.title(vars[i]+' ('+units[i]+')')
        plt.xticks(range(0,360,60))
        dy = 8
        plt.yticks(range(len(pictime))[::dy], [(time-inittime).days for time in pictime[::dy]])
        plt.axis([0,360,0,len(pictime)])

    plt.tight_layout()
    pfig = pcases+'../figures/'+casename+'/'
    os.system('mkdir '+pfig)
    ffig = pfig+casename+'_timelon.png'; print ffig
    plt.savefig(ffig)
    plt.close(1)

