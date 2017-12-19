from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
matplotlib.rcParams.update({'font.size':18})
import warnings;  warnings.filterwarnings('ignore')
import os

#import glob
#==================================================================
#for case in ['fp_T21L8_s303','fp_T21L8_s304']:
#for case in ['fp_T30L18_d_dp','fp_T30L18_d_sh']:
#for case in ['fp_T30L18_o_dp','fp_T30L18_o_sh']:
#for case in ['fp_T30L18_p_dp','fp_T30L18_p_sh_ptop850']:

cases = ['fp_T30L18_o_md','fp_T30L18_o_md_nophy','fp_T30L18_o_md_decay', 
    'fp_T30L18_o_md_p_dp_k3_rf5_nc8_maxconv2','fp_T30L18_o_md_p_sh_k3_rf5_nc8_maxconv2']

tstrs = ['a) ctrl', 'b) no_phy','c) decay','d) CISK_deep','e) CISK_shallow']

figid = plt.figure(1, figsize=(26,8))

for icase in range(len(cases)):

    case = cases[icase]
    pcases = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'
    #pcases = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'
    import sys;  sys.path.append(pcases+case+'/')
    from namelist import *
    import commands

    casename = case
    #vars = ['T','U', 'PS', 'OMEGA']
    #piclevs = [500,300, -1, 500]
    vars = ['U', 'PS', 'PREC']
    piclevs = [1000, -1, -1]
    facs = [1,0.01, 1]
    units = ['m/s', 'hPa', 'mm/day']
    cms = [cm.RdYlBu_r]*2 + [cm.hot_r]
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
    print range(len(cases)), icase
    
    plt.subplot(1,len(cases)+1,icase+1)
    clev = np.array([0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1, 2, 4, 7, 10])
    plotid = plt.contourf(lon, range(len(pictime)), picdata['PREC'], clev, 
            norm = mc.BoundaryNorm(clev, 256), cmap=cm.hot_r)
    if icase==len(cases)-1:
        position=figid.add_axes([0.84,0.11,0.01,0.82])  
        plt.colorbar(cax=position, ticks=clev)
        plt.subplot(1,len(cases)+1,icase+1)
        
    plt.plot([180,180], [0,400], 'b--', lw=2)
    
    for itime in range(len(pictime)):
        picdata['U'][itime,:] -= np.mean(picdata['U'][itime])
    clev = np.linspace(-20,20,41)
    plt.contour(lon, range(len(pictime)), picdata['U'], clev, colors='k')

    plt.title(tstrs[icase], loc='left')
    plt.xticks(range(0,360,60))
    plt.xlabel('longitude')
    if icase==0: plt.ylabel('time (days)')
    dy = 40
    plt.yticks(range(len(pictime))[::dy], [(time-inittime).days for time in pictime[::dy]])
    plt.axis([0,360,0,len(pictime)])


plt.tight_layout()
pfig = pcases+'../figures/'
os.system('mkdir '+pfig)
ffig = pfig+'mar565_fig2_timelon.png'; print ffig
plt.savefig(ffig)
plt.close(1)

