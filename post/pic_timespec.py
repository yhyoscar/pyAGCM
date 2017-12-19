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

from scipy.fftpack import fft, ifft    

#import glob
#==================================================================
#for case in ['fp_T21L8_s303','fp_T21L8_s304']:
#for case in ['fp_T30L18_s_dp','fp_T30L18_s_sh']:
#for case in ['fp_T30L18_d_dp','fp_T30L18_d_sh']:
#for case in ['fp_T30L18_o_dp','fp_T30L18_o_sh']:
#for case in ['fp_T30L18_p_dp','fp_T30L18_p_sh_ptop850']:
for case in ['fp_T30L18_d_spec']:
    
    pcases = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'
    #pcases = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'
    import sys;  sys.path.append(pcases+case+'/')
    from namelist import *
    import commands

    casename = case
    vars = ['T', 'U', 'V', 'VOR', 'DIV','OMEGA', 'PS']
    piclevs = [500, 300, 300, 300, 300, 300, -1]

    picdata = {var:[] for var in vars}
    pictime = []
    reallev = {}
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
    nlon = len(lon)
    plt.figure(1, figsize=(20*5, 40))
    for i in range(5):
        if i==0: 
            pic = fft(picdata['T'], axis=1)/len(lon)
            tstr = 'T @ 500'
        if i==1: 
            pic = fft(picdata['PS'], axis=1)/len(lon)
            tstr = 'PS'
        if i==2: 
            pic = fft(picdata['U']**2+picdata['V']**2, axis=1)/len(lon)
            tstr = 'KE @ 300'
        if i==3: 
            pic = fft(picdata['DIV']**2, axis=1)/len(lon)
            tstr = 'DIV**2 @ 300'
        if i==4: 
            pic = fft(picdata['OMEGA']**2, axis=1)/len(lon)
            tstr = 'OMEGA**2 @ 300'
        plt.subplot(1,5,i+1)
        ks = np.linspace(-nlon/2,nlon/2-1,nlon)
        datafft = np.log(np.sqrt(0.5*(pic.real**2 + pic.imag**2)))
        #datafft[:,0] = np.nan    
        picfft = datafft*0
        picfft[:, 0:nlon/2] = datafft[:, nlon/2:]
        picfft[:, nlon/2:] = datafft[:, 0:nlon/2]

        plt.contourf(ks, range(len(pictime)), picfft, 41, cmap=cm.jet)
        #plt.contourf(lon, range(len(pictime)), picdata[vars[i]] * facs[i], clevs[i], 
        #        norm = mc.BoundaryNorm(clevs[i], 256), cmap=cms[i])
        plt.colorbar(orientation='vertical')
        plt.plot([0,0], [0,300], 'k--', lw=2)
        #plt.contour(lon, range(len(pictime)), picdata[vars[i]], [0], colors='m')
        plt.title(tstr)
        plt.xticks(range(-nlon/2,nlon/2,1))
        plt.grid('on')
        #plt.axis([-20, 20, 0, len(pictime)])
        #plt.axis([-10, 10, 0,30])

    plt.tight_layout()
    pfig = pcases+'../figures/'+casename+'/'
    os.system('mkdir '+pfig)
    ffig = pfig+casename+'_timespec.png'; print ffig
    plt.savefig(ffig)
    plt.close(1)

