from netCDF4 import Dataset as netcdf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
matplotlib.rcParams.update({'font.size':16})
import warnings;  warnings.filterwarnings('ignore')
import os

from scipy.fftpack import fft2 




def fplotwk(he = 50, n=1):
    beta = 2.289e-11; g=9.8; R=6371000.0
    ks = np.linspace(-20,20,401)/R
    p = -ks*ks*g*he - (2*n+1)*beta*np.sqrt(g*he)
    q = -beta*ks*g*he
    d = (q/2)**2 + (p/3)**3 

    w1 = d*0
    w2 = d*0
    w3 = d*0

    for i in range(len(d)):
        if d[i]>0:
            w1[i] = (-q[i]/2+np.sqrt(d[i]))**(1.0/3) + (-q[i]/2-np.sqrt(d[i]))**(1.0/3)
            omg = (-1.0 + np.sqrt(3.0) * 1.0j )/2
            w2[i] = omg*(-q[i]/2+np.sqrt(d[i]))**(1.0/3) + \
                omg*omg*(-q[i]/2-np.sqrt(d[i]))**(1.0/3)
            w3[i] = omg*omg*(-q[i]/2+np.sqrt(d[i]))**(1.0/3) + \
                omg*(-q[i]/2-np.sqrt(d[i]))**(1.0/3)
        else:
            w1[i] = (-q[i]/2+np.sqrt(-d[i])*1j)**(1.0/3) + (-q[i]/2-np.sqrt(-d[i])*1j)**(1.0/3)
            omg = (-1.0 + np.sqrt(3.0) * 1.0j )/2
            w2[i] = omg*(-q[i]/2+np.sqrt(-d[i])*1j)**(1.0/3) + \
                omg*omg*(-q[i]/2-np.sqrt(-d[i])*1j)**(1.0/3)
            w3[i] = omg*omg*(-q[i]/2+np.sqrt(-d[i])*1j)**(1.0/3) + \
                omg*(-q[i]/2-np.sqrt(-d[i])*1j)**(1.0/3)

    fac = 90*86400.0/np.pi
    if n<0: 
        ks = np.linspace(0,20,21)
        ws = ks/R*np.sqrt(g*he)
        plt.plot(ks, ws*fac, 'b-', lw=2)
    else:
        plt.plot(ks[ks<0]*R, w3.real[ks<0]*fac, 'g', lw=2)
        plt.plot(ks*R, w1.real*fac, 'r-', lw=2)
        
    return




#import glob
#==================================================================
#for case in ['fp_T21L8_s303','fp_T21L8_s304']:
#for case in ['fp_T30L18_s_dp','fp_T30L18_s_sh']:
#for case in ['fp_T30L18_d_dp','fp_T30L18_d_sh']:
#for case in ['fp_T30L18_o_dp','fp_T30L18_o_sh']:
#for case in ['fp_T30L18_p_dp','fp_T30L18_p_sh_ptop850']:
#for case in ['fp_T30L18_d_spec_sh']:

cases = ['fp_T30L18_o_md_nophy','fp_T30L18_o_md_p_dp_k3_rf5_nc8_maxconv2','fp_T30L18_o_md_p_sh_k3_rf5_nc8_maxconv2']

tstrs = ['a) no_phy','b) CISK_deep', 'c) CISK_shallow']

figid = plt.figure(1, figsize=(18,7))

for icase in range(len(cases)):

    case = cases[icase]
    pcases = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'
    #pcases = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'
    import sys;  sys.path.append(pcases+case+'/')
    from namelist import *
    import commands

    casename = case
    vars = ['U', 'PS']
    piclevs = [1000, -1]

    picdata = {var:[] for var in vars}
    pictime = []
    reallev = {}

    foutfreq = -6
    if timeend*foutfreq < 0: ntime = (-timeend*tunit//timestep)//foutfreq
    else: ntime = timeend//foutfreq
    for itime in range(40, ntime):
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
                #data[0,:] -= np.mean(data[0,:])
                if len(pictime) == 1:
                    picdata[vars[ivar]] = data + 0
                    print time
                else:
                    picdata[vars[ivar]] = np.append(picdata[vars[ivar]], data, axis=0)

    #=====================================================================
    nlon = len(lon)
    print time, picdata['PS'].shape

    #for i in range(picdata['PS'].shape[0]):
    #    picdata['PS'][i,:] -= np.mean(picdata['PS'][i,:])
    #for i in range(picdata['PS'].shape[1]):
    #    picdata['PS'][:,i] -= np.mean(picdata['PS'][:,i])
    pic = fft2(picdata['PS'])
    tstr = 'PS'
    clev = np.linspace(4,14,21)

    plt.subplot(1,3,icase+1)
    ks = -1*np.linspace(-nlon/2,nlon/2-1,nlon)
    ws = np.linspace(0,len(pictime)/2-1,len(pictime)/2)
    datafft = np.log(np.sqrt(0.5*(pic.real**2 + pic.imag**2)))
    #datafft[:,0] = np.nan    
    picfft = datafft*0
    picfft[:, 0:nlon/2] = datafft[:, nlon/2:]
    picfft[:, nlon/2:] = datafft[:, 0:nlon/2]

    #clev = np.linspace(0,10,41)
    #plt.contourf(ks, ws, picfft[:len(pictime)/2,:], 41, cmap=cm.jet)
    plt.contourf(ks, ws, picfft[:len(pictime)/2,:], clev, 
            norm = mc.BoundaryNorm(clev, 256), cmap=cm.hot_r)
    #plt.colorbar(orientation='vertical')
    plt.plot([0,0], [0,len(pictime)/2], 'k--', lw=2)
    
    for he in [5,20,50]:
        fplotwk(he = he, n=1)
        fplotwk(he = he, n=-1)
    #plt.contour(lon, range(len(pictime)), picdata[vars[i]], [0], colors='m')
    if icase ==0:
        plt.ylabel('Period')
    plt.xlabel('[westward]    wave-number    [eastward]')
    #--------------------------
    
    #--------------------------

    plt.title(tstrs[icase], loc='left')
    plt.xticks(range(-16,16,4))
    nday = 90.0; ts=[90,20,10,5,4,3,2,1]
    plt.yticks([nday/t for t in ts], 
            [format(t)+'d' for t in ts])
    plt.grid('on')
    plt.axis([-15, 15, 0, 90.0/0.85])
    #plt.axis([-10, 10, 0,30])

plt.tight_layout()
pfig = pcases+'../figures/'
os.system('mkdir '+pfig)
ffig = pfig+'mar565_fig4_freqwavenum.png'; print ffig
plt.savefig(ffig)
plt.close(1)


