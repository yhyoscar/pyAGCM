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

from datetime import datetime, timedelta

import commands
#import glob
#==================================================================
#for case in ['fp_T21L8_s301','fp_T21L8_s302','fp_T21L8_s303', 
#    'fp_T21L8_s304','fp_T21L8_s305','fp_T21L8_s306']:

#for case in [ 'fp_T30L18_s_sh', 'fp_T30L18_s_dp']:
#for case in [ 'fp_T30L18_o_sh', 'fp_T30L18_o_dp']:
for case in [ 'fp_T30L18_p_sh_ptop850', 'fp_T30L18_p_dp']:

    #pcases = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'
    pcases = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'
    vars   = ['T', 'U', 'OMEGA', 'PS']
    #picday = [200,300]
    picday = [25, 25]
    inittime   = datetime(1, 1, 1)
    nlat = 48

    picdata = {var:[] for var in vars}
    pictime = []
    nday = 0
    for iday in range(picday[0], picday[-1]+1):
        time = inittime + timedelta(days=iday)
        tstr = format(time.year,'04')+format(time.month,'02')+ \
                format(time.day,'02')+'_'+format(time.hour,'02')+ \
                '-'+format(time.minute,'02')+'-'+format(time.second,'02')
        fn = pcases+case+'/'+case+'_'+tstr+'.nc'
        status, fn  = commands.getstatusoutput('ls '+fn)
        if status == 0:
            pictime += [time]
            print fn;  nday+=1
            fid = netcdf(fn, 'r')
            lev = fid.variables['lev'][:]
            lon = fid.variables['lon'][:]
            hyai = fid.variables['hyai'][:]
            hybi = fid.variables['hybi'][:]
            for ivar in range(len(vars)):
#                if piclevs[ivar]<0:
#                    data = np.mean(fid.variables[vars[ivar]][:,(nlat/2-1):(nlat/2+1),:], axis=1)
#                    reallev[vars[ivar]] = -1
#                else:
#                    ilev = np.argmin(np.abs(lev - piclevs[ivar]))
#                    reallev[vars[ivar]] = lev[ilev]
#                    data = np.mean(fid.variables[vars[ivar]][:,ilev,(nlat/2-1):(nlat/2+1),:], axis=1)
                if len(pictime) == 1:
                    picdata[vars[ivar]] = fid.variables[vars[ivar]][:] + 0
                else:
                    picdata[vars[ivar]] += fid.variables[vars[ivar]][:]
    print 'days = ',nday
    for ivar in range(len(vars)):
        picdata[vars[ivar]] /= nday
    picdata['dp'] = picdata['U'] * 0
    picdata['ua'] = picdata['U'] + 0
    picdata['mflx'] = picdata['U'] * 0

    um = np.mean(picdata['ua'], axis=3)
    for ilon in range(len(lon)):
        picdata['ua'][:,:,:,ilon] -= um 
    for k in range(len(lev)):
        picdata['dp'][:,k,:,:] = (hyai[k+1]-hyai[k])*1e5 + (hybi[k+1]-hybi[k])*picdata['PS']
    for k in range(len(lev)):
        picdata['mflx'][:,k,:,:] = np.sum(picdata['ua'][:,k:,:,:]*picdata['dp'][:,k:,:,:], 
            axis=1)/9.8

    Tsurf = 288.15  # surface temperature (K)
    Tstra = 216.65  # stratospheric temperature (K)
    dTdz  = 6.5     # lapse rate (K/km) -> tropopause at 226hPa
    Te = lev*0
    for k in range(len(lev)):
        Te[k] = max( Tstra, \
            np.exp(np.log(Tsurf) + (dTdz*0.001*287/9.8) * np.log(lev[k]/1000.0)) )
    #=====================================================================
    plt.figure(1, figsize=(8,8))
    plt.subplot(2,1,1)
    clev = np.linspace(-0.03,0.02,21) 
    pic = (picdata['OMEGA'][0,:,nlat/2-1,:]+picdata['OMEGA'][0,:,nlat/2,:])/2
    #plt.contourf(lon, lev, pic, 41, cmap=cm.jet)
    plt.contourf(lon, lev, pic, clev, 
        norm = mc.BoundaryNorm(clev, 256), cmap=cm.RdYlBu_r)
    plt.colorbar()
    plt.contour(lon, lev, pic, [0], colors='m')
    pic = (picdata['mflx'][0,:,nlat/2-1,:]+picdata['mflx'][0,:,nlat/2,:])/2
    clev = np.linspace(-1,1,41)*5e4
    plt.contour(lon, lev, pic, clev, colors='k', 
        norm = mc.BoundaryNorm(clev, 256))
    plt.xticks(range(0,360,30))
    
    dpicx = 3
    dpicy = 1
    picx = (picdata['U'][0,:,nlat/2-1,:]+picdata['U'][0,:,nlat/2,:])/2
    for k in range(len(lev)):
        picx[k,:] -= np.mean(picx[k,:])
    picy = -300*(picdata['OMEGA'][0,:,nlat/2-1,:]+picdata['OMEGA'][0,:,nlat/2,:])/2
    wind  = plt.quiver( lon[::dpicx],lev[::dpicy], 
        picx[::dpicy,::dpicx], picy[::dpicy,::dpicx], 
        color=[0.3,1,1],scale=1/0.008)
    #qk = plt.quiverkey(wind, 0.92, 1.02, 20, r'20 m/s', 
    #    fontproperties={'size':ft})

    plt.title('zonal mass flux, OMEGA')
    plt.gca().invert_yaxis()

    plt.subplot(2,1,2)
    pic = (picdata['T'][0,:,nlat/2-1,:]+picdata['T'][0,:,nlat/2,:])/2
    for k in range(len(lev)):
        pic[k,:] -= Te[k]
    clev = np.linspace(-5,11,33)
    plt.contourf(lon, lev, pic, 41, cmap=cm.jet)
    #plt.contourf(lon, lev, pic, clev, 
    #    norm = mc.BoundaryNorm(clev, 256), cmap=cm.jet)
    plt.colorbar()
    plt.contour(lon, lev, pic, [0], colors='k')
    plt.title('Temperature change (K)')
    plt.xticks(range(0,360,30))
    plt.gca().invert_yaxis()
    
    plt.tight_layout()

    pfig = pcases+'../figures/'+case+'/';  os.system('mkdir '+pfig)
    ffig = pfig+case+'_levlon_d'+format(picday[0])+'_'+format(picday[-1])+'.png'
    print ffig;  plt.savefig(ffig);  plt.close(1)

    
    
    
