import numpy as np

from namelist import *

from netCDF4 import Dataset as netcdf

# ============================================================================
# Physical parameterization schemes
# ============================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mc
import warnings;  warnings.filterwarnings('ignore')
import os
#import sys
#sys.path.append('/disk2/yhy/Work/tests/pymods')
#from mod_std_fplot import fplotcoast


def forcing(atm, coord, dtime, runtime, modeltime):

    if Rayleigh_friction:
        for k in range(nlev):
            sigma = (atm.p.g[k]/100.0 - pbltop) / (1000.0-pbltop)  # 2D
            sigma[sigma<0] = 0
            atm.dudt_phy.g[k] -= taufriction * sigma * atm.u.g[k]
            atm.dvdt_phy.g[k] -= taufriction * sigma * atm.v.g[k]
            
    if Newtonian_cooling:
        atm.dTdt_phy.g -= taucooling * (atm.T.g - atm.Te.g)

        #for k in range(nlev):
            #tempeq[k,:,:] = np.exp(np.log(288.0) + 287.0/1004/2 * np.log(coord.hyam[k] + coord.hybm[k]))
            #atm.dTdt_phy.g[k] -= taucooling * (atm.T.g[k]-tempeq[k])
    
    # Held-Suarez forcing
    if artforcing == 1:
        sigmab = 0.7;  kf = 1.0/86400; ka = 0.025/86400; ks = 0.25/86400
        dT = 60; dtheta = 10;
        for k in range(nlev):
            ps = np.exp(atm.lnps.g)
            p  = coord.hyam[k] * P0 + coord.hybm[k] * ps      # 2D
            sigma = (p/ps - sigmab) / (1-sigmab)  # 2D
            sigma[sigma<0] = 0             
            kv  = kf * sigma       # 2D
            atm.dudt_phy.g[k] -= kv * atm.u.g[k]
            atm.dvdt_phy.g[k] -= kv * atm.v.g[k]

            for j in range(nlat):
                kt  = ka + (ks-ka)*sigma[j,:]*(coord.coslat[j]**4)    # 1D
                Teq = (p[j,:]/P0)**(Rd/Cpd) * (315 - dT*(coord.sinlat[j]**2)-\
                    dtheta*np.log(p[j,:]/P0)*(coord.coslat[j]**2)) # 1D
                Teq[Teq<200] = 200
                atm.dTdt_phy.g[k,j,:] -= kt * (atm.T.g[k,j,:] - Teq)

    # elevated heating induced tropical waves, e.g. Matsuno-Gill modes
    if artforcing == 2:
        # heating center location
        latc = 0.0/180 * np.pi; lonc=np.pi; sigmac=0.75  
        R = 1.5e6;   H = 0.2    # heating size: horizontal (m); vertical (bar)
        heating = 10.0/86400    # heating rate (K s-1)
 
        r = np.zeros([nlat, nlon])
        for j in range(nlat):
            for i in range(nlon):
                r[j,i] = rad*np.arccos(np.sin(latc)*coord.sinlat[j] + \
                    np.cos(latc)*coord.coslat[j] * \
                    np.cos(coord.lon[i]*np.pi/180 - lonc))

        for k in range(nlev):
            ps = np.exp(atm.lnps.g)
            sigma  = (coord.hyam[k] * P0 + coord.hybm[k] * ps)/ps      # 2D
            atm.dTdt_phy.g[k,:,:] += heating * \
                np.exp(-(r/R)**2 - ((sigma-sigmac)/H)**2)
    
    # surface heating, e.g. SST
    if artforcing == 3:
        latc=0; lonc=np.pi;             # heating center location
        R = 1.5e6;                      # heating size: horizontal (m)
        heating = 10.0/86400            # heating rate (K s-1)
        for j in range(nlat):
            for i in range(nlon):
                r = rad*np.arccos(np.sin(latc)*coord.sinlat[j] + \
                    np.cos(latc)*coord.coslat[j]* \
                    np.cos(coord.lon[i]*np.pi/180 - lonc))
                # heating the bottom level only
                atm.T.g[-1,j,i] += 0.5*dtime* heating * np.exp(-(r/R)**2) 
    

    # oscillated heating
    if artforcing == 4:
        # heating center location and oscilation period
        latc = 0.0/180 * np.pi; lonc=np.pi; levc=500.0;  posci=44*86400.0
        dx = 4.0e6;  dy = 2.0e6;  dp = 300.0    # heating size        
        heating = 5.0/86400    # heating rate (K s-1)
 
        tmp = np.zeros([nlev, nlat, nlon]) + \
            np.sin((runtime.days*86400.0+runtime.seconds)/posci*2*np.pi)
        for k in range(nlev):
            ps = np.exp(atm.lnps.g)
            p  = coord.hyam[k] * P0 + coord.hybm[k] * ps      # 2D, Pa
            tmp[k,:,:] *= np.exp(-((p/100.0 - levc)/dp)**2)
##        for i in range(nlon):
##            if np.abs(coord.lon[i]*np.pi/180 - lonc) <= np.pi/kx:
##                tmp[:,:,i] *= np.sin((coord.lon[i]*np.pi/180 - lonc)*kx)
##            else:
##                tmp[:,:,i] *= 0.0
        for i in range(nlon):
            r = rad * np.abs(coord.lon[i]*np.pi/180 - lonc)
            tmp[:,:,i] *= np.exp(-(r/dx)**2)
        for j in range(nlat):
            r = rad * np.abs(coord.lat[j]*np.pi/180 - latc)
            tmp[:,j,:] *= np.exp(-(r/dy)**2)
        atm.dTdt_phy.g += heating * tmp


    # wave-CISK like heating
    if artforcing == 5:
        ctop = 200.0;  cbase = 850.0    # cloud(heating) top and base (hPa)
        #skew = 0 ~ 1: deep ~ shallow
        if heatingmode == 'shallow': skew = 0.8
        if heatingmode == 'deep': skew = 0.2
        heating = 5.0/86400  # heating rate (K/s)
        maxconverg = 1.0e-5  # maximum convergence (1/s)
        scalefac   = 4.0e-6  # D=1.0e-6 1/s ~ Q1=5 K/day
        ptop       = 850.0   # PBL top (hPa) 
        
        npbl = len(coord.lev[coord.lev > ptop])
        div3d = np.zeros([nlev,nlat,nlon])
        spfilter_3d(atm.div.s, div3d, coord, nlayer=nlev, kmax=3 )
        tmp2d = np.mean(-1*div3d[-npbl:,:,:], axis=0)  # mean PBL convergence
        tmp2d[tmp2d<0] = 0.0                    # mask divergence grids
        tmp2d[tmp2d>maxconverg] = maxconverg    # limit the max
        tmp2d /= scalefac                       # normalized

        # vertical: consine with skewness
        fac1 = np.cos((atm.p.g/100-(ctop+cbase)/2.0)/(cbase-ctop)*np.pi)
        fac2 = 2*(atm.p.g/100-cbase+skew*(cbase+ctop-2*atm.p.g/100)) / (ctop-cbase)
        fac1[fac1<0] = 0.0; fac2[fac2<0] = 0.0
        tmp = fac1 * fac2**3
        for k in range(nlev):
            tmp[k,:,:] *= tmp2d
        atm.dTdt_phy.g += heating * tmp


    # vertical cosine heating mimik deep/shallow mode
    if artforcing == 6:
        latc = 0.0/180*np.pi;  lonc = np.pi # heating location
        rx = 3.0e6;  ry = 1.5e6         # heating size: radius (m)
        ctop = 200.0;  cbase = 850.0    # cloud(heating) top and base (hPa)
        #skew = 0 ~ 1: deep ~ shallow
        if heatingmode == 'shallow': skew = 0.8
        if heatingmode == 'deep': skew = 0.2
        heating = 5.0/86400 # heating rate (K/s)

        # time oscillation
        tmp = np.zeros([nlev, nlat, nlon]) + \
            np.cos((runtime.days*86400.0+runtime.seconds)/heatingperiod*2*np.pi)
        # vertical: consine with skewness
        fac1 = np.cos((atm.p.g/100-(ctop+cbase)/2.0)/(cbase-ctop)*np.pi)
        fac2 = 2*(atm.p.g/100-cbase+skew*(cbase+ctop-2*atm.p.g/100)) / (ctop-cbase)
        fac1[fac1<0] = 0.0; fac2[fac2<0] = 0.0
        tmp *= fac1 * fac2**3
        # horizontal Gaussian distribution
        for i in range(nlon):
            r = rad * np.abs(coord.lon[i]*np.pi/180 - lonc)
            tmp[:,:,i] *= np.exp(-(r/rx)**2)
        for j in range(nlat):
            r = rad * np.abs(coord.lat[j]*np.pi/180 - latc)
            tmp[:,j,:] *= np.exp(-(r/ry)**2)
        atm.dTdt_phy.g += heating * tmp

    # vertical cosine heating mimik deep/shallow mode
    # x-direction: periodic cosine heating
    # y-direction: Gaussian distribution
    if artforcing == 7:
        latc = 0.0/180*np.pi;  lonc = np.pi # heating location
        rx = 3.0e6;  ry = 1.5e6         # heating size: radius (m)
        ctop = 200.0;  cbase = 850.0    # cloud(heating) top and base (hPa)
        #skew = 0 ~ 1: deep ~ shallow
        if heatingmode == 'shallow': skew = 0.8
        if heatingmode == 'deep': skew = 0.2
        heating = 5.0/86400 # heating rate (K/s)
        kx = 6

        # time oscillation
        tmp = np.zeros([nlev, nlat, nlon]) + \
            np.cos((runtime.days*86400.0+runtime.seconds)/heatingperiod*2*np.pi)
        # vertical: consine with skewness
        fac1 = np.cos((atm.p.g/100-(ctop+cbase)/2.0)/(cbase-ctop)*np.pi)
        fac2 = 2*(atm.p.g/100-cbase+skew*(cbase+ctop-2*atm.p.g/100)) / (ctop-cbase)
        fac1[fac1<0] = 0.0; fac2[fac2<0] = 0.0
        tmp *= fac1 * fac2**3
        # x: consine
        for i in range(nlon):
            tmp[:,:,i] *= np.cos(kx*(coord.lon[i] - 180.0)*np.pi/180)
        # y: Gaussian distribution
        for j in range(nlat):
            r = rad * np.abs(coord.lat[j]*np.pi/180 - latc)
            tmp[:,j,:] *= np.exp(-(r/ry)**2)
        atm.dTdt_phy.g += heating * tmp

    
    # heating from ERA interim 6hourly files
    if artforcing == 200:
        d0 = datetime(modeltime.year, modeltime.month, modeltime.day)
        h0 = (modeltime-d0).seconds//(6*3600) * 6
        era = ERA(d0 + timedelta(hours=h0))
        era.fheating(coord, atm)

##    ######################
##    if runtime.seconds == 0:
##        plt.figure(1, figsize=(12,6))
##        piclev = [200, 500, 800, 1000]
##        for i in range(4):
##            plt.subplot(2,2,i+1)
##            ilev = np.argmin(np.abs(coord.lev-piclev[i]))
##            if i == 0: 
##                plt.contourf(coord.lon, coord.lat, atm.u.g[ilev,:,:], 21, cmap=cm.jet)
##                plt.title('U @ '+format(coord.lev[ilev]))
##            if i == 1: 
##                plt.contourf(coord.lon, coord.lat, atm.dTdt_phy.g[ilev,:,:]*86400, 21, cmap=cm.jet)
##                plt.title('heating @ '+format(coord.lev[ilev]))
##            if i == 2: 
##                plt.contourf(coord.lon, coord.lat, atm.T.g[ilev,:,:], 21, cmap=cm.jet)
##                plt.title('T @ '+format(coord.lev[ilev]))
##            if i == 3: 
##                plt.contourf(coord.lon, coord.lat, np.exp(atm.lnps.g[:,:])/100, 21, cmap=cm.jet)
##                plt.title('ps')
##            plt.colorbar()
##
##        plt.tight_layout()
##        tstr = format(modeltime.year,'04')+format(modeltime.month,'02')+ \
##            format(modeltime.day,'02')+'_'+format(modeltime.hour,'02')+ \
##            '-'+format(modeltime.minute,'02')+'-'+format(modeltime.second,'02')
##        pfig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+casename+'/'
##        os.system('mkdir '+pfig)
##        ffig = pfig+casename+'_heating_'+\
##            tstr+'.png'
##        plt.savefig(ffig)
##        plt.close(1)
##
##    ######################

    # diagnostic precipitation
    atm.prec.g *= 0
    for k in range(nlev):
        atm.prec.g += precfrac * 86400 * Cpd / grav/Lref * atm.dTdt_phy.g[k] * atm.dp.g[k]
    atm.prec.g[atm.prec.g < 0] = 0.0
    
    return


# ============================================================================
#   ERA Interim 6hourly grib files
# ============================================================================
#import pygrib
class ERA:
    def __init__(self, time):
        self.time = time
        self.fp = '/disk2/yhy/Work/Data/ERA_Interim/6hourly_grib/'+ \
            format(self.time.year,'04') + format(self.time.month,'02') + '/'
        self.str  = format(self.time.year, '04') + \
            format(self.time.month, '02') + \
            format(self.time.day,'02') + \
            '_' + format(self.time.hour, '02') + '00_000'
        self.fn = self.fp + 'erainterim_'+self.str+'.grb'
        self.lon = np.linspace(0,358.5,240)
        self.lat = np.linspace(90,-90,121)
        self.lev = np.array([1,2,3,5,7,10,20,30,50,70,100,125,150,175, 
            200,225,250,300,350,400,450,500,550,600,650,700,\
            750,775,800,825,850,875,900,925,950,975,1000.0])

    def ftimeforward(self,dhour=6):
        self.__init__(self.time + timedelta(hours=dhour))
        
    def ftimeback(self,dhour=6):
        self.__init__(self.time - timedelta(hours=dhour))

    def freadvar(self, var):
        if var in ['Q1', 'OMEGA']:
            fn = self.fp + 'ERA_'+var+'_'+self.str+'.nc'
            print fn, var
            data = fncread(fn, var)
        else:
            print self.fn, var
            fid = pygrib.open(self.fn)
            if var=='Surface geopotential':
                dall = fid.select(name='Geopotential', level=0)
                data = np.array(dall[0].values)
            else:
                dall = fid.select(name=var)
                if len(dall) == 1:
                    data = np.array(dall[0].values)        
                else:
                    data = np.zeros([len(dall), len(self.lat), len(self.lon)])
                    for ilev in range(len(dall)):
                        data[ilev,:,:] = np.array(dall[ilev].values)

            fid.close()
        return data

    def fheating(self, coord, atm):
        nlon1 = len(self.lon)
        nlat1 = len(self.lat)
        nlon2 = len(coord.lon)
        nlat2 = len(coord.lat)
        
        data0 = self.freadvar('Q1')
        data0[0:12,:,:] *= 0

        data1 = np.zeros([len(self.lev), nlat2, nlon2])
        for k in range(len(self.lev)):
            temp = np.zeros([nlat1,nlon2])
            for ilat in range(nlat1):
                temp[ilat,:] = np.interp(coord.lon,self.lon,data0[k,ilat,:],period=360.0)
            for ilon in range(nlon2):
                data1[k,:,ilon] = np.interp(coord.lat,self.lat[::-1],temp[::-1,ilon])
        data2 = np.zeros([nlev, nlat2, nlon2])
        for ilat in range(nlat2):
            for ilon in range(nlon2):
                p = coord.hyam * 1e5 + coord.hybm*np.exp(atm.lnps.g[ilat,ilon])
                data2[:,ilat,ilon] = np.interp(p, self.lev*100.0, data1[:,ilat,ilon])

        atm.dTdt_phy.g = data2 + 0
        return
    

    def feratomodel(self, coord, atm):
        #import matplotlib
        #matplotlib.use('Agg')
        #import matplotlib.pyplot as plt
        #import matplotlib.cm as cm
        #import matplotlib.colors as mc
        #import warnings;  warnings.filterwarnings('ignore')
        #import sys
        #sys.path.append('/disk2/yhy/Work/tests/pymods')
        #from mod_std_fplot import fplotcoast


        nlon1 = len(self.lon)
        nlat1 = len(self.lat)
        nlon2 = len(coord.lon)
        nlat2 = len(coord.lat)

        for ivar in range(2):
            if ivar==0:  data0 = self.freadvar('Surface pressure')
            if ivar==1:  data0 = self.freadvar('Surface geopotential')
            data1 = np.zeros([nlat2, nlon2])
            temp = np.zeros([nlat1,nlon2])
            for ilat in range(nlat1):
                temp[ilat,:] = np.interp(coord.lon,self.lon,data0[ilat,:],period=360.0)
            for ilon in range(nlon2):
                data1[:,ilon] = np.interp(coord.lat, self.lat[::-1],temp[::-1,ilon])
            if ivar==0:  atm.lnps.g = np.log(data1)
            if ivar==1:  atm.gps.g  = data1 + 0

        #u = self.freadvar('U component of wind')
        #v = self.freadvar('V component of wind')
        #vor,div = self.fuvtovordiv(u, v, self.lon, self.lat)
        
        for ivar in range(3):
            if ivar==0:  data0 = self.freadvar('U component of wind')
            if ivar==1:  data0 = self.freadvar('V component of wind')
            if ivar==2:  data0 = self.freadvar('Temperature')
            #if ivar==3:
            #    fn = '/disk2/yhy/Work/Data/ERA_Interim/6hourly_deriv/'+format(self.year*100+self.mon,'06')+'/ERA_'+var+'_'+self.str+'.nc'
            #    data0 = fncread(fn, 'Q1')
            
            #data1 = np.zeros([nlev, nlat1, nlon1])
            #ps = np.exp(atm.lnps.g)
            #for ilat in range(nlat1):
            #    for ilon in range(nlon1):
            #        p = coord.hyam * 1e5 + coord.hybm*ps[ilat,ilon]
            #        data1[:,ilat,ilon] = np.interp(p, self.lev*100.0, data0[:,ilat,ilon])            
            #data2 = np.zeros([nlev, nlat2, nlon2])
            #for k in range(nlev):
            #    temp = np.zeros([nlat1,nlon2])
            #    for ilat in range(nlat1):
            #        temp[ilat,:] = np.interp(coord.lon,self.lon,data1[k,ilat,:],period=360.0)
            #    for ilon in range(nlon2):
            #        data2[k,:,ilon] = np.interp(coord.lat,self.lat[::-1],temp[::-1,ilon])


            data1 = np.zeros([len(self.lev), nlat2, nlon2])
            for k in range(len(self.lev)):
                temp = np.zeros([nlat1,nlon2])
                for ilat in range(nlat1):
                    temp[ilat,:] = np.interp(coord.lon,self.lon,data0[k,ilat,:],period=360.0)
                for ilon in range(nlon2):
                    data1[k,:,ilon] = np.interp(coord.lat,self.lat[::-1],temp[::-1,ilon])
            data2 = np.zeros([nlev, nlat2, nlon2])
            for ilat in range(nlat2):
                for ilon in range(nlon2):
                    p = coord.hyam * 1e5 + coord.hybm*np.exp(atm.lnps.g[ilat,ilon])
                    data2[:,ilat,ilon] = np.interp(p, self.lev*100.0, data1[:,ilat,ilon])

            if ivar==0:  atm.u.g = data2 + 0
            if ivar==1:  atm.v.g = data2 + 0
            if ivar==2:  atm.T.g = data2 + 0
            #if ivar==3:  atm.Q1.g = data2 + 0
        
        atm.vor.g, atm.div.g = self.fuvtovordiv(atm.u.g, atm.v.g, coord.lon, coord.lat)


##        plt.figure(1, figsize=(12,6))
##        plt.subplot(2, 2, 1)
##        plt.contourf(coord.lon, coord.lat, atm.vor.g[-1,:,:], 21, cmap=cm.jet)
##        plt.colorbar()
##        plt.subplot(2, 2, 2)
##        plt.contourf(coord.lon, coord.lat, atm.div.g[-1,:,:], 21, cmap=cm.jet)
##        plt.colorbar()
##        fplotcoast(coord.lon)
##
##        plt.tight_layout()
##        ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+casename+'_era.png'
##        plt.savefig(ffig)
##        plt.close(1)

        return

    def fuvtovordiv(self, u,v,lon,lat):
        clat = u*0
        for ilat in range(len(lat)):
            clat[:,ilat,:] = np.cos(lat[ilat]*np.pi/180.0)
        dudx,dudy = self.fdxdy(u*clat,lon,lat)
        dvdx,dvdy = self.fdxdy(v,lon,lat)
        vor = dvdx - dudy/clat
        vor[:,0,:] *= 0
        vor[:,-1,:] *= 0

        dudx,dudy = self.fdxdy(u,lon,lat)
        dvdx,dvdy = self.fdxdy(v*clat,lon,lat)
        div = dudx + dvdy/clat
        div[:,0,:] *= 0
        div[:,-1,:] *= 0

        return vor, div

    def fdxdy(self, ps,lon,lat):
        a = 6371000.0
        dpdx = ps*0; dpdy = ps*0    
        for ilon in range(1,len(lon)-1):
            dpdx[:,:,ilon] = (ps[:,:,ilon+1]-ps[:,:,ilon-1]) / ((lon[ilon+1]-lon[ilon-1])*np.pi/180)
        dpdx[:,:,0]  = (ps[:,:,1]-ps[:,:,-1]) / ((lon[1]-lon[-1]+360)*np.pi/180)
        dpdx[:,:,-1] = (ps[:,:,0]-ps[:,:,-2]) / ((lon[0]-lon[-2]+360)*np.pi/180)
        for ilat in range(1,len(lat)-1):
            dpdx[:,ilat,:] /= a*np.cos(lat[ilat]*np.pi/180)
            dpdy[:,ilat,:] = (ps[:,ilat+1,:]-ps[:,ilat-1,:])/((lat[ilat+1]-lat[ilat-1])*np.pi/180) / a
        dpdy[:,0,:]  = (ps[:,1,:]-ps[:,0,:])/((lat[1]-lat[0])*np.pi/180) / a
        dpdy[:,-1,:] = (ps[:,-1,:]-ps[:,-2,:])/((lat[-1]-lat[-2])*np.pi/180) / a
        return dpdx, dpdy


def fncread(fn,var):
    fid = netcdf(fn, 'r')
    data = fid.variables[var][:]
    fid.close()
    return data

# ============================================================================
# spectral filtering: 3 dimensions: spectra -> grid point
# ============================================================================
from scipy.fftpack import fft, ifft
def spfilter_3d(ss, gg, coord, nlayer=nlev, kmax=3):
    ff = np.zeros([nlayer, nlat,nlon]) + 0j
    for m in range(kmax+1):   ff[:, :, m] = np.dot(ss[0:nlayer, m,m:kmax+1], coord.lp[m,m:kmax+1,:])
    for m in range(1,kmax+1): ff[:, :,-m] = np.conj(ff[:, :,m])
    gg[:,:,:] = ifft(ff*nlon, axis=2).real
    return



