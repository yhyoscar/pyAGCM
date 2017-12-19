import numpy as np
import pygrib
import sys
sys.path.append('/disk2/yhy/Work/tests/pymods')
from mod_std_read import *
from mod_std_fdiag import *
from mod_std_fplot import *
from mod_std_OOP import *

# ============================================================================
#   ERA Interim 6hourly grib files
# ============================================================================
class ERA:
    def __init__(self, time):
        self.time = time
        self.year = time//1000000
        self.mon  = np.mod(time,1000000)//10000
        self.day  = np.mod(time,10000)//100
        self.hour = np.mod(time,100)
        self.str  = format(self.year*10000+self.mon*100+self.day,'08')+'_'+format(self.hour,'02')+'00_000'
        if np.mod(self.year,4)==0:        
            self.ndays = [31,29,31,30,31,30,31,31,30,31,30,31]
        else:
            self.ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
        if self.year<1979:
            self.__init__(1979010100)
        self.fn = '/disk2/yhy/Work/Data/ERA_Interim/6hourly_grib/'+format(self.year*100+self.mon)+'/erainterim_'+self.str+'.grb'
        self.lon = np.linspace(0,358.5,240)
        self.lat = np.linspace(90,-90,121)
        self.lev = np.array([1,2,3,5,7,10,20,30,50,70,100,125,150,175, 
            200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000.0])

    def ftimeforward(self,dhour=6):
        self.hour += dhour
        if self.hour >= 24:           self.day  += 1;  self.hour = 0
        if self.day  > self.ndays[self.mon-1]:    self.mon  += 1;  self.day  = 1
        if self.mon  > 12:           self.year += 1;  self.mon  = 1
        self.__init__(self.year*1000000+self.mon*10000+self.day*100+self.hour)
        
    def ftimeback(self,dhour=6):
        self.hour -= dhour
        if self.hour < 0:   self.day -= 1;  self.hour=24-dhour
        if self.day  < 1:   
            self.mon -= 1
            if self.mon < 1:
                self.year -= 1; self.mon = 12; self.day = 31
            else:
                self.day = self.ndays[self.mon-1]
        self.__init__(self.year*1000000+self.mon*10000+self.day*100+self.hour)

    def freadvar(self, var):
        if var in ['Q1', 'OMEGA']:
            fn = '/disk2/yhy/Work/Data/ERA_Interim/6hourly_deriv/'+format(self.year*100+self.mon,'06')+'/ERA_'+var+'_'+self.str+'.nc'
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

    def feratomodel(self, coord, atm):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mc
        import warnings;  warnings.filterwarnings('ignore')
        import sys
        sys.path.append('/disk2/yhy/Work/tests/pymods')
        from mod_std_fplot import fplotcoast


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
        
        atm.vor.g, atm.div.g = self.fuvtovordiv(atm.u.g, atm.v.g, coord.lon, coord.lat)

        plt.figure(1, figsize=(12,6))
        plt.subplot(2, 2, 1)
        plt.contourf(coord.lon, coord.lat, atm.vor.g[-1,:,:], 21, cmap=cm.jet)
        plt.colorbar()
        plt.subplot(2, 2, 2)
        plt.contourf(coord.lon, coord.lat, atm.div.g[-1,:,:], 21, cmap=cm.jet)
        plt.colorbar()
        fplotcoast(coord.lon)

        plt.tight_layout()
        ffig = '/disk2/yhy/ModelingGroup/XYZ_GCM/figures/'+casename+'_era.png'
        plt.savefig(ffig)
        plt.close(1)
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


class ModelFile:
    def __init__(self, casename, runtime):
        self.time = runtime
        self.case = casename
        iday  = runtime//86400
        ihour = np.mod(runtime,86400)//3600
        imin  = np.mod(runtime,3600)//60
        isec  = np.mod(runtime,60)
        tstr  = 'd'+format(iday,'04')+'_'+format(ihour,'02')+'-'+format(imin,'02')+'-'+format(isec,'02')
        self.fn = '/disk2/yhy/ModelingGroup/XYZ_cases/'+casename+'/'+casename+'_'+tstr+'.nc'
        self.lon = fncread(self.fn, 'lon')
        self.lat = fncread(self.fn, 'lat')
        self.lev = fncread(self.fn, 'lev')
        self.hyam = fncread(self.fn, 'hyam')
        self.hybm = fncread(self.fn, 'hybm')

    
    def freadvar(self, var, plev=-1):
        print self.fn, var

        if plev<0:
            data = fncread(self.fn, var)[0]
        else:
            data = np.zeros([len(self.lat), len(self.lon)]) + np.nan
            tmp = fncread(self.fn, var)[0]    
            self.ps  = fncread(self.fn, 'PS')[0]
            for ilev in range(len(self.lev)-1):
                flag1 = np.zeros([len(self.lat), len(self.lon)]).astype(int)
                flag2 = np.zeros([len(self.lat), len(self.lon)]).astype(int)
                p1 = self.hyam[ilev]*1e5 + self.hybm[ilev]*self.ps
                p2 = self.hyam[ilev+1]*1e5 + self.hybm[ilev+1]*self.ps
                flag1[p1 <= plev*100] = 1
                flag2[p2 >= plev*100] = 1
                data[flag1*flag2==1] = tmp[ilev][flag1*flag2==1] + (tmp[ilev+1][flag1*flag2==1] - tmp[ilev][flag1*flag2==1]) / (p2[flag1*flag2==1]-p1[flag1*flag2==1]) * (plev*100-p1[flag1*flag2==1])
        return data




