import numpy as np
from scipy.fftpack import fft, ifft    
from os import system

from namelist import *
from Physics import *

from time import time as timer
from netCDF4 import Dataset as netcdf


# ============================================================================
#   2d Fields: grid points, spherical harmonic coefficients
# ============================================================================
class Field2d:
    def __init__(self, grid=['g','s']):
        if 'g' in grid: self.g   = np.zeros([nlat,nlon])        # Grid point field [lat, lon]
        if 's' in grid: self.s   = np.zeros([M+1,N+1]) + 0j     # Spherical harmonic field [m, n]
        self.name = ''
        self.unit = ''


# ============================================================================
#   3d Fields: grid points, spherical harmonic coefficients
# ============================================================================
class Field3d:
    def __init__(self, nlayer=nlev, grid=['g','s']):
        if 'g' in grid: self.g   = np.zeros([nlayer,nlat,nlon])        # Grid point field [lev, lat, lon]
        if 's' in grid: self.s   = np.zeros([nlayer,M+1,N+1]) + 0j     # Spherical harmonic field [lev, m, n]
        self.name = ''
        self.unit = ''


# ============================================================================
#   Coordinate and normalized associated Legendre polynomial
# ============================================================================
class Coordinate:
    def __init__(self):
        self.lon  = np.linspace(0.0, 360.0-360.0/nlon, nlon)
        x = np.polynomial.legendre.leggauss(nlat)  # a tuple contains: sin(lat) and weights
        self.lat    = np.arcsin(x[0])*180.0/np.pi
        self.sinlat = x[0] + 0
        self.coslat = np.cos(self.lat*np.pi/180.0)
        self.wt     = x[1] + 0
        
        self.clat2d = np.zeros([nlat,nlon])
        self.slat2d = np.zeros([nlat,nlon])
        self.f2d    = Field2d()
        for ilat in range(nlat):
            self.clat2d[ilat,:] = self.coslat[ilat]
            self.slat2d[ilat,:] = self.sinlat[ilat]
            self.f2d.g[ilat,:] = 2*omg*self.sinlat[ilat]
        self.f2d.s[0,1] = omg/np.sqrt(0.375)

        self.ns   = np.arange(N+1)*np.arange(1,N+2)/rad/rad  # n*(n+1)/a/a
        self.wt2      = np.zeros([nlev, nlat])
        self.coslat2  = np.zeros([nlev, nlat])
        self.ns2      = np.zeros([nlev, N+1])
        for k in range(nlev):
            self.ns2[k] = self.ns + 0
            self.wt2[k] = self.wt + 0
            self.coslat2[k]  = self.coslat*self.coslat 
      

        self.lp   = np.zeros([M+1,N+1,nlat])
        self.hp   = np.zeros([M+1,N+1,nlat])
        self.lp2  = np.zeros([M+1,nlat,N+1])
        self.hp2  = np.zeros([M+1,nlat,N+1])
        #self.lpwt2  = np.zeros([M+1,nlat,N+1])
        #self.hpwt2  = np.zeros([M+1,nlat,N+1])

        for ilat in range(nlat):
            self.lp[:,:,ilat], self.hp[:,:,ilat] = self.nalf(self.sinlat[ilat],M,N)
            self.lp2[:,ilat,:] = self.lp[:,:,ilat] + 0
            self.hp2[:,ilat,:] = self.hp[:,:,ilat] + 0
            #self.lpwt2[:,ilat,:] = self.lp[:,:,ilat]*self.wt[ilat] + 0
            #self.hpwt2[:,ilat,:] = self.hp[:,:,ilat]*self.wt[ilat] + 0

        if nlev == 26:
            self.hyai = np.array([ 0.002194067, 0.004895209, 0.009882418, 0.01805201, 0.02983724, 
                    0.04462334, 0.06160587, 0.07851243, 0.07731271, 0.07590131, 
                    0.07424086, 0.07228744, 0.06998933, 0.06728574, 0.06410509, 
                    0.06036322, 0.05596111, 0.05078225, 0.04468960, 0.03752191, 
                    0.02908949, 0.02084739, 0.01334443, 0.00708499, 0.00252136, 
                    0, 0])
            self.hybi = np.array([0, 0, 0, 0, 0, 
                    0, 0, 0, 0.01505309, 0.03276228, 
                    0.05359622, 0.07810627, 0.1069411, 0.1408637, 0.180772, 
                    0.2277220, 0.2829562, 0.3479364, 0.4243822, 0.5143168, 
                    0.6201202, 0.7235355, 0.8176768, 0.8962153, 0.9534761, 
                    0.9851122, 1])
        if nlev == 18:
            self.hyai = np.array([0.00251499, 0.00710361, 0.01904260, 0.04607560, 0.08181860, 0.07869805, 
                    0.07463175, 0.06955308, 0.06339061, 0.05621774, 0.04815296, 0.03949230, 
                    0.03058456, 0.02193336, 0.01403670, 0.007458598, 0.002646866, 0., 0.])
            self.hybi = np.array([0., 0., 0., 0., 0., 0.03756984, 0.08652625, 0.1476709, 0.221864, 
                    0.308222, 0.4053179, 0.509588, 0.6168328, 0.7209891, 0.816061, 
                    0.8952581, 0.953189, 0.985056, 1.])
        if nlev == 8:
            self.hyai = np.zeros([nlev+1])
            self.hybi = np.array([0, 0.05, 0.14, 0.26, 0.42, 0.6, 0.77, 0.9, 1.0 ])
        
        self.dhybi = self.hybi[1:] - self.hybi[:-1]
        self.hyam  = (self.hyai[:-1] + self.hyai[1:])/2
        self.hybm  = (self.hybi[:-1] + self.hybi[1:])/2
        self.lev   = (self.hyam + self.hybm)*1000
        self.ilev  = (self.hyai + self.hybi)*1000

    #------------------------------------------------------------------
    # Return the Normalized Associated Legendre Functions and their 1st-order derivatives
    #------------------------------------------------------------------
    def nalf(self,slat,mmax,nmax):
        lp = np.zeros([mmax+1, nmax+2])
        hp = np.zeros([mmax+1, nmax+1])
        
        lp[0,0] = 1.0/np.sqrt(2.0)
        for m in range(1,mmax+1):
            lp[m,m] = -np.sqrt(1-slat*slat)*np.sqrt((2.0*m+1)/2.0/m) * lp[m-1,m-1]
            lp[m-1,m] = np.sqrt(2.0*m+1)*slat * lp[m-1,m-1]
            hp[m,m] = -m*slat*lp[m,m]
        
        for m in range(mmax+1):
            for n in range(m+1,nmax+1):
                fac  = np.sqrt((1.0*n*n-m*m)/(4.0*n*n-1.0))
                fac1 = np.sqrt(((n+1.0)*(n+1.0)-m*m)/(4.0*(n+1.0)*(n+1.0)-1.0))
                lp[m,n+1] = (slat*lp[m,n] - fac*lp[m,n-1]) / fac1
                hp[m,n] = (2.0*n+1)*fac*lp[m,n-1] - n*slat*lp[m,n]

        return lp[:,:-1], hp


# ============================================================================
#   AGCM fields
# ============================================================================
class Atmosphere:
    def __init__(self):
        # prognostic variables
        self.vor      = Field3d()        # relative vorticity (s^-1)
        self.vor.name = 'relative vorticity'
        self.vor.unit = 's^-1'

        self.div      = Field3d()        # divergence         (s^-1)
        self.div.name = 'divergence'
        self.div.unit = 's^-1'

        self.T      = Field3d()        # temperature        (K)
        self.T.name = 'temperature'
        self.T.unit = 'K'

        self.q      = Field3d()        # Specific humidity  (kg/kg)
        self.q.name = 'specific humidity'
        self.q.unit = 'kg kg^-1'

        self.lnps      = Field2d()        # log(surface pressure)   (log(Pa))
        self.lnps.name = 'log(surface pressure)'
        self.lnps.unit = 'log(Pa)'

        self.gps      = Field2d()        # surface geopotential (topography) (m2 s-2)
        self.gps.name = 'surface geopotential'
        self.gps.unit = 'm2 s^-2'

        # diagnostic variables
        self.u      = Field3d(grid=['g'])  # zonal wind (m/s)
        self.u.name = 'zonal wind'
        self.u.unit = 'm s^-1'

        self.v      = Field3d(grid=['g'])  # meridional wind (m/s)
        self.v.name = 'meridional wind'
        self.v.unit = 'm s^-1'

        self.w      = Field3d(grid=['g'])  # omega (pa/s)
        self.w.name = 'vertical pressure velocity'
        self.w.unit = 'Pa s^-1'

        self.z      = Field3d(grid=['g'])  # geopotential height  (m)
        self.z.name = 'geopotential height'
        self.z.unit = 'm'

        self.ps      = Field2d(grid=['g'])  # surface pressure (Pa)
        self.ps.name = 'surface pressure'
        self.ps.unit = 'Pa'

        self.dlnpsx      = Field2d(grid=['g'])        # dlog(surface pressure)/dlon   (log(Pa) rad-1)
        self.dlnpsx.name = 'dlog(surface pressure)/dlon'
        self.dlnpsx.unit = 'log(Pa) rad-1'

        self.dlnpsy      = Field2d(grid=['g'])        # dlog(surface pressure)/dlat   (log(Pa) rad-1)
        self.dlnpsy.name = 'dlog(surface pressure)/dlat'
        self.dlnpsy.unit = 'log(Pa) rad-1'

        self.p      = Field3d(grid=['g'])  # pressure (Pa)
        self.p.name = 'pressure'
        self.p.unit = 'Pa'

        self.dp      = Field3d(grid=['g'])  # delta pressure (Pa)
        self.dp.name = 'delta pressure'
        self.dp.unit = 'Pa'

        self.etap      = Field3d(nlayer=nlev+1, grid=['g'])  # vertical advection term ()
        self.etap.name = 'vertical advection term'
        self.etap.unit = ''

        self.KE      = Field3d(grid=['g'])  # kinetic energy (m2/s2)
        self.KE.name = 'kinetic energy'
        self.KE.unit = 'm2/s2'

        self.Tv      = Field3d()        # virtual temperature     (K)
        self.Tv.name = 'virtual temperature'
        self.Tv.unit = 'K'

        self.advlnps      = Field3d(grid=['g'])  # horizontal advection of lnps ()
        self.advlnps.name = 'horizontal advection of lnps'
        self.advlnps.unit = 'log(Pa) s-1 * '

        # physics tendencies
        self.dudt_phy      = Field3d(grid=['g'])
        self.dudt_phy.name = 'zonal wind tendency due to physics'
        self.dudt_phy.unit = 'm s^-2'

        self.dvdt_phy      = Field3d(grid=['g'])
        self.dvdt_phy.name = 'meridional wind tendency due to physics'
        self.dvdt_phy.unit = 'm s^-2'

        self.dTdt_phy      = Field3d(grid=['g']) 
        self.dTdt_phy.name = 'temperature tendency due to physics'
        self.dTdt_phy.unit = 'K s^-1'

        self.dqdt_phy      = Field3d(grid=['g']) 
        self.dqdt_phy.name = 'moisture tendency due to physics'
        self.dqdt_phy.unit = 'kg kg^-1 s^-1'

        # horizontal diffusion tendencies
        self.dudt_hd      = Field3d(grid=['g'])
        self.dudt_hd.name = 'zonal wind tendency due to horizontal diffusion'
        self.dudt_hd.unit = 'm s^-2'

        self.dvdt_hd      = Field3d(grid=['g'])
        self.dvdt_hd.name = 'meridional wind tendency due to horizontal diffusion'
        self.dvdt_hd.unit = 'm s^-2'

        self.dTdt_hd      = Field3d(grid=['g'])
        self.dTdt_hd.name = 'temperature tendency due to horizontal diffusion'
        self.dTdt_hd.unit = 'K s^-1'

        # horizontal diffusion tendencies
        self.prec      = Field2d(grid=['g'])
        self.prec.name = 'diagnostic precipitation'
        self.prec.unit = 'mm day^-1'

        # equilibrium temperature
        self.Te      = Field3d(grid=['g'])
        self.Te.name = 'equilibrium temperature'
        self.Te.unit = 'K'
    
    #------------------------------------------------------------------
    # initializing atmosphere
    #------------------------------------------------------------------
    def finitial(self, coord, modeltime, iflag=0):

        # constant T
        if iflag==0:
            self.T.g    += 288.0
            self.lnps.g += np.log(1.0e5)

        # US standard atmosphere:
        #  constant lapse rate in troposphere; constant temperature in stratosphere
        if iflag==10:
            Tsurf = 288.15  # surface temperature (K)
            Tstra = 216.65  # stratospheric temperature (K)
            dTdz  = 6.5     # lapse rate (K/km) -> tropopause at 226hPa
            self.lnps.g  += np.log(1.0e5)
            for k in range(nlev):
                self.T.g[k,:,:] = max( Tstra, \
                    np.exp(np.log(Tsurf) + (dTdz*0.001*Rd/grav) * np.log(coord.hyam[k] + coord.hybm[k])) )


        # Jablonowski and Williamson, 2006: 
        # iflag==1: steady-state test;
        # iflag==2: baroclinic waves test
        if iflag in [1,2]:
            eta0 = 0.252; u0 = 35; ps = 1.0e5; etat = 0.2; T0 = 288.0; gama = 0.005; dT = 4.8e5        

            self.lnps.g += np.log(ps)
            eta  = (coord.hyam*P0 + coord.hybm*ps)/ps
            etav = (eta - eta0) * np.pi/2
            Tm = np.zeros([nlev])
            
            for k in range(nlev):
                if eta[k]<etat:
                    Tm[k] = T0*(eta[k]**(Rd*gama/grav)) + dT*((etat-eta[k])**5)
                else:
                    Tm[k] = T0*(eta[k]**(Rd*gama/grav))
                
                for j in range(nlat):
                    self.vor.g[k,j,:] = -4*u0/rad*(np.cos(etav[k])**1.5)*coord.sinlat[j]*coord.coslat[j]*(2.0-5.0*coord.sinlat[j]*coord.sinlat[j]) 
                    self.T.g[k,j,:] = Tm[k] + 0.75*eta[k]*np.pi*u0/Rd*np.sin(etav[k])*np.sqrt(np.cos(etav[k]))*((-2*(coord.sinlat[j]**6)*(coord.coslat[j]**2+1.0/3)+10.0/63)*2*u0*(np.cos(etav[k])**1.5)  +  (1.6*(coord.coslat[j]**3)*(coord.sinlat[j]**2+2.0/3)-np.pi/4)*rad*omg)
            
            etavs = (1-eta0)*np.pi/2
            for j in range(nlat):
                self.gps.g[j,:] = u0*(np.cos(etavs)**1.5)*((-2*(coord.sinlat[j]**6)*(coord.coslat[j]**2+1.0/3)+10.0/63)*u0*(np.cos(etavs)**1.5)  +  (1.6*(coord.coslat[j]**3)*(coord.sinlat[j]**2+2.0/3)-np.pi/4)*rad*omg)

            # add perturbation on the steady-state            
            if iflag==2:
                R = rad/10.0; up = 1.0; lonc = np.pi/9; latc = np.pi*2/9
                for j in range(nlat):
                    for i in range(nlon):
                        x = np.sin(latc)*coord.sinlat[j]+np.cos(latc)*coord.coslat[j]*np.cos(coord.lon[i]*np.pi/180-lonc)
                        r = rad*np.arccos(x)
                        self.vor.g[:,j,i] += up/rad*np.exp(-(r/R)**2) * (coord.sinlat[j]/coord.coslat[j]-2*((rad/R)**2)*np.arccos(x)*(np.sin(latc)*coord.coslat[j]-np.cos(latc)*coord.sinlat[j]*np.cos(coord.lon[i]*np.pi/180-lonc))/np.sqrt(1-x*x))
                        self.div.g[:,j,i] -= 2*up*rad/R/R*np.exp(-(r/R)**2) * np.arccos(x) * np.cos(latc)*np.sin(coord.lon[i]*np.pi/180-lonc)/np.sqrt(1-x*x)
    
        # Rossby-Haurwitz 4-waves
        if iflag==3:
            MK = 1.962e-6;   n=4;  pref = 95500.0; gama=0.0065;  T0=288.0
            for j in range(nlat):
                A = MK*(2*omg+MK)/2*(coord.coslat[j]**2) + (MK**2)/4*(coord.coslat[j]**(2*n))*((n+1)*coord.coslat[j]**2 + 2*n*n-n-2) - n*n*MK*MK/2*(coord.coslat[j]**(2*(n-1)))
                B = 2*(omg+MK)*MK/(n+1)/(n+2)*(coord.coslat[j]**n)*(n*n+2*n+2 - ((n+1)**2)*(coord.coslat[j]**2))
                C = MK*MK/4*(coord.coslat[j]**(2*n)) * ((n+1)*coord.coslat[j]**2 -n-2)
                for i in range(nlon):
                    self.vor.g[:,j,i] = 2*MK*coord.sinlat[j]-MK*coord.sinlat[j]*(coord.coslat[j]**n)*np.cos(n*coord.lon[i]*np.pi/180)*(n*n+3*n+2) 
                    gpsp = rad*rad*(A + B*np.cos(n*coord.lon[i]*np.pi/180) + C*np.cos(2*n*coord.lon[i]*np.pi/180))
                    self.lnps.g[j,i] = np.log(pref*(1+gama/grav/T0*gpsp)**(grav/gama/Rd))

            for k in range(nlev):
                p = coord.hyam[k]*P0 + coord.hybm[k]*np.exp(self.lnps.g)
                self.T.g[k] = T0*(p/pref)**(gama*Rd/grav)

        # mountain-induced Rossby wave train
        if iflag==4:
            u0 = 20.0;  T0=288;  h0=2000.0;  d=1.5e6;  psp=9.3e4; latc=np.pi/6; lonc=np.pi/2; N2=0.0182**2
            self.T.g += T0
            for j in range(nlat):
                self.vor.g[:,j,:] = 2*u0*coord.sinlat[j]/rad
                for i in range(nlon):
                    r = rad*np.arccos(np.sin(latc)*coord.sinlat[j] + np.cos(latc)*coord.coslat[j]*np.cos(coord.lon[i]*np.pi/180 - lonc))
                    self.gps.g[j,i] = grav * h0 * np.exp(-(r/d)**2)
                    self.lnps.g[j,i] = np.log(psp) + (rad*u0*(u0/rad+2*omg)*(coord.coslat[j]**2)/2 - self.gps.g[j,i])*N2/grav/grav/(Rd/Cpd)

        # start from an atmosphere in netcdf file (can be used for restart run)
        if iflag==100:
            fid = netcdf(initfile, 'r')
            self.vor.g  = fid.variables['VOR'][0]
            self.div.g  = fid.variables['DIV'][0]
            self.T.g    = fid.variables['T'][0]
            self.lnps.g = np.log(fid.variables['PS'][0])
            self.gps.g  = fid.variables['PHIS'][0]
            fid.close()

        # start from an atmosphere in ERA Interim 6hourly grib file
        if iflag==200:
            era = ERA(modeltime)
            print 'Reading ERA interim 6hourly grib data: ', era.fn
            era.feratomodel(coord,atm=self)
        
        # initialize equilibrium temperature
        if Teinit == 4:
            self.Te.g = self.T.g + 0
        if Teinit == 0:
            self.Te.g += 288.0
        if Teinit == 1:
            Tsurf = 288.15  # surface temperature (K)
            Tstra = 216.65  # stratospheric temperature (K)
            dTdz  = 6.5     # lapse rate (K/km) -> tropopause at 226hPa
            for k in range(nlev):
                self.Te.g[k,:,:] = max( Tstra, \
                    np.exp(np.log(Tsurf) + (dTdz*0.001*Rd/grav) * np.log(coord.hyam[k] + coord.hybm[k])) )

        
        return

    #------------------------------------------------------------------
    # horizontal diffusion 
    #------------------------------------------------------------------
    def fhoridiff(self, coord, dtime):
        # Courant number limiter
        dcfl = np.zeros([nlev,N+1]) + 1.0
        if kmxhdc>0:
            u = np.zeros([kmxhdc, nlat,nlon]); v=np.zeros([kmxhdc, nlat,nlon])
            fvordivstouvg_3d(self.vor.s[0:kmxhdc], self.div.s[0:kmxhdc], coord, u, v, nlayer=kmxhdc)
            for k in range(kmxhdc):
                wmax = np.max(np.sqrt(u[k]**2+v[k]**2))
                if wmax > 1e-5:
                    nc = int( rad/dtime/wmax )
                    if nc < N:
                        dcfl[k,nc+1:] = 1000.0
                        print 'Courant number limiter at lev ',k+1,', nc=',nc

        # both K2 and K4
        Kwind = np.zeros([nlev,N+1])  # for VOR, DIV
        KT    = np.zeros([nlev,N+1])  # for T only

        for k in range(kmnhd4-1):
            fac = 2.0**max(kmxhd2-k, 0)
            Kwind[k] = 1.0/(1+2*dtime*dcfl[k]*K2*fac*(coord.ns-2/rad/rad) )
            KT[k]    = 1.0/(1+2*dtime*dcfl[k]*K2*fac*coord.ns)
        for k in range(kmnhd4-1,nlev):
            Kwind[k] = 1.0/(1+2*dtime*dcfl[k]*K4*(coord.ns*coord.ns-4/(rad**4)) )
            KT[k]    = 1.0/(1+2*dtime*dcfl[k]*K4*coord.ns*coord.ns)
        Kwind[:,0] = 1.0

        for m in range(M+1):
            self.vor.s[:,m,:] *= Kwind[:,:]
            self.div.s[:,m,:] *= Kwind[:,:]
            self.T.s[:,m,:]   *= KT[:,:]

        # horizontal diffusion tendency of T,u,v
        Ttend   = np.zeros([nlev,M+1,N+1]) + 0j
        vortend = np.zeros([nlev,M+1,N+1]) + 0j
        divtend = np.zeros([nlev,M+1,N+1]) + 0j
        for m in range(M+1):
            Ttend[:,m,:] = (1-1.0/KT)/2/dtime * self.T.s[:,m,:]
            vortend[:,m,:] = (1-1.0/Kwind)/2/dtime * self.vor.s[:,m,:]
            divtend[:,m,:] = (1-1.0/Kwind)/2/dtime * self.div.s[:,m,:]
        
        fstog_3d(Ttend, self.dTdt_hd.g, coord)
        fvordivstouvg_3d(vortend, divtend, coord, self.dudt_hd.g, self.dvdt_hd.g)

        return

    #------------------------------------------------------------------
    # temperature correction: part of horizontal diffusion 
    #------------------------------------------------------------------
    def ftempcorr(self, coord, dtime):
        # frictional heating rate
        Cpx  = (1 + (Cpv/Cpd-1)*self.q.g) * Cpd
        fric = -(self.u.g*self.dudt_hd.g + self.v.g*self.dvdt_hd.g)/Cpx

        # delt**4 (lnps)
        lnpsd4 = np.zeros([nlat, nlon])
        fstog_d4(self.lnps.s, lnpsd4, coord)

        # vertical derivatives of T
        Tadv = np.zeros([nlev,nlat,nlon])      # Pa
        
        for k in range(1,nlev-1):
            Tadv[k] = self.ps.g/2/self.dp.g[k]*(coord.hybi[k+1]*(self.T.g[k+1]-self.T.g[k]) + coord.hybi[k]*(self.T.g[k]-self.T.g[k-1]))
        Tadv[0]  = self.ps.g/2/self.dp.g[0]*coord.hybi[1]*(self.T.g[1] - self.T.g[0])
        Tadv[-1] = self.ps.g/2/self.dp.g[-1]*coord.hybi[-2]*(self.T.g[-1] - self.T.g[-2])

        # correction on T
        for k in range(nlev):
            dT = 2*dtime*fric[k]
            if k>=kmnhd4-1:
                dT += 2*dtime*Tadv[k]*K4*lnpsd4
            self.T.g[k] += dT
            self.dTdt_hd.g[k] += dT/2/dtime
        return

    # ============================================================================
    # calculate grid point variables:
    # u,v,KE,ps,dlnpsx,dlnpsy,p,dp,advlnps,Tv,Z,etap,omega
    # ============================================================================
    def fgridvars(self, coord):
        t0 = timer()
        # vors, divs -> u,v,KE
        fvordivstouvg_3d(self.vor.s, self.div.s, coord, self.u.g, self.v.g)        
        self.KE.g = (self.u.g*self.u.g + self.v.g*self.v.g)/2.0

        #t1 = timer(); print 'u,v,KE: ',t1-t0; t0=t1

        # ps, dlnps/dlon, dlnps/dmu *(1-mu**2)
        self.ps.g = np.exp(self.lnps.g)     # Pa
        fstog_dlon(self.lnps.s, self.dlnpsx.g, coord)
        fstog_dlat(self.lnps.s, self.dlnpsy.g, coord)

        #t1 = timer(); print 'ps, dlnps/dlon, dlnps/dlat: ',t1-t0; t0=t1

        # p, dp
        for k in range(nlev):
            self.p.g[k]  = coord.hyam[k] * P0 + coord.hybm[k]*self.ps.g   # Pa
            self.dp.g[k] = (coord.hyai[k+1] - coord.hyai[k]) * P0 + (coord.hybi[k+1] - coord.hybi[k]) * self.ps.g   # Pa

        #t1 = timer(); print 'p, dp: ',t1-t0; t0=t1

        # horizontal advection of lnps    
        for k in range(nlev):
            self.advlnps.g[k] = (self.u.g[k]*self.dlnpsx.g + self.v.g[k]*self.dlnpsy.g)/rad/coord.clat2d

        #t1 = timer(); print 'horizontal advection of lnps: ',t1-t0; t0=t1

        # Tv, Z
        self.Tv.g  = (1 + (Rv/Rd-1)*self.q.g) * self.T.g
        for k in range(nlev):
            self.z.g[k]  = self.gps.g + Rd*self.Tv.g[k]*self.dp.g[k]/2/self.p.g[k] 
            self.z.g[k] += np.sum(Rd*self.Tv.g[k+1:]*self.dp.g[k+1:]/self.p.g[k+1:], axis=0)
        self.z.g /= grav

        #t1 = timer(); print 'Tv, Z: ',t1-t0; t0=t1
              
        # vertical advection: etadot * dp/deta
        for k in range(1,nlev):
            self.etap.g[k]  = coord.hybi[k] * np.sum(self.div.g*self.dp.g, axis=0) 
            self.etap.g[k] += coord.hybi[k]*self.ps.g * np.tensordot(self.advlnps.g, coord.dhybi, axes=([0],[0]))
            self.etap.g[k] -= np.sum(self.div.g[0:k]*self.dp.g[0:k], axis=0)
            self.etap.g[k] -= self.ps.g * np.tensordot(self.advlnps.g[0:k], coord.dhybi[0:k], axes=([0],[0]))
            
        #t1 = timer(); print 'etap: ',t1-t0; t0=t1

        # omega
        for k in range(nlev):
            self.w.g[k]  = coord.hybm[k]*self.ps.g*self.advlnps.g[k]
            self.w.g[k] -= (self.div.g[k]*self.dp.g[k] + self.ps.g*self.advlnps.g[k]*coord.dhybi[k])/2
        for k in range(1,nlev):
            self.w.g[k] -= np.sum(self.div.g[0:k]*self.dp.g[0:k], axis=0) 
            self.w.g[k] -= self.ps.g * np.tensordot( self.advlnps.g[0:k], coord.dhybi[0:k], axes=([0],[0]) )

        #t1 = timer(); print 'omega: ',t1-t0; t0=t1
        return


    #------------------------------------------------------------------
    # output netcdf file
    #------------------------------------------------------------------
    def foutput(self, coord, runtime, modeltime):
        # output variables list
        allvars = ['VOR','DIV','T','Q','PS','U','V','OMEGA','Z3','PHIS','PREC','TE']   

        tstr = format(modeltime.year,'04')+format(modeltime.month,'02')+ \
            format(modeltime.day,'02')+'_'+format(modeltime.hour,'02')+ \
            '-'+format(modeltime.minute,'02')+'-'+format(modeltime.second,'02')

        fnout = casedir+casename+'_'+tstr+'.nc'
        print 'output: ', fnout
        system('rm -f '+fnout)
        fidout = netcdf(fnout, 'w')

        fidout.createDimension('time', 1)
        fidout.createDimension('lat',  nlat)
        fidout.createDimension('lon',  nlon)
        fidout.createDimension('lev',  nlev)
        fidout.createDimension('ilev', nlev+1)

        ftime = fidout.createVariable('time', 'f', ('time',))
        ftime[:] = runtime.days*86400 + runtime.seconds
        ftime.units = 'second'
        ftime.long_name = 'seconds from initialization'

        flat = fidout.createVariable('lat', 'd', ('lat',))
        flat[:] = coord.lat + 0
        flat.units = 'degrees_north'
        flat.long_name = 'latitude'

        fwt = fidout.createVariable('wt', 'd', ('lat',))
        fwt[:] = coord.wt + 0
        fwt.units = 'degrees_north'
        fwt.long_name = 'Gaussian weights'

        flon = fidout.createVariable('lon', 'f', ('lon',))
        flon[:] = coord.lon + 0
        flon.units = 'degrees_east'
        flon.long_name = 'longitude'

        flon = fidout.createVariable('lev', 'f', ('lev',))
        flon[:] = coord.lev + 0
        flon.units = 'hPa'
        flon.long_name = 'reference pressure at mid-point layer'

        flon = fidout.createVariable('hyam', 'f', ('lev',))
        flon[:] = coord.hyam + 0
        flon.units = '#'
        flon.long_name = 'A at mid-point layer'

        flon = fidout.createVariable('hybm', 'f', ('lev',))
        flon[:] = coord.hybm + 0
        flon.units = '#'
        flon.long_name = 'B at mid-point layer'

        flon = fidout.createVariable('ilev', 'f', ('ilev',))
        flon[:] = coord.ilev + 0
        flon.units = 'hPa'
        flon.long_name = 'reference pressure at interface'

        flon = fidout.createVariable('hyai', 'f', ('ilev',))
        flon[:] = coord.hyai + 0
        flon.units = '#'
        flon.long_name = 'A at interface'

        flon = fidout.createVariable('hybi', 'f', ('ilev',))
        flon[:] = coord.hybi + 0
        flon.units = '#'
        flon.long_name = 'B at interface'

        for i in range(len(foutvars)):
            varstr = foutvars[i]
            if not (varstr in allvars): 
                print 'ERROR: foutput: ', varstr, ' is not in our output list!'
                return
            if varstr == 'VOR':   fdata = self.vor;
            if varstr == 'DIV':   fdata = self.div;
            if varstr == 'T':     fdata = self.T;
            if varstr == 'TE':    fdata = self.Te;
            if varstr == 'Q':     fdata = self.q;
            if varstr == 'PS':    fdata = self.ps;
            if varstr == 'PHIS':  fdata = self.gps;
            if varstr == 'U':     fdata = self.u;
            if varstr == 'V':     fdata = self.v;
            if varstr == 'OMEGA': fdata = self.w;
            if varstr == 'Z3':    fdata = self.z;
            if varstr == 'PREC':  fdata = self.prec;

            if len(fdata.g.shape) == 2:
                fdout = fidout.createVariable(varstr, 'f', ('time','lat','lon',))
            if len(fdata.g.shape) == 3:
                fdout = fidout.createVariable(varstr, 'f', ('time','lev','lat','lon',))

            fdout[:] = fdata.g + 0
            fdout.units = fdata.unit + ''
            fdout.long_name = fdata.name + ''
        
        fidout.close()
        return


# ============================================================================
# Reference atmosphere
# ============================================================================
class ReferenceAtmosphere:
    # variables list:
    #   psr, Tr, pr, pri, dpr, Hr, Dr, br, hr, Aninv
    def __init__(self, coord, dtime):
        self.psr = psref            # from namelist: Pa
        self.Tr  = np.array(Tref)   # from namelist: K

        # derived variables based on psr and Tr
        self.pr  = coord.hyam*P0 + coord.hybm*self.psr  # Pa
        self.pri = coord.hyai*P0 + coord.hybi*self.psr  # Pa
        self.dpr = self.pri[1:] - self.pri[:-1]  # Pa

        # Hr
        self.Hr = np.zeros([nlev,nlev])
        for l in range(nlev):
            self.Hr[l,l] = self.dpr[l]/self.pr[l]/2.0
            self.Hr[0:l,l] = self.dpr[l]/self.pr[l]

        # br, hr 
        # NOTE: Our calculation here is based on Simmons and Burridge, 1981;
        # but in CAM3: br=Tr, hr=0. The difference between our calculation and 
        # CAM3 is very small and negligible.
        self.hr = self.Tr*self.psr*coord.hybm/self.pr
        self.br = np.zeros([nlev])
        for k in range(nlev-1):
            self.br[k] = self.psr * np.sum(self.Tr[k+1:]*(coord.hybi[k+2:]/self.pri[k+2:]-coord.hybi[k+1:-1]/self.pri[k+1:-1]))
        self.br[nlev-1] = self.Tr[nlev-1] - self.hr[nlev-1]

        # Dr
        ee  = np.zeros([nlev,nlev])
        ff  = np.zeros([nlev,nlev])
        self.Dr = np.zeros([nlev,nlev])
        for k in range(nlev):
            ee[k,0:k]   = 1
            ff[k,0:k+1] = 1         

        for l in range(nlev):
            self.Dr[l:,l] = Rd*self.Tr[l:]*self.dpr[l]/self.pr[l:]/Cpd
            self.Dr[l,l] /= 2
            for k in range(1,nlev-1):
                self.Dr[k,l] += 0.5*self.dpr[l]/self.dpr[k]*( (self.Tr[k]-self.Tr[k-1])*(coord.hybi[k]-ee[k,l]) + (self.Tr[k+1]-self.Tr[k])*(coord.hybi[k+1]-ff[k,l]) )
            self.Dr[0,l]  += 0.5*self.dpr[l]/self.dpr[0]* (self.Tr[1]-self.Tr[0])*(coord.hybi[1]-ff[0,l]) 
            self.Dr[-1,l] += 0.5*self.dpr[l]/self.dpr[-1]* (self.Tr[-1]-self.Tr[-2])*(coord.hybi[-1]-ee[-1,l]) 

        # Aninv for Helmholtz equation
        self.Aninv = np.zeros([N+1,nlev,nlev])
        temp1 = np.zeros([nlev,1]);  temp1[:,0] = self.br+self.hr
        temp2 = np.zeros([1,nlev]);  temp2[0,:] = self.dpr
        temp3 = np.dot(self.Hr, self.Dr) + np.dot(temp1,temp2)/self.psr
        for n in range(1, N+1):
            An = np.identity(nlev) + dtime*dtime*n*(n+1)/rad/rad*Rd*temp3
            self.Aninv[n] = np.linalg.inv(An)


# ============================================================================
# 3D: relative vorticity, divergence in spectral space -> ucos, vcos in grid points
# ============================================================================
def fvordivstouvg_3d(vors,divs,coord, u, v, nlayer=nlev):    
    ufft  = np.zeros([nlayer,nlat,nlon]) + 0j
    vfft  = np.zeros([nlayer,nlat,nlon]) + 0j
    tmp = np.dot(np.append(vors[:,0,1:N+1]/(coord.ns2[0:nlayer,1:N+1]*rad), divs[:,0,1:N+1]/(coord.ns2[0:nlayer,1:N+1]*rad), axis=0), 
            coord.hp[0,1:N+1,:])
    ufft[:,:,0] =  tmp[0:nlayer] + 0
    vfft[:,:,0] = -tmp[nlayer:]  + 0

    for m in range(1,M+1):
        tmp = np.dot( np.append(vors[:,m,m:N+1]/(coord.ns2[0:nlayer,m:N+1]*rad), divs[:,m,m:N+1]/(coord.ns2[0:nlayer,m:N+1]*rad), axis=0), 
            np.append(coord.lp[m,m:N+1,:], coord.hp[m,m:N+1,:], axis=1) )
        ufft[:,:, m] = -tmp[nlayer:, 0:nlat]*1j*m   + tmp[0:nlayer, nlat:]
        vfft[:,:, m] = -tmp[0:nlayer, 0:nlat]*1j*m - tmp[nlayer:, nlat:]
        ufft[:,:,-m]  =  np.conj(ufft[:,:,m])
        vfft[:,:,-m]  =  np.conj(vfft[:,:,m])

    for k in range(nlayer):
        u[k] = ifft(ufft[k]*nlon, axis=1).real / coord.clat2d  # u
        v[k] = ifft(vfft[k]*nlon, axis=1).real / coord.clat2d  # v
    return


# ============================================================================
# 2 dimensions: grid point -> spectra
# ============================================================================
def fgtos_2d(gg, ss, coord):
    ff = fft(gg, axis=1)/nlon    
    for m in range(M+1): ss[m, m:N+1] = np.dot( coord.lp[m, m:N+1,:], ff[:,m]*coord.wt )
    return


# ============================================================================
# 3 dimensions: grid point -> spectra
# ============================================================================
def fgtos_3d(gg, ss, coord):
    ff = fft(gg, axis=2)/nlon
    for m in range(M+1): ss[:, m, m:N+1] = np.dot( ff[:,:,m]*coord.wt2, coord.lp2[m, :, m:N+1])
    return


# ============================================================================
# 2 dimensions: spectra -> grid point
# ============================================================================
def fstog_2d(ss, gg, coord):
    ff = np.zeros([nlat,nlon]) + 0j
    for m in range(M+1):   ff[:, m] = np.dot(ss[m,m:N+1], coord.lp[m,m:N+1,:])
    for m in range(1,M+1): ff[:,-m] = np.conj(ff[:,m])
    gg[:,:] = ifft(ff*nlon, axis=1).real
    return


# ============================================================================
# 3 dimensions: spectra -> grid point
# ============================================================================
def fstog_3d(ss, gg, coord, nlayer=nlev):
    ff = np.zeros([nlayer, nlat,nlon]) + 0j
    for m in range(M+1):   ff[:, :, m] = np.dot(ss[0:nlayer, m,m:N+1], coord.lp[m,m:N+1,:])
    for m in range(1,M+1): ff[:, :,-m] = np.conj(ff[:, :,m])
    gg[:,:,:] = ifft(ff*nlon, axis=2).real
    return

# ============================================================================
# 2 dimensions: spectra -> grid point, but dlon
# ============================================================================
def fstog_dlon(ss, gg, coord):
    ff = np.zeros([nlat,nlon]) + 0j
    for m in range(1,M+1):
        ff[:, m] = np.dot(ss[m,m:N+1]*1j*m, coord.lp[m,m:N+1,:])
        ff[:,-m] = np.conj(ff[:,m])
    gg[:,:] = ifft(ff*nlon, axis=1).real
    return


# ============================================================================
# 2 dimensions: spectra -> grid point, but dlat
# ============================================================================
def fstog_dlat(ss, gg, coord):
    ff = np.zeros([nlat,nlon]) + 0j
    for m in range(M+1):   ff[:, m] = np.dot(ss[m,m:N+1], coord.hp[m,m:N+1,:])
    for m in range(1,M+1): ff[:,-m] = np.conj(ff[:,m])
    gg[:,:] = ifft(ff*nlon, axis=1).real
    return

# ============================================================================
# 2 dimensions: spectra -> grid point, but d4
# ============================================================================
def fstog_d4(ss, gg, coord):
    ff = np.zeros([nlat,nlon]) + 0j
    for m in range(M+1):   ff[:, m] = np.dot(ss[m,m:N+1]*coord.ns[m:N+1]*coord.ns[m:N+1], coord.lp[m,m:N+1,:])
    for m in range(1,M+1): ff[:,-m] = np.conj(ff[:,m])
    gg[:,:] = ifft(ff*nlon, axis=1).real
    return


# ============================================================================
# vertical advection of a 3-dimensional variable
# ============================================================================
def fvertadv(etap, dp, x3d):
    dout = np.zeros([nlev,nlat,nlon])
    for k in range(1,nlev-1):
        dout[k] = etap[k+1]*(x3d[k+1]-x3d[k]) + etap[k]*(x3d[k]-x3d[k-1])
    dout[0]  = (x3d[1]-x3d[0])*etap[1]
    dout[nlev-1] = (x3d[nlev-1]-x3d[nlev-2])*etap[nlev-1]
    return dout/dp/2


# ============================================================================
# semi-implicit time integration
# ============================================================================
def ftimeinte(atms, flags, coord, refatm, dtime, output, runtime, modeltime):
    # atms: 2 time layers
    # flags: 1 --- current time step; 0 --- t-dt backward time
    # coord: coordinate
    # refatm: reference atmosphere
    # dtime: time step (unit: s)
    # output: flag of output 
    # runtime: run time (sec), for initial time step, output and DIV damping    
    # modeltime: model time (yyyymmdd_hh:mm:ss), for forcing and output
    # Note: The time integration starts from the grid points of t and t-dt

    
    ic = flags.index(1); ib = flags.index(0)

    t0 = timer()

    # clean physics and diffusion tendency
    atms[ic].dudt_phy.g *= 0
    atms[ic].dvdt_phy.g *= 0
    atms[ic].dTdt_phy.g *= 0
    atms[ic].dqdt_phy.g *= 0
    atms[ic].dudt_hd.g *= 0
    atms[ic].dvdt_hd.g *= 0
    atms[ic].dTdt_hd.g *= 0

    #t1 = timer(); print 'clean physics and diffusion tendency: ',t1-t0; t0=t1
    

    # calculate grid point values for the first time step
    # vor.s,div.s -> ug,vg;  lnps.s -> dlnpsx,dlnpsy
    if runtime.days*86400+runtime.seconds == 0: 
        print '+++++++++++++++++ initial time step +++++++++++++++++++'
        fgtos_3d(atms[ic].vor.g, atms[ic].vor.s, coord)
        fgtos_3d(atms[ic].div.g, atms[ic].div.s, coord)
        fgtos_2d(atms[ic].lnps.g, atms[ic].lnps.s, coord)
        atms[ic].fgridvars(coord)

    
    ###############################################################
    # physics 
    forcing(atms[ic], coord, dtime, runtime, modeltime)
    ###############################################################

    # (nU, nV)
    nu = -fvertadv(atms[ic].etap.g, atms[ic].dp.g, atms[ic].u.g)
    nv = -fvertadv(atms[ic].etap.g, atms[ic].dp.g, atms[ic].v.g)
    for k in range(nlev):
        nu[k] *= coord.clat2d
        nv[k] *= coord.clat2d
        
    for k in range(nlev):
        nu[k] += (atms[ic].vor.g[k]+coord.f2d.g)*atms[ic].v.g[k]*coord.clat2d        
        nu[k] -= Rd*atms[ic].Tv.g[k]*(coord.hybm[k]/atms[ic].p.g[k])*atms[ic].ps.g*atms[ic].dlnpsx.g/rad
        nu[k] += atms[ic].dudt_phy.g[k] * coord.clat2d         # dudt -> ducosdt

        nv[k] -= (atms[ic].vor.g[k]+coord.f2d.g)*atms[ic].u.g[k]*coord.clat2d
        nv[k] -= Rd*atms[ic].Tv.g[k]*(coord.hybm[k]/atms[ic].p.g[k])*atms[ic].ps.g*atms[ic].dlnpsy.g/rad
        nv[k] += atms[ic].dvdt_phy.g[k] * coord.clat2d         # dvdt -> dvcosdt

    #t1 = timer(); print 'nU, nV: ',t1-t0; t0=t1
    
    # D_triangular
    Dtri = 2*(atms[ic].KE.g + atms[ic].z.g*grav)
    for k in range(nlev):
        # if Tv-Tr in nu,nv, then use [ib] only
        Dtri[k] += Rd*(refatm.br[k]+refatm.hr[k])*(atms[ib].lnps.g - 2*atms[ic].lnps.g)
        # both T and T-Tr are OK, because del2(Tr)=0
        Dtri[k] += Rd* np.tensordot(refatm.Hr[k,:], (atms[ib].T.g - 2*atms[ic].T.g), axes=([0],[0]) )
    Dtri *= dtime

    #t1 = timer(); print 'D_tri: ',t1-t0; t0=t1


    # VS(n,m) and DS(n,m)
    VS   = np.zeros([nlev,M+1,N+1]) + 0j    
    DS   = np.zeros([nlev,M+1,N+1]) + 0j
    VS0  = fft(atms[ib].vor.g, axis=2)/nlon
    DS0  = fft(atms[ib].div.g, axis=2)/nlon 
    nu   = fft(2*dtime*nu, axis=2)/nlon 
    nv   = fft(2*dtime*nv, axis=2)/nlon
    temp = 1.0/rad/coord.coslat2  # [nlev, nlat]

    dtris = np.zeros([nlev, M+1, N+1]) + 0j
    fgtos_3d(Dtri, dtris, coord)

    for m in range(M+1):
        VS[:,m,m:N+1]  = np.dot( coord.wt2*(VS0[:,:,m]+nv[:,:,m]*1j*m*temp),    coord.lp2[m, :, m:N+1] )
        VS[:,m,m:N+1] += np.dot( coord.wt2*nu[:,:,m]*temp,                      coord.hp2[m, :, m:N+1])
        DS[:,m,m:N+1]  = np.dot( coord.wt2*(DS0[:,:,m]+nu[:,:,m]*1j*m*temp),    coord.lp2[m, :, m:N+1] )
        DS[:,m,m:N+1] -= np.dot( coord.wt2*nv[:,:,m]*temp,                      coord.hp2[m, :, m:N+1] )
        DS[:,m,m:N+1] += coord.ns2[:,m:N+1] * dtris[:,m,m:N+1]


    #t1 = timer(); print 'VSnm, DSnm: ',t1-t0; t0=t1
        
    # TS(n,m)
    TS   = np.zeros([nlev,M+1,N+1]) + 0j
    TS0 = 2*dtime*( atms[ic].T.g*atms[ic].div.g + \
        Rd*atms[ic].Tv.g*atms[ic].w.g/atms[ic].p.g/((1+(Cpv/Cpd-1)*atms[ic].q.g)*Cpd) + \
        atms[ic].dTdt_phy.g - fvertadv(atms[ic].etap.g, atms[ic].dp.g, atms[ic].T.g) )
    for k in range(nlev):
        TS0[k] += atms[ib].T.g[k]
        TS0[k] -= dtime* np.tensordot(refatm.Dr[k,:], (atms[ib].div.g - 2*atms[ic].div.g), axes=([0],[0]) )
    TS0 = fft(TS0, axis=2)/nlon

    ut  = np.zeros([nlev, nlat, nlon]) + 0j
    vt  = np.zeros([nlev, nlat, nlon]) + 0j
    for k in range(nlev):
        ut[k]  = fft(2*dtime*atms[ic].u.g[k]*coord.clat2d*atms[ic].T.g[k], axis=1)/nlon
        vt[k]  = fft(2*dtime*atms[ic].v.g[k]*coord.clat2d*atms[ic].T.g[k], axis=1)/nlon

    for m in range(M+1):
        TS[:,m,m:N+1]  = np.dot( coord.wt2*(TS0[:,:,m]-ut[:,:,m]*1j*m*temp),    coord.lp2[m, :, m:N+1])
        TS[:,m,m:N+1] += np.dot( coord.wt2*vt[:,:,m]*temp,                      coord.hp2[m, :, m:N+1] )

    #t1 = timer(); print 'TSnm: ',t1-t0; t0=t1


    # PS(n,m)
    PS    = np.zeros([M+1,N+1]) + 0j
    psrhs = atms[ib].lnps.g + 0
    for l in range(nlev):
        psrhs -= 2*dtime*atms[ic].div.g[l]*atms[ic].dp.g[l]/atms[ic].ps.g
        psrhs -= 2*dtime*atms[ic].advlnps.g[l]*coord.dhybi[l]
        psrhs -= dtime*refatm.dpr[l]*(atms[ib].div.g[l] - 2*atms[ic].div.g[l])/refatm.psr
    fgtos_2d(psrhs, PS, coord)

    #t1 = timer(); print 'PSnm: ',t1-t0; t0=t1


    # Helmholtz equation -> DIV(n,m) at t+dt
    for m in range(M+1):
        for n in range(m, N+1):
            rhs  = DS[:,m,n] + dtime*coord.ns[n]*Rd*(np.dot(refatm.Hr,TS[:,m,n])+(refatm.br+refatm.hr)*PS[m,n])
            atms[ib].div.s[:,m,n] = np.dot(refatm.Aninv[n,:,:], rhs)
    atms[ib].div.s[:,0,0] = 0.0 

    #t1 = timer(); print 'Helmholtz -> DIVnm at t+dt: ',t1-t0; t0=t1

    # VOR(n,m), T(n,m), lnps(n,m) at t+dt
    atms[ib].vor.s = VS + 0
    for m in range(M+1):
        for n in range(m, N+1):
            atms[ib].T.s[:,m,n]  = TS[:,m,n] - dtime*np.dot(refatm.Dr, atms[ib].div.s[:,m,n])
            atms[ib].lnps.s[m,n] = PS[m,n] - dtime*np.dot(refatm.dpr, atms[ib].div.s[:,m,n])/refatm.psr


    #t1 = timer(); print '-> VORnm, Tnm, lnpsnm at t+dt: ',t1-t0; t0=t1
    
    # output current time step variables to ncfiles
    if output:   atms[ic].foutput(coord, runtime, modeltime)

    #t1 = timer(); print 'release memories and output ncfile: ',t1-t0; t0=t1

    # 1st step of time filter: backup
    if timefilter>0:
        atms[ic].vor.g = atms[ic].vor.g + timefilter * (atms[ib].vor.g - 2*atms[ic].vor.g)
        atms[ic].div.g = atms[ic].div.g + timefilter * (atms[ib].div.g - 2*atms[ic].div.g)
        atms[ic].T.g = atms[ic].T.g + timefilter * (atms[ib].T.g - 2*atms[ic].T.g)
        atms[ic].q.g = atms[ic].q.g + timefilter * (atms[ib].q.g - 2*atms[ic].q.g)
        atms[ic].lnps.g = atms[ic].lnps.g + timefilter * (atms[ib].lnps.g - 2*atms[ic].lnps.g)

    #t1 = timer(); print '1st step of time filter: ',t1-t0; t0=t1
    
    # horizontal diffusion of VOR, DIV, T, q
    if flaghd:   atms[ib].fhoridiff(coord, dtime)

    #t1 = timer(); print 'horizontal diffusion: ',t1-t0; t0=t1

    # initial divergence damping
    runsec = runtime.days*86400+runtime.seconds
    if runsec < tdamp*86400.0:
        damprate = max(1-runsec/(tdamp*86400.0), 0)
        atms[ib].div.s /= 1+2.0*damprate

    #t1 = timer(); print 'inital DIV damping: ',t1-t0; t0=t1

    tmpvar = Field3d(nlayer=4*nlev+1)
    
    tmpvar.s[0:nlev]        = atms[ib].vor.s 
    tmpvar.s[nlev:nlev*2]   = atms[ib].div.s 
    tmpvar.s[nlev*2:nlev*3] = atms[ib].T.s 
    tmpvar.s[nlev*3:nlev*4] = atms[ib].q.s 
    tmpvar.s[nlev*4]        = atms[ib].lnps.s 
    
    fstog_3d(tmpvar.s[0:4*nlev+1], tmpvar.g[0:4*nlev+1], coord, nlayer=4*nlev+1)

    atms[ib].vor.g  = tmpvar.g[0:nlev] + 0
    atms[ib].div.g  = tmpvar.g[nlev:nlev*2] + 0
    atms[ib].T.g    = tmpvar.g[nlev*2:nlev*3] + 0 
    atms[ib].q.g    = tmpvar.g[nlev*3:nlev*4] + 0 
    atms[ib].lnps.g = tmpvar.g[nlev*4] + 0 


    #t1 = timer(); print 'VOR, DIV, T, q, lnps: ss -> gg : ',t1-t0; t0=t1

    # grid point values
    atms[ib].fgridvars(coord)
    
    #t1 = timer(); print 'diagnostic grid points variables',t1-t0; t0=t1
    
    # temperature correction
    if flaghd:   atms[ib].ftempcorr(coord, dtime)

    #t1 = timer(); print 'temperature correction : ',t1-t0; t0=t1

    # 2nd step of time filter: smoothing 
    if timefilter > 0:
        atms[ic].vor.g += timefilter * atms[ib].vor.g
        atms[ic].div.g += timefilter * atms[ib].div.g
        atms[ic].T.g += timefilter * atms[ib].T.g
        atms[ic].q.g += timefilter * atms[ib].q.g
        atms[ic].lnps.g += timefilter * atms[ib].lnps.g

    #t1 = timer(); print '2nd step of time filter : ',t1-t0; t0=t1
  
    # time step moves forward
    flags[0] = 1-flags[0] 
    flags[1] = 1-flags[1]

    return
#=============== End of semi-implicit time integration function ===============


def fprintvar(var):
    print var.name, np.min(var.g), np.max(var.g)

