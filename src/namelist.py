from datetime import datetime, timedelta

# namelist

casename = 'fp_T30L18_spec_p_dp_initsh_rf10_nc5_nohd'
casedir  = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/'+casename+'/'
#casedir  = '/T3/yhy/ModelingGroup/XYZ_GCM/cases/'+casename+'/'
trun = ['T', 30]     # Truncation method and number: 'T' -- triangular
nlev = 18           # model levels; available nlev: 18, 26, 8

#=============================================================================
#   Model initialization
#=============================================================================
initialization = 100   # 1 steady-state; 2 barolinic wave; 3 Rossby-Hauritz wave
                       # 0 constant temperature; 4 mountain-induced Rossby wave
                       # 10 US standard atmospheric profile
                       # 100 start from netcdf file initfile
                       # 200 start from ERA Interim 6hourly grib files
Teinit = 1  # 0  constant temperature: 288K
            # 1  US standard atmospheric profile            
            # 3  ERA climatology
            # 4  the same as initial condition
inittime   = datetime(1, 1, 1)
initfile   = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/fp_T30L18_o_sh/fp_T30L18_o_sh_00010213_00-00-00.nc' 
#initfile   = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/fp_T30L18_o_dp/fp_T30L18_o_dp_00010213_00-00-00.nc' 
#initfile   = '/disk2/yhy/ModelingGroup/XYZ_GCM/cases/fp_T30L18_s_sh/fp_T30L18_s_sh_00011018_00-00-00.nc' 

#=============================================================================
#   Physics and forcing setting
#=============================================================================
artforcing = 5  # 0 no forcing; 1 Held-Suarez forcing;
                # 2 elevated heating; 3 bottom heating
                # 5 wave-CISK like (i.e. PBL convergence induced heating)
                # 6 vertical cosine heating with deep/shallow mode
                # 7 vertical and x: cosine; y: Gaussian
                # 200 using Q1 from ERA interim 6hourly data

heatingmode = 'deep'

heatingperiod = 10000*365 * 86400.0  # oscillation period (s)

Rayleigh_friction = True
taufriction = 1.0/(86400 * 10)  # relaxation time scale (s^-1)
pbltop = 700.0   # top of PBL (hPa)

Newtonian_cooling = True
taucooling = 1.0/(86400 * 5)  # relaxation time scale (s^-1)

precfrac = 1.0  # fraction of precipitation from heating

#=============================================================================
#   Time integration information
#=============================================================================
tunit    =  3600
timeend  = -24 * 200     # Ending time: + time steps; - per tunit second

#=============================================================================
#   Output NetCDF files   
# Output frequency: 0  no ouput; 
#                   +  time steps; 
#                   -  per tunit second
#=============================================================================
foutfreq = -6
#foutvars = ['T','PS','U','V','OMEGA','Z3','PREC']
foutvars = ['VOR','DIV','T','PS','U','V','OMEGA','Z3','PHIS','PREC']

#=============================================================================
#   Reference atmosphere for semi-implicity time integration
#=============================================================================
Tref  = [300.0 for i in range(nlev)]       # reference temperature (K)
psref = 1.0e5                              # reference surface pressure (Pa)

#=============================================================================
#   Smoothing, damping, diffusion
#=============================================================================
tdamp      = 0 # initial divergence damping days, recommanded value: 2 (days)
flaghd     = False # appling horizontal diffusion
timefilter = 0.06 # recommanded value: 0.06

# The following horizontal diffusion parameters are from CAM3
# kmnhd4: implicit del**2 form above level kmnhd4; 
#         implicit del**4 form at level kmnhd4 and below
# kmxhdc: courant number based truncation at level kmxhdc and above
# kmxhd2: increased del**2 coefficient at level kmxhd2 and above 
#         i.e. increased by 2**(kmxhd2-k+1)
kmnhd4 = 4      # CAM3, L26: 4
kmxhdc = 5      # CAM3, L26: 5
kmxhd2 = 2      # CAM3, L26: 2






#===============================================================================
#                                  constants
#===============================================================================
# Parameters
# ----------------------------------------------------------------------------
rad  = 6371229.0             # Earth's radius    (m)
omg  = 7.29212e-5            # Angular speed of Earth's rotation (1/s)
grav = 9.80616               # Gravity 


# These following parameters are from CESM codes
SHR_CONST_BOLTZ   = 1.38065e-23 # Boltzmann's constant ~ J/K/molecule
SHR_CONST_AVOGAD  = 6.02214e26  # Avogadro's number ~ molecules/kmole
SHR_CONST_MWDAIR  = 28.966      # molecular weight dry air ~ kg/kmole
SHR_CONST_MWWV    = 18.016      # molecular weight water vapor ~ kg/kmole
SHR_CONST_RGAS    = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ  # Universal gas constant ~ J/K/kmole
Rd   = SHR_CONST_RGAS/SHR_CONST_MWDAIR      # Dry air gas constant     ~ J/K/kg
Rv   = SHR_CONST_RGAS/SHR_CONST_MWWV        # Water vapor gas constant ~ J/K/kg
epsilo = SHR_CONST_MWWV/SHR_CONST_MWDAIR    # ratio of h2o to dry air molecular weights

Cpd  = 1004.64            # dry air heat capacity at const p (J/kg/K)
Cpv  = 1.810e3            # water vapor heat capacity at const p   (J/kg/K)
Cl   = 4.188e3            # liquid water heat capacity  (J/kg/K)
Lref = 2.501e6            # vaporization heat at Tref (J/kg)
P0   = 1.0e5



#   Size of grid, timestep, diffusion coefficient K4, K2
# ----------------------------------------------------------------------------
# Recommanded values, 09/26/2016
ntruns = [    21,     30,     42,     63,     85,    106,    170,    213,    340,    680,   1360]
nlons  = [    64,     96,    128,    192,    256,    320,    512,    640,   1024,   2048,   4096]
nlats  = [    32,     48,     64,     96,    128,    160,    256,    320,    512,   1024,   2048]
tsteps = [  3600,   2400,   1800,   1200,    900,    600,    450,    300,    200,    120,     60]
K4s    = [2.0e16, 1.5e16, 1.0e16, 2.5e15, 1.0e15, 5.0e14, 1.5e14, 7.0e13, 1.5e13, 1.0e12, 1.0e11]
K2 = 2.5e5  


#   Size of grid, timestep, diffusion coefficient K4, K2
# ----------------------------------------------------------------------------
M = trun[1]
N = trun[1]

itrun  = ntruns.index(trun[1])
nlon   = nlons[itrun]
nlat   = nlats[itrun]

timestep = tsteps[itrun]
K4 = K4s[itrun]


