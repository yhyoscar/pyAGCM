
from time import time as timer
import os

from namelist import *
from DynCore  import *

#=============================================================================
os.system('mkdir '+casedir);    os.system('cp ./namelist.py '+casedir)

coord  = Coordinate()
refatm = ReferenceAtmosphere(coord, dtime=timestep)
atms   = [Atmosphere(), Atmosphere()]

modeltime = inittime
atms[0].finitial(coord, modeltime, iflag=initialization)
atms[1].finitial(coord, modeltime, iflag=initialization)
flags  = [0,1]

if timeend>0: endtime = int(timeend*timestep)
if timeend<0: endtime = int(-timeend*tunit)

if abs(foutfreq)>0:
    if foutfreq>0: foutdtime = int(foutfreq*timestep)
    if foutfreq<0: foutdtime = int(-foutfreq*tunit)


# start time integration
t0 = timer()
ftimeinte(atms, flags, coord, refatm, dtime=timestep/2, output=True, runtime=modeltime-inittime, modeltime=modeltime)
print 'step = 0, run time = ',0,', model time = ',modeltime,', cpu time = ',timer()-t0
print '================================================================='

for istep in range(1, endtime/timestep+1):
    runtime   = timedelta(seconds=istep*timestep)
    modeltime = modeltime + timedelta(seconds=timestep)
    flagoutput = False
    if np.mod(istep*timestep, foutdtime)==0: flagoutput = True
    ftimeinte(atms, flags, coord, refatm, dtime=timestep, output=flagoutput, runtime=modeltime-inittime, modeltime=modeltime)
    print 'step = ',istep,', run time = ',modeltime-inittime,', model time = ',modeltime,', cpu time = ',timer()-t0
    #print '================================================================='




