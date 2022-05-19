#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])
print(sys.path)
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 

#Load Fdata
fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/relax/Fdata_HC_minimal"
fb.f2py_initbasics(fdatalocation)


#Load positions
din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
n_atomos=din.step[0].getNatoms()
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_natoms(n_atomos) 
fb.f2py_nucz(Zin) 
fb.f2py_ratom(pos) 

#Load options
fb.set_icluster(1)
fb.set_iquot(3)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)
fb.set_verbosity(1)


#run fireball 
fb.f2py_init()
fb.f2py_run()
fb.f2py_deallocate_all()
#fb.f2py_deallocate_charges()
#run fireball 
#fb.f2py_init()
#fb.f2py_run()

#Load new positions
din2=dinamic()
din2.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
n_atomos=din2.step[0].getNatoms()
pos=din2.step[0].getnumpy_pos()
Zin=np.array(din2.step[0].getZarray())
fb.f2py_natoms(n_atomos)
fb.f2py_nucz(Zin)
fb.f2py_ratom(pos)

#run fireball 
fb.f2py_init()
fb.f2py_run()


