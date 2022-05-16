#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])
print(sys.path)
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 


din=dinamic()
#din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/speedtest/serie/A/ini.bas")
din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
n_atomos=din.step[0].getNatoms()
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
#din.print()

fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/relax/Fdata_HC_minimal"
#fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/speedtest/Fdata"
fb.f2py_initbasics(fdatalocation)

fb.set_icluster(1)
fb.set_iquot(4)
fb.set_iquench(0)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)

fb.f2py_natoms(n_atomos) 
fb.f2py_nucz(Zin) 
fb.f2py_ratom(pos) 
#np.asfortranarray(pos)) #pos,Zin)
fb.f2py_init()
fb.f2py_run()
