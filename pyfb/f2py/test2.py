#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])
print(sys.path)
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 
import tarfile
from os.path import exists
from os import remove


#Load Fdata
fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/relax/Fdata_HC_minimal"
if not exists(fdatalocation):
  file = tarfile.open(os.environ["FIREBALLHOME"]+"/TESTS/Fdata.tar.gz")
  file.extractall(os.environ["FIREBALLHOME"]+"/TESTS/relax/")
  file.close()

fb.f2py_initbasics(fdatalocation)


#Load positions
if exists("CHARGES"):
  remove("CHARGES")
din=dinamic()
#din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
din.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_getbas(Zin,pos)

#Load options
fb.set_icluster(1)
fb.set_iqout(3)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)
fb.set_verbosity(1)


#run fireball 
fb.f2py_init()
fb.f2py_run()
fb.f2py_deallocate_all()

#Load new positions
if exists("CHARGES"):
  remove("CHARGES")
din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_getbas(Zin,pos)

#run fireball 
fb.f2py_init()
fb.f2py_run()
fb.f2py_deallocate_all()

#Load new positions
if exists("CHARGES"):
  remove("CHARGES")
din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_getbas(Zin,pos)

#run fireball 
fb.f2py_init()
fb.f2py_run()
fb.f2py_deallocate_all()

#Load new positions
if exists("CHARGES"):
  remove("CHARGES")
din=dinamic()
#din.loadbas("/home/dani/bases/H2O.bas")
din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_getbas(Zin,pos)

#run fireball 
fb.f2py_init()
fb.f2py_run()
fb.f2py_deallocate_all()
