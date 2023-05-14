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
#fdatalocation="/home/dani/Fdata_HC/"
if not exists(fdatalocation):
  file = tarfile.open(os.environ["FIREBALLHOME"]+"/TESTS/Fdata.tar.gz")
  file.extractall(os.environ["FIREBALLHOME"]+"/TESTS/relax/")
  file.close()
idipole=0
fb.f2py_initbasics_opt(fdatalocation,idipole)

#Load options
fb.set_icluster(1)
fb.set_iqout(3)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)
fb.set_idipole(idipole)

def delauxfiles():
  for i in ['CHARGES','ac.dat','answer.bas','Charges_and_Dipoles','dipole_Tot','dipole_Tot_proy','restart.xyz','xv.dat','dipole_Qout','PCHARGES']:
    if exists(i):
      os.remove(i)


def runFB(pos):
  delauxfiles()
  #Load new positions
  din=dinamic()
  din.loadbas(pos)
  pos=din.step[0].getnumpy_pos()
  Zin=np.array(din.step[0].getZarray())
  fb.f2py_getbas(Zin,pos)
  fb.f2py_init()
  fb.f2py_run()
  fb.f2py_print_charges()
  fb.f2py_deallocate_all()

runFB(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
runFB(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")

