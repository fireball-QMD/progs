import sys
import os
import tarfile
from os.path import exists
sys.path.append(os.environ["FIREBALLAPP"])
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 


def delauxfiles():
  for i in ['CHARGES','ac.dat','answer.bas','Charges_and_Dipoles','dipole_Tot','dipole_Tot_proy','restart.xyz','xv.dat','dipole_Qout','PCHARGES']:
    if exists(i):
      os.remove(i)

def runFB(pos):
  delauxfiles()
  fb.f2py_deallocate_all()
  #Load new positions
  din=dinamic()
  din.loadbas(pos)
  pos=din.step[0].getnumpy_pos()
  Zin=np.array(din.step[0].getZarray())
  fb.f2py_getbas(Zin,pos)
  fb.f2py_init()
  fb.f2py_run()
  fb.f2py_print_pcharges()
  fb.f2py_deallocate_all()
  delauxfiles()


#Load Fdata
fdatalocation=os.environ["FIREBALLAPP"]+"/../Fdata_HC"
fb.f2py_initbasics(fdatalocation)

#run Fireball
runFB(os.environ["FIREBALLAPP"]+"/TESTS/relax/input.bas")
runFB(os.environ["FIREBALLAPP"]+"/TESTS/relax/input.bas")
