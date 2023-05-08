#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
import tarfile
from os.path import exists
sys.path.append(os.environ["FIREBALLHOME"])
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 

#Load Fdata
fdatalocation="/home/dani/Fdata_HCNOS"
fb.f2py_initbasics(fdatalocation)

#Load bas
din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
#n_atomos=din.step[0].getNatoms()
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
#fb.f2py_natoms(n_atomos) 
#fb.f2py_nucz(Zin) 
fb.f2py_getbas(Zin,pos) 

#Load options
#PARA que funcione poner en READFILES/readparam.f90
#idipole = 1
#icluster = 1
fb.set_icluster(1)
fb.set_iqout(7)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)
fb.set_idipole(1)
fb.set_iks(1)
fb.set_imcweda(0)
fb.set_idogs(0)
fb.set_verbosity(10)
fb.set_iwrtcharges(1)
fb.set_iwrtdipole(0)


#Run fireball 
fb.f2py_init()
ETOT=fb.f2py_getenergy()
din.step[0].line2=('ETOT = {0:12.6f}  eV '.format(ETOT))
fb.f2py_run()
fb.f2py_print_pcharges()
fb.f2py_deallocate_all()

#Load new positions
if exists("CHARGES"):
  remove("CHARGES")
din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/f2py/CH4.bas")
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
fb.f2py_getbas(Zin,pos)

#Run fireball 
fb.f2py_init()
ETOT=fb.f2py_getenergy()
din.step[0].line2=('ETOT = {0:12.6f}  eV '.format(ETOT))
fb.f2py_run()
fb.f2py_print_pcharges()
fb.f2py_deallocate_all()

