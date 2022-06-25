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
fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/relax/Fdata_HC_minimal"
if not exists(fdatalocation):
  file = tarfile.open(os.environ["FIREBALLHOME"]+"/TESTS/Fdata.tar.gz")
  file.extractall(os.environ["FIREBALLHOME"]+"/TESTS/relax/")
  file.close()
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
fb.set_icluster(1)
fb.set_iquot(4)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)

#Run fireball 
fb.f2py_init()
ETOT=fb.f2py_getenergy()
din.step[0].line2=('ETOT = {0:12.6f}  eV '.format(ETOT))
fb.f2py_print_charges()
#print("CHARGES")
for i in range(1,din.step[0].getNatoms()+1):
  line=fb.f2py_charge(i).split()
  charge=0
  for j in range(len(line)):
    din.step[0].atom[i-1].q.append(line[j])
    charge=charge+float(line[j])
  din.step[0].atom[i-1].Q=charge
din.step[0].print_charges()

#fb.f2py_run()
#din=dinamic()
#din.laststep("restart.xyz",True)
#din.print_charges()


