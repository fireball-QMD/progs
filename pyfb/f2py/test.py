#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])
print(sys.path)
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np 

fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/speedtest/Fdata"
#fdatalocation="/home/dani/Fdata_HCNOS/"

print(fb.f2py_initbasics)
fb.f2py_initbasics(fdatalocation)
fb.f2py_options(1)

din=dinamic()
din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/speedtest/serie/A/ini.bas")
din.print()


n_atomos=din.step[0].getNatoms()
#fb.f2py_setnatoms(n_atomos)
#print(natoms)

pos=din.step[0].getnumpy_pos()
#print(pos)

Zin=np.array(din.step[0].getZarray())
#print(Zin)


fb.f2py_natoms(n_atomos) 
fb.f2py_nucz(Zin) 
print(pos)
fb.f2py_ratom(pos) 
#fb.f2py_loadbas(n_atomos,Zin,pos) 
#np.asfortranarray(pos)) #pos,Zin)
fb.f2py_init()
