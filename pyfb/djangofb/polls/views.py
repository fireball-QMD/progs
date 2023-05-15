#from django.shortcuts import render
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, redirect
import tarfile
import sys
import os
from os.path import exists
sys.path.append(os.environ["FIREBALLHOME"])
import libpyfb as fb
from pyfb.geometry.dinamic import *
import numpy as np
from openbabel import pybel

mol = pybel.readstring("smi", "C")
mol.make3D()
peticion=mol.write("xyz")
calculando=False

#Para probar con una base mas pequeña
BASE_minima=False
BASE_minima=True

fdatalocation=""
#Load Fdata
if BASE_minima:
  fdatalocation=os.environ["FIREBALLHOME"]+"/TESTS/relax/Fdata_HC_minimal"
else:
  fdatalocation="/home/dani/Fdata_HCNOS"
#  fdatalocation="/home/dani/Fdata_HC"

if not exists(fdatalocation):
  file = tarfile.open(os.environ["FIREBALLHOME"]+"/TESTS/Fdata.tar.gz")
  file.extractall(os.environ["FIREBALLHOME"]+"/TESTS/relax/")
  file.close()


idipole=0
fb.f2py_initbasics_opt(fdatalocation,idipole)
#fb.f2py_initbasics(fdatalocation)

#Load bas
din=dinamic()
if BASE_minima:
  din.loadbas(os.environ["FIREBALLHOME"]+"/TESTS/relax/input.bas")
else:
  din.loadbas(os.environ["FIREBALLHOME"]+"/pyfb/djangofb/polls/input.bas")
n_atomos=din.step[0].getNatoms()
pos=din.step[0].getnumpy_pos()
Zin=np.array(din.step[0].getZarray())
#fb.f2py_natoms(n_atomos) 
#fb.f2py_nucz(Zin) 
fb.f2py_getbas(Zin,pos) 


#Load options
fb.set_icluster(1)
fb.set_iqout(7)
fb.set_iquench(-1)
fb.set_dt(0.5)
fb.set_nstepf(1)
fb.set_iwrtxyz(1)
fb.set_idipole(idipole)
fb.set_iks(1)
fb.set_imcweda(0)
fb.set_idogs(0)
fb.set_verbosity(10)
fb.set_iwrtcharges(1)
fb.set_iwrtdipole(0)

#Run fireball 
fb.f2py_init()
fb.f2py_run()
#ETOT=fb.f2py_getenergy()
#print(ETOT)
#fb.f2py_print_pcharges()
#fb.f2py_print_charges()

atomo_infodat=[]
carga_infodat=[]
archivo=fdatalocation+"/info.dat"
text=open(archivo).readlines()
atomos_cargados=""
linea=0
for line in text:
  if(len(line.split())>2): 
    if line.split()[2] == 'Element':
      atomos_cargados=atomos_cargados+line.split()[0]+","
      atomo_infodat.append(int(text[linea+1].split()[0]))
      q=0.0
      for iq in text[linea+8].split():
        q=q+float(iq)
      carga_infodat.append(q)
  linea=linea+1

def index(request):
  #print('****************************')
  #print(request.POST) 
  #print('****************************')
  if request.POST.get('peticion_pybel') != "Load xyz-format to get charges":
    return render(request, 'polls/index.html',{ 'peticion': peticion, 'atomos_cargados': atomos_cargados })
  else:
    error=False
    aux=""
    
    try:
      aux=pybel.readstring(request.POST['format'],request.POST['input'])
    except (IOError):
      error=True

    if error:
      peticion_aux="Open Babel Error  in ReadMolecule"
    else:
      peticion_aux=aux.write("xyz")

    return render(request, 'polls/index.html',{ 'peticion': peticion_aux, 'atomos_cargados': atomos_cargados })

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def ERROR(peticion):
  nline=0
  error=False
  N=int(peticion.splitlines()[0].split()[0])
  for i in peticion.splitlines():
    nline+=1
    if nline > 2 and nline < (N+3):
      if (len(i.split()) > 3):
        if (((i.split()[0]) in atomos_cargados) == False) or (isfloat(i.split()[1]) == False) or (isfloat(i.split()[2]) == False) or (isfloat(i.split()[3]) == False) :
          error=True
      else:
          error=True
  return error

def delauxfiles():
  for i in ['CHARGES','ac.dat','answer.bas','Charges_and_Dipoles','dipole_Tot','dipole_Tot_proy','restart.xyz','xv.dat','dipole_Qout','PCHARGES']:
    if exists(i):
      os.remove(i)
   
def runFB(peticion):
  fb.f2py_deallocate_all()
  delauxfiles()
  din=dinamic()
  din.loadxyz_fromString(peticion)
  pos=din.step[0].getnumpy_pos()
  Zin=np.array(din.step[0].getZarray())
  fb.f2py_getbas(Zin,pos)
  fb.f2py_init()
  fb.f2py_run()

  aux=pybel.readstring("xyz",peticion)
  mol2=aux.write("mol2")
  for i in range(1,len(aux.atoms)+1):
    aux.atoms[i-1].OBAtom.SetPartialCharge(float(fb.f2py_pcharge(i)))
    mol2=aux.write("mol2")
  peticion=mol2.replace("GASTEIGER","Mulliken-dipole")
  pdb=aux.write("pdb")
  return pdb,peticion

def getPositions(request):
    #print('****************************************')
    #print(request.POST)
    #print(request.POST.get('fileSubmit'))
    #print('****************************************')
    peticion=request.POST['input']
    if request.POST.get('pybel') != "Use pybel to get xyz-format":
      if ERROR(peticion) == False:
        pdb,peticion = runFB(peticion)
        return render(request, 'polls/getpositions.html',{ 'peticion': peticion , 'pdb' : pdb })
      else:
        peticion=peticion+"\r There are some error in the imput"
        return render(request, 'polls/index.html',{ 'peticion': peticion, 'atomos_cargados': atomos_cargados  })
    else:
      mol = pybel.readstring("smi", "C=C")
      mol.make3D()
      peticion=mol.write("pdb")
      formatos=pybel.informats.keys() 
      return render(request, 'polls/pybel.html',{ 'formatos': formatos, 'peticion': peticion })
