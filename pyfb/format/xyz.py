from  pyfb.geometry.dinamic import dinamic
import pyfb

dinamic=dinamic()

def load_xyz(archivo,name=""):
  natoms = 0
  nstep = 0
  nlinea = 0
  bas=pyfb.geometry.structure.structure()
  if name != "":
    bas.name=name
  for linea in open(archivo):
    line = linea.split()
    nlinea+=1
    if nlinea == 1 :
      natoms = int(line[0])
      for i in range(0,natoms):
        bas.append(pyfb.geometry.atom.atom())
#    print(nlinea, line)
    mod=(nlinea-2) % (natoms + 2)
    if mod == natoms:
      nstep+=1
#      print ("nstep++")
    if mod == 0 :
      bas.line2=linea 
    if mod > 0 and mod < (natoms+1):
      r=[]
      r.append(line[1])
      r.append(line[2])
      r.append(line[3])
      bas.atom[mod-1].setR(r)
      bas.atom[mod-1].setZ(line[0])
      if mod == natoms:
        dinamic.append(bas)

def print_xyz():
  dinamic.print()

def print_2line():
  dinamic.print_2line()
