from  pyfb.geometry.structure import structure
from  pyfb.geometry.atom import atom

class dinamic:
  read_charges=False

  def __init__(self,read_charges=False):
    self.structure=[]
    self.read_charges=read_charges
    #lee las cargas de cada atomo despues de las posiciones:
    # x y z Qtot qs qp qd .....
    self.out=[]

  def append(self,structure):
    self.structure.append(structure)

  def print(self):
    for i in self.structure:
      i.print()

  def print_charges(self):
    for i in self.structure:
      i.print_charges()

  def print_2line(self):
    for i in self.structure:
      print(i.line2)

  def print_name(self):
    for i in self.structure:
      print(i.name)

  def get(self,col):
    count=0.0
    aux=[]
    for i in range(len(col)):
       a=0.0
       aux.append(a)
    for i in self.structure:
      salida=[]
      for c in range(len(col)): 
       if col[c][0] == 'x':
         salida.append(i.atom[col[c][1]-1].r[0])
       if col[c][0] == 'y':
         salida.append(i.atom[col[c][1]-1].r[1])
       if col[c][0] == 'z':
         salida.append(i.atom[col[c][1]-1].r[2])
       if col[c][0] == 'd':
         salida.append(str(i.atom[col[c][1]-1].distancia(i.atom[col[c][2]-1])))
       if col[c][0] == 'ang':
         salida.append(str(i.atom[col[c][1]-1].ang(i.atom[col[c][2]-1],i.atom[col[c][3]-1])))
       if col[c][0] == 'X':
         aux[c]=float(i.atom[col[c][1]-1].r[0])/(count+1)+float(aux[c])*count/(count+1)
         salida.append(str(aux[c]))
       if col[c][0] == 'Y':
         aux[c]=float(i.atom[col[c][1]-1].r[1])/(count+1)+float(aux[c])*count/(count+1)
         salida.append(str(aux[c]))
       if col[c][0] == 'Z':
         aux[c]=float(i.atom[col[c][1]-1].r[2])/(count+1)+float(aux[c])*count/(count+1)
         salida.append(str(aux[c]))
       if col[c][0] == 'D':
         aux[c]=float(i.atom[col[c][1]-1].distancia(i.atom[col[c][2]-1]))/(count+1)+float(aux[c])*count/(count+1)
         salida.append(str(aux[c]))
       if col[c][0] == 'ANG':
         aux[c]=float(i.atom[col[c][1]-1].ang(i.atom[col[c][2]-1],i.atom[col[c][3]-1]))/(count+1)+float(aux[c])*count/(count+1)
         salida.append(str(aux[c]))
      count=count+1
#      print(salida)
      self.out.append(salida)

  def load_xyz(self,archivo,name=""):
    natoms = 0
    nstep = 0
    nlinea = 0
    bas=structure()
    if name != "":
      bas.name=name
    for linea in open(archivo):
      line = linea.split()
      nlinea+=1
      if nlinea == 1 :
        natoms = int(line[0])
        for i in range(0,natoms):
          bas.append(atom())
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
        if self.read_charges:
          bas.atom[mod-1].Q=line[4]
        bas.atom[mod-1].setR(r)
        bas.atom[mod-1].setZ(line[0])
        if mod == natoms:
          self.append(bas)
    

