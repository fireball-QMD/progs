from  pyfb.geometry.step import step
from  pyfb.geometry.atom import atom

class dinamic:
  read_charges=False

  def __init__(self,read_charges=False):
    self.step=[]
    self.read_charges=read_charges
    #lee las cargas de cada atomo despues de las posiciones:
    # x y z Qtot qs qp qd .....
    self.out=[]

  def print_total_steps(self):
    print(len(self.step))
  
  def append(self,step):
    self.step.append(step)

  def print(self):
    for i in self.step:
      i.print()

  def print_charges(self):
    for i in self.step:
      i.print_charges()

  def print_2line(self):
    for i in self.step:
      print(i.line2)

  def print_name(self):
    for i in self.step:
      print(i.name)

  def get(self,col):
    count=0.0
    aux=[]
    for i in range(len(col)):
       a=0.0
       aux.append(a)
    for i in self.step:
      salida=[]
      i.print()
      print(i.atom[col[0][1]-1].r[0])
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

  def print_out(self):
    print(self.out)
#    for i in range(len(self.out)):
#      print(self.out[i])

  def load_xyz(self,archivo,name=""):
    natoms = 0
    nstep = 0
    text=open(archivo).readlines()
    nmaxlines=len(text)
    i=0
    while i < nmaxlines:
      line=text[i].split()
      bas=step()
      if name != "":
        bas.name=name
      else:
        bas.name="step = "+str(nstep)
      if i == 0 :
        natoms = int(line[0])
        i=i+1
      bas.line2=text[i]
      for j in range(natoms):
        i=i+1
        line=text[i].split()
        a=line[0]
        ra=[]
        ra.append(float(line[1]))
        ra.append(float(line[2]))
        ra.append(float(line[3]))
        if self.read_charges:
          bas.atom[mod-1].Q=line[4]
          self[istep].atom[mod-1].setR(r)
          self[istep].atom[mod-1].setZ(line[0])
        bas.append(atom(a,ra))
      i=i+1 #natom
      i=i+1 #line2
      self.append(bas)
    
