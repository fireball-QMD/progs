from  pyfb.geometry.step import step
from  pyfb.geometry.atom import atom
import pandas as pd
import numpy as np

class dinamic:
  def __init__(self):
    self.step=[]
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

  def print_bas_format(self):
    for i in self.step:
      i.print_bas_format()

  def print_charges(self):
    for i in self.step:
      i.print_charges()

  def print_2line(self):
    for i in self.step:
      print(i.line2)

  def print_name(self):
    for i in self.step:
      print(i.name)

  def print_out(self):
    #print(self.out)
    for i in range(len(self.out[0])):
      a=""
      for j in range(len(self.out)):
        a=a+'{0:12.6f}  '.format(self.out[j][i])
      print(a)

  def dintra_matrix(self):
    for bas in self.step:
      out=[]
      col=[]
      for iatom in bas.atom:
        col.append(iatom.Z)
        aux=[]
        for jatom in bas.atom:
          aux.append(iatom.distancia(jatom))
        out.append(aux)
      c=np.array(out)
      df2 = pd.DataFrame(c,index=col,columns=col)
      print(df2.to_string())

  def load_enlaces(self):
      for bas in self.step: #solo debería haber un step cargado
          bas.load_enlaces()
          #bas.print_enlaces()

  def get(self,info,col):
    count=0.0
    aux=[]
    salida=[]
    for i in self.step:
      if info == '-rescal':
        for j in i.atom:
          j.r[0]=j.r[0]*float(col[0])
          j.r[1]=j.r[1]*float(col[0])
          j.r[2]=j.r[2]*float(col[0])

      if info == '-x':
        salida.append(i.atom[col[0]-1].r[0])
      if info == '-y':
        salida.append(i.atom[col[0]-1].r[1])
      if info == '-z':
        salida.append(i.atom[col[0]-1].r[2])

      if info == '-dx':
        i.atom[col[0]-1].r[0]=i.atom[col[0]-1].r[0]+float(col[1])
      if info == '-dy':
        i.atom[col[0]-1].r[1]=i.atom[col[0]-1].r[1]+float(col[1])
      if info == '-dz':
        i.atom[col[0]-1].r[2]=i.atom[col[0]-1].r[2]+float(col[1])


      if info == '-X':
        if int(count)==0:
          salida.append(float(i.atom[col[0]-1].r[0]))
        else:
          salida.append(float(i.atom[col[0]-1].r[0])/(count+1)+float(salida[int(count-1)])*count/(count+1))
      if info == '-Y':
        if int(count)==0:
          salida.append(float(i.atom[col[0]-1].r[1]))
        else:
          salida.append(float(i.atom[col[0]-1].r[1])/(count+1)+float(salida[int(count-1)])*count/(count+1))
      if info == '-Z':
        if int(count)==0:
          salida.append(float(i.atom[col[0]-1].r[2]))
        else:
          salida.append(float(i.atom[col[0]-1].r[2])/(count+1)+float(salida[int(count-1)])*count/(count+1))

      if info == '-d':
        salida.append(float(i.atom[col[0]-1].distancia(i.atom[col[1]-1])))
      if info == '-D':
        if int(count)==0:
          salida.append(float(i.atom[col[0]-1].distancia(i.atom[col[1]-1])))
        else:
          salida.append(float(i.atom[col[0]-1].distancia(i.atom[col[1]-1]))/(count+1)+float(salida[int(count-1)])*count/(count+1))

      if info == '-ang':
        salida.append(float(i.atom[col[0]-1].ang(i.atom[col[1]-1],i.atom[col[2]-1])))
      if info == '-ANG':
        if int(count)==0:
          salida.append(float(i.atom[col[0]-1].ang(i.atom[col[1]-1],i.atom[col[2]-1])))
        else:
          salida.append(float(i.atom[col[0]-1].ang(i.atom[col[1]-1],i.atom[col[2]-1]))/(count+1)+float(salida[int(count-1)])*count/(count+1))
      count=count+1
#     print(salida)
    self.out.append(salida)


  def loadxyz(self,archivo,name=""):
    natoms = 0
    istep = 0
    text=open(archivo).readlines()
    nmaxlines=len(text)
    i=0
    while i < nmaxlines:
      line=text[i].split()
      bas=step()
      istep=istep+1
      if name != "":
        bas.name=name
      else:
        bas.name="step = "+str(istep)
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
        bas.append(atom(a,ra))
      i=i+1 #natom
      i=i+1 #line2
      self.append(bas)
    

  def loadstep(self,archivo,istep):
    natoms = 0
    text=open(archivo).readlines()
    nmaxlines=len(text)
    i=0
    while i < nmaxlines:
      line=text[i].split()
      bas=step()
      bas.name="step = "+str(istep)
      if i == 0 :
        natoms = int(line[0])
      i=(istep-1)*(natoms+2)+1
      bas.line2=text[i]
      for j in range(natoms):
        i=i+1
        line=text[i].split()
        a=line[0]
        ra=[]
        ra.append(float(line[1]))
        ra.append(float(line[2]))
        ra.append(float(line[3]))
        bas.append(atom(a,ra))
      i=nmaxlines
      self.append(bas)
 
  def laststep(self,archivo,read_charges):
    natoms = 0
    text = open(archivo).readlines()
    nmaxlines = len(text)
    natoms = int(text[0])
    bas=step()
    i=nmaxlines-(natoms+1)
    bas.line2=text[i]
    for j in range(natoms):
      i=i+1
      line=text[i].split()
      a=line[0]
      ra=[]
      ra.append(float(line[1]))
      ra.append(float(line[2]))
      ra.append(float(line[3]))
      bas.append(atom(a,ra))
    if read_charges: #load after read de step
      bas.loadcharges() 
    self.append(bas) 
