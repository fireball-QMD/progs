from  pyfb.geometry.step import step
from  pyfb.geometry.atom import atom
from pyfb.geometry.tabla import tabla
import pandas as pd
import numpy as np
import random


tabla=tabla()

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

  def print_enlaces(self):
    for bas in self.step:
      bas.print_enlaces()

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


  def loadxyz_fromString(self,cadena):
    natoms = 0
    text = cadena.split("\n")
    print(text)
    nmaxlines = len(text)
    natoms = int(text[0])
    bas=step()
    for j in range(natoms):
      line=text[j+2].split()
      print(line)
      a=line[0]
      ra=[]
      ra.append(float(line[1]))
      ra.append(float(line[2]))
      ra.append(float(line[3]))
      bas.append(atom(a,ra))
    self.append(bas)
    

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

  def rotX(self, istep,ang):
    self.step[istep-1].rotX(ang)

  def rotY(self, istep,ang):
    self.step[istep-1].rotY(ang)

  def rotZ(self, istep,ang):
    self.step[istep-1].rotZ(ang)
    
  def center(self, istep):
    self.step[istep-1].center()

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
  
  def merge(self,i,j): 
    """
    merge 2 steps of dinamic, -laststep2xyz also merge 2 files.xyz
    """
    step1=self.step[i-1].atom
    step2=self.step[j-1].atom
    N=(len(step1))+(len(step2))
    print(N)
    print()
    for s1 in step1:
      s1.print() 
    for s2 in step2:
      s2.print() 
  
  def distance(self,i,j): 
    """
    minimal distance 2 steps of dinamic
    -laststep2xyz distance 2 files.xyz
    """
    step1=self.step[i-1].atom
    step2=self.step[j-1].atom
    dmin=step1[0].distancia(step2[0])
    ri=step1[0].getRadio()
    rj=step2[0].getRadio()
    radiomin=(float(ri+rj))/2/self.step[0].tol/100
    for s1 in step1:
      for s2 in step2:
        d=s1.distancia(s2)
        ri=s1.getRadio()
        rj=s2.getRadio()
        radio=(float(ri+rj))/2/self.step[0].tol/100
        if dmin > d:
          dmin=d
          radiomin=radio
    #      print(dmin,radiomin) 
    return dmin,radiomin

  def join(self,i,j):
    """
    join 2 steps of dinamic, -laststep2xyz also join 2 files.xyz
    take in acount that dmin > radiomin, move random the j step
    and do a random rotate 
    """
    d,radio=self.distance(i,j)
    while d < radio:
      A=0.2*np.array([random.uniform(-1, 1),random.uniform(-1, 1),random.uniform(-1, 1)])
      self.step[j-1].rotX(np.pi*random.uniform(-1, 1))
      self.step[j-1].rotY(np.pi*random.uniform(-1, 1))
      self.step[j-1].rotZ(np.pi*random.uniform(-1, 1))
      for s in range(0,len(self.step[j-1].atom)):
        for k in range(0,3):
          self.step[j-1].atom[s].r[k]+=A[k]
    
      d,radio=self.distance(i,j)
    self.merge(i,j)
    #d,radio=self.distance(i,j)
    #print(d,radio)

  def joinrnd(self,i,j):
    """
    Use de max dimensión of the first step to put random the 2 step
    """
    self.step[i-1].setdim()
    #cambiamos RC dentro de los limites del step i, vaya tela, 
    #es para que no se nos vaya, ya que no tenemos encuenta lvs
    Desplaz=np.array([0.0,0.0,0.0])
    Desplaz=random.uniform(np.array(self.step[i-1].dim[0]),np.array(self.step[i-1].dim[1]))
    dentro=True
    self.step[j-1].center()
    for a in self.step[j-1].atom:
      a.r=a.r+Desplaz
    d,radio=self.distance(i,j)
    #en el caso de que ese muy cerca dejamos que se mueva pero si se sale del dim lo volvemos a meter
    while d < radio and not dentro:
      #rotamos de forma aleatoria
      self.step[j-1].rotX(np.pi*random.uniform(-1, 1))
      self.step[j-1].rotY(np.pi*random.uniform(-1, 1))
      self.step[j-1].rotZ(np.pi*random.uniform(-1, 1)) 
      #desplazamos 
      A=0.2*np.array([random.uniform(-1, 1),random.uniform(-1, 1),random.uniform(-1, 1)])
      for a in self.step[j-1].atom:
        a.r+=A
      #chequeamos si nos salimos
      dentro=True
      for k in range(0,3):
        if self.step[j-1].RC[k]+A[k] < self.step[i-1].dim[0][k]:
          dentro=False
        if self.step[j-1].RC[k]+A[k] > self.step[i-1].dim[1][k]:
          dentro=False
 
      if dentro: #aceptamos y actulizamos el centro
        self.step[j-1].center() 
      else: #media vuelta
        for a in self.step[j-1].atom:
          a.r-=A
      d,radio=self.distance(i,j) #vemos como estamos de cerca de algo

    self.merge(i,j)
    
 
  def loadbas(self,archivo,name=""):
    bas=step()
    text=open(archivo).readlines()
    nmaxlines=len(text)
    natoms = int(text[0].split()[0])
    for i in range(1,natoms+1):
      line=text[i].split()
      a=tabla.getZ(int(line[0]))
      ra=[]
      ra.append(float(line[1]))
      ra.append(float(line[2]))
      ra.append(float(line[3]))
      bas.append(atom(a,ra))
    self.append(bas)


 
