from pyfb.geometry.atom import atom
from pyfb.geometry.tabla import tabla
import numpy as np 
import random

tabla=tabla()


class step:

  def __init__(self):
    self.atom=[]
    self.charges=[]
    self.line2=""
    self.name=""
    self.tol=4
    self.enlaces=[] #conjunto enlaces primeros vecinos
    self.RC=np.array([0.0,0.0,0.0])
    self.dim=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]])

  def append(self,atomo):
    self.atom.append(atomo)

  def getEnergy(self):
    return self.line2.split()[2]

  def getNumberof(self,ele):
    aux=0
    for i in self.atom:
      if i.Z == ele:
        aux=aux+1
    return aux

  def load_enlaces(self):
    for i in range(len(self.atom)-1):
        for j in range(i+1,len(self.atom)):
            ri=self.atom[i].getRadio()
            rj=self.atom[j].getRadio()
            radio=(float(ri+rj))/2/self.tol
            if((100*self.atom[i].distancia(self.atom[j])) < radio ):
                 self.atom[i].neighbor.append(self.atom[j])
                 self.atom[j].neighbor.append(self.atom[i])
                 aux=[]
                 aux.append(self.atom[i].Z)
                 aux.append(i+1)
                 aux.append(self.atom[j].Z)
                 aux.append(j+1)
                 aux.append(self.atom[i].distancia(self.atom[j])) 
                 self.enlaces.append(aux)

  def print_enlaces(self):
    for i in self.enlaces:
      print(f'{i[0]} ({i[1]}) - {i[2]} ({i[3]}) {i[4]:.2f}')


  def print_enlaces_atomo(self):
    for i in self.atom:
      print(i.Z,self.atom.index(i)+1)
      for j in i.neighbor:
        radio=float(i.getRadio()+j.getRadio())/2/self.tol
        print("  - "+str(j.Z)+str(self.atom.index(j)+1)+" "+str(radio)+" "+str(100*i.distancia(j)))
      print("  ")
      print("  ")


  def print(self):
    linea=str(len(self.atom))+"\n"
    if self.line2 != "":
      linea=linea+self.line2
      linea = linea[:-1]
    print(linea)
    for i in self.atom:
      i.print()  


  def rotZ0(self,ang):
    for i in self.atom:
      aux=np.array([0.0,0.0,0.0])
      aux[0]= i.r[0]*np.cos(ang)+i.r[1]*np.sin(ang)
      aux[1]=-i.r[0]*np.sin(ang)+i.r[1]*np.cos(ang)
      aux[2]=i.r[2]
      i.r=aux

  def rotY0(self,ang):
    for i in self.atom:
      aux=np.array([0.0,0.0,0.0])
      aux[0]= i.r[0]*np.cos(ang)+i.r[2]*np.sin(ang)
      aux[1]=-i.r[0]*np.sin(ang)+i.r[2]*np.cos(ang)   
      aux[2]=i.r[1]      
      i.r=aux

  def rotX0(self,ang):
    """
    ang = radianes, rot from 0,0,0                      
    """
    for i in self.atom:
      aux=np.array([0.0,0.0,0.0])
      aux[0]= i.r[1]*np.cos(ang)+i.r[2]*np.sin(ang)
      aux[1]=-i.r[1]*np.sin(ang)+i.r[2]*np.cos(ang)   
      aux[2]=i.r[0]      
      i.r=aux
  
  def rotX(self,ang):
    self.setRC()
    for i in self.atom:
      i.r=i.r-self.RC
    self.rotX0(ang)
    for i in self.atom:
      i.r=i.r+self.RC    

  
  def rotY(self,ang):
    self.setRC()
    for i in self.atom:
      i.r=i.r-self.RC
    self.rotY0(ang)
    for i in self.atom:
      i.r=i.r+self.RC    

  
  def rotZ(self,ang):
    self.setRC()
    for i in self.atom:
      i.r=i.r-self.RC
    self.rotZ0(ang)
    for i in self.atom:
      i.r=i.r+self.RC    


  def setRC(self):
    for i in self.atom:
      self.RC=self.RC+i.r
    self.RC=self.RC/len(self.atom)

  def center(self):
    self.setRC()
    for i in self.atom:
      i.r=i.r-self.RC
 
  def setdim(self):
    s=0
    self.dim=([
      [self.atom[s].r[0],self.atom[s].r[1],self.atom[s].r[2]],
      [self.atom[s].r[0],self.atom[s].r[1],self.atom[s].r[2]]
      ])
    for s in self.atom:
      for k in range(0,3):
        if self.dim[0][k]>s.r[k]:
          self.dim[0][k]=s.r[k] 
        if self.dim[1][k]<s.r[k]:
          self.dim[1][k]=s.r[k]

  def print_charges(self):
    linea=str(len(self.atom))+"\n"
    if self.line2 != "":
      linea=linea+self.line2
      linea = linea[:-1]
    print(linea)
    for i in self.atom:
      i.print_charges()  

  def print_charges_to_string(self):
    out=""
    linea=str(len(self.atom))+"\n"
    if self.line2 != "":
      linea=linea+self.line2
      linea = linea[:-1]
    out=out+linea+"\n"
    for i in self.atom:
      out=out+i.print_charges_to_string()+"\n"
    return out


  def print_bas_format(self):
    print(len(self.atom))
    for i in self.atom:
      i.print_bas_format()

  def loadcharges(self):
    txt=open("CHARGES").readlines()
    for i in range(1,len(txt)):
      charge=0
      line=txt[i].split()
      for j in range(len(line)):
        self.atom[i-1].q.append(line[j])
        charge=charge+float(line[j])
      self.atom[i-1].Q=charge

  def getNatoms(self):
    print(len(self.atom),"N")
    return len(self.atom)

  def getnumpy_pos(self):
    aux=[]
    for i in self.atom:
      aux.append(i.r)
    return np.array(aux) 

  def getZarray(self):
     aux=[]
     for i in self.atom:
       aux.append(tabla.getN(i.Z))
     return aux
