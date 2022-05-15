from pyfb.geometry.atom import atom
from pyfb.geometry.tabla import tabla
import numpy as np 

tabla=tabla()

class step:

  def __init__(self):
    self.atom=[]
    self.charges=[]
    self.line2=""
    self.name=""
    self.tol=4

  def append(self,atomo):
    self.atom.append(atomo)

  def load_enlaces(self):
    for i in range(len(self.atom)-1):
        for j in range(i+1,len(self.atom)):
            ri=self.atom[i].getRadio()
            rj=self.atom[j].getRadio()
            radio=(float(ri+rj))/2/self.tol
            if((100*self.atom[i].distancia(self.atom[j])) < radio ):
                 self.atom[i].neighbor.append(self.atom[j])
                 self.atom[j].neighbor.append(self.atom[i])

  def print_enlaces(self):
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

  def print_charges(self):
    linea=str(len(self.atom))+"\n"
    if self.line2 != "":
      linea=linea+self.line2
      linea = linea[:-1]
    print(linea)
    for i in self.atom:
      i.print_charges()  

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
