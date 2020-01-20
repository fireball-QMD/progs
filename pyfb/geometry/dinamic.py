from  pyfb.geometry.structure import structure

class dinamic:
  read_charges=False

  def __init__(self,read_charges=False):
    self.structure=[]
    self.read_charges=read_charges
    #lee las cargas de cada atomo despues de las posiciones:
    # x y z Qtot qs qp qd .....

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
    for i in self.structure:
      salida=""
      for c in col: 
       if c[0] == 'x':
         salida=salida+" "+i.atom[c[1]-1].r[0]
       if c[0] == 'y':
         salida=salida+" "+i.atom[c[1]-1].r[1]
       if c[0] == 'z':
         salida=salida+" "+i.atom[c[1]-1].r[2]
       if c[0] == 'd':
         salida=salida+" "+str(i.atom[c[1]-1].distancia(i.atom[c[2]-1]))
      print(salida)

