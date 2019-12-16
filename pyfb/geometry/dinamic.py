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
