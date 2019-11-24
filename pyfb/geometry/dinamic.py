from  pyfb.geometry.structure import structure

class dinamic:

  def __init__(self):
    self.structure=[]

  def append(self,structure):
    self.structure.append(structure)

  def print(self):
    for i in self.structure:
      i.print()

