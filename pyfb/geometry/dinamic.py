from  pyfb.geometry.structure import structure

class dinamic:

  def __init__(self):
    self.dinamic=[]

  def append(self,structure):
    self.dinamic.append(structure)

  def print(self):
    for i in self.dinamic:
      i.print()
