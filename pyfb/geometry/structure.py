from  pyfb.geometry.atom import atom

class structure:

  def __init__(self):
    self.structure=[]
    self.line2=""
    self.name=""

  def append(self,atomo):
    self.structure.append(atomo)

  def print(self):
    linea=str(len(self.structure))+"\n"
    if self.line2 != "":
      linea=linea+self.line2
      linea = linea[:-1]
    print(linea)
    for i in self.structure:
      i.print()  

