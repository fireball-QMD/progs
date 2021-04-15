from  pyfb.geometry.atom import atom

class step:

  def __init__(self):
    self.atom=[]
    self.line2=""
    self.name=""

  def append(self,atomo):
    self.atom.append(atomo)

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

