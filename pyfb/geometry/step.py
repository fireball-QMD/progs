from  pyfb.geometry.atom import atom

class step:

  def __init__(self):
    self.atom=[]
    self.charges=[]
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

  def loadcharges(self):
    txt=open("CHARGES").readlines()
    for i in range(1,len(txt)):
      charge=0
      line=txt[i].split()
      for j in range(len(line)):
        self.atom[i-1].q.append(line[j])
        charge=charge+float(line[j])
      self.atom[i-1].Q=charge
     
