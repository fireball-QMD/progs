
class tabla:
  def __init__(self):
    self.ele=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc ', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', ' As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', ' La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Tb', 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']

    self.masa=['1.0079','4.0026','6.941','9.0122','10.811','12.0107','14.0067','15.9994','18.9984','20.1797','22.9897','24.305','26.9815','28.0855','30.9738','32.065','35.453','39.0983','39.948','40.078','44.9559','47.867','50.9415','51.9961','54.938','55.845','58.6934','58.9332','63.546','65.39','69.723','72.64','74.9216','78.96','79.904','83.8','85.4678','87.62','88.9059','91.224','92.9064','95.94','98','101.07','102.9055','106.42','107.8682','112.411','114.818','118.71','121.76','126.9045','127.6','131.293','132.9055','137.327','138.9055','140.116','140.9077','144.24','145','150.36','151.964','157.25','158.9253','162.5','164.9303','167.259','168.9342','173.04','174.967','178.49','180.9479','183.84','186.207','190.23','192.217','195.078','196.9665','200.59','204.3833','207.2','208.9804','209','210','222','223','226','227','231.0359','232.0381','237','238.0289','243','244','247','247','251','252','257','258','259','261','262','262','264','266','268','277']


  def getN(self,Z):
    s=1
    a=0
    for j in self.ele :
      if Z == j :
        a=s
      s=s+1
    return a
  
  def getM(self,Z):
    s=1
    for j in self.ele :
      if Z == j :
        a=s
      s=s+1
    return masa[int(a)]

