from pyfb.geometry.tabla import tabla
import math

tabla=tabla()
class atom:
    def __init__(self,Z='H',r=['0.00', '0.00', '0.00']):
        self.r = []
        self.v = []
        self.Z = Z
        self.r = r
        self.Q = "null"
        self.q = []

    def print(self):
        print ('{0:2} {1:12.6f} {2:12.6f} {3:12.6f}'.format(self.Z,float(self.r[0]),float(self.r[1]),float(self.r[2])))

    def print_charges(self):
        a='{0:2} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f}'.format(self.Z,float(self.r[0]),float(self.r[1]),float(self.r[2]),float(self.Q))
        for i in self.q:
           a=a+'{0:12.6f}'.format(float(i))
        print(a)

    def print_bas_format(self):
        print ('{0:2} {1:12.6f} {2:12.6f} {3:12.6f}'.format(tabla.getN('H'),float(self.r[0]),float(self.r[1]),float(self.r[2])))

    def setV(self , v):
        self.v = v

    def setZ(self , Z):
        self.Z = Z

    def setR(self , r):
        self.r = r

    def V2(self):
        return float(self.v[0]**2+self.v[1]**2+self.v[2]**2)

    def getM(self):
        return float(getM(self.Z))

    def getEc(self):
        return str(0.5*(self.getM()*self.V2()/fovermp))

    def distancia(self,atomo2):
        d=((float(self.r[0])-float(atomo2.r[0]))**2+(float(self.r[1])-float(atomo2.r[1]))**2+(float(self.r[2])-float(atomo2.r[2]))**2)**0.5
        return d

    def diffQ(self,atomo2):
        d=(float(self.Q)-float(atomo2.Q))
        return d

    def ang(self,atomo1,atomo2):
        x0 = float(self.r[0])
        y0 = float(self.r[1])
        z0 = float(self.r[2])
        x1 = float(atomo1.r[0])
        y1 = float(atomo1.r[1])
        z1 = float(atomo1.r[2])
        x2 = float(atomo2.r[0])
        y2 = float(atomo2.r[1])
        z2 = float(atomo2.r[2])
        ax=x0-x1;
        ay=y0-y1;
        az=z0-z1;
        bx=x2-x1;
        by=y2-y1;
        bz=z2-z1;
        cos=(ax*bx+ay*by+az*bz)/(ax*ax+ay*ay+az*az)**0.5/(bx*bx+by*by+bz*bz)**0.5;
        d=math.acos(cos)*180/math.pi;
        return d

