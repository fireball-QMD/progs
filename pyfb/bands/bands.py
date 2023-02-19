import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import *

class bands:
  def __init__(self):
    self.data=[]
    self.cristal="null"
    self.titulo="null"
    self.archivo="ek.dat"
    self.savefile=False
    self.fermi=0.00
    self.ajustargap=False
    self.ajustarmin=True
    self.ajustarmax=True
    self.y_min=0.00
    self.y_max=0.00
    self.min_read=0.00
    self.max_read=0.00
    self.ticks=[]
    aux1=[]
    aux2=[]
    self.ticks.append(aux1)
    self.ticks.append(aux2)
    self.printinfo=False
    self.label_e=[] #pintar label de energía
    self.ei=[] #pintar bandas
    self.charge_dens=False
    self.densfile=[]
    self.densfiletitle=[]
    self.up=[]
    self.x=[]
    self.down=[]
    self.paintinfo=False
    self.igapd=0
    self.igapi=0
    self.jgapi=0
    self.AE=0.00

    #def getinfo(self.data,self.fermi,self.AE,self.ticks):
  def getinfo(self):
    #primera banda vacia y ultima llena
    self.igapd
    self.igapi
    self.igapi
    self.data
    self.fermi
    self.AE
    self.ticks
    self.salida=""
    for i in range(len(self.data[:,0])):
      y1=self.data[i,1]
      y2=y1
      for j in range(1,len(self.data[i,:])):
        y=self.data[i][j]
        if y1 < 0:
          if y1 < y:
            y1=y
        if y2 < 0:
          if y2 < y:
            if y < 0:
              y2=y  
      #gap.append([i,y1,y2])
      self.x.append(i+1)
      self.up.append(y1-self.fermi-self.AE)
      self.down.append(y2-self.fermi-self.AE)
    max=0
    med=0
    c=0
    for i in range(1,len(self.down)):
      y1=self.down[i-1]
      y2=self.down[i]
      d=((y1-y2)**2)**0.5
      c=c+1
      med=med*(c-1)/c+d/c
      if max < d :
        max=d
  #  print (max,med)
    metal=0
    if max > med*10 :
      metal=1
    self.salida="-metal "+str(metal)+" "
    c=0
    for i in self.ticks[0]:
      #print(self.ticks[1][c]+" : "+str(self.up[i-1])+"-("+str(self.down[i-1])+")="+str(self.up[i-1]-self.down[i-1]))
      self.salida=self.salida+' -'+self.ticks[1][c]+" "+str(self.up[i-1]-self.down[i-1])+" "
      self.salida=self.salida+' -'+self.ticks[1][c]+"_"+self.ticks[1][2]+" = "+str(self.up[i-1]-self.up[self.ticks[0][2]-1])
      a=self.ticks[1][c]
      c=c+1
    #gap directo
    gap=self.up[0]-self.down[0]
    #self.igapd=0
    
 

    for i in range(len(self.up)):
      if gap > self.up[i]-self.down[i]: 
        gap=self.up[i]-self.down[i]
       # print (self.up[i],self.down[i],gap) 
        self.igapd=i
    self.salida=self.salida+" -gap "+str(gap)
  #  print("Direct gap = ",gap)
    
     #gap indirecto
    gap=self.up[0]-self.down[0]
    self.igapi=0
    self.jgapi=0
    for i in range(len(self.up)):
      for j in range(len(self.up)):
        if gap > self.up[i]-self.down[j]:
          gap=self.up[i]-self.down[j]
          self.igapi=i
          self.jgapi=j
  #  print("Indirect gap = ",gap)
    self.salida=self.salida+" -indirect = "+str(gap)
    print(self.salida)
 
  def getCarga(self,archivo,iemin,iemax):
    datachar=np.loadtxt(archivo)
    qmin=0.00
    qmax=0.00
    first=True
    for i in range(len(datachar)):
        if datachar[i][0] > iemin and datachar[i][0] < iemax:  
            #print(datachar[i][0],datachar[i][len(datachar[i])-1])
            qmax=datachar[i][len(datachar[i])-1]
            if first:
                first=False
                qmin=qmax
    return float(qmax-qmin)

 
  def ajustar(self):
    self.AE
    self.y_max
    self.y_min
    if self.ajustargap==True:
      first=False
      for j in range(1,len(self.data[0,:])):
        for i in self.data[:,j]:
          if i < 0:
            if first==False:
              first=True
              self.AE=i
            else:
              if self.AE < i:
                self.AE=i
  
    self.y_min=self.data[0,1]
    self.y_max=self.data[0,1]
    for i in range(1,len(self.data[0,:])):
      plt.plot(self.data[:,0], self.data[:,i]-self.fermi-self.AE , '-',color='black', linewidth=1.0)
    for i in range(1,len(self.data[0,:])):
      for j in self.data[:,i]:
        if self.y_min > j:
          self.y_min=j
        if self.y_max < j:
          self.y_max=j
    if self.ajustarmax == False :
      self.y_max=self.max_read
    if self.ajustarmin == False :
      self.y_min=self.min_read
    self.y_min=self.y_min-self.fermi-self.AE
    self.y_max=self.y_max-self.fermi-self.AE
    if self.ajustargap and not self.ajustarmax: 
      self.y_max=self.y_max+self.AE
    if self.ajustargap and not self.ajustarmin:
      self.y_min=self.y_min+self.AE
   
  def plot(self):
    fig= plt.figure(figsize=(4,6))
    self.ajustar()
    plt.axhline(y=0,color='grey', linewidth=0.50)
    if self.cristal != "null":
      for i in self.ticks[0]:
        plt.axvline(x=i,color='grey', linewidth=0.50)
    if self.printinfo:
      #getinfo(self.data,self.fermi,self.AE,self.ticks)
      self.getinfo() #self.data,self.fermi,self.AE,self.ticks)
      if self.paintinfo:
        plt.plot([self.x[self.igapd],self.x[self.igapd]],[self.up[self.igapd],self.down[self.igapd]], '-',color='blue', linewidth=1.0) 
        plt.plot([self.x[self.jgapi],self.x[self.igapi]],[self.down[self.jgapi],self.down[self.jgapi]], '-',color='green', linewidth=1.0) 
        plt.plot([self.x[self.igapi],self.x[self.igapi]],[self.down[self.jgapi],self.up[self.igapi]], '-',color='green', linewidth=1.0) 
        plt.plot(self.x,self.up, '-',color='red', linewidth=1.2)
        plt.plot(self.x,self.down, '-',color='orange', linewidth=1.2)

    if len(self.label_e) > 0:
      for inE in self.label_e:
        print("print Energy", inE )
        xl=inE
        #print(self.data[self.nE])
        for il in range(1,len(self.data[inE])):
          yl=self.data[inE][il]-self.fermi-self.AE
          label = "{:.2f}".format(yl)
          #print(xl,yl,label)
          plt.annotate(label, # this is the text
                 (xl,yl), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(2,0), # distance from text to points (x,y)
                 ha='left', # horizontal alignment can be left, right or center
                 size = 6)
    if len(self.ei)>0:
      ie=self.ei[0]
      print(self.ei)
      iemax=self.data[0][ie]-self.fermi-self.AE
      iemin=self.data[0][ie]-self.fermi-self.AE
      for ie in self.ei:
        plt.plot(self.data[:,0], self.data[:,ie]-self.fermi-self.AE , '-',color='blue', linewidth=1.0)
        for j in range(len(self.data[:,ie])):
          if iemax < self.data[j,ie]-self.fermi-self.AE:
            iemax=self.data[j][ie]-self.fermi-self.AE
          if iemin > self.data[j,ie]-self.fermi-self.AE:
            iemin=self.data[j][ie]-self.fermi-self.AE
      plt.axhline(y=iemin,color='grey', linewidth=0.50)
      plt.axhline(y=iemax,color='grey', linewidth=0.50)
      print("iemin = "+str(iemin))
      print("iemax = "+str(iemax))
   
      if self.charge_dens:
         med=(self.y_min+self.y_max)/2
         medh=(self.y_max-self.y_min)
         kk=-medh/24
         for c in self.densfile:
           e1=iemin+self.fermi+self.AE
           e2=iemax+self.fermi+self.AE 
           titulo=self.densfiletitle[self.densfile.index(c)]
           print(titulo,c,e1,e2,self.getCarga(c,e1,e2))
           label = titulo+' '+"{:.2f}".format(self.getCarga(c,e1,e2))
           print(label,((self.y_min+self.y_max)/2+kk))
           plt.annotate(label, # this is the text
                 (150,med+medh*2/4+kk), 
                 textcoords="offset points", # how to position the text
                 xytext=(2,0), # distance from text to points (x,y)
                 ha='left', # horizontal alignment can be left, right or center
                 size = 8)
           kk-=medh/24 
    if self.cristal != "null":
      print(self.ticks)
      plt.xticks(self.ticks[0],self.ticks[1])
    # plt.figure(figsize=(4,3))
    plt.ylim(top=self.y_max)
    print("self.y_min =",self.y_min)
    print("self.y_max =",self.y_max)
    plt.ylim(bottom=self.y_min)
    plt.xlim(left=0)
    plt.xlim(right=len(self.data[:,0]))
    if self.titulo != "null":
      plt.title(self.titulo)
    if self.savefile:
      outfile=self.archivo[0:len(self.archivo)-3]
      #print ("------------")
      print("save file "+outfile+"png")
      plt.savefig(outfile)
    else:
      plt.show()
  
