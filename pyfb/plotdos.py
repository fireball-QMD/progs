#!/usr/bin/env python

# head dens_001.dat 
#    Energía                                                                                             dens_sum  charge
#   -3.1000    0.0455    0.0123    0.4164    0.0122    0.0000    0.0001    0.0009    0.0001    0.0000    0.4876    0.0244
#   -3.0000    0.0372    0.0094    0.3892    0.0094    0.0000    0.0001    0.0008    0.0001    0.0000    0.4461    0.0711
#   -2.9000    0.0119    0.0050    0.1206    0.0051    0.0000    0.0001    0.0003    0.0001    0.0000    0.1431    0.1005

#plotdos -r dens_001.dat -l atm1 -fermi -3 -yrange -1.5 1.5 -save dos_atm1.png


import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import *

X=[]
Y=[]
Tx=""
Ty=""
L=[]
Legend=False
save=False
outfile="out"
fermi=0.00
ymin=0.0
ymax=0.0
yrange=False
xrange=False
onlyprint=False

for i in range(len(sys.argv)):
  if sys.argv[i] == "-fermi":
    fermi=float(sys.argv[i+1])
  if sys.argv[i] == "-yrange":
    ymin=float(sys.argv[i+1])
    ymax=float(sys.argv[i+2])
    yrange=True
  if sys.argv[i] == "-xrange":
    xmin=float(sys.argv[i+1])
    xmax=float(sys.argv[i+2])
    xrange=True

for i in range(len(sys.argv)):
  if sys.argv[i] == "-r":
    data = loadtxt(sys.argv[i+1])
#    print(len(data[0])) 
    x,y=data[:,len(data[0])-2], data[:,0] #pintamos la penultima columna DensTOT
    #print(x,y)
    X.append(x)  
    Y.append(y-fermi)  
  if sys.argv[i] == "-Tx":
    Tx=sys.argv[i+1]
  if sys.argv[i] == "-Ty":
    Ty=sys.argv[i+1]
  if sys.argv[i] == "-l":
    L.append(sys.argv[i+1])
  if sys.argv[i] == "-save":
    save=True
    outfile=sys.argv[i+1]
  if sys.argv[i] == "-print":
    onlyprint=True


for o in range(len(sys.argv)):
  if sys.argv[o] == "-Y0":
    for i in Y:
      e_min=i[0]
      for j in i:
        if float(e_min) > float(j):
          e_min=float(j)
      print(e_min)
      for j in range(len(i)): #acceder al valor
        i[j]=float(i[j])-float(e_min)
if onlyprint == False:
  fig=plt.figure(figsize=(2,6))
  if len(L) == len(Y):
    Legend=True
  else:
    Legend=False
     
  if yrange:
    plt.ylim((ymin,ymax))
  if xrange:
    plt.xlim((xmin,xmax))
  
  for l in range(len(Y)):
    if Legend:
  #     plt.plot(X[l], Y[l],'-', label=L[l])
      if l == 0:
        plt.plot(X[l], Y[l],'-', label=L[l],linewidth=1.5,color='black')
      else:
        plt.plot(X[l], Y[l],'-', label=L[l],linewidth=1.00)
    else:
      plt.plot(X[l], Y[l],'-', label=l)
  
  
  plt.xlabel(Tx)
  plt.ylabel(Ty)
  plt.xlim(left=0)
  plt.axhline(y=0,color='grey', linewidth=0.50)
  
  if Legend:
    plt.legend(loc='lower right')
  
  if save:
    plt.savefig(outfile)
  else:
    plt.show()
else:
  for i in data:
    a='{0:12.6f}'.format(i[0]-fermi)
    for j in range(1,len(i)):
      a=a+'{0:12.6f}'.format(i[j])
    print(a) 

