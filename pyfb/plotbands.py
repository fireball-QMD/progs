#!/usr/bin/env python

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.bands.bands import *

bands=bands()

def help():                                                                                                                 
  print(sys.argv[0]+" -c fcc")                                                                                              
  print('plotbands.py -c bcc -r ek.dat -save -plot')                                                                        
  print('plotbands.py -c bcc -r ek.dat -plot')                                                                              
  print('plotbands.py -c bcc -r ek.dat -titulo titulo -plot')                                                               
  print('plotbands.py -c bcc -r ek0_C.pbe-van_ak.UPF.dat -plot -ajustargap')                                                
  print('plotbands.py -r ek0_Si.pbesol-n-rrkjus_psl.0.1.UPF.dat -min -10 -max 2.5 -ajustargap -plot')                       
  print('plotbands.py -r ek0_Si.pbesol-n-rrkjus_psl.0.1.UPF.dat -min -10 -max 2.5 -info -ajustargap -plot')                 
  print('plotbands.py -r ek0_Si.pbesol-n-rrkjus_psl.0.1.UPF.dat -ajustargap  -info  -plot')                                 
  print('plotbands.py -r ek0_Si.pbesol-n-rrkjus_psl.0.1.UPF.dat -ajustargap  -info -paintinfo -plot')                       
  print('plotbands.py -c bcc -r ek0_As.pbe-n-van.UPF.dat -onlyinfo')                                                        
  print('plotbands.py -c bcc -r ek0_As.pbe-n-van.UPF.dat -print_E 0')                                                        
  print('plotbands.py -c bcc -r ek0_As.pbe-n-van.UPF.dat -ajustargap -getMAX')


if len(sys.argv) == 1 :
  help()

for i in range(1,len(sys.argv)):
  if sys.argv[i] == '-h' :
    help()
  if sys.argv[i] == '-titulo':
   bands.titulo=sys.argv[i+1]
  if sys.argv[i] == '-c' :
    bands.cristal=sys.argv[i+1]
    if bands.cristal == "fcc" or bands.cristal == "dia":
      bands.ticks=[[1,52,116,190,227,255],["W","L","G","X","W","K"]]
    if bands.cristal == "my":
      bands.ticks=[[1,52,82,142,202,255],["G","M","X","X","G","M"]]
    if bands.cristal == "my2":
      bands.ticks=[[1,107,160,255],["G","K","M","G"]]
    if bands.cristal == "bcc" :
      bands.ticks=[[1,77,132,185,255],["G","H","N","G","P"]]
    if bands.cristal == "sc" :
      bands.ticks=[[1,86,136,184,255],["R","G","X","M","G"]]
  if sys.argv[i] == '-r' :
    bands.archivo=sys.argv[i+1]
    bands.data=loadtxt(bands.archivo) 
  if sys.argv[i] == '-fermi':
    bands.fermi=float(sys.argv[i+1])
  if sys.argv[i] == '-min':
    bands.ajustarmin=False
    bands.min_read=float(sys.argv[i+1])
    print(min)
  if sys.argv[i] == '-max':
    bands.ajustarmax=False
    bands.max_read=float(sys.argv[i+1])
  if sys.argv[i] == '-info':
    bands.printinfo=True
  if sys.argv[i] == '-print_E':
    bands.print_E=True
    bands.nE=int(sys.argv[i+1])
  if sys.argv[i] == '-onlyinfo':
    bands.getinfo()
  if sys.argv[i] == '-paintinfo':
    bands.printinfo=True
    bands.paintinfo=True

for i in range(1,len(sys.argv)):
  if sys.argv[i] == '-plot' :
   bands. plot()
  if sys.argv[i] == '-save':
    bands.savefile=True
    bands.plot()
  if sys.argv[i] == '-print' :
    aux=""
    for k in range(len(bands.data[:,0])):
      aux=str(int(bands.data[k,0]))
      for j in range(1,len(bands.data[k,:])):
        aux=aux+" "+str(bands.data[k,j]-bands.fermi-bands.AE)
      print(aux)
  if sys.argv[i] == '-getMAX':
    print("-min "+str(bands.y_min)+" -max "+str(bands.y_max))


if len(sys.argv) == 2:
  if sys.argv[i] != '-h' :
    bands.cristal="null"
    bands.archivo=sys.argv[1]
    bands.data=loadtxt(bands.archivo)
    bands.plot()

