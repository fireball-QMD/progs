#!/usr/bin/env python

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.geometry.dinamic import *

din=dinamic()
acumular_salida=False
sep="                                                                                                            "
def print_help() :
  print (sys.argv[0], '-i <file xyz format>')
  print (sys.argv[0], '-i <file xyz format> -print')
  print (sys.argv[0], '-i <file xyz format> -op <atom position>')
  print (sep[:len(sys.argv[0])], '                     -op = -x -y -z -X -Y -X') 
  print (sys.argv[0], '-i <file xyz format> -op <atom1 position> <atom2 position>' ) 
  print (sep[:len(sys.argv[0])],'                     -op = -d -D' )
  print (sys.argv[0], '-i <file xyz format> -op <atom1 position> <float>' ) 
  print (sep[:len(sys.argv[0])], '                     -dx -dy -dz'  )
  print (sys.argv[0], '-i <file xyz format> -op <atom1 position> <atom2 position> <atom3 position>')
  print (sep[:len(sys.argv[0])],'                     -op = -ang -ANG ' )
  print (sys.argv[0], '-loadstep <file xyz format> <step> -print , -print_bas_format , -loadstep ...')
  print (sys.argv[0], '-laststep <file xyz format> -print')
  print (sys.argv[0], '-laststep_charges <file xyz format> -print_charges')
  print (sys.argv[0], '-i <file xyz format> -dintra_matrix')

if len(sys.argv) == 1 :
  print_help()

for i in range(1,len(sys.argv)):
  if sys.argv[i] == '-i' :
    din.loadxyz(sys.argv[i+1])

  if sys.argv[i] == '-loadstep' :
    din.loadstep(sys.argv[i+1],int(sys.argv[i+2]))

  if sys.argv[i] == '-laststep' :
    din.laststep(sys.argv[i+1],False)

  if sys.argv[i] == '-laststep_charges' :
    din.laststep(sys.argv[i+1],True)

  if sys.argv[i] == '-print' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.print()

  if sys.argv[i] == '-print_bas_format' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.print_bas_format()

  if sys.argv[i] == '-print_charges' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.print_charges()

  if sys.argv[i] == '-dintra_matrix' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.dintra_matrix()


  if sys.argv[i] == '-x' or sys.argv[i] == '-y' or sys.argv[i] == '-z' or sys.argv[i] == '-X' or sys.argv[i] == '-Y' or sys.argv[i] == '-Z':
    col=[int(sys.argv[i+1])]
    din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-d' or sys.argv[i] == '-D' or sys.argv[i] == '-dx' or sys.argv[i] == '-dy' or sys.argv[i] == '-dz':
    col=[int(sys.argv[i+1]),int(sys.argv[i+2])]
    din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-ang' or sys.argv[i] == '-ANG':
    col=[int(sys.argv[i+1]),int(sys.argv[i+2]),int(sys.argv[i+3])]
    din.get(sys.argv[i],col)
    acumular_salida=True

if acumular_salida:
  din.print_out()
