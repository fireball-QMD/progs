#!/usr/bin/env python3

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.geometry.dinamic import *

din=dinamic()
acumular_salida=False
sep="                                                                                                            "
def print_help() :
  print (sys.argv[0], '-i <file xyz format>')
  print (sys.argv[0], '-ibas <file bas format>')
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
  print (sys.argv[0], '-loadstep <file xyz format> <step> -enlaces -print_enlaces')
  print (sys.argv[0], '-loadstep <file xyz format> <step> -loadstep <file xyz format> <step> -diff_enlaces')
  print (sys.argv[0], '-i <file xyz format> -rescal <float>')
  print (sys.argv[0], '-i <file xyz format> -center <step>')
  print (sys.argv[0], '-i <file xyz format> -rotX <step> <float>')
  print (sep[:len(sys.argv[0])], '                     -rotX -rotY'  )
  print (sys.argv[0], '-laststep2xyz <file xyz format> <file xyz format> -merge')
  print (sys.argv[0], '-laststep2xyz <file xyz format> <file xyz format> -distance')
  print (sys.argv[0], '-laststep2xyz <file xyz format> <file xyz format> -join')
  print (sys.argv[0], '-laststep2xyz <file xyz format> <file xyz format> -joinrnd')
  print (sys.argv[0], '-i <file xyz format> -lvs <file lattice vectors> -R 2 1 1')
  print (sys.argv[0], '-i <file xyz format> -lvs <file lattice vectors> -unit')
  print (sys.argv[0], '-i <file xyz format> -lvs <file lattice vectors> -unit_index <iatom_max>')

if len(sys.argv) == 1 :
  print_help()

for i in range(1,len(sys.argv)):
  if sys.argv[i] == '-i' :
    din.loadxyz(sys.argv[i+1])

  if sys.argv[i] == '-ibas' :
    din.loadbas(sys.argv[i+1])

  if sys.argv[i] == '-laststep2xyz' :
    din.laststep(sys.argv[i+1],False)
    din.laststep(sys.argv[i+2],False)
    if sys.argv[i+3] == '-merge' :
      din.merge(1,2)
    if sys.argv[i+3] == '-distance' :
      print(din.distance(1,2))
    if sys.argv[i+3] == '-join' :
      din.join(1,2)
    if sys.argv[i+3] == '-joinrnd' :
      din.joinrnd(1,2)

  if sys.argv[i] == '-loadstep' :
    din.loadstep(sys.argv[i+1],int(sys.argv[i+2]))

  if sys.argv[i] == '-laststep' :
    din.laststep(sys.argv[i+1],False)

  if sys.argv[i] == '-laststep_charges' :
    din.laststep(sys.argv[i+1],True)

  if sys.argv[i] == '-center' :
    din.center(int(sys.argv[i+1]))
  
  if sys.argv[i] == '-rotX' :
    din.rotX(int(sys.argv[i+1]),float(sys.argv[i+2]))

  if sys.argv[i] == '-rotY' :
    din.rotY(int(sys.argv[i+1]),float(sys.argv[i+2]))

  if sys.argv[i] == '-rotZ' :
    din.rotZ(int(sys.argv[i+1]),float(sys.argv[i+2]))

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

  if sys.argv[i] == '-enlaces' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.load_enlaces()
  
  if sys.argv[i] == '-print_enlaces' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.print_enlaces()

  diff_enlaces=[]
  if sys.argv[i] == '-diff_enlaces' :
    if len(din.step) == 0:
      print("nothing is load")
    else:
      din.load_enlaces()
      for bas in din.step: #solo debería haber un step cargado
        diff_enlaces.append(bas.enlaces) 
    dis=0.0
    if (len(diff_enlaces[0]) == len(diff_enlaces[1])):
      for ienlace in range(len(diff_enlaces[0])):
        dis=dis+(diff_enlaces[0][ienlace]-diff_enlaces[1][ienlace])**2     
      print(dis**0.5/len(diff_enlaces[0]))
    else:
      print("El número de enlaces ha cambiado") 
  
  if sys.argv[i] == '-x' or sys.argv[i] == '-y' or sys.argv[i] == '-z' or sys.argv[i] == '-X' or sys.argv[i] == '-Y' or sys.argv[i] == '-Z' or sys.argv[i] == '-rescal':
    if sys.argv[i] == '-rescal':
       din.get(sys.argv[i],[sys.argv[i+1]])
    else:
      col=[int(sys.argv[i+1])]
      din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-d' or sys.argv[i] == '-D':
    col=[int(sys.argv[i+1]),int(sys.argv[i+2])]
    din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-dx' or sys.argv[i] == '-dy' or sys.argv[i] == '-dz':
    col=[int(sys.argv[i+1]),float(sys.argv[i+2])]
    din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-ang' or sys.argv[i] == '-ANG':
    col=[int(sys.argv[i+1]),int(sys.argv[i+2]),int(sys.argv[i+3])]
    din.get(sys.argv[i],col)
    acumular_salida=True

  if sys.argv[i] == '-lvs' :
    din.loadlvs(sys.argv[i+1])

  if sys.argv[i] == '-R' :
    din.repitelvs(int(sys.argv[i+1]),int(sys.argv[i+2]),int(sys.argv[i+3]))

  if sys.argv[i] == '-unit' :
    din.unitlvs_iatom_max(0)

  if sys.argv[i] == '-unit_index' :
    din.unitlvs_iatom_max(int(sys.argv[i+1]))

if acumular_salida:
  din.print_out()
