#!/usr/bin/env python

import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])

from pyfb.geometry.dinamic import *

din=dinamic()

def print_help() :
  print (sys.argv[0], '-i <file xyz format>')
  print (sys.argv[0], '-i <file xyz format> -print')
  print (sys.argv[0], '-i <file xyz format> -x <atom position>'  )

if len(sys.argv) == 1 :
  print_help()

for i in range(1,len(sys.argv)):
  if sys.argv[i] == '-i' :
    din.load_xyz(sys.argv[i+1])

  if sys.argv[i] == '-print' :
    if len(din.structure) == 0:
      print("nothing is load")
    else:
      din.print()

  if sys.argv[i] == '-x' :
    col=[["x",int(sys.argv[i+1])]]
    din.get(col)
    din.print_out()
