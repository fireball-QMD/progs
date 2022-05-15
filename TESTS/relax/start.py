#!/opt/intel/oneapi/intelpython/python3.7/bin/python3.7
import sys
import os
sys.path.append(os.environ["FIREBALLHOME"])
print(sys.path)
import libpyfb as fb

#Si no existe Fdata lo extraemos'

import tarfile
from os.path import exists 

if not exists('Fdata_HC_minimal'):
  file = tarfile.open('../Fdata.tar.gz')                                       
  file.extractall('.')   
  file.close()


if not exists('fireball.in'):
  f = open("fireball.in", "a")
  firein='&OPTION'+"\n"+\
  'fdatalocation = "Fdata_HC_minimal"'+"\n"+\
  'nstepf = 5000'+"\n"+\
  'iquench = -1'+"\n"+\
  'icluster = 1'+"\n"+\
  'iqout = 3'+"\n"+\
  'dt = 0.5'+"\n"+\
  '&END'+"\n"+\
  '&OUTPUT'+"\n"+\
  'iwrtxyz = 1'+"\n"+\
  '&END'

  f.write(firein) 
  f.close()


if not exists('input.bas'):
  f = open("input.bas", "a")
  bas='12'+"\n"+\
  '6     -4.197174     -3.920062     -0.001077'+"\n"+\
  '6     -3.132111     -3.058028      0.002942'+"\n"+\
  '6     -3.346311     -1.704680      0.002885'+"\n"+\
  '6     -4.625536     -1.213686      0.001660'+"\n"+\
  '6     -5.690602     -2.075731     -0.001805'+"\n"+\
  '6     -5.476402     -3.429078     -0.004659'+"\n"+\
  '1     -6.310509     -4.036755     -0.008219'+"\n"+\
  '1     -3.985726     -5.256165     -0.001792'+"\n"+\
  '1     -2.242726     -3.399675      0.005306'+"\n"+\
  '1     -2.295045     -0.853275      0.005357'+"\n"+\
  '1     -4.774460     -0.272773      0.003081'+"\n"+\
  '1     -6.953318     -1.590684     -0.003133'
  f.write(bas)
  f.close()

fb.start()
