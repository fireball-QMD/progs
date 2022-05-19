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

fb.start()
