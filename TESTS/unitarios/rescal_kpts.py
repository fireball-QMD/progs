#!/usr/bin/env python
import sys

archivo=sys.argv[1]
rescal=sys.argv[2]
j=0
for line in open(archivo):
   line = line.split()
   if len(line) > 0 :
     if j==0:
        print(line[0])
     else:
        x = float(line[0])/float(rescal)
        y = float(line[1])/float(rescal)
        z = float(line[2])/float(rescal)
        t = float(line[3])
        print("%12.6f %12.6f %12.6f %12.6f" % (x, y, z, t))
     j=j+1

