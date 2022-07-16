#!/usr/bin/env python
import sys

archivo=sys.argv[1]
rescal=sys.argv[2]

for line in open(archivo):
   line = line.split()
   x = float(rescal)*float(line[0])
   y = float(rescal)*float(line[1])
   z = float(rescal)*float(line[2])
   print("%12.6f %12.6f %12.6f" % (x, y, z))
