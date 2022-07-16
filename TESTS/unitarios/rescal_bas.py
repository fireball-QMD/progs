#!/usr/bin/env python
import sys


print_first_line=True
for i in sys.argv:
  if i == '-without_first_line':
    print_first_line=False

archivo=sys.argv[1]
rescal=sys.argv[2]

j=0

for line in open(archivo):
   line = line.split()
   if j==0:
      if print_first_line:
        print(line[0])
   else:
      zz= line[0]
      x = float(line[1])*float(rescal)
      y = float(line[2])*float(rescal)
      z = float(line[3])*float(rescal)
      print("%s %12.6f %12.6f %12.6f" % (zz, x, y, z))
   j=j+1

