#!/usr/bin/python

import sys

x0 = []
y0 = []
i=0
suma=0
for line in file(sys.argv[1]):
   line = line.split()
   if len(line)>2:
      #print len(line)
      x = line[1]
      y = line[2]
      z=abs(float(x)-float(y))
      i=i+1
      suma=suma+z
   #rint x,y,z,i,suma
   #print suma/i
print suma/i

