#!/usr/bin/python

import sys

y0 = []
suma=0
i=0
for line in file(sys.argv[1]):
   line = line.split()
   if len(line)>2:
      y = line[3]
      i=i+1
      suma=suma+float(y)
print suma/i

