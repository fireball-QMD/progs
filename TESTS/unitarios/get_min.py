#!/usr/bin/python3
import sys

archivo=sys.argv[1]

x0 = []
y0 = []
for line in open(archivo):
   line = line.split()
   if len(line) > 1 :
     x = line[0]
     y = line[1]
     x0.append(x)
     y0.append(y)  

j=0
for i in range(len(x0)):
  if float(y0[j]) > float(y0[i]):
    j=i

print (x0[j])

