### GPLv3 ###

 
##----- Parametros de control (el parametro de red tiene que encontrase entre ini fin  --------##
iqout=${1:-1} 
a0=${2:-1} 
b0=${3:-1} 
Z=${4:-1}
c0=${5:-1}
cristal=${6:-1}

N=12
ini=$(python -c "print '%.6f' % ($a0-0.60)")
fin=$(python -c "print '%.6f' % ($a0+0.60)")
echo $a0 $b0 $ini $fin > log.$$
echo "           1
  $Z      0.000000      0.000000      0.000000" > uno.bas

#sobre escribo el uno.bas en el caso de ser diamante
echo $cristal >> log.$$
if [ $cristal == dia ]
then
echo "           2
  $Z      0.000000      0.000000      0.000000
  $Z      0.250000      0.250000      0.250000" > uno.bas
fi


##----------Funcion analisis para un atomos/celda-----------------------------------------##
function analisis {
ETOT=$(grep 'ETOT' salida.out|cut -d'=' -f2)
sigma=$(grep sigma salida.out | cut -d'=' -f16 | tail -1)
charge=$(head -2 uno.bas | tail -1 | tr -s ' ' | cut -d' ' -f1)' -> '$(head -2 CHARGES | tail -1)
if [[ $((sigma)) != 200 ]]; 
then
if [[ $sigma != '' ]];
then
#echo $rescal$'\t'$ETOT$'\t'$sigma$'\t'$charge$'\t'>>salida
echo $rescal$'\t'$ETOT>>salida
fi
fi
}
##------------------------------------------------------------------------------------------##
function start {
rm -fr salida
for((i=0;i<=N;i++))
do
rescal=$(python -c "print '%.6f' % ($i*1.0*($fin-$ini)/$N+$ini)")
echo "&option
basisfile = uno.bas
lvsfile = uno.lvs
kptpreference = uno.kpts
rescal = $rescal
iqout = $iqout
sigmatol = 0.000001
nstepf = 1
&end
&output
iwrtxyz = 1
&end" > fireball.in 
./fireball.x > salida.out    

analisis
done
}
##----------------------------------1º start  ----------------------------------------------##
start
##-----------------------------buscamos el minimo ------------------------------------------##
min=$(python -c "
x0 = []
y0 = []
for line in file(\"salida\"):
   line = line.split()
   x = line[0]
   y = line[1]
   x0.append(x)
   y0.append(y)  

j=0
for i in range(len(x0)):
  if float(y0[j]) > float(y0[i]):
    j=i
print x0[j]")
cp salida borrar
rm -fr aux.py		        

##----------------------------------2º start -----------------------------------------------##
N=20
d=0.2
ini=$(python -c "print '%.6f' % (1.0*($min-$d))")
fin=$(python -c "print '%.6f' % (1.0*($min+$d))")
start
mv salida Vol.dat
##------------------------------------------------------------------------------------------## 


#      print \"%4.0f %12.6f %12.6f %12.6f' % (zz, x , y , z)

echo $(python -c "
x0 = []
y0 = []
for line in file(\"Vol.dat\"):
   line = line.split()
   x = line[0]
   y = line[1]
   x0.append(x)
   y0.append(y)  

j=0
for i in range(len(x0)):
  if float(y0[j]) > float(y0[i]):
    j=i
f=open(\"cero.dat\",\"w\")
for i in range(len(x0)):
  z=float(y0[i])-float(y0[j])
  f.write(x0[i]+\" \"+str(z)+\"\\n\")
")


echo $(python -c "
x0 = []
y0 = []
for line in file(\"Vol.dat\"):
   line = line.split()
   x = line[0]
   y = line[1]
   x0.append(x)
   y0.append(y)  

min=0
for i in range(len(x0)):
  if float(y0[min]) > float(y0[i]):
    min=i

j=0
f=open(\"min.bas\",\"w\")
for line in file(\"uno.bas\"):
   line = line.split()
   if j==0:
      f.write(line[0]+\"\\n\")
   else:
      zz= float(line[0])
      x = float(x0[min])*float(line[1])
      y = float(x0[min])*float(line[2])
      z = float(x0[min])*float(line[3])
      f.write(\"%4.0f %12.6f %12.6f %12.6f\" % (zz, x, y, z)+\"\\n\")
   j=j+1

f=open(\"min.lvs\",\"w\")
for line in file(\"uno.lvs\"):
   line = line.split()
   x = float(x0[min])*float(line[0])
   y = float(x0[min])*float(line[1])
   z = float(x0[min])*float(line[2])
   f.write(\"%12.6f %12.6f %12.6f\" % (x, y, z)+\"\\n\")

j=0
f=open(\"min.kpts\",\"w\")
for line in file(\"uno.kpts\"):
   line = line.split()
   if j==0:
      f.write(line[0]+\"\\n\")
   else:
      x = float(line[0])/float(x0[min])
      y = float(line[1])/float(x0[min])
      z = float(line[2])/float(x0[min])
      t = float(line[3])
      f.write(\"%12.6f %12.6f %12.6f %12.6f\" % (x, y, z, t)+\"\\n\")
   j=j+1
")




minimos=$(python -c "
x0 = []
y0 = []
for line in file(\"cero.dat\"):
   line = line.split()
   x = line[0]
   y = line[1]
   x0.append(x)
   y0.append(y)  

ymin=0
ymax=0
xmax=0
xmin=0
for i in range(len(x0)):
  if float(y0[ymin]) > float(y0[i]):
    ymin=i
  if float(y0[ymax]) < float(y0[i]):
    ymax=i
  if float(x0[xmin]) > float(x0[i]):
    xmin=i
  if float(x0[xmax]) < float(x0[i]):
    xmax=i
   

print str(x0[xmin])+', '+str(y0[ymin])+', '+str(x0[xmax])+', '+str(y0[ymax])")




salida=$(/home/dani/bin/xeo/xeo -Bulk Vol.dat file min.lvs)
if [[ salida != '' ]]
then
a0s=$(echo $salida | cut -d' ' -f1 | cut -c-5)
b0s=$(echo $salida | cut -d' ' -f3 | cut -c-5)
t1='a0 = '$a0s' ; b0 = '$b0s
t2='exp. a0 = '$a0' ; b0 = '$b0
sed "s/title\ \"title\"/title\ \"$t1\\\n $t2\"/" bulk.agr | sed "s/WORLD/world\ $minimos/"  > Vol.agr
fi

cat cero.dat >> Vol.agr
echo '&' >>  Vol.agr
 


