### GPLv3 ###
 
##----- Parametros de control (el parametro de red tiene que encontrase entre ini fin  --------##
iqout=${1:-1}
rescal=${2:-1} #min 
Z=${3:-1}


cp ../bulk/uno.bas .
cp ../bulk/uno.lvs .
cp ../bulk/uno.kpts .

#------------------ damos un paso -----------------
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

fermi=$(grep "Fermi Level" salida.out | tail -1 | cut -d'=' -f2)
f=$(echo "$fermi-3-0.1" | bc -l)

echo "1.0                   ! scale factor (leave 1.0)
1        1            ! list of atoms to analyze DOS
100                   ! number of energy steps
$f 0.1
0                     ! leave untouched
0.0     0.0           ! leave untouched
0.05                  ! imaginary part of Green function (controls energy level smearing)" > dos.optional


echo "&option
basisfile = uno.bas
lvsfile = uno.lvs
kptpreference = camino.kpts
rescal = $rescal
iqout = $iqout
sigmatol = 0.000001
nstepf = 1
ifixcharge = 1
&end
&output
iwrteigen = 1
iwrtxyz = 1
iwrtdos = 1
&end" > fireball.in 


./fireball.x >> salida.out    



cp bandas.agr ek.agr
archivo=ek.dat

ncol=$(python -c "
i=0
f = open("'"ekf.dat"'","'"w"'")
for line in file('$archivo'):
   line = line.split()
   out = []
   c=0
   for col in line:
      c=c+1
      if(c>1):
         out.append(float(col)-$fermi)
      else:
         out.append(col)
   f.write(str(out)+"'"\n"'")
   print c
")

ncol=$(echo $ncol | cut -d' ' -f1)

cat ekf.dat | sed 's/\[//g' | sed 's/\]//g' |sed 's/\,//g'| sed "s/'//g" | grep . > ek0.dat
rm -fr ekf.dat
for((i=0;i<$ncol;i++))
do
echo "@    s$i hidden false
@    s$i type xy
@    s$i symbol 0
@    s$i symbol size 1.000000
@    s$i symbol color 1
@    s$i symbol pattern 1
@    s$i symbol fill color 1
@    s$i symbol fill pattern 0
@    s$i symbol linewidth 1.0
@    s$i symbol linestyle 1
@    s$i symbol char 65
@    s$i symbol char font 0
@    s$i symbol skip 0
@    s$i line type 1
@    s$i line linestyle 1
@    s$i line linewidth 2.0
@    s$i line color 1
@    s$i line pattern 1
@    s$i baseline type 0
@    s$i baseline off
@    s$i dropline off
@    s$i fill type 0
@    s$i fill rule 0
@    s$i fill color 1
@    s$i fill pattern 1
@    s$i avalue off
@    s$i avalue type 2
@    s$i avalue char size 1.000000
@    s$i avalue font 0
@    s$i avalue color 1
@    s$i avalue rot 0
@    s$i avalue format general
@    s$i avalue prec 3
@    s$i avalue prepend \"\"
@    s$i avalue append \"\"
@    s$i avalue offset 0.000000 , 0.000000
@    s$i errorbar on
@    s$i errorbar place both
@    s$i errorbar color 1
@    s$i errorbar pattern 1
@    s$i errorbar size 1.000000
@    s$i errorbar linewidth 1.0
@    s$i errorbar linestyle 1
@    s$i errorbar riser linewidth 1.0
@    s$i errorbar riser linestyle 1
@    s$i errorbar riser clip off
@    s$i errorbar riser clip length 0.100000
@    s$i comment \"\"
@    s$i legend  \"\"
" >> ek.agr
done

for((i=2;i<=$ncol;i++))
do
echo "
@target G0.S$(($i-2))
@type xy" >> ek.agr
cut -d' ' -f1,$i ek0.dat >> ek.agr
echo "&" >> ek.agr
done


grep -v @ ek.agr | grep -v '#' | grep -v '&' | grep '.'  > max_min.dat

minimos=$(python -c "
x0 = []
y0 = []
for line in file(\"max_min.dat\"):
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
   

print str(0.0)+', '+str(y0[ymin])+', '+str(x0[xmax])+', '+str(y0[ymax])")


sed "s/WORLD/world\ $minimos/g" ek.agr > z
mv z ek.agr


