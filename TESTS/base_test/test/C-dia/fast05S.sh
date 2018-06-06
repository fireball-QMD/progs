


e=6 ;Z[${e}]=006; X[${e}]=C ; M[${e}]=12.011; orb[${e}]=sp ;  C[${e}]=dia ; a0[${e}]=3.567; B0[${e}]=442

Z=6
iqout=$1
create=$2
name=$3

cristal=${C[e]}
a0=${a0[e]}
b0=${B0[e]}
c0=1 #hcp




cd bulk
ln -s ../../Fdata .
ln -s ../../fireball.x .
./fast05S.sh $iqout $a0 $b0 $Z $c0 $cristal
if [[ $iqout == 1 ]]; then  d="../../../info-Lowdin/";fi
if [[ $iqout == 3 ]]; then  d="../../../info-NAP/";fi
cp min.lvs $d/${cristal}_${create}.lvs
cp min.bas $d/${cristal}_${create}.bas
cp min.kpts $d/${cristal}_${create}.kpts
cp Vol.dat $d/bulk_${cristal}_${create}.dat
cp Vol.agr $d/bulk_${cristal}_${create}.agr

if [[ $cristal == hcp ]]
then
min_c=$(python -c "
x0 = []
y0 = []
for line in file(\"Vol.dat\"):
   line = line.split()
   x = line[2]
   y = line[1]
   x0.append(x)
   y0.append(y)  
j=0
for i in range(len(x0)):
  if float(y0[j])< float(y0[i]):
    j=i
print x0[j]")
fi


min=$(python -c "
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
  if float(y0[j])< float(y0[i]):
    j=i
print x0[j]")

cd ..
cd dos-bandas
ln -s ../../Fdata .
ln -s ../../fireball.x .
./fast05S.sh $iqout $min $Z $min_c

cp ek0.dat  $d/ek_${create}.dat
cp ek.agr $d/ek_${create}.agr

cd ..

if ! test -d ${d}FIX; then mkdir ${d}FIX; fi
head -$((natoms+2))  answer.xyz  >${d}FIX/${name}_${create}.xyz
