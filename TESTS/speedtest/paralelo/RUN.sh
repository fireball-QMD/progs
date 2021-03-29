#MPI + SCALAPACK
for j in A B C D E F 
do
for i in 01 02 04 08 16 
do
cd $j$i
rm -fr ac.dat  answer.bas  answer.xyz NEIGHBORS  NEIGHBORS_PP  param.dat  salida.out  xv.dat
cp charges.bk CHARGES
mpirun  -np $i ../fireball.x > salida.out
n=$(head -1 answer.xyz )
ETOT=$(grep ETOT salida.out|tail -1 | tr -s ' ')
sec=$(grep RUN salida.out | cut -d':' -f2 | cut -d'[' -f1 | head -1)
echo 'np = '$i'  : n_atom = '$n' : time = '$sec' : '$ETOT >> ../out
cd ..
done
done
