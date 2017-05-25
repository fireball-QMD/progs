
for j in A B C D E F 
do
cd $j
rm -fr ac.dat  answer.bas  answer.xyz NEIGHBORS  NEIGHBORS_PP  param.dat  salida.out  xv.dat
cp charges.bk CHARGES
../fireball.x > salida.out
n=$(head -1 answer.xyz )
ETOT=$(grep ETOT salida.out|tail -1 | tr -s ' ')
sec=$(grep RUN salida.out | cut -d':' -f2 | cut -d'[' -f1)
echo 'n_atom = '$n' : time = '$sec' : '$ETOT >> ../out
cd ..
done
