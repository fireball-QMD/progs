#./fast05S.sh $iqout  ${basedir} ${create} ${name} ${idftd3} ${dftd3_func}

unitarios="$FIREBALLHOME/TESTS/unitarios/"
fdata=$FIREBALLHOME/TESTS/speedtest/Fdata/

Z=14
name=Si
cristal=dia
a0=5.43

iqout=4

sed "s/Z/${Z}/" ${unitarios}/${cristal}.bas > uno.bas
cp ${unitarios}/${cristal}.lvs uno.lvs
cp ${unitarios}/${cristal}.kpts uno.kpts
cp ${unitarios}/${cristal}_*.kpts camino.kpts


function start {
echo "&option
fdatalocation = '$fdata'
basisfile = 'uno.bas'
lvsfile = 'uno.lvs'
kptpreference = 'uno.kpts'
iqout = $iqout
rescal = $1
iquench = -1
sigmatol = 0.000001
nstepf = 1
&end
&output
iwrtxyz = 1
&end" > fireball.in
#$FIREBALLHOME/fireball.x > salida.out
/home/dani/FB/git/bin/fireball.21.09_intel-double.x > salida.out 
ETOT=$(grep 'ETOT' salida.out|cut -d'=' -f2)
echo $1$'\t'$ETOT>>salida
}



function analisis {
sigma=$(grep sigma salida.out | cut -d'=' -f16 | tail -1)
charge=$(head -2 uno.bas | tail -1 | tr -s ' ' | cut -d' ' -f1)' -> '$(head -2 CHARGES | tail -1)
if [[ $((sigma)) != 200 ]]; 
then
if [[ $sigma != '' ]];
then
echo $rescal$'\t'$ETOT$'\t'$sigma$'\t'$charge$'\t'>>salida
fi
fi
}

##----------------------------------1º start  ----------------------------------------------##
  N=12
  ini=$(python3 -c "print ('%.6f' % ($a0-0.60))")
  fin=$(python3 -c "print ('%.6f' % ($a0+0.60))")
  echo $a0 $b0 $ini $fin > log.$$
  echo $cristal >> log.$$

  rm -fr salida
  for((i=0;i<=N;i++))
  do
    rescal=$(python3 -c "print ('%.6f' % ($i*1.0*($fin-$ini)/$N+$ini))")
    start $rescal
  done

  min=$($FIREBALLHOME/TESTS/unitarios/get_min.py salida)
  cp salida borrar

  N=20
  d=0.2
  ini=$(python3 -c "print ('%.6f' % (1.0*($min-$d)))")
  fin=$(python3 -c "print ('%.6f' % (1.0*($min+$d)))")
  rm -fr salida
  for((i=0;i<=N;i++))
  do
    rescal=$(python3 -c "print ('%.6f' % ($i*1.0*($fin-$ini)/$N+$ini))")
    start $rescal
  done
  mv salida Vol.dat


##------------------------------------------------------------------------------------------## 

min=$(/home/dani/bin/get_min.py Vol.dat)
$FIREBALLHOME/TESTS/unitarios/./rescal_lvs.py uno.lvs $min > min.lvs
$FIREBALLHOME/TESTS/unitarios/./rescal_bas.py uno.bas $min > min.bas
$FIREBALLHOME/TESTS/unitarios/./rescal_kpts.py uno.kpts $min > min.kpts
$FIREBALLHOME/TESTS/unitarios/./rescal_kpts.py camino.kpts $min > camino_min.kpts

start 1.0


echo "&option
fdatalocation = '$fdata'
basisfile = 'min.bas'
lvsfile = 'min.lvs'
kptpreference = 'camino_min.kpts'
iqout = $iqout
iquench = -1
sigmatol = 0.000001
nstepf = 1
ifixcharge = 1
&end
&output
iwrteigen = 1
iwrtxyz = 1
&end" > fireball.in
#$FIREBALLHOME/fireball.x > salida.out 
/home/dani/FB/git/bin/fireball.21.09_intel-double.x > salida.out 

fermi=$(grep "Fermi Level" salida.out | tail -1 | cut -d'=' -f2)
$FIREBALLHOME/pyfb/plotbands.py -r ek.dat -fermi ${fermi} -print > ek0.dat
$FIREBALLHOME/pyfb/plotbands.py -r ek0.dat -c dia -save
cp Vol.dat Vol_${name}.dat
cp uno.lvs ${name}.lvs
cp uno.bas ${name}.bas
cp uno.kpts ${name}.kpts
cp CHARGES CHARGES_${name}
cp ek.dat ek_${name}.dat
cp ek0.dat ek0_${name}.dat


