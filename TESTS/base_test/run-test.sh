
basedir=$1
create=$2

h=$(pwd)
rm -fr aux
cp -r test aux
cd aux
rm -fr Fdata
ln -s ${basedir}/coutput Fdata

for iqout in 1 #3
do

for i in *  
do 
if [[ $i != 'run-info.sh' && $i != 'fireball.x' &&  $i != 'Fdata' && $i != 'run-all.sh' && $i != 'med.py' ]] 
then 
if test -d $i
then
cd $i
echo $i >> log
rm -fr Fdata
ln -s ../Fdata .
ln -s ../fireball.x .
./fast05S.sh $iqout ${create} $i

cd ..
fi
fi
done
done
cd ..
rm -fr ${h}/aux

