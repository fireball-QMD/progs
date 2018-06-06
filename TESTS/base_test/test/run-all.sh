

#for i in info*
#do
#echo run-info.sh $(echo $i| cut -d'_' -f2-)
#../test/run-info.sh $(echo $i| cut -d'_' -f2-) 
#done


Elemento=${1:-C}
entrada=${2:-info}

mkdir bio bio_FIX
rm -fr bio/* bio_FIX/*


if [[ $entrada == info ]]
then
for i in info*
do
if [[  $(echo $i| cut -d'_' -f2) != FIX ]]
then
echo run-info.sh $(echo $i| cut -d'_' -f2-) $Elemento
../test/run-info.sh $(echo $i| cut -d'_' -f2-) $Elemento 
fi
done
else
#no se han generado los info... ?
for i in H2*
do
create=$(echo $i| cut -d'_' -f2-|cut -d'.' -f1)
echo run-info.sh $create $Elemento
../test/run-info.sh $create $Elemento
done
fi


#test='D E_IOHB E_s66 PA P E_ALLs66_MP2-cc-pVTZ-CP E_ALLs66_CCSDT-CBShaTZ-CP'

#for j in $test
#do
#rm -fr INFO_$j
#rm -fr INFO_FIX_$j
#rm -fr INFO_suma_$j
#rm -fr INFO_FIX_suma_$j
#
#done
#
#
#cd bio
#
#for j in $test
#do
#for i in ${j}*
#do
#echo $i $(../../test/med.py $i)
#done  >> ../INFO_$j
#if [[ $j == 'E_s66' || $j == 'E_ALLs66_MP2-cc-pVTZ-CP' || $j == 'E_ALLs66_CCSDT-CBShaTZ-CP' ]]
#then
#for i in ${j}*
#do
#echo $i $(../../test/suma.py $i)
#done  >> ../INFO_suma_$j
#fi
#done
#
#cd ..
#
#mkdir bio_FIX
#
#cd bio_FIX
#
#for j in $test
#do
#for i in ${j}*
#do
#echo $i $(../../test/med.py $i)
#done  >> ../INFO_FIX_$j
#if [[ $j == 'E_s66' || $j == 'E_ALLs66_MP2-cc-pVTZ-CP' || $j == 'E_ALLs66_CCSDT-CBShaTZ-CP' ]]
#then
#for i in ${j}*
#do
#echo $i $(../../test/suma.py $i)
#done  >> ../INFO_FIX_suma_$j
#fi
#done

#cd ..
