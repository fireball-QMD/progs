#!/bin/bash

### GPLv3 ###
#autor : Daniel G. Trabada GPLv3
#
#El script necesita esas variables:
if [ -z "$FIREBALLHOME" ] || [ -z "$BEGINHOME" ] || [ -z "$CREATEHOME" ]; then
    echo "You should define the following variables in your .bashrc or export them:"
    echo "export FIREBALLHOME=/path/to/fireball"
    echo "export BEGINHOME=/path/to/begin"
    echo "export CREATEHOME=/path/to/create"
    exit 1  # Detener el script
fi

#El caso de que no tengas los ppfiles en el param tendras que tener definido tambien $PPFIREBALL#para ello:
#git clone git@github.com:fireball-QMD/fireball-qmd.github.io.git
#export PPFIREBALL=/path/fireball-qmd.github.io/tabla


here=$(pwd)
basis_file=${FIREBALLHOME}/TESTS/workflow_begin_create/basis_configuration
electronic_configuration=${FIREBALLHOME}/TESTS/workflow_begin_create/electronic_configuration

#preparamos create
if test -f create
then
  rm create
fi

cp -r $CREATEHOME .
basis_create=create/coutput/basis
mkdir $basis_create
mkdir ${basis_create}/param
if [ ! -e "create/coutput/basis" ]; then
  rm -fr create/cinput/*
fi

# Verificar si el archivo BASE existe
if [ ! -d BASE ]; then
    echo "la carpeta BASE no existe. Creándolo..."
    mkdir BASE
fi

#preparamos begin
if test -d begin_rcatms
then
  rm -fr begin_rcatms
fi
cp -r ${BEGINHOME}/begin_rcatms/ .

#ele='1 5 6 7 8 14 29'
ele='14'

#exchange-correlation model
#LDA = 3   : GGA = 9
EX=$(head param/EX)
cp param/EX ${basis_create}/param

if [ "$EX" -eq 3 ]; then
    label_ex="LDA"
else
    label_ex="GGA"
fi

N=0
for j in $ele
do
N=$((N+1))
n0s[$j]=0.00; n0p[$j]=0.00; n0d[$j]=0.00; 
ns[$j]=0.00;  np[$j]=0.00;  nd[$j]=0.00; 
rs[$j]=0.00;  rp[$j]=0.00;  rd[$j]=0.00; 
tols[$j]=FALSE; tolp[$j]=FALSE; told[$j]=FALSE; 
pots[$j]=FALSE; potp[$j]=FALSE; potd[$j]=FALSE; 
ion[$j]=FALSE #valores s p d
ions[$j]=FALSE; ionp[$j]=FALSE; iond[$j]=FALSE; 
v0s[$j]=0.0;  v0p[$j]=0.0;  v0d[$j]=0.0;
vrs[$j]=0.0;  vrp[$j]=0.0;  vrd[$j]=0.0;
Zpp[$j]=FALSE
ionPP[$j]=0; dionPP[$j]=0
mixns[$j]=0.00;  mixnp[$j]=0.00;  mixnd[$j]=0.00; 
mixs[$j]=0.00;  mixp[$j]=0.00;  mixd[$j]=0.00; 
opEx[$j]=0
excited[$j]=N
done
cd $here
#--------------------------------------#
#--------------------------------------#
#--------------------------------------#
create=""
first=0
for j in $ele
do
zzz=000$j
Z[${j}]=${zzz:$((${#zzz}-3))}
X[${j}]=$(cut -d' ' -f1,2 $basis_file | grep ${zzz:$((${#zzz}-3))} | cut -d' ' -f2) #Na Cl ....
cp param/${X[${j}]}_* ${basis_create}/param	

#Importante, en el caso de que haya un 00Z.pp en param lo copia de param, si no:
if [ ! -f param/${Z[${j}]}.pp ]; then
    echo "El archivo ${X[${j}]}.pp no existe lo copiamos de $PPFIREBALL"
    pp_path=$PPFIREBALL/${X[${j}]}/${X[${j}]}_00.more/${Z[${j}]}_${X[${j}]}_${label_ex}/${Z[${j}]}
    cp ${pp_path}.pp ${here}/begin_rcatms
else
    cp param/${Z[${j}]}.pp ${here}/begin_rcatms
fi


M[${j}]=$(grep \ ${X[${j}]}\  $basis_file | cut -d' ' -f3)

if [[ $(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | grep s -c) == 1 ]]
then
n0s[${j}]=$(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | cut -d's' -f2 | cut -d' ' -f2)
else
n0s[${j}]=0.00
fi

if [[ $(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | grep p -c) == 1 ]]
then
n0p[${j}]=$(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | cut -d'p' -f2 | cut -d' ' -f2)
else
n0p[${j}]=0.00
fi

if [[ $(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | grep d -c) == 1 ]]
then
n0d[${j}]=$(grep \ ${X[${j}]}\   $electronic_configuration | cut -d'+' -f2  | cut -d'|' -f1 | cut -d'd' -f2 | cut -d' ' -f2)
else
n0d[${j}]=0.00
fi

orb[${j}]=$(head -1 param/${X[j]}_orb)
rs[${j}]=$(head -1 param/${X[j]}_rs)
rp[${j}]=$(head -1 param/${X[j]}_rp)
rd[${j}]=$(head -1 param/${X[j]}_rd)
ns[${j}]=$(head -1 param/${X[j]}_ns)
np[${j}]=$(head -1 param/${X[j]}_np)
nd[${j}]=$(head -1 param/${X[j]}_nd) 
excited[${j}]=$(head -1 param/${X[j]}_exc)


#Importante, en el caso de que haya un 00Z++.pp en param lo copia de param, si no:
if [ ! -f param/${Z[${j}]}++.pp ]; then
    echo "El archivo ${X[${j}]}.pp no existe lo copiamos de $PPFIREBALL"
    pp_path=$PPFIREBALL/${X[${j}]}/${X[${j}]}_00.more/${Z[${j}]}_${X[${j}]}_${label_ex}/${Z[${j}]}
    if test -f ${pp_path}++.pp
    then
       cp ${pp_path}++.pp ${here}/begin_rcatms
    else #si no exsite copiamos el 00Z.pp a 00Z++.pp
       cp ${pp_path}.pp ${here}/begin_rcatms/${Z[${j}]}++.pp
    fi
else
    cp param/${Z[${j}]}++.pp ${here}/begin_rcatms/
fi


opEx[$j]=$(head -1 param/${X[j]}_opEx)
mixns[${j}]=$(head -1 param/${X[j]}_mixns); 
mixs[${j}]=$(head -1 param/${X[j]}_mixs)

ion[${j}]=$(head -1 param/${X[j]}_ion)
ions[${j}]=$(head -1 param/${X[j]}_ions)
ionp[${j}]=$(head -1 param/${X[j]}_ionp)
iond[${j}]=$(head -1 param/${X[j]}_iond)

ionPP[${j}]=$(head -1 param/${X[j]}_ionPP)
dionPP[${j}]=$(head -1 param/${X[j]}_dionPP)

if test -f param/${X[j]}_Zpp
then
Zpp[${j}]=$(head -1 param/${X[j]}_Zpp)
fi

if [[ ${mixs[$j]::1} == '.' ]] ; then mixs[$j]=0${mixs[$j]}; fi
if [[ ${ns[$j]::1} == '.' ]] ; then ns[$j]=0${ns[$j]}; fi
if [[ ${np[$j]::1} == '.' ]] ; then np[$j]=0${np[$j]}; fi
if [[ ${nd[$j]::1} == '.' ]] ; then nd[$j]=0${nd[$j]}; fi
if [[ ${rs[$j]::1} == '.' ]] ; then rs[$j]=0${rs[$j]}; fi
if [[ ${rp[$j]::1} == '.' ]] ; then rp[$j]=0${rp[$j]}; fi
if [[ ${rd[$j]::1} == '.' ]] ; then rd[$j]=0${rd[$j]}; fi
if [[ ${mixs[j]} ==  0 ]] ; then mixs[$j]=0.00; fi
if [[ ${ns[j]} ==  0 ]] ; then ns[$j]=0.00; fi
if [[ ${np[j]} ==  0 ]] ; then np[$j]=0.00; fi
if [[ ${nd[j]} ==  0 ]] ; then nd[$j]=0.00; fi

 if [[ $first == 0 ]]
 then
  create=${X[j]}
  first=1
 else
 create=${create}_${X[j]}
 fi


 if [[ ${Zpp[j]} != FALSE ]]
 then
  create=${create}zpp${Zpp[j]}
 fi

  for((o=1;o<=${#orb[j]};o++))
  do
   orb=${orb[j]:$((o-1)):1}
   if [[ $orb == s ]]
   then
    if [[ ${tols[j]} != FALSE ]]
    then
     create=${create}sT${tols[j]}
    else
     if [[ ${excited[j]} == Y ]]
     then
       if [[ ${opEx[j]} == 1 ]] ; then create=${create}s1 ; fi
       if [[ ${opEx[j]} == 2 ]] ; then create=${create}s ;  fi
       if [[ ${opEx[j]} == 4 ]] ; then create=${create}opEX4 ; fi
       if [[ ${opEx[j]} == 3 || ${opEx[j]} == 4 ]] ; then create=${create}mx${mixns[j]}${mixs[j]} ;  fi
     fi
     create=${create}s${ns[j]}${rs[j]}
    fi
    if [[ ${pots[j]} != FALSE ]]
    then
     create=${create}p${vrs[j]}${v0s[j]}
    fi
    if [[ ${ion[j]} == s ]]
    then
     create=${create}Is${ions[j]}${ionp[j]}${iond[j]}
    fi
   fi 
   
   if [[ $orb == p ]]
   then
    if [[ ${tolp[j]} != FALSE ]]
    then
     create=${create}pT${tolp[j]}
    else
     if [[ ${excited[j]} == Y ]]
     then
       if [[ ${opEx[j]} == 1 ]] ; then create=${create}p1 ; fi
       if [[ ${opEx[j]} == 2 ]] ; then create=${create}p ;  fi
       if [[ ${opEx[j]} == 4 ]] ; then create=${create}opEX4 ; fi
       if [[ ${opEx[j]} == 3 ||  ${opEx[j]} == 4 ]] ; then create=${create}mx${mixns[j]}${mixs[j]} ;  fi
     fi
     create=${create}p${np[j]}${rp[j]}
    fi
    if [[ ${potp[j]} != FALSE ]]
    then
     create=${create}p${vrp[j]}${v0p[j]}
    fi
    if [[ ${ion[j]} == p ]]
    then
     create=${create}Ip${ions[j]}${ionp[j]}${iond[j]}
    fi
   fi 
 
  if [[ $orb == d ]]
  then
   if [[ ${told[j]} != FALSE ]]
   then
    create=${create}dT${told[j]}
   else
    if [[ ${excited[j]} == Y ]]
    then
       if [[ ${opEx[j]} == 1 ]] ; then create=${create}d1 ; fi
       if [[ ${opEx[j]} == 2 ]] ; then create=${create}d ;  fi
       if [[ ${opEx[j]} == 4 ]] ; then create=${create}opEX4 ; fi
       if [[ ${opEx[j]} == 3 || ${opEx[j]} == 4 ]] ; then create=${create}mx${mixns[j]}${mixs[j]} ;  fi
    fi
    create=${create}d${nd[j]}${rd[j]}
   fi
   if [[ ${potd[j]} != FALSE ]]
   then
    create=${create}d${vrd[j]}${v0d[j]}
   fi
   if [[ ${ion[j]} == d ]]
   then
    create=${create}Id${ions[j]}${ionp[j]}${iond[j]}
   fi
  fi 
 done #orbitales
if [[ ${ionPP[j]} == 1 || ${ionPP[j]} == 2 || ${ionPP[j]} == 3 ]]
then
create=${create}ionPP${ionPP[j]}${dionPP[j]}
fi
#quitamos los puntos
create=$(echo $create | sed 's/\.//g')


done #atomos

#quitamos los puntos
create=$(echo $create | sed 's/\.//g')
basedir=${here}/create

cd $here
 cp ${FIREBALLHOME}/TESTS/workflow_begin_create/switch.input create/
 cd $here
 cd begin_rcatms
 for j in $ele 
 do 
  #Importante nos aseguramos que los *pp son los que hay en param
  cp ${basis_create}/param/basis/${Z[${j}]}*.pp .

  op="-ex $EX -z ${Z[j]} -ele ${X[j]} -mass ${M[j]} "
  if [[ ${Zpp[j]} != FALSE ]]
  then
    op=$op"-pp ${Zpp[j]} "
  fi

  for((o=1;o<=${#orb[j]};o++))
  do
   orb=${orb[j]:$((o-1)):1}
   op=$op"-orb $orb "

   if [[ $orb == s ]]
   then
    if [[ ${excited[j]} == Y ]] 
    then
     op=$op" -excited "${opEx[j]}" " 
     if [[ ${opEx[j]} = 3 || ${opEx[j]} == 4 ]]
     then
       op=$op" ${mixns[j]} ${mixs[j]} "
     fi
    fi
    op=$op"-n ${ns[j]} -n0 ${n0s[j]} "
    if [[ ${tols[j]} != FALSE ]]
    then
    op=$op"-tol ${tols[j]} "
    else
    op=$op"-r ${rs[j]} "
    fi
    if [[ ${pots[j]} != FALSE ]]
    then
     op=$op"-pot -vr ${vrs[j]} -v0 ${v0s[j]} "
    fi
    if [[ ${ion[j]} == s ]]
    then
     op=$op"-ion ${ions[j]} ${ionp[j]} ${iond[j]} "
    fi
   fi 
 
   if [[ $orb == p ]]
   then
    if [[ ${excited[j]} == Y ]] 
    then
     op=$op"-excited "${opEx[j]}" " 
     if [[ ${opEx[j]} = 3 || ${opEx[j]} == 4 ]]
     then
       op=$op" ${mixns[j]} ${mixnp[j]} ${mixs[j]} ${mixp[j]} "
     fi
    fi
    op=$op"-n ${np[j]} -n0 ${n0p[j]} "
    if [[ ${tolp[j]} != FALSE ]]
    then
     op=$op"-tol ${tolp[j]} "
    else
     op=$op"-r ${rp[j]} "
    fi
    if [[ ${potp[j]} != FALSE ]]
    then
     op=$op"-pot -vr ${vrp[j]} -v0 ${v0p[j]}  "
    fi
    if [[ ${ion[j]} == p ]]
    then
     op=$op"-ion ${ions[j]} ${ionp[j]} ${iond[j]} "
    fi
   fi 
 
   if [[ $orb == d ]]
   then
    if [[ ${excited[j]} == Y ]] 
    then
     op=$op"-excited "${opEx[j]}" " 
     if [[ ${opEx[j]} = 3 || ${opEx[j]} == 4 ]]
     then
       op=$op" ${mixns[j]} ${mixnp[j]} ${mixnd[j]} ${mixs[j]} ${mixp[j]} ${mixd[j]} "
     fi
    fi
    op=$op"-n ${nd[j]} -n0 ${n0d[j]} "
     if [[ ${told[j]} != FALSE ]]
    then
     op=$op"-tol ${told[j]} "
    else
     op=$op"-r ${rd[j]} "
    fi
    if [[ ${potd[j]} != FALSE ]]
    then
     op=$op"-pot -vr ${vrd[j]} -v0 ${v0d[j]} "
    fi
    if [[ ${ion[j]} == d ]]
    then
     op=$op"-ion ${ions[j]} ${ionp[j]} ${iond[j]} "
    fi
   fi
  done #orbitales 
  if [[ ${ionPP[j]} == 1 || ${ionPP[j]} == 2 || ${ionPP[j]} == 3 ]]
  then
    op=$op" -ionPP${ionPP[j]} ${dionPP[j]}"
  fi

  echo  ./begin.sh $op  >> ${basedir}/coutput/basis/begin.log
  ./begin.sh $op  > /dev/null
  cp -r cinput/* $basedir/cinput/
 done #en atomos

 cd $basedir
 cp cinput/*input .
 echo $N > create.input
 for j in $ele
 do
  echo ${X[j]}.input >> create.input
 done
 echo 
 echo start create $create $(date) >> ${here}/create/base.log
 cd $basedir
 ./create.com
 echo $create | sed 's/_/\n/g' >  $basedir/cinput/README
 cp -r $basedir/cinput/* ${basedir}/coutput/basis/

 cd $here
 
# rm -fr create begin_rcatms 


