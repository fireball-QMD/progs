### GPLv3 ###
#autor : Daniel G. Trabada GPLv3

here=$(pwd)

#--------------------------------------#
#--------------------------------------#
#--------------------------------------#

FB_HOME=/home/dani/FB/git/

ln -s $FB_HOME/progs/fireball.x fireball.x
cp -r $FB_HOME/create/ .
cp -r $FB_HOME/begin/begin_rcatms/ .
cp pp/* begin_rcatms/
mkdir info-NAP info-Lowdin BASES

ele='1 6 7 8' # 12 15 16'
#exchange-correlation model
#LDA = 3   : GGA = 9
EX=9


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
mixns[$j]=0.00;  mixnp[$j]=0.00;  mixnd[$j]=0.00; 
mixs[$j]=0.00;  mixp[$j]=0.00;  mixd[$j]=0.00; 
opEx[$j]=0
excited[$j]=N
done

e=1
#X[${e}]=H  ; Z[${e}]=001; M[${e}]=1.0079; orb[${e}]=s; 
X[${e}]=H  ; Z[${e}]=001; M[${e}]=1.0079; orb[${e}]=s; 
excited[${e}]=Y ; opEx[${e}]=3 
mixns[${e}]=1.00; mixs[${e}]=0.80;
n0s[${e}]=1.00; ns[${e}]=1.00; rs[${e}]=5.42


e=6
X[$e]=C ; Z[$e]=006 ; M[$e]=12.011 ; orb[$e]=sp; excited[${e}]=N
n0s[${e}]=1.00 ; n0p[${e}]=3.00 ; n0d[${e}]=0.00
ns[${e}]=2.00  ;  np[${e}]=1.00  ; nd[${e}]=0.00
rs[${e}]=5.95  ; rp[${e}]=5.95  ; rd[${e}]=5.95

e=7
X[${e}]=N ; Z[${e}]=007; M[${e}]=14.0067 ; orb[${e}]=sp ; excited[${e}]=N
n0s[${e}]=2.00; n0p[${e}]=3.00; n0d[${e}]=0.00
ns[${e}]=2.00 ; np[${e}]=2.50 ; nd[${e}]=0.00
rs[${e}]=5.42 ; rp[${e}]=5.42 ; rd[${e}]=5.42

e=8
X[${e}]=O; Z[${e}]=008; M[${e}]=15.9994 ; orb[${e}]=sp ; excited[${e}]=N
n0s[${e}]=2.00 ; n0p[${e}]=4.00 ; n0d[${e}]=0.00
ns[${e}]=2.00  ; np[${e}]=4.00  ; nd[${e}]=0.00
rs[${e}]=5.32  ; rp[${e}]=5.32  ; rd[${e}]=5.32

e=12
X[${e}]=Mg; Z[${e}]=012 ; M[${e}]=24.305; orb[${e}]=s
n0s[${e}]=2.00 ; n0p[${e}]=0.00 ; n0d[${e}]=0.00
ns[${e}]=0.00  ; np[${e}]=0.00  ; nd[${e}]=0.00
rs[${e}]=4.50  ; rp[${e}]=4.50  ; rd[${e}]=4.50

e=15
X[${e}]=P ; Z[${e}]=015 ; M[${e}]=30.974; orb[${e}]=spd
n0s[${e}]=2.00 ; n0p[${e}]=3.00 ; n0d[${e}]=0.00
ns[${e}]=1.20  ; np[${e}]=1.20  ; nd[${e}]=0.00
rs[${e}]=4.50  ; rp[${e}]=5.00  ; rd[${e}]=7.00

e=16
X[${e}]=S ; Z[${e}]=016 ; M[${e}]=32.065; orb[${e}]=sp
n0s[${e}]=2.00; n0p[${e}]=4.00 ; n0d[${e}]=0.00
ns[${e}]=2.00  ; np[${e}]=3.50  ; nd[${e}]=0.00
rs[${e}]=3.80  ; rp[${e}]=4.30  ; rd[${e}]=5.00

#------------------------------------------------
#------------------------------------------------
#------------------------------------------------
#------------------------------------------------


mixs[1]=$(echo ${mixs[1]} + $(head -1 param/H_mixs) | bc -l)
rs[1]=$(echo ${rs[1]} + $(head -1 param/H_rs) | bc -l)

np[6]=$(echo ${np[6]} + $(head -1 param/C_np) | bc -l)
rp[6]=$(echo ${rp[6]} + $(head -1 param/C_rp) | bc -l)
rs[6]=${rp[6]}

np[7]=$(echo ${np[7]} + $(head -1 param/N_np) | bc -l)
rp[7]=$(echo ${rp[7]} + $(head -1 param/N_rp) | bc -l)
rs[7]=${rp[7]} 

np[8]=$(echo ${np[8]} + $(head -1 param/O_np) | bc -l)
rp[8]=$(echo ${rp[8]} + $(head -1 param/O_rp) | bc -l)
rs[8]=${rp[8]} 

rs[12]=$(echo ${rs[12]} + $(head -1 param/Mg_rs) | bc -l)

np[15]=$(echo ${np[15]} + $(head -1 param/P_np) | bc -l)
ns[15]=${np[15]}
rp[15]=$(echo ${rp[15]} + $(head -1 param/P_rp) | bc -l)
rs[15]=$(echo ${rp[15]} - 0.50 | bc -l)

np[16]=$(echo ${np[16]} + $(head -1 param/S_np) | bc -l)
rp[16]=$(echo ${rp[16]} + $(head -1 param/S_rp) | bc -l)
rs[16]=$(echo ${rp[16]} - 0.50 | bc -l)


for r in 1 
#for((r=0 ; r<1 ; r++))  #0.1 0.5
# l='-1.50 -1.40 -1.30 -1.20 -1.10 -1.00 -0.90 -0.80 -0.70 -0.60 -0.50 -0.40 -0.30 -0.20 -0.10 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 1.10 1.20 1.30 1.40 1.50'
do
cd $here
#ns[$v]=$(echo ${ns[$v]} + 0.10 | bc -l)
#np[$v]=$(echo ${np[$v]} + 0.10 | bc -l)
#rs[$v]=TOL
#rp[$v]=TOL
#rd[$v]=TOL

#tols[$v]=0.01
#tolp[$v]=0.01
#told[$v]=0.01

#tols[$v]=$tol 
#tolp[$v]=$tol
#told[$v]=$tol

#tols[$v]=POT
#vrs[$v]=0.00
#v0s[$v]=0.00
#tolp[$v]=POT
#vrp[$v]=0.00
#v0p[$v]=0.00
#told[$v]=POT
#vrd[$v]=$vr
#v0d[$v]=$v0

#--------------------------------------#
#--------------------------------------#
#--------------------------------------#
create=""
first=0
for j in $ele
do
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
done #atomos

#quitamos los puntos
create=$(echo $create | sed 's/\.//g')
basedir=$(pwd)/BASES/$create

if test -s ${basedir}.cinput
then
 echo ${basedir}.cinput existe no lo hacemos
else
 echo ${basedir}.cinput no existe lo hacemos
 cp -r create $basedir
 cd $here
 cd begin_rcatms
 for j in $ele 
 do 
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
  echo  ./begin.sh $op  >> log
  ./begin.sh $op  > /dev/null
  mv cinput/* $basedir/cinput/
 done #en atomos

fi # para que ejecute el test si o si 
 cd $basedir
 cp cinput/*input .
 echo $N > create.input
 for j in $ele
 do
  echo ${X[j]}.input >> create.input
 done
 ./create.com
 cd -
 cd $here
 ./run-test.sh ${basedir} $create 
 cd $here
 cd info-Lowdin
 ../test/run-info.sh $create
 cd ../info-NAP
 ../test/run-info.sh $create
 cd $here
 mv $basedir/cinput ${basedir}.cinput
#rm -fr $basedir 

#fi # si existe no ejecuta el test
done
cd $here
rm -fr aux test
rm -fr create begin* fireball.x
