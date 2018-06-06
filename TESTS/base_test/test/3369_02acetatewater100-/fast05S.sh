### GPLv3 ###

iqout=${1:-1}
create=$2
name=$3

function start {
echo "&option
basisfile = uno.bas
lvsfile = uno.lvs
kptpreference = uno.kpts
idipole=1
idftd3=0
iqout = $iqout
icluster = 1
qstate = 1
sigmatol = 0.000001
nstepf = 400
&end
&output
iwrtxyz = 1
&end" > fireball.in 
./fireball.x > salida.out    
}
start
natoms=$(head -1 answer.xyz )

if [[ $iqout == 1 ]]; then  d="../../info-Lowdin/";fi
if [[ $iqout == 3 ]]; then  d="../../info-NAP/";fi


tail -$((natoms+2))  answer.xyz  >${d}${name}_${create}.xyz
cp CHARGES ${d}CHARGES_${name}_${create}

if ! test -d ${d}FIX; then mkdir ${d}FIX; fi
head -$((natoms+2))  answer.xyz  >${d}FIX/${name}_${create}.xyz
