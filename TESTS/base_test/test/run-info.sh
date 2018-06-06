#for i in ../BASES/Base_H_ns_0.10_*; do create=$(echo $i|sed 's/\.\.\/BASES\///'| sed 's/\.cinput//'); ../test/./run-info.sh $create;done
create=$1
Elemento=${2:-C}

mkdir bio
rm -fr info_${create}
rm -fr bio/E_s66_${create}
rm -fr bio/E_IOHB_${create}
rm -fr bio/D_${create}
rm -fr bio/P_${create}


rm -fr info_FIX_${create}
rm -fr bio/E_s66_FIX_${create}
rm -fr bio/E_IOHB_FIX_${create}
rm -fr bio/D_FIX_${create}
rm -fr bio/P_FIX_${create}




#function xeo { 
#/home/dani/bin/jdk1.8.0_05/bin/java -jar /home/dani/bin/xeo/xeo.jar $@ 2> /dev/null; 
#}


DEC=12
function Calc {
aux=$(echo ${@} | sed 's/X/\*/g' |  sed s/\ //)
python -c  "import random; import math ; print '%.${DEC}f' % (1.0*${aux})"
}

function D { #anser.xyz 1 3
f=$(python -c "
i=-1
for line in file(\"$1\"):
   line = line.split()
   if i == $2 :
      x1 = float(line[1])
      y1 = float(line[2])
      z1 = float(line[3])
   if i == $3 :
      x2 = float(line[1])
      y2 = float(line[2])
      z2 = float(line[3])
   i=i+1 
d=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
print d")

echo $f
}

function nETOT {
E=$(grep ETOT ${1:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
n=$(head -1  ${1:-0}_${create}.xyz )
echo $(Calc $E / $n |cut -c-$DEC)
} 
function uETOT {
E=$(grep ETOT ${1:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
n=$(head -1  ${1:-0}_${create}.xyz )
ele=$(grep -c $Elemento  ${1:-0}_${create}.xyz )
echo $(Calc $E X $ele / $n |cut -c-$DEC)
}
#function D { #archivo.xyz 1 2  
#echo $(rm -fr o.dat; xeo -xyz ${1} o.dat 1 d[$2][$3] && cat o.dat | cut -d' ' -f2-| cut -c-6  ; rm -fr o.dat)
#}

function dC { 
D ${1}_${create}.xyz $2 $3
#aux=$(rm -fr o.dat; xeo -xyz ${1}_${create}.xyz o.dat 1 d[$2][$3] && cat o.dat | cut -d' ' -f2-| cut -c-6  ; rm -fr o.dat) 
#echo $aux
}


function A { 
#echo $(rm -fr o.dat; xeo -xyz ${1} o.dat 1 a[$2][$3][$4] && cat o.dat | cut -d' ' -f2-| cut -c-6  ; rm -fr o.dat) 
echo 0
}

function completa15 { 
t=${1}............................................... ; echo ${t::15} 
}


function EV {
#1eV = 23.061 Kcal/mol
Calc ${1}/23.061
}

function ETOT {
echo $(grep ETOT ${1:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
}

function diffEH {
E1=$(grep ETOT ${1:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E2=$(grep ETOT ${2:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EH=$(grep ETOT H_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1 | cut -d'e' -f1)
Calc  $E1 - $E2
}

function E1_23 {
EA=$(grep ETOT ${2}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1 | cut -d'e' -f1)
EB=$(grep ETOT ${3}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1|  cut -d'e' -f1)
EAB=$(Calc  $EA + $EB |cut -c-10)
ETOT=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
diff=$(Calc $ETOT - $EAB | cut -c-10)
echo $(completa15 $1)$'\t'$'\t'$'\t'$EAB$'\t'$ETOT$'\t'$diff 
}


function diffE1_23 {
EA=$(grep ETOT ${2}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EB=$(grep ETOT ${3}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EAB=$(Calc  $EA + $EB | cut -c-10)
ETOT=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
diff=$(Calc $ETOT - $EAB |cut -c-10)
echo $diff 
}

function Eab {
EA=$(grep ETOT ${1}A_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EB=$(grep ETOT ${1}B_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EAB=$(Calc  $EA + $EB  |cut -c-10)
ETOT=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
diff=$(Calc $ETOT - $EAB  |cut -c-10)
echo $(completa15 $1)$'\t'$'\t'$'\t'$EAB$'\t'$ETOT$'\t'$diff 
}
function diffEab {
EA=$(grep ETOT ${1}A_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EB=$(grep ETOT ${1}B_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
EAB=$(Calc $EA + $EB | cut -c-10)
ETOT=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
diff=$( Calc $ETOT - $EAB | cut -c-10)
echo $diff 
}


function Barrera {
E1=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E2=$(grep ETOT ${2}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E3=$(grep ETOT ${3}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
echo $( Calc $E1 - $E2 - $E3 | cut -c-$DEC)
}

function Barrera2 {
E1=$(grep ETOT ${1}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E2=$(grep ETOT ${2}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E3=$(grep ETOT ${3}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E4=$(grep ETOT ${4}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
echo $(Calc $E1 + $E2- $E3 - $E4 | cut -c-$DEC)
}

function diffE {
E1=$(grep ETOT ${1:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
E2=$(grep ETOT ${2:-0}_${create}.xyz | cut -d'=' -f2 | cut -d'T' -f1| cut -d'e' -f1)
echo $(Calc $E1 - $E2 |cut -c-$DEC)$
}

function INTR {
echo INTR
echo "-------INTR. -------
ENLACE			l(exp)	l(fire)	Ang(exp)	Ang(fire)
H-H     H-H		0.74	$(D H2_${create}.xyz 1 2)
C-C	HC-CH		1.20	$(D CH-CH_${create}.xyz 1 2)	180º		$(A CH-CH_${create}.xyz 2 1 3) 
	CH2-CH2         1.33	$(D CH2-CH2_${create}.xyz 1 2)	116º 122º 	$(A CH3-CH3_${create}.xyz 5 2 6 )	$(A CH3-CH3_${create}.xyz 2 1 3)
	CH3-CH3		1.54	$(D CH3-CH3_${create}.xyz 1 2)
	Benceno		1.39    $(D benceno_${create}.xyz 1 3)
C-H	CH3CH2-H	1.07	$(D CH3-CH3_${create}.xyz 2 7)
C-O	CH2-OH		1.43	$(D CH2-OH_${create}.xyz 1 2)
	CH2-O		1.22	$(D COH2_${create}.xyz 1 2) 
C-N	CH3-NH2		1.47	$(D CH3-NH2_${create}.xyz 3 2) 
	CH3CH-NOH	1.28	$(D CH3CH-NOH_${create}.xyz 3 2) 
	CH3C-N		1.16	$(D CH3C-N_${create}.xyz 1 2 )
O-H	CH3O-H		0.96	$(D CH3O-H_${create}.xyz 1 2 )
	H2O		0.96	$(D water_${create}.xyz 1 2 )	$(ETOT water )
O-O	O2		1.21  	$(D O2_${create}.xyz 1 2) 
N-N	N2		1.112	$(D N2_${create}.xyz 1 2)
" >>  info_${create}

echo D

echo "HC-CH 1.20 $(D CH-CH_${create}.xyz 1 2)
CH2-CH2 1.33 $(D CH2-CH2_${create}.xyz 1 2)
CH3-CH3 1.54 $(D CH3-CH3_${create}.xyz 1 2)
Benceno 1.39 $(D benceno_${create}.xyz 1 3)
CH3CH2-H 1.07 $(D CH3-CH3_${create}.xyz 2 7)
CH2-OH 1.43 $(D CH2-OH_${create}.xyz 1 2)
CH2-O 1.22 $(D COH2_${create}.xyz 1 2) 
CH3-NH2 1.47 $(D CH3-NH2_${create}.xyz 3 2) 
CH3CH-NOH 1.28 $(D CH3CH-NOH_${create}.xyz 3 2) 
CH3C-N 1.16 $(D CH3C-N_${create}.xyz 1 2 )
CH3O-H 0.96 $(D CH3O-H_${create}.xyz 1 2 )
H2O 0.96 $(D water_${create}.xyz 1 2 ) 
acetate 1.36 $(D acetate_${create}.xyz 10 8)
acetic 1.43 $(D acetic_${create}.xyz 9 8) 
acetic 1.23 $(D acetic_${create}.xyz 10 8 )
acetone 1.23 $(D acetone_${create}.xyz 2 1)  
ethanol 1.43 $(D ethanol_${create}.xyz 9 8) 
diethyl_ether 1.43 $(D diethyl_ether_${create}.xyz 1 15) 
diethyl_peroxi. 1.49 $(D diethyl_peroxide_${create}.xyz 4 2) 
furan 1.43 $(D furan_${create}.xyz 1 2) 
peroxide 1.49 $(D peroxide_${create}.xyz 1 2) 
amoniaco 1.017 $(D amoniaco_${create}.xyz 1 2)
guanidina 1.47 $(D guanidina_${create}.xyz 7 9)
guanidina 1.35 $(D guanidina_${create}.xyz 6 7)
guanidinio 1.35 $(D guanidinio+_${create}.xyz 1 2)" >>  ../bio/D_${create}
}


function ANG {
echo ANG
echo "
-------ANGULOS-------

 H        H
  \      / 
   C----C  116º -> $(A CH3-CH3_${create}.xyz 5 2 6 )   
  /      \ 
 H        H 
    122º -> $(A CH3-CH3_${create}.xyz 2 1 3)

H--C----C--H  180º ->  $(A CH-CH_${create}.xyz 2 1 3)

 H   H
  \ /
   C
  / \ 
CH3 CH3
 112º -> $(A C3H8_${create}.xyz 1 2 3)

   O 
  / \ 
 H   H 
  105º -> $(A water_${create}.xyz 2 1 3)

   O  
  / \ 
CH3  CH3 
   111º -> $(A CH3-O-CH3_${create}.xyz 1 2 3)  

CH3
  \
   C--O  
  /      
CH3   120º -> $(A C2H6-O_${create}.xyz 1 2 3)
" >>  info_${create}
}

function TESTBIO {
echo TESTBIO
echo "-------TEST-BIO--------
HCO			exp.			fireball
	acetate		d(O,C) = 1.36		d(O,C)=$(D acetate_${create}.xyz 10 8) $(D acetate_${create}.xyz 9 8 )		a(C,O,C)=$(A acetate_${create}.xyz 10 8 9) 
	acetic 		d(O,C) = 1.43, 1.23	d(O,C)=$(D acetic_${create}.xyz 9 8) $(D acetic_${create}.xyz 10 8 )		a(C,O,C)=$(A acetic_${create}.xyz 10 8 9)	
	acetone 	d(O,C) = 1.23		d(C,O)=$(D acetone_${create}.xyz 2 1)			a(C,C,C)=$(A acetone_${create}.xyz 7 2 3);
	ethanol 	d(O,C) = 1.43 		d(C,O)=$(D ethanol_${create}.xyz 9 8) d(O,H)=$(D ethanol_${create}.xyz 9 5) "~H2O"
	diethyl_ester				d(C,O)=$(D diethyl_ester_${create}.xyz 4 1) d(C,O)=$(D diethyl_ester_${create}.xyz 1 1); d(C,O)=$(D diethyl_ester_${create}.xyz 2 5);
	diethyl_ether 	d(O,C) = ~1.43		d(C,O)=$(D diethyl_ether_${create}.xyz 1 15) d(C,O)=$(D diethyl_ether_${create}.xyz 14 15);  
	diethyl_peroxi. d(O,O) = ~1.49 		d(O,O)=$(D diethyl_peroxide_${create}.xyz 4 2) 
	furan		d(O,C) = 1.43		d(O,C)=$(D furan_${create}.xyz 1 2) d(C,O)=$(D furan_${create}.xyz 1 5) 
	peroxide 	d(O,O) = 1.49		d(O,O)=$(D peroxide_${create}.xyz 1 2) d(O,H)=$(D peroxide_${create}.xyz 1 3)
	water		H-O 0.98 Ang		$(D water_${create}.xyz 1 2)
			a(H,N,H) 108º		$(A water_${create}.xyz 2 1 3)
NC
        RET  					$(D RET_${create}.xyz 38 37) $(D RET_${create}.xyz 37 35) $(D RET_${create}.xyz 35 33) $(D RET_${create}.xyz 33 28) $(D RET_${create}.xyz 28 26) $(D RET_${create}.xyz 26 24) $(D RET_${create}.xyz 24 22) $(D RET_${create}.xyz 22 17) $(D RET_${create}.xyz 17 15) $(D RET_${create}.xyz 15 13) $(D RET_${create}.xyz 13 11) $(D RET_${create}.xyz 11 8)                          
        benceno                                 d(C,C)=$(D benceno_${create}.xyz 1 3)
HCN
        amoniaco        H-N 1.017 Ang           $(D amoniaco_${create}.xyz 1 2)	$(ETOT amoniaco ) 
                        a(H,N,H) 107.8º         $(A amoniaco_${create}.xyz 2 1 3)
        guanidina       360º			$(echo $(A guanidina_${create}.xyz 6 7 8)+$(A guanidina_${create}.xyz 8 7 9) + $(A guanidina_${create}.xyz 9 7 6 )|bc -l) 
                        d(N,C)=1.47, 1.35       d(N,C)=$(D guanidina_${create}.xyz 7 9), $(D guanidina_${create}.xyz 8 7); d(N,C)=$(D guanidina_${create}.xyz 6 7); a(N,C,N)=$(D guanidina_${create}.xyz 6 7 8)
        guanidinio	d(N,C) =  1.35; 120º	d(N,C)=$(D guanidinio+_${create}.xyz 1 2); d(N,C)=$(D guanidinio+_${create}.xyz 3 2); d(N,C)=$(D guanidinio+_${create}.xyz 4 2); a(N,C,N)=$(A guanidinio+_${create}.xyz 1 2 3)
        imidazol				d(C,N)=$(D imidazol_${create}.xyz 2 4); d(N,C)=$(D imidazol_${create}.xyz 4 5); d(C,N)=$(D imidazol_${create}.xyz 5 6);  d(N,C)=$(D imidazol_${create}.xyz 6 3);  d(C,C)=$(D imidazol_${create}.xyz 3 2);  
        imidazol+                             d(C,N)=$(D imidazol+_${create}.xyz 1 3); d(N,C)=$(D imidazol+_${create}.xyz 3 4); d(C,N)=$(D imidazol+_${create}.xyz 4 5);  d(N,C)=$(D imidazol+_${create}.xyz 5 2);  d(C,C)=$(D imidazol+_${create}.xyz 2 1);

" >>  info_${create}
}

function TESTCRC {
echo TESTCRC
echo "
#----- TEST - CRC Handbook of Chemistry and Physics Editor-in-Chief David R. Lide 
#Characteristic lengths of single bonds. (Å)																
	As	Br	C	Cl	F	Ge	H	I	N	O	P	S	Sb	Se	Si	Sn
As	2.10															
Br	2.32	2.28														
C	1.96	1.94	1.53													
Cl	2.17	2.14	1.79	1.99												
F	1.71	1.76	1.39	1.63	1.S66-41											
Ge		2.30	1.95	2.15	1.73	2.40										
H	1.51	1.S66-41	1.09	1.28	0.92	1.53	0.74									
I		2.47	2.13	2.32	1.91	2.51	1.61	2.67								
N			1.46	1.90	1.37		1.02		1.45 							
O			1.42	1.70	1.42		0.96		1.43 	1.48						
P		2.22	1.85	2.04	1.57		1.42		1.65 		2.25					
S		2.24	1.82	2.05	1.56		1.34					2.00				
Sb				2.33			1.70									
Se			1.95		1.71		1.47							2.33		
Si		2.21	1.87	2.05	1.58		1.48	2.44		1.63		2.14			2.33	
Sn			2.14	2.28			1.71	2.67								
Te					1.82		1.66									


---FIREBALL---        
	As	Br	C	Cl	F	Ge	H	I	N	O	P	S	Sb	Se	Si	Sn
As
Br
C			$(D CH3-CH3_${create}.xyz 1 2)
Cl
F
Ge
H			$(D CH3-CH3_${create}.xyz 2 7)				$(D H2_${create}.xyz 1 2)
I
N			$(D CH3-NH2_${create}.xyz 3 2)				$(D amoniaco_${create}.xyz 1 2)			$(D N2H4_${create}.xyz 1 2)		
O			$(D CH3-O-CH3_${create}.xyz 1 2)				$(D water_${create}.xyz 1 2)		$(D HNO2_${create}.xyz 1 2)	$(D H2O2_${create}.xyz 1 2)			
P			$(D C3H9-P_${create}.xyz 2 3)				$(D PH3_${create}.xyz 2 3) 		
S			$(D CH3SH_${create}.xyz 2 4)				$(D H2S_${create}.xyz 1 2)	$(D CH2N2S_${create}.xyz 2 4)			$(D H2S2_${create}.xyz 4 3)	
Sb
Se
Si
Sn
Te
																
Lengths of multiple bonds (non-ring molecules).												
	Handbook	FIREBALL				
C=C	1.34		$(D CH2-CH2_${create}.xyz 1 2)																
C÷C	1.20		$(D CH-CH_${create}.xyz 1 2)										
C=N	1.21		$(D CH3CH-NOH_${create}.xyz 3 2)													
C÷N	1.16		$(D CH3C-N_${create}.xyz 1 2 )		 											
C=O	1.21		$(D C2H6-O_${create}.xyz 1 2) 													
C=S	1.61		$(D CH2S_${create}.xyz 3 2) 																	
N=N	1.24		$(D CH3-N2-CH3_${create}.xyz 8 7)															
N÷N	1.13		$(D CH3-N2-CH3_${create}.xyz 7 6)															
N=O	1.18		$(D CH3NO2_${create}.xyz 5 4)															
O=O	1.21		$(D H2O2_${create}.xyz 1 2)													
											
Effect of environment on carbon-carbon single bonds (other single bonds not shown).																
Conf.		C÷C length	FIREBALL	Examples of molecules										
C÷C		1.526		$(D CH3-CH3_${create}.xyz 1 2)		H3C÷CH3										
C÷C=		1.501				H3C÷CH=CH2										
C÷C=		1.459				H3C÷C÷CH										
÷C÷C=		1.467				H2C=CH÷CH=CH2										
÷C÷C÷		1.445				HC÷C÷CH=CH2										
C÷C		1.378				HC÷C÷C÷CH										
																
Some metal-carbon bond lengths in gas-phase molecules. 																
Al÷C	1.96		Bi÷C	2.26		Pb÷C	2.24									
B÷C	1.58		Cd÷C	2.11		Sn÷C	2.14									
Be÷C	1.70		Hg÷C	2.08		Zn÷C	1.93
" >>  info_${create}
}


function S66 {
echo S66
echo '---------------- TEST S66 -------------------
ref. http://www.begdb.com/index.php?action=oneDataset&id=41&state=show&order=ASC&by=name_m&method= 
test              Ea+Eb      Etot	Etot-Ea-Eb	ref.E	Dc	 Dl	 ref.Dc	 ref.Dl'   >>  info_${create} 
a=S66-4113_01WaterWater
e=-0.2149
d1=$(dC $a 1 4)
d1r=2.924
d2=$(dC $a 3 4)
d2r=1.963
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)  >> ../bio/P_${create}
echo $a' '${d2r}' '${d2}' '$(uETOT $a)  >> ../bio/P_${create} 

a=S66-4115_03WaterMeNH2 
e=-0.3072
d1=$(dC $a 1 4)
d1r=2.929
d2=$(dC $a 3 4)
d2r=1.961
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=S66-4116_04WaterPeptide 
 e=-0.3500
d1=$(dC $a 1 9)
d1r=2.814
d2=$(dC $a 3 9)
d2r=1.861
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4117_05MeOHMeOH 
 e=-0.2532
d1=$(dC $a 1 7)
d1r=2.863
d2=$(dC $a 2 7)
d2r=1.905
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4118_06MeOHMeNH2 
 e=-0.3351
d1=$(dC $a 1 7)
d1r=2.909
d2=$(dC $a 2 7)
d2r=1.938
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4119_07MeOHPeptide 
 e=-0.3548
d1=$(dC $a 1 12)
d1r=2.814
d2=$(dC $a 2 12)
d2r=1.855
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4120_08MeOHWater 
 e=-0.2182
d1=$(dC $a 1 7)
d1r=2.917
d2=$(dC $a 7 2)
d2r=1.952
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4121_09MeNH2MeOH 
 e=-0.1327
d1=$(dC $a 1 8 )
d1r=3.178
d2=$(dC $a 2 8)
d2r=2.199
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4122_10MeNH2MeNH2 
 e=-0.1860
d1=$(dC $a 1 8 )
d1r=3.179
d2=$(dC $a 2 8 )
d2r=2.244
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4123_11MeNH2Peptide 
 e=-0.2398
d1=$(dC $a 1 13)
d1r=3.070
d2=$(dC $a 2 13)
d2r=2.203
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4124_12MeNH2Water 
 e=-0.3261
d1=$(dC $a 8 1)
d1r=2.888
d2=$(dC $a  9 1)
d2r=1.941
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4125_13PeptideMeOH 
 e=-0.2742
d1=$(dC $a 7 13 )
d1r=2.974
d2=$(dC $a 8  13)
d2r=1.989
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4126_14PeptideMeNH2 
 e=-0.3331
d1=$(dC $a 7 13 )
d1r=3.052
d2=$(dC $a 8 13 )
d2r=2.048
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4127_15PeptidePeptide 
 e=-0.3761
d1=$(dC $a 7 18 )
d1r=2.947
d2=$(dC $a 8 18)
d2r=1.947
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4128_16PeptideWater 
 e=-0.2232
d1=$(dC $a 7 13  )
d1r=3.058
d2=$(dC $a 8 13 )
d2r=2.052
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4129_17UracilUracilBP 
 e=-0.7445
d1=$(dC $a 1 24)
d1r=2.850
d2=$(dC $a 2 24)
d2r=1.824
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

#echo $(Eab $a)$'\t'-0.725$'\t'$(dC $a 1 24)$'\t'$(dC $a 2 24)$'\t'$(dC $a 4 21)$'\t'$(dC $a 4 22)$'\t'2.850   1.824   2.792   1.767  >> info_${create}
a=S66-4130_18WaterPyridine 
 e=-0.3065
d1=$(dC $a 1 4)
d1r=2.923
d2=$(dC $a 3 4)
d2r=1.953
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4131_19MeOHPyridine 
 e=-0.3328
d1=$(dC $a 1 7)
d1r=2.900
d2=$(dC $a 2 7 )
d2r=1.931
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

a=S66-4132_20AcOHAcOH 
 e=-0.8239
d1=$(dC $a 3 10 )
d1r=2.675
d2=$(dC $a 4 10  )
d2r=1.681
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

#echo $(Eab $a)$'\t'-0.808$'\t'$(dC $a 3 10 )$'\t'$(dC $a 4 10  )$'\t'$(dC $a 2 11 )$'\t'$(dC $a 2 12  )$'\t'2.675   1.681   2.675   1.681  >> info_${create}
a=S66-4133_21AcNH2AcNH2 
 e=-0.6991
d1=$(dC $a 3 11 )
d1r=2.864
d2=$(dC $a 2 12 )
d2r=1.844
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

#echo $(Eab $a)$'\t'-0.687$'\t'$(dC $a 3 11 )$'\t'$(dC $a 4 11  )$'\t'$(dC $a 2 12 )$'\t'$(dC $a 2 13  )$'\t'2.864   1.844   2.864   1.844  >> info_${create}
a=S66-4134_22AcOHUracil 
 e=-0.8413
d1=$(dC $a 2 17 )
d1r=2.788
d2=$(dC $a 2 18  )
d2r=1.767
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

#echo $(Eab $a)$'\t'-0.825$'\t'$(dC $a 2 17 )$'\t'$(dC $a 2 18  )$'\t'$(dC $a 3 20 )$'\t'$(dC $a 4 20  )$'\t'2.788   1.767   2.686   1.695  >> info_${create}
a=S66-4135_23AcNH2Uracil 
 e=-0.8283
d1=$(dC $a 2 18 )
d1r=2.763
d2=$(dC $a 2 19  )
d2r=1.730
echo $(Eab $a)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffEab $a)' '$(uETOT $a) >> ../bio/E_s66_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}

#echo $(Eab $a)$'\t'-0.813$'\t'$(dC $a 2 18 )$'\t'$(dC $a 2 19  )$'\t'$(dC $a 21 3 )$'\t'$(dC $a 21 4  )$'\t'2.763   1.730   2.896   1.883   >> info_${create}
#a=S66-4136_24BenzeneBenzenepipi 
# echo $(Eab $a) >> info_${create}
#a=S66-4137_25PyridinePyridinepipi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4138_26UracilUracilpipi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4139_27BenzenePyridinepipi 
# echo $(Eab $a) >> info_${create}
#a=S66-4140_28BenzeneUracilpipi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4141_29PyridineUracilpipi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4142_30BenzeneEthene 
# echo $(Eab $a)  >> info_${create}
#a=S66-4143_31UracilEthene 
# echo $(Eab $a)  >> info_${create}
#a=S66-4144_32UracilEthyne 
# echo $(Eab $a)  >> info_${create}
#a=S66-4145_33PyridineEthene 
# echo $(Eab $a)  >> info_${create}
#a=S66-4146_34PentanePentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4147_35NeopentanePentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4148_36NeopentaneNeopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4149_37CyclopentaneNeopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4150_38CyclopentaneCyclopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4151_39BenzeneCyclopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4152_40BenzeneNeopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4153_41UracilPentane
# echo $(Eab $a)  >> info_${create}
#a=S66-4154_42UracilCyclopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4155_43UracilNeopentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4156_44EthenePentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4157_45EthynePentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4158_46PeptidePentane 
# echo $(Eab $a)  >> info_${create}
#a=S66-4159_47BenzeneBenzeneTS 
# echo $(Eab $a)  >> info_${create}
#a=S66-4160_48PyridinePyridineTS 
# echo $(Eab $a)  >> info_${create}
#a=S66-4161_49BenzenePyridineTS 
# echo $(Eab $a)  >> info_${create}
#a=S66-4162_50BenzeneEthyneCHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4163_51EthyneEthyneTS 
# echo $(Eab $a)  >> info_${create}
#a=S66-4164_52BenzeneAcOHOHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4165_53BenzeneAcNH2NHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4166_54BenzeneWaterOHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4167_55BenzeneMeOHOHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4168_56BenzeneMeNH2NHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4169_57BenzenePeptideNHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4170_58PyridinePyridineCHN 
# echo $(Eab $a)  >> info_${create}
#a=S66-4171_59EthyneWaterCHO 
# echo $(Eab $a)  >> info_${create}
#a=S66-4172_60EthyneAcOHOHpi 
# echo $(Eab $a)  >> info_${create}
#a=S66-4173_61PentaneAcOH 
# echo $(Eab $a)  >> info_${create}
#a=S66-4174_62PentaneAcNH2 
# echo $(Eab $a)  >> info_${create}
#a=S66-4175_63BenzeneAcOH 
# echo $(Eab $a)  >> info_${create}
#a=S66-4176_64PeptideEthene 
# echo $(Eab $a)  >> info_${create}
#a=S66-4177_65PyridineEthyne 
# echo $(Eab $a)  >> info_${create}
#a=S66-4178_66MeNH2Pyridine 
# echo $(Eab $a)  >> info_${create}

}



function ALLS66 {
echo ALLS66
#CCSD(T) /CBS(haTZ) CP              MP2 /cc-pVTZ CP              
A[1]=S66-4113_01WaterWater;F[1]=-0.2149;E[1]=-0.1955
A[2]=S66-4114_02WaterMeOH;F[2]=-0.2469;E[2]=-0.2209
A[3]=S66-4115_03WaterMeNH2;F[3]=-0.3072;E[3]=-0.2771
A[4]=S66-4116_04WaterPeptide;F[4]=-0.3500;E[4]=-0.3148
A[5]=S66-4117_05MeOHMeOH;F[5]=-0.2532;E[5]=-0.2233
A[6]=S66-4118_06MeOHMeNH2;F[6]=-0.3351;E[6]=-0.2940
A[7]=S66-4119_07MeOHPeptide;F[7]=-0.3548;E[7]=-0.3120
A[8]=S66-4120_08MeOHWater;F[8]=-0.2182;E[8]=-0.1958
A[9]=S66-4121_09MeNH2MeOH;F[9]=-0.1327;E[9]=-0.1086
A[10]=S66-4122_10MeNH2MeNH2;F[10]=-0.1860;E[10]=-0.1506
A[11]=S66-4123_11MeNH2Peptide;F[11]=-0.2398;E[11]=-0.1931
A[12]=S66-4124_12MeNH2Water;F[12]=-0.3261;E[12]=-0.2907
A[13]=S66-4125_13PeptideMeOH;F[13]=-0.2742;E[13]=-0.2367
A[14]=S66-4126_14PeptideMeNH2;F[14]=-0.3331;E[14]=-0.2868
A[15]=S66-4127_15PeptidePeptide;F[15]=-0.3761;E[15]=-0.3257
A[16]=S66-4128_16PeptideWater;F[16]=-0.2232;E[16]=-0.1999
A[17]=S66-4129_17UracilUracilBP;F[17]=-0.7445;E[17]=-0.6626
A[18]=S66-4130_18WaterPyridine;F[18]=-0.3065;E[18]=-0.2768
A[19]=S66-4131_19MeOHPyridine;F[19]=-0.3328;E[19]=-0.2941
A[20]=S66-4132_20AcOHAcOH;F[20]=-0.8239;E[20]=-0.7506
A[21]=S66-4133_21AcNH2AcNH2;F[21]=-0.6991;E[21]=-0.6377
A[22]=S66-4134_22AcOHUracil;F[22]=-0.8413;E[22]=-0.7643
A[23]=S66-4135_23AcNH2Uracil;F[23]=-0.8283;E[23]=-0.7542
A[24]=S66-4136_24BenzeneBenzenepipi;F[24]=-0.2039;E[24]=-0.1260
A[25]=S66-4137_25PyridinePyridinepipi;F[25]=-0.2604;E[25]=-0.1756
A[26]=S66-4138_26UracilUracilpipi;F[26]=-0.4830;E[26]=-0.3578
A[27]=S66-4139_27BenzenePyridinepipi;F[27]=-0.2356;E[27]=-0.1535
A[28]=S66-4140_28BenzeneUracilpipi;F[28]=-0.3267;E[28]=-0.2217
A[29]=S66-4141_29PyridineUracilpipi;F[29]=-0.3741;E[29]=-0.2691
A[30]=S66-4142_30BenzeneEthene;F[30]=-0.1008;E[30]=-0.0555
A[31]=S66-4143_31UracilEthene;F[31]=-0.1739;E[31]=-0.1201
A[32]=S66-4144_32UracilEthyne;F[32]=-0.1912;E[32]=-0.1412
A[33]=S66-4145_33PyridineEthene;F[33]=-0.1225;E[33]=-0.0742
A[34]=S66-4146_34PentanePentane;F[34]=-0.1721;E[34]=-0.1018
A[35]=S66-4147_35NeopentanePentane;F[35]=-0.1160;E[35]=-0.0681
A[36]=S66-4148_36NeopentaneNeopentane;F[36]=-0.0753;E[36]=-0.0431
A[37]=S66-4149_37CyclopentaneNeopentane;F[37]=-0.1078;E[37]=-0.0618
A[38]=S66-4150_38CyclopentaneCyclopentane;F[38]=-0.1360;E[38]=-0.0807
A[39]=S66-4151_39BenzeneCyclopentane;F[39]=-0.1986;E[39]=-0.1282
A[40]=S66-4152_40BenzeneNeopentane;F[40]=-0.1561;E[40]=-0.1018
A[41]=S66-4153_41UracilPentane;F[41]=-0.2359;E[41]=-0.1489
A[42]=S66-4154_42UracilCyclopentane;F[42]=-0.2037;E[42]=-0.1275
A[43]=S66-4155_43UracilNeopentane;F[43]=-0.1756;E[43]=-0.1151
A[44]=S66-4156_44EthenePentane;F[44]=-0.0933;E[44]=-0.0555
A[45]=S66-4157_45EthynePentane;F[45]=-0.0912;E[45]=-0.0581
A[46]=S66-4158_46PeptidePentane;F[46]=-0.1957;E[46]=-0.1238
A[47]=S66-4159_47BenzeneBenzeneTS;F[47]=-0.1624;E[47]=-0.1117
A[48]=S66-4160_48PyridinePyridineTS;F[48]=-0.1904;E[48]=-0.1379
A[49]=S66-4161_49BenzenePyridineTS;F[49]=-0.1810;E[49]=-0.1300
A[50]=S66-4162_50BenzeneEthyneCHpi;F[50]=-0.1502;E[50]=-0.1148
A[51]=S66-4163_51EthyneEthyneTS;F[51]=-0.0720;E[51]=-0.0581
A[52]=S66-4164_52BenzeneAcOHOHpi;F[52]=-0.2277;E[52]=-0.1798
A[53]=S66-4165_53BenzeneAcNH2NHpi;F[53]=-0.2048;E[53]=-0.1627
A[54]=S66-4166_54BenzeneWaterOHpi;F[54]=-0.1547;E[54]=-0.1249
A[55]=S66-4167_55BenzeneMeOHOHpi;F[55]=-0.2065;E[55]=-0.1576
A[56]=S66-4168_56BenzeneMeNH2NHpi;F[56]=-0.1666;E[56]=-0.1188
A[57]=S66-4169_57BenzenePeptideNHpi;F[57]=-0.2688;E[57]=-0.2027
A[58]=S66-4170_58PyridinePyridineCHN;F[58]=-0.1894;E[58]=-0.1551
A[59]=S66-4171_59EthyneWaterCHO;F[59]=-0.1246;E[59]=-0.1107
A[60]=S66-4172_60EthyneAcOHOHpi;F[60]=-0.2181;E[60]=-0.1849
A[61]=S66-4173_61PentaneAcOH;F[61]=-0.1313;E[61]=-0.0809
A[62]=S66-4174_62PentaneAcNH2;F[62]=-0.1586;E[62]=-0.1006
A[63]=S66-4175_63BenzeneAcOH;F[63]=-0.1976;E[63]=-0.1392
A[64]=S66-4176_64PeptideEthene;F[64]=-0.1374;E[64]=-0.0979
A[65]=S66-4177_65PyridineEthyne;F[65]=-0.1824;E[65]=-0.1594
A[66]=S66-4178_66MeNH2Pyridine;F[66]=-0.1974;E[66]=-0.1483

#A[1]=S66-4113_01WaterWater;E[1]=-0.2060
#A[2]=S66-4114_02WaterMeOH;E[2]=-0.2320
#A[3]=S66-4115_03WaterMeNH2;E[3]=-0.2860
#A[4]=S66-4116_04WaterPeptide;E[4]=-0.3383
#A[5]=S66-4117_05MeOHMeOH;E[5]=-0.2363
#A[6]=S66-4118_06MeOHMeNH2;E[6]=-0.3073
#A[7]=S66-4119_07MeOHPeptide;E[7]=-0.3378
#A[8]=S66-4120_08MeOHWater;E[8]=-0.2078
#A[9]=S66-4121_09MeNH2MeOH;E[9]=-0.1204
#A[10]=S66-4122_10MeNH2MeNH2;E[10]=-0.1626
#A[11]=S66-4123_11MeNH2Peptide;E[11]=-0.2121
#A[12]=S66-4124_12MeNH2Water;E[12]=-0.2999
#A[13]=S66-4125_13PeptideMeOH;E[13]=-0.2506
#A[14]=S66-4126_14PeptideMeNH2;E[14]=-0.2987
#A[15]=S66-4127_15PeptidePeptide;E[15]=-0.3482
#A[16]=S66-4128_16PeptideWater;E[16]=-0.2119
#A[17]=S66-4129_17UracilUracilBP;E[17]=-0.7063
#A[18]=S66-4130_18WaterPyridine;E[18]=-0.2823
#A[19]=S66-4131_19MeOHPyridine;E[19]=-0.3004
#A[20]=S66-4132_20AcOHAcOH;E[20]=-0.7972
#A[21]=S66-4133_21AcNH2AcNH2;E[21]=-0.6764
#A[22]=S66-4134_22AcOHUracil;E[22]=-0.8121
#A[23]=S66-4135_23AcNH2Uracil;E[23]=-0.7997
#A[24]=S66-4136_24BenzeneBenzenepipi;E[24]=-0.0764
#A[25]=S66-4137_25PyridinePyridinepipi;E[25]=-0.1201
#A[26]=S66-4138_26UracilUracilpipi;E[26]=-0.3522
#A[27]=S66-4139_27BenzenePyridinepipi;E[27]=-0.1013
#A[28]=S66-4140_28BenzeneUracilpipi;E[28]=-0.1856
#A[29]=S66-4141_29PyridineUracilpipi;E[29]=-0.2344
#A[30]=S66-4142_30BenzeneEthene;E[30]=-0.0359
#A[31]=S66-4143_31UracilEthene;E[31]=-0.1163
#A[32]=S66-4144_32UracilEthyne;E[32]=-0.1354
#A[33]=S66-4145_33PyridineEthene;E[33]=-0.0535
#A[34]=S66-4146_34PentanePentane;E[34]=-0.1243
#A[35]=S66-4147_35NeopentanePentane;E[35]=-0.0854
#A[36]=S66-4148_36NeopentaneNeopentane;E[36]=-0.0569
#A[37]=S66-4149_37CyclopentaneNeopentane;E[37]=-0.0773
#A[38]=S66-4150_38CyclopentaneCyclopentane;E[38]=-0.0980
#A[39]=S66-4151_39BenzeneCyclopentane;E[39]=-0.1157
#A[40]=S66-4152_40BenzeneNeopentane;E[40]=-0.0947
#A[41]=S66-4153_41UracilPentane;E[41]=-0.1585
#A[42]=S66-4154_42UracilCyclopentane;E[42]=-0.1334
#A[43]=S66-4155_43UracilNeopentane;E[43]=-0.1244
#A[44]=S66-4156_44EthenePentane;E[44]=-0.0661
#A[45]=S66-4157_45EthynePentane;E[45]=-0.0582
#A[46]=S66-4158_46PeptidePentane;E[46]=-0.1440
#A[47]=S66-4159_47BenzeneBenzeneTS;E[47]=-0.0971
#A[48]=S66-4160_48PyridinePyridineTS;E[48]=-0.1251
#A[49]=S66-4161_49BenzenePyridineTS;E[49]=-0.1173
#A[50]=S66-4162_50BenzeneEthyneCHpi;E[50]=-0.1081
#A[51]=S66-4163_51EthyneEthyneTS;E[51]=-0.0595
#A[52]=S66-4164_52BenzeneAcOHOHpi;E[52]=-0.1787
#A[53]=S66-4165_53BenzeneAcNH2NHpi;E[53]=-0.1666
#A[54]=S66-4166_54BenzeneWaterOHpi;E[54]=-0.1278
#A[55]=S66-4167_55BenzeneMeOHOHpi;E[55]=-0.1563
#A[56]=S66-4168_56BenzeneMeNH2NHpi;E[56]=-0.1150
#A[57]=S66-4169_57BenzenePeptideNHpi;E[57]=-0.1934
#A[58]=S66-4170_58PyridinePyridineCHN;E[58]=-0.1599
#A[59]=S66-4171_59EthyneWaterCHO;E[59]=-0.1184
#A[60]=S66-4172_60EthyneAcOHOHpi;E[60]=-0.1958
#A[61]=S66-4173_61PentaneAcOH;E[61]=-0.0972
#A[62]=S66-4174_62PentaneAcNH2;E[62]=-0.1193
#A[63]=S66-4175_63BenzeneAcOH;E[63]=-0.1343
#A[64]=S66-4176_64PeptideEthene;E[64]=-0.1087
#A[65]=S66-4177_65PyridineEthyne;E[65]=-0.1613
#A[66]=S66-4178_66MeNH2Pyridine;E[66]=-0.1458
for((i=1;i<67;i++))
do
#F = CCSD(T) /CBS(haTZ) CP             E = MP2 /cc-pVTZ CP              
echo ${A[i]}' '${F[i]}' '$(diffEab ${A[i]})' '$(uETOT ${A[i]}) >> ../bio/E_ALLs66_CCSDT-CBShaTZ-CP_${create}
echo ${A[i]}' '${E[i]}' '$(diffEab ${A[i]})' '$(uETOT ${A[i]}) >> ../bio/E_ALLs66_MP2-cc-pVTZ-CP_${create}
done

}

function TS {
echo TS
echo "
http://comp.chem.umn.edu/db/dbs/htbh38.html
http://comp.chem.umn.edu/db/dbs/nhtbh38.html
					REF1		REF2 (B.E.)
					eV		eV		Fireball
TS1	H + HCl -> H2 + Cl	Vf  	$(EV 5.70)	$(EV 5.70)
				Vr  	$(EV 7.86)	$(EV 8.70)
TS2	OH + H2 -> H2O + H	Vf  	$(EV 4.90)	$(EV 5.10)	$(Barrera2 H2O HH H2)
				Vr  	$(EV 21.20)	$(EV 21.20)     $(Barrera TS2 OH H2)
TS3	CH3 + H2 -> CH4 + H	Vf  	$(EV 12.10)	$(EV 12.10)	$(Barrera TS3 CH3 H2)
				Vr  	$(EV 15.30)	$(EV 15.30)	$(Barrera TS3 CH4 H)
TS4	OH + CH4 -> H2O + CH3	Vf  	$(EV 6.50)	$(EV 6.70)	$(Barrera TS4 OH CH4)
				Vr  	$(EV 19.60)	$(EV 19.60)	$(Barrera TS4 H2O CH3)
TS5	H + H2 -> H2 + H	Vf  	$(EV 9.60)	$(EV 9.60)	$(Barrera TS5 H H2)
				Vr  	$(EV 9.60)	$(EV 9.60)	$(Barrera TS5 H2 H)
TS6	OH + NH3 -> H2O + NH2	Vf  	$(EV 3.00)	$(EV 3.20)	$(Barrera TS6 OH NH3)
				Vr  	$(EV 12.70)	$(EV 12.70)	$(Barrera TS6 H2O NH2)
TS7	HCl + CH3 -> CH4 + Cl	Vf  	$(EV 1.70)	$(EV 1.70)
				Vr  	$(EV 7.06)	$(EV 7.90)
TS8	OH + C2H6 -> H2O + C2H5	Vf  	$(EV 3.20)	$(EV 3.40)	$(Barrera TS8 OH C2H6)
				Vr  	$(EV 19.90)	$(EV 19.90)	$(Barrera TS8 H2O C2H5)
TS9	F + H2 -> HF + H	Vf  	$(EV 1.42)	$(EV 1.80)	
				Vr  	$(EV 33.40)	$(EV 33.40)
TS10	O + CH4 -> OH + CH3	Vf  	$(EV 13.47)	$(EV 13.70)	$(Barrera TS10 O CH4)
				Vr  	$(EV 7.90)	$(EV 8.10)	$(Barrera TS10 OH CH3)
TS11	H + PH3 -> H2 + PH2	Vf  	$(EV 3.10)	$(EV 3.10)
				Vr  	$(EV 23.20)	$(EV 23.20)
TS12	H + HO -> H2 + O	Vf  	$(EV 10.50)	$(EV 10.70)	$(Barrera TS12 H OH)
				Vr  	$(EV 12.87)	$(EV 13.10)	$(Barrera TS12 H2 O)
TS13	H + H2S -> H2 + HS	Vf  	$(EV 3.50)	$(EV 3.50)
				Vr  	$(EV 16.76)	$(EV 17.30)
TS14	O + HCl -> OH + Cl	Vf  	$(EV 9.57)	$(EV 9.80)
				Vr  	$(EV 9.36)	$(EV 10.40)
TS15	CH3 + NH2 -> CH4 + NH	Vf  	$(EV 8.00)	$(EV 8.00)	$(Barrera TS15 CH3 NH2)
				Vr  	$(EV 22.40)	$(EV 22.40)	$(Barrera TS15 CH4 NH)
TS16	C2H5 + NH2 -> C2H6 + NH	Vf  	$(EV 7.50)	$(EV 7.50)	$(Barrera TS16 C2H5 NH2)
				Vr  	$(EV 18.30)	$(EV 18.30)	$(Barrera TS16 C2H6 NH)
TS17	NH2 + C2H6 -> NH3 + C2H5Vf  	$(EV 10.40)	$(EV 10.40)	$(Barrera TS17 NH2 C2H6)
				Vr  	$(EV 17.40)	$(EV 17.40)	$(Barrera TS17 NH3 C2H5)
TS18	NH2 + CH4 -> NH3 + CH3	Vf  	$(EV 14.50)	$(EV 14.50)	$(Barrera TS18 NH2 CH4)
				Vr  	$(EV 17.80)	$(EV 17.80)	$(Barrera TS18 NH3 CH3)
TS19	stranscis-C5H8->s-transcis-C5H8 $(EV 38.40)	$(EV 38.40)
				Vr  	$(EV 38.40)	$(EV 38.40)
TS21	H + N2O -> OH + N2	Vf 	$(EV 17.13)	$(EV 17.13)	$(Barrera TS21 H N2O)
				Vr 	$(EV 82.27)	$(EV 82.47)	$(Barrera TS21 OH N2)
TS22	H + FH -> HF + H	Vf 	$(EV 42.18)	$(EV 42.18)
				Vr 	$(EV 42.18)	$(EV 42.18)
TS23	H + ClH -> HCl + H	Vf 	$(EV 18.00)	$(EV 18.00)
				Vr 	$(EV 18.00)	$(EV 18.00)
TS24	H + FCH3 -> HF + CH3	Vf 	$(EV 30.38)	$(EV 30.38)
				Vr 	$(EV 57.02)	$(EV 57.02)
TS25	H + F2 -> HF + F	Vf 	$(EV 2.27)	$(EV 2.27)
				Vr 	$(EV 105.80)	$(EV 105.80)
TS26	CH3 + FCl -> CH3F + Cl	Vf 	$(EV 6.75)	$(EV 6.75)
				Vr 	$(EV 59.16)	$(EV 60.00)
TS27	F- + CH3F -> FCH3 + F-	Vf 	$(EV -0.34)	$(EV -0.34)
				Vr 	$(EV -0.34)	$(EV -0.34)
TS28	F-...CH3F -> FCH3...F-	Vf 	$(EV 13.38)	$(EV 13.38)
				Vr 	$(EV 13.38)	$(EV 13.38)
TS29	Cl- + CH3Cl->ClCH3+Cl-	Vf 	$(EV 3.10)	$(EV 3.10)
				Vr 	$(EV 3.10)	$(EV 3.10)
TS30	Cl-...CH3Cl->ClCH3...Cl-Vf 	$(EV 13.41)	$(EV 13.41)
				Vr 	$(EV 13.41)	$(EV 13.41)
TS31	F- + CH3Cl -> FCH3 + Cl-Vf 	$(EV -12.54)	$(EV -12.54)
				Vr 	$(EV 20.11)	$(EV 20.11)
TS32	F-...CH3Cl -> FCH3...Cl-Vf 	$(EV 3.44)	$(EV 3.44)
				Vr 	$(EV 29.42)	$(EV 29.42)
TS33	OH- + CH3F -> HOCH3 + F-Vf 	$(EV -2.44)	$(EV -2.44)
				Vr 	$(EV 17.66)	$(EV 17.66)
TS34	OH-...CH3F -> HOCH3...F-Vf 	$(EV 10.96)	$(EV 10.96)
				Vr 	$(EV 47.20)	$(EV 47.20)
TS35	H + N2 -> HN2		Vf 	$(EV 14.36)	$(EV 14.36)	$(Barrera TS35 H N2)
				Vr 	$(EV 10.61)	$(EV 10.61)	$(Barrera TS35 HN2)
TS36	H + CO -> HCO		Vf 	$(EV 3.17)	$(EV 3.17)	$(Barrera TS36 H CO)
				Vr 	$(EV 22.68)	$(EV 22.68)	$(Barrera TS36 HCO)
TS37	H + C2H4 -> CH3CH2	Vf 	$(EV 1.72)	$(EV 1.72)	$(Barrera TS37 H C2H4)
				Vr 	$(EV 41.75)	$(EV 41.75)	$(Barrera TS37 CH3CH2)
TS38	CH3 + C2H4 -> CH3CH2CH2	Vf 	$(EV 6.85)	$(EV 6.85)	$(Barrera TS38 CH3 C2H4)
				Vr 	$(EV 32.97)	$(EV 32.97)	$(Barrera TS38 CH3CH2CH2)
TS39	HCN -> HNC		Vf 	$(EV 48.07)	$(EV 48.07)	$(Barrera TS39 HCN FFF)
		  		Vr 	$(EV 32.82)	$(EV 32.82)	$(Barrera TS39 HNC FFF)
">> info_${create}
}

function PA {
echo PA
echo "
http://comp.chem.umn.edu/db/dbs/pa8.html
Proton Affinities
		REF1		REF2 (B.E.)
NH3	 $(EV 211.90)	 $(EV 211.90)	$(diffEH NH3 NH4+)
H2O	 $(EV 171.80)	 $(EV 171.80)	$(diffEH H2O H3O+)
C2H2	 $(EV 156.60)	 $(EV 156.60)	$(diffEH C2H2 C2H3+)
SiH4	 $(EV 156.50)	 $(EV 156.50)	
PH3	 $(EV 193.10)	 $(EV 193.10)	 
H2S	 $(EV 173.70)	 $(EV 173.70)	 
HCl	 $(EV 137.10)	 $(EV 137.10)	
H2	 $(EV 105.90)	 $(EV 105.90)	$(diffEH H2 H3+)
">> info_${create}

echo NH3 $(EV 211.90) $(diffEH NH3 NH4+) >../bio/PA_${create}
echo H2O $(EV 171.80) $(diffEH H2O H3O+) >>../bio/PA_${create} 
echo C2H2 $(EV 156.60) $(diffEH C2H2 C2H3+) >>../bio/PA_${create}
echo H2 $(EV 105.90) $(diffEH H2 H3+) >>../bio/PA_${create}

}


function IOHB {
echo IOHB
echo '
ref : http://www.begdb.com/index.php?action=oneDataset&id=29&state=show&order=ASC&by=name_m&method=    (1.00)
------------------ TEST IOHB  ----------------
test              Ea+Eb      Etot	Etot-Ea-Eb	ref.E	Dc	 Dl	 ef. Dc Dl'   >>  info_${create} 
a=3361_01acetatemethanol100- ; b=C2O2H3- ; c=COH4
e=-0.8564
d1=$(dC $a 4 8)
d1r=1.692
d2=$(dC $a 3 12)
d2r=3.149
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


#echo $(E1_23 $a $b $c)$'\t'-0.8564$'\t'$(dC $a 4 8)$'\t'$(dC $a 3 12)$'\t'$'\t'$'\t'1.692	3.149 >> info_${create}
a=3369_02acetatewater100- ; b=C2O2H3- ; c=H2O
#echo $(E1_23 $a $b $c)$'\t'-0.9132$'\t'$(dC $a 4 8)$'\t'$(dC $a 3 10)$'\t'$'\t'$'\t'2.042	1.999  >> info_${create}
e=-0.9132
d1=$(dC $a 4 8)
d1r=2.042
d2=$(dC $a 3 10)
d2r=1.999
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3377_03acetatemethylamine100- ; b=C2O2H3- ; c=CNH5
#echo $(E1_23 $a $b $c)$'\t'-0.4970$'\t'$(dC $a 4 8)$'\t'$(dC $a 3 12)$'\t'$'\t'$'\t'1.995	3.164  >> info_${create}
e=-0.4970
d1=$(dC $a 4 8)
d1r=1.995
d2=$(dC $a 3 12)
d2r=3.164
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3385_04methylammoniumformaldehyde100+ ; b=CNH6+ ; c=COH2
#echo $(E1_23 $a $b $c)$'\t'-0.8283$'\t'$(dC $a 9 5)$'\t'$'\t'$'\t'$'\t'1.729  >> info_${create}
e=-0.8283
d1=$(dC $a 9 5)
d1r=1.729
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3393_05methylammoniummethylamine100+  ; b=CNH6+ ; c=CNH5
#echo $(E1_23 $a $b $c)$'\t'-1.2385$'\t'$(dC $a 2 5)$'\t'$'\t'$'\t'$'\t'1.106  >> info_${create}
e=-1.2385
d1=$(dC $a 2 5)
d1r=1.106
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3401_06methylammoniummethanol100+ ;  b=CNH6+ ; c=COH4 
##echo $(E1_23 $a $b $c)$'\t'-0.9206$'\t'$(dC $a 5 9)$'\t'$'\t'$'\t'$'\t'1.642  >> info_${create}
e=-0.9206
d1=$(dC $a 5 9)
d1r=1.642
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3409_07methylammoniumwater100+ ;  b=CNH6+ ; c=H2O
#echo $(E1_23 $a $b $c)$'\t'-0.8027$'\t'$(dC $a 5 9)$'\t'$'\t'$'\t'$'\t'1.705  >> info_${create}
e=-0.8027
d1=$(dC $a 5 9)
d1r=1.705
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3417_08guanidiniumformaldehyde100+ ;  b=CN3H6+ ; c=COH2
#echo $(E1_23 $a $b $c)$'\t'-0.7845$'\t'$(dC $a 11 7)$'\t'$(dC $a 11 8)$'\t'$'\t'$'\t'2.012	2.012 >> info_${create}
e=-0.7845
d1=$(dC $a 11 7)
d1r=2.012
d2=$(dC $a 11 8)
d2r=2.012
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3425_09guanidiniummethylamine100+ ;  b=CN3H6+ ; c=CNH5
#echo $(E1_23 $a $b $c)$'\t'-0.8760$'\t'$(dC $a 11 9)$'\t'$'\t'$'\t'$'\t'1.787 >> info_${create}
e=-0.8760
d1=$(dC $a 11 9)
d1r=1.787
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3433_10guanidiniummethanol100+ ; b=CN3H6+ ; c=COH4
#echo $(E1_23 $a $b $c)$'\t'-0.8582$'\t'$(dC $a 11 9)$'\t'$(dC $a 11 7)$'\t'$'\t'$'\t'1.977	1.979 >> info_${create}
e=-0.8582
d1=$(dC $a 11 9)
d1r=1.977
d2=$(dC $a 11 7)
d2r=1.979
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3441_11guanidiniumwater100+ ; b=CN3H6+ ; c=H2O
#echo $(E1_23 $a $b $c)$'\t'-0.7576$'\t'$(dC $a 11 10)$'\t'$(dC $a 11 7)$'\t'$'\t'$'\t'2.034	2.033 >> info_${create}
e=-0.7576
d1=$(dC $a 11 10)
d1r=2.034
d2=$(dC $a 11 7)
d2r=2.033
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3449_12imidazoliumformaldehyde100+ ; b=C3N2H5+ ; c=COH2
#echo $(E1_23 $a $b $c)$'\t'-0.7116$'\t'$(dC $a 11 7)$'\t'$'\t'$'\t'$'\t'1.745  >> info_${create}
e=-0.7116
d1=$(dC $a 11 7)
d1r=1.745
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3457_13imidazoliummethylamine100+ ; b=C3N2H5+ ; c=CNH5
#echo $(E1_23 $a $b $c)$'\t'-1.1266$'\t'$(dC $a 11 7)$'\t'$'\t'$'\t'$'\t'1.596 >> info_${create}
e=-1.1266
d1=$(dC $a 11 7)
d1r=1.596
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3465_14imidazoliummethanol100+ ; b=C3N2H5+ ; c=COH4
#echo $(E1_23 $a $b $c)$'\t'-0.8200$'\t'$(dC $a 11 7)$'\t'$'\t'$'\t'$'\t'1.639 >> info_${create}
e=-0.8200
d1=$(dC $a 11 7)
d1r=1.639
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#echo $a' '${d2r}' '${d2} >> ../bio/P_${create}


a=3473_15imidazoliumwater100+ ; b=C3N2H5+ ; c=H2O
#echo $(E1_23 $a $b $c)$'\t'-0.7151$'\t'$(dC $a 11 10)$'\t'$'\t'$'\t'$'\t'1.700 >> info_${create}
e=-0.7151
d1=$(dC $a 11 10)
d1r=1.700
d2=
d2r=
echo $(E1_23 $a $b $c)$'\t'${e}$'\t'${d1}$'\t'${d2}$'\t'$'\t'$'\t'${d1r}$'\t'${d2r} >> info_${create}
echo $a' '${e}' '$(diffE1_23 $a $b $c) >> ../bio/E_IOHB_${create}
echo $a' '${d1r}' '${d1}' '$(uETOT $a)>> ../bio/P_${create}
#cho $a' '${d2r}' '${d2} >> ../bio/P_${create}


echo "
">>  info_${create}
}


########## SALIDA ###############


if test -f ${create}.tar.gz
then
rm -fr info_${create}
gunzip ${create}.tar.gz > /dev/null
tar -xvf ${create}.tar > /dev/null
rm -fr ${create}.tar
fi


if ! test -d ${create}
then
mkdir $create
fi


#si no existe lo crea ...
if ! test -d ${d}FIX; then mkdir ${d}FIX; fi


mv FIX $create
mv FIX/* ${create}/FIX
mv *${create}.xyz $create
mv *${create}.bas $create
mv *${create}.lvs $create
mv *${create}.dat $create
mv *${create}.agr $create
mv *${create}.kpts $create
mv CHARG*${create} $create




#%%%%%%%%%%%%%%% ANALISIS %%%%%%%%%%%%%%%%%%%
cd $create
rm -fr info_${create}
INTR
TESTBIO
ALLS66
S66
#TS
IOHB
PA
mkdir bio
cd FIX
rm -fr info_FIX_${create}
INTR
TESTBIO
S66
ALLS66
#TS
IOHB
PA
cd ..
if test -f ../bio_FIX
then
mkdir ../bio_FIX
fi
mv bio/* ../bio_FIX
rm -fr bio

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cp info_${create} ../
cp FIX/info_${create} ../info_FIX_${create}
cd ..

tar -cf ${create}.tar $create
gzip ${create}.tar
rm -fr $create



