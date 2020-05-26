
if test -d Fdata_HC_minimal
then
echo 'Fdata existe'
else
tar -xvf ../Fdata.tar.gz
fi

echo '
&OPTION
fdatalocation="'"Fdata_HC_minimal"'"
nstepf = 5000
iquench = -1
icluster = 1
iqout = 3
dt = 0.5
&END
&OUTPUT
iwrtxyz = 1
&END
' > fireball.in

echo '12
   6     -4.197174     -3.920062     -0.001077
   6     -3.132111     -3.058028      0.002942
   6     -3.346311     -1.704680      0.002885
   6     -4.625536     -1.213686      0.001660
   6     -5.690602     -2.075731     -0.001805
   6     -5.476402     -3.429078     -0.004659
   1     -6.310509     -4.036755     -0.008219
   1     -3.985726     -5.256165     -0.001792
   1     -2.242726     -3.399675      0.005306
   1     -2.295045     -0.853275      0.005357
   1     -4.774460     -0.272773      0.003081
   1     -6.953318     -1.590684     -0.003133
'> input.bas


../../fireball.x > out.$$
