
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

../../fireball.x > out.$$
