#!/bin/bash

#export OMP_NUM_THREADS=1 #important for some calculations on super-computers
tar -zxf Fdata.tar.gz

cd distorted_benzene
ln -s ../Fdata_HC_minimal Fdata

cp fireball.in-opt fireball.in
echo "running optimization of Benzene with FIRE & McWEDA"
../../fireball.x > relaxation.out
echo "optimization done"