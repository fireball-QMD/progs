#!/bin/bash

tar -zxf Fdata.tar.gz
cd distorted_benzene
ln -s ../Fdata_HC_minimal Fdata
cp fireball.in-opt fireball.in
echo "running optimization of Benzene with FIRE & McWEDA"
../../fireball.x > relaxation.out
echo "optimization done"
rm CHARGES
cp fireball.in-el fireball.in
echo "running calculation of electronic structure with McWEDA"
../../fireball.x > el_McWEDA.out
echo "electronic structure with McWEDA calculated"
cp fireball.in-dos fireball.in
echo "running calculation of DOS with McWEDA"
../../fireball.x > dos_McWEDA.out
echo "DOS structure with McWEDA calculated"
cp fireball.in-orb fireball.in
echo "running calculation for ploting the real-space orbitals with McWEDA"
../../fireball.x > dos_McWEDA.out
echo "real-space orbitals with McWEDA plotted"
cp fireball.in-vib fireball.in
echo "running calculation of molecular vibrations with McWEDA"
../../fireball.x > vibrations.out
echo "vibrations with McWEDA calculated"
rm CHARGES
cp fireball.in-grid fireball.in
echo "running calculation of electronic structure on a grid"
../../fireball.x > on_grid.out
echo "electronic structure on a grid calculated"
rm CHARGES
mv answer.xyz answer_opt.xyz
cp answer.bas answer_opt.bas
cp fireball.in-dyn fireball.in
echo "running dft-molecular dynamics with McWEDA"
../../fireball.x > dynamics.out
mv answer.xyz answer_dyn.xyz
echo "All test with Benzene performed"
