#!/bin/bash

#export OMP_NUM_THREADS=1 #important for some calculations on super-computers
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
../../fireball.x > electronic.out
echo "electronic structure with McWEDA calculated"
cp fireball.in-dos fireball.in
echo "running calculation of DOS with McWEDA"
../../fireball.x > dos.out
echo "DOS structure with McWEDA calculated"
cp fireball.in-orb fireball.in
echo "running calculation for ploting the real-space orbitals with McWEDA"
echo "degenerated HOMO - orbitals 14,15"
echo "degenerated LUMO - orbitals 16,17"
../../fireball.x > orbitals.out
echo "real-space orbitals with McWEDA plotted"
cp fireball.in-coefs fireball.in
echo "running calculation of LCAO coefficients McWEDA"
../../fireball.x > coefficients.out
echo "The coefficients with McWEDA calculated"
cp fireball.in-vib fireball.in
echo "running calculation of molecular vibrations with McWEDA"
../../fireball.x > vibrations.out
echo "vibrations with McWEDA calculated"
mv answer.xyz answer_opt.xyz
cp answer.bas answer_opt.bas
cp fireball.in-dyn fireball.in
echo "running dft-molecular dynamics with McWEDA"
../../fireball.x > dynamics.out
mv answer.xyz answer_dyn.xyz

cd ../benzene_on_grid
ln -s ../Fdata_HC_minimal Fdata

cp fireball.in-coefs fireball.in
echo "running calculation of LCAO coefficietns on a grid"
../../fireball.x > coefficients.out
echo "The coefficients on a grid calculated"
cp fireball.in-dos fireball.in
echo "running calculation of DOS on a grid"
../../fireball.x > dos_grid.out
echo "DOS on a grid calculated"
cp fireball.in-orbg fireball.in
echo "running calculation for ploting the real-space orbitals on a grid"
echo "degenerated HOMO - orbitals 14,15"
echo "degenerated LUMO - orbitals 16,17"
../../fireball.x > orb_grid.out
echo "real-space orbitals  on a grid plotted"

echo "All tests with benzene performed"

cd ../graphene_2x2
ln -s ../Fdata_HC_minimal Fdata

cp fireball.in-opt fireball.in
echo "running optimization of graphene with FIRE & McWEDA"
../../fireball.x > relaxation.out
echo "optimization done"
rm CHARGES
cp fireball.in-el fireball.in
echo "running calculation of electronic structure with McWEDA"
../../fireball.x > electronic.out
echo "electronic structure with McWEDA calculated"
grep "Fermi Level" electronic.out |tail -1 > fermi.dat
cp fireball.in-dos fireball.in
echo "running calculation of DOS with McWEDA"
../../fireball.x > dos.out
echo "DOS structure with McWEDA calculated"

cp fireball.in-orb fireball.in
echo "running calculation for ploting the real-space density nearby the Fermi Level with McWEDA"
../../fireball.x > orbitals.out
echo "real-space density with McWEDA plotted"
cp fireball.in-coefs fireball.in
echo "running calculation of LCAO coefficients McWEDA"
../../fireball.x > coefficients.out
echo "The coefficients with McWEDA calculated"
cp fireball.in-vib fireball.in
echo "running calculation of molecular vibrations with McWEDA"
../../fireball.x > vibrations.out
echo "vibrations with McWEDA calculated"
mv answer.xyz answer_opt.xyz
cp answer.bas answer_opt.bas
cp fireball.in-dyn fireball.in
echo "running dft-molecular dynamics with McWEDA"
../../fireball.x > dynamics.out
mv answer.xyz answer_dyn.xyz

cd ../graphene_grid
ln -s ../Fdata_HC_minimal Fdata

cp fireball.in-coefs fireball.in
echo "running calculation of LCAO coefficietns on a grid"
../../fireball.x > coefficients.out
echo "The coefficients on a grid calculated"
grep "Fermi Level" coefficients.out |tail -1 > fermi.dat
cp fireball.in-dos fireball.in
echo "running calculation of DOS on a grid"
../../fireball.x > dos_grid.out
echo "DOS on a grid calculated"
rm CHARGES denmat.dat
cp fireball.in-orbg fireball.in
echo "running calculation for ploting the real-space density nearby the Fermi Level on a grid"
../../fireball.x > orb_grid.out
echo "real-space density on a grid plotted"

echo "All tests with graphene performed"

echo "tests done"
