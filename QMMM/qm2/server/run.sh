
 mpiexec -n 1 -env I_MPI_DEBUG 2 -env I_MPI_DEVICE sock ./fireball_server &

rm -fr outmd ; mpiexec -n 1  -env I_MPI_DEBUG 2 -env I_MPI_DEVICE sock /home/jesus/amber16/bin/sander.MPI -O -i amber.in -o outmd -p input.tpp -c ini.rst -x answer.mdcrd -r answer.rst -ref ref.rst ; tail outmd 
