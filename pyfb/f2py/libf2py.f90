program libf2py
print*, 'Wrapping Fortran Python using F2PY'
end 

subroutine start()
  call cpu_time (time_begin)
  call initbasics ()
  call readdata ()
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end

subroutine f2py_initbasics()
   ! more or less initbasics ()
   use options
   use configuration
   use interactions
   use scf
   use integrals
   use outputs
   use kpoints
   use optimization
   use md
   use charges
   use barrier
   use transport
   use energy
   use neighbor_map
   use forces
   use mpi_main
   use hartree_fock

   implicit none
! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer in1
   integer icount
   integer counter
   integer counter_ini
   integer isorp
   integer ideriv
   integer issh
   integer numorbPP_max
   integer numorb
   integer l
   integer imu

   real distance
   real, dimension (3) :: vector

   logical file_exists
   call welcome
   call initconstants (sigma, sigmaold, scf_achieved)
   call diagnostics (ioff2c, ioff3c, itestrange, testrange)
   call readparam ()
   call readinfo ()
   if( iclassicMD > 0 ) call readdata_classicMD ()

end

subroutine do_readdata()
  !reads the different data file.
  call readdata ()
  ! call readdata_xczw ()
  ! call read_1c
  ! call read_2c 
  ! call read_3c
end

subroutine do_main_loop()
!call scf_loop_harris (1)
!integer itime_step
!logical scf_achieved 
!itime_step  = 1
!call assemble_h ()
!call diag_k ()
!call build_rho (itime_step)
!scf_achieved = .false.
!call scf_bcast(scf_achieved)
  call main_loop ()
!   call initcDFT()
end
