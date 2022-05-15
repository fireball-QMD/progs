program libf2py
print*, 'Wrapping Fortran Python using F2PY'
!call load_base()
!call load_data()
end 

subroutine start()
  call cpu_time (time_begin)
  call initbasics ()
  call readdata ()
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end

subroutine do_initbasics()
!   call initbasics ()
   call welcome
   call initconstants (sigma, sigmaold, scf_achieved)
   call diagnostics (ioff2c, ioff3c, itestrange, testrange)
   call readparam ()
   call readinfo ()
   ! call initconstants (sigma, sigmaold, scf_achieved) 
   ! Read the parameter file - fireball.param
   ! call readparam () 
   ! Read the info.dat
   ! call readinfo ()  
   ! Read data from the basis file - XXX.bas.
   ! call readbasis 
   ! call readlvs 
   ! call readfragments ()
   ! call initmasses 
   ! call readquench 
   ! call readbarrier 
   ! call make_munu (nspecies)
   ! call make_munuPP (nspecies)
   ! ....

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
