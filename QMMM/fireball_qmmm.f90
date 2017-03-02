        subroutine fireball_qmmm (itime_step,escf)
  
        use mpi_main
        use interactions
        use MD
        use options

        implicit none
! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
	integer, intent (in) :: itime_step
	integer k
        real , intent(out) :: escf

! Procedure
! ===========================================================================
        call cpu_time (time_begin)
        

! ***************************************************************************
!                           MPI Specific Stuff
! ***************************************************************************
! Call the initializer to see if we are doing MPI or not.  If we did not
! compile with the MPI flag, then this routine calls a dummy subroutine.
        call init_MPI (iammaster, iammpi, my_proc, nprocs)
        wrtout = .true.
        if (.not. iammaster) then
          write (*,*) '------ SLAVE:  CALL KSPACE_SLAVE ------'
          call kspace_slave (nprocs, my_proc)
        endif


        write(*,*) 'llamar a coord_to_fireball'
        call coords_to_fireball

        if (imdet .eq. 1 ) then
           call main_loop_MDET_qmmm (itime_step)
        else
           call main_loop_MD_qmmm (itime_step)
        endif

        write(*,*) 'coords_forces_charges_to_amber'
        call coords_forces_charges_to_amber(escf)
! timer
        call cpu_time (time_end)
!       write (*,1001) omp_get_max_threads()
        write (*,*) ' FIREBALL RUNTIME FOR STEP: ',itime_step ,'=',time_end-time_begin,'[sec]'
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
1001    format (2x, '    Number OpenMP threads: ', i3)
700     format (i2, 3(2x,f8.4))
300     format (2x, g10.10)
        
      end subroutine fireball_qmmm
