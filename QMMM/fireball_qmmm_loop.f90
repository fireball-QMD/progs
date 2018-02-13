        subroutine fireball_qmmm_loop(itime_step,qmcoords_recv,nclatoms_recv,clcoords_recv,escf_send,dxyzqm_send,dxyzcl_send)
  
        use mpi_main
        use interactions
        use MD
        use options
        use configuration
        use energy
        use forces

        implicit none
! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        integer, intent (in) :: itime_step
        real, dimension(3,natoms), intent (in) :: qmcoords_recv
        integer, intent (in) :: nclatoms_recv
        real, dimension(4,nclatoms_recv), intent (in) :: clcoords_recv
        integer k
        real , intent(out) :: escf_send
        real , dimension(3,natoms), intent (out) :: dxyzqm_send
        real , dimension(3,nclatoms_recv), intent(out) :: dxyzcl_send

! Procedure
! ===========================================================================
       
        ratom = qmcoords_recv
        qmmm_struct%qm_mm_pairs = nclatoms_recv


        if (itime_step .eq. 1) then
         allocate(qmmm_struct%qm_xcrd(4,nclatoms_recv))
         qmmm_struct%qm_xcrd = clcoords_recv
         allocate(qmmm_struct%dxyzcl(3,nclatoms_recv))
        else
         write(*,*) 'else allocate'
         deallocate(qmmm_struct%dxyzcl)
         deallocate(qmmm_struct%qm_xcrd)
         allocate(qmmm_struct%dxyzcl(3,nclatoms_recv))
         allocate(qmmm_struct%qm_xcrd(4,nclatoms_recv))
         qmmm_struct%qm_xcrd = clcoords_recv
        end if

        

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

        if (imdet .eq. 1 ) then
           call main_loop_MDET_qmmm (itime_step)
        else
           call main_loop_MD_qmmm (itime_step)
        endif

        escf_send = etot*23.061d0
        dxyzqm_send = -ftot*23.061d0 
        dxyzcl_send = qmmm_struct%dxyzcl

        !deallocate(qmmm_struct%qm_xcrd)
        !deallocate(qmmm_struct%dxyzcl)
        !deallocate(dxyzcl_send)

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
1001    format (2x, '    Number OpenMP threads: ', i3)
700     format (i2, 3(2x,f8.4))
300     format (2x, g10.10)
        
      end subroutine fireball_qmmm_loop
