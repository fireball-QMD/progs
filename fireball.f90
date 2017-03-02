! info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! fireball.f90
! Program Description
! ===========================================================================
!       This is the main driver program for the FIREBALL package.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        program fireball
  
        use mpi_main
        use interactions


!        use omp_lib

        implicit none

! Parameters and Data Declaration
! ===========================================================================
 
! Variable Declaration and Description
! ===========================================================================
 
! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end

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
        end if
 
! Read basic informations about the task
        call initbasics ()
        
! Read data
        call readdata ()
            
! Executing time loop
        call main_loop ()

! timer
        call cpu_time (time_end)
!       write (*,1001) omp_get_max_threads()
        write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
1001    format (2x, '    Number OpenMP threads: ', i3)
 
        stop
      end program fireball
