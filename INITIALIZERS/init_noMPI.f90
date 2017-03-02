! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! init_noMPI.f90
! Program Description
! ===========================================================================
! Initializes the MPI space
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine init_MPI (iammaster, iammpi, my_proc, nprocs)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output
        integer, intent (out) :: my_proc
        integer, intent (out) :: nprocs
 
        logical, intent (out) :: iammaster
        logical, intent (out) :: iammpi
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
        iammaster = .true.
        iammpi = .false.
        my_proc = 0
        nprocs = 1
 
! Format Statements
! ===========================================================================
 
        return
        end

        subroutine MPI_COMM_RANK ()
        write (*,*) '  '
        write (*,*) ' In MPI_COMM_RANK '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine MPI_COMM_RANK
 
        subroutine bcast ()
        write (*,*) '  '
        write (*,*) ' In bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine bcast

