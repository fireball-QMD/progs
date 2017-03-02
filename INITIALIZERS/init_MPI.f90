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

! init_MPI.f90 
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
!
! Program Declaration
! ===========================================================================
        subroutine init_MPI (iammaster, iammpi, my_proc, nprocs)
        use mpi_declarations
        implicit none
 
        include 'mpif.h'
 
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
        integer ierr_mpi
        real dummy
        double precision d_dummy

        integer mpic_real, mpic_double, mpic_max, mpic_sum
        common /mpiconsts/ mpic_real, mpic_double, mpic_max, mpic_sum
 
! Procedure
! ===========================================================================
! TaoInitialize will call this
        call MPI_INIT (ierr_mpi)
        if (ierr_mpi .ne. 0) then
         write (*,*) ' mpi error MPI_INIT '
         stop
        end if
        call MPI_COMM_RANK (MPI_COMM_WORLD, my_proc, ierr_mpi)
        if (ierr_mpi .ne. 0) then
         write (*,*) ' mpi error MPI_COMM_RANK '
         stop
        end if
        call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr_mpi)
        if (ierr_mpi .ne. 0) then
         write (*,*) ' mpi error MPI_COMM_SIZE '
         stop
        end if
        iammaster = .false.
        if (my_proc .eq. 0) iammaster = .true.
        iammpi = .true.

! set up mpi_whatever_real so that we can use real*4 or real*8 without
! affecting the MPI calls
        if (kind(dummy) .eq. 8) then
         mpi_whatever_real = MPI_DOUBLE_PRECISION
        else
         mpi_whatever_real = MPI_REAL
        end if
        if (kind(d_dummy) .eq. 8) then
         mpi_whatever_double = MPI_DOUBLE_PRECISION
        else
         write (*,*) ' ERROR: No support currently for 16-byte reals.'
!        mpi_whatever_double = MPI_LONG_DOUBLE
        end if

! initialize common variables for use by f77 modules
        mpic_real = mpi_whatever_real
        mpic_double = mpi_whatever_double
        mpic_max = mpi_max
        mpic_sum = mpi_sum

! MPI init flag 
        mpi_on = 0
! Format Statements
! ===========================================================================
 
        return
        end
