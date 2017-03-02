! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
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

! kspace_MPI_slave.f90
! Program Description
! ===========================================================================
!       This is part of the MPI version of kspace.f
!
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
        subroutine kspace_slave (nprocs, my_proc)
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: desc_length = 10

        real, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iceil
        integer icontext
        integer ierror
        integer ikpoint
        integer imu
        integer info
        integer inu
        integer jatom
        integer Kscf
        integer locc
        integer locr
        integer mineig
        integer mq0
        integer mycol
        integer myrow
        integer nb
        integer nkpoints
        integer np0
        integer npcol
        integer nprow
        integer norbitals
        integer norbitals_new
        integer numroc
        integer wrk1
        integer wrk2
        integer wrk3

        integer, dimension (desc_length) :: desc_x
        integer, dimension (desc_length) :: desc_y
        integer, dimension (desc_length) :: desc_z
        integer, dimension (:), allocatable :: iwork
        integer, dimension (7) :: mybuffer

        real abstol
        real pslamch
        real sqlami

        real, dimension (:), allocatable :: eigen
        real*8, dimension (:), allocatable :: rwork
        real, dimension (:), allocatable :: slam

        complex a0
        complex a1

        complex, dimension (:), allocatable :: cwork
        complex, dimension (:, :), allocatable :: xxxx
        complex, dimension (:, :), allocatable :: yyyy
        complex, dimension (:, :), allocatable :: zzzz
        complex, dimension (:, :, :), allocatable, save :: sm12_save

        external iceil
        external numroc

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)

100     continue

! Initialize BLACS
        call MPI_BCAST (mybuffer, 7, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        norbitals = mybuffer(1)
        kscf = mybuffer(2)
        ikpoint = mybuffer(3)
        nprow = mybuffer(4)
        npcol = mybuffer(5)
        nb = mybuffer(6)
        nkpoints = mybuffer(7)

        call blacs_pinfo (my_proc, nprocs)
        call blacs_get (0, 0, icontext)
        call blacs_gridinit (icontext, 'R', nprow, npcol)
        call blacs_gridinfo (icontext, nprow, npcol, myrow, mycoL)
        if (myrow .eq. -1) then
         call MPI_FINALIZE (ierror)
         stop
        end if

! Allocate memory
        locr = max(1, numroc(norbitals, nb, myrow, 0, nprow))
        locc = max(1, numroc(norbitals, nb, mycol, 0, npcol))
        wrk1 = 30*norbitals
        np0 = numroc(norbitals, nb, 0, 0, nprow)
        mq0 = numroc(norbitals, nb, 0, 0, npcol)
        wrk2 = 4*norbitals + max(5*norbitals, np0*mq0) +                     &
     &         iceil(norbitals, nprow*npcol)*norbitals + 2000*norbitals
        wrk3 = norbitals + (np0 + mq0 + nb)*nb + 2000*norbitals

! Allocate some arrays
        allocate (cwork(wrk3))
        allocate (eigen(norbitals))
        allocate (slam(norbitals))
        allocate (iwork(wrk1))
        allocate (rwork(wrk2))
        allocate (xxxx(locr, locc))
        allocate (yyyy(locr, locc))
        allocate (zzzz(locr, locc))
        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(1:locr, 1:locc, 1:nkpoints))
        end if

        call descinit (desc_x, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 locr, info)
        call descinit (desc_y, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 locr, info)
        call descinit (desc_z, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 locr, info)

!
! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
        if (Kscf .eq. 1) then
         call pclaputter (yyyy, desc_y, xxxx, norbitals)

! Call the diagonalizer
!        call pcheevx ('V', 'A', 'U', norbitals, yyyy, 1, 1, desc_y, 0, 0, 0,&
!    &                 0, abstol, ijunk1, ijunk2, slam, -1.0, xxxx, 1, 1,    &
!    &                 desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1, ijunk3,&
!    &                 ijunk4, xjunk5, info)
         call pcheevd ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx,  &
     &                 1, 1, desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1,  &
     &                 info)

! xxxx = eigenvectors of overlap
! yyyy = nothing
! zzzz = nothing

! Overlap is the eigenvectors, slam is the eigenvalues, mr is destroyed
         if (info .ne. 0) then
          call MPI_FINALIZE (ierror)
          stop
         end if

! CHECK THE LINEAR DEPENDENCE
! ****************************************************************************
         mineig = 0
         do imu = 1, norbitals
          if (slam(imu) .lt. overtol) mineig = imu
         end do
         mineig = mineig + 1
         norbitals_new = norbitals + 1 - mineig

         if (norbitals_new .ne. norbitals) then
          slam(1:mineig-1) = 0.0d0
         end if

! CALCULATE (S^-1/2) --> sm1
! ****************************************************************************
! Note: We do S^-1/4 here, because the sqlami contribution get squared
! after it is combined with overlap.
         do imu = 1, norbitals
          if (slam(imu) .lt. overtol) then
           sqlami = 0.0d0
          else
           sqlami = slam(imu)**(-0.25d0)
          end if
          do inu = 1, norbitals
           call blacsaba (xxxx, inu, imu, desc_x, sqlami, mycol, myrow,      &
     &                    npcol, nprow)
          end do
         end do

         call pcgemm ('N', 'C', norbitals, norbitals, norbitals, a1, xxxx, 1,&
     &                1, desc_x, xxxx, 1, 1, desc_x, a0, yyyy, 1, 1, desc_y)

! Now put S^-1/2 into sk^-1/2, (remembered for the duration of the scf cycle)
         sm12_save(1:locr,1:locc,ikpoint) = yyyy(1:locr,1:locc)
        else

! Now if not first iteration, then just restore S^-1/2
         yyyy(1:locr,1:locc) = sm12_save(1:locr,1:locc,ikpoint)
        end if

! xxxx = nothing
! yyyy = S^1/2 in AO basis
! zzzz = nothing

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        call pclaputter (zzzz, desc_z, xxxx, norbitals)

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = Hamiltonian in AO basis

! Set M=H*(S^-.5)
        call pchemm ('R', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               zzzz, 1, 1, desc_z, a0, xxxx, 1, 1, desc_x)

! Set Z=(S^-.5)*M
        call pchemm ('L', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               xxxx, 1, 1, desc_x, a0, zzzz, 1, 1, desc_z)


! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = Hamiltonian in MO basis

! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
! Eigenvectors are needed to calculate the charges and for forces!
!       call pcheevx ('V', 'A', 'U', norbitals, zzzz, 1, 1, desc_z, 0, 0, 0, &
!    &                0, abstol, ijunk1, ijunk2, eigen, -1.0, xxxx, 1, 1,    &
!    &                desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1, ijunk3, &
!    &                ijunk4, xjunk5, info)
        call pcheevd ('V', 'U', norbitals, zzzz, 1, 1, desc_z, eigen, xxxx,  &
     &                1, 1, desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1,   &
     &                info)

! xxxx = Eigenvectors of H in MO basis
! yyyy = S^-1/2 in AO basis
! zzzz = Nothing
        if (info .ne. 0) then
         call MPI_FINALIZE (ierror)
         stop
        end if

! Save the answer
        call pclagetter (xxxx, desc_x, zzzz, norbitals)

! ****************************************************************************
        call pchemm ('L', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               xxxx, 1, 1, desc_x, a0, zzzz, 1, 1, desc_z)

        call pclagetter (zzzz, desc_z, xxxx, norbitals)

! Exit BLACS/PBLAS
        call blacs_gridexit (icontext)

! Deallocate Arrays
! ===========================================================================
        deallocate (cwork)
        deallocate (eigen)
        deallocate (iwork)
        deallocate (rwork)
        deallocate (slam)
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

        goto 100

! Format Statements
! ===========================================================================

        end
