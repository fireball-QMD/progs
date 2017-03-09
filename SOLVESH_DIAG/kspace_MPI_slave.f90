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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


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

        real*8, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer, external :: iceil
        integer icontext
        integer ierror
        integer imu
        integer info
        integer inu
        integer jatom
        integer Kscf
        integer liwork
        integer locc
        integer locr
        integer lwork
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
        integer, external :: numroc

        integer, dimension (desc_length) :: desc_x
        integer, dimension (desc_length) :: desc_y
        integer, dimension (desc_length) :: desc_z
        integer, dimension (:), allocatable :: iwork
        integer, dimension (5) :: mybuffer

        real*8 sqlami

        real*8, dimension (:), allocatable :: eigen
        real*8, dimension (:), allocatable :: slam
        real*8, dimension (:, :), allocatable, save :: sm12_save
        real*8, dimension (:, :), allocatable :: xxxx
        real*8, dimension (:, :), allocatable :: yyyy
        real*8, dimension (:, :), allocatable :: zzzz
        real*8, dimension (:), allocatable :: work

! Procedure
! ===========================================================================
100     continue

! Initialize BLACS
        call MPI_BCAST (mybuffer, 5, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        norbitals = mybuffer(1)
        kscf = mybuffer(2)
        nprow = mybuffer(3)
        npcol = mybuffer(4)
        nb = mybuffer(5)

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
        np0 = numroc(norbitals, nb, 0, 0, nprow)
        mq0 = numroc(norbitals, nb, 0, 0, npcol)
        lwork = 1
        lwork = 5*norbitals + (np0 + mq0 + nb)*nb + 2000*norbitals
!        liwork = 30*norbitals

! Allocate some arrays
        allocate (eigen (norbitals))
!        allocate (iwork (liwork))
        allocate (slam (norbitals))
        allocate (work (lwork))
        allocate (xxxx (locr, locc))
        allocate (yyyy (locr, locc))
        allocate (zzzz (locr, locc))
        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(1:locr, 1:locc))
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

         call pdsyev ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx, 1,&
     &                1, desc_x, work, -1, info)
         lwork = work(1)
! reallocate working array
         deallocate (work)
         allocate (work(lwork))
         write (*,*) 'lwork= ',lwork
         write (*,*) 'norbitals =',norbitals
! Call the diagonalizer
         call pdsyev ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx, 1,&
     &                1, desc_x, work, lwork, info)
!        call pdsyevd ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx,  &
!     &                1, 1, desc_x, work, lwork, iwork, liwork, info)

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

         call pdgemm ('N', 'C', norbitals, norbitals, norbitals, 1.0d0,      &
     &                xxxx,1, 1, desc_x, xxxx, 1, 1, desc_x, 0.0d0, yyyy,    &
     &                1, 1, desc_y)

! Now put S^-1/2 into sk^-1/2, (remembered for the duration of the scf cycle)
         sm12_save(1:locr,1:locc) = yyyy(1:locr,1:locc)
        else

! Now if not first iteration, then just restore S^-1/2
         yyyy(1:locr,1:locc) = sm12_save(1:locr,1:locc)
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
        call pdsymm ('R', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,      &
     &               desc_y, zzzz, 1, 1, desc_z, 0.0d0, xxxx, 1, 1, desc_x)

! Set Z=(S^-.5)*M
        call pdsymm ('L', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,      &
     &               desc_y, xxxx, 1, 1, desc_x, 0.0d0, zzzz, 1, 1, desc_z)

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = Hamiltonian in MO basis

! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
! Eigenvectors are needed to calculate the charges and for forces!
! first find optimal working space
        call pdsyev ('V', 'U', norbitals, zzzz, 1, 1, desc_z, eigen, xxxx,   &
     &               1, 1, desc_x, work, -1, info)
        lwork = work(1)
! reallocate working array
        deallocate (work)
        allocate (work(lwork))
        call pdsyev ('V', 'U', norbitals, zzzz, 1, 1, desc_z, eigen, xxxx,   &
     &                1, 1, desc_x, work, lwork, info)
!        call pdsyevd ('V', 'U', norbitals, zzzz, 1, 1, desc_z, eigen, xxxx,  &
!     &                1, 1, desc_x, work, lwork, iwork, liwork, info)

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
        call pdsymm ('L', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,      &
     &               desc_y, xxxx, 1, 1, desc_x, 0.0d0, zzzz, 1, 1, desc_z)
        call pclagetter (zzzz, desc_z, xxxx, norbitals)

! Exit BLACS/PBLAS
        call blacs_gridexit (icontext)

! Deallocate Arrays
! ===========================================================================
        deallocate (eigen)
!        deallocate (iwork)
        deallocate (slam)
        deallocate (work)
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

        go to 100

! Format Statements
! ===========================================================================

        end
