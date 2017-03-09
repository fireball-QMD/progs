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


! kspace_MPI.f90
! Program Description
! ===========================================================================
!       This is part of the MPI version of kspace.f
! ===========================================================================
! Original code written by Otto F. Sankey with modification by Alex A. Demkov
! and Jose Ortega

! Code rewritten by:
! James P. Lewis
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
        subroutine kspace (nprocs, my_proc, Kscf, iqout, icluster,   &
     &                     iwrteigen, ikpoint, sks,            &
     &                     nkpoints, iwrtdos, iwrthop, iwrtatom, itrans)
        use configuration
        use density
        use dimensions
        use interactions
        use neighbor_map
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrteigen
        integer, intent (in) :: Kscf
        integer, intent (in) :: my_proc
        integer, intent (in) :: nkpoints
        integer, intent (in) :: nprocs
! CGP
        integer, intent (in) :: iwrtdos
        integer, intent (in) :: iwrthop
        integer, intent (in) :: iwrtatom
        integer, intent (in) :: itrans
! end CGP

        real, intent (in), dimension (3) :: sks

! Output
 

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: desc_length = 10

        real*8, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer, external :: iceil
        integer icontext
        integer ierror
        integer imu
        integer ineigh
        integer info
        integer inu
        integer in1
        integer in2
        integer ishort
        integer jatom
        integer jmu
        integer jnu
        integer liwork
        integer locc
        integer locr
        integer lwork
        integer mbeta
        integer mineig
        integer mq0
        integer mycol
        integer myrow
        integer nb
        integer np0
        integer npcol
        integer nprow
        integer, external :: numroc

        integer, dimension (desc_length) :: desc_x
        integer, dimension (desc_length) :: desc_y
        integer, dimension (desc_length) :: desc_z
!        integer, dimension (:), allocatable :: iwork   ! iwork vector for pssyev

        integer, dimension (5) :: mybuffer

        real*8 sqlami
        real*8 magnitude
        real*8, dimension (norbitals) :: eigen
        real*8, dimension (norbitals) :: slam
        real*8, dimension (:, :), allocatable, save :: sm12_save
        real*8, dimension (:), allocatable :: work   ! work vector for pssyev
        real*8, dimension (:, :), allocatable :: xxxx
        real*8, dimension (:, :), allocatable :: yyyy
        real*8, dimension (:, :), allocatable :: zzzz

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! For this version of kspace, we assume that the user is calculating only the
! gamma point.   It usually does not make sense to calculate something with
! a lot of kpoints if the system size is quite big.  Of course, if one were
! interested in metals, then many kpoints would be needed and for a large
! system it would still be quite useful to use the ScaLAPACK routines.
        magnitude = sqrt(sks(1)**2 + sks(2)**2 + sks(3)**2)
        if (magnitude .gt. 1.0d-3) then
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' You have chosen to do a calculation with k-points other '
         write (*,*) ' than just the gamma point.  Unfotunately, this version '
         write (*,*) ' of kspace calls the ScaLAPACK subroutine that assumes '
         write (*,*) ' only real matrices (i.e. for kpoints other than the '
         write (*,*) ' gamma point, the complex version of ScaLAPACK is needed.'
         write (*,*) ' We must stop here! '
         stop
        end if

        ishort = 1
        if (iwrteigen .eq. 1) ishort = 0

        write (*,*) '  '
        write (*,*) ' ******************************************************** '
        write (*,*) '  '
        write (*,*) '         Welcome to kspace -- ikpoint = ', ikpoint
        write (*,*) '  '
        write (*,*) ' ******************************************************** '

! Initialize BLACS
        if (nprocs .le. 8) then
         nprow = 1
         npcol = nprocs
        else if (nprocs .eq. 16) then
         nprow = 4
         npcol = 4
        else if (nprocs .eq. 32) then
         nprow = 4
         npcol = 8
        else if (nprocs .eq. 64) then
         nprow = 8
         npcol = 8
        else if (nprocs .eq. 128) then
         nprow = 8
         npcol = 16
        else if (nprocs .eq. 256) then
         nprow = 16
         npcol = 16
        else if (nprocs .eq. 512) then
         nprow = 16
         npcol = 32
        else if (nprocs .eq. 1024) then
         nprow = 32
         npcol = 32
        else if (nprocs .eq. 2048) then
         nprow = 32
         npcol = 64
        else
         nprow = int(sqrt(real(nprocs)))
         npcol = nprocs/nprow
        end if
        nb = 64 ! how the columns and rows are split

        mybuffer(1) = norbitals
        mybuffer(2) = kscf
        mybuffer(3) = nprow
        mybuffer(4) = npcol
        mybuffer(5) = nb

        call MPI_BCAST (mybuffer, 5, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        call blacs_pinfo (my_proc, nprocs)
        call blacs_get (0, 0, icontext)
        call blacs_gridinit (icontext, 'R', nprow, npcol)
        call blacs_gridinfo (Icontext, nprow, npcol, myrow, mycol)

        if (myrow .eq. -1) then
         write (*,*) 'kspace_MPI died in BLACS initing '
         stop
        end if

! Allocate memory
        locr = max(1, numroc(norbitals, nb, myrow, 0, nprow))
        locc = max(1, numroc(norbitals, nb, mycol, 0, npcol))
        np0 = numroc(norbitals, nb, 0, 0, nprow)
        mq0 = numroc(norbitals, nb, 0, 0, npcol)
!        lwork = norbitals + (np0 + mq0 + nb)*nb + 2000*norbitals
!        write (*,*) 'lwork= ',lwork
!        liwork = 30*norbitals

!        allocate (iwork (liwork))
        lwork = 1
        allocate (work (lwork))
        allocate (xxxx (norbitals,norbitals))
        allocate (yyyy (norbitals,norbitals))
        allocate (zzzz (norbitals,norbitals))
        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(1:locr, 1:locc))
        end if
        call descinit (desc_x, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 norbitals, info)
        call descinit (desc_y, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 norbitals, info)
        call descinit (desc_z, norbitals, norbitals, nb, nb, 0, 0, icontext, &
     &                 norbitals, info)

! COMPUTE S(K) AND H(K)
! ****************************************************************************
! Find the overlap and Hamiltonian matrices s(mu,nu,i,j) and h(mu,nu,i,j)
! in k space. Here iatom is an atom in the central cell and jatom a
! neighboring atom.
! Initialize to zero first
        xxxx = 0.0d0
        zzzz = 0.0d0

! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            xxxx(jmu,jnu) = xxxx(jmu,jnu) + s_mat(imu,inu,ineigh,iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + h_mat(imu,inu,ineigh,iatom)
           end do
          end do
         end do
!        end do

! Now loop over all neighbors jatom of iatom VNL
        do ineigh = 1, neighPPn(iatom)
          mbeta = neighPP_b(ineigh,iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + vnl(imu,inu,ineigh,iatom)
           end do ! do imu
          end do ! do inu
         end do ! do inegh

        end do ! do iatom


! xxxx = overlap (complete matrix) in AO basis
! yyyy = nothing
! zzzz = H (complete matrix) in AO basis

! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
        if (Kscf .eq. 1) then
         call pclaputter (yyyy, desc_y, xxxx, norbitals)

! Call the diagonalizer
         write (*,*) '            Call MPI_diagonalizer for overlap.           '
         write (*,*) '                The overlap eigenvalues:                 '
         write (*,*) ' ******************************************************* '

! first find optimal working space
         call pdsyev ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx, 1,&
     &                1, desc_x, work, -1, info)
         lwork = work(1)
! reallocate working array
         deallocate (work)
         allocate (work(lwork))
         write (*,*) 'lwork= ',lwork
         call pdsyev ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx, 1,&
     &                1, desc_x, work, lwork, info)
!         call pdsyevd ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx,  &
!     &                 1, 1, desc_x, work, lwork, iwork, liwork, info)
         write (*,*) '  leaving MPI_diagonalizer '
         if (info .ne. 0) call diag_error (info, 1)

         if (ishort .eq. 1) then
          write (*,100) slam(1), slam(norbitals)
         else
          write (*,200) (slam(imu), imu = 1, norbitals)
         end if

! xxxx = Overlap eigenvectors
! yyyy = nothing
! zzzz = H (complete matrix) in AO basis
!
! CHECK THE LINEAR DEPENDENCE
! ****************************************************************************
! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
! set by overtol. This parameter is now in input.inc

! Determine the smallest active eigenvector
         mineig = 0
         do imu = 1, norbitals
          if (slam(imu) .lt. overtol) mineig = imu
         end do
         mineig = mineig + 1
         norbitals_new = norbitals + 1 - mineig

         if (norbitals_new .ne. norbitals) then
          write (*,*) '  '
          write (*,*) ' WARNING. ############################ '
          write (*,*) ' Linear dependence encountered in basis set. An overlap '
          write (*,*) ' eigenvalue is very small. '
          write (*,*) norbitals - norbitals_new, ' vectors removed. '
          write (*,*) ' Spurious orbital energies near zero will '
          write (*,*) ' appear as a result of dropping these orbitals'
          write (*,*) ' You can change this by adjusting overtol in kspace: '
          write (*,*) ' It gets in via the input.inc file. '
          write (*,*) '  '
          write (*,*) '                The overlap eigenvalues: '
          write (*,*) ' ****************************************************** '
          write (*,200) (slam(imu), imu = 1, norbitals)
          slam(1:mineig-1) = 0.0d0
         end if

!
! CALCULATE (S^-1/2) --> sm1
! ****************************************************************************
! In a diagonal reperesentation (Udagger*S*U = s, s is a diagonal matrix)
! We just take the inverse of the square roots of the eigenvalues to get
! s^-1/2. Then we 'undiagonalize' the s^-1/2 matrix back to get
! S^-1/2 = U*s^-1/2*Udagger.
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

         call pdgemm ('N', 'C', norbitals, norbitals, norbitals, 1.0d0, xxxx,  &
     &             1, 1, desc_x, xxxx, 1, 1, desc_x, 0.0d0, yyyy, 1, 1, desc_y)

! Now put S^-1/2 into s(k)^-1/2, this will be remembered for the duration of
! the scf cycle.
         sm12_save(1:locr,1:locc) = yyyy(1:locr,1:locc)
        else

! Now if not first iteration
! Restore S^-1/2 from s(k)^-1/2,
         yyyy(1:locr,1:locc) = sm12_save(1:locr,1:locc)
        end if

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = H (complete matrix) in AO basis

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        call pclaputter (xxxx, desc_x, zzzz, norbitals)

! Set M=H*(S^-.5)
        call pdsymm ('R', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,    & 
    &                desc_y, xxxx, 1, 1, desc_x, 0.0d0, zzzz, 1, 1, desc_z)

! Set Z=(S^-.5)*M
        call pdsymm ('L', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,    &
    &                desc_y, zzzz, 1, 1, desc_z, 0.0d0, xxxx, 1, 1, desc_x)


! xxxx = H in the MO basis
! yyyy = S^-1/2 in AO basis
! zzzz = nothing

!
! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
        write (*,*) '  '
        write (*,*) ' Call diagonalizer for hamiltonian. '
        write (*,*) '                The energy eigenvalues: '
        write (*,*) ' ******************************************************* '

! Eigenvectors are needed to calculate the charges and for forces!
! first find optimal working space
        call pdsyev ('V', 'U', norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,   &
     &               1, 1, desc_z, work, -1, info)
        lwork = work(1)
! reallocate working array
        deallocate (work)
        allocate (work(lwork))
        write (*,*) 'lwork= ',lwork
        call pdsyev ('V', 'U', norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,   &
     &                1, 1, desc_z, work, lwork, info)
!        call pdsyevd ('V', 'U', norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,  &
!     &                1, 1, desc_z, work, lwork, iwork, liwork, info)

        if (info .ne. 0) call diag_error (info, 1)

! FIXME-Should only go up to norbitals_new, but we do not know what eigenvalues
! and eigenvectors correspond to the removed MO's.  Their eigenvalues will be
! very close to zero, but not exactly.  Also, we do not know if a real
! eigen value is near zero.

        if (ishort .eq. 1) then
         write (*,100) eigen(1), eigen(norbitals)
        else
         write (*,200) (eigen(imu), imu = 1, norbitals)
        end if

!
! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! Save the answer.
        eigen_k(1:norbitals,ikpoint) = eigen(1:norbitals)

        call pclagetter (zzzz, desc_z, xxxx, norbitals)
        if (iqout .ne. 2) then
         blowre(1:norbitals,1:norbitals,ikpoint) =                           &
     &    real(xxxx(1:norbitals,1:norbitals))
         if (icluster .ne. 1) blowim(1:norbitals,1:norbitals,ikpoint) = 0.0d0
        end if

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = Eigenvectors of H in MO basis
        call pdsymm ('L', 'U', norbitals, norbitals, 1.0d0, yyyy, 1, 1,     &
     &               desc_y, zzzz, 1, 1, desc_z, 0.0d0, xxxx, 1, 1, desc_x)

        call pclagetter (xxxx, desc_x, zzzz, norbitals)
        bbnkre(1:norbitals,1:norbitals,ikpoint) =                            &
     &   real(zzzz(1:norbitals,1:norbitals))
        if (icluster .ne. 1) bbnkim(1:norbitals,1:norbitals,ikpoint) = 0.0d0


! xxxx = bbnkre/im stuff
! yyyy = S^-1/2 in AO basis
! zzzz = Eigenvectors of H in MO basis

! Exit BLACS/PBLAS
        call blacs_gridexit (icontext)

! Deallocate Arrays
! ===========================================================================
!        deallocate (iwork)
        deallocate (work)
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f10.6,                              &
     &              ' eigenvalue(norbitals) = ', f10.6)
200     format (4(2x, f12.4))

        return
        end
