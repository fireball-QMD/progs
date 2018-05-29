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
        subroutine kspace (natoms, nprocs, my_proc, Kscf, iqout, icluster,   &
     &                     iwrteigen, ikpoint, sks,  eigen_k,  &
     &                     nkpoints)
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
        integer, intent (in) :: natoms
        integer, intent (in) :: nkpoints
        integer, intent (in) :: nprocs

        real, intent (in), dimension (3) :: sks

! Output
       

        real, intent (inout), dimension (norbitals, nkpoints) :: eigen_k

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: desc_length = 10

        real, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iceil
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
        integer locc
        integer locr
        integer wrk1
        integer wrk2
        integer wrk3
        integer mbeta
        integer mineig
        integer mq0
        integer mycol
        integer myrow
        integer nb
        integer np0
        integer npcol
        integer nprow
        integer numroc

        integer, dimension (desc_length) :: desc_x
        integer, dimension (desc_length) :: desc_y
        integer, dimension (desc_length) :: desc_z
        integer, dimension (:), allocatable :: iwork
        integer, dimension (7) :: mybuffer

        real dot
        real sqlami

        real, dimension (norbitals) :: eigen
        real*8, dimension (:), allocatable :: rwork
        real, dimension (norbitals) :: slam
        real, dimension (3) :: vec

        complex a0
        complex a1
        complex phase

        complex, dimension (:), allocatable :: cwork
        complex, dimension (:,:), allocatable :: xxxx
        complex, dimension (:,:), allocatable :: yyyy
        complex, dimension (:,:), allocatable :: zzzz
        complex, dimension (:, :, :), allocatable, save :: sm12_save

        external iceil
        external numroc

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)
        ishort = 1
        if (iwrteigen .eq. 1) ishort = 0

        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        write (*,*) '  '
        write (*,*) '         Welcome to kspace -- ikpoint = ', ikpoint
        write (*,*) '  '
        write (*,*) ' ****************************************************** '

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
        mybuffer(3) = ikpoint
        mybuffer(4) = nprow
        mybuffer(5) = npcol
        mybuffer(6) = nb
        mybuffer(7) = nkpoints

        call MPI_BCAST (mybuffer, 7, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

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
        wrk1 = 30*norbitals
        np0 = numroc(norbitals, nb, 0, 0, nprow)
        mq0 = numroc(norbitals, nb, 0, 0, npcol)
        wrk2 = 4*norbitals + max(5*norbitals, np0*mq0) +                     &
     &         iceil(norbitals, nprow*npcol)*norbitals + 2000*norbitals
        wrk3 = norbitals + (np0 + mq0 + nb)*nb + 2000*norbitals

        allocate (cwork(wrk3))
        allocate (iwork(wrk1))
        allocate (rwork(wrk2))
        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))
        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(1:locr, 1:locc, 1:nkpoints))
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
! neighboring atom. The k-space matrix is found from the real space matrix by:
! s(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) s((0,iatom), (l,jatom))
! h(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) h((0,iatom), (l,jatom))

! Initialize to zero first
        xxxx = a0
        zzzz = a0

! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
          dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
          phase = cmplx(cos(dot),sin(dot))

          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            xxxx(jmu,jnu) = xxxx(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
           end do ! do inu
          end do ! do imu 
         end do ! do ineigh


! Now loop over all neighbors jatom of iatom VNL
        do ineigh = 1, neighPPn(iatom)
          mbeta = neighPP_b(ineigh,iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
          dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
          phase = cmplx(cos(dot),sin(dot))

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
!         write (*,*) ' Call diagonalizer for overlap. '
!         write (*,*) '                The overlap eigenvalues: '
!         write (*,*) ' ******************************************************* '
!
!        call pcheevx ('V', 'A', 'U', norbitals, yyyy, 1, 1, desc_y, 0, 0, 0,&
!    &                 0, abstol, ijunk1, ijunk2, slam, -1.0, xxxx, 1, 1,    &
!    &                 desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1, ijunk3,&
!    &                 ijunk4, xjunk5, info)
         call pcheevd ('V', 'U', norbitals, yyyy, 1, 1, desc_y, slam, xxxx,  &
     &                 1, 1, desc_x, cwork, wrk3, rwork, wrk2, iwork, wrk1,  &
     &                 info)
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

         call pcgemm ('N', 'C', norbitals, norbitals, norbitals, a1, xxxx, 1,&
     &                1, desc_x, xxxx, 1, 1, desc_x, a0, yyyy, 1, 1, desc_y)

! Now put S^-1/2 into s(k)^-1/2, this will be remembered for the duration of
! the scf cycle.
         sm12_save(1:locr,1:locc,ikpoint) = yyyy(1:locr,1:locc)
        else

! Now if not first iteration
! Restore S^-1/2 from s(k)^-1/2,
         yyyy(1:locr,1:locc) = sm12_save(1:locr,1:locc,ikpoint)
        end if

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = H (complete matrix) in AO basis

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        call pclaputter (xxxx, desc_x, zzzz, norbitals)
! Set M=H*(S^-.5)
        call pchemm ('R', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               xxxx, 1, 1, desc_x, a0, zzzz, 1, 1, desc_z)

! Set Z=(S^-.5)*M
        call pchemm ('L', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               zzzz, 1, 1, desc_z, a0, xxxx, 1, 1, desc_x)

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
!       call pcheevx ('V', 'A', 'U', norbitals, xxxx, 1, 1, desc_x, 0, 0, 0, &
!    &                0, abstol, ijunk1, ijunk2, eigen, -1.0, zzzz, 1, 1,    &
!    &                desc_z, cwork, wrk3, rwork, wrk2, iwork, wrk1, ijunk3, &
!    &                ijunk4, xjunk5, info)
        call pcheevd ('V', 'U', norbitals, xxxx, 1, 1, desc_x, eigen, zzzz,  &
     &                1, 1, desc_z, cwork, wrk3, rwork, wrk2, iwork, wrk1,  &
     &                info)

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
         if (icluster .ne. 1) blowim(1:norbitals,1:norbitals,ikpoint) =      &
     &    aimag(xxxx(1:norbitals,1:norbitals))
        end if

! xxxx = nothing
! yyyy = S^-1/2 in AO basis
! zzzz = Eigenvectors of H in MO basis

        call pchemm ('L', 'U', norbitals, norbitals, a1, yyyy, 1, 1, desc_y, &
     &               zzzz, 1, 1, desc_z, a0, xxxx, 1, 1, desc_x)

        call pclagetter (xxxx, desc_x, zzzz, norbitals)
        bbnkre(1:norbitals,1:norbitals,ikpoint) =                            &
     &   real(zzzz(1:norbitals,1:norbitals))
        if (icluster .ne. 1) bbnkim(1:norbitals,1:norbitals,ikpoint) =       &
     &   aimag(zzzz(1:norbitals,1:norbitals))

! xxxx = bbnkre/im stuff
! yyyy = S^-1/2 in AO basis
! zzzz = Eigenvectors of H in MO basis

! Exit BLACS/PBLAS
        call blacs_gridexit (icontext)

! Deallocate Arrays
! ===========================================================================
        deallocate (cwork)
        deallocate (iwork)
        deallocate (rwork)
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
