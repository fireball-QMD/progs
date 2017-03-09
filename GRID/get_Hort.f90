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


! ortho_Hk.f90
! Program Description
! ===========================================================================
!       This subroutine orthogonalize H via Lowdin transformation
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_Hort (icluster, ikpoint, sks, iwrtdos,  &
     &                      iwrthop, iwrtatom, itrans, igap)

        use configuration
        use dimensions
        use interactions
        use neighbor_map
        use charges
        use density
        use transport
        use hartree_fock

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: iwrtdos
        integer, intent (in) :: iwrthop
        integer, intent (in) :: iwrtatom
        integer, intent (in) :: itrans
        integer, intent (in) :: igap
        real, intent (in), dimension (3) :: sks


! Local Parameters and Data Declaration
! ===========================================================================
        real*8, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer info
        integer inu
        integer in1, in2
        integer ineigh
        integer jatom
        integer jmu
        integer jnu
        integer mbeta
        integer lm
        integer issh
        integer mineig

        real dot
        real*8, dimension (norbitals) :: eigen
        real, dimension (3) :: vec
        real*8 sqlami
        real*8, dimension (norbitals) :: slam

        complex*16 a0
        complex*16 a1
        complex*16 phase

! A bunch of memory to be used in many ways
        complex*16, dimension (:, :), allocatable :: xxxx
        complex*16, dimension (:, :), allocatable :: yyyy
        complex*16, dimension (:, :), allocatable :: zzzz


! work vector for cheev/cheevd
        complex*16, allocatable, dimension (:) :: work
        real*8, allocatable, dimension (:) :: rwork
        integer, allocatable, dimension (:) :: iwork
        integer lwork, lrwork, liwork

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)

        if (wrtout) then
          write (*,*) ''
          write (*,*) '   ---       Orthogonalize H(k)     ---'
        end if

        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))
        lwork = 1
        allocate (work(lwork))
        lrwork = 3*norbitals - 2
        allocate (rwork(lrwork))

! We built H/S(K) going through the list of atoms and their neighbors.
! How do we know where the particular (i,j) spot belongs in the grand
! matrix?  This should help you understand the degelec shifting business:
!
!                  atom 1         atom2            atom  3
!
!              nu1 nu2 nu3   nu1 nu2 nu3 nu4     nu1 nu2 nu3
!
!                           _________________
!         mu1               |               |
!  atom 1 mu2               |    H(1,2)     |
!         mu3               |               |
!                           -----------------
!         mu1
!  atom 2 mu2
!         mu3
!         mu4
!
!         mu1
!  atom3  mu2
!         mu3
!
! to the (1,2) portion at the right place we use degelec(iatom), which is
! passed, it remembers how many orbitals one must skip to get to the spot
! reserved for iatom, e.g. in this example degelec(1)=0, degelc(2)=3.

!
! COMPUTE S(K) AND H(K)
! ****************************************************************************
! Find the overlap and Hamiltonian matrices s(mu,nu,i,j) and h(mu,nu,i,j)
! in k space. Here iatom is an atom in the central cell and jatom a
! neighboring atom. The k-space matrix is found from the real space matrix by:
! s(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) s((0,iatom), (l,jatom))
! h(k) = sum(l) exp(iatom*k - dot - (l + bj - bi)) h((0,iatom), (l,jatom))

! Initialize to zero first
        zzzz = a0
        yyyy = a0
        xxxx = a0

! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
          dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
          phase = cmplx(cos(dot),sin(dot))
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
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
            yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*vnl(imu,inu,ineigh,iatom)
           end do ! do imu
          end do ! do inu
         end do ! do inegh

        end do ! do iatom


        if  ( igap .eq. 3 ) then
          do iatom = 1, natoms
            in1 = imass(iatom)
! Now loop over all atoms, because there are non-neighbours important contributions
            do jatom = 1, natoms
              in2 = imass(jatom)
! So this matrix element goes in the i,j slot.
              vec(:) = ratom(:,jatom) - ratom(:,iatom)
              dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
              phase = cmplx(cos(dot),sin(dot))
              do inu = 1, num_orb(in2)
                jnu = inu + degelec(jatom)
                do imu = 1, num_orb(in1)
                  jmu = imu + degelec(iatom)
                  yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*hs_mat(jmu,jnu)
                end do
              end do
            end do
          end do
        end if


! xxxx = unused (used as complex workspace in cheev call below)
! yyyy = Hamiltonian in AO basis

! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
! NPA


! Call the diagonalizer
        if (wrtout) then
           write (*,*) ' Call diagonalizer for overlap. '
           write (*,*) '                  The overlap eigenvalues: '
           write (*,*) ' ******************************************************* '
        end if

! first find optimal working space
        call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work,      &
                       -1, rwork, info)
! resize working space
        lwork = work(1)
        deallocate (work)
        allocate (work(lwork))
! diagonalize the overlap matrix with the new working space
        call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, work,       &
     &               lwork, rwork , info)

        if (info .ne. 0) call diag_error (info, 0)

! writeout eigenvalues (min,max)
        write (*,100) slam(1), slam(norbitals)

! xxxx = unused
! zzzz = Overlap eigenvectors
! yyyy = Hamiltonian


! CHECK THE LINEAR DEPENDENCE
! ****************************************************************************
! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
! set by overtol.

! Determine the smallest active eigenvector
        mineig = 0
        do imu = 1, norbitals
         if (slam(imu) .lt. overtol) mineig = imu
        end do

! You can specify a specific number of orbitals to drop with this
! next line, by uncommenting it.
! mineig = 0  {Don't drop any}

        mineig = mineig + 1
        norbitals_new = norbitals + 1 - mineig

        if (norbitals_new .ne. norbitals) then
         write (*,*) '  '
         write (*,*) ' WARNING. ############################ '
         write (*,*) ' Linear dependence encountered in basis set. '
         write (*,*) ' An overlap eigenvalue is very small. '
         write (*,*) norbitals - norbitals_new, ' vectors removed. '
         write (*,*) ' Spurious orbital energies near zero will '
         write (*,*) ' appear as a result of dropping these orbitals'
         write (*,*) ' You can change this by adjusting overtol in '
         write (*,*) ' kspace.f '
         write (*,*) '  '
         write (*,*) '            The overlap eigenvalues: '
         write (*,*) ' ********************************************** '
         write (*,200) (slam(imu), imu = 1, norbitals)
         write (*,*) ' '

         do imu = mineig, norbitals
          jmu = imu - mineig + 1
          zzzz(:,jmu) = zzzz(:,imu)
          slam(jmu) = slam(imu)
         end do
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
        do imu = 1, norbitals_new
         sqlami = slam(imu)**(-0.25d0)
         zzzz(:,imu) = zzzz(:,imu)*sqlami
        end do

        call zgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz,&
     &              norbitals, zzzz, norbitals, a0, xxxx, norbitals)



! xxxx = S^-1/2 in AO basis
! zzzz = Unused (used as temporary work area below)
! yyyy = Hamiltonian in AO basis


! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
! Set M=H*(S^-.5)
        call zhemm ( 'R', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )

! Set Z=(S^-.5)*M
        call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, zzzz, norbitals, a0, yyyy, norbitals )


! xxxx = S^-1/2 in AO basis
! zzzz = Unused (used as complex workspace in cheev call below)
! yyyy = Hamiltonian in the MO basis set


! CGP
        if(iwrtdos.ge.1 .or. iwrthop.ge.1 .or. iwrtatom.ge.1) then
            hamk(1:norbitals,1:norbitals) = yyyy(1:norbitals,1:norbitals)
        end if
! end CGP
       if (itrans .eq. 1) then
         write (*,*)
         write (*,*)  '  -- Store H(k) for transport calculation --'
         H_k(1:norbitals,1:norbitals,ikpoint) = yyyy(1:norbitals,1:norbitals)
       endif


! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

        deallocate (rwork)
        deallocate (work)

! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f16.8, &
     &              ' eigenvalue(norbitals) = ', f16.8)
200     format (8f14.6)
!300     format (2x, <norbitals>f12.4)
301     format (2x, 2i4,2f12.4)
302     format (2x, 4i4,3f12.4)

        return
      end subroutine  get_Hort


