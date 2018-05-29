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


! kspace_blas.f90
! Program Description
! ===========================================================================
!       This is a version of kspace.f that uses the blas library
! ===========================================================================
! Original code written by Otto F. Sankey with modification by Alex A. Demkov
! and Jose Ortega

! Code rewritten by:
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
        subroutine diag_Hk (iqout, icluster, iwrteigen, ikpoint, sks, nkpoints)

        use configuration
        use density
        use dimensions
        use interactions
        use neighbor_map
        use charges
        use tdse

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrteigen
        integer, intent (in) :: nkpoints

!        real, intent (in), dimension (3, natoms) :: ratom
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
        integer ishort
        integer jatom
        integer jmu
        integer jnu
        integer mbeta
        integer mineig
        integer lm
        integer issh

        real dot
        real*8, dimension (norbitals) :: eigen
        real, dimension (3) :: vec

        complex*16 a0
        complex*16 a1
        complex*16 phase

! A bunch of memory to be used in many ways
        complex*16, dimension (:, :), allocatable :: xxxx
        complex*16, dimension (:, :), allocatable :: yyyy
        complex*16, dimension (:, :), allocatable :: zzzz
!NPA
        complex*16, dimension (:, :), allocatable :: ssss
        real*8, dimension (:), allocatable :: ww

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

        ishort = 1
        if (iwrteigen .eq. 1) ishort = 0

!        if (wrtout) then
!          write (*,*) '  '
!          write (*,*) ' ****************************************************** '
!          write (*,*) '  '
!          write (*,*) '         Diagonalize H(k) -- ikpoint = ', ikpoint
!          write (*,*) '  '
!          write (*,*) ' ****************************************************** '
!        end if
!
        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))
!NPA
        if (iqout .eq. 3) then
         allocate (ssss(norbitals,norbitals))
         allocate (ww(norbitals))
        endif

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
!NPA
        if (iqout .eq. 3) then
         do inu = 1, norbitals
          imu = getissh(inu)
          iatom = getiatom(inu)
          lm = getlssh(inu)
          in1 = imass(iatom)
!          ww(inu) = ( Qin(imu,iatom)/(2.0d0*lm+1))*qaux(imu,in1) + 1.0d0
          ww(inu) = ( Qin(imu,iatom)/(2.0d0*lm+1))*100.0 + 1.0d0
         enddo
        endif

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
!NPA
        if (iqout .eq. 3) then

         do inu = 1, norbitals
          do imu = 1, norbitals
           ssss(inu,imu) = zzzz(inu,imu)
          end do
         end do
         do inu = 1, norbitals
          do imu = 1, norbitals
           zzzz(inu,imu) = zzzz(inu,imu)*ww(inu)*ww(imu)
          end do
         end do

        endif  ! end if (iqout .eq. 3)

! xxxx = unused (used as complex workspace in cheev call below)
! yyyy = Hamiltonian in AO basis

! write out blow
!         do inu = 1, norbitals
!          write (5001,300) ( real(yyyy(inu,imu)), imu=1,norbitals)
!          if (icluster .ne. 1) write (5002,300) (imag(yyyy(inu,imu)), imu=1,norbitals)
!         enddo



! Restore S^-1/2 from s(k)^-1/2,
        do inu = 1, norbitals
         do imu = 1, norbitals
          xxxx(inu,imu) = sm12(inu,imu,ikpoint)
         end do
        end do

! xxxx = S^-1/2 in AO basis
! zzzz = Unused (used as temporary work area below)
! yyyy = Hamiltonian in AO basis

! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        if (iqout .ne. 3) then

! Set M=H*(S^-.5)
         call zhemm ( 'R', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )

! Set Z=(S^-.5)*M
         call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, zzzz, norbitals, a0, yyyy, norbitals )

        else

! FIXME: I think these two calls we don't need them!!
         call zgemm ('C', 'N', norbitals, norbitals, norbitals, a1, xxxx,&
     &               norbitals, ssss, norbitals, a0, zzzz, norbitals)

         call zgemm ('N', 'N', norbitals, norbitals, norbitals, a1, zzzz,&
     &               norbitals, xxxx, norbitals, a0, zzzz, norbitals)
! FIXME
! Set conjg((W(WSW)^-1/2)T)*H
         call zgemm ( 'C', 'N', norbitals, norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
! Set M*(W(WSW)^-1/2)
         call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, zzzz, &
     &               norbitals, xxxx, norbitals, a0, yyyy, norbitals )

! so we have conjg((W(WSW)^-1/2)T)*H*(W(WSW)^-1/2) now

        endif

! xxxx = S^-1/2 in AO basis
! zzzz = Unused (used as complex workspace in cheev call below)
! yyyy = Hamiltonian in the MO basis set



! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
!        if (wrtout) then
!          write (*,*) '  '
!          write (*,*) ' Call diagonalizer for Hamiltonian. '
!          write (*,*) '            The energy eigenvalues: '
!          write (*,*) ' *********************************************** '
!!        end if
!

! set default size of working space
        lwork = 1
        deallocate (work)
        allocate (work(lwork))
! first find optimal working space
        call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work,      &
     &                 -1, rwork, info)
! resize working space
        lwork = work(1)
!        write (*,*) 'lwork =',lwork
        deallocate (work)
        allocate (work(lwork))
! diagonalize the overlap matrix with the new working space
        call zheev ('V', 'U', norbitals, yyyy, norbitals, eigen, work,       &
     &                lwork, rwork, info)

        if (info .ne. 0) call diag_error (info, 0)

! FIXME - Should only go up to norbitals_new, but we do not know what
! eigenvalues and eigenvectors correspond to the removed MO's.  Their
! eigenvalues will be very close to zero, but not exactly.  Also, we do not
! know if a real eigen value is near zero.

!        if (ishort .eq. 1 .and. wrtout) then
!         write (*,100) eigen(1), eigen(norbitals)
!        else if (wrtout) then
!         write (*,200) (eigen(imu), imu = 1, norbitals)
!        end if

!
! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! xxxx = S^-1/2 in AO basis
! zzzz = Unused
! yyyy = Hamiltonian eigenvectors in the MO basis
        eigen_k(1:norbitals,ikpoint) = eigen(:)
        if (iqout .ne. 2) blowre(:,:,ikpoint) = real(yyyy(:,:))
        if (iqout .ne. 2 .and. icluster .ne. 1)                              &
     &   blowim(:,:,ikpoint) = aimag(yyyy(:,:))

        if (iqout .ne. 3) then
         call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
!NPA
        else
         call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx,   &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
        end if

        bbnkre(:,:,ikpoint) = real(zzzz(:,:))
        if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(zzzz(:,:))

! write out blow
!         do inu = 1, norbitals
!          write (3001,300) (blowre(inu,imu,ikpoint), imu=1,norbitals)
!          if(icluster .ne. 1) write (3002,300) (blowim(inu,imu,ikpoint), imu=1,norbitals)
!         enddo
! write out bbnk
!         do inu = 1, norbitals
!          write (4001,300) (bbnkre(inu,imu,ikpoint), imu=1,norbitals)
!          if(icluster .ne. 1) write (4002,300) (bbnkim(inu,imu,ikpoint), imu=1,norbitals)
!         enddo

! check orthogonality
!         zzzz = Conjg (Transpose (yyyy))
!         call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, zzzz,   &
!     &               norbitals, yyyy, norbitals, a0, xxxx, norbitals )

! write out blow
!         do inu = 1, norbitals
!          write (5001,300) ( real(xxxx(inu,imu)), imu=1,norbitals)
!          if (icluster .ne. 1) write (5002,300) (imag(xxxx(inu,imu)), imu=1,norbitals)
!         enddo


! We did a symmetric orthogonalization followed by a diagnalization
! of the Hamiltonian in this "MO" basis set. This yields a net
! canonical diagnolization with matrix bbnk.

! xxxx = S^-1/2 in AO basis
! zzzz = S^-1/2 * yyyy
! yyyy = Hamiltonian eigenvectors in the MO basis

! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)
! NPA
        if (iqout .eq. 3) then
	     deallocate (ww)
	     deallocate (ssss)
        endif

        deallocate (rwork)


! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f16.8, &
     &              ' eigenvalue(norbitals) = ', f16.8)
200     format (4(2x, f12.4))
!300     format (2x, <norbitals>f12.4)
301     format (2x, 2i4,2f12.4)
302     format (2x, 4i4,3f12.4)

        return
      end subroutine  diag_Hk


