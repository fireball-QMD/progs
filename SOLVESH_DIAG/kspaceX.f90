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
!
! Code rewritten by:

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine kspace  (natoms, nprocs, my_proc, Kscf, iqout, icluster,  &
     &                     iwrteigen, ikpoint, sks, norbitals_new,          &
     &                     nkpoints, iwrtdos, iwrthop, iwrtatom)
        use configuration  
        use density
        use dimensions
        use interactions
        use neighbor_map
        use charges
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrteigen
        integer, intent (in) :: Kscf
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
        integer, intent (in) :: nkpoints
! CGP
        integer, intent (in) :: iwrtdos
        integer, intent (in) :: iwrthop
        integer, intent (in) :: iwrtatom
! end CGP

        real, intent (in), dimension (3) :: sks
! Output
        integer, intent (out) :: norbitals_new


! Local Parameters and Data Declaration
! ===========================================================================
        real*16, parameter :: overtol = 1.0d-4

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

        real dot
        real sqlami

        real, dimension (norbitals) :: eigen

        real, dimension (norbitals) :: slam
        real, dimension (3) :: vec

        complex a0
        complex a1
        complex phase

! A bunch of memory to be used in many ways
        complex, dimension (:, :), allocatable :: xxxx
        complex, dimension (:, :), allocatable :: yyyy
        complex, dimension (:, :), allocatable :: zzzz
        complex, dimension (:, :), allocatable :: ssss 
        complex, dimension (:, :, :), allocatable, save :: sm12_save
!NPA
        real, dimension (:), allocatable :: ww

! work vector for cheev/cheevd
        complex, allocatable, dimension (:) :: work
        real, allocatable, dimension (:) :: rwork
        integer, allocatable, dimension (:) :: iwork
        integer lwork, lrwork, liwork
        logical, parameter :: divide = .false.  ! do we divide and conquer?

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)

! CGP
        if(iwrtdos.ge.1 .or. iwrthop.ge.1 .or. iwrtatom.ge.1) hamk = a0
!end CGP

        ishort = 1
        if (iwrteigen .eq. 1) ishort = 0

!        if (wrtout) then
!          write (*,*) '  '
!          write (*,*) ' ****************************************************** '
!          write (*,*) '  '
!          write (*,*) '         Welcome to kspaceG -- ikpoint = ', ikpoint
!          write (*,*) '  '
!          write (*,*) ' ****************************************************** '
!        end if

        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))
        allocate (ssss(norbitals,norbitals))

        if (divide) then
          lwork = 100*norbitals + norbitals*norbitals
          lrwork = 100*norbitals + 3*norbitals*norbitals ! for old versions of cheevd
          liwork = 10*norbitals
          allocate (work(lwork))
          allocate (iwork(liwork))
          allocate (rwork(lrwork))
        else
          lwork = norbitals*norbitals ! Use xxxx, yyyy and zzzz for work area
          lrwork = 3*norbitals
          allocate (rwork(lrwork))
        end if

        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(norbitals,norbitals,nkpoints))
        end if

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
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
            yyyy(jmu,jnu) = yyyy(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
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

! save the Overlap matrix for S^1/2
        do inu = 1, norbitals
         do imu = 1, norbitals
          ssss(inu,imu) = zzzz(inu,imu)
         enddo
        enddo

!        write (*,*) 'H hamiltonian'
!        do inu = 1, norbitals
!          write (*,300) (real(yyyy(inu,imu)), imu=1,norbitals) 
!        enddo

!        write (*,*) ''
!        write (*,*) 'S hamiltonian'
!        do inu = 1, norbitals
!          write (*,300) (real(zzzz(inu,imu)), imu=1,norbitals) 
!        enddo

        norbitals_new = norbitals

! xxxx = unused (used as complex workspace in cheev call below)
! zzzz = Overlap in AO basis
! yyyy = Hamiltonian in AO basis

! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
! 
! ZHEGV  -  compute all the eigenvalues, and optionally, the eigenvectors
!   A*x=(lambda)*B*x

!        if (wrtout) then
!          write (*,*) '  '
!          write (*,*) ' Call diagonalizer for Hamiltonian. '
!          write (*,*) '            The energy eigenvalues: '
!          write (*,*) ' *********************************************** '
!        end if

! Eigenvectors are needed to calculate the charges and for forces!
! Compute the eigenvalues and the eigenvectors of a  complex  generalized  
! Hermitian-definite  eigenproblem,  of the form A*x=(lambda)*B*x,
        call zhegv (1, 'V', 'U', norbitals, yyyy, norbitals, zzzz,         &
     &               norbitals, eigen, xxxx, lwork, rwork , info)

        if (info .ne. 0) call diag_error (info, 0)


        if (ishort .eq. 1 .and. wrtout) then
         write (*,100) eigen(1), eigen(norbitals)
        else if (wrtout) then
         write (*,200) (eigen(imu), imu = 1, norbitals)
        end if

!
! INFORMATION FOR THE DENSITY MATRIX
! ****************************************************************************
! xxxx = unused
! zzzz = triangular factor U from the Cholesky factorization B =  U**H*U
! yyyy = Hamiltonian eigenvectors in the AO basis

        eigen_k(1:norbitals,ikpoint) = eigen(:)
        bbnkre(:,:,ikpoint) = real(yyyy(:,:))
        if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(yyyy(:,:))

! *****************************************************************
!         Peform S^1/2 transformation for LOWDIN CHARGES
! *****************************************************************
        if (iqout .ne. 2) then

! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
! NPA
         if (Kscf .eq. 1 .or. iqout .eq. 3) then

! copy the Overlap matrix
          do inu = 1, norbitals
           do imu = 1, norbitals
            zzzz(inu,imu) = ssss(inu,imu)
           enddo
          enddo

!NPA
          if (iqout .eq. 3) then
           allocate (ww(norbitals))
           do inu = 1, norbitals
            imu = getissh(inu)
            iatom = getiatom(inu)
            lm = getlssh(inu)
            in1 = imass(iatom)
            ww(inu) = ( Qin(imu,iatom)/(2.0d0*lm+1))*qaux(imu,in1) + 1.0d0
           enddo
           do inu = 1, norbitals
            do imu = 1, norbitals
             zzzz(inu,imu) = zzzz(inu,imu)*ww(inu)*ww(imu)
            end do
           end do
          endif ! (iqout .eq. 3)

! Call the diagonalizer
          if (wrtout) then
           write (*,*) ' Call diagonalizer for overlap. '
           write (*,*) '                  The overlap eigenvalues: '
           write (*,*) ' ******************************************************* '
          end if

          if (divide) then
           call zheevd('V', 'U', norbitals, zzzz, norbitals, slam, work,      &
     &               lwork, rwork , lrwork, iwork, liwork, info )
          else
           call zheev ('V', 'U', norbitals, zzzz, norbitals, slam, xxxx,       &
     &               lwork, rwork , info)
          end if

          if (info .ne. 0) call diag_error (info, 0)

          if (ishort .eq. 1 .and. wrtout) then
           write (*,100) slam(1), slam(norbitals)
          else if (wrtout) then
           write (*,200) (slam(imu), imu = 1, norbitals)
          end if


! xxxx = unused
! zzzz = Overlap eigenvectors
! yyyy = Hamiltonian eigenvectors in the AO basis


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
           if(ishort .eq. 1) then    ! Don't print out again if done above
            write (*,*) '            The overlap eigenvalues: '
            write (*,*) ' ********************************************** '
            write (*,200) (slam(imu), imu = 1, norbitals)
           else                      ! They asked for extra printout
            write(*,*) ' '
            write(*,*) ' Eigenvectors that correspond to eigenvalues'
            write(*,*) ' that were eliminated.  These might provide'
            write(*,*) ' insight into what atomic orbitals are causing'
            write(*,*) ' the problem.'
            write(*,*) ' '
            do imu = 1, mineig - 1
             write(*,*) ' eigenvector',imu
             do jmu = 1, norbitals
              write(*,*) jmu,' ',zzzz(jmu,imu)
             end do ! do jmu
            end do ! do imu
           end if
           write (*,*) ' '


           do imu = mineig, norbitals
            jmu = imu - mineig + 1
            zzzz(:,jmu) = zzzz(:,imu)
            slam(jmu) = slam(imu)
           end do
          end if

!
! CALCULATE (S^1/2) --> sm1
! ****************************************************************************
! In a diagonal reperesentation (Udagger*S*U = s, s is a diagonal matrix)
! We just take the inverse of the square roots of the eigenvalues to get
! s^1/2. Then we 'undiagonalize' the s^1/2 matrix back to get
! S^1/2 = U*s^1/2*Udagger.
! Note: We do S^1/4 here, because the sqlami contribution get squared
! after it is combined with overlap.
          do imu = 1, norbitals_new
           sqlami = slam(imu)**(0.25d0)
           zzzz(:,imu) = zzzz(:,imu)*sqlami
          end do

          call zgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz,&
     &               norbitals, zzzz, norbitals, a0, xxxx, norbitals)

!NPA  now we do X=W(WSW)^1/2, before X=S^1/2
          if (iqout .eq. 3) then
           do imu=1, norbitals
            xxxx(imu,:)=xxxx(imu,:)*ww(imu)
           end do
          endif

! Now put S^1/2 into s(k)^1/2, this will be remembered for the duration of
! the scf cycle.
          do inu = 1, norbitals
           do imu = 1, norbitals
            sm12_save(imu,inu,ikpoint) = xxxx(imu,inu)
           end do
          end do

         else ! (if Kscf .eq. 1 .and iqout .ne. 3)

! Now if not first iteration
! Restore S^1/2 from s(k)^1/2,
          xxxx(:,:) = sm12_save(:,:,ikpoint)
         end if

! xxxx = S^1/2 in AO basis
! zzzz = Unused
! yyyy = Hamiltonian eigenvectors in the AO basis

         if (iqout .ne. 3) then
          call zhemm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
!NPA
         else
          call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx,   &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
         end if

! xxxx = S^1/2 in AO basis
! zzzz = Hamiltonian eigenvectors in the MO basis
! yyyy = Hamiltonian eigenvectors in the AO basis

         if (iqout .eq. 1) blowre(:,:,ikpoint) = real(zzzz(:,:))
         if (iqout .eq. 1 .and. icluster .ne. 1)                              &
     &    blowim(:,:,ikpoint) = aimag(zzzz(:,:))

! write out blow 
         do inu = 1, norbitals
          write (301,300) (blowre(inu,imu,ikpoint), imu=1,norbitals)
          if (icluster .ne. 1) write (302,300) (blowim(inu,imu,ikpoint), imu=1,norbitals)
         enddo
! write out bbnk
         do inu = 1, norbitals
          write (401,300) (bbnkre(inu,imu,ikpoint), imu=1,norbitals)
          if (icluster .ne. 1) write (402,300) (bbnkim(inu,imu,ikpoint), imu=1,norbitals)
         enddo

! check orthogonality
         yyyy(:,:) = Conjg (Transpose (zzzz(:,:)))
         call zgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, yyyy,   &
     &               norbitals, zzzz, norbitals, a0, xxxx, norbitals )

! write out blow 
         do inu = 1, norbitals
          write (501,300) ( real(xxxx(inu,imu)), imu=1,norbitals)
          if (icluster .ne. 1) write (502,300) (imag(xxxx(inu,imu)), imu=1,norbitals)
         enddo

        endif ! if (iqout .ne. 2)

! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)
        deallocate (ssss)
        if (iqout .eq. 3) deallocate (ww)

        deallocate (rwork)
        if (divide) then
          deallocate (work)
          deallocate (iwork)
        end if

! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f10.6, &
     &              ' eigenvalue(norbitals) = ', f10.6)
200     format (4(2x, f12.4))
300     format (<norbitals>f12.4)

        return
      end subroutine kspace
!
      subroutine kspace_slave (nprocs, my_proc)
        implicit none
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
        stop
      end subroutine kspace_slave

