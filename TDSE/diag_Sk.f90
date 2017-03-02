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

! diag_Sk.f90
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
        subroutine diag_Sk (iqout, icluster, iwrteigen, ikpoint, sks, nkpoints)

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
        real*8 sqlami
        real*8, dimension (norbitals) :: slam
        real, dimension (3) :: vec

        complex*16 a0
        complex*16 a1
        complex*16 phase

! A bunch of memory to be used in many ways
        complex*16, dimension (:, :), allocatable :: xxxx
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

        if (wrtout) then
          write (*,*) '  '
          write (*,*) ' ****************************************************** '
          write (*,*) '  '
          write (*,*) '         Diagonalize S(k) -- ikpoint = ', ikpoint
          write (*,*) '  '
          write (*,*) ' ****************************************************** '
        end if

        allocate (xxxx(norbitals,norbitals))
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
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + phase*s_mat(imu,inu,ineigh,iatom)
           end do ! do inu
          end do ! do imu
         end do ! do ineigh
        end do ! do iatom

! xxxx = unused (used as complex workspace in cheev call below)
! zzzz = Overlap in AO basis


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

        if (ishort .eq. 1 .and. wrtout) then
         write (*,100) slam(1), slam(norbitals)
        else if (wrtout) then
         write (*,200) (slam(imu), imu = 1, norbitals)
        end if

! xxxx = unused
! zzzz = Overlap eigenvectors



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
           end do
          end do
         end if
         write (*,*) ' '

         do imu = mineig, norbitals
          jmu = imu - mineig + 1
          zzzz(:,jmu) = zzzz(:,imu)
          slam(jmu) = slam(imu)
         end do
        end if ! if (norbitals_new .ne. norbitals)

!
! CALCULATE (S^-1/2) --> sm12
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
     &               norbitals, zzzz, norbitals, a0, xxxx, norbitals)

!NPA  now we do X=W(WSW)^-1/2, before X=S^-1/2
        if (iqout .eq. 3) then
         do imu=1, norbitals
          xxxx(imu,:)=xxxx(imu,:)*ww(imu)
         end do
        endif

! Now put S^-1/2 into s(k)^-1/2, this will be remembered for the duration of
! the ionic time step
        do inu = 1, norbitals
         do imu = 1, norbitals
          sm12(imu,inu,ikpoint) = xxxx(imu,inu)
         end do ! do imu
        end do ! do inu

!        write (*,*) 'S^-1/2'
!        do imu = 1, norbitals
!         write (*,600) (real(sm12(imu,inu,ikpoint)), inu=1,norbitals)
!        enddo
!       write (*,*)

! ****************************************************************************
! CALCULATE (S^1/2) --> s12
! ****************************************************************************
!        do imu = 1, norbitals_new
!         sqlami = slam(imu)**(0.25d0)
!         zzzz(:,imu) = zzzz(:,imu)*sqlami
!        end do
!        call zgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz,&
!    &               norbitals, zzzz, norbitals, a0, xxxx, norbitals)

!NPA  now we do X=W(WSW)^-1/2, before X=S^-1/2
!        if (iqout .eq. 3) then
!         do imu=1, norbitals
!         xxxx(imu,:)=xxxx(imu,:)*ww(imu)
!        end do
!        endif

! Now put S^1/2 into s(k)^1/2, this will be remembered for the duration of
! the ionic time step
!        do inu = 1, norbitals
!         do imu = 1, norbitals
!          s12(imu,inu,ikpoint) = xxxx(imu,inu)
!        end do ! do imu
!        end do ! do inu

! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
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
600     format (2x, 8f12.6)

        return
      end subroutine diag_Sk
