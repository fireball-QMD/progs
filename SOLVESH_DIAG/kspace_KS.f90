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

! kspace_KS.f90
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
        subroutine kspace_KS (nprocs, my_proc, Kscf, iqout, icluster,  &
     &                     iwrteigen, ikpoint, sks, igap)

        use configuration
        use density
        use dimensions
        use interactions
        use neighbor_map
        use hartree_fock

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrteigen
        integer, intent (in) :: Kscf
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
        integer, intent (in) :: igap

        real, intent (in), dimension (3) :: sks
! Output


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

        real dot
        real*8 sqlami

        real*8, dimension (norbitals) :: eigen

        real*8, dimension (norbitals) :: slam
        real, dimension (3) :: vec

        complex*16 a0
        complex*16 a1
        complex*16 phase

! A bunch of memory to be used in many ways
        complex*16, dimension (:, :), allocatable :: xxxx
        complex*16, dimension (:, :), allocatable :: yyyy
        complex*16, dimension (:, :), allocatable :: zzzz
        complex*16, dimension (:, :, :), allocatable, save :: sm12_save

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
          write (*,*) '         Welcome to kspaceG -- ikpoint = ', ikpoint
          write (*,*) '  '
          write (*,*) ' ****************************************************** '
        end if

        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))

        lwork = norbitals*norbitals ! Use xxxx, yyyy and zzzz for work area
        lrwork = 3*norbitals
        allocate (rwork(lrwork))
        allocate (work(lwork))

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


! optionally, save Hamiltonian & Overlap for DOS
        xxxx = yyyy


        norbitals_new = norbitals

! xxxx = unused (used as complex workspace in cheev call below)
! zzzz = Overlap in AO basis
! yyyy = Hamiltonian in AO basis

! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
!
! ZHEGV  -  compute all the eigenvalues, and optionally, the eigenvectors
!   A*x=(lambda)*B*x

        if (wrtout) then
          write (*,*) '  '
          write (*,*) ' Call diagonalizer for Hamiltonian. '
          write (*,*) '            The energy eigenvalues: '
          write (*,*) ' *********************************************** '
        end if

! Eigenvectors are needed to calculate the charges and for forces!
! Compute the eigenvalues and the eigenvectors of a  complex  generalized
! Hermitian-definite  eigenproblem,  of the form A*x=(lambda)*B*x,
        call zhegv (1, 'V', 'U', norbitals, yyyy, norbitals, zzzz,         &
     &               norbitals, eigen, work, lwork, rwork , info)

        if (info .ne. 0) call diag_error (info, 0)


        if (ishort .eq. 1 .and. wrtout) then
         write (*,100) eigen(1), eigen(norbitals)
        else if (wrtout) then
         write (*,200) (eigen(imu), imu = 1, norbitals)
        end if

!
! INFORMATION FOR THE DENSITY MATRIX
! ****************************************************************************
! xxxx =
! zzzz = the triangular factor U or L from the Cholesky factorization
! yyyy = matrix of eigenvectors

        eigen_k(1:norbitals,ikpoint) = eigen(:)
        bbnkre(:,:,ikpoint) = real(yyyy(:,:))


        if (icluster .ne. 1) bbnkim(:,:,ikpoint) = aimag(yyyy(:,:))

! Deallocate Arrays
! ===========================================================================
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

        deallocate (rwork)
        deallocate (work)


! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f10.6, &
     &              ' eigenvalue(norbitals) = ', f10.6)
200     format (4(2x, f12.4))
!300     format (<norbitals>f12.4)

        return
      end subroutine kspace_KS
!
      subroutine kspace_KS_slave (nprocs, my_proc)
        implicit none
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
        stop
      end subroutine kspace_KS_slave

