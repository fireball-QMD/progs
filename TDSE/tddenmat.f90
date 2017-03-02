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

! tddenmat.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy for TDSE.
!
! ===========================================================================
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
        subroutine tddenmat ()

        use charges
        use configuration
        use constants_fireball
        use options
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use energy
        use tdse

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iband
        integer ielec
        integer ikpoint
        integer imu
        integer inu
        integer ineigh
        integer in1
        integer in2
        integer issh
        integer jatom
        integer jneigh
        integer mbeta
        integer mmu
        integer nnu

        real dot
        real gutr
        real pcharge
        real rnorm
        real cnorm
        real, dimension (3) :: vec

        complex a0
        complex a1
        complex phase, phasex
        complex step1
        complex step2
        complex, dimension (:, :), allocatable :: xxxx
        complex, dimension (:, :), allocatable :: yyyy
        complex, dimension (:, :), allocatable :: zzzz

        real, dimension (nelec, nkpoints) :: pek

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

        rhoPP = 0.0d0
        rho = 0.0d0
        cape = 0.0d0

        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        write (*,*) '  '
        write (*,*) '                   Welcome to TD-denmat --              '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '

! allocate aux arrays
	allocate (xxxx(norbitals,norbitals))
	allocate (yyyy(norbitals,nelec))
        allocate (zzzz(norbitals,nelec))

        xxxx = a0
        yyyy = a0
        zzzz = a0

!***************************************************
!
! Project MO basis set to AO basis set
!
! **************************************************
        do ikpoint = 1,nkpoints

! restore s-1/2 matrix
         do imu = 1, norbitals
          do inu = 1, norbitals
           xxxx(inu,imu) = sm12(inu,imu,ikpoint)
          enddo
         enddo
! restore phi
         do ielec = 1, nelec
          do inu = 1, norbitals
           yyyy(inu,ielec) = psi(inu,ielec,ikpoint)
          enddo
         enddo

! xxxx = S^-1/2 in AO basis
! yyyy = eigenvectors in the MO basis
        if (iqout .ne. 3) then
         call zhemm ( 'L', 'U', norbitals, nelec, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
!NPA
        else
         call zgemm ( 'N', 'N', norbitals, nelec, norbitals, a1, xxxx,   &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
        end if

! save psi in AO basis set
         do ielec = 1, nelec
          do inu = 1, norbitals
           psiAO(inu,ielec,ikpoint) = zzzz(inu,ielec)
          enddo ! do inu
         enddo ! do imu
! end k-loop
        enddo ! do ikpoint


! deallocate xxxx
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)


! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
        pek = 0.0d0
! Loop over electrons
        do ielec = 1, nelec
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
         do iatom = 1, natoms
          in1 = imass(iatom)
          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

! Loop over the special k points.
           do ikpoint = 1, nkpoints
            dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
            phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
            do imu = 1, num_orb(in1)
             mmu = imu + degelec(iatom)
             step1 = phase*psiAO(mmu,ielec,ikpoint)
             do inu = 1, num_orb(in2)
              nnu = inu + degelec(jatom)
              step2 = step1*conjg(psiAO(nnu,ielec,ikpoint))
              gutr = real(step2)
! Finally the expressions.........
              rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
              pek(ielec,ikpoint) = pek(ielec,ikpoint)                          &
      &         + h_mat(imu,inu,ineigh,iatom)*gutr
             end do ! inu
            end do ! imu
           end do ! ikpoint
          end do ! ineigh
         end do ! iatom
        end do ! ielec

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   M A T R I X     D E N S I T Y
!                          PP-neighbors
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! Loop over paired electrons
        do ielec = 1, nelec
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
         do iatom = 1, natoms
          in1 = imass(iatom)
          do ineigh = 1, neighPPn(iatom)
           mbeta = neighPP_b(ineigh,iatom)
           jatom = neighPP_j(ineigh,iatom)
           in2 = imass(jatom)
           vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

! Loop over the special k points.
           do ikpoint = 1, nkpoints
            dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
            phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
            do imu = 1, num_orb(in1)
             mmu = imu + degelec(iatom)
             step1 = phase*psiAO(mmu,ielec,ikpoint)
             do inu = 1, num_orb(in2)
              nnu = inu + degelec(jatom)
              step2 = step1*conjg(psiAO(nnu,ielec,ikpoint))
              gutr = real(step2)
! Finally the expressions.........
              rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom) + gutr
              pek(ielec,ikpoint) = pek(ielec,ikpoint)                         &
     &           + vnl(imu,inu,ineigh,iatom)*gutr
             end do ! inu
            end do ! imu
           end do ! ikpoint
          end do ! ineigh
         end do ! iatom
        end do ! ielec

! Now calculate Ebs term using <c|H|c>, we need that to calculate forces
        ebs = 0.0d0
! 1. we sum hamiltonian with generic neighbors
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
        do iatom = 1, natoms
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          do imu = 1, num_orb(in1)
           do inu = 1, num_orb(in2)
! Finally the expressions.........
            ebs =  ebs + rho(imu,inu,ineigh,iatom)*h_mat(imu,inu,ineigh,iatom)
           end do ! inu
          end do ! imu
         end do ! ineigh
        end do ! iatom

! 2. we add Vnl part

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
        do iatom = 1, natoms
         in1 = imass(iatom)
         do ineigh = 1, neighPPn(iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)
          do imu = 1, num_orb(in1)
           do inu = 1, num_orb(in2)
! Finally the expressions.........
            ebs =  ebs + rhoPP(imu,inu,ineigh,iatom)*vnl(imu,inu,ineigh,iatom)
           end do ! inu
          end do ! imu
         end do ! ineigh
        end do ! iatom
        write (*,100) ebs

! Last we calculate cape term, part of H-F theorem
! ****************************************************************************
!
!                      C O M P U T E    F O R C E S
! ****************************************************************************

! Loop over electrons
        do ielec = 1, nelec
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
         do iatom = 1, natoms
          in1 = imass(iatom)
          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

! Loop over the special k points.
           do ikpoint = 1, nkpoints
            dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
            phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
            do imu = 1, num_orb(in1)
             mmu = imu + degelec(iatom)
             step1 = phase*psiAO(mmu,ielec,ikpoint)
             do inu = 1, num_orb(in2)
              nnu = inu + degelec(jatom)
              step2 = step1*conjg(psiAO(nnu,ielec,ikpoint))
              gutr = real(step2)
! Finally the expressions.........
              cape(imu,inu,ineigh,iatom) =                                 &
     &           cape(imu,inu,ineigh,iatom) + pek(ielec,ikpoint)*gutr
             end do ! inu
            end do ! imu
           end do ! ikpoint
          end do ! ineigh
         end do ! iatom
        end do ! ielec

! writeout the density matrix into file
! NOTE: we need the density matrix to restart TDSE
        open (unit = 20, file = 'denmat.dat', status = 'unknown')
! loop over atoms
        do iatom = 1, natoms
         in1 = imass(iatom)
         write (20,1100) iatom,neighn(iatom),num_orb(in1)
! loop over neighbors
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          write (20,1101) ineigh, jatom, num_orb(in2)
          do imu = 1, num_orb(in1)
           do inu = 1, num_orb(in2)
            write (20,1102) rho(imu,inu,ineigh,iatom)
           end do ! do inu
          end do ! do imu
         end do ! do ineigh
        end do ! do iatom
! close file
        close (20)

! ****************************************************************************
!
! W R I T E    O U T    C H A R G E S    F I L E
! ****************************************************************************
! Open the file CHARGES which contain the Lowdin charges for restart purposes.
        if (ifixcharge .eq. 0 .and. wrtout) then
         open (unit = 21, file = 'CHARGES', status = 'unknown')
         write (21,1000) natoms, basisfile, iqout

! Write the output charges to the CHARGES file.
         do iatom = 1, natoms
          in1 = imass(iatom)
          write (21,1001) (Qin(issh,iatom), issh = 1, nssh(in1))
         end do
         close (unit = 21)
        end if

! Format Statements
! ===========================================================================
100     format (2x, '  TD EBS = ',f16.8,' [eV]')
400     format (2x, 2f12.6)
500     format (3x, 'norm =',3f12.6)
600     format (2x,8f12.6)
800     format (2x,4i3,f12.6)
!900     format (<norbitals>f12.6)
1000    format (2x, i5, 2x, a40, 2x, i2)
1001    format (2x, 18f14.8)
1100    format (2x, 3i5)
1101    format (3x, 3i5)
1102    format (f16.8)

        return
      end subroutine tddenmat

