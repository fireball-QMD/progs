!	 copyright info:
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


! initpsi.f90
! Program Description
! ===========================================================================
!       This routine initialize wf coeff's and other matrices for TDSE run
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
        subroutine initpsi ()

        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use options
        use outputs
        use tdse
        use energy

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
        integer imu, inu
        integer ineigh
        integer in1, in2
        integer iorbital
        integer issh
        integer jatom
        integer jneigh
        integer mqn
        integer mbeta
        integer mmu
        integer noccupy
        integer nnu
        integer nelec_pair
        integer iexc

        real aux1, aux2
        real dot
        real gutr
        real pcharge
        real ztest
        real norm
        complex cnorm

        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec

        complex cphi
        complex fe
        complex fh
        complex a1
        complex a0
        complex phase, phasex
        complex step1, step2
        complex, dimension (:, :), allocatable :: xxxx
        complex, dimension (:, :), allocatable :: yyyy
        complex, dimension (:, :), allocatable :: zzzz

! Procedure
! ===========================================================================
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)

! allocate aux arrays
	    allocate (xxxx(norbitals,norbitals))
	    allocate (yyyy(norbitals,nelec))
        allocate (zzzz(norbitals,nelec))

        xxxx = a0
        yyyy = a0
        zzzz = a0


        aux1 = 0.0d0
        aux2 = 0.0d0
        rhoPP = 0.0d0
        rho = 0.0d0

        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        write (*,*) '  '
        write (*,*) '              Initialize TDSE run --              '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '

! setup projector from psi to phi
! here we assume the ground state configuration,
! after we can create some excitation of individual electrons
        psi2es = a0
        do ikpoint = 1, nkpoints
         do ielec = 1, nelec
          iband = (ielec-1)/2+1
!          write (*,*) 'elec, iband',ielec, iband
          psi2es(iband,ielec,ikpoint) = a1
         end do ! do ielec
        end do ! do ikpoint

        do ikpoint = 1, nkpoints
         do iexc = 1, nexcite
          ielec = idelec(iexc)
          iband = (ielec-1)/2+1
          fh = cmplx(hoccup(iexc),0.0d0)
          fe = cmplx(1.0d0-hoccup(iexc),0.0d0)
! create hole
          psi2es(iband,ielec,ikpoint) = sqrt(fh)
! add electron
          iband = eband(iexc)
          psi2es(iband,ielec,ikpoint) = sqrt(fe)
         enddo ! iexc
         do ielec = 1, nelec
          write (*,*)  '  Electron no.',ielec
          do iband = 1, norbitals
           write (*,200) iband, real(psi2es(iband,ielec,ikpoint))
          enddo ! do iband
         enddo ! do ielec
         write (*,100)
        enddo ! ikpoint


! setup psi TD-wavefunction
        psi = 0.0d0
        if (icluster .eq. 1) then
         do ikpoint = 1, nkpoints
          do ielec = 1, nelec
           do iband = 1, norbitals
            do imu = 1, norbitals
             cphi = cmplx(blowre(imu,iband,ikpoint),0.0d0)
             psi(imu,ielec,ikpoint) = psi(imu,ielec,ikpoint)                &
                  +  cphi*psi2es(iband,ielec,ikpoint)
            enddo ! do imu
           enddo ! do iband
          enddo ! do ielec
         enddo ! do ikpoint
        else
         do ikpoint = 1,nkpoints
          do ielec = 1, nelec
           do iband = 1, norbitals
            do imu = 1, norbitals
             cphi = cmplx(blowre(imu,iband,ikpoint),blowim(imu,iband,ikpoint))
             psi(imu,ielec,ikpoint) = psi(imu,ielec,ikpoint)                &
       &          +  cphi*psi2es(iband,ielec,ikpoint)
            enddo ! do imu
           enddo ! do ielec
          enddo ! do iband
         enddo ! do ikpoint
        endif

! write out GS wavefunction
       write (*,*) ' Write out GS wf-coeficients'
       do ikpoint = 1, nkpoints
        write (*,*) '  kpoint no. =',ikpoint
        do iband = 1, norbitals
         write (*,*) '   band no.', iband
         write (*,*) ' ---------------------------'
         norm = 0.0d0
         do imu = 1, norbitals
          if (icluster .eq. 1) then
           write (*,400) blowre(imu,iband,ikpoint)
           norm = norm + blowre(imu,iband,ikpoint)**2
          else
           write (*,400) blowre(imu,iband,ikpoint),blowim(imu,iband,ikpoint)
           norm = norm + blowre(imu,iband,ikpoint)**2 + blowim(imu,iband,ikpoint)**2
          endif
         enddo
         write (*,500) norm
         write (*,*) ' ---------------------------'
        enddo
       enddo

! write out initil wavefunction
       write (*,*) ' Write out initial wf-coeficients'
       do ikpoint = 1, nkpoints
        write (*,*) '  kpoint no. =',ikpoint
        do ielec = 1, nelec
         write (*,*) '   electron no.', ielec
         write (*,*) ' ---------------------------'
         cnorm = a0
         do imu = 1, norbitals
          write (*,400) psi(imu,ielec,ikpoint)
          cnorm = cnorm + (psi(imu,ielec,ikpoint))**2
         enddo
         write (*,500) cnorm
         write (*,*) ' ---------------------------'
        enddo
       enddo

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

        write (*,*) 'S^-1/2'
        do imu = 1, norbitals
         write (*,600) (real(xxxx(imu,inu)), inu=1,norbitals)
        enddo
        write (*,*)

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
        psiAO = a0
        do ielec = 1, nelec
         do inu = 1, norbitals
          psiAO(inu,ielec,ikpoint) = zzzz(inu,ielec)
         enddo ! do inu
        enddo ! do ielec
! end k-loop
       enddo ! do ikpoint

! write out psiAO & real AO eigenvectors
! write out initil wavefunction
       write (*,*) ' Write out psi AO'
       do ikpoint = 1, nkpoints
        write (*,*) '  kpoint no. =',ikpoint
        do ielec = 1, nelec
         write (*,*) '   electron no.', ielec
         write (*,*) ' ---------------------------'
         cnorm = a0
         do imu = 1, norbitals
          write (*,400) psiAO(imu,ielec,ikpoint)
          cnorm = cnorm + (psiAO(imu,ielec,ikpoint))**2
         enddo
         write (*,500) cnorm
         write (*,*) ' ---------------------------'
        enddo
       enddo
! write out original eigenvectors in AO
       write (*,*) ' Write out eigenvectors in AO'
       do ikpoint = 1, nkpoints
        write (*,*) '  kpoint no. =',ikpoint
        do iband = 1, norbitals
         write (*,*) '   band no.', iband
         write (*,*) ' ---------------------------'
         cnorm = a0
         do imu = 1, norbitals
          if (icluster .eq. 1) then
           write (*,400) bbnkre(imu,iband,ikpoint)
           cnorm = cnorm + cmplx(bbnkre(imu,iband,ikpoint),0.0d0)**2
          else
           write (*,400) bbnkre(imu,iband,ikpoint),bbnkim(imu,iband,ikpoint)
           cnorm = cnorm + cmplx(bbnkre(imu,iband,ikpoint),bbnkim(imu,iband,ikpoint))**2
          endif
         enddo
         write (*,500) cnorm
         write (*,*) ' ---------------------------'
        enddo
       enddo
! deallocate xxxx
       deallocate (xxxx)
       deallocate (yyyy)
       deallocate (zzzz)

! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
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
              rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
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
       write (*,*) 'ebs_1 =',ebs
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
        write (*,700) ebs

! Last we calculate cape term, part of HF theorem

! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
        do iatom = 1, natoms
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          in2 = imass(jatom)
          do imu = 1, num_orb(in1)
           do inu = 1, num_orb(in2)
! Finally the expressions.........
            cape(imu,inu,ineigh,iatom) = cape(imu,inu,ineigh,iatom)         &
       &           +        ebs*rho(imu,inu,ineigh,iatom)*0.5d0
           end do ! inu
          end do ! imu
         end do ! ineigh
        end do ! iatom

       if (iwrtdensity .eq. 1) then
         do iatom = 1, natoms
          write (*,*) '  '
          write (*,*) ' Write out the density matrix for iatom = ', iatom
          write (*,*) ' Number of neighbors = ', neighn(iatom)
          in1 = imass(iatom)
          write (*,*) ' Number of orbitals on iatom, num_orb(in1) = ',       &
     &     num_orb(in1)
          do ineigh = 1, neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           write (*,*) '  '
           write (*,*) ' iatom, ineigh = ', iatom, ineigh
           write (*,*) ' Number of orbitals on jatom, num_orb(in2) = ',      &
     &      num_orb(in2)
           write (*,*) ' ----------------------------------------------------- '
           do imu = 1, num_orb(in1)
            write (*,400) (rho(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           do imu = 1, num_orb(in1)
            write (*,400) (h_mat(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          end do
         end do
        end if ! iwrtdensity


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (' Band n = ', i4, ' occupy = ', f6.3)
300     format (2x, ' This is band number: ',2x, i6)
301     format (2x, i4, f10.6)
400     format (2x, 2f12.6)
500     format (3x, 'norm =',2f12.6)
600     format (2x,8f12.6)
700     format (2x, '  TD EBS = ',f16.8,' [eV]')
800     format (2x,4i3,f12.6)

        return
      end subroutine initpsi

