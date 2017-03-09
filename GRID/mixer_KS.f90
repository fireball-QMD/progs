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


! Program Description
! ===========================================================================

! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine mixer_KS (natoms, ifixcharge)

   use scf
   use grid
   use interactions
   use charges
   use neighbor_map
   use density
   use configuration , only : ratom, xl
   implicit none

! Argument Declaration and Description
! ===========================================================================
! input

! input and output
   integer, intent (in) :: natoms
   integer, intent (in) :: ifixcharge   ! fixed charges

! output

! Local Parameters and Data Declaration
! ===========================================================================
  real, parameter :: A1 = 1.0d0  ! a scale factor for occupied states of DM
  real, parameter :: A2 = 0.01d0 ! a scale factor for unoccupied states of DM

! Local Variable Declaration and Description
! ===========================================================================

   integer ipoint
   integer imix
   integer iatom
   integer jatom
   integer matom
   integer ineigh
   integer jneigh
   integer imu
   integer inu
   integer in1
   integer in2
   integer issh
   integer i

   real dqrms
   real dqmax
   real renorm
   real dens
   real zouttot
   real zcheck
! jel-MN
   integer mbeta
   real, dimension (3) :: r
   real dr

   real, dimension (numorb_max, natoms) :: QMulliken_in
   real, dimension (numorb_max, natoms) :: QMulliken_out

! Procedure
! ===========================================================================

! Check to see if input charges and output charges are within tolerance.
! If they are within tolerance, then perform a motional time step. If
! they are not within tolerence, then perform another iteration to make
! charges self-consistent.

! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
  write (*,*) 'mixer_KS'
   QMulliken_in = 0.0d0

   do iatom = 1, natoms
      in1 = imass(iatom)

! Loop over neighbors
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,ineigh)

        do imu = 1, num_orb(in1)
         do inu = 1, num_orb(in2)
          QMulliken_in(imu,iatom) = QMulliken_in(imu,iatom)                &
     &        + 0.5d0*(rho_old(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho_old(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
         end do ! do imu
        end do ! do inu
      enddo ! do ineigh
    enddo ! do iatom
! Also, only do this step if ifixcharge not equal to 1.
   if (ifixcharge .ne. 1) then

! density matrix mixing
     imix = 0
     do iatom = 1, natoms
      in1 = imass(iatom)
      matom = neigh_self(iatom)
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass (jatom)
        mbeta = neigh_b(ineigh,iatom)
        r(:) = ratom(:,iatom) - (ratom(:,jatom) + xl(:,mbeta))
        dr = sqrt (r(1)**2 + r(2)**2 + r(3)**2)
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            imix = imix + 1
! density matrix in step i (after diagonalization)
            rho_in(imix) = rho(imu,inu,ineigh,iatom)
! denisty matrix in step i-1
            rho_out(imix) = rho_old(imu,inu,ineigh,iatom)
          end do ! issh
        end do ! inu
      end do ! ineigh
     end do ! iatom

! aux vectors for Pulay mixing
     if (ialgmix .eq. 4) then
      imix = 0
      do iatom = 1, natoms
       in1 = imass(iatom)
       matom = neigh_self(iatom)
       do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass (jatom)
        mbeta = neigh_b(ineigh,iatom)
        r(:) = ratom(:,iatom) - (ratom(:,jatom) + xl(:,mbeta))
        dr = sqrt (r(1)**2 + r(2)**2 + r(3)**2)
        do inu = 1, num_orb(in2)
           imu = 1
           do issh = 1, nssh(in1)
! evaluate neutral charge per orbital
            do i = 1, (2*lssh(issh,in1)+1)
             imix = imix + 1
! get overlap as a norm
!             mwe(imix) = s_mat(imu,inu,ineigh,iatom)
! get distance as a norm
!              mwe(imix) = exp(A1*(1.0d0-dr)) + 0.03d0
             mwe(imix) = 1.0d0
!             if(Qneutral(issh,in1) .gt. 0.05)then
!               mwe(imix) = exp(A1*(1.0d0-dr)) + 0.03d0
!             else
!               mwe(imix) = A2*exp(A1*(1.0d0-dr)) + 0.01d0
!             endif
             drwe(imix) = dr
             imu = imu + 1
            enddo
!            mwe(imix) = 1.0d0
          end do ! inu
        end do ! imu
       end do ! ineigh
      end do ! iatom
     endif ! if (ialgmix)
 
     dqrms = 0
     dqmax = -99
     do imu = 1,imix
        dqmax = amax1(abs(rho_in(imu) - rho_out(imu)),dqmax)
        dqrms = dqrms + (rho_in(imu) - rho_out(imu))**2
     end do
     dqrms = sqrt(dqrms/imix)

     write (*,300) dqrms
     write (*,301) dqmax


! ===========================================================================
!                                   mixing
! ===========================================================================


! changed by honza


! call mixing procedure
     select case (ialgmix)
     case (1)
        write(*,501) 'anderson'
        call anderson (rho_in, rho_out, bmix, sigma, Kscf, idmix,    &
             &               imix , max_scf_iterations)
     case (2)
        write(*,501) 'broyden'
        call broyden (rho_in, rho_out, bmix, sigma, Kscf, idmix,    &
             &               imix , max_scf_iterations)
     case (3)
        write(*,501) 'louie'
        call louie (rho_in, rho_out, bmix, sigma, Kscf, idmix,    &
             &               imix , max_scf_iterations)
     case (4)
        write(*,501) 'pulay'
        call pulay (rho_in, rho_out, bmix, sigma, Kscf, idmix,    &
             &               imix , max_scf_iterations)
     end select !ialgmix

! end changed by honza

! restore the density matrix
     imix = 0
     do iatom = 1, natoms
      in1 = imass(iatom)
      matom = neigh_self(iatom)
      do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        mbeta = neigh_b(ineigh,iatom)
        in2 = imass (jatom)
        r(:) = ratom(:,iatom) - (ratom(:,jatom) + xl(:,mbeta))
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
            imix = imix + 1
            rho(imu,inu,ineigh,iatom) = rho_out(imix)
! add the atomic density matrix
!            if ((imu .eq. inu) .and. (ineigh .eq. matom)) then
!             rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + rhoA(imu,iatom)
!            endif
          end do ! inu
        end do ! imu
      end do ! ineigh
     end do ! iatom

     if (Kscf .gt. 1) then
        if (sigma .lt. sigmaold) then
           sigmaold = sigma
        end if
     else
        sigmaold = sigma
     end if


! ===========================================================================
!                                 end mixing
! ===========================================================================

    zouttot = 0.0d0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
       jatom = neigh_j(ineigh,iatom)
       in2 = imass(jatom)
       do imu = 1, num_orb(in1)
        do inu = 1, num_orb(in2)
         zouttot = zouttot +                                              &
  &        rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom)
        end do ! do imu
       end do ! do inu
      end do ! do ineigh
    end do ! do iatom
    renorm = ztot / zouttot
    write (*,303) renorm

    zcheck = 0.0d0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do ineigh = 1, neighn(iatom)
       jatom = neigh_j(ineigh,iatom)
       in2 = imass(jatom)
       do imu = 1, num_orb(in1)
        do inu = 1, num_orb(in2)
         rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom)*renorm
         zcheck = zcheck +                                                &
  &        rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom)
        end do ! do imu
       end do ! do inu
      end do ! do ineigh
    end do ! iatom

! write out resume
    write (*,*) ' (Before renormalization) zouttot = ', zouttot
    write (*,*) ' (After  renormalization)  zcheck = ', zcheck
    write (*,*) ' (What it must be)           ztot = ', ztot
    write (*,*) '  '
    write (*,304) sigma, sigmatol, Kscf

! Check convergence of charge; sigmatol is in scf.optional
     if (sigma .lt. sigmatol) then
      scf_achieved = .true.
! Restore true density matrix
!      rho(:,:,:,:) = rho_old(:,:,:,:)
     else
      rho_old(:,:,:,:) = rho(:,:,:,:)
     endif

! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
    QMulliken_out = 0.0d0

    do iatom = 1, natoms
       in1 = imass(iatom)

! Loop over neighbors
       do ineigh = 1, neighn(iatom)
        jatom = neigh_j(ineigh,iatom)
        in2 = imass(jatom)
        jneigh = neigh_back(iatom,ineigh)

        do imu = 1, num_orb(in1)
         do inu = 1, num_orb(in2)
          QMulliken_out(imu,iatom) = QMulliken_out(imu,iatom)              &
     &        + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
         end do
        end do

! End loop over neighbors
       end do
! End loop over atoms
    end do

! write out Mulliken Charges IN & OUT of mixerG
!    write (*,*) '  Mulliken Charges BEFORE mixing:'
!    do iatom = 1, natoms
!      in1 = imass(iatom)
!      write (*,401) (QMulliken_in(imu,iatom),imu=1,num_orb(in1))
!    enddo ! iatom
!    write (*,*) '  ====================================='
!    write (*,*) '  Mulliken Charges AFTER mixing:'
!    do iatom = 1, natoms
!      in1 = imass(iatom)
!      write (*,402) (QMulliken_out(imu,iatom),imu=1,num_orb(in1))
!    enddo ! iatom
!    write (*,*) '  ====================================='

! Skip to here if charges are fixed
   end if ! end if (ifixcharge)



! Format Statements
! ===========================================================================
300     format (2x, ' Deviation (rms) of input/output charges = ', f9.4)
301     format (2x, ' Deviation (max) of input/output charges = ', f9.4)
303     format (2x, ' Renormalization of density:   renorm = ', 4x, f14.8)
304     format (' =======> sigma = ', e14.7,                              &
     &          ' Must be less than', e14.7, ' SCF step = ', i3)
401     format (2x,'QmIN ',18(2x,f8.4))
402     format (2x,'QmOUT',18(2x,f8.4))
501     format (2x,'Your chosen mixing routine is: ',a)
600     format (4i4, f8.4)

      return
    end subroutine mixer_KS
