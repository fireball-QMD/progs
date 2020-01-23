! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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

 
! denmat.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! (modified by P. Jelinek; May 2005 Utah)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine compute_all_charges  (ifixcharge, icluster)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ifixcharge
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iband
        integer jband
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
 
        integer, dimension (norbitals) :: ioccupy
        integer, dimension (norbitals, nkpoints) :: ioccupy_k
        
        real y
        real aux1, aux2, aux3
        real deltae
        real dot
        real gutr
        real pcharge
        real ztest
        real norm 
        real, dimension (norbitals, nkpoints) :: foccupy
        real, dimension (3) :: r1, r2, r21
        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec
 
        !complex ai
        !complex phase, phasex
        !complex step1, step2
 
        logical read_occupy

! A bunch of memory to be used in many ways
        integer jnu,jmu
        !complex*16, dimension (:, :), allocatable :: xxxx
        !complex*16, dimension (:, :), allocatable :: yyyy
        !complex*16 a0
        !complex*16 a1

! Procedure
! ===========================================================================

 ! allocate (eigen_k (norbitals, nkpoints))

 ! call fermie (norbitals, ztot, eigen_k, efermi, ioccupy_k, foccupy)

 ! deallocate (eigen_k)
! ****************************************************************************
!
!  C O M P U T E    L O W D I N    C H A R G E S
! ****************************************************************************
! Initialize
! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
        !if (iqout .eq. 2) then
         open(unit = 11, file = 'Mulliken_Charges', status = 'unknown')
         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom) 
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(imu,iatom) 
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)
 
! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            jneigh = neigh_back(iatom,ineigh)
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                    &
     &        + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
             end do
            end do
! End loop over neighbors
           end do
! Finally the imu loop.
           imu = 0
           do imu = 1, num_orb(in1)
             Qout(imu,iatom) = Qout(imu,iatom) + QMulliken(imu,iatom)
             QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(imu,iatom) 
           end do
           write(11,*) iatom, (Qout(imu,iatom), imu=1, num_orb(in1))
! End loop over atoms
          end do
         end if     ! endif of ifixcharges
         close(11)
        !end if      ! endif of iqout .eq. 2
 


! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N - D I P O L E    C H A R G E S
! ****************************************************************************
! Compute Mulliken-dipole charges.
        !if (iqout .eq. 4) then
         open(unit = 22, file = 'MullDip_Charges', status = 'unknown')
         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom)
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qin(imu,iatom)
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)
           r1(:) = ratom(:,iatom)
! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            r2(:) = ratom(:,jatom)


! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the
! bondcharge, and the bc distance, y
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

            jneigh = neigh_back(iatom,ineigh)
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                   &
     &        +0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
             end do
            end do
             
! dipole correction. Only if the two atoms are different
          if (y .gt. 1.0d-05) then
           

            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)

              QMulliken(imu,iatom) = QMulliken(imu,iatom)+                  &
     &        (-rho(imu,inu,ineigh,iatom)*dip(imu,inu,ineigh,iatom)         &
     &        + rho(inu,imu,jneigh,jatom)*dip(inu,imu,jneigh,jatom))/y
             
 
           !         write(*,*) 'DIPOLE when iatom,jatom,imu,inu = ',  &
     ! &        iatom,jatom,imu,inu,' is',dip(imu,inu,ineigh,iatom),    &
     ! &        dip(inu,imu,jneigh,jatom) 

     !          write(*,*) 'Qa(',imu,',',iatom,')= ',QMulliken(imu,iatom)

             end do
            end do
          end if !end if y .gt. 1.0d-05)

         

! End loop over neighbors
           end do

! Finally the imu loop.
           imu = 0
           do imu = 1, num_orb(in1)               
               Qout(imu,iatom) = Qout(imu,iatom)+QMulliken(imu,iatom)
               QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qout(imu,iatom)
           end do

!Check whether there are negative charges and correct
!If there's more than one shell whose charge is negative, more work is
!needed, but that'd be quite pathological a situation...
            do imu = 1, num_orb(in1)

               if( Qout(imu,iatom) .lt. 0 .and. num_orb(in1) .gt. 1 ) then
           
                  do inu = 1,num_orb(in1)

                     if ( inu .ne. imu ) then

                        Qout(inu,iatom) = Qout(inu,iatom)+            &
                 &                         Qout(imu,iatom)/(num_orb(in1)-1)

                     end if !end if inu .ne. imu

                  end do !end if inu = 1,num_orb(in1)

                  Qout(imu,iatom) = 0.0d0               

               end if !end if  Qout(imu,iatom) .lt. 0
            end do !end do imu = 1, num_orb(in1)

            write(22,*) iatom, (Qout(imu,iatom), imu = 1, num_orb(in1))
! End loop over atoms
          end do
         end if     ! endif of ifixcharges
         close(22)
        !end if      ! endif of iqout .eq. 4




! Format Statements
! ===========================================================================
100     format (2x, 2i4, f8.4)
200     format (' Band n = ', i4, ' k-points: ioccupy = ', i2)
201     format (' Band n = ', i4, ' foccupy = ', f6.3)
300     format (2x, ' This is band number: ',2x, i6)
301     format (2x, i4, f10.6)
110     format (2x, 4f10.6)
!120     format (2x, <norbitals>f10.6) 
400     format (2x, 'Qmull =',10f10.6)

        return
      end subroutine compute_all_charges
        
