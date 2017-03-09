! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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
! University of Regensburg - Juergen Fritsch

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


! num_neigh_tot.f90
! Program Description
! ===========================================================================
! The subroutine gets total number of neighbors (normal+PP) of each atom
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! ==========================================================================
!
! Program Declaration
! ===========================================================================
         subroutine num_neigh_tot (numorb_max)

         use neighbor_map
         use module_dos
         use configuration

         implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
         integer, intent(in)          :: numorb_max

! Local Parameters and Data Declaration
! ============================================================================

! Local Variable Declaration and Description
! ===========================================================================
         integer,allocatable   :: neighj_aux(:,:)
         integer,allocatable   :: neighb_aux(:,:)
         integer               :: iatom
         integer               :: jatom
         integer               :: jatomPP
         integer               :: num_neigh
         integer               :: num_neighPP
         integer               :: num_neig_tot
         integer               :: ineighPP
         integer               :: ineigh
         integer               :: count_neig
         integer               :: mbeta
         integer               :: mbetaPP

! ============================================================================
! Initializate

! test allocated arrays
         if ( allocated (neighj_tot)) then
         deallocate (neighj_tot)
         deallocate (neighb_tot)
         endif
         if ( allocated (hr_box)) deallocate (hr_box)
! JOM-warning : these allocations seem arbitrary. we should
! improve
         allocate (neighj_aux(neigh_max+neighPP_max**2,natoms))
         allocate (neighb_aux(neigh_max+neighPP_max**2,natoms))
         neighj_aux = 0
         neighb_aux = 0

! Loop over atoms
         do iatom = 1,natoms
           num_neigh = neighn(iatom)
           num_neighPP = neighPPn(iatom)
           num_neig_tot = num_neigh 
           neighj_aux(1:num_neigh,iatom) = neigh_j(1:num_neigh,iatom)
           neighb_aux(1:num_neigh,iatom) = neigh_b(1:num_neigh,iatom)

! Loop over PP neighbors
           do ineighPP = 1,num_neighPP
             count_neig = 0
             jatomPP = neighPP_j(ineighPP,iatom)
             mbetaPP = neighPP_b(ineighPP,iatom)
! Loop over neighbors
             do ineigh = 1, num_neigh
               jatom = neigh_j(ineigh,iatom)
               mbeta = neigh_b(ineigh,iatom)
               if ((jatomPP .eq. jatom .and. mbetaPP .eq. mbeta)) then
               count_neig = 1
               end if
             end do ! ineigh 
             if (count_neig .eq. 0) then
             num_neig_tot = num_neig_tot + 1
             neighj_aux(num_neig_tot,iatom) = jatomPP 
             neighb_aux(num_neig_tot,iatom) = mbetaPP
             end if
           end do ! ineighPP
         neighn_tot(iatom) = num_neig_tot
         end do ! iatom

         num_neig_maxtot = maxval(neighn_tot(1:natoms))
         allocate (neighj_tot(num_neig_maxtot,natoms))
         allocate (neighb_tot(num_neig_maxtot,natoms))
         allocate (hr_box(numorb_max,numorb_max,natoms,0:num_neig_maxtot))
         neighj_tot = 0
         neighb_tot = 0
         hr_box = 0.0d0

! Loop over atoms
         do iatom = 1,natoms
            neighj_tot(1:neighn_tot(iatom),iatom) =                     &
     &               neighj_aux(1:neighn_tot(iatom),iatom)
            neighb_tot(1:neighn_tot(iatom),iatom) =                     &
     &               neighb_aux(1:neighn_tot(iatom),iatom)
         end do

         deallocate (neighb_aux)
         deallocate (neighj_aux)

         return

! Format Statements
! ===========================================================================

          end subroutine num_neigh_tot
