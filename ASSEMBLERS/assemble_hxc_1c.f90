! copyright info:
!
!                             @Copyright 2002
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

 
! assemble_1c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the one-center exchange-correlation
! interactions. The results are stored in vxc_1c and etotxc_1c.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_hxc_1c (natoms, itheory, iforce)
        use charges
        use dimensions
        use interactions
        use neighbor_map
        use energy
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iforce
        integer, intent(in) :: itheory
        integer, intent(in) :: natoms
 
! Output

 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu,inu
        integer in1
        integer matom
        integer ixc
 
        real dccexc_1c
        real exc_1c
        real muexc_1c
        real, dimension (numorb_max,numorb_max) :: mu1xc
 
! Procedure
! ===========================================================================
! Initialize
        etotxc_1c = 0.0d0
        Uexc_1c = 0.0d0
        Umuxc_1c = 0.0d0
        vxc_1c = 0.0d0

! Loop over all atoms.
        do iatom = 1, natoms
         matom = neigh_self(iatom)
         in1 = imass(iatom)
         ! set horsfield flag
         ixc = 0
         call unocentros (in1, iatom, iforce, itheory, ixc, exc_1c, muexc_1c, &
      &                   dccexc_1c, mu1xc)
         
  
! Now add Integral n*(exc-muxc) d^3r to etotxc_1c - The atom "overcounting"
! correction. We set derivative of this EQUAL TO ZERO. Certainly the neutral
! atom part has no derivative, but the charge part does have a derivative.
!        write (*,*) ' No derivative of overcounting xc correction. '
!        write (*,*) ' OK only for neutral atoms.'
!        write (*,*) ' Maybe we should fix this later! '
 
         etotxc_1c = etotxc_1c + dccexc_1c

! Note vxc is initialized in initnanlxc above.
         do imu = 1, num_orb(in1)
          do inu = 1, num_orb(in1)
           vxc_1c(imu,inu,matom,iatom) =                                &
     &     vxc_1c(imu,inu,matom,iatom) + mu1xc(imu,inu) 
          end do
         end do
 
! End loop over atoms
        end do
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine assemble_hxc_1c
 
