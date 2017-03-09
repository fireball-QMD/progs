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


! mpairnay.f90
! Program Description
! ======================================================================
! This is a neighbor routine, that we need for the nonlocal
! pseudopotential. We need this for Dassemble_3c.
!
! < iatom | jatom >. Given iatom, jatom and rdiff=r(iatom)-r(jatom),
!                  then what neighbor to iatom is jatom.
! In other words, what is its "mvalue". We need this
! because the interaction is stored in sVNL(mu,nu,iatom,m).
! So what the heck is m?
!
! This routine is a function call.
!
! Variables are iatom,jatom,rdiff(3). Other input are "constants_fireball".
! The output is mpairnay, which is the m value for iatom, jatom.
! 
! ===========================================================================
! Original code written by Otto F. Sankey.

! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        integer function mpairnay (iatom, jatom, rdiff)

        use configuration  
        use dimensions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: iatom
        integer, intent(in) :: jatom

        real, intent(in), dimension (3) :: rdiff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imatch
        integer ineigh
        integer jjatom
        integer mbeta

        real diff

        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
 
! Procedure
! ===========================================================================
! Initialize
        mpairnay = 0
        r1(:) = ratom(:,iatom)

! Loop over the neighbors of  iatom.
        imatch = 0
        do ineigh = 1, nPPn(iatom)     ! <==== loop 2 over i's neighbors
         jjatom = nPP_j(ineigh,iatom)

! Eliminate the obvious.
         if (jjatom .eq. jatom) then

! OK we have a candidate. If the vector connceting them is
! rdiff(3). Then we have a match. if there is no match, then
! we should not be here. Then there is a problem.
          mbeta = nPP_b(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)

! Vector from 1 to 2 is r21
          r21(:) = r2(:) - r1(:)

! Now compare rdiff to r21.
          diff = (r21(1) - rdiff(1))**2 + (r21(2) - rdiff(2))**2             &
                + (r21(3) - rdiff(3))**2

! We should find only one match.
          if (diff .lt. 0.0001d0) then
           imatch = imatch + 1
           mpairnay = ineigh
          end if
         end if
        end do
 
! Sanity checks
        if (imatch .ne. 1) then
         write (*,*) ' imatch = ', imatch
         write (*,*) ' The variable imatch MUST be ONE! NO EXCEPTIONS '
         write (*,*) ' Bad imatch value in mpairnay.f90; must abort! '
         write (*,*) ' iatom = ', iatom, r1(:)
         write (*,*) ' jatom = ', jatom, rdiff(:) 
         stop
        end if
        if (mpairnay .le. 0) then
         write (*,*) ' Huh? The variable, mpairnay < 0, mpairnay = ', mpairnay
         write (*,*) ' It MUST be > 1. NO EXCEPTIONS '
         write (*,*) ' Bad mpairnay value in mpairnay.f90; must abort! '
         stop
        end if
 
! Format Statements
! ===========================================================================
 
        return
      end function mpairnay
