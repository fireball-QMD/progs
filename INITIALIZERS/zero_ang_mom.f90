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


! zero_ang_mom.f90
! Program Description
! ===========================================================================
!       Contains subroutine zero_ang_mom which takes a set of random
! velocity and adjusts them to get total angular momentum = 0.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
! 
! Program Declaration
! ===========================================================================
        subroutine zero_ang_mom ()
        use dimensions
        use configuration
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input: 
     
! Output:

! Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix

        real, dimension (3) :: crossa 
        real, dimension (3, 3) :: xinertia 

        real, dimension (3, 3) :: xinvert 
        real, dimension (3) :: xlcm 
        real, dimension (3) :: wvec 

! Procedure
! ===========================================================================
        xlcm = 0.0d0
        do iatom = 1, natoms 
         xlcm(1) = xlcm(1) + (ratom(2,iatom)*xmass(iatom)*vatom(3,iatom) -  & 
     &                        ratom(3,iatom)*xmass(iatom)*vatom(2,iatom))
         xlcm(2) = xlcm(2) + (ratom(3,iatom)*xmass(iatom)*vatom(1,iatom) -  &
     &                        ratom(1,iatom)*xmass(iatom)*vatom(3,iatom))
         xlcm(3) = xlcm(3) + (ratom(1,iatom)*xmass(iatom)*vatom(2,iatom) -  &
     &                        ratom(2,iatom)*xmass(iatom)*vatom(1,iatom))
        end do

        write (*,*) '  '
        write (*,100) xlcm

! Calculate the inertia tensor
! I(i,j) = sumoverk m(k) * ( r**2(k) delk(i,j)  -  ri(k) * rj(k) )
        xinertia = 0.0d0
        do iatom = 1, natoms 
         xinertia(1,1) = xinertia(1,1)                                      &
     &    + xmass(iatom)*(ratom(2,iatom)**2 + ratom(3,iatom)**2)
         xinertia(1,2) = xinertia(1,2)                                      &
     &    - xmass(iatom)*(ratom(1,iatom)*ratom(2,iatom))
         xinertia(1,3) = xinertia(1,3)                                      &
     &    - xmass(iatom)*(ratom(1,iatom)*ratom(3,iatom))
         xinertia(2,1) = xinertia(2,1)                                      &
     &    - xmass(iatom)*(ratom(2,iatom)*ratom(1,iatom))
         xinertia(2,2) = xinertia(2,2)                                      &
     &    + xmass(iatom)*(ratom(1,iatom)**2 + ratom(3,iatom)**2)
         xinertia(2,3) = xinertia(2,3)                                      &
     &    - xmass(iatom)*(ratom(2,iatom)*ratom(3,iatom))
         xinertia(3,1) = xinertia(3,1)                                      &
     &    - xmass(iatom)*(ratom(3,iatom)*ratom(1,iatom))
         xinertia(3,2) = xinertia(3,2)                                      &
     &    - xmass(iatom)*(ratom(3,iatom)*ratom(2,iatom))
         xinertia(3,3) = xinertia(3,3)                                      &
     &    + xmass(iatom)*(ratom(1,iatom)**2 + ratom(2,iatom)**2)
        end do
 
! Here the inertia tensor has units of mp*A**2. This is okay
! because the units will work out in the end.
        call invert3x3 (xinertia, xinvert) 

! L = I - dot - omega  so   omega = I(inverse) - dot - L
        wvec = 0.0d0
        do ix = 1, 3
         wvec(:) = wvec(:) + xinvert(:,ix)*xlcm(ix)
        end do
 
        do iatom = 1, natoms
         crossa(1) = wvec(2)*ratom(3,iatom) - wvec(3)*ratom(2,iatom)
         crossa(2) = wvec(3)*ratom(1,iatom) - wvec(1)*ratom(3,iatom)
         crossa(3) = wvec(1)*ratom(2,iatom) - wvec(2)*ratom(1,iatom)
         vatom(:,iatom) = vatom(:,iatom) - crossa
        end do

        xlcm = 0.0d0
        do iatom = 1, natoms 
         xlcm(1) = xlcm(1) + (ratom(2,iatom)*xmass(iatom)*vatom(3,iatom) -   &
     &                        ratom(3,iatom)*xmass(iatom)*vatom(2,iatom))
         xlcm(2) = xlcm(2) + (ratom(3,iatom)*xmass(iatom)*vatom(1,iatom) -   &
     &                        ratom(1,iatom)*xmass(iatom)*vatom(3,iatom))
         xlcm(3) = xlcm(3) + (ratom(1,iatom)*xmass(iatom)*vatom(2,iatom) -   &
     &                        ratom(2,iatom)*xmass(iatom)*vatom(1,iatom))
        end do

        write (*,200) xlcm

! Format Statements
! ===========================================================================
100     format (2x, ' initial lcm = ', 3d16.7)
200     format (2x, '   final lcm = ', 3d16.7)

        return
        end
