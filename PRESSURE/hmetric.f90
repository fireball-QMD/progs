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


! hmetric.f90
! Program Description
! ===========================================================================
!       This routine calculates the metric g = hdagger h, and its derivative
! The subroutine sets an h matrix for a1vec, a2vec and a3vec, and then makes
! a metric tensor: g=h'*h, where h' - is a transpose of h.
!
! input hmet = [ a1(x)  a2(x)  a3(x)  ]
!              [ a1(y)  a2(y)  a3(y)  ]
!              [ a1(z)  a2(z)  a3(z)  ]
!
! input hp = d hmet / dt = hmetdot
!
! The algorithm is from:
! S. Nose, M.L. Klein, Molec. Phys. 50, 1055-1076 (1983)
!
! ===========================================================================
! Original code written by Alex A. Demkov.
 
! Code rewritten by:
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
        subroutine hmetric (hp, g, gdot, hmet, volume)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: hp
        real, intent (in), dimension (3, 3) :: hmet
 
! Output:
        real, intent (out) :: volume
        real, intent (out), dimension (3, 3) :: g
        real, intent (out), dimension (3, 3) :: gdot
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ix
 
        real determinant
 
        real, dimension (3, 3) :: hmett                   ! transpose of hmet
        real, dimension (3, 3) :: hpt
 
! Procedure
! ===========================================================================
! Set up hmett matrix (transpose of hmet)
        hmett(1,1) = hmet(1,1)
        hmett(1,2) = hmet(2,1)
        hmett(1,3) = hmet(3,1)
        hmett(2,1) = hmet(1,2)
        hmett(2,2) = hmet(2,2)
        hmett(2,3) = hmet(3,2)
        hmett(3,1) = hmet(1,3)
        hmett(3,2) = hmet(2,3)
        hmett(3,3) = hmet(3,3)
 
! Set up the g matrix= hdagger*hmet
        g = 0.0d0
        do ix = 1, 3
         g(1,:) = g(1,:) + hmett(1,ix)*hmet(ix,:)
         g(2,:) = g(2,:) + hmett(2,ix)*hmet(ix,:)
         g(3,:) = g(3,:) + hmett(3,ix)*hmet(ix,:)
        end do
 
! Get the determinant of h in a fancy way!
        determinant =                                                        &
     &     hmet(1,1)*(hmet(2,2)*hmet(3,3) - hmet(2,3)*hmet(3,2))             &
     &   - hmet(1,2)*(hmet(2,1)*hmet(3,3) - hmet(2,3)*hmet(3,1))             &
     &   + hmet(1,3)*(hmet(2,1)*hmet(3,2) - hmet(2,2)*hmet(3,1))
 
! The volume is sqrt(determinant**2) for volume. This is because determinant
! may be negative. Also, find that
! dV/dh(alpha,beta) = (determinant**2/V)*(h**-1)((beta,alpha)
!                   = V(h**-1)(a,b)A
! which is in agreement with formula 2.9 of Nose (Molec. Phys. 50, 1055 (1983).)
        volume = sqrt(determinant*determinant)
        write (*,*) '  '
        write (*,*) ' volume = ', volume
 
! Set up hpt matrix (hpt = transpose of hp where hp= d hmet /dt )
        hpt(1,1) = hp(1,1)
        hpt(1,2) = hp(2,1)
        hpt(1,3) = hp(3,1)
        hpt(2,1) = hp(1,2)
        hpt(2,2) = hp(2,2)
        hpt(2,3) = hp(3,2)
        hpt(3,1) = hp(1,3)
        hpt(3,2) = hp(2,3)
        hpt(3,3) = hp(3,3)
 
! Need two products: hpt*hmet and at*hpt
        gdot = 0.0d0
        do ix = 1, 3
         gdot(1,:) = gdot(1,:) + hpt(1,ix)*hmet(ix,:) + hmett(1,ix)*hp(ix,:)
         gdot(2,:) = gdot(2,:) + hpt(2,ix)*hmet(ix,:) + hmett(2,ix)*hp(ix,:)
         gdot(3,:) = gdot(3,:) + hpt(3,ix)*hmet(ix,:) + hmett(3,ix)*hp(ix,:)
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
