! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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


! chebft.f90
! Program Description
! ===========================================================================
!	This routine calculates the Chebyshev coefficients for a given 
! function. In particular, it is used to obtain the coefficients for the
! functions x^1/2 and x^-1/2, which in turn is used for the calculation of
! S^1/2 and S-1/2. Also given a range of [a,b], this routine rescales the
! range to [-1,1].
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
      	subroutine chebft (a, b, c, n, func)
      	implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
	integer, intent(in) :: n		! number of coefficients

	real, intent(in) :: a, b                ! range of function
	real, external :: func

! Output
	real, intent(out), dimension (1:n) :: c                ! coefficients

! Local Parameters and Data Declaration
! ===========================================================================
	integer, parameter :: nmax = 50

! Local Variable Declaration and Description
! ===========================================================================
	integer j, k

	real bma
	real bpa
	real fac
	real pi
	real sum
	real y

	real, dimension (1: nmax) :: f

! Procedure
! ===========================================================================
! Initialize pi
        pi = 4.0d0*atan(1.0d0)

        bma = 0.5d0*(b - a)
        bpa = 0.5d0*(b + a)

        do k = 1, n
         y = cos(pi*(k - 0.5d0)/n)
         f(k) = func(y*bma + bpa)
        end do
        fac = 2.0d0/n
        do j = 1, n
         sum = 0.0d0
         do k = 1, n
          sum = sum + f(k)*cos((pi*(j - 1.0d0))*((k - 0.5d0)/n))
         end do
         c(j) = fac*sum
        end do

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end


        real function func12 (x)
        implicit none
        real x
        func12 = x**(0.5d0)
        return
        end


        real function funcm12 (x)
        implicit none
        real x
        funcm12 = 1.0d0/x**(0.5d0)
        return
        end
 
