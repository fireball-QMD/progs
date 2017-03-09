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


! buildspline_1d.f90
! Program Description
! ===========================================================================
!       This code takes a 1D integral array and generates a single natural 
! cubic spline to fit the entire curve.  Values will later be interpolated 
! based on this "super spline".  The spline and its derivatives are continuous.
! See reference in Numerical Analysis by R. L. Burden and J. D. Faires, 
! section 3.6
! There are several types of splines: not-a-knot, natural spline,
! clamped spline, extrapolated spline, parabolically terminated spline,
! endpoint curvature adjusted spline
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine buildspline_1d (integral, numz, itype, in1, in2, xmax,   &
     &                             interaction)
        use dimensions
        use constants_fireball
        use integrals

        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: integral
        integer, intent(in) :: numz
        integer, intent(in) :: itype
        integer, intent(in) :: in1
        integer, intent(in) :: in2
        integer, intent(in) :: interaction

        real, intent(in) :: xmax
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imid
        integer iorder
        integer norder
        integer numz_used

        real fn
        real fnm1
        real fo
        real fop1
        real fpo
        real fpn
        real h
        real xmin

        real, dimension (0: nfofx) :: a
        real, dimension (0: nfofx) :: b
        real, dimension (0: nfofx) :: c
        real, dimension (0: nfofx) :: d
        real, dimension (0: nfofx) :: alpha
        real, dimension (0: nfofx) :: L
        real, dimension (0: nfofx) :: mu
        real, dimension (0: nfofx) :: Z

! Procedure
! ===========================================================================
! Note : the points must be equally spaced!!!!
        xmin = 0.0d0
        h = (xmax - xmin)/(numz - 1)

! Almost everything goes to zero
!        if ((interaction.ne.12).and.(itype.ne.ind2c(4,0))) xintegral_2c(integral,numz,itype,in1,in2) = 0  ! JOM-ENRIQUE

        numz_used = numz
        checker: do iorder = numz, 3, -1
          if (xintegral_2c(integral,iorder,itype,in1,in2) .eq. 0) then
              numz_used = iorder
          else
              exit checker
          end if
        end do checker
          
! Initialize some variables
        fn = xintegral_2c(integral,numz_used,itype,in1,in2)
        fnm1 = xintegral_2c(integral,numz_used-1,itype,in1,in2)
        fo = xintegral_2c(integral,1,itype,in1,in2)
        fop1 = xintegral_2c(integral,2,itype,in1,in2)
! Get estimates of derivatives
        fpn = (fn - fnm1)/h
        fpo = (fop1 - fo)/h 
! Slight correction for 1/r derivative JOM-ENRIQUE
        if ((interaction.eq.12).or.(interaction.eq.4)) then
	  if (itype.ne.ind2c(4,0)) then
       	    fpn = fpn * (xmax - h)/xmax
	  end if
	end if

! We are not doing natural splines anymore, but rather we now do clamped
! Cubic splices: "clamped" splines with f'(x) given at both endpoints.
        norder = numz_used - 1
        do iorder = 0, norder
         a(iorder) = xintegral_2c(integral,iorder+1,itype,in1,in2)
        end do

        alpha(0) = 3.0d0*(a(1) - a(0))/h - 3.0d0*fpo
        do iorder = 1, norder - 1
         alpha(iorder) = 3.0d0*(a(iorder+1) - 2.0d0*a(iorder) + a(iorder-1))/h
        end do
        alpha(norder) = 3.0d0*fpn - 3.0d0*(a(norder) - a(norder-1))/h

        L(0) = 2.0d0*h
        mu(0) = 0.5d0
        Z(0) = alpha(0)/L(0)
        do iorder = 1, norder - 1
         L(iorder) = (4.0d0 - mu(iorder-1))*h
         mu(iorder) = h/L(iorder)
         Z(iorder) = (alpha(iorder) - h*Z(iorder-1))/L(iorder)
        end do
        L(norder) = (2.0d0 - mu(norder-1))*h
        mu(norder) = 0.0d0
        Z(norder) = (alpha(norder) - h*Z(norder-1))/L(norder) 
        c(norder) = Z(norder)

        do iorder = norder - 1, 0, -1
         c(iorder) = z(iorder) - mu(iorder)*c(iorder+1)
         b(iorder) = (a(iorder+1) - a(iorder))/h                             &
     &              - h*(c(iorder+1) + 2.0d0*c(iorder))/3.0d0
         d(iorder) = (c(iorder+1) - c(iorder))/(3.0d0*h)
        end do
        b(norder) = 0.0d0
        d(norder) = 0.0d0

! We now have a, b, c, d (the coefficients to the spline segments)
! Now copy these into splineint_2c
        do iorder = 1, numz_used
         splineint_2c(1,integral,iorder,itype,in1,in2) = a(iorder-1)
         splineint_2c(2,integral,iorder,itype,in1,in2) = b(iorder-1)
         splineint_2c(3,integral,iorder,itype,in1,in2) = c(iorder-1)
         splineint_2c(4,integral,iorder,itype,in1,in2) = d(iorder-1)
        end do
! Set end to zero if necessary
        if (numz_used .ne. numz) then
         do iorder = numz_used, numz
          splineint_2c(:,integral,iorder,itype,in1,in2) = 0
         end do
        else if (splineint_2c(1,integral,numz,itype,in1,in2) .eq. 0) then
         splineint_2c(2,integral,numz,itype,in1,in2) = 0.0d0
         splineint_2c(3,integral,numz,itype,in1,in2) = 0.0d0
         splineint_2c(4,integral,numz,itype,in1,in2) = 0.0d0
        end if

! Format Statements
! ===========================================================================

        return
        end
