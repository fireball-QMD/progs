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


! interpolate_2d.f90
! Program Description
! ===========================================================================
!       This routine is a two-dimensional interpolater on a 4x4 sub-grid.
!
! ===========================================================================
! Code rewritten by:
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
        subroutine interpolate_2d (xin, yin, iforce, nx, ny, hx, hy,  &
     &                             xintegral, Q_L, dQ_Ldx, dQ_Ldy)
        use dimensions
        use interactions
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce

! Number of points in 3-center integrals
        integer, intent(in) :: nx
        integer, intent(in) :: ny

        real, intent (in) :: xin
        real, intent (in) :: yin

! Difference between points along x and y
        real, intent (in) :: hx
        real, intent (in) :: hy

        real, intent (in), dimension (numXmax, numYmax) :: xintegral
 
! Output
        real, intent (out) :: Q_L    ! the contibutions for a matrix element
        real, intent (out) :: dQ_Ldx ! d/dx Q_L (Computed only if iforce = 1)
        real, intent (out) :: dQ_Ldy ! d/dy Q_L (Computed only if iforce = 1)

! Local Parameters and Data Declaration
! ===========================================================================
! If you think that lower order interpolation for slowly changing parts
! of the surface are bad, set these two to zero.
        real, parameter :: tiny = 1.0d-5
        real, parameter :: small= 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer imidx, imidy
        integer k

        real f0p3, f0p6, f1m1, f1m2, f1m3, f1p3, f1p6, f2p1, flm1
        real gradtest
        real gradx, grady
        real px, py
        real, parameter :: xmin = 0
        real, parameter :: ymin = 0

        real bb0,bb1,bb2,bb3
! -1 to 2 since we do third order
        real, dimension (-1:2,-1:2) :: fun
        real, dimension (-1:2) :: g, gp

! Procedure
! ===========================================================================
! We assume xmin = ymin = 0.0d0 always in our interpolations.
! We assume that x<xmin and x>xmax has been checked for
! ===========================================================================
! We need to find what point of the grid to use
        imidx = int((xin - xmin)/hx) + 1
        if (imidx .lt. 2) then
          imidx = 2
        else if (imidx .gt. nx - 2) then
          imidx = nx - 2
        end if

        imidy = int((yin - ymin)/hy) + 1
        if (imidy .lt. 2) then
          imidy = 2
        else if (imidy .gt. ny - 2) then
          imidy = ny - 2
        end if

        px = xin/hx - (imidx - 1)
        py = yin/hy - (imidy - 1)
 
! Now copy grid values to fun array for third order
! Note:  F90 syntax is slower f(-1:2,-1:2) = xintegral(imidx - 1:imidx + 2,imidy - 1:imidy + 2)
        fun(-1,-1) = xintegral(imidx - 1,imidy - 1)
        fun(0,-1) = xintegral(imidx,imidy - 1) 
        fun(1,-1) = xintegral(imidx + 1,imidy - 1)
        fun(2,-1) = xintegral(imidx + 2,imidy - 1)

        fun(-1, 0) = xintegral(imidx - 1,imidy)
        fun(0, 0) = xintegral(imidx,imidy) 
        fun(1, 0) = xintegral(imidx + 1,imidy) 
        fun(2, 0) = xintegral(imidx + 2,imidy)

        fun(-1, 1) = xintegral(imidx - 1,imidy + 1)
        fun(0, 1) = xintegral(imidx,imidy + 1)
        fun(1, 1) = xintegral(imidx + 1,imidy + 1)
        fun(2, 1) = xintegral(imidx + 2,imidy + 1)

        fun(-1, 2) = xintegral(imidx - 1,imidy + 2)
        fun(0, 2) = xintegral(imidx,imidy + 2)
        fun(1, 2) = xintegral(imidx + 1,imidy + 2)
        fun(2, 2) = xintegral(imidx + 2,imidy + 2)
 
! ===========================================================================
! Adaptive interpolation - estimate gradient:
        gradx = (fun(1,0) - fun(0,0))/hx
        grady = (fun(0,1) - fun(0,0))/hy
        gradtest = abs(gradx) + abs(grady)

        if (gradtest .lt. tiny) then

!
! METHOD 1
!

! Do three point linear bivariate interpolation
! Handbook of Mathematical Functions..., edited by M. Abramowitz
! and I.A. Stegun, Dover edition, Pg. 882, Eq. 25.2.65
         Q_L = (1.0d0 - px - py)*fun(0,0) + px*fun(1,0) + py*fun(0,1)
         if (iforce .eq. 1) then
          dQ_Ldx = gradx
          dQ_Ldy = grady
         end if
        else if (gradtest .lt. small) then

!
! METHOD 2
!

! Do quadratic bivariate interpolation (six point formula, Eq. 25.2.67)
         Q_L = py*(py-1)*0.5d0*fun(0,-1) + px*(px-1)*0.5d0*fun(-1,0)         &
     &        + (1 + px*py - px*px - py*py)*fun(0,0)                         &
     &        + px*(px - 2*py + 1)*fun(1,0)*0.5d0                            &
     &        + py*(py - 2*px + 1)*fun(0,1)*0.5d0                            &
     &        + px*py*fun(1,1)
         if (iforce .eq. 1) then
          dQ_Ldx = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*py           &
     &            + (fun(-1,0) + fun(1,0) - 2.0d0*fun(0,0))*px               &
     &            - 0.5d0*(fun(-1,0) - fun(1,0)))/hx
          dQ_Ldy = ((fun(1,1) - fun(1,0) - fun(0,1) + fun(0,0))*px           &
     &            + (fun(0,-1) + fun(0,1) - 2.0d0*fun(0,0))*py               &
     &            - 0.5d0*(fun(0,-1) - fun(0,1)))/hy
         end if

!
! METHOD 3
!

! Interpolate one direction, then interpolate using these values to get
! the final value.  Use Eq. 25.2.13 (4 point intepolater).
! Total number of points used is 16
        else if (D2intMeth .eq. 1) then
         do k = -1, 2
          f1m1 = fun(k,-1)
          f1m2 = 2*f1m1
          f1m3 = 3*f1m1

          f0p3 = 3*fun(k,0)
          f0p6 = 2*f0p3

          f1p3 = 3*fun(k,1)
          f1p6 = 2*f1p3

          f2p1 = fun(k,2)

          bb3 = - f1m1 + f0p3 - f1p3 + f2p1
          bb2 = f1m3 - f0p6 + f1p3
          bb1 = - f1m2 - f0p3 + f1p6 - f2p1
          bb0 = f0p6

          g(k) = ((bb3*py + bb2)*py + bb1)*py + bb0
          if (iforce .eq. 1) gp(k) = ((3*bb3*py + 2*bb2)*py + bb1)
         end do

         f1m1 = g(-1)
         f1m2 = 2*f1m1
         f1m3 = 3*f1m1
    
         f0p3 = 3*g(0)
         f0p6 = 2*f0p3
 
         f1p3 = 3*g(1)
         f1p6 = 2*f1p3

         f2p1 = g(2) 

         bb3 = - f1m1 + f0p3 - f1p3 + f2p1
         bb2 = f1m3 - f0p6 + f1p3
         bb1 = - f1m2 - f0p3 + f1p6 - f2p1
         bb0 = f0p6

         Q_L = (((bb3*px + bb2)*px + bb1)*px + bb0)/36.0d0

         if (iforce .eq. 1) then
           dQ_Ldx = ((3*bb3*px + 2*bb2)*px + bb1)/(36.0d0*hx)

           f1m1 = gp(-1)
           f1m2 = 2*f1m1
           f1m3 = 3*f1m1
    
           f0p3 = 3*gp(0)
           f0p6 = 2*f0p3
    
           f1p3 = 3*gp(1)
           f1p6 = 2*f1p3

           f2p1 = gp(2) 

           bb3 = - f1m1 + f0p3 - f1p3 + f2p1
           bb2 = f1m3 - f0p6 + f1p3
           bb1 = - f1m2 - f0p3 + f1p6 - f2p1
           bb0 = f0p6

           dQ_Ldy = (((bb3*px + bb2)*px + bb1)*px + bb0)/(36.0d0*hy)

         end if
        else
          write(*,*) ' Invalid interpolation in interpolate_2d '
          write(*,*) ' Some methods that sucked have been removed '
          stop
        end if

        return
! A final note, if you are interested in splines, you should start with
! "Handbook on SPLINES for the User"
! by Eugene V. Shikin and Alexander I. Plis.  1995, CRC Press.  Most other
! books on slines and bivariate interpolation are nothing but proofs and
! abstract math, but this one gives the real equations you need to get
! the actual work done.

! Format Statements
! ===========================================================================
 
        return
        end
