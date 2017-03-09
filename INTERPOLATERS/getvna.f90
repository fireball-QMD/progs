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


! getvna.f90
! Program Description
! ===========================================================================
! getpsi returns neutral atomic potential at given point r 
!
! ===========================================================================
! Code written by:
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
 subroutine getvna (ispec, xin, vnn, dvnn)

   use dimensions
   use constants_fireball
   use vnneutral

   implicit none

! Argument Declaration and Description
! ===========================================================================

   integer, intent(in)               :: ispec    ! ispec 
   real, intent(in)                  :: xin      ! r(:)
   real, intent(out)                 :: vnn      ! vna(x)
   real, intent(out)                 :: dvnn     ! dvna/dx
 
! Local Parameters and Data Declaration
! ===========================================================================
! tol=tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
   real, parameter :: tol = 1.0d-5
   real, parameter :: trash = 1.0d-7

! Local Variable Declaration and Description
! ===========================================================================
!
! -------------------------------------------------------------
   

   integer imid
   integer nnum
   integer ilo
   integer ihi
 
   real h
   real xmax
   real, parameter :: xmin = 0.0d0
   real xxp
   real xlo
   real xhi

! Superspline variables
   real a
   real b
   real vlo
   real vhi
   real d2vlo
   real d2vhi
 
! Procedure
! ===========================================================================

! --------------------------------------------------------------------
!                   -------- Radial part -------
! --------------------------------------------------------------------


! get number of mesh points
   nnum = mesh_na(ispec)
! get xmax of R(l)
   xmax = rmax_na(ispec)
! note : the points must be equally spaced and start at 0
   h = (xmax - xmin)/(nnum - 1)

! check interpolation range
   if (xin .lt. xmin) then
      write (*,*) ' xin, xmin = ', xin, xmin
      write (*,*) '  error in getvna : xin < xmin'
      stop ! Negative value is very bad - should never get there
   else if (xin .gt. xmax) then
      vnn = vnna(nnum,ispec)
      dvnn = 0.0d0
      return  
   else if (xin .lt. (xmin+h)) then
      vnn = vnna(1,ispec)
      dvnn = 0.0d0
      return  
   else
      xxp = xin - xmin
! note : imid is the point to the right (code assumes this)
      imid = int((xxp+trash)/h) + 1
      xhi = real(imid-1)*h  
      xlo = real(imid-2)*h
      ihi = imid
      ilo = imid-1
   end if

! interpolation formulae
   a = (xhi - xxp)/h   
   b = (xxp - xlo)/h
   vhi = vnna(ihi,ispec)
   vlo = vnna(ilo,ispec)
   d2vhi = vnna_spline(ihi,ispec)
   d2vlo = vnna_spline(ilo,ispec)

! value
   vnn = a*vlo + b*vhi + ((a**3-a)*d2vlo + (b**3-b)*d2vhi)*(h**2)/6.0d0
! 1-th derivative
   dvnn = (vhi-vlo)/h - (3.0d0*a**2-1.0d0)*h*d2vlo                       &
 &        + (3.0d0*b**2-1.0d0)*h*d2vhi
! 2-nd derivative
!   d2vnn = a*d2vlo + b*d2vhi


! Format Statements
! ===========================================================================

   return
 end subroutine getvna
