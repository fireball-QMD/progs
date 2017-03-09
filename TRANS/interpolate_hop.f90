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


! interpolate_jop.f90
! Program Description
! ===========================================================================
! perform spline interpolation of nonzero elements of hopping matrix  
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
 subroutine interpolate_hop (index, in1, in2, xin, val)

   use dimensions
   use constants_fireball
   use transport

   implicit none

! Argument Declaration and Description
! ===========================================================================

   integer, intent(in)               :: index    ! index 
   integer, intent(in)               :: in1            ! in1
   integer, intent(in)               :: in2            ! in2
   real, intent(in)                  :: xin            ! r
   real, intent(out)                 :: val            ! val(x)
 
! Local Parameters and Data Declaration
! ===========================================================================
! tol=tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
   real, parameter :: tol = 3.0d-1
   real, parameter :: trash = 0.000001d0

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
   real xmin 
   real xxp
   real xlo
   real xhi
   
! Superspline variables
   real a
   real b
   real flo
   real fhi
   real d2flo
   real d2fhi

 
! Procedure
! ===========================================================================

! get number of mesh points
   nnum = nzh(in1,in2)
! get xmax
   xmax = zh_max(in1,in2)
   xmin = zh_min(in1,in2)

! note : the points must be equally spaced and start at 0 (??)
   h = (xmax - xmin)/(nnum - 1)

! get distance

! check interpolation range
   if (xin .lt. xmin) then
      write (*,*) ' xin, xmin = ', xin, xmin
      write (*,*) '  error in gethop : xin < xmin'
      stop ! Negative value is very bad - should never get there
   else if (xin .ge. xmax) then
      val = hops(nnum,index,in1,in2)
      return
   else if (xin .lt. (xmin+h)) then
      val = hops(1,index,in1,in2)
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
   fhi = hops(ihi,index,in1,in2)
   flo = hops(ilo,index,in1,in2)
   d2fhi = hop_spline(ihi,index,in1,in2)
   d2flo = hop_spline(ilo,index,in1,in2)

! value
   val = a*flo + b*fhi + ((a**3-a)*d2flo + (b**3-b)*d2fhi)*(h**2)/6.0d0

! Format Statements
! ===========================================================================

   return
 end subroutine interpolate_hop
