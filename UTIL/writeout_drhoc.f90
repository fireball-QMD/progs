! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek
!
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


! writeout_xsf.f90
! Program Description
! ===========================================================================
!       The subroutine writes out complex grid densities and does other
! processing
!
! ===========================================================================
! Code written by:
! B. Hess
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine writeout_drhoc (filename, dgridc, maxgrid, gridtol)

   use configuration
   use grid
   use charges
   use interactions
   use options
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

   complex,   dimension (0:nrm-1), intent (in) :: dgridc
!   character (len=40), intent(in) :: filename
   character (len=50) :: filename
   real :: maxgrid, gridtol
!, intent(in)

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer i
   integer j
   integer k
   integer i0
   integer j0
   integer k0
   integer index, inonzero
   real xx, yy, zz
   real xavg, yavg, zavg

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! write out dden into *.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org

   write (*,*) '  Write out complex imported density ',filename
!   write (*,*) 'test grid ',dgridc(1:100)
   write (*,100)
!Find center of grid to subtract off. Determine which points to save
   inonzero = 0
   xavg = 0.0d0
   yavg = 0.0d0
   zavg = 0.0d0
index = 0
   do k = 0, rm3-1
    do j = 0, rm2-1
     do i = 0, rm1-1
	  if (abs(dgridc(index)) .gt. gridtol * maxgrid) inonzero = inonzero + 1
      index = index + 1
      xx = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1)
      yy = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2)
      zz = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3)
      xavg = xavg + xx
   	  yavg = yavg + yy
   	  zavg = zavg + zz
     enddo ! i
    enddo ! j
  enddo ! k
  xavg = xavg/index
  yavg = yavg/index
  zavg = zavg/index

! open file
  open ( unit = 302, file = filename, status = 'unknown' )

  write (302,*) inonzero !, 'nonzero points on density grid'
  write (*,*) 'Number of nonzero gridpoints saved', inonzero

!  inonzero = 0
  index = 0
  do k = 0, rm3-1
   do j = 0, rm2-1
    do i = 0, rm1-1
     if (abs(dgridc(index)) .gt. gridtol * maxgrid) then
!     inonzero = inonzero + 1
      xx = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1) - xavg - elvec(1,1)
      yy = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2) - yavg - elvec(2,2)
      zz = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3) - zavg - elvec(3,3)
	  write (302,'(5f14.8)') xx,yy,zz,dgridc(index)
	 end if
	 index = index + 1
    enddo ! i
   enddo ! j
  enddo ! k
!  write(*,*), 'inonzero2 written', inonzero

   write (302,*) 'ATOMS'
   do iatom = 1,natoms
!   write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
    write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom(i,iatom),i=1,3)
   enddo


  write (302,*)
! close file
  close (302)
! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (70('='))
200     format (e14.6)

        return
      end subroutine writeout_drhoc
