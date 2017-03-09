! Copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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


! readgrid.f90
! Program Description
! ===========================================================================
!       This routine reads information about a numerical grid used to solve
! Poisson equation.
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine readgrid ( iwrtewf )

   use grid
   use charges
   use interactions
   use configuration
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

   integer, intent (in) :: iwrtewf

! Local Parameters and Data Declaration
! ===========================================================================
!   integer, parameter :: max_mesh = 5000

! Local Variable Declaration and Description
! ===========================================================================
   integer ispec
   integer imu
   integer issh
   integer in1
   character (20) :: tmp
   character (len = 30) str
   logical isfile
   logical isstr
   Namelist /mesh/ Ecut, iewform, npbands, pbands, ewfewin_max, ewfewin_min, &
      &            ifixg0, g0

! Allocate Arrays
! ===========================================================================


! Procedure
! ===========================================================================
! default settings
   Ecut = 100.0d0
   if (iwrtewf .eq. 1) then
    iewform = 2
    npbands = 0
    ewfewin_min =-30.0d0
    ewfewin_max = 10.0d0
   endif
   ifixg0 = 0
   g0 (:) = 0.0d0

   write (*,*) '  '
   write (*,*) ' ---------------------  DEFINE GRID  ----------------------*'
   inquire (file = initfile, exist = isfile)
! file fireball.in exists so let's read it
   if (isfile) then
! open input file
    write (*,*) '  Openning fireball.in file'
    open (unit = 17, file = initfile, status = 'old')

! section GRID
! check if the section QUENCH is listed
    str = '&MESH'
    call issection (str,isstr)
    if (isstr) read (60, NML=mesh)
   endif

! read data
! Energy cutoff
   write (*,200) Ecut

   if (iwrtewf .eq. 1) then
! write format 1 .. individual bands; 2.. energy window
    if (iewform .eq. 1) then
! write list of the bands
     write (*,*) ' Number of plotted bands: ',npbands
     if (npbands .gt. npbands_max) then
      write (*,*) 'Number of the bands exceeded max. value:',npbands_max
      write (*,*) 'Please decrease it and try again'
      stop
     endif
     write (*,*) ' List of the selected bands: ',(pbands(imu), imu=1,npbands)
    else
! read energy window
     write (*,*) ' The selected energy windows (emin,emax): '
     write (*,500) ewfewin_min, ewfewin_max
! allocate aux arrays
    endif ! if (iewform)
   endif

! read flag to fix the initial point of the grid
   write (*,600) ifixg0
   if (ifixg0 .eq. 1) then
! read the initial point of the grid
     write (*,610),(g0(imu),imu=1,3)
! avoid move atoms if some of them is 0,0,0
     ishiftO  = 0
   endif

   write (*,*) '  '
   write (*,*) ' ---------------------  END READGRID  ----------------------*'

! close the section GRID
   close (17)

! write out information into param.dat
   open (unit = 50, file = 'param.dat', position='append' ,status='old')
   write (50, *) ' MESH:'
   write (50, *) '  Ecut [Ry]         : ',Ecut
   write (50, *) '  ifixg0            : ',ifixg0
   if (iwrtewf .eq. 1) then
    write (50, *) '  iewform           : ',iewform
    if (iewform .eq. 1) then
     write (50, *) '  npbands           : ',npbands
     do imu = 1,npbands
      write (50, *) '    band no.        : ',pbands(imu)
     enddo
    else
      write (50,*) '  ewfewin_min       : ',ewfewin_min
      write (50,*) '  ewfewin_max       : ',ewfewin_max
    endif
   endif

   close (50)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (a25)
200     format (4x,'Ecut = ', f12.6)
300     format (4x,i3,2x,9f6.2)
400     format (a50)
500     format (4x,'Range :',2f12.6)
600     format (4x,'ifixg0 =',i2)
610     format (4x,'g0 :',3f12.6)

  return
end subroutine readgrid

