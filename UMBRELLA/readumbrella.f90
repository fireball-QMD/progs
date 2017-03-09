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


! readumbrella.f90
! Program Description
! ===========================================================================
!       This routine reads in the information from the umbrella.optional file.
!
! ===========================================================================
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
        subroutine readumbrella (natoms, ratom, nstepf)
        use umbrella
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nstepf
 
        real, intent (in), dimension (3, natoms) :: ratom      

! Local Parameters and Data Declaration
! ===========================================================================
        integer ip
 
        real dist_XY

! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
         write (*,*) '  '
         write (*,*) ' Reading information from the umbrella.optional file. '
         write (*,100)
         write (*,*) ' You are doing umbrella sampling MD calcvulations so we '
         write (*,*) ' read in the information from the umbrella.optional file.'
         write (*,*) ' The umbrella analysis must be made with the iquench = 0 '
         write (*,*) ' option and with T .ne. 0 : in script.input and '
         write (*,*) ' quench.optional. '
         write (*,*) '  '  
 
         open (unit = 78, file = 'umbrella.optional', status = 'old')
         write(*,*) ' Read in name of output umbrella sampling analysis file'
         read (78,200) umb_output_file
         write (*,200) umb_output_file
         write (*,*) '  '

         write (*,*) ' Read in : time before starting the umbrella analysis'
         read (78,300) umb_time_start 
         write (*,*)  umb_time_start
         write (*,*) '  '

         write (*,*) ' Read in the number of atoms pairs'
         read (78,400) umb_pair
         write (*,*) ' Read in from umbrella.optional, umb_pair = ', umb_pair
         call allocate_umb (natoms) 
 
         write (*,*) '  '
         write (*,*) ' Read in atom A number, atom B number, force constant '
         do ip = 1, umb_pair
          read (78,*) umb_p1(ip), umb_p2(ip) 
          read (78,*) dist_XY
          read (78,*) umb_CFd(ip)
          write (*,*) ' write pair ', ip
          write (*,*) ' atom1 = ', umb_p1(ip), ' atom2', umb_p2(ip)
          write (*,*) ' d0 ', dist_XY, ' Constant force ', umb_CFd(ip)
          write (*,*) '  '

! Initialization of the reference distance used in the constrained potential...
! Thus is the distance found in the input geom...

          umb_d0(ip) = dist_XY               
         end do

! Initialisation of the initial reaction coordinate used for the definition of
! umbrella sampling window.
!
! Here this is simply the distance between A and B.
! (other definition may be coded, when many atomic pairs are of concerns.)
         umb_react_coord_init = dist_XY     
         close (unit = 78)  
         write (*,100)
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a30)  
300     format (f15.6)
400     format (i2)

        return
        end
