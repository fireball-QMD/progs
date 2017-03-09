! copyright info:
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


! getsections.f90
! Program Description
! ===========================================================================
! A subroutine contains auxiliar subroutines to handle string operations
! ===========================================================================
! Code rewritten by:
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
        subroutine issection (string, status)

        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        character (len = 30), intent(in) :: string
        logical, intent (out) :: status
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        character (len=30) line
        integer eof
        logical eqstr

! Procedure
! ===========================================================================
! open in file  
!        open(unit=50,file=initfile,status='unknown') it is already open - unit=60

! initialize status        
        status = .false.  

		rewind 60
! loop through the file        
        do 
         read (60,*, iostat=eof) line
! compare two strings (no matter on upper & lower case)
         call compare_strings (line, string, eqstr) 
         if (eqstr) then
          status = .true.
          exit
         endif
! end of file
         if (eof < 0) exit
       enddo
   
! rewind file   
		rewind 60
  
! Format Statements
! ===========================================================================
        return
      end subroutine issection
      
! aux subroutine
      subroutine compare_strings (str1, str2, flag)
      
      implicit none

! Argument Declaration and Description
! ===========================================================================
! Input     
      character (len = 30), intent(in) :: str1
      character (len = 30), intent(in) :: str2
! Output
      logical, intent (out) :: flag

!  Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! =========================================================================== 
      character (len=30) :: str1a
      character (len=30) :: str2a
      integer i
      
! Procedure
! ===========================================================================  
      flag = .false.
  
! copy originals  
      str1a = str1
      str2a = str2
  
! now shift lower case to upper case
      do i = 1, len(str1a)
       if ( str1a(i:i) >= 'a' .and. str1a(i:i) <= 'z') then
        str1a(i:i) = achar (iachar (str1a(i:i))-32)
       endif   
      enddo  
      do i = 1, len(str2a)
       if ( str2a(i:i) >= 'a' .and. str2a(i:i) <= 'z') then
        str2a(i:i) = achar (iachar (str2a(i:i))-32)
       endif   
      enddo    

      if ( str1a .eq. str2a ) flag = .true.
      return   
      end subroutine compare_strings
      
