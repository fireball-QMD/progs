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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! append_string.f90
! Program Description
! ===========================================================================
!       This function takes a string and appends a given extension to it.
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
        function append_string (filename, extension)
        implicit none
        character (len = 200) append_string
 
! Argument Declaration and Description
! ===========================================================================
        character (len = 200) filename
        character (len = 200) extension
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer lenx
        integer leny
        integer len_trim
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
      lenx = len_trim(filename)
      if (lenx .eq. 0) then
       append_string = ' '
      else
       leny = len_trim(extension)
       if (leny .eq. 0) then
        append_string = filename(1:lenx)
       else
        !lenx = min(lenx,200)
        append_string = filename(1:lenx)//extension(1:leny)
       end if
      end if
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
