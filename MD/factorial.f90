! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! factorial.f90
! Program Description
! ===========================================================================
!       This computes the factorial and returns it as a real*4
!
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
        real function factorial (ifac)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: ifac
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer jj
        integer countit
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if (ifac .lt. 0) then
         write (*,*) ' Error in factorial  ----  ifac < 0 '
         stop
        else
         countit = 1
! Note: if ifac=0,1 then this do loop is skipped, and factorial=1
         do jj = 2, ifac
          countit = countit*jj
         end do
        end if
        factorial = real(countit)
 
! Format Statements
! ===========================================================================
        return
        end
