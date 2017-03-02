! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
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

! fireball.f90
! Program Description
! ===========================================================================
!       This finds the Fdata if it does not exist in cwd.
!
! ===========================================================================
! Code written by:
! J.B. Keith
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine findFdata (fdataLocation)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        character(len = 200), intent(out) :: fdataLocation
        logical qexist

! Procedure
! ===========================================================================

        inquire(file='Fdata/info.dat',exist=qexist)
        if(qexist) then
          fdataLocation='Fdata'
        else
          open(unit=412,file='Fdata.optional',status='unknown',action='read')
          read(412,'(a)')fdataLocation
          close(412)
        end if  

        end subroutine
