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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! Stubs for linking against, when we do not have xmgr libraries installed
        subroutine graceinit (time, etotper, getotper)
        implicit none
 
        real, intent (in) :: etotper
        real, intent (in) :: getotper
        real, intent (in) :: time
 
        stop
        end


        subroutine graceupdate (time, etotper, getotper)
        implicit none
 
        real, intent (in) :: etotper
        real, intent (in) :: getotper
        real, intent (in) :: time
 
        stop
        end
