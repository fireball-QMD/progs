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

! xmgrupdate.f90
! Program Description
! ===========================================================================
!       This routine updates the xmgr plots that have been previously
! opened and initialized.
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
        subroutine xmgrupdate (time, etotper, getotper)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: etotper
        real, intent (in) :: getotper
        real, intent (in) :: time
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        character (len=64) buffer
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (buffer, '("g0.s0 point", f10.4, ",", f12.4)') time, etotper
        call ACEgrCommandf (buffer)

        write (buffer, '("g0.s1 point", f10.4, ",", f12.4)') time, getotper
        call ACEgrCommandf (buffer)

        call ACEgrCommandf ('autoscale')
        call ACEgrCommandf ('redraw')

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
