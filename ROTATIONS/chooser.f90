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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! chooser.f90
! Program Description
! ===========================================================================
!       This routine prepares the L (left) and R (right) matrices for the
! rotation of L*Matrix*R.
!
! ===========================================================================
! Original code written by Alex Demkov.
!
! Code rewritten by:
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
        subroutine chooser (l, dmat, pmat, rmatrix)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        integer, intent(in) :: l
 
        real, intent(in), dimension(5, 5) :: dmat
        real, intent(in), dimension(3, 3) :: pmat
 
! Output:
        real, intent(out), dimension(5, 5) :: rmatrix
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the rotation matrix
        rmatrix = 0.0d0
 
        if (l .eq. 0) then
         rmatrix(1,1) = 1.0d0
        else if (l .eq. 1) then
         rmatrix(1:3,1:3) = pmat
        else if (l .eq. 2) then
         rmatrix = dmat
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
