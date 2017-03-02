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

! sparse_getpacksize.f90
! Program Description
! ===========================================================================
!       This routine gets the size of the packet sent to this processor.
! FIXME This routine can return very, very large numbers for very, very
! sparse matrices.  There must be a better way.
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
        subroutine sparse_getpacksize (nAmax, nprows, packsize)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
! input byte array
        integer, intent (in) :: nAmax
        integer, intent (in) :: nprows

! Output
! packsize = on output, minimum size in bytes necessary for the pack array
        integer, intent (out) :: packsize

! Local Parameters and Data Declaration
! ===========================================================================
        real cplx_dummy
        integer, parameter :: sinteger = kind(packsize)
        integer, parameter :: sreal = kind(cplx_dummy)*2

! Local Variable Declaration and Description
! ===========================================================================

! Procedure
! ===========================================================================
        packsize = sinteger*(3 + nprows*(1 + nAmax)) + 2*sreal*nprows*nAmax

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end

