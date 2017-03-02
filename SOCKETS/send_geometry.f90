! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
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
 
! send_geometry.f90
! Program Description
! ===========================================================================
! This sends energy and its derivatives back to the program that is 
! driving fireball via a socket.
! Real work is done in C code.
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine send_geometry(natoms, energy)
        use forces
!       Add something here for Hessian, once we have it
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: natoms
        real, intent(in)    :: energy 
! Local Variable Declaration and Description
! ===========================================================================
        integer realsize
 
! Procedure
! ===========================================================================
        realsize = 4 ! Real(4)
        if (precision(energy) .ge. 10) realsize = 8 ! Real(8)
        call soc_send(energy,realsize)
        call soc_send(ftot,3*natoms*realsize)
!       Hessian is ((3*natoms)**2) * realsize
 
        return
        end

