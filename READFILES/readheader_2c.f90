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

! readheader_2c.f90
! Program Description
! ===========================================================================
!       Read the header of the 2-center data files.
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
        subroutine readheader_2c (interaction, iounit, nsh_max, numz, rc1,   &
     &                            rc2, zmin, zmax, npseudo, cl_pseudo)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: interaction
        integer, intent (in) :: iounit
        integer, intent (in) :: nsh_max
 
! Output
        integer, intent (out) :: npseudo
        integer, intent (out) :: numz
 
        real, intent (out) :: rc1, rc2
        real, intent (out) :: zmin, zmax
 
        real, intent (out), dimension (nsh_max) :: cl_pseudo
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iline
        integer issh
        integer nucz1, nucz2
 
        character (len = 70) message
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        do iline = 1, 9
         read (iounit,100) message
        end do
 
        read (iounit,*) nucz1, rc1
        read (iounit,*) nucz2, rc2
 
! The pseudopotential has an extra 2 lines.
        if (interaction .eq. 5) then
         read (iounit,*) npseudo
         read (iounit,*) (cl_pseudo(issh), issh = 1, npseudo)
        end if
 
        read (iounit,*) zmax, numz
        zmin = 0.0d0
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (a70)
 
        return
        end
