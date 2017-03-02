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
 
! umbrella_energy.f90 
! Program Description
! ===========================================================================
!       This is a subroutine for calculating energy contribution relative 
! to the umbrella sampling potential
!
! ===========================================================================
! Code written by:
! Hilaire Chevreau 
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_umbrella (umb_e, coord_reac_inter, natoms, ratom)
        use umbrella
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

        real, intent (in), dimension (3, natoms) :: ratom

! Ouput
        real, intent (out) :: umb_e          ! Energy for reaction coordinate
        real, intent (out) :: coord_reac_inter  ! Reaction coordinate MD 
 
! Local Parameters and Data Declaration
! ===========================================================================
        
        
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ip

        real distance
 
! Procedure
! ===========================================================================
! Initialize 
        umb_e = 0.0d0
        do ip = 1, umb_pair
         distance = sqrt((ratom(1,umb_p2(ip)) - ratom(1,umb_p1(ip)))**2      &
     &                   + (ratom(2,umb_p2(ip)) - ratom(2,umb_p1(ip)))**2    &
     &                   + (ratom(3,umb_p2(ip)) - ratom(3,umb_p1(ip)))**2)
         umb_e = umb_e + (umb_CFd(ip)*(distance - umb_d0(ip))**2)/2.0d0
        end do

! Reaction coordinate value, at this moment only a distance, corresponding 
! to only one pair...
        coord_reac_inter = distance
 
! Format Statements
! ===========================================================================
        return
        end
