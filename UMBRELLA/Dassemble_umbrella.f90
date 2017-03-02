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
 
! get_Dumbrella.f90
! Program Description
! ===========================================================================
!       This is a subroutine for calculating forces contribution relative 
! relative to the umbrella sampling potential
!
! ===========================================================================
! Original code written by Hilaire Chevreau.

! Code rewritten by:
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
        subroutine Dassemble_umbrella (natoms, iwrtfpieces, ratom, ftot)
        use umbrella
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iwrtfpieces
        integer, intent (in) :: natoms

        real, intent (in), dimension (3, natoms) :: ratom

! Ouput
        real, intent(inout), dimension (3, natoms) :: ftot      ! total force         
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ix
        integer jx 

        real distance
        real norm 
 
! Procedure
! ===========================================================================
! Initialize 
        do ix = 1, umb_pair
         distance = sqrt((ratom(1,umb_p2(ix)) - ratom(1,umb_p1(ix)))**2      &
     &                   + (ratom(2,umb_p2(ix)) - ratom(2,umb_p1(ix)))**2    &
     &                   + (ratom(3,umb_p2(ix)) - ratom(3,umb_p1(ix)))**2)
         do jx = 1, 3
          fumb(jx,umb_p1(ix)) =                                              &
     &     - umb_CFd(ix)*(distance - umb_d0(ix))                             &
     &                  *(ratom(jx,umb_p1(ix)) - ratom(jx,umb_p2(ix)))/distance 
          fumb(jx,umb_p2(ix)) = - fumb(jx,umb_p1(ix))

          ftot(jx,umb_p1(ix)) = ftot(jx,umb_p1(ix)) + fumb(jx,umb_p1(ix)) 
          ftot(jx,umb_p2(ix)) = ftot(jx,umb_p2(ix)) + fumb(jx,umb_p2(ix)) 
         end do

         if (iwrtfpieces .eq. 1) then 
          write (*,*) ' At A : ', umb_p1(ix), ' at B ', umb_p2(ix), ' Kd ',  &
     &                            umb_CFd(ix)
          write (*,*) ' Umbrella force: '
          write (*,*) fumb(1,umb_p1(ix)), fumb(2,umb_p1(ix)), fumb(3,umb_p1(ix))
          write (*,*) '   '
         end if
        end do
 
! Format Statements
! ===========================================================================
        return
        end
