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
! Motorola, Physical Sciences Research Labs - Jun Wang
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! setgear.f90
! Program Description
! ===========================================================================
! This sets up the gear algorithm constants_fireball
!
! Gear proceedures use the Gear algorithm to move the atoms
!
! Gear algorithm --- use a Taylor series expansion up to fifth-order
!                   in the derivatives.
!
!                   Prediction step:
!                     x(t+dt)=x(t)+dx(t)*dt+(1/2)*ddx(t)*dt**2+.....
!                    dx(t+dt)=dx(t)+ddx(t)*dt+(1/2)*dddx(t)*dt**2+....
!                   ddx(t+dt)=ddx(t)+dddx(t)*dt+(1/2)*ddddx(t)*dt**2+...
!                  dddx(t+dt)=dddx(t)+ddddx(t)*dt+(1/2)*dddddx(t)*dt**2
!                 ddddx(t+dt)=ddddx(t)+dddddx(t)*dt
!                dddddx(t+dt)=dddddx(t)
!
!                    (here dx(t) means first derivative of x wrt t,
!                       ddx(t) means second derivative of x wrt t,etc.)
!
!                   Correction step:
!                      First calculate the force at the predicted
!                        position x(t+dt) from above.  From the force
!                        calculate the acceleration a.  Then take
!                        the difference between the predicted and
!                        actual accelerations:
!
!                                D = a - ddx(t+dt)
!
!                      Now correct all the derivatives:
!                        xnew(t+dt) = x(t+dt)      + c0*D*(0.5*dt**2)
!                       dxnew(t+dt) = dx(t+dt)     + c1*D*(0.5*dt)
!                      ddxnew(t+dt) = ddx(t+dt)    + c2*D
!                     dddxnew(t+dt) = dddx(t+dt)   + c3*D*(3/dt)
!                    ddddxnew(t+dt) = ddddx(t+dt)  + c4*D*(12/dt**2)
!                   dddddxnew(t+dt) = dddddx(t+dt) + c5*D*(60/dt**3)
!
!                     References:
!                        C.W. Gear, ANL Report No. ANL-7126, 1966.
!                        Young Hee Lee, Ph.D. Dissertation, Kent
!                            State University, August 1986.
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
        subroutine setgear
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Procedure
! ===========================================================================
! Gear coeficients for second-order equation (thus cfac(2) = 1.0d0 always)
        cfac = 0.0d0
        write (*,*)
        write (*,*) ' In setgear, using gear order of ', gear_order

        if (gear_order .eq. 2) then   ! velocity verlet style
         cfac(0) = 0.0d0
         cfac(1) = 1.0d0
         cfac(2) = 1.0d0
        else if (gear_order .eq. 3) then
         cfac(0) = 1.0d0/6.0d0
         cfac(1) = 5.0d0/6.0d0
         cfac(2) = 1.0d0
         cfac(3) = 1.0d0/3.0d0
        else if (gear_order .eq. 4) then
         cfac(0) = 19.0d0/120.0d0
         cfac(1) = 3.0d0/4.0d0
         cfac(2) = 1.0d0
         cfac(3) = 1.0d0/2.0d0
         cfac(4) = 1.0d0/12.0d0
        else if (gear_order .eq. 5) then
         cfac(0) = 3.0d0/20.0d0
         cfac(1) = 251.0d0/360.0d0
         cfac(2) = 1.0d0
         cfac(3) = 11.0d0/18.0d0
         cfac(4) = 1.0d0/6.0d0
         cfac(5) = 1.0d0/60.0d0
        else if (gear_order .eq. 6) then
         cfac(0) = 863.0d0/6048.0d0
         cfac(1) = 665.0d0/1008.0d0
         cfac(2) = 1.0d0
         cfac(3) = 25.0d0/36.0d0
         cfac(4) = 35.0d0/144.0d0
         cfac(5) = 1.0d0/24.0d0
!        cfac(6) = 1.0d0/360.0d0
        else if (gear_order .eq. 7) then
         cfac(0) = 1925.0d0/14112.0d0
         cfac(1) = 19087.0d0/30240.0d0
         cfac(2) = 1.0d0
         cfac(3) = 137.0d0/180.0d0
         cfac(4) = 5.0d0/16.0d0
         cfac(5) = 17.0d0/240.0d0
!        cfac(6) = 1.0d0/120.0d0
!        cfac(7) = 1.0d0/2520.0d0
        else
         write (*,*) ' in setgear, gear_order = ', gear_order
         write (*,*) ' only gears of 2 through 7 are allowed'
         stop
        endif
        write (*,*)

! Format Statements
! ===========================================================================
 
        return
        end
