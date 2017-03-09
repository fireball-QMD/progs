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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! DgelementsG.f90
! Program Description
! ===========================================================================
! This program gets derivatives of three-center matrix elements.
! Very similar to gauss_matel but with derivatives
! ---> See gauss_matel subroutine for more details.

! ===========================================================================
! Code rewritten by:
! Hao Wang
! Department of Physics and Astronomy
! Brigham Young University
! N311 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-9270
! ===========================================================================
!
! Program Declaration
! =====================================================================
        subroutine Dgelements_VXC(l1,l2,ish,jsh,imu,jmu,in1,in2,indna,  &
     &                            x,y,cost,rmatel,iparam)
! =====================================================================

        use dimensions
        use gaussG
        use constants_fireball
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna

        integer, intent (in) :: l1
        integer, intent (in) :: l2

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y

        integer ish
        integer jsh
        integer imu
        integer jmu

        integer iparam 

! ===========================================================================
! output
        real rmatel

! Local Parameters and Data Declaration
! ===========================================================================
        integer isorp
        integer ial1
        integer ial2
        integer ial3

        real rmu
        real lilm
        real bigM
        real d12
        real dcm3
! d12 is the distance between centers 1 and 2
! dcm3 is the distance from the "center of mass" of centers
! 1 and 2 to center 3.
        real bigQ
        real gM
        real xcm
        real ycm
        real zcm

        real al1
        real al2
        real al3

        real coef1
        real coef2
        real coef3

        real dQ
        real factor
        real dfactor

        real c2
        real c4

        real x1
        real y1
        real z1

        real x2
        real y2
        real z2

        real x3
        real y3
        real z3

! derivatives of xcm,ycm,zcm,x,y,and z...
        real dxcm
        real dycm
        real dzcm
 
        real dx1
        real dy1
        real dz1

        real dx2
        real dy2
        real dz2

        real dx3
        real dy3
        real dz3

! Proceedure
! =====================================================================
! Coordinates and Arrangement...
!
!
!           X  ^                                * (NA)
!              !                             +
! X-Z PLANE    !                          +
!              !                      +   THETA               -----> Z
!              O ---------------------------------------------O
!            1 END OF BC                                OTHER END OF BC
!               MU                                          NU

        d12=y

        x1 = 0.
        y1 = 0.
        z1 = -y/2.
        x2 = 0.
        y2 = 0.
        z2 = y/2.
        x3 = x*sqrt(abs(1. - cost**2))
        y3 = 0.
        z3 = x*cost

! Derivatives of coordinates...
! The iparam says whether the derivatives
! are w/respect to the y (the distance between centers 1 and 2) and
! x3 or z3 (coordinates of na potential...)
!    iparam = 1 ---> y
!    iparam = 2 ---> x3
!    iparam = 3 ---> z3

! Initialize Derivatives
        dx1=0.
        dy1=0.
        dz1=0.
        
        dx2=0.
        dy2=0.
        dz2=0.
        
        dx3=0.
        dy3=0.
        dz3=0.
        
        if (iparam.eq.1) then
           dz1=-1./2.
           dz2=1./2.
        else if (iparam.eq.2) then
           dx3=1.
        else if (iparam.eq.3) then
           dz3=1.
        end if

        rmatel=0.

! ---------------------------------------------------------------------
! MHL (Sep. 29. 2004)
! VNA doesn't have shell information
!        real gcoefficientsVNA(max_alphas,nspec_max)
!        real alphaVNA(max_alphas,nspec_max)
!        integer nalphaVNA(nspec_max)

! Electron density (0: total, 1: first shell, etc.)
!        real gcoefficientsN(max_alphas,0:nsh_max,nspec_max)
!        real alphaN(max_alphas,0:nsh_max,nspec_max)
!        integer nalphaN(0:nsh_max,nspec_max)

! Wavefunction (begins from 1 - meaning first shell)
!        real gcoefficientsPSI(max_alphas,1:nsh_max,nspec_max)
!        real alphaPSI(max_alphas,1:nsh_max,nspec_max)
!        integer nalphaPSI(1:nsh_max,nspec_max)

! Wavefunction for shperical approximation (begins from 1)
!        real gcoefficientsPSIS(max_alphas,1:nsh_max,nspec_max)
!        real alphaPSIS(max_alphas,1:nsh_max,nspec_max)
!        integer nalphaPSIS(1:nsh_max,nspec_max)

! ----------------------------------------------------------------------
! alpha1 (wavefunction)
        do ial1=1,nalphaPSI(ish,in1)
           al1=alphaPSI(ial1,ish,in1)
           coef1=gcoefficientsPSI(ial1,ish,in1)

! alpha2 (wavefunction)
        do ial2=1,nalphaPSI(jsh,in2)
        al2=alphaPSI(ial2,jsh,in2)
        coef2=gcoefficientsPSI(ial2,jsh,in2)

        lilm=al1+al2
        rmu=al1*al2/lilm
        dcm3=(z3 - (al1*z1 + al2*z2)/lilm)**2 + x3**2
        dcm3=sqrt(abs(dcm3))

! alpha3
        isorp = 0
        do ial3=1,nalphaN(isorp,indna)
        al3=alphaN(ial3,isorp,indna)
        coef3=gcoefficientsN(ial3,isorp,indna)

           bigM=al1+al2+al3

! The "center of mass" is
!   (al1*r1+al2*r2+al3*r3)/(al1+al2+al3)
! This is location of the gaussian that you get when
! the three gaussians are multiplied together.
        xcm=x3*al3/bigM
        ycm=0.
        zcm=(al1*z1 + al2*z2 + al3*z3)/bigM

! initialize derivatives of coordinates...
           dxcm=0.
           dycm=0.
           dzcm=0.

           if (iparam.eq.1) then
           dzcm=(al1*dz1 + al2*dz2 + al3*dz3)/bigM
           else if (iparam.eq.2) then
           dxcm=al3/bigM
           else if (iparam.eq.3) then
           dzcm=al3/bigM
           end if

! bigQ is a factor that you get in front of the product gaussian...
           bigQ=exp( -rmu*d12**2 - lilm*al3*dcm3**2/bigM)

! gM = integral from -inf to +inf of exp^(-M*x**2)
           gM=sqrt(pi/bigM)

! dcm3**2 = (z3 - (al1*z1 + al2*z2)/lilm)**2 + x3**2

! Derivative of bigQ...
           if (iparam .eq. 1) then
            dQ= - 2.*rmu*y - lilm*al3*(2.*(z3 - (al1*z1 + al2*z2)/lilm)*&
     &            (dz3 - (al1*dz1 + al2*dz2)/lilm) + 2.*x3*dx3)/bigM
            else if (iparam .eq. 2) then
            dQ= - 2.*lilm*al3*(x3)/bigM
            else if (iparam .eq. 3) then
            dQ= - 2.*lilm*al3*(z3 - (al1*z1 + al2*z2)/lilm)/bigM
            end if
            dQ=dQ*bigQ

! Arrangement of matrix elements in fireball format
!
!           orbital     S        l= 0,   m=  0
!
!           orbital     PY       l= 1,   m= -1
!           orbital     PZ       l= 1,   m=  0
!           orbital     PX       l= 1,   m=  1
!
!           orbital     xy       l= 1,   m= -2
!           orbital     yz       l= 1,   m= -1
!           orbital  3z^2-r^2    l= 1,   m=  0
!           orbital     xz       l= 1,   m=  1
!           orbital   x^2-y^2    l= 1,   m=  2

! Examples:
! p_z state centered on atom 2 is basically:
!    (z-z2) = (z-zcm)+(zcm-z2)
!
! d_xy on 1st center --->
! (x-x1)*(y-y1) = (x-xcm)*(y-ycm) + (ycm-y1)*(x-xcm) +
!                                   (xcm-x1)*(y-ycm) + (ycm-y1)*(xcm-x1)

! A lot of the coordinates below are zero:
!    y1,y2,ycm = x1,x2 = z1 = 0.

! ==========================
!        get factor
! ==========================
! This is what you get when you multiply the two polynomials (states) together
! and integrate against the product gaussian...

          factor = 0.d0
         dfactor = 0.d0

! It's an extra factor that comes up if you integrate
! x**2 e^(-Mx**2) versus e^(-M*x**2).
        c2=1./2./bigM
! or x**4 e^(-Mx**2)...
        c4=3./4./bigM**2

! SS
        if (l1 .eq. 0 .and. l2 .eq. 0) then
           factor = 1.d0
          dfactor = 0.d0
        end if
! SP
        if (l1 .eq. 0 .and. l2 .eq. 1) then
          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq.  0) then
           factor = zcm-z2
          dfactor = dzcm-dz2
          end if
          if (imu .eq. 0 .and. jmu .eq.  1) then
           factor = xcm-x2
          dfactor = dxcm-dx2
          end if
        end if ! SP

! PS
        if (l1 .eq. 1 .and. l2 .eq. 0) then
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq.  0 .and. jmu .eq. 0) then
           factor = zcm-z1
          dfactor = dzcm-dz1
          end if
          if (imu .eq.  1 .and. jmu .eq. 0) then
           factor = xcm-x1
          dfactor = dxcm-dx1
          end if
        end if ! PS

! PP
        if (l1 .eq. 1 .and. l2 .eq. 1) then
          if (imu .eq. -1 .and. jmu .eq. -1) then
           factor = c2
          dfactor = 0.d0
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if

          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = c2 + (zcm-z1)*(zcm-z2)
          dfactor = (dzcm-dz1)*(zcm-z2) + (zcm-z1)*(dzcm-dz2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
           factor = (zcm-z1)*(xcm-x2)
          dfactor = (dzcm-dz1)*(xcm-x2) + (zcm-z1)*(dxcm-dx2)
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
           factor = (xcm-x1)*(zcm-z2)
          dfactor = (dxcm-dx1)*(zcm-z2) + (xcm-x1)*(dzcm-dz2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
           factor = c2 + (xcm-x1)*(xcm-x2)
          dfactor = (dxcm-dx1)*(xcm-x2) + (xcm-x1)*(dxcm-dx2)
          end if

        end if ! PP

! SD
        if (l1 .eq. 0 .and. l2 .eq. 2) then
          if (imu .eq. 0 .and. jmu .eq. -2) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = 2.*(zcm-z2)**2-(xcm-x2)**2
          dfactor = 4.*(zcm-z2)*(dzcm-dz2) - 2.*(xcm-x2)*(dxcm-dx2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
           factor = (xcm-x2)*(zcm-z2)
          dfactor = (dxcm-dx2)*(zcm-z2) + (xcm-x2)*(dzcm-dz2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
           factor = (xcm-x2)**2
          dfactor = 2.*(xcm-x2)*(dxcm-dx2)
          end if
        end if ! SD

! DS
        if (l1 .eq. 2 .and. l2 .eq. 0) then
          if (imu .eq. -2 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = 2.*(zcm-z1)**2-(xcm-x1)**2
          dfactor = 4.*(zcm-z1)*(dzcm-dz1) - 2.*(xcm-x1)*(dxcm-dx1)
          end if
          if (imu .eq. 1 .and. jmu .eq. 0) then
           factor = (xcm-x1)*(zcm-z1)
          dfactor = (dxcm-dx1)*(zcm-z1) + (xcm-x1)*(dzcm-dz1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
           factor = (xcm-x1)**2
          dfactor = 2.*(xcm-x1)*(dxcm-dx1)
          end if
        end if ! DS

! PD
        if (l1 .eq. 1 .and. l2 .eq. 2) then
          if (imu .eq. -1 .and. jmu .eq. -2) then
           factor = c2*(xcm-x2)
          dfactor = c2*(dxcm-dx2)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
           factor = c2*(zcm-z2)
          dfactor = c2*(dzcm-dz2)
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if


          if (imu .eq. 0 .and. jmu .eq. -2) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = 4.*c2*(zcm-z2)+(zcm-z1)*(2.*(zcm-z2)**2-(xcm-x2)**2)
          dfactor = 4.*c2*(dzcm-dz2) +                                  &
     &     (dzcm-dz1)*(2.*(zcm-z2)**2-(xcm-x2)**2) +                    &
     &       (zcm-z1)*(4.*(zcm-z2)*(dzcm-dz2) - 2.*(xcm-x2)*(dxcm-dx2))
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
           factor = c2*(xcm-x2) + (zcm-z1)*(xcm-x2)*(zcm-z2)
          dfactor = c2*(dxcm-dx2) + (dzcm-dz1)*(xcm-x2)*(zcm-z2) +      &
     &                              (zcm-z1)*(dxcm-dx2)*(zcm-z2) +      &
     &                              (zcm-z1)*(xcm-x2)*(dzcm-dz2) 
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
           factor = (zcm-z1)*(xcm-x2)**2
          dfactor = (dzcm-dz1)*(xcm-x2)**2 +                            &
     &              2.*(zcm-z1)*(xcm-x2)*(dxcm-dx2)
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
           factor =-2.*c2*(xcm-x2)+(xcm-x1)*(2.*(zcm-z2)**2-(xcm-x2)**2)
          dfactor =-2.*c2*(dxcm-dx2) +                                  &
     &    (dxcm-dx1)*(2.*(zcm-z2)**2 - (xcm-x2)**2) +                   &
     &      (xcm-x1)*(4.*(zcm-z2)*(dzcm-dz2) - 2.*(xcm-x2)*(dxcm-dx2))
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
           factor = c2*(zcm-z2) + (xcm-x1)*(xcm-x2)*(zcm-z2)
          dfactor = c2*(dzcm-dz2) + (dxcm-dx1)*(xcm-x2)*(zcm-z2) +      &
     &                              (xcm-x1)*(dxcm-dx2)*(zcm-z2) +      &
     &                              (xcm-x1)*(xcm-x2)*(dzcm-dz2) 
          end if
          if (imu .eq. 1 .and. jmu .eq. 2) then
           factor = 2.*c2*(xcm-x2) + (xcm-x1)*(xcm-x2)**2
          dfactor = 2.*c2*(dxcm-dx2) + (dxcm-dx1)*(xcm-x2)**2 +         &
     &                              2.*(xcm-x1)*(xcm-x2)*(dxcm-dx2)
          end if

        end if ! PD

! DP
        if (l1 .eq. 2 .and. l2 .eq. 1) then
          if (imu .eq. -2 .and. jmu .eq. -1) then
           factor = c2*(xcm-x1)
          dfactor = c2*(dxcm-dx1)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
           factor = c2*(zcm-z1)
          dfactor = c2*(dzcm-dz1)
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.d0
          dfactor = 0.d0
          end if


          if (imu .eq. -2 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.d0
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = 4.*c2*(zcm-z1)+(zcm-z2)*(2.*(zcm-z1)**2-(xcm-x1)**2)
          dfactor = 4.*c2*(dzcm-dz1) +                                  &
     &    (dzcm-dz2)*(2.*(zcm-z1)**2-(xcm-x1)**2) +                     &
     &    (zcm-z2)*(4.*(zcm-z1)*(dzcm-dz1)-2.*(xcm-x1)*(dxcm-dx1))
          end if
          if (imu .eq. 1 .and. jmu .eq. 0) then
           factor = c2*(xcm-x1) + (zcm-z2)*(xcm-x1)*(zcm-z1)
          dfactor = c2*(dxcm-dx1) + (dzcm-dz2)*(xcm-x1)*(zcm-z1) +      &
     &                              (zcm-z2)*(dxcm-dx1)*(zcm-z1) +      &
     &                              (zcm-z2)*(xcm-x1)*(dzcm-dz1)  
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
           factor = (zcm-z2)*(xcm-x1)**2
          dfactor = (dzcm-dz2)*(xcm-x1)**2 +                            &
     &                            2.*(zcm-z2)*(xcm-x1)*(dxcm-dx1)
          end if

          if (imu .eq. 0 .and. jmu .eq. 1) then
           factor =-2.*c2*(xcm-x1)+(xcm-x2)*(2.*(zcm-z1)**2-(xcm-x1)**2)
          dfactor = -2.*c2*(dxcm-dx1) +                                 &
     &    (dxcm-dx2)*(2.*(zcm-z1)**2-(xcm-x1)**2) +                     &
     &      (xcm-x2)*(4.*(zcm-z1)*(dzcm-dz1)-2.*(xcm-x1)*(dxcm-dx1))
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
           factor = c2*(zcm-z1) + (xcm-x2)*(xcm-x1)*(zcm-z1)
          dfactor = c2*(dzcm-dz1) + (dxcm-dx2)*(xcm-x1)*(zcm-z1) +      &
     &                              (xcm-x2)*(dxcm-dx1)*(zcm-z1) +      &
     &                              (xcm-x2)*(xcm-x1)*(dzcm-dz1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 1) then
           factor = 2.*c2*(xcm-x1) + (xcm-x2)*(xcm-x1)**2
          dfactor = 2.*c2*(dxcm-dx1) + (dxcm-dx2)*(xcm-x1)**2 +         &
     &                                2.*(xcm-x2)*(xcm-x1)*(dxcm-dx1)
          end if

        end if ! DP


! DD
        if (l1 .eq. 2 .and. l2 .eq. 2) then

          if (imu .eq. -2 .and. jmu .eq. -2) then
           factor = c2*c2 + c2*(xcm-x1)*(xcm-x2)
          dfactor = c2*(dxcm-dx1)*(xcm-x2) +  c2*(xcm-x1)*(dxcm-dx2)
          end if
          if (imu .eq. -2 .and. jmu .eq. -1) then
           factor = c2*(xcm-x1)*(zcm-z2)
          dfactor = c2*(dxcm-dx1)*(zcm-z2) + c2*(xcm-x1)*(dzcm-dz2)
          end if
          if (imu .eq. -2 .and. jmu .eq. 0) then
           factor = 0.
          dfactor = 0.d0
          end if

          if (imu .eq. -1 .and. jmu .eq. -2) then
           factor = c2*(xcm-x2)*(zcm-z1)
          dfactor = c2*(dxcm-dx2)*(zcm-z1) + c2*(xcm-x2)*(dzcm-dz1)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
           factor = c2*c2 + c2*(zcm-z1)*(zcm-z2)
          dfactor = c2*(dzcm-dz1)*(zcm-z2) + c2*(zcm-z1)*(dzcm-dz2)
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
           factor = 0.
          dfactor = 0.d0
          end if

          if (imu .eq. 0 .and. jmu .eq. -2) then
           factor = 0.
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
           factor = 0.
          dfactor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
           factor = (2.*(zcm-z1)**2 - (xcm-x1)**2)*                     &
     &              (2.*(zcm-z2)**2 - (xcm-x2)**2) +                    &
     &              16.*c2*(zcm-z1)*(zcm-z2) -                          &
     &              4.*c2*(xcm-x1)*(xcm-x2) +6.*c4 - 6.*c2*c2
          dfactor = (4.*(zcm-z1)*(dzcm-dz1) - 2.*(xcm-x1)*(dxcm-dx1))*  &
     &              (2.*(zcm-z2)**2 - (xcm-x2)**2) +                    &
     &              (4.*(zcm-z2)*(dzcm-dz2) - 2.*(xcm-x2)*(dxcm-dx2))*  &
     &              (2.*(zcm-z1)**2 - (xcm-x1)**2) +                    &
     &              16.*c2*(dzcm-dz1)*(zcm-z2) +                        &
     &              16.*c2*(zcm-z1)*(dzcm-dz2) -                        &
     &               4.*c2*(dxcm-dx1)*(xcm-x2) -                        &
     &               4.*c2*(xcm-x1)*(dxcm-dx2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
           factor = (xcm-x2)*(zcm-z2)*(2.*(zcm-z1)**2 - (xcm-x1)**2) +  &
     &              4.*c2*(xcm-x2)*(zcm-z1) - 2.*c2*(zcm-z2)*(xcm-x1)
          dfactor = (dxcm-dx2)*(zcm-z2)*(2.*(zcm-z1)**2 - (xcm-x1)**2) + &
     &             (xcm-x2)*(dzcm-dz2)*(2.*(zcm-z1)**2 - (xcm-x1)**2) + &
     &             (xcm-x2)*(zcm-z2)*(4.*(zcm-z1)*(dzcm-dz1) -          &
     &             2.*(xcm-x1)*(dxcm-dx1)) +                            &
     &             4.*c2*(dxcm-dx2)*(zcm-z1) + 4.*c2*(xcm-x2)*(dzcm-dz1)&
     &           - 2.*c2*(dzcm-dz2)*(xcm-x1) - 2.*c2*(zcm-z2)*(dxcm-dx1)
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
           factor = (2.*(zcm-z1)**2 - (xcm-x1)**2)*(xcm-x2)**2 -        &
     &              4.*c2*(xcm-x1)*(xcm-x2)
          dfactor = (4.*(zcm-z1)*(dzcm-dz1) - 2.*(xcm-x1)*(dxcm-dx1))*  &
     &              (xcm-x2)**2 + (2.*(zcm-z1)**2 - (xcm-x1)**2)*       &
     &              2.*(xcm-x2)*(dxcm-dx2) -                            &
     &              4.*c2*(dxcm-dx1)*(xcm-x2) -                         &
     &              4.*c2*(xcm-x1)*(dxcm-dx2)
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
           factor = (xcm-x1)*(zcm-z1)*(2.*(zcm-z2)**2 - (xcm-x2)**2) +  &
     &              4.*c2*(xcm-x1)*(zcm-z2) - 2.*c2*(zcm-z1)*(xcm-x2)
          dfactor = (dxcm-dx1)*(zcm-z1)*(2.*(zcm-z2)**2 - (xcm-x2)**2) +&
     &              (xcm-x1)*(dzcm-dz1)*(2.*(zcm-z2)**2 - (xcm-x2)**2) +&
     &              (xcm-x1)*(zcm-z1)*(4.*(zcm-z2)*(dzcm-dz2) -         &
     &              2.*(xcm-x2)*(dxcm-dx2)) +                           &
     &              4.*c2*(dxcm-dx1)*(zcm-z2) +                         &
     &              4.*c2*(xcm-x1)*(dzcm-dz2) -                         &
     &              2.*c2*(dzcm-dz1)*(xcm-x2) -                         &
     &              2.*c2*(zcm-z1)*(dxcm-dx2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
           factor = (xcm-x1)*(zcm-z1)*(xcm-x2)*(zcm-z2) + c2*c2 +       &
     &              c2*(zcm-z1)*(zcm-z2) + c2*(xcm-x1)*(xcm-x2)
          dfactor = (dxcm-dx1)*(zcm-z1)*(xcm-x2)*(zcm-z2) +             &
     &              (xcm-x1)*(dzcm-dz1)*(xcm-x2)*(zcm-z2) +             &
     &              (xcm-x1)*(zcm-z1)*(dxcm-dx2)*(zcm-z2) +             &
     &              (xcm-x1)*(zcm-z1)*(xcm-x2)*(dzcm-dz2) +             &
     &              c2*(dzcm-dz1)*(zcm-z2) + c2*(zcm-z1)*(dzcm-dz2) +   &
     &              c2*(dxcm-dx1)*(xcm-x2) + c2*(xcm-x1)*(dxcm-dx2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 2) then
           factor = (xcm-x1)*(zcm-z1)*(xcm-x2)**2 +                     &
     &              2.*c2*(zcm-z1)*(xcm-x2)
          dfactor = (dxcm-dx1)*(zcm-z1)*(xcm-x2)**2 +                   &
     &              (xcm-x1)*(dzcm-dz1)*(xcm-x2)**2 +                   &
     &              2.*(xcm-x1)*(zcm-z1)*(xcm-x2)*(dxcm-dx2) +          &
     &              2.*c2*(dzcm-dz1)*(xcm-x2) +                         &
     &              2.*c2*(zcm-z1)*(dxcm-dx2)
          end if

          if (imu .eq. 2 .and. jmu .eq. 2) then
           factor = (xcm-x1)**2*(xcm-x2)**2 + 2.*c4 -2.*c2*c2 +         &
     &              4.*c2*(xcm-x1)*(xcm-x2)
          dfactor = 2.*(xcm-x1)*(dxcm-dx1)*(xcm-x2)**2 +                &
     &              2.*(xcm-x1)**2*(xcm-x2)*(dxcm-dx2) +                &
     &              4.*c2*(dxcm-dx1)*(xcm-x2) + 4.*c2*(xcm-x1)*(dxcm-dx2)
          end if
          if (imu .eq. 2 .and. jmu .eq. 1) then
           factor = (xcm-x2)*(zcm-z2)*(xcm-x1)**2 +                     &
     &             2.*c2*(zcm-z2)*(xcm-x1)
          dfactor = (dxcm-dx2)*(zcm-z2)*(xcm-x1)**2 +                    &
     &             (xcm-x2)*(dzcm-dz2)*(xcm-x1)**2 +                    &
     &             (xcm-x2)*(zcm-z2)*2.*(xcm-x1)*(dxcm-dx1) +           &
     &             2.*c2*(dzcm-dz2)*(xcm-x1) + 2.*c2*(zcm-z2)*(dxcm-dx1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
           factor = (2.*(zcm-z2)**2 - (xcm-x2)**2)*(xcm-x1)**2 -        &
     &              4.*c2*(xcm-x2)*(xcm-x1)
          dfactor =                                                     &
     &    (4.*(zcm-z2)*(dzcm-dz2) - 2.*(xcm-x2)*(dxcm-dx2))*(xcm-x1)**2 &
     &    + (2.*(zcm-z2)**2 - (xcm-x2)**2)*2.*(xcm-x1)*(dxcm-dx1) -     &
     &      4.*c2*(dxcm-dx2)*(xcm-x1) - 4.*c2*(xcm-x2)*(dxcm-dx1)
          end if

        end if ! end if of DD

        rmatel = rmatel + coef1*coef2*coef3*gM**3*                      &
     &                  (bigQ*dfactor + factor*dQ)
        end do
        end do
        end do

        rmatel = rmatel*gfactor(l1,imu)*gfactor(l2,jmu)

        return
        end
! =======================================================================
