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


! DgelementsGS.f90
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
! ===========================================================================
       subroutine DgelementsGS_VXC (l1,l2,ish,jsh,im,jm,in1,in2,indna,  &
     &                              x,y,cost,rmatel,iparam,isorp)
! ===========================================================================
!
        use dimensions
        use interactions
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
        integer im
        integer jm

        integer iparam 
!MHL
        integer, intent (in) :: isorp
!        integer, intent (in) :: interaction

! ===========================================================================
! output
        real rmatel

! Local Parameters and Data Declaration
! ===========================================================================
! MHL
!        integer malpha

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
           dz1 = -1./2.
           dz2 = 1./2.
        else if (iparam.eq.2) then
           dx3 = 1.
        else if (iparam.eq.3) then
           dz3 = 1.
        end if

        rmatel=0.d0

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
! alpha1 (spherical approximation, nalphaPSI -> nalphaPSIS)
        do ial1=1,nalphaPSIS(ish,in1)
        al1=alphaPSIS(ial1,ish,in1)
        coef1=gcoefficientsPSIS(ial1,ish,in1)

! alpha2 (spherical approximation, nalphaPSI -> nalphaPSIS)
        do ial2=1,nalphaPSIS(jsh,in2)
        al2=alphaPSIS(ial2,jsh,in2)
        coef2=gcoefficientsPSIS(ial2,jsh,in2)

        lilm=al1+al2
        rmu=al1*al2/lilm
        dcm3=(z3 - (al1*z1 + al2*z2)/lilm)**2 + x3**2
        dcm3=sqrt(abs(dcm3))

! alpha3
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

! derivatives of coordinates...
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

! BigQ is a factor that you get in front of the product gaussian...
        bigQ=exp(-rmu*d12**2-lilm*al3*dcm3**2/bigM)

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
! p_z state centered on atom 2 is
!    (z-z2) = (z-zcm)+(zcm-z2)
!
! d_xy on 1st center is
!    (x-x1)*(y-y1) = (x-xcm)*(y-ycm) + (ycm-y1)*(x-xcm) +
!                                      (xcm-x1)*(y-ycm) + (ycm-y1)*(xcm-x1)

! A lot of the coordinates below are zero:
!    y1,y2,ycm = x1,x2 = z1 = 0.

! ==========================
!        get factor
! ==========================
! This is what you get when you multiply the two polynomials (states) together
! and integrate against the product gaussian...
        factor=0.
        dfactor=0.

! For the spherical approxiamtion, only l1 = 0 and l2 = 0 considered,
! so m1 = 0 and m2 =0.
! IT IS ONLY ONE TERM NEED TO BE EVALUATED
         factor = 1.d0
        dfactor = 0.d0

        rmatel = rmatel + coef1*coef2*coef3*gM**3*                      &
     &                (bigQ*dfactor + factor*dQ)
        end do
        end do
        end do

        rmatel = rmatel*gfactor(0,0)*gfactor(0,0)

        return
        end
! =======================================================================
