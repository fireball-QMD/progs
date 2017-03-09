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

 
! gelementsG.f90
! Program Description
! ===========================================================================
! This routine gets the matrix element - <psi|V|psi>, where V is the potential
! of the neutral atom or exchange-correlation.
! Example: l1 = 1, imu = 1 and l2 = 1, jmu = 0 means get <p_x|Vna|p_z>.
! (p_x on center 1 and p_z on center 2...)
!
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
        subroutine gelementsG_overlap (in1, in2, l1,l2, issh, jssh,     &
     &                                 imu, jmu, y, rmatel)
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
        integer, intent (in) :: imu
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: issh
        integer, intent (in) :: jmu
        integer, intent (in) :: jssh
        integer, intent (in) :: l1
        integer, intent (in) :: l2
 
        real, intent (in) :: y

! Output
        real, intent (out) :: rmatel

! Local Parameters and Data Declaration
! ===========================================================================
        integer ial1
        integer ial2

        real rmu
        real lilm
        real bigQ
        real d12
        real gM

        real al1
        real al2

        real coef1
        real coef2
        real factor

        real c2
        real c4

        real x1
        real y1
        real z1

        real x2
        real y2
        real z2

        real xcm
        real ycm
        real zcm

! Proceedure
! ===========================================================================

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

        rmatel=0.d0
! MHL 
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
        do ial1=1,nalphaPSI(issh,in1)
        al1=alphaPSI(ial1,issh,in1)
        coef1=gcoefficientsPSI(ial1,issh,in1)

! alpha2 (wavefunction)
        do ial2=1,nalphaPSI(jssh,in2)
        al2=alphaPSI(ial2,jssh,in2)
        coef2=gcoefficientsPSI(ial2,jssh,in2)

        lilm = al1+al2
         rmu = al1*al2/lilm

        xcm = 0.d0
        ycm = 0.d0
        zcm = (al1*z1 + al2*z2)/lilm

! bigQ is a factor that you get in front of the product gaussian...
        bigQ = exp(-rmu*d12**2)

! gM = integral from -inf to +inf of exp^(-M*x**2)
        gM=sqrt(pi/lilm)

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

        factor=0.

! It's an extra factor that comes up if you integrate
! x**2 e^(-Mx**2) versus e^(-M*x**2).
        c2=1./2./lilm
! or x**4 e^(-Mx**2)...
        c4=3./4./lilm**2

! SS
        if (l1 .eq. 0 .and. l2 .eq. 0) then
          factor = 1.d0
        end if
! SP
        if (l1 .eq. 0 .and. l2 .eq. 1) then
          if (imu .eq. 0 .and. jmu .eq. -1) then 
          factor = 0.d0
          end if 
          if (imu .eq. 0 .and. jmu .eq.  0) then
          factor = zcm-z2
          end if
          if (imu .eq. 0 .and. jmu .eq.  1) then
          factor = xcm-x2
          end if
        end if ! SP

! PS
        if (l1 .eq. 1 .and. l2 .eq. 0) then
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.d0
          end if
          if (imu .eq.  0 .and. jmu .eq. 0) then
          factor = zcm-z1
          end if
          if (imu .eq.  1 .and. jmu .eq. 0) then
          factor = xcm-x1
          end if
        end if ! PS

! PP    
        if (l1 .eq. 1 .and. l2 .eq. 1) then
          if (imu .eq. -1 .and. jmu .eq. -1) then
          factor = c2
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.d0
          end if

          if (imu .eq. 0 .and. jmu .eq. -1) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = c2 + (zcm-z1)*(zcm-z2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
          factor = (zcm-z1)*(xcm-x2)
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
          factor = (xcm-x1)*(zcm-z2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
          factor = c2 + (xcm-x1)*(xcm-x2)
          end if

        end if ! PP
           
! SD
        if (l1 .eq. 0 .and. l2 .eq. 2) then
          if (imu .eq. 0 .and. jmu .eq. -2) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = 2.*(zcm-z2)**2-(xcm-x2)**2
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
          factor = (xcm-x2)*(zcm-z2) 
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
          factor = (xcm-x2)**2
          end if
        end if ! SD

! DS
        if (l1 .eq. 2 .and. l2 .eq. 0) then
          if (imu .eq. -2 .and. jmu .eq. 0) then
          factor = 0.d0
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = 2.*(zcm-z1)**2-(xcm-x1)**2
          end if
          if (imu .eq. 1 .and. jmu .eq. 0) then
          factor = (xcm-x1)*(zcm-z1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
          factor = (xcm-x1)**2
          end if
        end if ! DS

! PD
        if (l1 .eq. 1 .and. l2 .eq. 2) then
          if (imu .eq. -1 .and. jmu .eq. -2) then
          factor = c2*(xcm-x2)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
          factor = c2*(zcm-z2)
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.d0
          end if


          if (imu .eq. 0 .and. jmu .eq. -2) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = 4.*c2*(zcm-z2)+(zcm-z1)*(2.*(zcm-z2)**2-(xcm-x2)**2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
          factor = c2*(xcm-x2) + (zcm-z1)*(xcm-x2)*(zcm-z2)
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
          factor = (zcm-z1)*(xcm-x2)**2
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
          factor =-2.*c2*(xcm-x2)+(xcm-x1)*(2.*(zcm-z2)**2-(xcm-x2)**2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
          factor = c2*(zcm-z2) + (xcm-x1)*(xcm-x2)*(zcm-z2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 2) then
          factor = 2.*c2*(xcm-x2) + (xcm-x1)*(xcm-x2)**2
          end if

        end if ! PD

! DP
        if (l1 .eq. 2 .and. l2 .eq. 1) then
          if (imu .eq. -2 .and. jmu .eq. -1) then
          factor = c2*(xcm-x1)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
          factor = c2*(zcm-z1)
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
          factor = 0.d0
          end if


          if (imu .eq. -2 .and. jmu .eq. 0) then
          factor = 0.d0
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.d0
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = 4.*c2*(zcm-z1)+(zcm-z2)*(2.*(zcm-z1)**2-(xcm-x1)**2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 0) then
          factor = c2*(xcm-x1) + (zcm-z2)*(xcm-x1)*(zcm-z1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
          factor = (zcm-z2)*(xcm-x1)**2
          end if

          if (imu .eq. 0 .and. jmu .eq. 1) then
          factor =-2.*c2*(xcm-x1)+(xcm-x2)*(2.*(zcm-z1)**2-(xcm-x1)**2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
          factor = c2*(zcm-z1) + (xcm-x2)*(xcm-x1)*(zcm-z1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 1) then
          factor = 2.*c2*(xcm-x1) + (xcm-x2)*(xcm-x1)**2
          end if

        end if ! DP

! DD
        if (l1 .eq. 2 .and. l2 .eq. 2) then

          if (imu .eq. -2 .and. jmu .eq. -2) then
          factor = c2*c2 + c2*(xcm-x1)*(xcm-x2)
          end if
          if (imu .eq. -2 .and. jmu .eq. -1) then
          factor = c2*(xcm-x1)*(zcm-z2)
          end if
          if (imu .eq. -2 .and. jmu .eq. 0) then
          factor = 0.
          end if

          if (imu .eq. -1 .and. jmu .eq. -2) then
          factor = c2*(xcm-x2)*(zcm-z1)
          end if
          if (imu .eq. -1 .and. jmu .eq. -1) then
          factor = c2*c2 + c2*(zcm-z1)*(zcm-z2)
          end if
          if (imu .eq. -1 .and. jmu .eq. 0) then
          factor = 0.
          end if

          if (imu .eq. 0 .and. jmu .eq. -2) then
          factor = 0.
          end if
          if (imu .eq. 0 .and. jmu .eq. -1) then
          factor = 0.
          end if
          if (imu .eq. 0 .and. jmu .eq. 0) then
          factor = (2.*(zcm-z1)**2 - (xcm-x1)**2)*                      &
     &             (2.*(zcm-z2)**2 - (xcm-x2)**2) +                     &
     &             16.*c2*(zcm-z1)*(zcm-z2) -                           &
     &             4.*c2*(xcm-x1)*(xcm-x2) +6.*c4 - 6.*c2*c2
          end if
          if (imu .eq. 0 .and. jmu .eq. 1) then
          factor = (xcm-x2)*(zcm-z2)*(2.*(zcm-z1)**2 - (xcm-x1)**2) +   &
     &             4.*c2*(xcm-x2)*(zcm-z1) - 2.*c2*(zcm-z2)*(xcm-x1)
          end if
          if (imu .eq. 0 .and. jmu .eq. 2) then
          factor = (2.*(zcm-z1)**2 - (xcm-x1)**2)*(xcm-x2)**2 -         &
     &             4.*c2*(xcm-x1)*(xcm-x2)
          end if

          if (imu .eq. 1 .and. jmu .eq. 0) then
          factor = (xcm-x1)*(zcm-z1)*(2.*(zcm-z2)**2 - (xcm-x2)**2) +   &
     &             4.*c2*(xcm-x1)*(zcm-z2) - 2.*c2*(zcm-z1)*(xcm-x2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 1) then
          factor = (xcm-x1)*(zcm-z1)*(xcm-x2)*(zcm-z2) + c2*c2 +        &
     &             c2*(zcm-z1)*(zcm-z2) + c2*(xcm-x1)*(xcm-x2)
          end if
          if (imu .eq. 1 .and. jmu .eq. 2) then
          factor = (xcm-x1)*(zcm-z1)*(xcm-x2)**2 +                      &
     &             2.*c2*(zcm-z1)*(xcm-x2)
          end if

          if (imu .eq. 2 .and. jmu .eq. 2) then
          factor = (xcm-x1)**2*(xcm-x2)**2 + 2.*c4 -2.*c2*c2 +          &
     &             4.*c2*(xcm-x1)*(xcm-x2)
          end if
          if (imu .eq. 2 .and. jmu .eq. 1) then
          factor = (xcm-x2)*(zcm-z2)*(xcm-x1)**2 +                      &
     &             2.*c2*(zcm-z2)*(xcm-x1)
          end if
          if (imu .eq. 2 .and. jmu .eq. 0) then
          factor = (2.*(zcm-z2)**2 - (xcm-x2)**2)*(xcm-x1)**2 -         &
     &             4.*c2*(xcm-x2)*(xcm-x1)
          end if

        end if ! end if of DD

        rmatel=rmatel+factor*coef1*coef2*bigQ*gM**3
        end do
        end do

! Put in spherical harmonic type factors...
        rmatel=rmatel*gfactor(l1,imu)*gfactor(l2,jmu)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
