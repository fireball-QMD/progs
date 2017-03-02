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

! DgelementsGS_overlap.f90
! Program Description
! ===========================================================================
! This routine gets drivatives of the overlap matrix element - using gaussian wavefactions
!
! ===========================================================================
! Code rewritten by:
! Hao Wang
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-9270
! ===========================================================================
!
! Program Declaration
! =====================================================================
       subroutine DgelementsGS_overlap(l1,l2,in1,in2,ish,jsh,y,rmatel)
! =====================================================================
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
        integer, intent (in) :: l1
        integer, intent (in) :: l2
        integer, intent (in) :: ish
        integer, intent (in) :: jsh

        real, intent (in) :: y

! ===========================================================================
! output
        real rmatel

! Local Parameters and Data Declaration
! ===========================================================================
! MHL
!        integer malpha

        integer ial1
        integer ial2

        real rmu
        real lilm
        real d12

        real bigQ
        real gM

        real al1
        real al2

        real coef1
        real coef2

        real dQ
        real factor
        real dfactor


        real x1
        real y1
        real z1

        real x2
        real y2
        real z2

! Proceedure
! =====================================================================

        d12 = y

        x1 = 0.
        y1 = 0.
        z1 = -y/2.
        x2 = 0.
        y2 = 0.
        z2 = y/2.

        rmatel = 0.d0
! ---------------------------------------------------------------------
! MHL (Sep. 29. 2004)

! Wavefunction for shperical approximation (begins from 1)
!        real gcoefficientsPSIS(max_alphas,1:nsh_max,nspec_max)
!        real alphaPSIS(max_alphas,1:nsh_max,nspec_max)
!        integer nalphaPSIS(1:nsh_max,nspec_max)

! ----------------------------------------------------------------------
! MHL (nalpha --> nalphaPSI, etc.)
! alpha1 (spherical approximation, nalphaPSI -> nalphaPSIS)
        do ial1=1,nalphaPSIS(ish,in1)
        al1=alphaPSIS(ial1,ish,in1)
        coef1=gcoefficientsPSIS(ial1,ish,in1)

! alpha2 (spherical approximation, nalphaPSI -> nalphaPSIS)
        do ial2=1,nalphaPSIS(jsh,in2)
        al2=alphaPSIS(ial2,jsh,in2)
        coef2=gcoefficientsPSIS(ial2,jsh,in2)

        rmu=al1*al2/(al1+al2)
        lilm=al1+al2

! Big Q
        bigQ=exp(-rmu*d12**2)
        gM=sqrt(pi/lilm)

! Derivative of bigQ...
            dQ=-2.*rmu*y
            dQ=dQ*bigQ

! For the spherical approxiamtion, only l1 = 0 and l2 = 0 considered,
! so m1 = 0 and m2 =0.
! IT IS ONLY ONE TERM NEED TO BE EVALUATED
         factor = 1.d0
        dfactor = 0.d0

        rmatel = rmatel + coef1*coef2*gM**3*(bigQ*dfactor + factor*dQ)
        end do
        end do

        rmatel = rmatel*gfactor(0,0)*gfactor(0,0)

        return
        end
! =======================================================================
