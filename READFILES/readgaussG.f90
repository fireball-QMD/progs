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
 
! readgaussG.f90 
! Program Description
! ===========================================================================
!       This routine reads in the alpha's and coef's from gauss.dat file.
! It also makes nalpha(ish,ispec)= number of alphas (i.e., gaussians) in
! the expansion of the ish'th shell of the ispec'th species. If ish = 0
! then the expansion is of the neutral atom potential...
!
! WANG begin, gauss, Oct. 1, 2001
! If ish = -1, then the expansion is of the electron density
! WANG end
!
!         e.g.:  ish = 1   --->  expansion of wavefunction of 1st shell
!              ish = 0   --->  expansion of neutral atom potential...
! WANG gauss,  ish = -1, --->  expansion of the electron density
! ===========================================================================
! Code written by John Tomfohr
! Department of Physics and Astronomy
! Arizona State University 
! Tempe, AZ 85287-1504 
! Office Telephone (480) 965-0667
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine  readgaussG (nspecies, nzx, fdataLocation)

        use dimensions
        use interactions
        use gaussG

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer nspecies
        integer nzx(nspec_max)
        character (len=200), intent(in) ::fdataLocation
 
!! Output
! Beginning of old statement
! WANG gauss Oct. 1, 2001
! Change the array index 0:nsh_max to -1:nsh_max  
! -1 corresponding to gauss fit of electron density
! New definations of old statement
!        real gcoefficients(max_alphas,-1:nsh_max,nspec_max)
!        real alpha(max_alphas,-1:nsh_max,nspec_max)
!        integer nalpha(-1:nsh_max,nspec_max)
!
! Old definations 
!        real gcoefficients(max_alphas,0:nsh_max,nspec_max)
!        real alpha(max_alphas,0:nsh_max,nspec_max)
!        integer nalpha(0:nsh_max,nspec_max)
! End of old statement 
! ==========================================================================
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

! Neutral atom potential by shell
!        real gcoefficientsVNA_SH(max_alphas,0:nsh_max,nspec_max)
!        real alphaVNA_SH(max_alphas,0:nsh_max,nspec_max)
!        integer nalphaVNA_SH(0:nsh_max,nspec_max)

! Local Variable Declaration and Description
! ===========================================================================
        integer i,ispec,ish,icounter
        integer nspec,ial,natomic,numshells

! WANG gauss, Oct. 1, 2001, introduce ish00 = 0 OR -1
! ish00 = 0,  neutral atom potential       
! ish00 = -1, electron density        
        integer ish00

        real alpha_junkS(max_alphas,1:nsh_max,nspec_max)
        real coef_junkS(max_alphas,1:nsh_max,nspec_max)
        real alpha_junk(max_alphas,nspec_max)
        real coef_junk(max_alphas,nspec_max)
        real R_na_junk(1:nsh_max,nspec_max)
        integer num_als_junk        
 
! Procedure
! ===========================================================================

        write(*,*)
        write(*,*)' Welome to gauss_read'
        write(*,*)

        open(unit=44, file=trim(fdataLocation)//'/gauss.dat', status='unknown')

        do i=1,4
        read(44,*)
        end do
        read(44,*)nspec
        write(*,*)' Number of species in gauss.dat  = ',nspec
        write(*,*)' Number of species to be read in = ',nspecies

! WANG gauss, Oct. 1, 2001.
! Read in ish00
        read(44,*)
        read(44,*)ish00

        icounter=1

        do ispec=1,nspec
        read(44,*)
        read(44,*)
        read(44,*)natomic
        read(44,*)
        read(44,*)numshells

        if (natomic.eq.nzx(icounter)) then
           if (numshells.ne.nssh(icounter)) then
              write(*,*)' numshells in gauss.dat and info.dat differ.'
              write(*,*)numshells,nssh(icounter)
              write(*,*)' for species with atomic number = ',natomic
              stop
           end if 
! Read gaussians for total electron density. 
! ish=0 means total electron density.
        ish=0
        read(44,*)
        read(44,*)nalphaN(ish,icounter)
        read(44,*)
        read(44,*)(alphaN(ial,ish,icounter),ial=1,nalphaN(ish,icounter))
        read(44,*)
        read(44,*)(gcoefficientsN(ial,ish,icounter),ial=1,nalphaN(ish,icounter))

! Read gaussians for neutral atomic potential. 
        read(44,*)
        read(44,*)nalphaVNA(icounter)
        read(44,*)
        read(44,*)(alphaVNA(ial,icounter),ial=1,nalphaVNA(icounter))
        read(44,*)
        read(44,*)(gcoefficientsVNA(ial,icounter),ial=1,nalphaVNA(icounter))

! Read gaussians for wavefunction
        do ish=1,numshells
        read(44,*)
        read(44,*)nalphaPSI(ish,icounter)
        read(44,*)
        read(44,*)(alphaPSI(ial,ish,icounter),ial=1,nalphaPSI(ish,icounter))
        read(44,*)
        read(44,*)(gcoefficientsPSI(ial,ish,icounter),ial=1,nalphaPSI(ish,icounter))
        end do
 
! Read gaussians for wavefunction (sqrt(psi**2))
        do ish=1,numshells
        read(44,*)
        read(44,*)nalphaPSIS(ish,icounter)
        read(44,*)
        read(44,*)(alphaPSIS(ial,ish,icounter),ial=1,nalphaPSIS(ish,icounter))
        read(44,*)
        read(44,*)(gcoefficientsPSIS(ial,ish,icounter),ial=1,nalphaPSIS(ish,icounter))
        end do

! Read gaussians for electron density of each shell
        do ish=1,numshells
        read(44,*)
        read(44,*)nalphaN(ish,icounter)
        read(44,*)
        read(44,*)(alphaN(ial,ish,icounter),ial=1,nalphaN(ish,icounter))
        read(44,*)
        read(44,*)(gcoefficientsN(ial,ish,icounter),ial=1,nalphaN(ish,icounter))
        end do

! Read gaussians for electron density of each shell
        do ish=1,numshells
        read(44,*)
        read(44,*)nalphaVNA_SH(ish,icounter)
        read(44,*)
        read(44,*)(alphaVNA_SH(ial,ish,icounter),ial=1,nalphaVNA_SH(ish,icounter))
        read(44,*)
        read(44,*)(gcoefficientsVNA_SH(ial,ish,icounter),ial=1,nalphaVNA_SH(ish,icounter))
        read(44,*)
        read(44,*)R_na(ish,icounter)
        end do

        icounter=icounter+1

        else
        do i  = 1, 2
        read(44,*)
        read(44,*)num_als_junk
        read(44,*)
        read(44,*)(alpha_junk(ial,icounter),ial=1,num_als_junk)
        read(44,*)
        read(44,*)(coef_junk(ial,icounter),ial=1,num_als_junk)
        end do

        do  i = 1, 3
         do ish=1,numshells
         read(44,*)
         read(44,*)num_als_junk
         read(44,*)
         read(44,*)(alpha_junkS(ial,ish,icounter),ial=1,num_als_junk)
         read(44,*)
         read(44,*)(coef_junkS(ial,ish,icounter),ial=1,num_als_junk)
         end do
        end do

         do ish=1,numshells
         read(44,*)
         read(44,*)num_als_junk
         read(44,*)
         read(44,*)(alpha_junkS(ial,ish,icounter),ial=1,num_als_junk)
         read(44,*)
         read(44,*)(coef_junkS(ial,ish,icounter),ial=1,num_als_junk)
         read(44,*)
         read(44,*)R_na_junk(ish,icounter)
         end do

        end if
! ispec
        end do

        if (icounter-1.ne.nspecies) then
        write(*,*)' You are missing one or more species in your gauss.dat file...'
        stop
        end if

        close (unit=44)

        write(*,*)' Done reading alphas and coefs for gaussian expansions.'
        write(*,*)
        write(*,*)' Bye from gauss_read!!!'
        write(*,*)

        return
        end
