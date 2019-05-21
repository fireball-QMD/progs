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


! readdata_mcweda.f90
! Program Description
! ===========================================================================
!       This routine reads the data of McWeda approximation
!
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readdata_mcweda ()

        use configuration
        use constants_fireball
        use interactions 
        use integrals
        use options
        use gaussG
        use charges
        
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
 
! Local Variable Declaration and Description
! ===========================================================================
        integer interaction
        integer interaction_start
 
! Procedure
! ===========================================================================
! ============================================================================
!                              read data files
! ============================================================================

! ****************************************************************************
! About ioff2c and ioff3c - the interactions are turned on or off according to 
! the ioff2c and ioff3c arrays which are read in from file diagnostics.input. 

! Two Center Interactions:
! Note the last three short range Ewald, long range Ewald, and coulomb
! (extended Hubbard) are not interactions from datafiles like the rest.
! However, we give the option here to turn off these interactions anyways.
!                ioff - 2c overlap
!                ioff - 2c vna_ontopl
!                ioff - 2c vna_ontopr
!                ioff - 2c vna_atom-atom
!                ioff - 2c non-local
!                ioff - 2c xc_ontop
!                ioff - 2c xc_atom-atom
!                ioff - 2c xc_correction
!                ioff - 2c z-dipole
!                ioff - 2c y-dipole
!                ioff - 2c x-dipole
!                ioff - 2c coulomb
!                ioff - 2c kinetic
!                ioff - 2c extended Hubbard
!                ioff - 2c den_ontopl
!                ioff - 2c den_ontopr
!                ioff - 2c den_atom
!                ioff - 2c denS_ontopl
!                ioff - 2c denS_ontopr
!                ioff - 2c denS_atom
!                ioff - 2c_overlapS
!                ioff - 2c coulomb (extended Hubbard)
!                ioff - 2c short range Ewald
!                ioff - 2c long range Ewald
!
! Three Center Interactions:
!                ioff - 3c neutral-atom
!                ioff - 3c exchange-correlation
!                ioff - 3c average density OLSXC
!                ioff - 3c average density OLSXC (spheric)
! ****************************************************************************


        write (*,*) '  '
        write (*,100)
        write (*,*) ' REEEEEADING DATA FROM DATA FILES IN FDATA DIRECTORY: '
        write (*,100)
! one-center
        write (*,*) '  '
        write (*,*) ' Calling   for 1-center exchange-correlation '
        call read_1c (nspecies, itheory, itheory_xc, ispin, ioff2c(7))
! two-center
        interaction_start = 1
        do interaction = interaction_start, 5
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx)
        end do
! horsfield-like interaction (not needed for SNXC) 
        if (itheory_xc .ne. 1) then
         !  xc_ontop
         interaction = 6
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx)
         ! xc_atom-atom
         interaction = 7
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx)
         ! xc_correction
         interaction = 8
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx)
        end if
        
        if (itheory .eq. 1) then
         interaction = 9
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx)
        end if

! JIMM
!       interaction = 10, 11 are y and x dipole
        if (idipole .eq. 1) then

         if (itheory .eq. 1) then
          interaction = 10
          write (*,*) '  '
          write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
          write (*,100)
          call read_2c (interaction, nspecies, itheory,                  &
     &                  ioff2c(interaction), nzx)
         end if

         if (itheory .eq. 1) then
          interaction = 11
          write (*,*) '  '
          write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
          write (*,100)
          call read_2c (interaction, nspecies, itheory,                  &
     &                  ioff2c(interaction), nzx)
         end if

        end if ! idipole .eq. 1
 
        interaction = 12
        write (*,*) '  '
        write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        write (*,100)
        call read_2c (interaction, nspecies, itheory,            &
     &                ioff2c(interaction), nzx)
        interaction = 13
        write (*,*) '  '
        write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        write (*,100)
        call read_2c (interaction, nspecies, itheory,            &
     &                ioff2c(interaction), nzx)
 
! Spherical OLSXC exchange-correlation
        do interaction = 15, 23
         write (*,*) '  '
         write (*,*) ' Calling read2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory, ioff2c(interaction), &
     &                  nzx)
        end do

! Do not need xintegral_2c if doing superspline (use splineint_2c)
        if (superspline) deallocate (xintegral_2c)


        if (igauss .eq. 0) then
         interaction = 1   ! bcna
         write (*,*) '  '
         write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx)
                   
         interaction = 3   ! den3 (3c - OLSXC) - average density
         write (*,*) '  '
         write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, &
     &                  nzx)
         interaction = 4   ! den3 (3c - OLSXC) - spherical average density
         write (*,*) '  '
         write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, &
     &                  nzx)
 
         write (*,100)

! Set up some tables for the 2d interpolator
         call setterp_2d (nspecies, itheory_xc, itheory)
        end if ! end if igauss

        if (igauss .eq. 1) then
         gfactor(0,0)  = sqrt(1.0d0/4.0d0/pi)
         gfactor(1,-1) = sqrt(3.0d0/4.0d0/pi)
         gfactor(1,0)  = sqrt(3.0d0/4.0d0/pi)
         gfactor(1,1)  = sqrt(3.0d0/4.0d0/pi)
         gfactor(2,-2) = sqrt(15.0d0/4.0d0/pi)
         gfactor(2,-1) = sqrt(15.0d0/4.0d0/pi)
         gfactor(2,0)  = sqrt(15.0d0/12.0d0/4./pi)
         gfactor(2,1)  = sqrt(15.0d0/4.0d0/pi)
         gfactor(2,2)  = sqrt(15.0d0/4.0d0/pi)/2.0d0
! MHL
! NA potential
         allocate (gcoefficientsVNA (max_alphas, nspec_max))
         allocate (alphaVNA (max_alphas, nspec_max))
         allocate (nalphaVNA (nspec_max))
! Electron density
         allocate (gcoefficientsN (max_alphas, 0:nsh_max, nspec_max))
         allocate (alphaN (max_alphas, 0:nsh_max, nspec_max))
         allocate (nalphaN (0:nsh_max, nspec_max))
! Wavefunction
         allocate (gcoefficientsPSI (max_alphas, 1:nsh_max, nspec_max))
         allocate (alphaPSI (max_alphas, 1:nsh_max, nspec_max))
         allocate (nalphaPSI (1:nsh_max, nspec_max))
! Wavefunction for spherical approximation
         allocate (gcoefficientsPSIS (max_alphas, 1:nsh_max, nspec_max))
         allocate (alphaPSIS (max_alphas, 1:nsh_max, nspec_max))
         allocate (nalphaPSIS (1:nsh_max, nspec_max))
! NA potential by shell
         allocate (gcoefficientsVNA_SH (max_alphas, 1:nsh_max, nspec_max))
         allocate (alphaVNA_SH (max_alphas, 1:nsh_max, nspec_max))
         allocate (nalphaVNA_SH (1:nsh_max, nspec_max))
         allocate (R_na (1:nsh_max, nspec_max))

         call readgaussG (nspecies, nzx )
        end if
 
! Broadcast some variables to other processors
        if (iordern .eq. 1) then
         call ME_max_bcast
         call readdata_ordern_init (nspecies, ioff2c)
        end if



! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine readdata_mcweda
 
