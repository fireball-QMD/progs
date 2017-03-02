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

! readdata.f90
! Program Description
! ===========================================================================
!       This routine reads the Extended Hubbard data 
!
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readdata_eh ()

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
        write (*,*) '  '
        write (*,100)
        write (*,*) ' READING DATA FROM DATA FILES IN FDATA DIRECTORY: '
        write (*,100)
! one-center
        write (*,*) '  '
        write (*,*) ' Calling read_1c for 1-center exchange-correlation '
        call read_1c (nspecies, itheory, itheory_xc, ispin, ioff2c(7), nsup, &
     &                nsu, fdataLocation)
! two-center
        interaction_start = 1
        do interaction = interaction_start, 5
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx, fdataLocation)
        end do
! horsfield interaction 
        if (itheory_xc .ne. 1) then
         !  xc_ontop
         interaction = 6
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx, fdataLocation)
         ! xc_atom-atom
         interaction = 7
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx, fdataLocation)
         ! xc_correction
         interaction = 8
         write (*,*) '  '
         write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
         write (*,100)
         call read_2c (interaction, nspecies, itheory,                  &
     &                 ioff2c(interaction), nzx, fdataLocation)

        end if
 
!       interaction = 10, 11 are x and y dipole
        interaction = 12
        write (*,*) '  '
        write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        write (*,100)
        call read_2c (interaction, nspecies, itheory,            &
     &                ioff2c(interaction), nzx, fdataLocation)
        interaction = 13
        write (*,*) '  '
        write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        write (*,100)
        call read_2c (interaction, nspecies, itheory,            &
     &                ioff2c(interaction), nzx, fdataLocation)
 
! Extended-hubbard addition
        
        interaction = 14
        write (*,*) '  '
        write (*,*) ' Calling read_2c for 2-Center Interaction # ', interaction
        write (*,100)
        call read_2c (interaction, nspecies, itheory,           &
     &                 ioff2c(interaction), nzx, fdataLocation)

! Spherical OLSXC exchange-correlation
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2) then 
         do interaction = 15, 23
          write (*,*) '  '
          write (*,*) ' Calling read2c for 2-Center Interaction # ', interaction
          write (*,100)
          call read_2c (interaction, nspecies, itheory, ioff2c(interaction), &
     &                  nzx, fdataLocation)
         end do
        end if
! Do not need xintegral_2c if doing superspline (use splineint_2c)
        if (superspline) deallocate (xintegral_2c)



        if (igauss .eq. 0) then
         interaction = 1   ! bcna
         write (*,*) '  '
         write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
         write (*,100)
         call read_3c (interaction, nspecies, ioff3c(interaction), itheory, nzx, &
     &                  fdataLocation)
         if (itheory_xc .eq. 0) then 
          interaction = 2  ! xc3c (Horsfield) 
          write (*,*) '  '
          write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
          write (*,100)
          call read_3c (interaction, nspecies, ioff3c(interaction), itheory, &
     &                  nzx, fdataLocation)
         end if

         if ((itheory_xc .eq. 1) .or. (itheory_xc .eq. 2)) then 
          interaction = 3   ! den3 (3c - OLSXC) - average density
          write (*,*) '  '
          write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
          write (*,100)
          call read_3c (interaction, nspecies, ioff3c(interaction), itheory, &
     &                  nzx, fdataLocation)
          interaction = 4   ! den3 (3c - OLSXC) - spherical average density
          write (*,*) '  '
          write (*,*) ' Calling read3c for 3-Center Interaction # ', interaction
          write (*,100)
          call read_3c (interaction, nspecies, ioff3c(interaction), itheory, &
     &                  nzx, fdataLocation)
         end if
         write (*,100)

! Set up some tables for the 2d interpolator
         call setterp_2d (nspecies, itheory_xc, itheory)
        end if

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

         call readgaussG (nspecies, nzx, fdataLocation)
        end if
 


! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine readdata_eh
 
