         module gaussG 
         use dimensions

         integer, parameter :: max_alphas = 10

! Note: in the 0:nsh_max slot  1 ---> wavefunction of 1st shell
!                              2 ---> wavefunction of 2nd shell
!                              etc... BUT:
!                              0 ---> neutral atom potential
        real, dimension (:,:,:), allocatable :: gcoefficients
        real, dimension (:,:,:), allocatable :: alpha
        integer, dimension (:,:), allocatable :: nalpha
! MHL
! For neutral atomic potential
         real, dimension (:,:), allocatable :: gcoefficientsVNA
         real, dimension (:,:), allocatable :: alphaVNA
         integer, dimension (:), allocatable :: nalphaVNA
! For electron density
         real, dimension (:,:,:), allocatable :: gcoefficientsN
         real, dimension (:,:,:), allocatable :: alphaN
         integer, dimension (:,:), allocatable :: nalphaN
! For wavefunction
         real, dimension (:,:,:), allocatable :: gcoefficientsPSI
         real, dimension (:,:,:), allocatable :: alphaPSI
         integer, dimension (:,:), allocatable :: nalphaPSI
! For wavefunction using spherical approximation
         real, dimension (:,:,:), allocatable :: gcoefficientsPSIS
         real, dimension (:,:,:), allocatable :: alphaPSIS
         integer, dimension (:,:), allocatable :: nalphaPSIS

! For electron density
         real, dimension (:,:), allocatable :: R_na
         real, dimension (:,:,:), allocatable :: gcoefficientsVNA_SH
         real, dimension (:,:,:), allocatable :: alphaVNA_SH
         integer, dimension (:,:), allocatable :: nalphaVNA_SH

         real, dimension (0:2,-2:2) :: gfactor

! The variable ptemp is "pseudo-temperature". It must NOT be zero!
! This is for making the bcna 3-centers go continuously to zero at RC.
! ptemp = 0.01 should be about right... Maybe smaller...
! Basically, the smaller ptemp is, the sharper the drop to zero at RC.
! But don't make ptemp too big because the bcna integrals won't be accurate.
! See tresgaussians and Dtresgaussians for details.
         real, parameter :: ptemp = 0.01d0

         end module
