        module dimensions

! =============================================================================
! These should not have to be changed, because they have no effect on
! memory use within fireball.

! Maximum number of species
         integer, parameter :: nspec_max = 20

! =============================================================================
! These control the integration quadrature.

! Maximum number of 2-center grid points for interpolation
         integer, parameter :: nfofx = 207        

! Maximum x,y grids for 3 center interactions. (nbcba,nnaba in create)
         integer, parameter :: numXmax = 31
         integer, parameter :: numYmax = 31

! =============================================================================
! These should not be changed.  They have to do with the inner workings of 
! the fireball program
! Maximum number of 2-center interactions
         integer interactions2c_max

! The number of integration theta
         integer, parameter :: ithetamax = 5
         integer, parameter :: ntheta = 5
        end module
