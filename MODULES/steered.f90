        module steered
         integer umb_pair                    ! number of atoms pairs of concern

! First geometry atom number pair x
         integer, dimension (:), allocatable :: umb_p1

! Second geometry atom number pair x
         integer, dimension (:), allocatable :: umb_p2

! Umbrella potential parameters 
         real, dimension (:), allocatable :: umb_CFd_ini     ! Initial force constant
         real, dimension (:), allocatable :: umb_CFd_end     ! Final force constant
         real, dimension (:), allocatable :: umb_d0_ini      ! Initial equilibrium distance
         real, dimension (:), allocatable :: umb_d0_end      ! Final equilibrium distance
         real, dimension (:), allocatable :: umb_d0_step     ! Step
         real, dimension (:), allocatable :: umb_CFd_step    ! Step 

         real, dimension (:, :), allocatable :: fumb     ! Force

! Initial reaction coordinate
! During the MD simulation, this value is always determined with the input 
! geometry during the first step.
         real umb_react_coord_init 

! This is the time that we have to wait for before making the umbrella analysis
! (i.e. time of equilibration required for the MD trajectory)
        ! real umb_time_start 
                ! Shouldn't this be an integer??
        integer umb_time_start
                ! We want it to be a time step, not a real time

! -------------------------------------------------------
!  Analysis of the data and Preparation of the output         
! -------------------------------------------------------

! Summary file of the umbrella sampling analysis over one MD simulation, 
! i.e. one window along the reaction coordinate. This file will be used after
! in the WHAM analysis    
         character (len=30)  umb_output_file
                                                           
! Value of the reaction coordinate which is stored and for which we count the 
! number of occurences along the MD simulation in one window.  This value is 
! average each time we find another one very close. (comparison with umb_epsi)
         real, dimension (:,:), allocatable :: val_coord

! Counter for the different tabulars, number of different val_coord really 
! stored after comparison using umb_epsi
         integer ival_coord
        end module
