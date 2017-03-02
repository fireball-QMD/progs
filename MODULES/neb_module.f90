 module neb
! This module contains variables definition for Nugged Elastic Band method 
! ===========================================================================

!----------------------------------
! name of files
!   character (100), dimension(:,:), allocatable :: name_wf
!   character (100), dimension (:), allocatable :: name_na


! number of neb images
   integer :: nimg_neb

! the spring constant
   real :: k_neb

! displacement tolerance controling the NEB convergence
   real :: tol_displ_neb

! displacement tolerance controling the NEB convergence
   real :: tol_ftot_neb

! displacement tolerance controling the NEB convergence
   real :: tol_etot_neb

! maximal of NEB iterations
   integer :: niter_neb_max

! NEB iteration  
   integer :: iter_neb

!
! total energy NEB of the images
   real, dimension (:), allocatable :: etot_neb

! difference total energy 
   real, dimension (:), allocatable :: detot_neb

! define time step used for the minimisation verlet algorithm
   real  :: dt_neb

! comment; we should think also to store charges for each 
! image to save time for next iterations 
 

! neutral atomic potential
   real, dimension (:,:,:), allocatable :: ratom_neb
   real, dimension (:,:,:), allocatable :: ftot_neb
   real, dimension (:,:,:), allocatable :: vatom_neb

! the tangent vector
   real, dimension (:,:), allocatable :: tang
! the spring force
   real, dimension (:,:), allocatable :: Fs
! the projected (rectangular) total force 
   real, dimension (:,:), allocatable :: Frec
! the true NEB force acting on an imgae 
   real, dimension (:,:,:), allocatable :: Fneb


 end module neb
