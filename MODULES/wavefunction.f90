 module wavefunction
! This module provide basic information about wavefunctions of each specie 
! ===========================================================================


! max. available w.f. mesh points
   integer, parameter  :: wfmax_points = 5000

! real w.f. mesh size
   integer, dimension (:,:), allocatable :: mesh_wf 
! interval of w.f. mesh of each specie
   real, dimension (:,:), allocatable ::  drr_wf
! mesh of distance 
   real, dimension (:,:,:), allocatable :: rr_wf
! rmax
   real, dimension (:,:), allocatable :: rmax_wf

! w.f. spline array
   real, dimension (:,:,:), allocatable :: wf_spline 
! w.f. array
   real, dimension (:,:,:), allocatable :: wf

 end module wavefunction
