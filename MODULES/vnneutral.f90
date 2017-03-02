 module vnneutral
! This module provide basic information about neutral atomic potential 
! of each specie 
! ===========================================================================


! max. available w.f. mesh points
   integer, parameter  :: max_vna_points = 5000

! real w.f. mesh size
   integer, dimension (:), allocatable :: mesh_na 
! r_max
   real, dimension (:), allocatable :: rmax_na 
! interval of w.f. mesh of each specie
   real, dimension (:), allocatable ::  drr_na
! mesh of distance 
   real, dimension (:,:), allocatable :: rr_na

! w.f. spline array
   real, dimension (:,:), allocatable :: vnna_spline 
! w.f. array
   real, dimension (:,:), allocatable :: vnna

 end module vnneutral
