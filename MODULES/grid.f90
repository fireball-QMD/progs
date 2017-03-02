 module grid
! This module provide information about numerical grid used to solve Poisson
! equation and Vxc stuff
! ===========================================================================

! wave function
!----------------------------------
! name of files
!   character (100), dimension(:,:), allocatable :: name_wf
!   character (100), dimension (:), allocatable :: name_na
! mesh size
   integer, parameter :: max_mesh = 5000
! energy cutoff
   real :: Ecut

! Rc_max
   real :: Rc_max

! the regular mesh (spread out over the unit cell)
   integer :: rm1
   integer :: rm2
   integer :: rm3
   integer :: nrm

! local fine grid (spread out over elementary mesh cell, used for the wf. interpolation)
!   integer, parameter :: mloc = 1
!   real, parameter :: dloc = 1.0d0
!   real, dimension (3,9) :: xmloc

! the extended mesh (defined by the overlap of the atoms in neighbors periodic cells)
   integer :: em1
   integer :: em2
   integer :: em3
   integer :: nem

! the atomic mesh
   integer :: am1
   integer :: am2
   integer :: am3
   integer :: nam

! offset of extended mesh:  emi = cmi+ 2*mexi
   integer :: emx1
   integer :: emx2
   integer :: emx3

! Map the extended mesh to the regular mesh
   integer, dimension (:), allocatable :: e2r

   integer :: noff

! number of atomic mesh point within Rc_max
   integer, dimension (:), allocatable :: am2rc
   real, dimension (:,:), allocatable :: ram2rc

   real, dimension (:,:), allocatable :: ratom2g

! mesh for finite difference method
   integer :: mfd1
   integer :: mfd2
   integer :: mfd3
   integer :: nmfd

! maping array between extended and normal mesh
   integer, dimension (:), allocatable :: e2n
   integer, dimension (:), allocatable :: n2e

! modified by honza
! auxiliary mesh of local neighbors map
! (1,:),(2,:),(3,:) are x,y,z; (4,:) is auxiliary (e.g. for the atom itself)
! (:,1) is forward, (:,2) is backward  
   integer :: nneighij
   integer, dimension (4,2) :: neighij
! coefficients for the laplace operator (2nd derivative)
   real, dimension (4,2) :: d2f
! coefficients for 1st derivative
   real, dimension (4,2) :: d1f
! aux elementary vector
   real, dimension (3,3) :: ervec

! end - modified by honza   
   
! initial point of the grid
   real, dimension (3) :: g0

! flag to fix the initial point of the grid
   integer :: ifixg0

! unit (elementary) lattice vector of the submesh
   real, dimension (3,3) :: elvec

! unit (elementary) lattice vector of the submesh in the reciprocal space
!   real, dimension (3,3) :: relvec

! the reciprocal lattice vector
   real, dimension (3,3) :: rlvec

! grid spacing
   real :: dr1
   real :: dr2
   real :: dr3
   real :: dr12
   real :: dr22
   real :: dr32



! elementary volume
   real  :: dvol
! neutral atomic potential
   real, target, dimension (:), allocatable :: vnaG
! Hartree potential (related to density0)
   real, target, dimension (:), allocatable :: vcaG
! xc potential
   real, target, dimension (:), allocatable :: vxcG
! variation of density
   real, target, dimension (:), allocatable :: drhoG
! atomic density
   real, target, dimension (:), allocatable :: rhoG0

! mapped wavefunction onto the mesh
!   real, target, dimension (:,:,:), allocatable :: psi2m
!   real, target, dimension (:,:), allocatable :: ipsi2m
!   real, target, dimension (:,:,:,:,:), allocatable :: psi22m
!   real, target, dimension (:,:,:), allocatable ::ipsi22m
!   real, target, dimension (:,:), allocatable :: npsi22m

! XC energy
!   real :: Etot_XC
! Hartree energy
!   real :: Etot_H


! ++++++++++ Projection of the bands ++++++++++++
   integer npbands
   integer, parameter :: npbands_max = 200
   integer, dimension(npbands_max) :: pbands
   integer iewform
   real  ewfewin_max
   real  ewfewin_min


 end module grid
