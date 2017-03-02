module dynamo

! displacement 
  real        u0disp
! dimension
  integer     ndx
  integer, dimension (3) :: ndvec

! Dynamical matrix
  real, dimension (:,:), allocatable :: phidm

! Forces vector from first (+) displacement
  real, dimension (:,:), allocatable :: ftot1

! Displacement vector 
  real, dimension (:), allocatable :: u0vec

! Displacement vector 
  real, dimension (:,:), allocatable :: ratom0

! number of atoms to move
  integer     natoms_dm

! list of moved atoms 
  integer, dimension (:), allocatable :: jatoms_dm

! name of a file where dynamical matrix will be stored
  character(len=30)  :: filephi

! flag determines which dispalcement is doing +(-)u0disp
  logical   ldynamo

! local time clock
  integer   ltime

! *** e-ph couplings
! number of e-ph coupling to be evaluated
  integer  nephc

! name of a file where e-ph couplings will be stored
  character(len=30)  :: fileephc

! saving eigenavalues
  real, dimension (:,:), allocatable :: eigsave
 
! Derivatives of eigenvalues
  real, dimension (:,:), allocatable :: deigen

! Reference eigenvalues (zero displacement)
  real, dimension (:), allocatable :: eigref

! List of eigenstate from which e-ph coupling will be evaluated
  integer, dimension (:), allocatable :: eiglist

! number of vibrational modes be considered for e-ph coupling
  integer  nvmodes

! list of normal modes 
  integer, dimension (:), allocatable :: idvmode

! temperature for which e-ph coupling will be calculated
  real temp_ephc

! effective mass for e-ph coupling
  real mass_ephc


end module dynamo
