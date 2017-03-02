module transport
! Description 
! module contains variables related to transport calculation 
!

! basic variables

  logical ichannel
  integer ifithop
  logical ieta
  logical iwrt_trans

! lower bound energy
  real :: Elow
! upper bound energy
  real :: Eup
! energy step
  complex :: dE
! number of energy steps
  integer :: nE
! imag part of energy
  real :: eta
! optional imag part of energy
  real :: eta0
! list of atoms on which the optional eta0 will be applied
  integer, dimension (:), allocatable   :: ideta

  type system    
! sample
     integer                              :: natom
     integer                              :: norb
     integer                              :: ncell
     integer, dimension(:), pointer       :: atom
! tip
     integer                              :: natom_tip
     integer                              :: norb_tip
     integer, dimension(:), pointer       :: atom_tip
     real, dimension(:,:), pointer        :: ratom_tip
     integer                              :: nspec
     integer, dimension(:), pointer       :: spec
     integer, dimension(:), pointer       :: ispec
     integer, dimension(:), pointer       :: t2s
  end type system

! poniters mapping atoms into tip's arrays
  integer, dimension (:), allocatable  :: pointer1
  integer, dimension (:), allocatable  :: pointer2

! Hamiltonian in k-space of sample_x
  complex,    dimension (:,:), allocatable    :: Hsam1_k
  complex,    dimension (:,:), allocatable    :: Hsam2_k
! auxiliar Green's function in k-space of sample_1 
!  complex,    dimension (:,:), allocatable    :: Gsam1_k
!  complex,    dimension (:,:), allocatable    :: Idn1
! auxiliar Green's function in k-space of sample_2 
!  complex,    dimension (:,:), allocatable    :: Gsam2_k
!  complex,    dimension (:,:), allocatable    :: Idn2
! Orthogonal unperturbed Hamiltonian in kspace
  complex,    dimension (:,:,:), allocatable    :: H_k

! Dos of sample_x
!  real,  dimension (:,:), allocatable  :: dos1
!  real,  dimension (:,:), allocatable  :: dos2

! Green's function of tip
! retarded G.f.
  complex,  dimension (:,:,:), allocatable  :: Gr_tip1
  complex,  dimension (:,:,:), allocatable  :: Gr_tip2
! advanced G.f.
  complex,  dimension (:,:,:), allocatable  :: Ga_tip1
  complex,  dimension (:,:,:), allocatable  :: Ga_tip2

! complex current matrix
  complex,  dimension (:,:,:), allocatable  :: Jc

! Denominators D_xx
! retarded 
!  complex,  dimension (:,:,:), allocatable  :: Dr_11
!  complex,  dimension (:,:,:), allocatable  :: Dr_22
! advanced 
!  complex,  dimension (:,:,:), allocatable  :: Da_11
!  complex,  dimension (:,:,:), allocatable  :: Da_22

! t_12
! retarded 
  complex,  dimension (:,:), allocatable  :: t_12
  complex,  dimension (:,:), allocatable  :: t_21

! mapping matrices
  integer,  dimension (:), allocatable  :: degelec1
  integer,  dimension (:), allocatable  :: degelec2

  type (system)   :: sample1 
  type (system)   :: sample2 

! hoppings
  real,  dimension (:,:,:,:), allocatable   :: hops
  real, dimension (:,:), allocatable        :: zh_min
  real, dimension (:,:), allocatable        :: zh_max
  integer, dimension (:,:), allocatable     :: nzh
  real, dimension (:,:,:), allocatable      :: rc_hfit
  integer, dimension (:,:,:,:), allocatable :: lsh_hop
  real,  dimension (:,:,:,:), allocatable   :: hop_spline

  integer :: max_int

end module transport
