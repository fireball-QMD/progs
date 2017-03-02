 module optimization 
!=======================================================================
! Conjugate gradient optimization
!=======================================================================
   real                 :: cg_drmax
   real                 :: cg_dummy
   real                 :: energy_tol
   real                 :: force_tol
   integer              :: cg_iter
   integer              :: cg_maxstep
   integer              :: cg_minint
   integer              :: initer
   integer              :: cgiter
   integer              :: icg2md
   integer              :: istatus 
   integer		:: freeParamsCount
  
   real,dimension(:,:),allocatable    :: x0
   real,dimension(:,:),allocatable    :: f0
   real,dimension(:,:),allocatable    :: Q0
   real,dimension(:,:),allocatable    :: g
   real,dimension(:,:),allocatable    :: h
   real,dimension(:,:),allocatable      :: mask

! VARIABLES:
!   drmax      .. maximal permited displacement per atom of a trial step
!   dummmy     .. read the scale to reduce the search step if e1<e2 
!   energy_tol .. convergence criteria of total energy
!   force_tol  .. convergence criteria of forces
!   cg_maxstep ..
!   cg_minint  ..
 
! SAVE variables cgo
   integer      :: minflag
   integer      :: iflag2
   real         :: etot0
   real         :: etot1
   real         :: etot2
   real         :: alpha_cg
   real         :: gamma
! end JEL-CG ----------------------------------------------------


! FIRE variables (  Prokop Hapala for main_loop_FIRE.f90 )
   real         :: FIRE_finc
   real         :: FIRE_fdec
   real         :: FIRE_falpha
   integer      :: FIRE_Nmin
   real         :: FIRE_dtmax
   real         :: FIRE_dt
   real         :: FIRE_acoef0
   real         :: FIRE_acoef
   real         :: FIRE_mass
   real         :: FIRE_Ftot


        end module
