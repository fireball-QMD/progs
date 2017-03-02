        module tdse

         integer :: netime        ! number of electron time steps within one ionic step
         integer :: nelec         ! number of electrons in the system
         real    :: dte           ! electron time step

         integer :: nexcite
         integer, parameter :: nexcite_max = 10     ! max number of excited electrons
         integer, dimension (nexcite_max) :: idelec ! an excited electron
         integer, dimension (nexcite_max) :: eband  ! a band where electron was excited
         real, dimension (nexcite_max) :: hoccup    ! hole occupancy

         complex, dimension (:,:,:), allocatable :: psi     ! TD-wf in MO basis set
         complex, dimension (:,:,:), allocatable :: psiAO   ! TD-wf in AO basis set
         complex, dimension (:,:,:), allocatable :: HLow    ! H in MO basis set
         complex, dimension (:,:,:), allocatable :: sm12    ! S^(-1/2)
         real, dimension (:,:), allocatable :: Enev         ! Energy expectation value of electron
         complex, dimension (:,:,:), allocatable :: psi2es  ! projection psi to eigenfuction

! writeout wavefunction
         integer  :: np2es       ! number of ionic steps after psi is projected on eigenstate
         logical  :: isp2es      ! indicator of firts attempt two writeout data


       end module tdse
