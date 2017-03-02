 module scf
! This module provide information about SCF cycle
! ===========================================================================

! max. number of SCF iterations
 integer  max_scf_iterations
! SCF step
 integer  Kscf

! cDFT vlada
 integer id_hole
 integer id_elec
 real occup_hole
 real occup_elec
 logical  cDFT_active
 real, dimension (:,:), allocatable :: wf_elec
 real, dimension (:,:), allocatable :: wf_hole
! cDFT vlada

! idmix: Defines the mixing procedure: idmix=1 means simple mixing
! For larger idmix values, the choice of bmix becomes less important
 integer, parameter ::  idmix = 6

 real tempfe                 ! Fermi temperature (Smearing)

 real   bmix
 real   sigmatol
 real   sigmaold
 real   sigma

! mixing algorithm
 integer ialgmix

 logical  scf_achieved

 end module scf
