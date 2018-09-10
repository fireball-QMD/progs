!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Module to define the variables we use in the DOS (dos,   !!!
!!! hopping or atom) calculation (except hamk)               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module module_dos

     complex, dimension (:,:,:), allocatable :: green
     real, dimension(:,:,:,:), allocatable   :: hr_box 
! variable bellow are defined in MODULES/neighbor_map.f90 now
!    integer, dimension(:), allocatable      :: neighn_tot
!    integer, dimension(:,:), allocatable    :: neighj_tot
!    integer, dimension(:,:), allocatable    :: neighb_tot
!    integer                                 :: num_neig_maxtot

     integer natom_beg           ! number of beggining atom for the DOS
     integer natom_end           ! number of final atom for the DOS
     integer norb_act            ! total # of active orbitals in STM
     integer nener               ! # of enegies for the dos calculation
     integer iwrttip             ! write out the tip_e_str.inp for the STM
     real ener_beg               ! initial energy for the DOS
     real ener_step              ! energy step
     real ener_min               ! energy minimum for tip_e_str.inp
     real ener_max               ! energy maximum for tip_e_str.inp
     real eta                    ! imaginary part in the DOS calculation
     real lattice                ! lattice parameter in atomic units
!dosng (option for computing DOS and electronci states without computing the GF)
     real, dimension (: ,:, :), allocatable :: dngcof
     real, dimension (:, :), allocatable    :: E_KS
     real, dimension (:), allocatable       :: DOS_total
     real, dimension (:,:), allocatable     :: States_total
     integer                                :: dstep !to count time inside dosng
!End of dosng
   
    end module

