        module interactions

!$ volatile norbitals, nbands, numorb_max, ME2c_max, ME3c_max, index_max2c
!$ volatile index_max3c, lssh, mu, mvalue, nssh, nu, num_orb, degelec
!$ volatile mu2shell, index_maxPP, lsshPP, muPP, nuPP, nsshPP, num_orbPP, imass 
!$ volatile cl_PP, h_mat, s_mat, sVNL, t_mat, vna, vnl, vxc, vxc_1c, dip
!$ volatile ewald, ewaldlr, ewaldsr, vca, vxc_ca, Vcoulomb, Vewaldsr, Vxcnu
!$ volatile ideriv_max, isorpmax
 
! Do we print out everything?
         logical wrtout

! Smoother parameters
         real, parameter :: smt_elect = 0.8d0 ! Ewald and electrostatic
         real, parameter :: smt_vnl = 0.8d0   ! VNL terms
! Total number of orbitals in the Hamiltonian and overlap matrices.
         integer norbitals

! Actual total number of orbitals in the Hamiltonian and overlap matrices.
         integer norbitals_new

! Total number of occupied bands in the solved problem.
         integer nbands

! Maximum number of shells in one atom (No effect on memory, can leave big)
! Examples: s ==> 1, sp^3 ==>  2, ss*p^3p*^3 ==>  4, sp^3d^5 ==>  3
!                                 ss*p^3p*^3d^5d*^5  ==>  6
! Cannot be greater than 8 (or else ind2c is wrong!)
! Calculated in fireball.f90
         integer, parameter :: nsh_max = 6

! Maximum number of shells and derivatives - needed for reading in the 
! data files.
         integer ideriv_max
         integer isorpmax
         integer isorpmax_xc

! Maximum number of orbitals in one atom.
! Examples: s ==> 1, sp^3 ==>  4, ss*p^3p*^3 ==>  8, sp^3d^5 ==>  9
         integer :: numorb_max

! Maximum number of two-center matrix elements: (calculated in make_munu.f90)
! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 39
         integer :: ME2c_max
         integer :: ME2cPP_max
         integer :: ME2cDipY_max
         integer :: ME2cDipX_max


! Maximum number of three-center matrix elements: (calculated in make_munu.f90)
! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
         integer :: ME3c_max

! Maximum number of two and three-center matrix elements in spherical density
! approximation (OLSXC) (calculated in make_munuS.f90)
! Examples: s ==> 1, sp^3 ==> 4, sp^3d^5 ==> 9
         integer :: MES_max

         integer, dimension (:, :), allocatable :: index_max2c
         integer, dimension (:, :), allocatable :: index_max3c
         integer, dimension (:, :), allocatable :: index_max2cDipY
         integer, dimension (:, :), allocatable :: index_max2cDipX
         integer, dimension (:, :), allocatable :: lssh
         integer, dimension (:, :, :), allocatable :: mu
         integer, dimension (:, :, :), allocatable :: mvalue
         integer, dimension (:), allocatable :: nssh
         integer  nssh_tot
         integer, dimension (:, :, :), allocatable :: nu
         integer, dimension (:), allocatable :: num_orb
         integer, dimension (:, :, :), allocatable :: muDipY
         integer, dimension (:, :, :), allocatable :: nuDipY
         integer, dimension (:, :, :), allocatable :: muDipX
         integer, dimension (:, :, :), allocatable :: nuDipX

! To get information of the orbitals, re loaded in getinfo_orbital
         integer, dimension (:), allocatable :: getmssh
         integer, dimension (:), allocatable :: getlssh
         integer, dimension (:), allocatable :: getissh
         integer, dimension (:), allocatable :: getiatom
 
! Placement of matrix element in Hamiltonian and overlap matrices.  
! The placement is based on which atom number is accessed and the number
! of orbitals of the atom.
         integer, dimension (:), allocatable :: degelec

! Needed for extended hubbard interactions.
         integer, dimension (:, :), allocatable :: mu2shell
 
! These variables are specifically for the Kleinmann-Bylander pseudo-potentials
         integer, dimension (:, :), allocatable :: index_maxPP
         integer, dimension (:, :), allocatable :: lsshPP
         integer, dimension (:, :, :), allocatable :: muPP
         integer, dimension (:, :, :), allocatable :: nuPP
         integer, dimension (:), allocatable :: nsshPP
         integer, dimension (:), allocatable :: num_orbPP


! These variables are specifically for spherical density approximation 
! used in OLSXC method
         integer, dimension (:, :), allocatable :: index_maxS
         integer, dimension (:, :, :), allocatable :: muS
         integer, dimension (:, :, :), allocatable :: nuS
         integer, dimension (:, :, :), allocatable :: mvalueS

! This array stores the species information for each atom.
         integer, dimension (:), allocatable :: imass
 
         real, dimension (:, :), allocatable :: cl_PP

! These arrays store the interactions for the Hamiltonian matrix. 
! These will be dimensioned according to (mu, nu, natoms, neigh_max),
! where mu, nu are the orbitals for iatom and its neighbor.
         real, dimension (:, :, :, :), allocatable :: h_mat
         real, dimension (:, :, :, :), allocatable :: s_mat
         real, dimension (:, :, :, :), allocatable :: sm_mat

         real, dimension (:, :, :, :), allocatable :: sVNL
         real, dimension (:, :, :, :), allocatable :: t_mat
         real, dimension (:, :, :, :), allocatable :: vna
         real, dimension (:, :, :, :), allocatable :: vnl
         real, dimension (:, :, :, :), allocatable :: vnl2c
         real, dimension (:, :, :, :), allocatable :: vnl3c
         real, dimension (:, :, :, :), allocatable :: vxc
         real, dimension (:, :, :, :), allocatable :: vxc_1c

! These arrays store interactions which are needed for the DOGS 
! contributions to the Hamiltonian matrix. 
         real, dimension (:, :, :, :), allocatable :: dip
         real, dimension (:, :), allocatable :: ewald
         real, dimension (:, :, :, :), allocatable :: ewaldlr
         real, dimension (:, :, :, :), allocatable :: ewaldsr
         real, dimension (:, :, :, :), allocatable :: vca
         real, dimension (:, :, :, :), allocatable :: vxc_ca
         real, dimension (:, :, :, :), allocatable :: ewaldqmmm
! Dipole with XYZ components
         real, dimension (:, :, :), allocatable :: dipcm
         real, dimension (:, :, :, :, :), allocatable :: dipc

! These arrays store interactions which are needed for the extended-Hubbard
! contributions to the Hamiltonian matrix. 
         real, dimension (:, :), allocatable :: Vcoulomb
         real, dimension (:, :), allocatable :: Vewaldsr
         real, dimension (:, :), allocatable :: Vxcnu

! van der Waals interactions
         real vdw
         real, dimension (:), allocatable :: C6
         real, dimension (:), allocatable :: p_alpha
         real, dimension (:), allocatable :: R0

! harmonic oscillator interactions
         real enHarmonic

! These arrays are for evaluating the gaussian approximation to the 
! three-center exchange-correlation contributions.
         real, dimension (:, :, :, :), allocatable :: bar_density_2c
         real, dimension (:, :, :, :), allocatable :: bar_density_3c
         real, dimension (:, :, :, :), allocatable :: density_2c
         real, dimension (:, :, :, :), allocatable :: density_3c
         real, dimension (:, :, :, :), allocatable :: nuxc_3c
         real, dimension (:, :, :, :), allocatable :: nuxc_total
         real, dimension (:, :, :, :), allocatable :: vxc_3c

! CGP
! These arrays for the dos calculation
         complex, dimension (:,:), allocatable   :: hamk

! Define for the gaussian
        real, parameter ::  xc_overtol = 5.0d-5
        real, parameter ::  xc_overtolG = 5.0d-5

! jel-grid
! These variables contains names of wavefunction and potential files
        character (len=25), dimension (:,:), allocatable :: wavefxn
        character (len=25), dimension (:,:), allocatable :: napot
! end jel-grid


! =============================================================================

        end module
