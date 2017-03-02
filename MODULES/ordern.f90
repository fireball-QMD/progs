        module ordern
! FIXME replace these volatile attributes with threadprivate when supported
!$ volatile mpi_whatever_real, nactualprocs, nhmax, ncmax, nctmax, nFmax
!$ volatile mpi_whatever_double
!$ volatile nFtmax, ncrowsmax, nprowsmax, numh, listh, h_s_compact, h_compact
!$ volatile s_compact, nums12_local, lists12_local, s12_compact_local
!$ volatile numcape_local, listcape_local, cape_compact_local, numrho_local
!$ volatile listrho_local, rho_compact_local, numc_local, listc_local
!$ volatile c_compact_local, numG_local, listG_local, G_compact_local
!$ volatile numRG_local, listRG_local, RG_compact_local
!$ volatile hs_bit_matrix, Xmax
!$ volatile numX_local, listX_local, X_compact_local

! Maximum number of ordern iterations and acceptable relative error tolerance.
         integer, parameter :: max_ordern_iterations = 300
         real, parameter :: ordern_tolerance = 1.0d-6
         real, parameter :: ordern_grad_tolerance = 5.0d-4
! If the energy is less than this then we call it a blowup and restart.
         real, parameter :: ebs_blowup = 1.0d+6

! Minimum stepsize (OK, actually maximum, but you get the idea.)
         real, parameter :: min_stepsize = -0.000001d0

! Actual number of processors to be used. 
         integer nactualprocs

! Maximum number of non-zero matrix elements in the Hamiltonian and overlap.
         integer nhmax

! Maximum number of non-zero matrix elements in the operator S^-1/2.
         integer nS12max

! Maximum number of non-zero matrix elements for the wavefunctions coefficients
! (matrix and transpose) 
         integer ncmax
         integer nctmax
         integer Xmax

! Maximum number of non-zero matrix elements for the F and Fs matrices.
! (matrix and transpose) 
         integer nFmax
         integer nFtmax

! Maximum number of rows in the sparse matrices.
         integer ncrowsmax
         integer nprowsmax

! Maximum cutoff radius for the LWF's
         real, parameter :: rcutoff_lwf = 5.5d0

! ****************************************************************************
! The variable arrays are organized as follows: 
! numA : Control vector of A matrix
!        (number of nonzero elements of each row of A)
! listA : Control vector of A matrix
!         (list of nonzero elements of each row of A)
! A_compact: Contains actual nonzero values of the A matrix
! ****************************************************************************
! Representation of the Hamiltonian and overlap matrices in a compact format.
         integer, dimension (:), allocatable :: numh
         integer, dimension (:, :), allocatable :: listh

         real, dimension (:, :, :), allocatable, target :: h_s_compact
         real, dimension (:, :), pointer :: h_compact
         real, dimension (:, :), pointer :: s_compact

! Representation of S^-1/2 for Lowdin transformation in a compact form. 
         integer, dimension (:), allocatable :: nums12_local
         integer, dimension (:, :), allocatable :: lists12_local

         real, dimension (:, :), allocatable :: s12_compact_local

! Representation of the density matrices in a compact form.  The indexing of
! these are the same as the indexing for h_compact.
         integer, dimension (:), allocatable :: numcape_local
         integer, dimension (:), allocatable :: numrho_local
         integer, dimension (:, :), allocatable :: listcape_local
         integer, dimension (:, :), allocatable :: listrho_local

         real, dimension (:, :, :), allocatable, target ::                   &
     &    cape_rho_compact_local
         real, dimension (:, :), pointer :: cape_compact_local
         real, dimension (:, :), pointer :: rho_compact_local

! Representation of the wavefunction coefficients in a compact format.
         integer, dimension (:), allocatable :: numc_local
         integer, dimension (:), allocatable :: numct_local
         integer, dimension (:, :), allocatable :: listc_local
         integer, dimension (:, :), allocatable :: listct_local

         real, dimension (:, :), allocatable, target :: c_compact_local
         real, dimension (:, :), allocatable :: ct_compact_local

! Representation of the X matrix in a compact format.
         integer, dimension (:), allocatable, target :: numX_local
         integer, dimension (:, :), allocatable, target :: listX_local

         real, dimension (:, :), allocatable, target :: X_compact_local

! Representation of the gradient in a compact format.
         integer, dimension (:), allocatable, target :: numG_local
         integer, dimension (:, :), allocatable, target :: listG_local

         real, dimension (:, :), allocatable, target :: G_compact_local

! ****************************************************************************
! Step directions (R) matrix declaration - local pieces.
! A step direction is determined for the CG methods.
! ****************************************************************************
         integer, dimension (:), allocatable, target :: numRG_local
         integer, dimension (:, :), allocatable, target :: listRG_local

         real, dimension (:, :), allocatable, target :: RG_compact_local

! Bit matrix indicating which blocks of the h/s matrices are nonzero.
! FIXME: Make this into a packed bit vector to save space.
         logical, dimension (:, :), allocatable :: hs_bit_matrix

        end module
