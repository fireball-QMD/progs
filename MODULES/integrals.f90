        module integrals
         use dimensions 

!$ volatile icon3c, numx3c_bcna, numy3c_bcna, hx_bcna, hy_bcna, x3cmax_bcna
!$ volatile y3cmax_bcna, bcna_01, bcna_02, bcna_03, bcna_04, bcna_05
!$ volatile numx3c_xc3c, numy3c_xc3c, hx_xc3c, hy_xc3c, x3cmax_xc3c, y3cmax_xc3c
!$ volatile xc3c_01, xc3c_02, xc3c_03, xc3c_04, xc3c_05, ind2c, numz2c
!$ volatile xintegral_2c, splineint_2c, z2cmax, exc1c_0, exc1c, xcnu1c

! Fdata location
        character (len = 200) fdataLocation


! ****************************************************************************
! One-center integrals
! ****************************************************************************

! Skipping some atoms in the Fdata files
        integer nsup
        integer, dimension (nspec_max) :: nsu

! ****************************************************************************
! Three-center integrals
! ****************************************************************************
         integer, dimension (:, :, :), allocatable :: icon3c

! Neutral (charged) atom interactions; there are five bcna arrays - one for
! each theta
         integer, dimension (:,:), allocatable :: numx3c_bcna
         integer, dimension (:,:), allocatable :: numy3c_bcna

         real, dimension (:,:), allocatable :: hx_bcna
         real, dimension (:,:), allocatable :: hy_bcna
         real, dimension (:,:), allocatable :: x3cmax_bcna
         real, dimension (:,:), allocatable :: y3cmax_bcna

         real, dimension (:, :, :, :, :), allocatable :: bcna_01
         real, dimension (:, :, :, :, :), allocatable :: bcna_02
         real, dimension (:, :, :, :, :), allocatable :: bcna_03
         real, dimension (:, :, :, :, :), allocatable :: bcna_04
         real, dimension (:, :, :, :, :), allocatable :: bcna_05

! XC interactions; 7 implies different derivative types; there are five xc3c
! arrays - one for each theta
         integer, dimension (:,:), allocatable :: numx3c_xc3c
         integer, dimension (:,:), allocatable :: numy3c_xc3c

         real, dimension (:,:), allocatable :: hx_xc3c
         real, dimension (:,:), allocatable :: hy_xc3c
         real, dimension (:,:), allocatable :: x3cmax_xc3c
         real, dimension (:,:), allocatable :: y3cmax_xc3c

         real, dimension (:, :, :, :, :), allocatable :: xc3c_01
         real, dimension (:, :, :, :, :), allocatable :: xc3c_02
         real, dimension (:, :, :, :, :), allocatable :: xc3c_03
         real, dimension (:, :, :, :, :), allocatable :: xc3c_04
         real, dimension (:, :, :, :, :), allocatable :: xc3c_05
!xc3c_SN
! XC interactions; 7 implies different derivative types; there are five den3
! arrays - one for each theta, used only for SNXC method
         integer, dimension (:,:), allocatable :: numx3c_den3
         integer, dimension (:,:), allocatable :: numy3c_den3

         real, dimension (:,:), allocatable :: hx_den3
         real, dimension (:,:), allocatable :: hy_den3
         real, dimension (:,:), allocatable :: x3cmax_den3
         real, dimension (:,:), allocatable :: y3cmax_den3

         real, dimension (:, :, :, :, :), allocatable :: den3_01
         real, dimension (:, :, :, :, :), allocatable :: den3_02
         real, dimension (:, :, :, :, :), allocatable :: den3_03
         real, dimension (:, :, :, :, :), allocatable :: den3_04
         real, dimension (:, :, :, :, :), allocatable :: den3_05

         real, dimension (:, :, :, :, :), allocatable :: den3S_01
         real, dimension (:, :, :, :, :), allocatable :: den3S_02
         real, dimension (:, :, :, :, :), allocatable :: den3S_03
         real, dimension (:, :, :, :, :), allocatable :: den3S_04
         real, dimension (:, :, :, :, :), allocatable :: den3S_05

!end xc3c_SN

! ****************************************************************************
! Two center integrals
! ****************************************************************************
! jel-F2c
!         integer, dimension (1:21, 0:8) :: ind2c
         integer, dimension (1:23, 0:8) :: ind2c
! end jel-F2c
         integer, dimension (:, :, :), allocatable :: numz2c

         real, dimension (:, :, :, :, :), allocatable :: xintegral_2c
         real, dimension (:, :, :, :, :, :), allocatable :: splineint_2c
         real, dimension (:, :, :), allocatable :: z2cmax

! One center integrals
         real, dimension (:, :), allocatable :: exc1c_0
         real, dimension (:, :, :, :), allocatable :: exc1c
         real, dimension (:, :, :), allocatable :: xcnu1c
         real, dimension (:, :, :), allocatable :: xcnu1cs
! jel-der
!         real, dimension (: ,), allocatable :: exc1c0
         real, dimension (:, :, :), allocatable :: exc1c0
         real, dimension (:, :, :), allocatable :: nuxc1c
         real, dimension (:, :, :, :), allocatable :: dexc1c
         real, dimension (:, :, :), allocatable :: d2exc1c
         real, dimension (:, :, :, :), allocatable :: dnuxc1c
         real, dimension (:, :, :, :, :), allocatable :: d2nuxc1c
! end jel-der


        end module
