        module nonadiabatic
! this module contains the variables requiered for the calculation of
! non-adiabatic couplings. 
!$ volatile gover, gover1c, gH_2c, gh_atm, gh_3c
!$ volatile gh_pp_2c, gh_pp_atm, gh_pp_3c, gks

 
! Matrix elements involving the Gradients
         real, dimension (:, :, :, :, :), allocatable :: gover
         real, dimension (:, :, :),    allocatable :: gover1c
         real, dimension (:, :, :),    allocatable :: f1nac1c
         real, dimension (:, :, :),    allocatable :: f2nac1c
         real, dimension (:, :, :, :, :), allocatable :: gh_2c 
!        real, dimension (:, :, :, :, :), allocatable :: gh_2c_ca    ! smear large arrays VLADA 
         real, dimension (:, :, :, :, :), allocatable :: gh_atm
!        real, dimension (:, :, :, :, :), allocatable :: gh_atm_ca   ! smear large arrays VLADA
         real, dimension (:,:, :, :, :, :), allocatable :: gh_3c
!        real, dimension (:,:, :, :, :, :), allocatable :: gh_3c_ca  ! smear large arrays VLADA
!        real, dimension (:,:, :, :, :, :), allocatable :: gh_lrew   ! smear large arrays VLADA
!        real, dimension (:,:, :, :, :, :), allocatable :: gh_xc_3c  ! smear large arrays VLADA
         real, dimension (:, :, :, :, :), allocatable :: gh_pp_otl
         real, dimension (:, :, :, :, :), allocatable :: gh_pp_otr
         real, dimension (:, :, :, :, :), allocatable :: gh_pp_atm
         real, dimension (:,:, :, :, :, :), allocatable :: gh_lrew_qmmm
         real, dimension (:,:, :, :, :, :), allocatable :: gh_pp_3c

         real, dimension (:, :, :, :),    allocatable :: gks
! coefficients for the TD-wfs
         complex, dimension (:, :, :),    allocatable :: c_na
! variables for the integration of TD-wfs
         integer nddt
         real, dimension (:, :), allocatable :: ratom_old
         real, dimension (:, :), allocatable :: vatom_old
         real, dimension (:, :), allocatable :: eigen_old
         real, dimension (:, :), allocatable :: eigen_0
         real, dimension (:, :), allocatable :: eigen_1
         real, dimension (:, :, :, :),    allocatable :: gks_old
         real, dimension (:, :, :, :),    allocatable :: gks_1
         real, dimension (:, :, :, :),    allocatable :: gks_0
         real, dimension (:, :, :),    allocatable :: bbnkre_old
         real, dimension (:, :, :),    allocatable :: blowre_old
! maps of lists of evolving KS-states:
         integer nele                          ! number of electrons
         integer nalph
!        integer, dimension (: ),    allocatable :: map_td
         integer, dimension (: ),    allocatable :: map_ks
         integer, dimension (: ),    allocatable :: map_proj
! numeric derivative (g*v)
         real, dimension (:,:), allocatable :: dnac
         real, dimension (:,:), allocatable :: dnac_old
         real, dimension (:,:), allocatable :: sumb
!jel-nac
	 real, dimension (:,:), allocatable :: ratom_opt
	 real, dimension (:,:), allocatable :: dnac_opt

         integer trans ! VLADA

        end module
