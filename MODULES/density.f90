        module density

!$ volatile bbnkre, bbnkim, blowim, blowre, cape, rho

! These arrays contain the coefficients to calculate the density matrix
! and the Lowdin charges:
         real, dimension (:, :, :), allocatable :: bbnkre
         real, dimension (:, :, :), allocatable :: bbnkim
         real, dimension (:, :, :), allocatable :: blowim
         real, dimension (:, :, :), allocatable :: blowre
         real, dimension (:, :), allocatable :: eigen_k

! These arrays store the densities.
         real, dimension (:, :, :, :), allocatable :: cape
         real, dimension (:, :, :, :), allocatable :: rho
         real, dimension (:, :, :, :), allocatable :: rhoPP
! jel-grid
         real, dimension (:, :), allocatable :: rhoA
         real, dimension (:, :, :, :), allocatable :: rho_old
! end jel-grid

! These arrays store stuff related to the average density.
! Used in the OLSXC exchange-correlation interactions.
         real, dimension (:, :, :), allocatable :: rho_on
         real, dimension (:, :, :), allocatable :: rhoi_on
         real, dimension (:, :, :, :), allocatable :: rho_off
         real, dimension (:, :, :, :), allocatable :: rhoij_off
         real, dimension (:, :, :), allocatable :: arho_on
         real, dimension (:, :, :), allocatable :: arhoi_on
         real, dimension (:, :, :, :), allocatable :: arho_off
         real, dimension (:, :, :, :), allocatable :: arhoij_off
! JOM-nonadiabatic
         real, dimension ( :, :), allocatable :: foccupy_na
         integer, dimension ( :, :), allocatable :: ioccupy_na
! VLADA-nonadiabatic
         real, dimension ( :, :), allocatable :: foccupy_na_TS    
         integer, dimension ( :, :), allocatable :: ioccupy_na_TS   
         real, dimension (:, :, :), allocatable :: bbnkre_o 
         real, dimension (:, :, :), allocatable :: blowre_o  
       end module density
