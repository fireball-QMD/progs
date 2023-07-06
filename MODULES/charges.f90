        module charges
         use dimensions

!$ volatile nelectron, nzx, Qin, Qinmixer, Qneutral, Qout, Qoutmixer, dq
!$ volatile QLowdin_TOT, QMulliken_TOT

         real    efermi
         real ztot
         real qstate

         integer, dimension (:), allocatable :: nelectron
         integer, dimension (:), allocatable :: nzx

         real, dimension (:, :), allocatable :: Qin
         real, dimension (:), allocatable :: Qinmixer
         real, dimension (:, :), allocatable :: Qneutral
         real, dimension (:, :), allocatable :: Qout
         real, dimension (:), allocatable :: Qoutmixer
         real, dimension (:), allocatable :: dq

         real, dimension (:), allocatable :: QLowdin_TOT
         real, dimension (:), allocatable :: QMulliken_TOT

! vlada-cdft-mdet
         real, dimension (:, :), allocatable :: Qin_es
         real, dimension (:, :), allocatable :: Qout_es      !vlada cdft
         real, dimension (:), allocatable :: QLowdin_TOT_es  !vlada cdft
         real, dimension(:), allocatable  :: Q0_TOT
         real, dimension(:), allocatable  :: Q_partial


! anderson iteration procedure
        real, allocatable, dimension(:,:) :: Fv   ! x_try-x_old 
        real, allocatable, dimension(:,:) :: Xv   ! x_old 
        real, allocatable, dimension(:,:) :: delF! F(m+1)-F(m) 
        real, allocatable, dimension(:,:) :: delX! X(m+1)-X(m) 
        real, allocatable, dimension(:)   :: r2_sav
! Broyden mixing
        real, allocatable, dimension(:,:) :: RJac
! Louie mixing
	real, allocatable, dimension(:,:) :: betaInvH
	real, allocatable, dimension(:,:) :: gamaH
! NPA
        real, allocatable, dimension(:,:) :: qaux
! Kohn-Sham
        real, dimension (:, :), allocatable :: Qxneutral
        real, dimension (:), allocatable :: rho_in
        real, dimension (:), allocatable :: rho_out
        real, dimension (:), allocatable :: mwe
        real, dimension (:), allocatable :: drwe
! Populations analysis (iwrtpop = 1). energy range
        real  Epop_L
        real  Epop_U
!  Dipoles for iqout = 7
        real    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x, dipRes_x
        real    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y, dipRes_y
        real    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z, dipRes_z
        real    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot, dipRes_tot 
        real, dimension (:), allocatable :: dq_DP
        end module
