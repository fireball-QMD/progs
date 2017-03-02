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

        end module
