         module neighbor_map
         use dimensions

!$ volatile neigh_max, neigh_max_vdw, neigh_b, neigh_b_vdw, neigh_j
!$ volatile neigh_j_vdw, neighn, neighn_vdw, neigh_comb, neigh_comj, neigh_comm
!$ volatile neigh_comn, neigh_back, neigh_self, neigh_vdw_self, range_vdw

! Maximum number of neighbors 
         integer neigh_max 
! Maximum number of neighbors (Pseudopotential)
! Only for simple overlap <phi_i|Psi_j>
         integer neighPP_max

! Maximum number of neighbors - for van der Waals interactions. 
         integer neigh_max_vdw 

! Neighbor mapping and storage.
         integer, dimension (:, :), allocatable :: neigh_b
         integer, dimension (:, :), allocatable :: neigh_j
         integer, dimension (:), allocatable :: neighn

! Neighbor (Pseudopotential) mapping and storage.
! nPP
         integer, dimension (:, :), allocatable :: nPP_b
         integer, dimension (:, :), allocatable :: nPP_j
         integer, dimension (:, :), allocatable :: nPP_map
         integer, dimension (:), allocatable :: nPPn
         integer, dimension (:), allocatable :: nPP_self
! nPPx
         integer, dimension (:, :), allocatable :: nPPx_b
         integer, dimension (:, :), allocatable :: nPPx_j 
         integer, dimension (:, :), allocatable :: nPPx_map
         integer, dimension (:, :), allocatable :: nPPx_point
         integer, dimension (:), allocatable :: nPPxn
         integer, dimension (:), allocatable :: nPPx_self
! neighPP
         integer, dimension (:, :), allocatable :: neighPP_b
         integer, dimension (:, :), allocatable :: neighPP_j
         integer, dimension (:), allocatable :: neighPPn

! Neighbor mapping and storage - for van der Waals interactions.
         integer, dimension (:, :), allocatable :: neigh_b_vdw
         integer, dimension (:, :), allocatable :: neigh_j_vdw
         integer, dimension (:), allocatable :: neighn_vdw

! Common neighbor mapping and storage.
         integer, dimension (:, :, :), allocatable :: neigh_comb 
         integer, dimension (:, :, :), allocatable :: neigh_comj
         integer, dimension (:, :, :), allocatable :: neigh_com_ng
         integer, dimension (:, :), allocatable :: neigh_comm 
         integer, dimension (:), allocatable :: neigh_comn 

! Common neighbor (Pseudopotential) mapping and storage.
! 3. party PP  common pairs
         integer, dimension (:, :, :), allocatable :: neighPP_comb 
         integer, dimension (:, :, :), allocatable :: neighPP_comj 
         integer, dimension (:, :), allocatable :: neighPP_comm 
         integer, dimension (:), allocatable :: neighPP_comn 

! Back neighbor mapping 
         integer, dimension (:,:), allocatable :: neigh_back

! neigh_self is the m value for the "self m"
         integer, dimension (:), allocatable :: neigh_self
         integer, dimension (:), allocatable :: neighPP_self
         integer, dimension (:), allocatable :: neigh_vdw_self

! total neighbor mapping (mapping neigh and neighPP together)
         integer, dimension(:), allocatable      :: neighn_tot
         integer, dimension(:,:), allocatable    :: neighj_tot
         integer, dimension(:,:), allocatable    :: neighb_tot
         integer                                 :: num_neig_maxtot

!CHROM neighbor for classical MD simulation
		 integer :: neigh_max_class = 0 !size of arrays down
         integer, dimension(:,:),   allocatable     ::  neigh_classic
         integer, dimension(:),     allocatable     ::  neighn_classic
         integer, dimension(:,:),   allocatable     ::  neigh_b_classic
!END CHROM                     

!SYMMETRIC FIREBALL
         integer :: tot_pairs
         integer, dimension(:), allocatable :: neigh_pair_a1
         integer, dimension(:), allocatable :: neigh_pair_a2
         integer, dimension(:), allocatable :: neigh_pair_n1
         integer, dimension(:), allocatable :: neigh_pair_n2 
!END SYMMETRIC FIREBALL


! Cutoff range for vdw interactions
         real range_vdw
        end module
