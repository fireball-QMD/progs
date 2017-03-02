        module kpoints
         
!$ volatile special_k,special_k_orig,weight_k,weight_k_orig,scale_k

! Special k-points
         integer nkpoints
         real, dimension (:, :), allocatable :: special_k
         real, dimension (:, :), allocatable :: special_k_orig
         real, dimension (:), allocatable :: weight_k
         real, dimension (:), allocatable :: weight_k_orig
         real, dimension (:, :), allocatable :: scale_k

! Input files         
         character (len = 40) kptpreference 

        end module
