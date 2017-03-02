!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Module to define the variables we use in the hartree-fock!!!
!!! stuff                                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module hartree_fock

    real, dimension(:,:,:), allocatable     :: Uisigma
    real, dimension(:,:,:,:), allocatable   :: Jijsigma
    real, dimension(:,:,:,:), allocatable   :: hf_mat
    real, dimension(:,:), allocatable       :: hs_mat
    real, dimension(:,:), allocatable       :: Jialpha
    complex, dimension (:,:,:,:), allocatable :: nij
    integer                                 :: natomshf
    integer                                 :: natomhf_beg
    integer                                 :: natomhf_end
    real                                    :: betha

    end module


