module configuration
!$ volatile ratom, ximage, xdot, xl

       integer :: natoms
       integer :: nspecies

       real, dimension (:, :), allocatable :: ratom
       real, dimension (:, :), allocatable :: ximage
       real, dimension (:, :, :), allocatable :: xdot
       real, dimension (:, :), allocatable :: nowMinusInitialPos
       real, dimension (:, :), allocatable :: initialPosition
!         real, dimension (:, :), allocatable :: original_shift
! PP change
!        real, dimension (3, 0:124) :: xl
!        real, dimension (3, 0:342) :: xl
!        integer, parameter :: mbeta_max = 728
!!         real, dimension (3, 0:728) :: xl
! now we do neighbors shells dynamically    jel(2/2/2005)
       integer :: mbeta_max
       real, dimension (:,:), allocatable :: xl

! atomixc mass
	   real, dimension (:), allocatable :: xmass
	   
! Stuff for shifting atoms - if any are at the origin-->problem 
! later.  Shift all atoms away from origin.
        integer ishiftO                ! Shift atoms by constant amount
        real, dimension (3) :: shifter 
! Lattice vectors
        real, dimension (3) :: a1vec
        real, dimension (3) :: a2vec
        real, dimension (3) :: a3vec
! Volume of unit cell
        real Vouc

! d@ni rescal = 1.0d0 default option, to rescalate de positions, lvs and kpts, 
! to do a volumen calculations, bulkmodulus etc ...
        real :: rescal
        integer :: xyz2line
! Positions, velocities of the atoms
        real, dimension (:, :), allocatable :: vatom
        character (len = 2), dimension (:), allocatable :: symbolA
        character (len = 2), dimension (:), allocatable:: symbol
! These declarations are for the charge transfer contributions:
        real, dimension (:), allocatable :: etotatom
        real, dimension (:,:), allocatable :: rcutoff
        real, dimension (:), allocatable :: rc_PP
!        real, dimension (:), allocatable :: xmass
        real, dimension (:), allocatable :: smass

! ---------------------------------------------------------------------------
! Input files 
! ---------------------------------------------------------------------------
        character (len = 40) basisfile            
        character (len = 40) lvsfile  
        character (len = 40), parameter :: initfile = 'fireball.in'

              
end module
