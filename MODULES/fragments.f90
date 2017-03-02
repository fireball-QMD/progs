        module fragments

!$ volatile numfrags,fragtemp,fraggots,ratom_frag,ratom_frag_save,fragsize,fragatm
! jel-fr
         integer :: ifrags                     ! type of fragments 
! end jel-fr
! dani.JOM
         integer :: nfragments                 ! number of fixed atoms
         integer :: numfrags                   ! number of framents
         integer :: fragtemp                   ! do we project forces or not?
         integer, dimension(:), allocatable :: fraggots ! what atoms are in fragments
         integer, dimension(:,:), allocatable :: fragxyz ! which axis-direction is fixed
         real, dimension (:,:), allocatable :: ratom_frag ! work space
         real, dimension (:,:), allocatable :: ratom_frag_save ! first geometry
         integer, dimension (:), allocatable :: fragsize ! size of fragment
         integer, dimension (: ,:), allocatable :: fragatm ! atom i of frag j is this atom in the cell

        end module
