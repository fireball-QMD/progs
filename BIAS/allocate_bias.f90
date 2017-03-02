! allocate_sys.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine allocate_bias ( natoms )
   
   use dimensions
   use interactions
   use bias
   use neighbor_map

! Argument Declaration and Description
! ===========================================================================
! Input
  integer, intent (in)             :: natoms

! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

   
   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to allocate_bias subroutine! '
   write (*,*) ' Allocate arrays of systems for the option ibias  '


   allocate (Vbias_mat (numorb_max, numorb_max, neigh_max, natoms))
   allocate (Vbiasp_mat (3,numorb_max, numorb_max, neigh_max, natoms))

   Vbias_mat = 0.0d0
   Vbiassp_mat = 0.0d0

! Format Statements
! ===========================================================================

   return
   
 end subroutine allocate_bias
