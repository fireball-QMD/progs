! allocate_sys.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine allocate_trans (nkpoints)
   
   use dimensions
   use transport
   use interactions
   use neighbor_map

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in)   ::  nkpoints

! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

   integer numorb1
   integer numorb2
   integer inu
   integer ncell1
   integer ncell2

   complex a0
   complex a1

! Procedure
! ===========================================================================

   
   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to allocate_trans subroutine! '
   write (*,*) ' Allocate arrays of systems for transport calculation '

   a0 = (0.0d0, 0.0d0)
   a1 = (1.0d0, 0.0d0)

   numorb1 = sample1%norb
   numorb2 = sample2%norb

! Hamiltonian of sample_x
   allocate (Hsam1_k (numorb1,numorb1))
   allocate (Hsam2_k (numorb2,numorb2))


! Global Hamiltonian in kspace
   allocate (H_k (norbitals,norbitals,nkpoints))
   H_k = (0.0d0, 0.0d0)

! considering only xy-plane therefore **2
   ncell1 = (2*sample1%ncell+1)**2
   ncell2 = (2*sample2%ncell+1)**2
   numorb1 = sample1%norb_tip
   numorb2 = sample2%norb_tip
!   numorb1 = sample1%norb_tip*ncell1
!   numorb2 = sample2%norb_tip*ncell2


! Green's function of tip
! retarded G.f.
  allocate (Gr_tip1 (numorb1,numorb1,nE))
  allocate (Gr_tip2 (numorb2,numorb2,nE))
! advanced G.f.
  allocate (Ga_tip1 (numorb1,numorb1,nE))
  allocate (Ga_tip2 (numorb2,numorb2,nE))

! complex current matrix
  allocate (Jc (numorb1,numorb1,nE))

! mapping matrices 
  allocate (t_12 (numorb1,numorb2))
  allocate (t_21 (numorb2,numorb1))

  Gr_tip1 = a0 
  Gr_tip2 = a0 
  Ga_tip1 = a0 
  Ga_tip2 = a0 
  t_12 = a0
  t_21 = a0
  Jc = a0


! Format Statements
! ===========================================================================

   return
   
 end subroutine allocate_trans
