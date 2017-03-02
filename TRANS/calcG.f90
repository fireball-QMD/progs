! calcG.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine calcG ()
   
   use kpoints
   use interactions
   use transport
   use neighbor_map
   use configuration

! Argument Declaration and Description
! ===========================================================================
! Input


! Output


! Local Parameters and Data Declaration
! ===========================================================================
! NOTE: unit conversion from quantum conductance unit to A/V (or Siemens)
  real, parameter :: go2siemens = 7.7480917330d-05
 
! Local Variable Declaration and Description
! ===========================================================================

   integer ikpoints
   real, dimension (3) :: k_temp
   
! Procedure
! ===========================================================================

    write (*,*)
    write (*,*)  ' Calculate Transport ...'
    write (*,*)

! calculate t_12 matrix from dimer approximation
!    call assemble_t12_bare ( )
    call assemble_t12_fit ( )


    do ikpoint = 1, nkpoints

      write (*,*) '  '
      write (*,*) ' ****************************************************** '
      write (*,*) '  '
      write (*,*) '         Doing ikpoint = ', ikpoint
      write (*,*) '  '
      write (*,*) ' ****************************************************** '

      k_temp(:) = special_k(:,ikpoint)

! Assemble Hamiltonians of samples
      call assemble_Hsam (k_temp, ikpoint)
! Assemble Green func
      call assemble_Gsam (ikpoint)
  
    enddo ! enddo ikpoint

! write DOS
!     call wrt_dostip (iwrtout)

! calc current
    call assemble_Dxx ()

   
! Format Statements
! ===========================================================================
100     format (2i5,2f14.6)

   return
   
 end subroutine calcG
