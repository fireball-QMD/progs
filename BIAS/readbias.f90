! readbias.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine readbias ( natoms )

   use interactions
   use neighbor_map
   use bias
   implicit none

! Argument Declaration and Description
! ===========================================================================
  integer, intent (in)             :: natoms


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

   integer i
   integer j
   integer iatom
   integer jatom
   integer nint
   integer imu
   integer inu
   integer index

   real delta

! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to readbias subtoutine! '
   write (*,*) ' Read in the data from the scriptfile - bias.optional '

   open (unit = 12, file = 'bias.optional', status = 'old')

   write (*,*) ' Input applied bias voltage [eV]: '
   read (12,*) Vbias
   write (*,600) Vbias
   write (*,*) ' Input z-axis value zb1, upper lead (-eV/2) :  '
   read (12,*) zb1
   write (*,*) ' Input z-axis value zb0, lower lead (eV/2) :  '
   read (12,*) zb0
   write (*,700) zb1,zb0

! close the input file
   close (unit = 12)

   return

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a30)
600     format (2x,' Vbias            = ',f16.8,' [eV]')
700     format (2x,' zb1  = ',f16.8,' zb0   = ',f16.8)

 end subroutine readbias
