! wrt_dos.f90
! Program Description
! ===========================================================================
! subroutine writes LDOS and DOS of tips in output files  
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine wrt_dostip (iwrtout)
   
   use dimensions
   use interactions
   use structure
   use transport
   use kpoints

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   logical, intent (in)   ::  iwrtout

! Output


! Local Parameters and Data Declaration
! ===========================================================================

   real, parameter :: pi = 3.141592653589793238462643

! Local Variable Declaration and Description
! ===========================================================================

   integer ie
   integer inu
   integer jnu
   integer trace
   integer trace1
   integer iatom
   integer in1
   integer itip

   real energy
   real, allocatable, dimension (:) :: dos1
   real, allocatable, dimension (:) :: dos2

   real dos1_tot
   real dos2_tot

   character (3)    :: ind_aux

! Procedure
! ===========================================================================
 

! dimension
   norb1 = sample1%norb_tip
   norb2 = sample2%norb_tip

   allocate ( dos1(norb1) )
   allocate ( dos2(norb2) )

   dos1 = 0.0d0
   dos2 = 0.0d0
   dos1_tot = 0.0d0
   dos2_tot = 0.0d0

   write (*,*) '  Write  atomic DOS into files  '
   open ( unit = 18, file = 'dos1.dat', status = 'unknown')
   open ( unit = 19, file = 'dos2.dat', status = 'unknown')

! set initial energy
   energy = Elow

! Loop of energy
   do ie = 1,nE

! SAMPLE_1
      dos1_tot = 0.0d0
! Lopp over atoms
      do itip = 1, sample1%natom_tip

         iatom = sample1%atom_tip(itip)
         in1 = imass(iatom) 
! Assign proper extension to each data file
         write (ind_aux,'(i3.3)') iatom
! open file
         open ( unit = 17, file = 'dos1_'//ind_aux//'.dat', status = 'unknown',position = 'append')

         trace = pointer1(itip)
         trace1 =trace + 1

! calc dos
         do jnu = 1,num_orb(in1)
            inu = trace + jnu 
!            if (ie .eq. 1) then 
            dos1(jnu) = (-1.0d0/PI)*Imag(Gr_tip1(inu,inu,ie)) 
            dos1_tot = dos1_tot + dos1(jnu)
!            else
!               dos1(jnu) = (-1.0d0/PI)*Imag(Gr_tip1(inu,inu,ie)) +         &
!   &            eta*(1.0d0/PI)*( real(Gr_tip1(inu,inu,ie))                 &
!   &            - real(Gr_tip1(inu,inu,ie-1)) )/real(dE)
!            endif
         enddo
         write (17,100) energy, dos1(1:num_orb(in1)),sum(dos1(1:num_orb(in1)))
! close DOS file of given atom
         close (17)
      enddo ! do iatom
      write (18,200) energy, dos1_tot

! SAMPLE_2
      dos2_tot = 0.0d0
! Lopp over atoms
      do itip = 1, sample2%natom_tip

         iatom = sample2%atom_tip(itip)
         in1 = imass(iatom) 
! Assign proper extension to each data file
         write (ind_aux,'(i3.3)') iatom
! open file
         open ( unit = 17, file = 'dos2_'//ind_aux//'.dat', status = 'unknown',position = 'append')
         trace = pointer2(itip)
         trace1 =trace + 1
! calc dos
         do jnu = 1,num_orb(in1)
            inu = trace + jnu 
!            if (ie .eq. 1) then 
            dos2(jnu) = (-1.0d0/PI)*Imag(Gr_tip2(inu,inu,ie)) 
            dos2_tot = dos2_tot + dos2(jnu)
!            dos2(jnu) = (-1.0d0/PI)*Imag(Gr_tip2(inu,inu,ie)) +            &
!   &            eta*(1.0d0/PI)*( real(Gr_tip2(inu,inu,ie))                 &
!   &            - real(Gr_tip2(inu,inu,ie-1)) )/real(dE)
!            endif
         enddo
         write (17,100) energy, dos2(1:num_orb(in1)),sum(dos2(1:num_orb(in1)))
! close DOS file of given atom
         close (17)
      enddo ! do iatom
      write (19,200) energy, dos2_tot
! increment the energy 
      energy = energy + dE
   enddo ! do ie

! close file total DOS
   close (18)
   close (19)

! deallocate
   deallocate (dos1)
   deallocate (dos2)

! Format Statements
! ===========================================================================
100     format ( f12.5, 9f10.4, f10.4 )
200     format ( f12.5, f14.6 )
300     format ( 2i5, 2f14.6 )

   return
 end subroutine wrt_dostip
