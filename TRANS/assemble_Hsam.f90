! assemble_Hsam.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_Hsam ( kpt, ikpt)
   
   use dimensions
   use interactions
   use transport
   use neighbor_map
   use configuration

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   real, intent (in), dimension (3)  :: kpt
   integer, intent (in)              :: ikpt

! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

   integer iatom
   integer jatom
   integer ineigh
   integer in1
   integer in2
   integer imu
   integer inu
   integer jmu
   integer jnu
   integer i
   integer j
   integer itip
   integer jtip
   integer iorb
   integer jorb
   integer ignu
   integer igmu
   integer mbeta
   
   real dot
   real, dimension (3) :: vec

   complex a0
   complex a1
   complex phase


! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) '  '

! set aux variables
   a0 = cmplx(0.0d0,0.0d0)
   a1 = cmplx(1.0d0,0.0d0)

   write (*,*) ' Assemble Hamiltonians of samples  ....'

! Now, create matrices H_1; H_2 of the samples H1 and H2.  

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          SAMPLE_1
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! first, copy unperturbed Hamiltonian htemp into HsamX_k      


   inu = 0
   do i = 1,sample1%natom
     iatom = sample1%atom(i)
     in1 = imass(iatom)
     imu = 0
     do j = 1,sample1%natom
       jatom = sample1%atom(j)
       in2 = imass(jatom)
       do iorb = 1, num_orb(in1)
         ignu = iorb + degelec(iatom)
         jnu = inu + iorb
         do jorb = 1, num_orb(in2)
           igmu = jorb + degelec(jatom)
           jmu = imu + jorb
           Hsam1_k(jnu,jmu) = H_k(ignu,igmu,ikpt)
         enddo
       enddo
       imu = imu + num_orb(in2)
     enddo ! do j
     inu = inu + num_orb(in1)
   enddo ! do i

!! loop over atoms in tip1
!   i = 0
!   do itip = 1,sample1%natom_tip
!      
!      iatom = sample1%atom_tip(itip)
!      in1 = imass(iatom)
!      j = 0
!! loop over atoms in tip2
!      do jtip = 1,sample2%natom_tip
!         
!         jatom = sample2%atom_tip(jtip)
!         in2 = imass(jatom)
!         
!! vec = r_i - r_j; we consider only interactions in the basic unit cell 
!!         vec(:) = ratom(:,jatom) - ratom(:,iatom)
!!@         vec(:) = 0.0d0
!         vec(:) = ratom(:,jatom) - ratom(:,iatom)
!         dot = kpt(1)*vec(1) + kpt(2)*vec(2) + kpt(3)*vec(3)
!         phase = cmplx(cos(dot),sin(dot))
!! sample_1
!! store piece of hopping between iatom and jatom       
!         do imu = 1, num_orb(in1)
!            jmu = imu + degelec1(itip)
!            do inu = 1, num_orb(in2)
!               jnu = inu + degelec2(jtip)
!               Hsam1_k(jmu,jnu) = Hsam1_k(jmu,jnu) - phase*t_12_orig(i+imu,j+inu)
!               Hsam1_k(jnu,jmu) = Hsam1_k(jnu,jmu) - phase*t_21_orig(j+inu,i+imu)
!            end do ! do inu
!         end do ! do imu
!         j = j + num_orb(in2)
!      enddo ! enddo jtip
!      i = i + num_orb(in1)
!   enddo ! enddo itip


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          SAMPLE_2
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! first, copy unperturbed Hamiltonian htemp into HsamX_k      

   inu = 0
   do i = 1,sample2%natom
     iatom = sample2%atom(i)
     in1 = imass(iatom)
     imu = 0
     do j = 1,sample2%natom
       jatom = sample2%atom(j)
       in2 = imass(jatom)
       do iorb = 1, num_orb(in1)
         ignu = iorb + degelec(iatom)
         jnu = inu + iorb
         do jorb = 1, num_orb(in2)
           igmu = jorb + degelec(jatom)
           jmu = imu + jorb
           Hsam2_k(jnu,jmu) = H_k(ignu,igmu,ikpt)
         enddo
       enddo
       imu = imu + num_orb(in2)
     enddo ! do j
     inu = inu + num_orb(in1)
   enddo ! do i


!! loop over atoms in tip2
!   i = 0
!   do itip = 1,sample2%natom_tip
!      
!      iatom = sample2%atom_tip(itip)
!      in1 = imass(iatom)
!      j = 0
!! loop over atoms in tip1
!      do jtip = 1,sample1%natom_tip
!         
!         jatom = sample1%atom_tip(jtip)
!         in2 = imass(jatom)
!         
!! vec = r_i - r_j; we consider only interactions in the basic unit cell 
!         vec(:) = ratom(:,jatom) - ratom(:,iatom)
!         dot = kpt(1)*vec(1) + kpt(2)*vec(2) + kpt(3)*vec(3)
!         phase = cmplx(cos(dot),sin(dot))
!
!! sample_2
!! store piece of hopping between iatom and jatom       
!         do imu = 1, num_orb(in1)
!            jmu = imu + degelec2(itip)
!            do inu = 1, num_orb(in2)
!               jnu = inu + degelec1(jtip)
!               Hsam2_k(jmu,jnu) = Hsam2_k(jmu,jnu) - phase*t_21_orig(i+imu,j+inu)
!               Hsam2_k(jnu,jmu) = Hsam2_k(jnu,jmu) - phase*t_12_orig(j+inu,i+imu)
!            end do ! do inu
!         end do ! do imu
!         j = j + num_orb(in2)
!      enddo ! enddo jtip
!      i = i + num_orb(in1)
!   enddo ! enddo itip

! Format Statements
! ===========================================================================
100     format (2x,2i5,f14.6)
200     format (9f12.5)
300     format (2i5,2f14.6)
400     format (4i5,f14.6)
450     format (5i5,2f14.6)
600     format ('ERR_SR:',2i5,f18.8)
610     format ('ERR_SI:',2i5,f18.8)
620     format ('ERR_HR:',2i5,f18.8)
630     format ('ERR_HI:',2i5,f18.8)

   return
   
 end subroutine assemble_Hsam
