! assemble_t12.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_t12_bare (kpt, ikpt, weightk)
   
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
   real, intent (in)                 :: weightk

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
   integer mbeta
 
   real dot
   real y
   real, dimension (3)                       :: vec

   complex a0
   complex a1
   complex phase

! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) ' Assemble T12 ....'

! set aux variables
   a0 = cmplx(0.0d0,0.0d0)
   a1 = cmplx(1.0d0,0.0d0)

! Now, create matrix T_12 coupling the samples H1 and H2.  
! Remeber, we want to treat only the hopping in basic unit cell.        
!   icelln = 0
!   do icellx = -sample1%ncell, sample1%ncell
!    do icelly = -sample1%ncell, sample1%ncell
!     norbtipx = icelln * system1%norb_tip
!     icelln = icelln + 1
!     r_icell(:) = icellx*a1vec(:) + icelly*a2vec(:)
!     jcelln = 0
!     do jcellx = -sample2%ncell, sample2%ncell
!      do jcelly = -sample2%ncell, sample2%ncell
!       norbtipy = jcelln * sample2%norb_tip
!       jcelln = jcelln + 1
!       r_jcell(:) = jcellx*a1vec(:) + jcelly*a2vec(:)

! loop over atoms in tip1
   i = 0
   do itip = 1,sample1%natom_tip
      
      iatom = sample1%atom_tip(itip)
      in1 = imass(iatom)
      j = 0
! loop over atoms in tip2
      do jtip = 1,sample2%natom_tip
         
         jatom = sample2%atom_tip(jtip)
         in2 = imass(jatom)

! --------------------------------------
!     NO FITTING T_12
! --------------------------------------
! vec = r_i - r_j; we consider only interactions in the basic unit cell 
! FIREABALL kspace2.f90
!         vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
!         dot = sks(1)*vec(1) + sks(2)*vec(2) + sks(3)*vec(3)
!         phase = cmplx(cos(dot),sin(dot))
!         vec(:) = 0.0d0
         vec(:) = ratom(:,jatom) - ratom(:,iatom)
         dot = kpt(1)*vec(1) + kpt(2)*vec(2) + kpt(3)*vec(3)
         phase = Conjg(cmplx(cos(dot),sin(dot)))

! store piece of hopping between iatom and jatom       
         do imu = 1, num_orb(in1)
          jmu = imu + degelec(iatom)
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
! hopping matrix
           t_12(i+imu,j+inu) = t_12(i+imu,j+inu)                    &
                               + real(phase*H_k(jmu,jnu,ikpt)*weightk)
          end do ! do inu
         end do ! do imu

! increment orbital pointer tip2
         j = j + num_orb(in2)
      enddo ! enddo jtip
! increment orbital pointer tip1
      i = i + num_orb(in1)
   enddo ! enddo itip

! get t_21
   t_21 = Conjg (Transpose (t_12))

!   write (222,*) ' ikpt =',ikpt,weightk,kpt(:)
!   do i = 1, sample1%norb_tip
!    write (222,300) (t_12_orig(i,j), j = 1, sample2%norb_tip)
!   end do
!   write (222,*) ''
!   do i = 1, sample1%norb_tip
!    write (222,300) (real(t_12(i,j)), j = 1, sample2%norb_tip)
!   end do
!   write (222,*) ''
!   do i = 1, sample1%norb_tip
!    write (222,300) (imag(t_12(i,j)), j = 1, sample2%norb_tip)
!   end do
!   write (222,*) ''
!   do i = 1, sample2%norb_tip
!    write (222,300) (real(t_21(i,j)), j = 1, sample1%norb_tip)
!   end do

! Format Statements
! ===========================================================================
300         format(9f14.5)
   return
   
 end subroutine assemble_t12_bare
