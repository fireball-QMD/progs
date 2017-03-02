! assemble_t12.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_t12_fit ( )
   
   use dimensions
   use interactions
   use transport
   use configuration

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

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
   integer ind1 
   integer ind2
 
   real y

   logical :: bind1
   logical :: bind2

   real, dimension (3, 3)                    :: eps
   real, dimension (3)                       :: sighat
   real, dimension (3)                       :: r1
   real, dimension (3)                       :: r2
   real, dimension (3)                       :: r21
   real, dimension (3)                       :: vec
   real, dimension (numorb_max, numorb_max)  :: hopx

! Procedure
! ===========================================================================

   write (*,*) '  '

! allocate auxiliar variables

   write (*,*) ' Assemble T_12 with fitting ....'


! Now, create matrix T_12 coupling the samples H1 and H2.  
! Remeber, we want to treat only the hopping in basic unit cell.        

   if (ifithop .eq. 1) then
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
     r1(:) = ratom(:,iatom)
     j = 0
! loop over atoms in tip2
     do jtip = 1,sample2%natom_tip
         
      jatom = sample2%atom_tip(jtip)
      in2 = imass(jatom)
      r2(:) = ratom(:,jatom)

! --------------------------------------
!     FITTING T_12
! --------------------------------------

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
      r21(:) = r2(:) - r1(:)
      y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

! Find the unit vector in sigma direction.
      if (y .lt. 1.0d-05) then
       sighat(1) = 0.0d0
       sighat(2) = 0.0d0
       sighat(3) = 1.0d0
      else
       sighat(:) = r21(:)/y
      end if

      call epsilon (r2, sighat, eps) 

      ind1 = sample1%ispec(itip) 
      ind2 = sample2%ispec(jtip)

      call gethop (in1, in2, ind1, ind2, y, eps, hopx)

! store piece of hopping between iatom and jatom
      write (*,*) 'hop ',itip,jtip,in1,in2
      do imu = 1, num_orb(in1)
!       write (*,300) (hopx(imu,inu),inu=1,num_orb(in2))
       do inu = 1, num_orb(in2)
! hopping matrix
        t_12(i+imu,j+inu) = (1.0d0, 0.0d0)*hopx(imu,inu)
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

   endif ! if (ifithop .eq. 1)


 
! Format Statements
! ===========================================================================
300         format(9e18.6)
   return
   
 end subroutine assemble_t12_fit
