! assemble_Gsam.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_Gsam (ikpoint)

   use dimensions
   use interactions
   use transport
   use kpoints
   use neighbor_map
   use configuration
   use charges

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in)   ::  ikpoint

   complex,    dimension (:,:), allocatable    :: Gsam1_k
   complex,    dimension (:,:), allocatable    :: Gsam2_k

! Output


! Local Parameters and Data Declaration
! ===========================================================================

   real, parameter :: PI = 3.141592653589793238462643

! Local Variable Declaration and Description
! ===========================================================================

   integer ie
   integer i
   integer iatom
   integer jatom
   integer itip
   integer jtip
   integer in1
   integer in2
   integer inu
   integer inu0
   integer imu
   integer imu0
   integer jnu
   integer jmu
   integer ineigh
   integer mbeta

   integer norb1
   integer norb2

   real dot
   real dz
   real, dimension (3) :: vec
   real, dimension (3) :: k_temp

   complex a0
   complex a1
   complex phase
   complex omega
   complex omega0
   complex omega1
   complex omega10
   complex omega2
   complex omega20

   logical bind1
   logical bind2

! Lapack variables
   integer            :: info
   integer            :: lwork1
   integer            :: lwork2
   complex, dimension(:),  allocatable               :: work1
   complex, dimension(:),  allocatable               :: work2
   integer, dimension(:),  allocatable               :: ipiv1
   integer, dimension(:),  allocatable               :: ipiv2

! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) ' Assemble tip Green`s functions  '

   a1 = (1.0d0, 0.0d0)
   a0 = (0.0d0, 0.0d0)

! set the kpoint
   k_temp(:) = special_k(:,ikpoint)

! number of orbitals of samples
   norb1 = sample1%norb
   norb2 = sample2%norb

! aux Gr.f.
   allocate (Gsam1_k (norb1,norb1))
   allocate (Gsam2_k (norb2,norb2))

! Lapack variables
   lwork1 = norb1*norb1
   lwork2 = norb2*norb2
   allocate(work1(lwork1))
   allocate(work2(lwork2))
   allocate(ipiv1(norb1))
   allocate(ipiv2(norb2))

! shift Fermi level to zero
   do inu = 1,norb1
      Hsam1_k(inu,inu) = Hsam1_k(inu,inu) - efermi*a1
   enddo ! do inu
   do inu = 1,norb2
      Hsam2_k(inu,inu) = Hsam2_k(inu,inu) - efermi*a1
   enddo ! do inu

! First calculate auxilliar energy step needed to calculate correction to DOS
! comming from real part of Green's function
! usually we take:
!   n(E;k) = -1/Pi * Im{Gr[E+i*eta;k]}
! with real part correction  we have
!   n(E;k) = -1/Pi * Im{Gr[E+i*eta;k]} - 1/Pi*eta*d/dE( Re{Gr[E+i*eta;k]} )
! where Eo is previous step in integration

! set intial energy
   omega = (Elow + 0.5d0*real(dE))*(1.0d0,0.0d0) + eta*(0.0d0,1.0d0)
   omega0 = (Elow + 0.5d0*real(dE))*(1.0d0,0.0d0) + eta0*(0.0d0,1.0d0)

! Loop of energy
   do ie = 1,nE

     write (*,100) ie, real(omega)

! assign energy due to the bias voltage to the system 1
     omega1 = omega
     omega10 = omega0
! assign energy due to the bias voltage to the system 1
     if(real(omega) .lt. 0.0d0) then
      omega2 = omega - a1*Elow
      omega20 = omega0 - a1*Elow
     else
      omega2 = omega - a1*Eup
      omega20 = omega0 - a1*Eup
     endif

! ---------------------------------
!  SAMPLE_1
! ---------------------------------
! calc Green function for given energy (in k-space)
     Gsam1_k(:,:) = a0 - Hsam1_k(:,:)

! diagonal terms
     inu = 1
     do i = 1,sample1%natom
       iatom =  sample1%atom(i)
       in1 = imass(iatom)

! application of optional eta0
       if ( ideta(iatom) .eq. 1 ) then
         do imu = 1,num_orb(in1)
           Gsam1_k(inu,inu) = omega10 - Hsam1_k(inu,inu)
           inu = inu + 1
         enddo
       else
         do imu = 1,num_orb(in1)
            Gsam1_k(inu,inu) = omega1 - Hsam1_k(inu,inu)
            inu = inu + 1
          enddo
       endif
     enddo ! do i

! inversion
     write (*,*) 'Doing inversion for Gsam_1'
! old fashion
!!      call inv (Gsam1_k, Gsam1_k, norb1, norb1)
! LU matrix factorization of general matrix
     call zgetrf (norb1, norb1, Gsam1_k, norb1, ipiv1, info )
     if (info .ne. 0)  then
      write (*,*) ' ***  error in zgetrf  '
      write (*,*) ' ***  info = ',info
      stop
     endif
! matrix inversion of general matrix
     call zgetri (norb1, Gsam1_k, norb1, ipiv1, work1, lwork1, info)
     if (info .ne. 0) then
      write (*,*) ' ***  error in zgetri  '
      write (*,*) ' ***  info = ',info
      stop
     endif



! Restore tip_1 Green's function in real space
! loop over atoms of tip1
     do itip = 1,sample1%natom_tip
       iatom = sample1%atom_tip(itip)
       in1 = imass(iatom)

       do jtip = 1,sample1%natom_tip
         jatom = sample1%atom_tip(jtip)
         in2 = imass(jatom)

	     vec(:) = ratom(:,jatom) - ratom(:,iatom)
         dot = k_temp(1)*vec(1) + k_temp(2)*vec(2) + k_temp(3)*vec(3)
         phase = Conjg (cmplx(cos(dot),sin(dot)))
         do jnu = 1,num_orb(in1)
          inu = pointer1(itip) + jnu
          inu0 = degelec1(sample1%t2s(itip)) + jnu
          do jmu = 1,num_orb(in2)
            imu = pointer1(jtip) + jmu
            imu0 = degelec1(sample1%t2s(jtip)) + jmu
! G.f. tip1
            Gr_tip1(inu,imu,ie) = Gr_tip1(inu,imu,ie)                 &
                       + phase*Gsam1_k(inu0,imu0)*weight_k(ikpoint)
!@                       + Gsam1_k(inu0,imu0)*weight_k(ikpoint)
          enddo ! do jmu
!               if (ie .eq. 1) then
!                  dos1(inu,ie) = dos1(inu,ie) -                              &
!     &              (1.0d0/PI)*weight_k(ikpoint)*Imag(Gsam1_k(inu0,inu0)) +
!               else
!                  dos1(inu,ie) = dos1(inu,ie) -                              &
!     &              (1.0d0/PI)*weight_k(ikpoint)*Imag(Gsam1_k(inu0,inu0)) +      &
!     &               eta*(1.0d0/PI)*weight_k(ikpoint)*( real(Gsam1_k(inu0,inu0)) &
!     &               - real(Gaux1_k(inu0,inu0)) )/real(dE)
!               endif

         enddo ! do imu
       enddo ! do jtip
     enddo ! do itip

! ---------------------------------
!  SAMPLE_2
! ---------------------------------
! calc Green function for given energy (in k-space)
! energy omega
     Gsam2_k(:,:) = a0 - Hsam2_k(:,:)
     inu = 1
     do i = 1,sample2%natom
       iatom =  sample2%atom(i)
       in1 = imass(iatom)
! application of optional eta0
       if ( ideta(iatom) .eq. 1) then
          do imu = 1,num_orb(in1)
            Gsam2_k(inu,inu) = omega20 - Hsam2_k(inu,inu)
            inu = inu + 1
          enddo
       else
          do imu = 1,num_orb(in1)
            Gsam2_k(inu,inu) = omega2 - Hsam2_k(inu,inu)
            inu = inu + 1
          enddo
       endif
     enddo
!      Gsam2_k = omega*Idn2 - Hsam2_k

! inversion
     write (*,*) 'Doing inversion for Gsam_2'
! old fashion
!!      call inv (Gsam2_k, Gsam2_k, norb2, norb2)
! LU matrix factorization of general matrix
     call zgetrf (norb2, norb2, Gsam2_k, norb2, ipiv2, info )
     if (info .ne. 0)  then
       write (*,*) ' ***  error in zgetrf  '
       write (*,*) ' ***  info = ',info
       stop
     endif
! matrix inversion of general matrix
     call zgetri (norb2, Gsam2_k, norb2, ipiv2, work2, lwork2, info)
     if (info .ne. 0) then
      write (*,*) ' ***  error in zgetri  '
      write (*,*) ' ***  info = ',info
      stop
     endif

! Restore tip_1 Green's function in real space
! loop over atoms of tip2
     do itip = 1,sample2%natom_tip
       iatom = sample2%atom_tip(itip)
       in1 = imass(iatom)

       do jtip = 1,sample2%natom_tip
         jatom = sample2%atom_tip(jtip)
         in2 = imass(jatom)

!@            vec(:) = 0.0d0
         vec(:) = ratom(:,jatom) - ratom(:,iatom)
         dot = k_temp(1)*vec(1) + k_temp(2)*vec(2) + k_temp(3)*vec(3)
         phase = Conjg (cmplx(cos(dot),sin(dot)))
         do jnu = 1,num_orb(in1)
           inu = pointer2(itip) + jnu
           inu0 = degelec2(sample2%t2s(itip)) + jnu
           do jmu = 1,num_orb(in2)
             imu = pointer2(jtip) + jmu
             imu0 = degelec2(sample2%t2s(jtip)) + jmu
! G.f. tip2
             Gr_tip2(inu,imu,ie) = Gr_tip2(inu,imu,ie)          &
                  + phase*Gsam2_k(inu0,imu0)*weight_k(ikpoint)
!@                       + Gsam2_k(inu0,imu0)*weight_k(ikpoint)
           enddo ! do jmu
         enddo ! do imu
! sum DOS2
!            if (ie .eq. 1) then
!               dos2(inu,ie) = dos2(inu,ie) -                              &
!     &              (1.0d0/PI)*weight_k(ikpoint)*Imag(Gsam2_k(inu0,inu0)) +
!            else
!               dos2(inu,ie) = dos2(inu,ie) -                              &
!     &              (1.0d0/PI)*weight_k(ikpoint)*Imag(Gsam2_k(inu0,inu0)) +      &
!     &               eta*(1.0d0/PI)*weight_k(ikpoint)*( real(Gsam2_k(inu0,inu0)) &
!     &               - real(Gaux2_k(inu0,inu0)) )/real(dE)
!            endif

       enddo ! do jtip
     enddo ! do itip

! increment the energy
     omega = omega + dE
     omega0 = omega0 + dE
   enddo ! do iE

   deallocate (Gsam1_k)
   deallocate (Gsam2_k)
   deallocate (work1)
   deallocate (work2)
   deallocate (ipiv1)
   deallocate (ipiv2)

! Format Statements
! ===========================================================================
100     format (2x,'Step = ', i6,'  Energy = ',f16.8)
300     format (2i5,2f14.6)
   return
 end subroutine assemble_Gsam
