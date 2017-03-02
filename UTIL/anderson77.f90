! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! anderson.f90
! Program Description
! ===========================================================================
!       Donald G. Anderson, J. Assoc. Computing Machinery, 12, 547 (1965)
!       Generalized by:  V. Eyert, J. Comp. Phys. 124, 271 (1996)
!       Use Eq. 7.7, 5.31 of Eyert
!       We set all betas to be equal for simplicity
!       Computes input vector for the next iteration
!       Order is equal to M+1 in Eyert.
!       This is the same as the best Broyden method (see Eyert)
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
      subroutine anderson(x_try,x_old,beta,r2,iter,max_order,nmsh,max_scf_iterations)

      implicit none

! Argument Declaration and Description
! ===========================================================================
! input
      integer, intent(in) :: nmsh      ! Size of vectors being optimized
      real, intent(in) :: beta         ! Mixing factor
      integer, intent(in) :: iter      ! iteration number
      integer, intent(in) :: max_order ! How far back do we go to extrapolate?
      integer, intent(in) :: max_scf_iterations ! Max SCF iterations
      real, intent(in), dimension(nmsh) :: x_try ! potential new vector on input

! input and output
      real, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real new vector on output

! output
      real, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2

! Local Parameters and Data Declaration
! ===========================================================================
      real, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged

! Local Variable Declaration and Description
! ===========================================================================
      real, save, allocatable, dimension(:,:) :: F   ! x_try-x_old saved
      real, save, allocatable, dimension(:,:) :: X   ! x_old saved
      real, save, allocatable, dimension(:,:) :: delF! F(m+1)-F(m) saved
      real, save, allocatable, dimension(:,:) :: delX! X(m+1)-X(m) saved
      real, save, allocatable, dimension(:)   :: r2_sav
      real, allocatable, dimension(:,:) :: a_matrix  ! Eq. 5.17
      real, allocatable, dimension(:) :: delF_F      ! <delF|F> in Eq. 5.31
      real, allocatable, dimension(:) :: contribution

      integer iloop
      integer jloop
      integer mix_order ! Actual order used min(iter,max_order)

! For ssysv call
      integer lwork
      real, allocatable, dimension(:) :: work
      integer, allocatable, dimension(:) :: ipiv
      integer info
 
! Procedure
! ===========================================================================
      if(.not. allocated(F))then
         allocate(F(nmsh,max_scf_iterations))
         allocate(X(nmsh,max_scf_iterations))
         allocate(delF(nmsh,max_scf_iterations))
         allocate(delX(nmsh,max_scf_iterations))
         allocate(r2_sav(max_scf_iterations))
      end if
      if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
         write(*,*) ' Stop in Anderson '
         write(*,*) ' max_order=', max_order
         write(*,*) ' iter=',iter
         write(*,*) ' max_scf_iterations=',max_scf_iterations
         stop
      end if
 
      X(:,iter) = x_old(:)
      F(:,iter) = x_try(:) - x_old(:)

      r2 = dot_product(F(:,iter),F(:,iter))
      r2 = r2/nmsh

! What order interpolation do we use this time
      mix_order = min(iter,max_order)

! Only doing simple extrapolation, or converged
      if(mix_order .eq. 1 .or. r2 .lt. tr2) then
         x_old(:)=x_old(:)+beta*F(:,iter)
         return
      end if

! Add to r2_sav
      r2_sav(iter) = r2

! Add to delF
      delF(:,iter-1) = F(:,iter)-F(:,iter-1) ! Eq. 5.6

! Add to delX
      delX(:,iter-1) = X(:,iter)-X(:,iter-1) ! Eq. 5.5

! Make sure step with lowest r2 value is always used (ie. not lost)
      if (iter .gt. max_order .and. max_order .ge. 6) then
        if (r2_sav(iter-max_order) .lt. minval(r2_sav(iter-max_order+1:iter))) then
!         Throw away second oldest step instead
          r2_sav(iter-max_order+1) = r2_sav(iter-max_order)
          F(:,iter-max_order+1) = F(:,iter-max_order)
          X(:,iter-max_order+1) = X(:,iter-max_order)
          delX(:,iter-max_order+1) = X(:,iter-max_order+2)-X(:,iter-max_order+1)
          delF(:,iter-max_order+1) = F(:,iter-max_order+2)-F(:,iter-max_order+1)
        end if
      end if

! Build a_matrix Eq. 5.17
 888  allocate(a_matrix(iter-mix_order+1:iter-1,iter-mix_order+1:iter-1))
      do iloop = iter-mix_order+1, iter-1
        do jloop = iter-mix_order+1, iter-1
          a_matrix(iloop,jloop) = dot_product(delF(:,iloop),delF(:,jloop))
        end do
      end do

! Build delF_F Eq. 5.31
      allocate(delF_F(iter-mix_order+1:iter-1))
      do iloop = iter-mix_order+1, iter-1
        delF_F(iloop) = dot_product(delF(:,iloop),F(:,iter))  
      end do

! Solve for gammas Eq. 5.31, 7.4 (We move a-inverse to other side: a * gamma = <|>)
      lwork = (mix_order-1)**2
      allocate (work(lwork))
      allocate (ipiv(mix_order-1))
      info=0
      call ssysv('U', mix_order-1, 1, a_matrix, mix_order-1, ipiv, delF_F, &
     &                  mix_order-1, work, lwork, info )
      if(info .ne. 0) then
         write(*,*) ' Error in Anderson, info =',info
         if(mix_order .le. 2)stop ! if you can't solve a 2x2 something is wrong
         mix_order=mix_order-1
         deallocate(work,ipiv,delF_F,a_matrix)
         goto 888 ! Try again with lower order
      end if

! Generate new guess at charges Eq. 7.7  (delF_F is now gamma)
      x_old(:)=x_old(:)+beta*F(:,iter)  ! First-order term
      do iloop = iter-mix_order+1, iter-1
        x_old(:)=x_old(:) - delF_F(iloop)*(delX(:,iloop) + beta*delF(:,iloop))
      end do

!     This code explicitly calculates the contribution of each
!     SCF iteration.  Simply a reworking of Eq. 7.7.  Useful
!     for debugging.
!     allocate(contribution(iter-mix_order+1:iter))
!     contribution(iter-mix_order+1) = (delF_F(iter-mix_order+1) - 0)
!     do iloop = iter-mix_order+2, iter-1
!       contribution(iloop) = delF_F(iloop) - delF_F(iloop-1)
!     end do
!     contribution(iter) = 1 - delF_F(iter-1)
!     Now use contributions (same answer as above).
!     x_old(:) = 0
!     do iloop = iter-mix_order+1, iter
!       x_old(:)=x_old(:) + contribution(iloop)*(X(:,iloop) + beta*F(:,iloop))
!     end do
!     deallocate(contribution)

      deallocate(delF_F,a_matrix,work,ipiv)

! Format Statements
! ===========================================================================
      return
      end
