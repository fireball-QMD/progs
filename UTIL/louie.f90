


subroutine louie ( x_try, x_old, beta, r2, iter, max_order, nmsh,      &
     &                   max_scf_iterations)

  use charges
  implicit none

! Argument Declaration and Description
! ===========================================================================
! input
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real, intent(in) :: beta         ! Mixing factor - not used
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate? - not used
  integer, intent(in) :: max_scf_iterations ! Max SCF iterations - not really used
  real, intent(in), dimension(nmsh) :: x_try ! potential new vector on input

! input and output
  real, intent(inout), dimension(nmsh) :: x_old ! old vector in input, real new vector on output

! output
  real, intent(out) :: r2 ! mean-square of (x_try(i)-x_old(i))**2

! Local Parameters and Data Declaration
! ===========================================================================
  real, parameter :: tr2=2.0e-15  ! convergence factor, if r2<tr2, assume converged


! Local Variables
! ===========================================================================  

! weigths for Vanderbilt-Louie method
! wgt0 is the weight of the initial guess
! wgtl is for subsequent vectors
! It is also possible to have wgtl as a vector, so that different iterations have different weights.
  real wgt0
  real wgrad
  real wgt1

  real, allocatable, dimension(:) :: deltaF
  real, allocatable, dimension(:) :: deltaX

! auxiliary matrices and vectors
  real, allocatable, dimension(:)   :: auxvec
  real, allocatable, dimension(:,:) :: amat
  real, allocatable, dimension(:,:) :: bmat

  real renorm
  real aux
  real aux2
  real norm
  integer i,j,k
  integer mix_order ! Actual order used min(iter,max_order)

! For lapack
  integer, allocatable, dimension(:) :: ipiv
  integer info
  integer nrhs
  integer lda
  integer ldb
  integer ldx


! ===================================================================
!   *                                                              *
!   *   Modified Broyden's method for convergence acceleration     *
!   *      in solving systems of non linear equations              *
!   *                                                              *
!       ( D.Vanderbilt & S.G.Louie, Phys.Rev.B 30 (1984) 6118 )


!  Adaptation of Fortran77 procedure hipl.f.
!  This file implements the Vanderbilt and Louie's modification of the original Broyden's method.
!  For the original Broyden's method (as described in the article), see the procedure broyden.f90
!
! ===========================================================================


  write(*,*)'  '
  write(*,*)' Welcome to louie; mix charges to self-consistency.'


! These should somehow depend on beta.
! However, I haven't been able to discover, how to set the weights reasonably,
! so the current settings is taken literally from hipl.f
! Note that the article suggests weights, however those seem unsuitable for charge mixing.
  wgt0 = 0.05d0
  wgrad = 0.01d0
  wgt1 = 1.0d0


! allocate static arrays
  if(.not. allocated(Fv))then
     ! Why 2? We need the current and the one before.
     ! I guess we could do with just one, but it would make the code less legible.
     allocate (Fv(nmsh,2))
     allocate (Xv(nmsh,2))

     allocate (r2_sav(max_scf_iterations))

     allocate (betaInvH(nmsh,nmsh))
     allocate (gamaH(nmsh,nmsh))


!!$!for testing if the inverse of beta goes well
!!$      allocate (betaH(nmsh,nmsh))       
!!$      betaH(:,:) = 0.0d0
!!$   end if
  end if
  
     
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in Louie '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        stop
     end if


     Xv(:,1) = Xv(:,2)
     Fv(:,1) = Fv(:,2)

     Xv(:,2) = x_old(:)
     Fv(:,2) = x_try(:) - x_old(:)

     r2 = dot_product(Fv(:,2),Fv(:,2))
     r2 = r2 / nmsh




! What order interpolation do we use this time - this is used just for recognizing the first iteration.
     mix_order = min(iter,max_order)

! Converged
     if(r2 .lt. tr2) then
        x_old(:) = x_old(:) + beta*Fv(:,2)
        return
     end if

! First step - has to initialize betaH and gamaH
     if((mix_order .eq. 1)) then ! .or. (itipo .eq. 2 .and. mix_order .le. 3)

        x_old(:) = x_old(:) + beta*Fv(:,2)

        betaInvH(:,:) = 0.0d0
        gamaH(:,:) = 0.0d0
        aux = 1.0d0 / wgrad**2
        aux2 = wgrad**2 / wgt0
        
        do i = 1, nmsh
           betaInvH(i,i) = aux
           gamaH(i,i) = aux2
        end do
        
        return
     end if !mix_order .eq. 1

! Allocate, what is needed.
! Ideally, just one auxiliary matrix (amat) should be used.
     allocate (deltaF(nmsh))
     allocate (deltaX(nmsh))
     allocate (auxvec(nmsh))
     allocate (amat(nmsh,nmsh))
     allocate (bmat(nmsh,nmsh))

     allocate (ipiv(nmsh))

     deltaX = Xv(:,2) - Xv(:,1)
     renorm = sqrt(dot_product( deltaX, deltaX ))

     deltaF = (Fv(:,2) - Fv(:,1)) / renorm !Eq A8
     deltaX = deltaX / renorm                      !Eq A7

! comment to realization of Eq A16:
! u = v = wgt1 * deltaX
! auxvec = betaInvH . deltaX
! aux = sigma
        aux = 1.0d0
        aux2 = wgt1**2

        do i = 1,nmsh
           do j = 1, nmsh
              gamaH(i,j) = gamaH(i,j) -  aux2 * deltaF(i) * deltaX(j)  !Eq A15
              amat(i,j) = aux2 * deltaX(i) * deltaX(j)                 !Eq A16 first part
              aux = aux + aux2 * deltaX(i) * betaInvH(i,j) * deltaX(j) !Eq A16 second part
           end do
        end do

        bmat(:,:) = 0.0d0

        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,amat,nmsh,betaInvH,nmsh,0.0d0,bmat,nmsh)
        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,betaInvH,nmsh,bmat,nmsh,0.0d0,amat,nmsh)

        betaInvH(:,:) = betaInvH(:,:) - (1/aux) * amat(:,:)  !Eq A16 finish

! get jacobian via updated betaH^-1
!        write(*,*) 'get Jacobian'
        call dgemm('n','n',nmsh,nmsh,nmsh,1.0d0,gamaH,nmsh,betaInvH,nmsh,0.0d0,amat,nmsh) !Eq A13, efficiently

! Invert the matrix to get the correction vector - Eq. A6
     auxvec(:) = Fv(:,2)
     info = 0
     nrhs = 1
     lda = nmsh
     ldb = nmsh

     call dgesv(nmsh,nrhs,amat,lda,ipiv,auxvec,ldb,info) !get x**(m+1) - x**m


     
     x_old(:) = x_old(:) + auxvec(:)

! deallocate all
     deallocate(deltaF,deltaX,auxvec)
     deallocate(amat,bmat)
     deallocate(ipiv)


     return


! Format Statements
! ===========================================================================
300  format (2x, ' norm =  ', f20.4)
305  format (2x, ' norm of difference between auxvecs =  ', f20.4)
301  format (2x, ' norm of the old vector =  ', f20.4)
302  format (2x, ' norm of the new guess =  ', f20.4)
308  format (2x, ' norm of difference between betaInvH, betaH^-1 =  ', f20.4)
242  format (2x, ' Selected method is:  ', a17)



!!$   contains
!!$   
!!$     real function vecnorm(vec)
!!$       real, intent(in), dimension(nmsh) :: vec
!!$       integer i,j
!!$       real r
!!$       r = 0.0d0
!!$       do i = 1,nmsh
!!$             r = r + vec(i)**2
!!$       end do
!!$       r = sqrt(r)
!!$       vecnorm = r     
!!$     end function  vecnorm
!!$
!!$     real function maxnorm(mat)
!!$       real, intent(in), dimension(nmsh,nmsh) :: mat
!!$       integer i,j
!!$       real m
!!$       real r
!!$       r = 0.0d0
!!$       m = 0.0d0
!!$       do i = 1,nmsh
!!$          do j = 1,nmsh
!!$             r = r + abs(mat(i,j))
!!$          end do
!!$          if (r .gt. m)  m = r
!!$          r = 0.0d0
!!$       end do
!!$       maxnorm = m
!!$   end function  maxnorm



   end subroutine louie
 
