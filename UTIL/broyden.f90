
subroutine broyden ( x_try, x_old, beta, r2, iter, max_order, nmsh,      &
     &                   max_scf_iterations)

  use charges
  implicit none

! Argument Declaration and Description
! ===========================================================================
! input
  integer, intent(in) :: nmsh      ! Size of vectors being optimized
  real, intent(in) :: beta         ! Mixing factor - used just for starting guess and simple mixing
  integer, intent(in) :: iter      ! iteration number
  integer, intent(in) :: max_order ! How far back do we go to extrapolate? - not used
  integer, intent(in) :: max_scf_iterations ! Max SCF iterations
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
  real, allocatable, dimension(:) :: deltaF
  real, allocatable, dimension(:) :: deltaX

! auxiliary matrices and vectors
  real, allocatable, dimension(:)   :: auxvec
  real, allocatable, dimension(:,:) :: amat

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

! Program Description
! ===================================================================
!   *                                                              *
!   *   Modified Broyden's method for convergence acceleration     *
!   *      in solving systems of non linear equations              *
!   *                                                              *
!       ( D.Vanderbilt & S.G.Louie, Phys.Rev.B 30 (1984) 6118 )


!  Adaptation of Fortran77 procedure hipl.f.
!  This file implements the original Broyden's method (as described in the article.
!  For the Vanderbilt and Louie's modification (as described in the article), see the procedure louie.f90
!
! ===========================================================================  

  write(*,*)'  '
  write(*,*)' Welcome to broyden; mix charges to self-consistency.'


! Allocation of arrays
! The important one for Broyden is RJac.
  if(.not. allocated(Fv))then
     ! Why 2? We need the current and the one before.
     ! I guess we could do with just one, but it would make the code less legible.
     allocate (Fv(nmsh,2))
     allocate (Xv(nmsh,2))
     allocate (r2_sav(max_scf_iterations))
     allocate (RJac(nmsh,nmsh))
  end if
  
     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in Broyden '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        stop
     end if


     Xv(:,1) = Xv(:,2)
     Fv(:,1) = Fv(:,2)
     Xv(:,2) = x_old(:)
     Fv(:,2) = x_try(:) - x_old(:)
     r2 = dot_product(Fv(:,iter),Fv(:,iter))
     r2 = r2 / nmsh




! What order interpolation do we use this time - this is used just for recognizing the very first iteration
     mix_order = min(iter,max_order)

! Converged
     if(r2 .lt. tr2) then
        x_old(:) = x_old(:) + beta*Fv(:,2)
        return
     end if

! First iteration - guess the starting Jacobian.
     if((mix_order .eq. 1)) then
        x_old(:) = x_old(:) + beta*Fv(:,2)

        RJac(:,:) = 0.0d0
        aux = 1/beta
        do i = 1, nmsh
           RJac(i,i) = aux
        end do

        return
     end if ! mix_order .eq. 1

! Allocate, what is needed.
     allocate (deltaF(nmsh))
     allocate (deltaX(nmsh))
     allocate (auxvec(nmsh))
     allocate (amat(nmsh,nmsh))
     allocate (ipiv(nmsh))


     deltaX = Xv(:,2) - Xv(:,1)
     renorm = sqrt(dot_product( deltaX, deltaX ))

     deltaF = (Fv(:,2) - Fv(:,1)) / renorm !Eq A8
     deltaX = deltaX / renorm                      !Eq A7

! This could be probably optimized further.
     auxvec(:) = 0.0d0
     call dgemv('n',nmsh,nmsh,1.0d0,RJac,nmsh,deltaX,1,0.0d0,auxvec,1)

     do i = 1,nmsh
        do j = 1, nmsh
           RJac(i,j) = RJac(i,j) - ( deltaF(i) + auxvec(i) ) * deltaX(j) !Eq A10
        end do
     end do
     
! Get the correction - Eq. A6
     amat(:,:) = RJac(:,:)
     auxvec(:) = Fv(:,2)
     ipiv(:) = 0
     info = 0
     nrhs = 1
     lda = nmsh
     ldb = nmsh

     call dgesv(nmsh,nrhs,amat,lda,ipiv,auxvec,ldb,info) !get x^(m+1) - x^(m) from Eq A6

     x_old(:) = x_old(:) + auxvec(:)


! deallocate all
     deallocate(deltaF,deltaX,auxvec)
     deallocate(amat)
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
!!$          r = r + vec(i)**2
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
!!$     end function  maxnorm



   end subroutine broyden

