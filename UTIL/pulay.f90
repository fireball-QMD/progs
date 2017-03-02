subroutine pulay ( x_try, x_old, beta, r2, iter, max_order, nmsh,      &
     &                   max_scf_iterations)

  use charges
!  use density, only: mwe, drwe
  implicit none



! pulay.f90
! Program Description
! ===========================================================================
! G.Kresse, J. Furthmuller: Efficient iterative schemes for ab initio
! total-energy calculations using a plane wave basis set.
! Physical Review B, vol.54, no 16, 1996
! paragraph IV A,B
!
! See also:
! Pulay: Convergence acceleration of iterative sequences. The case of SCF iteration
! Chemical Physics Letters, vol.73, no 2, 1980
!
!
!
!
!
!
!
!
! ===========================================================================



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

! initial mixing parameters - works fine according to the article
  real, parameter :: A0 = 0.2d0  ! 0.2 can sometimes work better
  real, parameter :: q0 = 1.5d0  ! in angstrom^-1


! Local Variables
! ===========================================================================
  real, allocatable, dimension(:) :: deltaF
  real, allocatable, dimension(:) :: deltaX

  real, allocatable, dimension(:) :: metric   !for Eq.64
  real, allocatable, dimension(:) :: mixcoeff !alphas in Eq.52

! auxiliary matrices and vectors
  real, allocatable, dimension(:)     :: auxvec
  real, allocatable, dimension(:)     :: auxvec2
  real, allocatable, dimension(:,:) :: amat
  real, allocatable, dimension(:,:) :: bmat
  real, allocatable, dimension(:,:) :: cmat
  real, allocatable, dimension(:,:) :: idmat

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

!=======================================================================================

  write(*,*)'  '
  write(*,*)' Welcome to pulay; mix charges to self-consistency.'





! Allocation of arrays
  if(.not. allocated(Fv))then
     allocate (Fv(nmsh,max_scf_iterations))
     allocate (Xv(nmsh,max_scf_iterations))
     allocate (r2_sav(max_scf_iterations))
!     allocate (pulmetric(nmsh))

  end if



     if(max_order .le. 0 .or. iter .gt. max_scf_iterations) then
        write (*,*) ' Stop in pulay '
        write (*,*) ' max_order=', max_order
        write (*,*) ' iter=',iter
        write (*,*) ' max_scf_iterations=',max_scf_iterations
        stop
     end if

! save previous step
     Xv(:,iter) = x_old(:)
 !residual vector
!     Fv(:,iter) = (x_try(:) - x_old(:))
     Fv(:,iter) = (x_try(:) - x_old(:))*mwe(:)        !residual vector
     r2 = dot_product(Fv(:,iter),Fv(:,iter))
     r2 = r2 / nmsh
! dump filtering
     do i = 1,nmsh
!       if (mwe(i)**2 .gt. 0.001d0) write (1000+iter,*) Fv(i,iter), mwe(i)
       write (1000+iter,*) Fv(i,iter), drwe(i)
       write (2000+iter,*) Fv(i,iter), mwe(i)
     enddo

! What order interpolation do we use this time
! Not used. All iterations are counted.
! The article seems to assume this - is it OK?
     mix_order = min(iter,max_order)

! Converged
     if(r2 .lt. tr2) then
        x_old(:) = x_old(:) + beta*Fv(:,iter)/mwe(:)
        return
     end if

! First 4 iterations - use initial approximation
     if(mix_order .lt. 4) then

!        allocate (auxvec(nmsh))
!        do i = 1, nmsh
!           auxvec(i) = A0
!           auxvec(i) = A0  x_old(i) / (x_old(i) * q0) !Eq.61
!        end do

!        x_old(:) = x_old(:) + auxvec(:) * Fv(:,iter)
!         x_old(:) = x_old(:) + bmix * Fv(:,mix_order)/mwe(:)
         write (*,*) ' Doing simple mixing', beta
         x_old(:) = (1.0d0-beta)*x_old(:) + beta*x_try(:)
!        deallocate (auxvec)
        return
     end if


     allocate (mixcoeff(mix_order + 1)) !the last coefficient stands for the lagrangian multiplier for constraint of Eq.55
     allocate (amat(mix_order + 1, mix_order + 1))
     allocate(ipiv(mix_order + 1))

! Determine a new metric
! "The choice of q1 is relatively unimportant and we set q1 in a way,
!  that the shortest wave vector is weighted 20 times stronger than the longest wave vector."
!     do i=1,nmsh
!        pulmetric(i) = x_old(i) + q1 !Eq.64
!        pulmetric(i) = 1.0d0 ! provizorni
!     end do

! Build the linear equations to be solved
     amat = 0.0d0
     do i=1,mix_order
        do j=1,mix_order
           do k=1,nmsh
              amat(i,j) = amat(i,j) +  Fv(k,iter - mix_order + i) * Fv(k,iter - mix_order + j)
           end do
        end do
     end do
     amat(mix_order + 1,:) = -1.0d0
     amat(:,mix_order + 1) = -1.0d0
     amat(mix_order + 1, mix_order + 1) = 0.0d0
     mixcoeff(:) = 0.0d0
     mixcoeff(mix_order + 1) = -1.0d0


! Solve linear equations using LAPACK
     ipiv = 0
     nrhs = 1
     lda = mix_order + 1
     ldb = mix_order + 1
     info = 0
!     write(*,*) 'Calling dgesv'
     call dgesv(mix_order + 1,nrhs,amat,lda,ipiv,mixcoeff,ldb,info)
!     write(*,*) info
!     write(*,*) sqrt(mixcoeff(mix_order+1))
     write (*,*) 'sum_alpha =', sum(mixcoeff(1:mix_order))

! Get new input vector
! Adding a little bit of Fv as not to stick in the space spanned by the first iteration vectors
     x_old(:) = 0.0d0
     do i = 1, mix_order
!        x_old(:) = x_old(:) + mixcoeff(i) * Xv(:,i) !Eq.52
        x_old(:) = x_old(:) + mixcoeff(i) * Xv(:,i) !Eq.52
!        x_old(:) = x_old(:) + mixcoeff(i) * (Xv(:,iter - mix_order + i) + beta * Fv(:,iter - mix_order + i))
     end do

! deallocate all
     deallocate(mixcoeff)
     deallocate(amat,ipiv)

!     write(*,*) 'ending'

     return


! Format Statements
! ===========================================================================
300  format (2x, ' norm =  ', f20.4)
305  format (2x, ' norm of difference between auxvecs =  ', f20.4)
301  format (2x, ' norm of the old vector =  ', f20.4)
302  format (2x, ' norm of the new guess =  ', f20.4)
308  format (2x, ' norm of difference between betaInvH, betaH^-1 =  ', f20.4)
242  format (2x, ' Selected method is:  ', a17)



   contains

!!$     real function dot_metr(x,y) ! scalar product with metric
!!$       real, intent(in), dimension(nmsh) :: x
!!$       real, intent(in), dimension(nmsh) :: y
!!$       real sum
!!$       sum = 0.0d0
!!$       do i=1,nmsh
!!$          sum = sum + metric(i) * x(i) * y(i)
!!$       end do
!!$       dot_metr = sum
!!$     end function dot_metr

! debugging...

     real function vecnorm(vec)
       real, intent(in), dimension(nmsh) :: vec
       integer i,j
       real r
       r = 0.0d0
       do i = 1,nmsh
          r = r + vec(i)**2
       end do
       r = sqrt(r)
       vecnorm = r
     end function  vecnorm

     real function maxnorm(mat)
       real, intent(in), dimension(nmsh,nmsh) :: mat
       integer i,j
       real m
       real r
       r = 0.0d0
       m = 0.0d0
       do i = 1,nmsh
          do j = 1,nmsh
             r = r + abs(mat(i,j))
          end do
          if (r .gt. m)  m = r
          r = 0.0d0
       end do
       maxnorm = m
     end function  maxnorm



   end subroutine pulay
