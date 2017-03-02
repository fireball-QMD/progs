! sqrt_mat.f90
! Program Description
! ===========================================================================
! A .. OUTPUT MATRIX [A^(1/2)]
! B .. INPUT MATRIX
! Formula:
!   A^(1/2) = (z^(-1) * D * z)^(1/2)  where D..diagonal matrix
!  if A is symmetric A = A_t (transpose matrix) then 
!  (z^(-1) * D * z)^(1/2) = z^(-1) * D^(1/2) * z
!   z ..matrix of eigenvectors
!   D .. matrix of eignevalues (off-diag = 0., trace = eigenvalues)
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================

subroutine c_sqrt_m (b, a, ndim, iprn)

  use matmult

! Argument Declaration and Description
! ===========================================================================
! Input

  integer, intent (in)                :: ndim
  logical, intent (in)                :: iprn
  complex, dimension (ndim,ndim), intent(in)   :: b
  complex, dimension (ndim,ndim), intent(out)  :: a

! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

  integer imu

  complex sqlami
  complex a0
  complex a1

  complex, allocatable, dimension (:,:) :: work
  real, allocatable, dimension (:) :: rwork
  real, allocatable, dimension (:) :: eigv
  integer lwork
  integer info

! Procedure
! ===========================================================================

! set aux variables
  a0 = cmplx(0.0d0,0.0d0)
  a1 = cmplx(1.0d0,0.0d0)

! allocate 
  lwork = ndim*ndim
  allocate (rwork(lwork))
  allocate (eigv(ndim))
  allocate (work(ndim, ndim))

! copy to output matrix
  a(:,:) = b(:,:)

  if ( iprn ) write (*,*) 'Calling zheev subroutine ...'
  call zheev ('V', 'U', ndim, a, ndim, eigv, work, lwork, rwork, info)

  if ( info .gt. 0) then
     write (*,*) ' Error occured in eigenvalue solver  !!!'
     write (*,*) ' Program will be aborted                '
     stop
  else
      
     if ( iprn ) then
        write (*,*) '                      eigenvalues:                    '
        write (*,*) '======================================================'
        write (*,200) (eigv(imu), imu=1,ndim)
!        write (*,*) '--------------------'
!        write (*,*) 'EigenVectors:'
!        do imu = 1, ndim
!           write(*,100) (a(inu,imu), inu=1,ndim)
!        enddo
     endif ! if(iprn)

     do imu = 1, ndim
        sqlami = eigv(imu)**(-0.25d0)
        a(:,imu) = a(:,imu)*sqlami
     end do

     if ( iprn ) write (*,*) 'Calling zgemm subroutine ...'
     call zgemm ('N', 'C', ndim, ndim, ndim, a1, a, ndim, a, ndim, a0,   &
     &           work, ndim)

! copy A^(-1/2)
     a(:,:) = work(:,:)

  endif ! if (info .gt. 0)
      

! Format Statements
! ===========================================================================
!100   format ( <ndim>f8.3 )
! 100   format ( 50f8.3 )
200   format (4(2x, f12.4))

  return

end subroutine c_sqrt_m








subroutine c_sqrt (b, a, ndim, iprn)

  use matmult

! Argument Declaration and Description
! ===========================================================================
! Input

  integer, intent (in)                :: ndim
  logical, intent (in)                :: iprn
  complex, dimension (ndim,ndim), intent(in)   :: b
  complex, dimension (ndim,ndim), intent(out)  :: a

! Output


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

  integer imu
  integer inu

  complex sqlami
  complex a0
  complex a1

  complex, allocatable, dimension (:,:) :: work
  real, allocatable, dimension (:) :: rwork
  real, allocatable, dimension (:) :: eigv
  complex, allocatable, dimension (:,:) :: bra
  complex, allocatable, dimension (:,:) :: ket

  integer lwork
  integer info

! Procedure
! ===========================================================================

! set aux variables
  a0 = cmplx(0.0d0,0.0d0)
  a1 = cmplx(1.0d0,0.0d0)

! allocate 
  lwork = ndim*ndim
  allocate (rwork (lwork))
  allocate (eigv (ndim))
  allocate (work (ndim, ndim))
  allocate (bra (ndim, ndim))
  allocate (ket (ndim, ndim))

! copy to output matrix
  a(:,:) = b(:,:)
!  do imu =1,ndim
!     do inu = 1,ndim
!        if (abs(a(imu,inu) - a(inu,imu)) .gt. 0.0001d0) then
!           write (*,*) ' ERROR  imu =',imu,' inu =', inu
!           write (*,*) ' diff =', abs(a(imu,inu)-a(inu,imu)),a(imu,inu)
!        endif
!     enddo
!  enddo

  if ( iprn ) write (*,*) 'Calling zheev subroutine ...'
  call zheev ('V', 'U', ndim, a, ndim, eigv, work, lwork, rwork, info)

  if ( info .gt. 0) then
     write (*,*) ' Error in zgeev subroutine'
     write (*,*) '   info = ',info   
     stop
  endif
      
  if (iprn) then
!     write (*,*) '                      eigenvalues:                    '
!     write (*,*) '======================================================'
!     write (*,200) (eigv(imu), imu=1,ndim)

! save eigenvectors u
! ket |u>
     ket = a
! bra <u|
     bra = Transpose(a)
     bra = Conjg(bra)

     work = matmul(ket,bra)
     write (*,*) ' Check orthonormality of eigenvectors ...'
     do i = 1,ndim
        if (abs(1.0d0 - real(work(i,i))) .gt. 0.001d0 ) then 
           write (*,*) '         ******  WARNNING  ******        '
           write (*,*) '  ',i,'-th vector is not orthonormal !!! '
           write (*,*) '  norm value should be 1.0 ', real(work(i,i))
        endif
     enddo
  endif

  do imu = 1, ndim
     if (eigv(imu) .lt. 0.0d0) then
        write (*,*) '         ******  WARNNING  ******        '
        write (*,300) imu,eigv(imu)
        write (*,*) '             It will be set zero         '
        eigv(imu) = 0.0d0
     endif
     sqlami = a1*eigv(imu)**(0.25d0)
     a(:,imu) = a(:,imu)*sqlami
  end do

  if (iprn) write (*,*) 'Calling zgemm subroutine ...'
  call zgemm ('N', 'C', ndim, ndim, ndim, a1, a, ndim, a, ndim, a0,   &
     &           work, ndim)

! copy A^(1/2)
  a(:,:) = work(:,:)

! deallocate
  deallocate (rwork)
  deallocate (eigv)
  deallocate (work)
  deallocate (bra)
  deallocate (ket)

! Format Statements
! ===========================================================================
!100   format ( <ndim>f8.3 )
!100   format ( 50f8.3 )
200   format (4(2x, f12.4))
300   format ('  ',i4,'-th eigenvalue is negative = ',f14.6)
  return

end subroutine c_sqrt

