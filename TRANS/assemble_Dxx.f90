! assemble_Dxx.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine assemble_Dxx ()
   
   use dimensions
   use transport
   use matmult
   use interactions

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output


! Local Parameters and Data Declaration
! ===========================================================================

   real, parameter :: PI = 3.141592653589793238462643d0
   real, parameter :: units = 4.0d0 * PI * PI
   real, parameter :: Go2nAmp = 7.7480162E-5

! Local Variable Declaration and Description
! ===========================================================================

   integer ie
   integer norb1
   integer norb2
   integer natm1
   integer natm2
   integer inu
   integer jnu
   integer imu
   integer jmu
   integer ix
   integer i
   integer j
   integer in1
   integer norb

! identity
   complex, dimension (:,:), allocatable :: ident1
   complex, dimension (:,:), allocatable :: ident2
! aux matrices
   complex, dimension (:,:), allocatable :: mat1
   complex, dimension (:,:), allocatable :: mat2
! tip's DOS
   complex, dimension (:,:), allocatable :: rho_11
   complex, dimension (:,:), allocatable :: rho_22
   complex, dimension (:,:), allocatable :: sq_rho_11
   complex, dimension (:,:), allocatable :: sq_rho_22
! denominators
   complex, dimension (:,:), allocatable :: Da_11
   complex, dimension (:,:), allocatable :: Dr_22
   complex, dimension (:,:), allocatable :: Dr_11
   complex a0
   complex a1

   real omega
   real vali
   real valr
   real cond
   real cond_tot
   real cur
   real cur_tot

! Lapack variables
   integer            :: info
   integer            :: lwork1
   integer            :: lwork2
   complex, dimension(:),  allocatable               :: work1
   complex, dimension(:),  allocatable               :: work2
   integer, dimension(:),  allocatable               :: ipiv1
   integer, dimension(:),  allocatable               :: ipiv2


! CALCULATION OF THE EIGENVECTORS AND EIGENVALUES OF THE CUURENT MATRIX
! Note: we will use the subroutines from Lapack library () to solve 
! eigenproblem, use only real part of the final current matrix. 
! For details about Lapack library see http://nacphy.physics.orst.edu/lapack/
! or http://www.netlib.org/lapack


! Temporary vector (dim >= 4*NOTIP1)
  real, dimension (:), allocatable  :: rwork
! Vector of eigenvalues
  complex, dimension (:), allocatable  :: eign
! Vector of eigenvalues real
  real, dimension (:), allocatable  :: reign
! Left eigenvector                 
  complex, dimension (:,:), allocatable  :: bra 
! Right eigenvector 
  complex, dimension (:,:), allocatable  :: ket 
! Temprary vector
  complex, dimension (:), allocatable  :: work 
! Matrix of current for given energy, after contains eigenvectors
  complex, dimension (:,:), allocatable  :: curm 
! Transmission matrix
  complex, dimension (:,:), allocatable  :: tm 
  complex, dimension (:,:), allocatable  :: tma
! Auxiliar vector of eigenvalues
  complex, dimension (:,:), allocatable  :: alpha
! Auxiliar vector of eigenvalues
  complex, dimension (:,:), allocatable  :: beta
  complex, dimension (:,:), allocatable  :: cbeta
  complex, dimension (:,:), allocatable  :: beta2


  integer itmp
  integer icurform

! Procedure
! ===========================================================================
 
   write (*,*) '  '
   write (*,*) ' Assemble Denominators D_xx  '
   icurform = 2

! dimension
   norb1 = sample1%norb_tip
   norb2 = sample2%norb_tip
   natm1 = sample1%natom_tip
   natm2 = sample2%natom_tip

!   norb1 = sample1%norb_tip*(2*sample1%ncell+1)**2
!   norb2 = sample2%norb_tip*(2*sample2%ncell+1)**2
!   natm1 = sample1%natom_tip*(2*sample1%ncell+1)**2
!   natm2 = sample2%natom_tip*(2*sample2%ncell+1)**2

! aux variable
   a1 = (1.0d0, 0.0d0)
   a0 = (0.0d0, 0.0d0)

   cur_tot = 0.0d0
   cond_tot = 0.0d0

! Identity matrices
   allocate ( ident1(norb1,norb1) )
   allocate ( ident2(norb2,norb2) )

! Lapack variables
   lwork1 = norb1*norb1
   lwork2 = norb2*norb2
   allocate(work1(lwork1))
   allocate(work2(lwork2))
   allocate(ipiv1(norb1))
   allocate(ipiv2(norb2))

! set up ident matrices
   ident1 = a0
   ident2 = a0
   do inu = 1, norb1
      ident1(inu,inu)  = a1
   enddo
   do inu = 1, norb2
      ident2(inu,inu)  = a1
   enddo

! auxiliar matrices
   allocate ( mat1(norb1,norb1) )
   allocate ( mat2(norb2,norb2) )
! tip dos
   allocate ( rho_11(norb1,norb1) )
   allocate ( rho_22(norb2,norb2) )
   allocate ( sq_rho_11(norb1,norb1) )
   allocate ( sq_rho_22(norb2,norb2) )
! denominators
   allocate ( Da_11(norb1,norb1) )
   allocate ( Dr_22(norb2,norb2) )
   allocate ( Dr_11(norb1,norb1) )


! arrays for EIGENproblem
  itmp = 6*norb1                       ! auxillary dimension for temporary vector WORK 
  allocate ( eign(norb1) )             ! eigenvalues
  allocate ( reign(norb1) )            ! eigenvalues (real)
  allocate ( bra(norb1,norb1) )      ! left eigenvectors
  allocate ( ket(norb1,norb1) )      ! right eigenvectors
  allocate ( rwork(itmp) )             ! tmp vector
  allocate ( curm(norb1,norb1) )       ! current matrix for given energy level
  allocate ( tm(norb1,norb2) )         ! transmission matrix 
  allocate ( tma(norb2,norb1) )         ! transmission matrix 
  allocate ( work(itmp) )              ! tmp vector
  allocate ( alpha(norb1,norb1) )      ! aux eigenvectors
  allocate ( beta(norb2,norb1) )       ! aux eigenvectors
  allocate ( cbeta(norb1,norb2) )       ! aux eigenvectors
  allocate ( beta2(norb2,norb2) )       ! aux eigenvectors

! open file eigenchannels
  open ( unit = 77, file = 'channels.dat', status = 'unknown')
! open file eigenvectors
  open ( unit = 78, file = 'ket.dat', status = 'unknown')
  open ( unit = 79, file = 'bra.dat', status = 'unknown')
  open ( unit = 80, file = 'analysis.dat', status = 'unknown')
! open current file
  open ( unit = 7, file = 'current.dat', status = 'unknown')
! open conductance file
  open ( unit = 8, file = 'conductance.dat', status = 'unknown')

! set intial energy
   omega = Elow + 0.5d0*real(dE)

   rho_11 = a0
   rho_22 = a0

   if (icurform .eq. 1) then 

      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) '+  CALCULATION of the CURRENT '
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) ' Initially we use following expression:'
      write(*,*) '' 
      write(*,*) 'G = (4*pi*e)/h sum_E ( t_12 rho_22 Dr_22 t_21 rho_11 Da_11 )'
      write(*,*) ''
      write(*,*) 'where'
      write(*,*) ' Dr_22 = 1 / [ I - t_21 Gr_11 t_12 Gr_22 ] '
      write(*,*) ' Da_11 = 1 / [ I - t_12 Ga_22 t_21 Ga_11 ] '
      write(*,*) ''
      write (*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
   else 

      write (*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      write (*,*) '+  CALCULATION of the CURRENT '
      write (*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
      write (*,*) ' Farther we use following expression:'
      write (*,*) '' 
      write (*,*) ' J = (4*pi*e)/h sum_E{ rho_11^(1/2)}*Da_11*T_12*    '
      write (*,*) '                  rho_22*T_21*Dr_11*rho_11^(1/2)}   '
      write (*,*) ''
      write (*,*) ' where'
      write (*,*) '  Da_11 = 1 / [1 - T_12*Ga_22*T_21*Ga_11]'
      write (*,*) '  Dr_11 = 1 / [1 - Gr_11*T_12*Gr_22*T_21]'
      write (*,*) ''
      write (*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++'
   endif

! energy loop
   do ie = 1,nE

      write (*,900) ie, omega

! advanced Green's function of tips
      Ga_tip1(:,:,ie) = Conjg (Transpose (Gr_tip1(:,:,ie)))
      Ga_tip2(:,:,ie) = Conjg (Transpose (Gr_tip2(:,:,ie)))

! DOS tip1
      do inu = 1,norb1
         do imu = 1,norb1
            rho_11 (inu,imu) = a1*(-1.0d0/PI)*Imag(Gr_tip1(inu,imu,ie))
!            rho_11 (inu,imu) = (-1.0d0/PI)*Imag(Gr_tip1(inu,imu,ie)) -    &
!   &            eta*(1.0d0/PI)*( real(Gr_tip1(inu,imu,ie))                 &
!   &            - real(Gr_tip1(inu,imu,ie-1)) )/real(dE)
         enddo
      enddo
! DOS tip2
      do inu = 1,norb2
         do imu = 1,norb2 
            rho_22 (inu,imu) = a1*(-1.0d0/PI)*Imag(Gr_tip2(inu,imu,ie)) 
!            rho_22 (inu,imu) = (-1.0d0/PI)*Imag(Gr_tip2(inu,imu,ie)) -      &
!   &            (eta/PI)*( real(Gr_tip2(inu,imu,ie))                       &
!   &            - real(Gr_tip2(inu,imu,ie-1)) )/real(dE)
         enddo
      enddo


! write down hopping matrix
      if (iwrt_trans) then

         write (*,*) ' rho_11 '
         do inu = 1,norb1
            do imu = 1,norb1
               write (110,200) inu, imu, real(rho_11(inu,imu))
            enddo ! imu
         enddo ! inu
         
         write (*,*) ' rho_22 '
         do inu = 1,norb2
            do imu = 1,norb2
               write (210,200) inu, imu, real(rho_22(inu,imu))
            enddo ! imu
         enddo ! inu
         
         write (*,*) ' Gr_tip1 '
         do inu = 1,norb1
            do imu = 1,norb1
               write (120,200) inu, imu, real(Gr_tip1(inu,imu,ie)),imag(Gr_tip1(inu,imu,ie))   
            enddo ! imu
         enddo ! inu
      
         write (*,*) ' Ga_tip1 '
         do inu = 1,norb1
            do imu = 1,norb1
               write (220,200) inu, imu, real(Ga_tip1(inu,imu,ie)),imag(Ga_tip1(inu,imu,ie))   
            enddo ! imu
         enddo ! inu
         
         write (*,*) ' Gr_tip2 '
         do inu = 1,norb2
            do imu = 1,norb2
               write (130,200) inu, imu, real(Gr_tip2(inu,imu,ie)),imag(Gr_tip2(inu,imu,ie))   
            enddo ! imu
         enddo ! inu
         write (*,*) ' Ga_tip2 '
         do inu = 1,norb2
            do imu = 1,norb2
               write (230,200) inu, imu, real(Ga_tip2(inu,imu,ie)),imag(Ga_tip2(inu,imu,ie))   
            enddo ! imu
         enddo ! inu
         
      endif

! ---------------------------------------
!                Da_11
! ---------------------------------------
      mat1(1:norb1,1:norb1) =                                               &
     &            ( (t_12(1:norb1,1:norb2) .x. Ga_tip2(1:norb2,1:norb2,ie)) &
     &          .x. (t_21(1:norb2,1:norb1) .x. Ga_tip1(1:norb1,1:norb1,ie)) ) 
! inversion
      mat1 = ident1 - mat1
! old fashion
!!      call inv(mat1, mat1, norb1, norb1)
! LU matrix factorization of general matrix
      call zgetrf (norb1, norb1, mat1, norb1, ipiv1, info )
      if (info .ne. 0)  then
       write (*,*) ' ***  error in zgetrf  '
       write (*,*) ' ***  info = ',info
       stop
      endif
! matrix inversion of general matrix
      call zgetri (norb1, mat1, norb1, ipiv1, work1, lwork1, info)
      if (info .ne. 0) then
       write (*,*) ' ***  error in zgetri  '
       write (*,*) ' ***  info = ',info
       stop
      endif

      Da_11(1:norb1,1:norb1) = mat1(1:norb1,1:norb1)

! ---------------------------------------
!                Dr_11
! ---------------------------------------
      mat1(1:norb1,1:norb1) =                                               &
     &            ( (Gr_tip1(1:norb1,1:norb1,ie) .x. t_12(1:norb1,1:norb2)) &
     &          .x. (Gr_tip2(1:norb2,1:norb2,ie) .x. t_21(1:norb2,1:norb1)) ) 
!! inversion
      mat1 = ident1 - mat1
! old fashion
!!      call inv(mat1, mat1, norb1, norb1)
! LU matrix factorization of general matrix
      call zgetrf (norb1, norb1, mat1, norb1, ipiv1, info )
      if (info .ne. 0)  then
       write (*,*) ' ***  error in zgetrf  '
       write (*,*) ' ***  info = ',info
       stop
      endif
! matrix inversion of general matrix
      call zgetri (norb1, mat1, norb1, ipiv1, work1, lwork1, info)
      if (info .ne. 0) then
       write (*,*) ' ***  error in zgetri  '
       write (*,*) ' ***  info = ',info
       stop
      endif

      Dr_11(1:norb1,1:norb1) = mat1(1:norb1,1:norb1)

! ---------------------------------------
!                Dr_22
! ---------------------------------------
      mat2(1:norb2,1:norb2) =                                                &
     &            ( (t_21(1:norb2,1:norb1) .x. Gr_tip1(1:norb1,1:norb1,ie)) &
     &          .x. (t_12(1:norb1,1:norb2) .x. Gr_tip2(1:norb2,1:norb2,ie)) ) 
! inversion
      mat2 = ident2 - mat2
! old fashion
!!      call inv(mat2, mat2, norb2, norb2)
! LU matrix factorization of general matrix
      call zgetrf (norb2, norb2, mat2, norb2, ipiv2, info )
      if (info .ne. 0)  then
       write (*,*) ' ***  error in zgetrf  '
       write (*,*) ' ***  info = ',info
       stop
      endif
! matrix inversion of general matrix
      call zgetri (norb2, mat2, norb2, ipiv2, work2, lwork2, info)
      if (info .ne. 0) then
       write (*,*) ' ***  error in zgetri  '
       write (*,*) ' ***  info = ',info
       stop
      endif

      Dr_22(1:norb2,1:norb2) = mat2(1:norb2,1:norb2)
   

      if (iwrt_trans) then
         write (*,*) ' Da_11 ',norb1
         do inu = 1,norb1
            do imu = 1,norb1
               write (140,200) inu, imu, real(Da_11(inu,imu)),imag(Da_11(inu,imu))
            enddo ! imu
         enddo ! inu
         write (*,*) ' Dr_11 ',norb1
         do inu = 1,norb1
            do imu = 1,norb1
               write (270,200) inu, imu, real(Dr_11(inu,imu)),imag(Dr_11(inu,imu))
            enddo ! imu
         enddo ! inu
         
         write (*,*) ' Dr_22 ',norb2
         do inu = 1,norb2
            do imu = 1,norb2
               write (240,200) inu, imu, real(Dr_22(inu,imu)),imag(Dr_22(inu,imu))
            enddo ! imu
         enddo ! inu
         
         write (*,*) ' t_12 ',norb1
         do inu = 1,norb1
            do imu = 1,norb2
               write (150,200) inu, imu, real(t_12(inu,imu)),imag(t_12(inu,imu))
            enddo ! imu
         enddo ! inu
      

         write (*,*) ' t_21 ',norb2
         do inu = 1,norb2
            do imu = 1,norb1
               write (250,200) inu, imu, real(t_21(inu,imu)),imag(t_21(inu,imu))
            enddo ! imu
         enddo ! inu
         
      endif


      if (icurform .eq. 1) then 
!!!!!!!!!!!!!!!!!!!!!!
! And now finally, we sum up in energies
! J = (4*pi*e)/h sum_E ( t_12 rho_22 Dr_22 t_21 rho_11 Da_11  )
! We would have to include the denominators
!     Dr_22 = 1 / [ I - t_21 Gr_11 t_12 Gr_22 ] 
!     Da_11 = 1 / [ I - t_12 Ga_22 t_21 Ga_11 ] 
! For that we would need to store the advanced and retarded Greens func.
! We multiply by 4 * pi * pi to get the conductance in units of the
! quantum of conductance (2 e^2 / h)
!
!         UNITS = 4. * 3.1415 * 3.1415 

! t_12 rho_22 Dr_22 t_21 rho_11 Da_11
         Jc(1:norb1,1:norb1,ie) = ( t_12(1:norb1,1:norb2) .x.              &
     &    ( rho_22(1:norb2,1:norb2) .x.  ( Dr_22(1:norb2,1:norb2) .x.      &
     &    ( t_21(1:norb2,1:norb1) .x. ( rho_11(1:norb1,1:norb1) .x.        &
     &      Da_11(1:norb1,1:norb1) ) ) ) ) )


         curm(1:norb1,1:norb1) = Jc(1:norb1,1:norb1,ie)
!  CGEEV computes for an N-by-N complex nonsymmetric matrix  A,  the  
!  eigenvalues and, optionally, the  left and/or right eigenvectors.  
!  The right eigenvector v(j) of A satisfies
!        A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!        u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm equal 
!  to 1 and largest component real.
! 
!  VARIABLES:
!  EIGN   ... eigenvalues
!  BRA  ... u(j) left eigenvectors [u(j) = BRA(:,j), 
!                  the j-th column of BRA]
!  KET  ... v(j) right eigenvectors [v(j) = KET(:,j), 
!                  the j-th column of KET]

         if (iwrt_trans) write (*,*) 'Calculate eigenchannels ...'
         call zgeev ('V', 'V', norb1, curm, norb1, eign, bra, norb1,    &
     &                ket, norb1, work, itmp, rwork, info)

!??? check it !!
         bra = Transpose(bra)
         bra = Conjg(bra)

      else if (icurform .eq. 2) then 

!============================================================================
!  SYMMETRIC COMPLEX CURRENT MATRIX
!============================================================================
     
! Reorder Formula for current to obtain EigenChanels
! Original Version:
!
! J = (4*pi*e)/h sum_E ( T_12 rho_22 D_r T_21 rho_11 D_a  )
!
!  where
!     D_r = 1 / [ I - T Gr_11 Tt Gr_22 ] 
!     D_a = 1 / [ I - T Ga_22 Tt Ga_11 ]
!
! NEW REORDERED version
!  
! sum_E{ rho_11^(1/2)}*Dr_11*T_12*rho_22*T_21*Da_11*rho_11^(1/2)}
!
!     Da_11 = 1 / [1 - T_12*Ga_22*T_21*Ga_11]
!     Dr_11 = 1 / [1 - Gr_11*T_12*Gr_22*T_21]
!

! calc rho_XX^(1/2)     
         write (*,*) ' Calculate rho_11^(1/2) matrix ....'
         call c_sqrt(rho_11, sq_rho_11, norb1, iwrt_trans)
         
         write (*,*) ' Calculate rho_22^(1/2) matrix ....'
         call c_sqrt(rho_22, sq_rho_22, norb2, iwrt_trans)

         curm(1:norb1,1:norb1) = ( sq_rho_11(1:norb1,1:norb1) .x.          &
     &    ( Da_11(1:norb1,1:norb1) .x.  ( t_12(1:norb1,1:norb2) .x.         &
     &    ( rho_22(1:norb2,1:norb2) .x. ( t_21(1:norb2,1:norb1) .x.         &
     &      (Dr_11(1:norb1,1:norb1) .x. sq_rho_11(1:norb1,1:norb1) ) ) ) ) ) )


         Jc(1:norb1,1:norb1,ie) = curm(1:norb1,1:norb1)

!  ZHEEV computes for an N-by-N complex  Hermitian  matrix  A,  the  
!  eigenvalues and, optionally, the  left and/or right eigenvectors.  
!  The right eigenvector v(j) of A satisfies
!        A * v(j) = lambda(j) * v(j)
!  where lambda(j) is its eigenvalue.
!  The left eigenvector u(j) of A satisfies
!        u(j)**H * A = lambda(j) * u(j)**H
!  where u(j)**H denotes the conjugate transpose of u(j).
!
!  The computed eigenvectors are normalized to have Euclidean norm equal 
!  to 1 and largest component real.
! 
!  VARIABLES:
!  EIGN   ... eigenvalues
!  BRA  ... u(j) left eigenvectors [u(j) = BRA(:,j), 
!                  the j-th column of BRA]
!  KET  ... v(j) right eigenvectors [v(j) = KET(:,j), 
!                  the j-th column of KET]


         if (iwrt_trans) write (*,*) ' Diagonalize current matrix ...'
         call zheev ('V', 'U', norb1, curm, norb1, reign, ket, itmp,    &
     &                rwork, info)
         
         if (info .eq. 0) then 
! save eigenvectors u
! ket |u>
            ket = curm
! bra <u|
            bra = Transpose(curm)
            bra = Conjg(bra)
            eign(:) = a1*reign(:)
         else
            write (*,*) ' Error in zgeev subroutine'
            write (*,*) '   info = ',info 
            stop
         endif
         
      endif ! if (icurform .eq. 2)


!++++++++++++++++++++++++++++++++++++++++++++++++
! CALCULATION OF CHANNELS TRANSMISSION
!++++++++++++++++++++++++++++++++++++++++++++++++


! Check orthonormality of eigenvectors: tr{|u><u|} = 1 
      curm = matmul(ket,bra)
      if (iwrt_trans) write (*,*) ' Check orthonormality of eigenvectors ...'
      do i = 1,norb1
         if (abs(1.0d0 - real(curm(i,i))) .gt. 0.001d0 ) then 
            write (*,*) '         ******  WARNNING  ******        '
            write (*,*) '  ',i,'-th vector is not orthonormal !!! '
            write (*,*) '  norm value should be 1.0 ', real(curm(i,i))
         endif
      enddo

! write eigenvalues and orbitals weights into file 
     do inu = 1,norb1
          write (77,300, advance='no')  real(eign(inu))*units
     enddo
          write(77,*) ' '
      if ( ichannel ) then 
         write (*,*) ' Store eigenvectors into file ...'
! <u|
         do inu = 1,norb1
            do imu = 1,norb1
               write (79,400, advance='no') bra(imu,inu)
            enddo
               write(79,*) ' '
         enddo
! |u>
         do inu = 1,norb1
            do imu = 1,norb1
               write (78,400,advance='no') ket(inu,imu)
            enddo
               write (78,*) ' '
         enddo  

         write (80,*) '   --------   TIP_1  --------'
         do inu = 1,norb1
            write (80,500) inu,real(eign(inu))*units
            do imu = 1,norb1
               rwork(imu) = real(ket(imu,inu))**2 + imag(ket(imu,inu))**2 
            enddo
            imu = 1
            do i = 1,natm1
               in1 = imass(sample1%atom_tip(i))
               norb = num_orb(in1)
                 write (80,351,advance="no") sample1%atom_tip(i)
              do  ix=imu,imu+norb-1
                 write (80,350,advance="no") rwork(ix)
              enddo    
                 write (80,*) ' '
               imu = imu + norb 
            enddo ! enddo i
         enddo ! enddo inu

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                  Transmission Matrix t
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! transmission matrix t
         tm(1:norb1,1:norb2) = ( sq_rho_11(1:norb1,1:norb1) .x.              &
     &          ( Da_11(1:norb1,1:norb1) .x. ( t_12(1:norb1,1:norb2) .x.      &
     &          sq_rho_22(1:norb2,1:norb2) ) ) )

! transmission matrix t+
         tma(1:norb2,1:norb1) = ( sq_rho_22(1:norb2,1:norb2) .x.             &
     &          ( t_21(1:norb2,1:norb1) .x. (Dr_11(1:norb1,1:norb1) .x.       &
     &            sq_rho_11(1:norb1,1:norb1) ) ) )


! |alpha> is eigenvector of current matrix Jc(iE) 
         alpha(1:norb1,1:norb1) = ket(1:norb1,1:norb1)
      
! {t+}*|alpha> = a1*|beta>
         beta(1:norb2,1:norb1) =                                           &
     &          ( tma(1:norb2,1:norb1) .x. alpha(1:norb1,1:norb1) ) 

         cbeta = Transpose(beta)
         cbeta = Conjg(cbeta)
! a1*{t}*|beta> = |a1|^2*|alpha> 
         alpha(1:norb1,1:norb1) =                                          &
     &          ( tm(1:norb1,1:norb2) .x. beta(1:norb2,1:norb1) ) 

         if (norb1 .eq. norb2) then 

! |a1|^2|beta><beta|
            beta2(1:norb2,1:norb2) =                                       &
    &            ( cbeta(1:norb2,1:norb1) .x. beta(1:norb1,1:norb2) )

! write into 'analysis.dat' file
            write (80,*) ''
            write (80,*) '   --------   TIP_2  --------'
            inu = norb2
            jnu = norb1
            do ix = 1,norb1
               write (80,500) inu,real(beta2(inu,inu))*units
               beta(:,jnu) = beta(:,jnu) / sqrt(real(beta2(inu,inu)))     
               do imu = 1,norb2
                  rwork(imu) = real(beta(imu,jnu))**2 + imag(beta(imu,jnu))**2
               enddo
               imu = 1
               do i = 1,natm2
                  in1 = imass(sample2%atom_tip(i))
                  norb = num_orb(in1)
                 write (80,351,advance="no") sample2%atom_tip(i)
                  do jmu=imu,imu+norb-1  
                      write (80,350,advance="no") rwork(jmu)
                  enddo
                      write(80,*) ' ' 
                  imu = imu + norb 
               enddo ! enddo i
               inu = inu - 1
               jnu = jnu -1
            enddo ! enddo ix
         else
            write (*,*) '         ---------   WARNNING   ---------         '
            write (*,*) ' Number of orbitals of tip1 and tip2 is different,'
            write (*,*) ' skip analysis of tip2'
         endif

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      endif ! if(ichannel)

! write current for given energy
      cur = 0.0d0
      cond = 0.0d0
      do inu = 1,norb1
         cond = cond + real(Jc(inu,inu,iE))
      enddo

! adjust units
      cond = cond * units 
      cur = cond * Go2nAmp * real(dE)

! sum up total values
      cur_tot = cur_tot + cur
      cond_tot = cond_tot + cond
! write out results into files
      write (7,601) iE, omega, cur, cur_tot
      write (8,600) iE, omega, cond 
   
! write total current 
      write (*,700) iE,cond
      write (*,701) iE,cur

! increment the energy 
      omega = omega + real(dE)
   enddo ! do iE

   write (8,730) cond_tot/nE 
   write (7,705)
   write (7,710) Elow, Eup
   write (7,720) cur_tot

! close files 
    close (77)
    close (78)
    close (79)
    close (80)
    close (7)
    

! Deallocate matrices

    deallocate ( mat1 )
    deallocate ( mat2 )
    deallocate ( ident1 )
    deallocate ( ident2 )
    deallocate ( rho_11 )
    deallocate ( rho_22 )
    deallocate ( Da_11 )
    deallocate ( Dr_22 )
    
    deallocate ( eign )
    deallocate ( bra )
    deallocate ( ket )
    deallocate ( rwork )
    deallocate ( curm )
    deallocate ( work )


    deallocate ( alpha )
    deallocate ( beta )
    deallocate ( cbeta )
    deallocate ( beta2 )
    deallocate ( tm )
    deallocate ( tma )

    deallocate ( sq_rho_11 )
    deallocate ( sq_rho_22 )

   deallocate (work1)
   deallocate (work2)
   deallocate (ipiv1)
   deallocate (ipiv2)

! Format Statements
! ===========================================================================
100     format (5x,i3,2f14.6)
200     format (2i5,2f14.6)
300     format (f14.6)
350     format (f12.6)
351     format (i4)
!301    format (<norb2>f14.6)
400     format (2f10.6,3x)
!401    format (<norb2>(2f10.6,3x))
500     format (i4,'-th channel: ',f10.6)
501     format (i4,'-th channel: ',3f10.6)
600     format (i4,2f16.8)
601     format (i4,3e18.8)
700     format (' Step :  ',i4,'  Go = ',f14.6,' [2*e^2/h]')
701     format (' Step :  ',i4,'  I = ',e18.8,' [A]')
705     format (' ======================================')
710     format (' Bias volatge :',2f12.6,' [eV]')
720     format (' Ic =',e18.8,'  [A]')
730     format (' Go =',e18.8,'  [2*e^2/h]')
800     format (9f12.6)
900     format (2x,'Step = ', i6,'  Energy = ',f16.8)
   return
 end subroutine assemble_Dxx
