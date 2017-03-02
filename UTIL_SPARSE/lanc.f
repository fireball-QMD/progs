      subroutine lanc1(opt,hbar,numh,listh,
     .                   nhmax,ncols,ncolsmax,
     .                   norbitals,ener)
! *********************************************************************
! Routine to calculate the mimimum or maximum eigenvalues of
! a given sparse Hamiltonian, by (2nd order) the Lanczos Method.
!
! Written by Maider Machado and P.Ordejon,  June'98
! ******************************** INPUT ******************************
! integer opt                     : 1 = compute minimun eigenval.
!                                   2 = compute maximum eigenval.
! real hbar(nhmax,ncolsmax)       : Sparse Hamiltonian
! integer numh(nhmax)             : control vector of hbar
! integer listh(nhmax,ncolsmax)   : control vector of hbar
! integer nhmax                   : Maximum number of nonzero elements of
!                                   each column of the hamiltonian
! integer ncols                   : number of columns on this processor
! integer ncolsmax                : maximum number of columns
! integer norbitals               : number of basis orbitals
! ******************************* OUTPUT *******************************
! real ener                     : Eigenvalue
! **********************************************************************

      implicit none

! BTN communication domain
       integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
       common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

      integer nhmax, ncols, ncolsmax, norbitals, opt

      integer listh(nhmax,ncolsmax), numh(ncolsmax)

      real ener, hbar(nhmax,ncolsmax)

!  Internal variables ...

      integer itmax
      parameter (itmax=200)

      real tol
      parameter (tol=0.0001d0)

      integer i, ii, j, k, ierror

      real
     .  a0, as0, b1, b2, bs2, c1, diff, ecomp, eivec(2), 
     .  hv0(norbitals), hvs0(norbitals), mod0, mod, norm, norms, ran3,
     .  u0(ncols), v0(norbitals), vs0(norbitals), wr(2)

      integer mpic_real, mpic_double, mpic_max, mpic_sum
      common /mpiconsts/ mpic_real, mpic_double, mpic_max, mpic_sum

! ...

!  An unlikely number ...
      ecomp=4321.0987d0
! ...

!  Generate an initial normalized random vector ........
      call random_seed
      mod0=0.
      do i=1,ncols
        call random_number(ran3)
        u0(i)=2.*ran3-1.
        mod0=mod0+u0(i)**2 
      enddo
      call mpi_allreduce (mod0, mod, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
      mod=sqrt(mod)
      do i=1,ncols
        u0(i)=u0(i)/mod
      enddo
! ...

      do i=1,norbitals
        v0(i)=0.
      enddo

!  Lanczos loop ................
      do k=1,itmax
        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              v0(ii)=v0(ii)+hbar(i,j)*u0(j)
           enddo
        enddo
        call mpi_allreduce (v0, vs0, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        a0=0.
        do i=1,ncols
           a0=a0+u0(i)*vs0(i)
        enddo
        call mpi_allreduce (a0, as0, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        norm=0.
        do i=1,ncols
           vs0(i)=vs0(i)-as0*u0(i)
           norm=norm+vs0(i)**2
        enddo
        call mpi_allreduce (norm, norms, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        b1=sqrt(norms)

        do i=1,norbitals
           hv0(i)=0.0
        enddo

        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              hv0(ii)=hv0(ii)+hbar(i,j)*vs0(j)
           enddo
        enddo
        call mpi_allreduce (hv0, hvs0, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        b2=0.
        do j=1,ncols
           b2=b2+u0(j)*hvs0(j)/b1
        enddo
        call mpi_allreduce (b2, bs2, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        c1=0.
        do i=1,norbitals
           c1=c1+vs0(i)*hvs0(i)/norms
        enddo


! eigenvalues and eigenvectors ...

        wr(1) = 0.5d0*((as0+c1) 
     .          + sqrt((as0+c1)**2 - 4.0d0*(as0*c1-b1*bs2)))
        wr(2) = 0.5d0*((as0+c1) 
     .          - sqrt((as0+c1)**2 - 4.0d0*(as0*c1-b1*bs2)))

        eivec(1)=1/sqrt(1+((as0+b1-wr(opt))/(bs2+c1-wr(opt)))**2)
        eivec(2)=-eivec(1)*(as0+b1-wr(opt))/(bs2+c1-wr(opt))

        ener=wr(opt)
! ...

 
        diff=abs((wr(opt)-ecomp)/wr(opt))
        if (diff.gt.tol) then
          norm=0.
          do i=1,ncols
             u0(i)=eivec(1)*u0(i)+eivec(2)*vs0(i)/b1
             norm=norm+u0(i)**2
          enddo
          call mpi_allreduce (norm, norms, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
          do i=1,norbitals
             v0(i)=0.
          enddo
          do j=1,ncols
            u0(j)=u0(j)/sqrt(norms)
          enddo
          ecomp=wr(opt)
        else 
          goto 20
        endif
      enddo

      write(*,*) 'WARNING: lanc1 not converged after ',itmax,
     .           ' iterations'

20    return
      end





      subroutine lanc2(hbar,numh,listh,nhmax,
     .                 ncols,ncolsmax,norbitals,icstart,ener)
! *********************************************************************
! Routine to calculate the eigenvalue closest to zero for a given
! sparse Hamiltonian H, using the the Folded Spectrum Method
! (by (2nd order) the Lanczos Method).
!
! Written by Maider Machado and P.Ordejon,  June'98
! ******************************** INPUT ******************************
! real hbar(nhmax,ncolsmax)       : Sparse Hamiltonian
! integer numh(nhmax)             : control vector of hbar
! integer listh(nhmax,ncolsmax)   : control vector of hbar
! integer nhmax                   : Maximum number of nonzero elements of
!                                   each column of the hamiltonian
! integer ncols                   : number of columns on this processor
! integer ncolsmax                : maximum number of columns
! integer norbitals               : number of basis orbitals
! integer icstart                 : starting column on this processor
! ******************************* OUTPUT *******************************
! real ener                     : Eigenvalue of H
! **********************************************************************

      implicit none

! BTN communication domain
       integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
       common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

      integer nhmax,ncols,ncolsmax,norbitals,icstart

      integer listh(nhmax,ncolsmax), numh(ncolsmax)

      real ener, hbar(nhmax,ncolsmax)

!  Internal variables ...

      integer itmax
      parameter (itmax=500)

      real tol
      parameter (tol=0.000001d0)

      integer i, ii, ij, j, k

      integer ierror

      real
     .  a0, b1, b2, bs2, c1, cs1, diff, ecomp, eivec(2), 
     .  hv0(ncols), mod0, mod, norm, u0(norbitals),
     .  us0(norbitals), v0(norbitals), vs0(norbitals),
     .  v00(norbitals), vs00(norbitals), wr, ran3

      integer mpic_real, mpic_double, mpic_max, mpic_sum
      common /mpiconsts/ mpic_real, mpic_double, mpic_max, mpic_sum

! ...


!  An unlikely number ...
      ecomp=4321.0987d0
! ...

      do i=1,norbitals
        u0(i)=0.
      enddo

!  Generate an initial normalized random vector ........
      mod0=0.
      call random_seed
      do i=icstart,icstart+ncols-1
        call random_number(ran3)
        u0(i)=2.*ran3-1.
        mod0=mod0+u0(i)**2 
      enddo
      call mpi_allreduce (mod0, mod, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
      call mpi_allreduce (u0, us0, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
      mod=sqrt(mod)
      do i=1,norbitals
        us0(i)=us0(i)/mod
      enddo
! ...

      do i=1,norbitals
        v0(i)=0.
        v00(i)=0.
      enddo

!  Lanczos loop ................
      do k=1,itmax
        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              v00(icstart+j-1)=v00(icstart+j-1)+hbar(i,j)*us0(ii)
           enddo
        enddo
        call mpi_allreduce (v00, vs00, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              v0(icstart+j-1)=v0(icstart+j-1)+hbar(i,j)*vs00(ii)
           enddo
        enddo
        call mpi_allreduce (v0, vs0, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        a0=0.
        do i=1,norbitals
           a0=a0+us0(i)*vs0(i)
        enddo
        norm=0.
        do i=1,norbitals
           vs0(i)=vs0(i)-a0*us0(i)
           norm=norm+vs0(i)**2
        enddo
        b1=sqrt(norm)

        do i=1,ncols
           hv0(i)=0.0
        enddo
        do i=1,norbitals
          v00(i)=0.0
        enddo

        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              v00(icstart+j-1)=v00(icstart+j-1)+hbar(i,j)*vs0(ii)
           enddo
        enddo
        call mpi_allreduce (v00, vs00, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        do j=1,ncols
           do i=1,numh(j)
              ii=listh(i,j)
              hv0(j)=hv0(j)+hbar(i,j)*vs00(ii)
           enddo
        enddo
        b2=0.
        do j=1,ncols
           b2=b2+us0(icstart+j-1)*hv0(j)/b1
        enddo
        call mpi_allreduce (b2, bs2, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
        c1=0.
        do i=1,ncols
           c1=c1+vs0(icstart+i-1)*hv0(i)/norm
        enddo
        call mpi_allreduce (c1, cs1, 1,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)


!  minimum eigenvalue ...
        wr = 0.5d0*((a0+cs1) 
     .                 - sqrt((a0+cs1)**2 - 4.0d0*(a0*cs1-b1*bs2)))
! ...

! eigenvector ...
        eivec(1)=1/sqrt(1+((a0+b1-wr)/(bs2+cs1-wr))**2)
        eivec(2)=-eivec(1)*(a0+b1-wr)/(bs2+cs1-wr)
! ...

 
        norm=0.

        do i=1,norbitals
          u0(i)=eivec(1)*us0(i)+eivec(2)*vs0(i)/b1
          norm=norm+u0(i)**2
        enddo
        do i=1,norbitals
          v0(i)=0.
          v00(i)=0.
          u0(i)=u0(i)/sqrt(norm)
        enddo

        diff=abs(wr-ecomp)
        ecomp=wr
        if (diff.lt.tol) goto 20
      enddo
! .................

      write(*,*) 'WARNING: lanc2 not converged after ',itmax,
     .           ' iterations'

20    continue

! Calculate eigehvalue of H ...
      ener=0.0
      do i=1,norbitals
        v0(i)=0.0
      enddo
      do i=1,ncols
        do j=1,numh(i)
          ij=listh(j,i)
          v0(icstart+i-1)=v0(icstart+i-1)+hbar(j,i)*u0(ij)
        enddo
      enddo
      call mpi_allreduce (v0, vs0, norbitals,
     .                    mpic_real, mpic_sum,
     .                    mpi_btn_world, ierror)
      do i=1,norbitals
        ener=ener+u0(i)*vs0(i)
      enddo
      
      return
      end

