! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

!
! readhop.f90
! Program Description
! ===========================================================================
!
! author
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine readhop ( nspecies )

   use interactions
   use transport
   use charges
   use constants_fireball

   implicit none

! Argument Declaration and Description
! ===========================================================================

   integer, intent(in) :: nspecies

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


   character (len=80)     :: title
   character (len=2)      :: wx1
   character (len=2)      :: wx2

   integer, dimension(3)  :: aux1
   integer                :: idist
   integer                :: ndist
   integer                :: in1
   integer                :: in2
   integer                :: num_nonzero
   integer                :: i
   integer                :: j
   integer                :: inter
   integer                :: ii

   real, dimension(3)     :: aux2
   real                   :: dmin
   real                   :: dmax
   real                   :: tmp
   real                   :: rcmax
   real                   :: dist
   real                   :: dx

   real, dimension (:,:), allocatable    :: hop_lr
   real, dimension (:), allocatable      :: rr


! nsh_tip(iat_tip) : # of shells (s,p,d) in tip (# of orbitals=nst_tip(iat_tip) **2
! nineq       : # of interactions =/ 0
! lsh_ts(1-3,I,J1,I1): l1 l2 m characterizing the interaction (eg 0 2 0 == sd sigma)
! Hoppings:
!    (1) hop_ts(I,J,J1,I1)  == short range hopping
!    (2) ts_hop_lr(1-3,I,J1,I1) == prefactor, power, exponent for long-range hopping
!

   write(*,*) ' read hoppings from files ...'
   in1 = sample1%nspec
   in2 = sample2%nspec
   allocate ( nzh (in1,in2) )
   allocate ( zh_min (in1,in2) )
   allocate ( zh_max (in1,in2) )

! loop over species of tip1
   do i = 1, sample1%nspec
    in1 = sample1%spec(i)
    write (wx1,'(i2.2)') nzx(in1)

! loop over species of tip2
    do j = 1, sample2%nspec
     in2 = sample2%spec(j)
     write (wx2,'(i2.2)') nzx(in2)

! open data file
     open(3,file='interaction_'//wx1//'_'//wx2//'.dat')
! read header
     read(3,'(a80)') title
! read number of points
     read(3,*) nzh(i,j)
! read min and max distance
     read(3,*) zh_min(i,j), zh_max(i,j)
! read number of nonzero elements
     read(3,*) num_nonzero

     if ( num_nonzero .ne. index_max2c(in1,in2) ) then
      write (*,*) ' number of non zero interactions is different'
      write (*,*) ' index_max2c = ',index_max2c(in1,in2)
      write (*,*) ' num_nonzero = ',num_nonzero
      write (*,*) ' file = ','interaction_'//wx1//'_'//wx2//'.dat'
      stop
     endif

! close data file
     close(3)

    end do ! do j
   end do ! do i

! maximal number of nonzero interactions
   max_int = int( maxval(index_max2c(:,:)))
! maximal number of distances
   ndist = int( maxval(nzh(:,:)))

   in1 = sample1%nspec
   in2 = sample2%nspec
! allocate arrays
   allocate ( hops ( ndist, max_int, in1, in2) )
   allocate ( hop_lr (3, max_int) )
   allocate ( lsh_hop (3, max_int, in1, in2) )
   allocate ( rc_hfit (max_int, in1, in2) )
! allocate spline table
   allocate ( hop_spline (ndist, max_int, in1, in2))

   hops = 0.0d0
   lsh_hop = 0
   hop_lr  = 0.0d0
   rc_hfit = 0.0d0


   do i = 1, sample1%nspec
    in1 = sample1%spec(i)
    write (wx1,'(i2.2)') nzx(in1)
    do j = 1, sample2%nspec
     in2 = sample2%spec(j)
     write (wx2,'(i2.2)') nzx(in2)

     write (*,*) ' Opening data file: ','interaction_'//wx1//'_'//wx2//'.dat'
! open data file
     open(3,file='interaction_'//wx1//'_'//wx2//'.dat')

! read header
     read(3,'(a80)') title
! read number of points
     read(3,*) ndist
! read min and max distance
     read(3,*) dmin, dmax
! read number of nonzero elements
     read(3,*) num_nonzero

! read fitting part
     do inter = 1, num_nonzero

! read line
      read (3,*) aux1(1:3), aux2(1:3), ii, ii, rcmax

      lsh_hop (1:3, inter, i, j) = aux1 (1:3)
      hop_lr (1:3, inter) = aux2 (1:3)
      rc_hfit (inter, i, j) = rcmax

     end do ! do inter

! read nonzero interactions
     do idist = 1, ndist
      read(3,*) dist, (hops(idist,inter,i,j), inter = 1,num_nonzero)
! convert hoppings from a.u. to eV
      hops(idist,:,i,j) = hops(idist,:,i,j)*Hartree
     end do

! close data file
     close(3)

! matching long range part into hoppings
     dist = zh_min(i,j)
     dx = real((zh_max(i,j) - zh_min(i,j)) / (nzh(i,j)-1))
! allocate rr arrays (store distances)
     allocate (rr(ndist))

     do idist = 1,ndist
      do inter = 1, num_nonzero
       if (dist .ge. rc_hfit(inter,i,j)) then
        aux2(1:3) = hop_lr(1:3,inter)
        tmp = aux2(1) * exp ( - aux2(3) * dist - aux2(2) * log(dist) )
! convert to eV
        hops(idist,inter,i,j) = tmp*Hartree
       endif
      end do ! do inter

! store & update actual distance
      rr(idist) = dist
      dist = dist + dx
     end do ! do idist

! convert distances into Angstrom units
     rr = rr*abohr
     zh_min(i,j) = zh_min(i,j)*abohr
     zh_max(i,j) = zh_max(i,j)*abohr

! loop over nonzero values
     do inter = 1, num_nonzero

! call spline builder
      call buildspline2_1d (rr(1), hops(1,inter,i,j), ndist,         &
                            hop_spline(1,inter,i,j))

     end do ! do inter

! test
     do idist = 1,ndist
      write (101,*) rr(idist), (hops(idist,inter,i,j), inter = 1,num_nonzero)
     end do

! deallocate dr
     deallocate ( rr )

    end do ! do j
   end do ! do i

   deallocate ( hop_lr )

! Format Statements
! ===========================================================================

   return

   end subroutine readhop
