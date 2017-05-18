!
!
!      This code project molecular orbitals on plane waves in 1D line
!      this is usefull to plot "bandstructure" of finite systems
!         made by Prokop Hapala       ProkopHapala@gmail.com
!

 subroutine ew2mesh_ARPES (icluster)

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints

   implicit none
   INCLUDE 'fftw3.f'  


! Argument Declaration and Description
! ===========================================================================
! Input
   integer, intent (in) :: icluster

!Output
! Local Parameters and Data Declaration
! ===========================================================================

   INTEGER, PARAMETER :: DPC = 8
!   COMPLEX, PARAMETER :: ai = (0.0,1.0)

   interface
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  end interface

   interface
    subroutine  project_orb_complex(iband,ikpoint,ewfaux)
     integer iband
     integer ikpoint
     complex, dimension (:), pointer, intent (out) :: ewfaux
    end
  end interface

! Local Variable Declaration and Description
! ===========================================================================

! ========== Grid varibales
   complex, target, dimension (:), allocatable :: ewfaux
   complex, dimension (:), pointer   :: pewf

   integer index
   integer i, j, k
   real x,y,z

 real Integral1,Integral0
 complex Integral2,kwave
 real dot
 COMPLEX ai,a0

 real, dimension    (:,:,:), allocatable :: results
 real, dimension    (:,:), allocatable :: results_depth
 real,    dimension (:), allocatable :: results_E
 real,    dimension (:), allocatable :: norm_i

 integer iband, ikpoint
 integer norbin, iorbin
 
! file name variables
   character(16)   :: namewf
   character(5)    :: name2
!   real, dimension (:), pointer   :: pmat
!   character (len=30) mssg


! ========= kpoints varibles

 integer nkps,    ikp

 real alat
 real kstep

 real, dimension (:,:), allocatable  :: kline

 real time_start,time_end

 integer byEnerg 
 integer, ALLOCATABLE, DIMENSION (:) :: myBands
 real Emax,Emin
 real valMax


! Procedure
! ===========================================================================

! ========= start READ INPUTS ============

write (*,*) " YOU ARE IN ew2mesh_ARPES "

  open (unit = 30, file = 'ARPES.optional', status = 'unknown')      
    read (30,*) Emin, Emax
    read (30,*) alat                           ! read lattice constant
    read (30,*) nkps                   ! read 
    write (*,*) "DEBUG: ", nkps, alat
    allocate(kline(3,nkps))

    read (30,*) kline(:,1), kline(:,nkps)   
    write (*,'(A,6f10.5)') "DEBUG: ",kline(:,1), kline(:,nkps) 
    kline(:,1)    =  kline(:,1 )*    2*3.14159265/alat
    kline(:,nkps) =  kline(:,nkps)*  2*3.14159265/alat
    close(30)

  ! interpolate k-points in lines

write (*,*) "k-points to project on: " 
    kstep = 1.0/float(nkps)
       do ikp = 2,nkps-1
         kline(:,ikp) =   kstep*(nkps-ikp+1)*kline(:,1) + kstep*(ikp-1)*kline(:,nkps)  
         write (*,'(i10,3f10.5)')   ikp,kline(:,ikp) 
       end do ! ipoint

write (*,*) " elvec :"
write (*,'(3f10.5)') elvec(1,:)
write (*,'(3f10.5)') elvec(2,:)
write (*,'(3f10.5)') elvec(3,:)

! ========= end READ INPUTS ============

   ai = cmplx(0.0d0, 1.0d0)
   a0 = cmplx(0.0d0, 0.0d0)

   allocate ( ewfaux(0:nrm-1))
   Integral0 = rm2*rm1    ! norm of grid

   Write (*,*)   " Grid norm =  ", Integral0

 ! =========  Band projection loop  ===========

open (unit = 33, file = "ARPES.dat", status = 'unknown')

write (*,*) " DEBUG: 1 "

norbin = 0
do iband = 1, norbitals
	do ikpoint = 1,nkpoints
       if( (eigen_k(iband,ikpoint) .lt. Emax) .AND. ( eigen_k(iband,ikpoint) .gt. Emin) ) norbin = norbin+1
	end do
end do

write (*,*)   " DEBUG: norbin   ",norbin

allocate (  results         ( norbin, nkps, rm3  )  )
allocate (  results_depth   ( norbin, rm3  )  )
allocate (  results_E       ( norbin) )
allocate (  norm_i          ( norbin) )

iorbin = 0

Write (*,*) "projected states:"
Write (*,*) " i   k    #    e(i,k) "
do iband = 1, norbitals
   do ikpoint = 1,nkpoints
       if( (eigen_k(iband,ikpoint) .lt. Emax) .AND. ( eigen_k(iband,ikpoint) .gt. Emin) ) then
       iorbin = iorbin+1
	   results_E ( iorbin) = eigen_k(iband,ikpoint)
       Write (*,'(3i10,f10.5)') iband,ikpoint,iorbin, results_E ( iorbin)

     pewf => ewfaux
     call project_orb_complex(iband,ikpoint,pewf)    

! =======================================
!              Start SCAN  
! =======================================

 call cpu_time (time_start)

! write (30,'(2x,f20.10)',advance='no')   eigen_k(iband,1)

! evaluate z-density
index = 0
norm_i (iorbin) = 0
do k = 1, rm3
  Integral1 = 0.0
   do j = 1, rm2
    do i = 1, rm1
      Integral1 = Integral1 + ewfaux(index)*conjg(ewfaux(index))
      index = index + 1
    enddo ! i
   enddo ! j
 results_depth ( iorbin, k ) = Integral1
 norm_i (iorbin) = norm_i(iorbin) + Integral1
enddo ! k
results_depth ( iorbin, : ) = results_depth ( iorbin, : ) /  norm_i (iorbin)

 do ikp = 1,nkps 

! Evaluate projection
   index = 0
   do k = 1, rm3
    Integral2 = a0
    !Integral1 = 0.0
    !Integral0 = 0.0
    do j = 1, rm2
     do i = 1, rm1
      x = (i-1)*elvec(1,1) + (j-1)*elvec(2,1)
      y = (i-1)*elvec(1,2) + (j-1)*elvec(2,2)
      z = (k-1)*elvec(3,3)
      dot = ( kline(1,ikp)*x + kline(2,ikp)*y )  
      kwave =  cos(dot) + ai*sin(dot) 
      Integral2 = Integral2 + (ewfaux(index)*kwave)
      Integral1 = Integral1 +  ewfaux(index)*conjg(ewfaux(index))
      Integral0 = Integral0 +  kwave*conjg(kwave)
      index = index + 1
     enddo ! i
    enddo ! j
     results ( iorbin, ikp, k  )  = abs(Integral2) / sqrt( Integral0 * norm_i (iorbin) )
    ! write (*,'(A,i5,3f16.8)')   "DEBUG: Integral1, Integral2: ", k, Integral1, abs(Integral2), results ( iorbin, ikp, k  )

    !write (*,'(A,3i5,5f16.8)')   "DEBUG: iorbin,ikp,k, |kw|, |wf| , Im(P), Re(P, |P|)  ", iorbin, ikp, k,  sqrt(Integral0), sqrt(Integral1),  aimag(Integral2), real(Integral2), results ( iorbin, ikp, k  )
   enddo ! k

end do ! ikp

call cpu_time (time_end)
Write (*,*) "         time : ", (time_end - time_start), " [sec]"

!  write (30,*)

  end if ! in Energy range?
 end do ! ikpoint
enddo ! iband

! ========  Write out  =========

write (*,*)   " DEBUG: norbin *  ",norbin

do iorbin = 1,norbin
write ( 33, '(A, f10.5)' ) "               ",results_E ( iorbin)
 do k = 1, rm3
  write  (33,'(f16.10)', advance='no') results_depth ( iorbin, k )
  do ikp = 1, nkps
    write  (33,'(f16.8)', advance='no') results ( iorbin, ikp, k )
  end do
  write (33,*)
  !write  (33,'(f16.10,6x,<nkps>f14.10)') results_depth ( iorbin, k ), ( results ( iorbin, ikp, k ), ikp=1,nkps )
 end do
end do


!do iorbin = 1,norbin
!write ( 33, '(A, f10.5)' ) "               ",results_E ( iorbin)
!do ikp = 1, nkps
! do k = 1, rm3
!  write  (33,'(3i5,f16.8)') iorbin, ikp, k, results ( iorbin, ikp, k )
!  end do
! write (33,*)
! end do
!end do


! =============== END OF BIG CYCLE ============


   deALLOCATE ( results_E )
   deALLOCATE ( results )
   deALLOCATE ( kline   )
   deallocate (ewfaux)

  close (33)
  close (30)

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ew2mesh_ARPES

