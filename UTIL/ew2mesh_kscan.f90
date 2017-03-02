!
!
!      This code project molecular orbitals on plane waves in 1D line
!      this is usefull to plot "bandstructure" of finite systems
!         made by Prokop Hapala       ProkopHapala@gmail.com
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                         !!!
!!!   WARRNING:  This code works only for rectangular box and icluster=1    !!!
!!!                                                                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 subroutine ew2mesh_kscan (icluster)

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
    subroutine  project_orb(iband,ewfaux)
     integer iband
     real, dimension (:), pointer, intent (out) :: ewfaux
    end
  end interface

! Local Variable Declaration and Description
! ===========================================================================

! ========== Grid varibales
   real, target, dimension (:), allocatable :: ewfaux
   real, dimension (:), pointer   :: pewf

   integer index
   integer i, j, k
   real x,y,z

! lookup tables  
 complex     ,dimension (:,:,:), allocatable :: eikxi 
 complex     ,dimension (:,:,:), allocatable :: eikyj
 complex     ,dimension (:,:,:), allocatable :: eikzk

 complex     ,dimension (:), allocatable :: eikxiii 
 complex eikyjj,eikzkk

 complex Integral, Integral2,kwave
 real dot
 COMPLEX ai

 integer iband,ipband
 
! file name variables
   character(16)   :: namewf
   character(5)    :: name2
!   real, dimension (:), pointer   :: pmat
!   character (len=30) mssg


! ========= kpoints varibles

 integer nkps,    ikp
 integer nlines,  iline

 real alat
 real kstep

 real, dimension (:,:,:), allocatable  :: klines
 complex, dimension (:,:)  , allocatable  :: linesVal
 !complex, dimension (:,:)  , allocatable  :: linesVal2


 real time_start,time_end

 integer byEnerg 
 integer, ALLOCATABLE, DIMENSION (:) :: myBands
 real Emax,Emin
 real valMax


! Procedure
! ===========================================================================

! ========= start READ INPUTS ============

  open (unit = 30, file = 'kscan.optional', status = 'unknown')
    read (30,*) byEnerg                        
    read (30,*) Emin, Emax
    read (30,*) alat                           ! read lattice constant
    read (30,*) nlines, nkps                   ! read 
    write (*,*) "DEBUG: ",nlines, nkps, alat
    allocate(klines(3,nkps,nlines))
    allocate(linesVal(nkps,nlines))
  !  allocate(linesVal2(nkps,nlines))    

    do iline = 1,nlines
    write (*,*) iline
    read (30,*) klines(:,1,iline), klines(:,nkps,iline)   
    write (*,'(A,6f16.8)') "DEBUG: ",klines(:,1,iline), klines(:,nkps,iline) 
    klines(:,1  ,iline)  =  klines(:,1  ,iline)*    2*3.14159265/alat
    klines(:,nkps,iline) =  klines(:,nkps,iline)*   2*3.14159265/alat
    end do
    close(30)

  ! interpolate k-points in lines
    kstep = 1.0/float(nkps)
    do iline = 1,nlines 
       do ikp = 2,nkps-1
         klines(:,ikp,iline) =   kstep*(nkps-ikp+1)*klines(:,1,iline) + kstep*(ikp-1)*klines(:,nkps,iline)  
         write (*,'(2i10,3f10.5)')    iline,ikp,klines(:,ikp,iline) 
       end do ! ipoint
    end do ! iline

write (*,*) " elvec :"
write (*,'(3f10.5)') elvec(1,:)
write (*,'(3f10.5)') elvec(2,:)
write (*,'(3f10.5)') elvec(3,:)

! ========= end READ INPUTS ============

Write (*,*) "Building lookup table ...."

   ai = cmplx(0.0d0, 1.0d0)
ALLOCATE(  eikxi(rm1,nkps,nlines)  ) 
ALLOCATE(  eikyj(rm2,nkps,nlines)  ) 
ALLOCATE(  eikzk(rm3,nkps,nlines)  ) 

ALLOCATE(  eikxiii(rm1)            ) 

do iline = 1,nlines
 do ikp = 1,nkps
   do i = 1, rm1
     x = klines(1,ikp,iline)  *(i-1)*elvec(1,1) 
     eikxi(i,ikp,iline) =  cos(x) + ai*sin(x)  
   end do ! i
   do j = 1, rm2
     y = klines(2,ikp,iline)  *(j-1)*elvec(2,2)
     eikyj(j,ikp,iline) =  cos(y) + ai*sin(y) 
   end do
   do k = 1, rm3
     z = klines(3,ikp,iline)  *(k-1)*elvec(3,3)
     eikzk(k,ikp,iline) =  cos(z) + ai*sin(z)
   end do
 end do ! ikp
end do ! iline

Write (*,*) "..... lookup table done"


! ======== END build loockup table =========

! allocate aux arrays
   allocate ( ewfaux(0:nrm-1))

 ! =========  Band projection loop  ===========



allocate (myBands(norbitals))
myBands=0
if ( byEnerg .eq. 1 ) then
    do iband = 1, norbitals
       if( (eigen_k(iband,1) .lt. Emax) .AND. ( eigen_k(iband,1) .gt. Emin) ) myBands(iband)=1
    end do
else
    do i = 1, npbands
       myBands(pbands (i))=1
    end do
end if 

!
Write (*,*) "DEBUG >>>"
    do iband = 1, norbitals
       if( myBands(iband) .eq. 1 ) write (*,*) iband 
    end do
Write (*,*) "<<<< DEBUG"
!



  open (unit = 33, file = "klines_TOT.dat", status = 'unknown')

do iband = 1, norbitals
  if ( myBands(iband) .eq. 1) then
       Write (*,*) iband, eigen_k(iband,1) 

     pewf => ewfaux
     call project_orb(iband,pewf)    

 ! =======================================
 !              Start SCAN  
 ! =======================================


 call cpu_time (time_start)

! write (30,'(2x,f20.10)',advance='no')   eigen_k(iband,1)
do iline = 1,nlines
 do ikp = 1,nkps

  eikxiii(:) = eikxi(:,ikp,iline)
  
! Evaluate projection
   index = 0
   Integral  = 0.0d0
   Integral2 = 0.0d0
   do k = 1, rm3
      eikzkk = eikzk(k,ikp,iline) 
    do j = 1, rm2
       eikyjj = eikyj(j,ikp,iline)
     do i = 1, rm1

      ! This code is ~6x faster but work just for rectangular lattice      
      kwave =  eikxiii(i) * eikyjj * eikzkk    
      Integral = Integral + ewfaux(index)*kwave

   !   ! This code should work for nonrectangular lattice  
   !   x = (i-1)*elvec(1,1) + (j-1)*elvec(2,1) + (k-1)*elvec(3,1)
   !   y = (i-1)*elvec(1,2) + (j-1)*elvec(2,2) + (k-1)*elvec(3,2)
   !   z = (i-1)*elvec(1,3) + (j-1)*elvec(2,3) + (k-1)*elvec(3,3)
   !   dot = (klines(1,ikp,iline)*x + klines(2,ikp,iline)*y + klines(3,ikp,iline)*z )  
   !   kwave =  cos(dot) + ai*sin(dot) 
   !   Integral2 = Integral2 + ewfaux(index)*kwave

      index = index + 1
     enddo ! i
    enddo ! j
   enddo ! k


  linesVal(ikp,iline)   = Integral
 !  linesVal2(ikp,iline)  = Integral2

  ! write (30,'(2x,f20.10)',advance='no')  abs(Integral)

 end do ! ikp
end do ! ilines

 call cpu_time (time_end)
Write (*,*) "         time : ", (time_end - time_start), "[sec]"

!  write (30,*)

 ! =======================================
 !              END SCAN  
 ! =======================================

! ========  Write out  =========

  Write (name2,'(i5.5)') iband
  namewf = 'klines_'//name2//'.dat'
  open (unit = 31, file = namewf, status = 'unknown')
  do ikp = 1,nkps
    write (31,'(i10)',advance='no') ikp
    do iline = 1,nlines 
        write (31,'(f16.8)',advance='no') abs(linesVal(ikp,iline))
    end do ! iline
   write (31,*)
 end do ! ipoint
 close (31)

  write (33,'(f10.5)',advance='no') eigen_k(iband,1) 
  do ikp = 1,nkps
    valMax=0
    do iline = 1,nlines 
      if(valMax .lt. abs(linesVal(ikp,iline)) ) valMax=abs(linesVal(ikp,iline)) 
    end do ! iline
    write (33,'(f16.8)',advance='no') valMax
 end do ! ipoint
 write (33,*)
 

!  Write (name2,'(i5.5)') iband
!  namewf = 'klines2'//name2//'.dat'
!  open (unit = 31, file = namewf, status = 'unknown')
!  do ikp = 1,nkps
!    write (31,'(i10)',advance='no') ikp
!    do iline = 1,nlines 
!        write (31,'(f16.8)',advance='no') abs(linesVal2(ikp,iline))
!    end do ! iline
!   write (31,*)
! end do ! ipoint
! close (31)

end if ! myBands(iband) 
enddo ! iband
! =============== END OF BIG CYCLE ============


   deALLOCATE ( linesVal )
   deALLOCATE ( klines   )

   deALLOCATE ( eikxi ) 
   deALLOCATE ( eikyj )
   deALLOCATE ( eikzk )

   deallocate (ewfaux)


  close (33)
  close (30)


! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ew2mesh_kscan

