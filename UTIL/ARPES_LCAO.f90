!
!
!      This code project molecular orbitals on plane waves in 1D line
!      this is usefull to plot "bandstructure" of finite systems
!         made by Prokop Hapala       ProkopHapala@gmail.com
!

 subroutine ARPES_LCAO ()

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use density
   use charges
   use kpoints

   implicit none

!Output
! Local Parameters and Data Declaration
! ===========================================================================

   integer iatom
   integer in1
   integer imu
   integer mmu

 real    absI
 complex Integral

 real,    dimension (:), allocatable :: norm_i

 real, dimension (3) :: vec

 integer iband, ikpoint
 integer norbin, iorbin

  character(5)    :: filenstr
 
! file name variables
   character(16)   :: namewf
   character(5)    :: name2
!   real, dimension (:), pointer   :: pmat
!   character (len=30) mssg


! ========= kpoints varibles

 integer decompose
 integer nCell !, iCell 
 integer nE, iE
 real Estep, E
 real Emax,Emin
 real, dimension    (:,:,:), allocatable :: Emap
 real, dimension    (:,:  ), allocatable :: EmapTot

 integer itype,ntypes
 integer, dimension (:), allocatable ::  type_z,type_orb


! Procedure
! ===========================================================================

! ========= start READ INPUTS ============

write (*,*) " you are in ARPES_LCAO  "

  open (unit = 30, file = 'ARPES_LCAO.optional', status = 'unknown')
    
    read (30,*) Emin, Emax, nE      
	Estep = (Emax-Emin)/float(nE)
    read (30,*) ntypes, decompose
	allocate (type_z(ntypes))
	allocate (type_orb(ntypes))
    do itype = 1,ntypes
    	read (30,*) type_z(itype),type_orb(itype)
    end do
    close(30)


 ! =========  Band projection loop  ===========

norbin = 0
do iband = 1, norbitals
	do ikpoint = 1,nkpoints
       if( (eigen_k(iband,ikpoint) .lt. Emax) .AND. ( eigen_k(iband,ikpoint) .gt. Emin) ) norbin = norbin+1
	end do
end do

write (*,*)   " DEBUG: norbin   ",norbin

iorbin = 0
allocate (  EmapTot        ( nE, nkpoints          ) )
EmapTot = 0.0
if (decompose .gt. 0) then
	allocate (  Emap ( nE, nkpoints, ntypes  ) )
	Emap = 0.0
end if


Write (*,*) "projected states:"
Write (*,*) " i   k    #    e(i,k) "
do iband = 1, norbitals
	do ikpoint = 1,nkpoints
   		if( (eigen_k(iband,ikpoint) .lt. Emax) .AND. ( eigen_k(iband,ikpoint) .gt. Emin) ) then
       		iorbin = iorbin+1

   			iE = 0
			E = Emin
    		do while ( E .lt. eigen_k(iband,ikpoint) )
    			E = E  + Estep
				iE = iE + 1
    		end do

    		Write (*,'(A, 4i10,2f10.5)') "iband,ikpoint,iorbin,iE, E:   ", iband,ikpoint,iorbin,iE, eigen_k(iband,ikpoint), (Estep*iE + Emin)

			do itype = 1,ntypes
				Integral  = 0.0
				do iatom = 1, natoms
					in1 = imass(iatom)
					if ( nzx(in1) .eq. type_z(itype) ) then
    	        		imu = type_orb(itype)
    	        		mmu = imu + degelec(iatom)
    	        		Integral  = Integral +  cmplx( bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint) ) 	  
					end if ! type_z(itype) 
				end do  ! iatom
				absI = abs(Integral)
        		if (decompose .gt. 0) Emap    ( iE, ikpoint , itype) =  Emap    ( iE, ikpoint , itype) +  absI
        		EmapTot ( iE, ikpoint )        =  EmapTot ( iE, ikpoint )        +  absI
    		end do ! itype
		end if ! in Energy range?
	end do ! ikpoint
enddo ! iband


! ========  Write out  =========

 write (*,*)   " DEBUG: norbin *  ",norbin

 open (unit = 34, file = 'ARPES_Emap_all.dat', status = 'unknown')
 do iE = 1,nE
	!write  (34,'(f10.5, 6x, <nkpoints>f18.10)')  (Estep*iE + Emin),   ( EmapTot ( iE, ikpoint ), ikpoint=1,nkpoints )
	write  (34,'(f10.5)',advance='no') (Estep*iE + Emin)
	do ikpoint=1,nkpoints
	    write  (34,'(f18.10)',advance='no') EmapTot ( iE, ikpoint )
	end do
	write (34,*)
 end do ! norbin
 close (34)

if (decompose .gt. 0) then
 do itype = 1, ntypes
	write (filenstr,'(i2.2,A,i2.2)') type_z(itype),"_",type_orb(itype)
	Write (*,*)  filenstr

	open (unit = 34, file = 'ARPES_Emap_'//filenstr//'.dat', status = 'unknown')
	do iE = 1,nE
  		!write  (34,'(f10.5, 6x, <nkpoints>f18.10)')  (Estep*iE + Emin),   ( Emap ( iE, ikpoint, itype ), ikpoint=1,nkpoints )
  	    write  (34,'(f10.5)',advance='no') (Estep*iE + Emin)
	    do ikpoint=1,nkpoints
	        write  (34,'(f18.10)',advance='no') Emap ( iE, ikpoint, itype )
	    end do
	    write (34,*)
	end do ! norbin
	close (34)
 end do ! itype
end if ! decompose

! =============== END OF BIG CYCLE ============

  if (decompose .gt. 0) deallocate ( Emap )
  deallocate ( EmapTot )

! Format Statements
! ===========================================================================
100 format (3f14.6,4x,f14.6)
200 format (e14.6)

   return
 end subroutine ARPES_LCAO

