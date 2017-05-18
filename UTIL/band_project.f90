
        subroutine band_project ()
        use charges
        use configuration
!        use constants
        use density
        use dimensions
        use interactions
!        use neighbor_map
        use kpoints
        implicit none

! Local Parameters and Data Declaration
! ===========================================================================

   COMPLEX, PARAMETER :: ai = (0.0,1.0) ! imaginary unit

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ik
        integer imu,mmu
        integer in1
        integer issh
        integer i

  complex projection,kwave
  real dens

 character(5)    :: filenstr

 real Emax, Emin, E
 integer Z_choosen

 integer iband,nfound, ifound
 integer iFirst, iLast 

 integer itype,ntypes
 integer, dimension (:), allocatable :: type_z,type_orb
 real, dimension (:,:,:), allocatable :: wk
 real                                 :: wk_tot


! Procedure
! ===========================================================================
! Initialize some things


        write (*,*) ' ****************************************************** '
        write (*,*) '                   band_project                         '
        write (*,*) '  = Project all bands on selected atoms and orbitals    '    
        write (*,*) ' ****************************************************** '

! ========= start READ INPUTS ============
 
  open (unit = 27, file = 'band_project.optional', status = 'unknown')
! read energy range
   read (27,*) Emin,Emax
	read (27,*) ntypes
	allocate (type_z(ntypes))
	allocate (type_orb(ntypes))
    do itype = 1,ntypes
    	read (27,*) type_z(itype),type_orb(itype)
    end do
    close(27)

! find range of idexes iFirst, iLast, nfound
iFirst = 1
iband = 0
write (*,*) " norbitals, nkpoints ",norbitals, nkpoints 
do 
	iband = iband + 1
if ( iband .eq. norbitals) exit
	E = -1000000.0
	do ik = 1, nkpoints
		E = max(eigen_k(iband,ik),  E )
	end do
	! write (*,*) " iFirst",iband, E
if ( E .gt. Emin) exit
end do
iFirst = iband
do
	iband = iband + 1
if ( iband .eq. norbitals) exit
	E = +1000000.0
	do ik = 1, nkpoints
		E = min(eigen_k(iband,ik),  E )
	end do
	! write (*,*) " iLast",iband, E
if ( E .gt. Emax) exit
end do
iLast  = iband - 1 
nfound = iband - iFirst

Write (* ,'(3i10,4f18.10)')  nfound, iFirst, iLast, eigen_k(iFirst,1),eigen_k(iLast,1), Emin, Emax

allocate ( wk     (nfound,nkpoints,ntypes) )

  open (unit = 31, file = 'Projected_ek.dat', status = 'unknown')
  open (unit = 32, file = 'Projected_wk.dat', status = 'unknown')

Write (31,'(3i10,4f18.10)')  nfound, iFirst, iLast, eigen_k(iFirst,1),eigen_k(iLast,1), Emin, Emax
Write (32,'(3i10,4f18.10)')  nfound, iFirst, iLast, eigen_k(iFirst,1),eigen_k(iLast,1), Emin, Emax

ifound = 0
do iband = iFirst, iLast
	Write (*,'(i10,f18.10)') iband, eigen_k(iband,1)
    ifound = ifound + 1
	do ik = 1, nkpoints
		Write (31,'(f18.10)', advance='no') eigen_k(iband,ik)
		wk_tot = 0.0
		do itype = 1, ntypes
    		dens = 0.0
			do iatom = 1, natoms
				in1 = imass(iatom)
				if ( nzx(in1) .eq. type_z(itype) ) then
    	   			imu = type_orb(itype)
           			mmu = imu + degelec(iatom)
           			dens  = dens + (bbnkre(mmu,iband,ik)**2) + (bbnkim(mmu,iband,ik)**2)
				end if ! Z_choosen
			end do  ! iatom
    		wk(ifound,ik,itype) = dens
    		wk_tot = wk_tot+dens 
		end do ! itype
		Write (32,'(f18.10)', advance='no') wk_tot
	end do ! ik
    Write (31,*)
    Write (32,*)
end do ! iband
 close(31) 
 close(32)

write (*,*) " Write out Projected_wk_XX_YY.dat : "

do itype = 1, ntypes
	write (filenstr,'(i2.2,A,i2.2)') type_z(itype),"_",type_orb(itype)
	Write (*,*)  filenstr

	open (unit = 33, file = 'Projected_wk_'//filenstr//'.dat', status = 'unknown')
    Write (33,'(3i10,4f18.10)')  nfound, iFirst, iLast, eigen_k(iFirst,1),eigen_k(iLast,1), Emin, Emax
	do ifound = 1,nfound
  		!write  (33,'(<nkpoints>f18.10)') ( wk ( ifound, ik, itype ), ik=1,nkpoints )
  		do ik=1,nkpoints
  		    write  (33,'(f18.10)',advance='no')  wk ( ifound, ik, itype )
  		end do
  		write  (33,*)
	end do ! iband
	close (33)
end do ! itype

deallocate ( wk      )

  return
end subroutine band_project

