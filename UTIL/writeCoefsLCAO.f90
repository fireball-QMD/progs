!
!
!      This code project molecular orbitals on plane waves in 1D line
!      this is usefull to plot "bandstructure" of finite systems
!         made by Prokop Hapala       ProkopHapala@gmail.com
!

 subroutine writeCoefsLCAO ()

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

!  COMPLEX, PARAMETER :: ai = (0.0,1.0)


 integer iatom
 integer in1, issh, lmu, l, l1
 integer imu, mmu, j

 real dot
 COMPLEX eik
 real, dimension (3) :: vec
 integer iband, ikpoint

 COMPLEX, DIMENSION(9) :: Ci,Ci2
 integer,dimension(3)  :: lstart=(/0,1,4/)

 character(4)    :: kstr

! Procedure
! ===========================================================================

write (*,*) " you are in writeCoefsLCAO  "

 ! =========  Band projection loop  ===========

do ikpoint = 1,nkpoints
	write (kstr,'(i4.4)') ikpoint
    open (unit = (70+1), file = 'phik_'//kstr//'_s.dat'    , status = 'unknown')
    open (unit = (70+2), file = 'phik_'//kstr//'_py.dat'   , status = 'unknown')
    open (unit = (70+3), file = 'phik_'//kstr//'_pz.dat'   , status = 'unknown')
    open (unit = (70+4), file = 'phik_'//kstr//'_px.dat'   , status = 'unknown')
    open (unit = (70+5), file = 'phik_'//kstr//'_dxy.dat'  , status = 'unknown')
    open (unit = (70+6), file = 'phik_'//kstr//'_dyz.dat'  , status = 'unknown')
    open (unit = (70+7), file = 'phik_'//kstr//'_dz2.dat'  , status = 'unknown')
    open (unit = (70+8), file = 'phik_'//kstr//'_dxz.dat'  , status = 'unknown')
    open (unit = (70+9), file = 'phik_'//kstr//'_dx2y2.dat', status = 'unknown')

    open (unit = (90+1), file = 'phie_'//kstr//'_s.dat'    , status = 'unknown')
    open (unit = (90+2), file = 'phie_'//kstr//'_py.dat'   , status = 'unknown')
    open (unit = (90+3), file = 'phie_'//kstr//'_pz.dat'   , status = 'unknown')
    open (unit = (90+4), file = 'phie_'//kstr//'_px.dat'   , status = 'unknown')
    open (unit = (90+5), file = 'phie_'//kstr//'_dxy.dat'  , status = 'unknown')
    open (unit = (90+6), file = 'phie_'//kstr//'_dyz.dat'  , status = 'unknown')
    open (unit = (90+7), file = 'phie_'//kstr//'_dz2.dat'  , status = 'unknown')
    open (unit = (90+8), file = 'phie_'//kstr//'_dxz.dat'  , status = 'unknown')
    open (unit = (90+9), file = 'phie_'//kstr//'_dx2y2.dat', status = 'unknown')

	do j =1,9
		write (70+j,*) natoms,norbitals, efermi
		write (90+j,*) natoms,norbitals, efermi
	end do
	
	do iband = 1, norbitals
		do j =1,9
			write (70+j,'(f10.5)',advance='no') eigen_k(iband,ikpoint)
			write (90+j,'(f10.5)',advance='no') eigen_k(iband,ikpoint)
		end do
		do iatom = 1, natoms
    		    in1 = imass(iatom)
        	    imu = 0
		    Ci (:) = 0
		    Ci2(:) = 0
        	    dot  = dot_product ( special_k(:,ikpoint) , ratom(:,iatom) )
		    !eik  = cmplx( cos(dot), sin(dot) )  
        	    l1 = -1
        		do issh = 1,nssh(in1)
         		    l = lssh(issh,in1)
         			do lmu = 1, (2*l+1)
          			    imu     = imu + 1
           			    mmu     = imu + degelec(iatom)
				    j       = lstart(l+1) + lmu
				    !Ci  (j) = cmplx( bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint) ) 
				    !Ci2 (j) = Ci(j) *eik   
					If (l .gt.  l1) then
					    Ci  (j) = cmplx( bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint) ) 
					else
					    Ci2 (j) = cmplx( bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint) ) 
					end if
				!write (*,'(A,5i5,2f10.5)') " DEBUG : ", iatom, issh, l, lmu, j, real(Ci(j)), aimag(Ci(j))
				end do ! lmu
			    l1 = max(l,l1)
			    !l1 = l
			end do ! l
			do j =1,9
				write (70+j,'(f14.5, f10.5)', advance='no')  real( Ci (j) ),aimag( Ci (j) )  
				write (90+j,'(f14.5, f10.5)', advance='no')  real( Ci2(j) ),aimag( Ci2(j) )  
			end do	  
		end do  ! iatom
		do j =1,9
			write (70+j,*)
			write (90+j,*)
		end do
 	end do ! iband
	do j =0,9
		close ( 70 + j )
		close ( 90 + j )
	end do  
end do ! ikpoint

   return
 end subroutine writeCoefsLCAO

!  for clusters:
!
!      This code project molecular orbitals on plane waves in 1D line
!      this is usefull to plot "bandstructure" of finite systems
!         made by Prokop Hapala       ProkopHapala@gmail.com
!

 subroutine writeCoefsLCAO_G ()

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

!  COMPLEX, PARAMETER :: ai = (0.0,1.0)


 integer iatom
 integer in1, issh, lmu, l, l1
 integer imu, mmu, j

 real dot
 COMPLEX eik
 real, dimension (3) :: vec
 integer iband, ikpoint

 real, DIMENSION(9) :: Ci,Ci2
 integer,dimension(3)  :: lstart=(/0,1,4/)
 real nichts
 character(4)    :: kstr

! Procedure
! ===========================================================================

write (*,*) " you are in writeCoefsLCAO_G  "

 ! =========  Band projection loop  ===========

ikpoint = 1
nichts = 0.0
	!write (kstr,'(i4.4)') ikpoint
    open (unit = (70+1), file = 'phik_s.dat'    , status = 'unknown')
    open (unit = (70+2), file = 'phik_py.dat'   , status = 'unknown')
    open (unit = (70+3), file = 'phik_pz.dat'   , status = 'unknown')
    open (unit = (70+4), file = 'phik_px.dat'   , status = 'unknown')
    open (unit = (70+5), file = 'phik_dxy.dat'  , status = 'unknown')
    open (unit = (70+6), file = 'phik_dyz.dat'  , status = 'unknown')
    open (unit = (70+7), file = 'phik_dz2.dat'  , status = 'unknown')
    open (unit = (70+8), file = 'phik_dxz.dat'  , status = 'unknown')
    open (unit = (70+9), file = 'phik_dx2y2.dat', status = 'unknown')

    open (unit = (90+1), file = 'phie_s.dat'    , status = 'unknown')
    open (unit = (90+2), file = 'phie_py.dat'   , status = 'unknown')
    open (unit = (90+3), file = 'phie_pz.dat'   , status = 'unknown')
    open (unit = (90+4), file = 'phie_px.dat'   , status = 'unknown')
    open (unit = (90+5), file = 'phie_dxy.dat'  , status = 'unknown')
    open (unit = (90+6), file = 'phie_dyz.dat'  , status = 'unknown')
    open (unit = (90+7), file = 'phie_dz2.dat'  , status = 'unknown')
    open (unit = (90+8), file = 'phie_dxz.dat'  , status = 'unknown')
    open (unit = (90+9), file = 'phie_dx2y2.dat', status = 'unknown')

	do j =1,9
		write (70+j,*) natoms,norbitals, efermi
		write (90+j,*) natoms,norbitals, efermi
	end do
	
	do iband = 1, norbitals
		do j =1,9
			write (70+j,'(f10.5)',advance='no') eigen_k(iband,ikpoint)
			write (90+j,'(f10.5)',advance='no') eigen_k(iband,ikpoint)
		end do
		do iatom = 1, natoms
    		    in1 = imass(iatom)
        	    imu = 0
		    Ci (:) = 0
		    Ci2(:) = 0
        	    dot  = dot_product ( special_k(:,ikpoint) , ratom(:,iatom) )
		    !eik  = cmplx( cos(dot), sin(dot) )  
        	    l1 = -1
        		do issh = 1,nssh(in1)
         		    l = lssh(issh,in1)
         			do lmu = 1, (2*l+1)
          			    imu     = imu + 1
           			    mmu     = imu + degelec(iatom)
				    j       = lstart(l+1) + lmu
				    !Ci  (j) = cmplx( bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint) ) 
				    !Ci2 (j) = Ci(j) *eik   
					If (l .gt.  l1) then
					    Ci  (j) = bbnkre(mmu,iband,ikpoint)
					else
					    Ci2 (j) = bbnkre(mmu,iband,ikpoint)
					end if
				!write (*,'(A,5i5,2f10.5)') " DEBUG : ", iatom, issh, l, lmu, j, real(Ci(j)), aimag(Ci(j))
				end do ! lmu
			    l1 = max(l,l1)
			    !l1 = l
			end do ! l
			do j =1,9
				!write (*,*) "Cij=", Ci (j)
				write (70+j,'(f14.5, f5.1)', advance='no')  Ci (j), 0.0
				write (90+j,'(f14.5, f5.1)', advance='no')  Ci2(j), 0.0
				!write (90+j,'(f14.5)', advance='no')  Ci2(j)
			end do	  
		end do  ! iatom
		do j =1,9
			write (70+j,*)
			write (90+j,*)
		end do
 	end do ! iband
	do j =0,9
		close ( 70 + j )
		close ( 90 + j )
	end do  

   return
 end subroutine writeCoefsLCAO_G

