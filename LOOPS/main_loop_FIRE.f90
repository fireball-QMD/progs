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

! main_loop_FIRE.f90
! Program Description
! ===========================================================================
!       This routine performs FIRE minimization
!       autor: Prokop Hapala (  hapala@fzu.cz ,  ProkopHapala@gmail.com  )
!   References:
!  Bitzek, E., Koskinen, P., Gähler, F., Moseler, M. & Gumbsch, P. Structural relaxation made simple. Phys. Rev. Lett. 97, 170201 (2006).  
!  Eidel, B., Stukowski, A. & Schröder, J. Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods. Pamm 11, 509–510 (2011).
!  DOI: 10.1002/pamm.201110246        
! ===========================================================================
! Code written by:
! P. Hapala
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: hapala@fzu.cz

! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine main_loop_FIRE ()
	use options
	use configuration
	use options
	use MD
	use forces
	use constants_fireball
	use energy
	use optimization
	use scf

	implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
	integer itime_step

! Procedure
! ===========================================================================
! in case of constrain DFT, first we perform SCF loop for the ground state;
! This will be used later as a reference state
      if(icDFT .eq. 1 ) then
       call scf_loop (0)
       call getenergy (0)
       cDFT_active = .true.
      endif  ! end icDFT 

      write (*,*) ' ======================================================= '
      write (*,*) '             Begin FIRE minimization                     '
      write (*,*) ' ======================================================= '

      call init_FIRE( )

      do itime_step = nstepi, nstepf
!  pre-step
       if (iimage .ge. 1) call imaged (icluster, iimage, itime_step, nstepi)
       if (iclassicMD == 0) then
        call scf_loop (itime_step)
        call postscf ()                 ! optionally perform post-processing (DOS etc.)
        call getenergy (itime_step)    ! calculate the total energy
        write (*,*) '  '
        write (*,*) ' ================= Move on to Forces ================= '
       end if
!  FIRE step
       call getforces ()   ! Assemble forces
       call move_ions_FIRE (itime_step)   ! Move ions now
! call move_ions (itime_step)   ! Move ions now
!  output
!       call writeOutRFE    ( ratom, ftot, etot, 'coordinates and forces        ', itime_step )
!       call writeAnswerBas ( ratom, natoms )
       call write_bas(  )
! convergence criteria
	write (*,'(A,i6,2f16.8)') ' ++++ i, Fmax, force_tol ', itime_step, deltaFmax, force_tol
!	if ( FIRE_Ftot .lt. force_tol ) then
	if ( deltaFmax .lt. force_tol ) then
	  write (*,*) ' +++++ FIRE.optionalimization converged +++++ '
	  write (*,*) 'That`sall for now, bye ..'
	  stop
	endif
      end do ! itime_step

      return
end subroutine main_loop_FIRE

! ======================================
! =============== init_FIRE( )
! ======================================

subroutine init_FIRE( )
	use optimization
	use MD
	use fragments
	use configuration
	logical file_exists
	inquire(file="FIRE.optional", exist=file_exists ) 
	! set FIRE parameters
	if ( file_exists ) then
		write (*,*) " loading parameters from FIRE.optional "
		open (unit= 33, file="FIRE.optional", status='unknown')
		read (33,*)  FIRE_finc 
		read (33,*)  FIRE_fdec 
		read (33,*)  FIRE_acoef0  
		read (33,*)  FIRE_falpha 
		read (33,*)  FIRE_Nmin 
		read (33,*)  FIRE_mass  
		read (33,*)  FIRE_Fclamp  
		read (33,*)  FIRE_dtmax
		read (33,*)  FIRE_dtmin
		close (33)   
	else
		write (*,*) " No FIRE.optional => setting default parameters "
		FIRE_finc     = 1.1D0
		FIRE_fdec     = 0.5D0
		FIRE_acoef0   = 0.1D0
		FIRE_falpha   = 0.99D0
		FIRE_Nmin     = 5        ! currently not used
		FIRE_mass     = 4.0D0
		FIRE_Fclamp   = 10.0D0   ! too big force
		FIRE_dtmax    = dt
		FIRE_dtmin    = dt*0.1
	end if
	! set FIRE varaibles
	FIRE_dt       = FIRE_dtmax * 0.25
	FIRE_acoef    = FIRE_acoef0

	if ( .not. allocated( fragxyz ) ) then
		allocate( fragxyz(3,natoms) ) 
		fragxyz( :,:) = 0
	end if

end subroutine init_FIRE

! ======================================
! =============== move_ions_FIRE( )
! ======================================

subroutine move_ions_FIRE( istep )
	use outputs, only: iwrtxyz
	use optimization
	use configuration
	use forces
	use fragments
        use energy, only: deltaFmax

 implicit none
! == arguments 
	integer, intent(in) :: istep
! == Local variables
	integer iatom, k
	real ff,vv,vf, cF, cV, dtv
 	character (200 ) xyz_headline
! == Procedure
	
! projection of force to velocity <v|f>
	ff = 0.D0
	vv = 0.D0
	vf = 0.D0
    deltaFmax = 0.0d0
	do iatom = 1, natoms
		do k=1,3
			if ( fragxyz(k,iatom) .eq. 0 ) then 
				deltaFmax = max(deltaFmax, abs(ftot(k,iatom)))
				ff = ff + ftot(k,iatom)**2
				vv = vv + vatom(k,iatom)**2
				vf = vf + vatom(k,iatom) * ftot(k,iatom)
			end if ! fragxyz(k,iatom)
		end do ! k
	end do ! iatom
	FIRE_Ftot = sqrt( ff )

! FIRE update depending of projection <v|f>
	if ( vf .lt. 0 ) then
		write (*,*) " DEBUG FIRE: <v|f>  < 0 "
		vatom(:,:)   = 0 
		!FIRE_dt      = FIRE_dt * FIRE_fdec
		FIRE_dt      = max( FIRE_dt * FIRE_fdec, FIRE_dtmin )
		FIRE_acoef   = FIRE_acoef0
    else
		cF           =     FIRE_acoef * sqrt(vv/ff)
		cV           = 1 - FIRE_acoef
		do iatom = 1, natoms
			do k=1,3
				if ( fragxyz(k,iatom) .eq. 0 ) then 
					vatom(k,iatom) = cV * vatom(k,iatom)    +    cF *  ftot(k,iatom)
				end if ! fragxyz(k,iatom)
			end do ! k
		end do ! iatom
		FIRE_dt     = min( FIRE_dt * FIRE_finc, FIRE_dtmax ) 
		FIRE_acoef  = FIRE_acoef   * FIRE_falpha
    end if

!  normal MD step using leap-frog 
	dtv  = FIRE_dt / FIRE_mass 
	if( deltaFmax .gt. FIRE_Fclamp ) then  ! if force is too big
		dtv        = dtv * FIRE_Fclamp / deltaFmax
		vatom(:,:) = vatom(:,:) * 0.5 
	end if
	do iatom = 1, natoms
		do k=1,3
			if ( fragxyz(k,iatom) .eq. 0 ) then 
				vatom(k,iatom) = vatom(k,iatom)   +  dtv     * ftot (k,iatom)
				ratom(k,iatom) = ratom(k,iatom)   +  FIRE_dt * vatom(k,iatom)
			end if
		enddo
	enddo

	write ( xyz_headline, '(A, i6, 5f16.8)' ) " #### FIRE: i,Fmax,|F|,v,<v|f>,dt: ", istep, deltaFmax, FIRE_Ftot, sqrt(vv), vf, FIRE_dt  
	write ( *, '(A)' )  xyz_headline
	if ( iwrtxyz .eq. 1 ) call write_to_xyz( xyz_headline, istep )

end subroutine move_ions_FIRE


! ======================================
!        OUTPUTS ROUTINES
! ======================================

! =============== write_bas(  )
subroutine write_bas(  )
	use configuration, only: natoms,ratom
	use interactions, only: imass
	use charges, only: nzx
 implicit none
! == variables 
	 integer iatom, in1
! == body
	open (unit = 17, file = 'answer.bas', status = 'unknown')
	write (17,*) natoms
	do iatom = 1, natoms
		in1 = imass(iatom)
		write (17, '(i5, 3f16.8)'  ) nzx(in1), ratom(:,iatom)
	enddo
	close (unit = 17)
end subroutine write_bas


! ===============  write_to_xyz(  )
subroutine write_to_xyz( headline, istep )
	use configuration, only: natoms,ratom,symbol
	use interactions, only: imass
 implicit none
! == arguments 
	character (*), intent(in) :: headline
	integer      , intent(in) :: istep
! == variables 
	 integer iatom, in1
! == body
	if ( istep .eq. 1 ) then
		open (unit = 18, file = 'answer.xyz', status = 'unknown')
	else
		open (unit = 18, file = 'answer.xyz', status = 'unknown', position = 'append')
	endif
	write (18,*) natoms
	write (18,'(A)') headline
	do iatom = 1, natoms
		in1 = imass(iatom)
		write (18, '(A, 3f16.8)' ) symbol(iatom), ratom(:,iatom)
	enddo
	close (unit = 18)
end subroutine write_to_xyz











