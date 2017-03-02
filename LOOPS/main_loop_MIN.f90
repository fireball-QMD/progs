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

! main_loop_CG.f90
! Program Description
! ===========================================================================
!       This routine performs Conjugated Gradient optimization loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_MIN ()
 
        use configuration, only: natoms,ratom,shifter,ishiftO
        use forces, only: ftot
        use energy, only: etot
        use optimization, only: cg_maxstep
        use charges, only:nzx
        use interactions, only: imass
        use options, only: iquench
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
     integer :: iatom
     integer :: istatus
         
! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                 MINIMIZATION   L O O P
! ---------------------------------------------------------------------------
! ===========================================================================
     write (*,*) '++   '
     write (*,*) '++  ======================================================= '
     write (*,*) '++              Begin BFGS minimization. '
     write (*,*) '++  ======================================================= '
     write (*,*) '++   '

!  L-BFGS-B optimization         
     call bfgs (istatus)
          
     if (istatus == 1) then
       	write(*,*)'++ =============================='
       	write(*,*)'++ ==     MINIMUM REACHED      =='
       	write(*,*)'++ =============================='
     elseif( istatus == 2 ) then 
       	write(*,*)'++ =============================='
       	write(*,*)'++ ==   MINIMIZATION FAILED    =='
       	write(*,*)'++ =============================='
       	write(*,*)'++ The BFGS algorithm did not find the next step in minimization.'
     elseif( istatus == 3 ) then
        write(*,*)'++ =============================='
        write(*,*)'++ ==   MINIMIZATION STOPPED   =='
        write(*,*)'++ =============================='
        write(*,*)'++ The minimization in the progress but you want to stop minimization after',cg_maxstep,'iteration'
    endif

! WRITE OUT OUTPUTS
	if( iquench /= -6 )then
	    open (unit= 86, file='answer.bas', status='unknown')
   		write (86,*) natoms
		write(*,*)'++ ==== ENERGY ===='
		write(*,*)'++ etot=',etot

		write(*,*)'++ ==== COORDINATES: ===='
		do iatom=1,natoms
			write(*,701)nzx(imass(iatom)),ratom(:,iatom)
        	if (ishiftO == 1) then
            	write (86,700) nzx(imass(iatom)), ratom(:,iatom) - shifter(:)
	        else
    	        write (86,700) nzx(imass(iatom)), ratom(:,iatom)
        	endif
		enddo
    	close(unit=86)
    endif
    
    700 format (2x, i2, 3(2x,f18.8))
    701 format ('++ ',2x, i2, 3(2x,f18.8))
    return
end subroutine main_loop_MIN
 
