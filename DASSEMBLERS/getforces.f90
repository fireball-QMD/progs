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

! getforces.f90
! Program Description
! ===========================================================================
!       This routine assemble forces
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getforces ()

        use options

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

!CHROM: classical forces given by a potential
		if ( iclassicMD > 0 ) then
			call getforces_classic ()
			return
		else
	        write (*,*) '  '
    	    write (*,100)
        	write (*,*) ' Now we are assembling forces. '
	        write (*,*) '  '
		endif
!END CHROM

! doing extended Hubbard
        if (ihubbard .eq. 1) then
         call getforces_eh ()
         return
        endif

! doing Kohn-Sham
        if (iKS .eq. 1) then
         call getforces_KS ()
         return
        endif

! doing Horsfield data
        if (itheory_xc .eq. 0) then
          call getforces_hxc ()
          return
        endif

! doing McWeda data
        if (itheory_xc .ne. 0) then
          call getforces_mcweda ()
          return
        endif



! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getforces

