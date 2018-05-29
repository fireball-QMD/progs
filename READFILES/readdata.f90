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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! readdata.f90
! Program Description
! ===========================================================================
!       This routine reads the different data file.
!
! ===========================================================================
! Code originally written by Juergen Fritsch and Otto F. Sankey
 
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readdata ()
 
        use options 
        
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
    
 
! Procedure
! ===========================================================================
! Initialize rcutoff array - we loop over all of it elsewhere
! ============================================================================
!                              read data files
! ============================================================================

! ****************************************************************************
! About ioff2c and ioff3c - the interactions are turned on or off according to 
! the ioff2c and ioff3c arrays which are read in from file diagnostics.input. 

! Two Center Interactions:
! Note the last three short range Ewald, long range Ewald, and coulomb
! (extended Hubbard) are not interactions from datafiles like the rest.
! However, we give the option here to turn off these interactions anyways.
!                ioff - 2c overlap
!                ioff - 2c vna_ontopl
!                ioff - 2c vna_ontopr
!                ioff - 2c vna_atom-atom
!                ioff - 2c non-local
!                ioff - 2c xc_ontop
!                ioff - 2c xc_atom-atom
!                ioff - 2c xc_correction
!                ioff - 2c z-dipole
!                ioff - 2c y-dipole
!                ioff - 2c x-dipole
!                ioff - 2c coulomb
!                ioff - 2c kinetic
!                ioff - 2c extended Hubbard
!                ioff - 2c den_ontopl
!                ioff - 2c den_ontopr
!                ioff - 2c den_atom
!                ioff - 2c denS_ontopl
!                ioff - 2c denS_ontopr
!                ioff - 2c denS_atom
!                ioff - 2c_overlapS
!                ioff - 2c coulomb (extended Hubbard)
!                ioff - 2c short range Ewald
!                ioff - 2c long range Ewald
!
! Three Center Interactions:
!                ioff - 3c neutral-atom
!                ioff - 3c exchange-correlation
!                ioff - 3c average density OLSXC
!                ioff - 3c average density OLSXC (spheric)
! ****************************************************************************

!        write (*,*) '  '
!        write (*,100)
!        write (*,*) ' Now we are reading Fdata file. '
!        write (*,*) '  '
!============OLD=================
! read extended Hubbard data
!        if (ihubbard .eq. 1) then
!         call readdata_eh ()
!         return
!        endif

! read Kohn_Sham data
!        if (iks .eq. 1) then
!         call readdata_KS ()
!         return
!        endif

! read Horsfield data
!        if (itheory_xc .eq. 0) then
!          call readdata_hxc ()
!          return 
!        endif

! read McWeda data        
!        if (itheory_xc .ne. 0) then
!          call readdata_mcweda ()
!          return 
!        endif 

!CHROM
!============NEW==============

		if(iclassicMD > 0 )then
!classical solution only
		elseif (ihubbard .eq. 1) then
! read extended Hubbard data
			call readdata_eh ()
		elseif (iks .eq. 1) then
! read Kohn_Sham data
		 	call readdata_KS ()
		elseif (itheory_xc .eq. 0) then
! read Horsfield data
			call readdata_hxc ()
		elseif (itheory_xc .ne. 0) then
! read McWeda data        
			call readdata_mcweda ()
		endif

! assign the right pointer on subroutine force
!		call init_getforces()
!END CHROM

! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine readdata
 
