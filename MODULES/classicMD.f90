!
!                             @Copyright 2010
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

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

! getforces_classic.f90
! Program Description
! ===========================================================================
! Definition of global variable (allocatable) Potential, user-defined types 
! and interface for subroutines and functions. 
! ===========================================================================
! Code written by:
! Zdenka Chromcova
! Institute of Physics of  the  AS CR,  v. v. i.
! Cukrovarnicka 10
! CZ-162 00 Praha 6
! Czech Republic
! email: chrom@fzu.cz
! webpage with description: http://nanosurf.fzu.cz/wiki/doku.php?id=classical_md
! ===========================================================================

module classicMD

	use configuration
!	procedure(getLJforce), pointer :: getClassForcePt
	integer :: freq_of_outputs

	type POTENTIALPARAM
		integer :: ptype
		character(10) :: type
		integer :: nparam
		real, dimension(:), allocatable :: params
		real :: cutoff
	end type POTENTIALPARAM
    
    type(POTENTIALPARAM), dimension(:,:),allocatable,target :: Potential

	interface
	subroutine readPotenialParams(path,path_len,param,locimass,nParam)
		integer, intent(in) :: path_len,nParam
		real, dimension (nParam), intent(out):: param
		character, dimension(path_len), intent(in) :: path
		integer, dimension(2), intent(in) :: locimass
	end subroutine readPotenialParams

	real function distance2(a,b)
		real, intent(in) :: a(3),b(3)
	end function

	real function distanceNeigh(i,ineigh,ratom,N)
		integer, intent(in) :: N,i,ineigh
		real,intent(in) :: ratom(3,N)
	end function distanceNeigh

	function atom(i,ineigh,ratom,N)
		real, dimension(3,N),intent(in) :: ratom
		real :: atom(3)
		integer, intent(in) ::i,ineigh,N
	end function atom
	
	integer function find_neigh_max_class(param)
		integer :: param
	end

	subroutine getLJforce (ratom,f,e,distance)
		use configuration, only: natoms
		real, intent(in) :: ratom(3,natoms)
		real,intent(out) :: f(3,natoms)
		double precision, intent(out) :: e
		real, intent(inout) :: distance(natoms,natoms)
	end subroutine

	subroutine getforce_vdw (ratom,f,e)
		use configuration, only: natoms
		real, intent(in) :: ratom(3,natoms)
		real,intent(out) :: f(3,natoms)
		double precision, intent(out) :: e
	end subroutine

	subroutine getRGLforce (ratom,f,e,distance)
		use configuration, only: natoms
		real, intent(in) :: ratom(3,natoms)
		real,intent(out) :: f(3,natoms)
		double precision, intent(out) :: e
		real, intent(inout) :: distance(natoms,natoms)
	end subroutine

	subroutine getTersoffforce (ratom,f,e,distance)
		use configuration, only: natoms
		real, intent(in) :: ratom(3,natoms)
		real,intent(out) :: f(3,natoms)
		double precision, intent(out) :: e
		real, intent(inout) :: distance(natoms,natoms)
	end subroutine
	
	end interface

end module classicMD

