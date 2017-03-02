! copyright info:
!
! @Copyright 2010
! Fireball Committee
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
!used in RGL and MEAM potentials 
! ===========================================================================
! Code written by:
! Zdenka Chromcova
! Institute of Physics oftheAS CR,v. v. i.
! Cukrovarnicka 10
! CZ-162 00 Praha 6
! Czech Republic
! email: chrom@fzu.cz
! webpage with description: http://nanosurf.fzu.cz/wiki/doku.php?id=classical_md
! ===========================================================================


module classicMEAMforces
interface
subroutine getForceMEAM(i,potential,nspecies,ratom,imass,N,msRho,Rho0,distance,Energ,Force)
	use classicMD, only: potentialParam
	type(PotentialParam),intent(in), dimension(nspecies,nspecies) :: potential
	integer, intent(in) :: N,i,nspecies,imass(N)
	real, dimension(3,N),intent(in) :: ratom
	real, intent(in) :: distance(N,N)
	real, intent(out) :: msRho(N,4),Energ,Force(3)
	real,intent(inout) :: Rho0(N,4)
end subroutine getForceMEAM

function distance2(a,b)
	real, intent(in) :: a(3), b(3)
	real :: distance2
end function distance2

function distanceNeigh(i,ineigh,ratom,N)
    integer, intent(in) :: N,i,ineigh
    real,intent(in) :: ratom(3,N)
    real :: distanceNeigh
end function distanceNeigh

function Cijk(ri,rj,rk,Cmin,Cmax)
    real, intent(in) :: ri(3),rj(3),rk(3),Cmin,Cmax
    real :: Cijk
end function Cijk			    

subroutine dC(ri,rj,rk,Cmin,Cmax,dCijk)
	real, intent(out) :: dCijk(3)
	real, intent(in) :: ri(3),rj(3),rk(3),Cmin,Cmax
end subroutine dC	

subroutine SdS(i,jneigh,Cmin,Cmax,ratom,N,dSij,Sij)
	real, intent(out) :: dSij(3),Sij
	real, intent(in) :: Cmin,Cmax
	integer, intent(in) :: i,jneigh,N	
	real, dimension(3,N),intent(in) :: ratom
end subroutine SdS

function x(i,jneigh,a,ratom,distance,N)
	integer, intent(in) :: i,jneigh,a,N
	real, dimension(3,N),intent(in) :: ratom
	real,intent(in) :: distance(N,N)
	real :: x
end function x

function dx(i,jneigh,a,ratom,distance,N,di)
	integer, intent(in) :: i,jneigh,a,di,N
	real, dimension(3,N),intent(in) :: ratom
	real, intent(in) :: distance(N,N)
	real :: dx
end function dx

function dR(di,ineigh,a,ratom,R,N)
	integer, intent(in) :: di,ineigh,N,a
	real, dimension(3,N),intent(in) :: ratom
	real, intent(in) :: R
	real :: dR
end function dR

subroutine Rho(i,Cmin,Cmax,ratom,imass,distance,N,Rho2ij,dRho2ij,beta,rho0,R0)
	real, intent(in) :: Cmin,Cmax,beta(0:3),rho0,R0, distance(N,N)
	real, intent(out) :: Rho2ij(0:3),dRho2ij(0:3,3)
	integer, intent(in) :: N,i,imass(N)
	real, dimension(3,N),intent(in) :: ratom
end subroutine Rho

subroutine gammad(i,Cmin,Cmax,ratom,imass,distance,N,beta,t,gamma,dgamma,Rho0,dRho0,rho0Param,R0)
	real, intent(in) :: Cmin,Cmax,beta(4),t(3),Rho0Param,R0,distance(N,N)
	integer, intent(in) :: N,i,imass(N)
	real, dimension(3,N),intent(in) :: ratom
	real, intent(out) :: gamma,dgamma(3),Rho0,dRho0(3)
end subroutine gammad

subroutine sRho(i,imass,Cmin,Cmax,ratom,distance,N,beta,Rho0Param,R0,tParam,msRho,mRho0)
	real, intent(in) :: Cmin,Cmax,beta(4),Rho0Param,tParam(3),R0,distance(N,N)
	integer, intent(in) :: N,i,imass
	real, dimension(3,N),intent(in) :: ratom
	real, intent(out) :: msRho(1:4),mRho0(1:4)
end subroutine sRho

subroutine backgroundDensity(beta,tt,R,R0,Z,dR,rho0,backRho,dBackRho)
	real,intent(in) :: beta(0:3),tt(3),R,dR(3),R0,rho0
	integer, intent(in) :: Z
	real, intent(out) :: backRho, dbackRho(3)
end subroutine backgroundDensity

function countNeigh(i,distance,natom,R0)
	real, intent(in) :: distance(natom,natom),R0
	integer, intent(in) :: i, natom
	integer :: countNeigh
end function countNeigh
end interface
end module classicMEAMforces
