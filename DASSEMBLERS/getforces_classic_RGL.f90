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


! getforces_classic.f90
! Program Description
! ===========================================================================
!RGL potential according:
!Philos. Mag. A 59 (1989) 321
!J.Chem.Phys 122 (2005) 194308 R.Ferrando Ag-Cu,Ag-Ni,Au-Cu,Ag-Au,Ag-Pd,Pd-Pt,Ag-Ni
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


subroutine getRGLforce (ratom,df,etotal,distance)
    use classicMD,only: distance2, Potential,getLJforce, distanceNeigh, atom
	use neighbor_map, only: neigh_classic,neighn_classic,neigh_b_classic
	use configuration, only: xl, natoms, nspecies
	use interactions, only: imass

	implicit none
	integer :: idx
	real, dimension(3,natoms),intent(in) :: ratom
	real, intent(inout) :: distance(natoms,natoms)
	real, intent(out) :: df(3,natoms)
	double precision, intent(out) :: etotal
	integer :: i,j
	double precision :: alpha, p, q, A, dzeta, ro, Z , r1, de
	double precision :: sr1(3),dr1(3),dEr(3),Eb(natoms),Er,dEb(3),cEb,Er_ij,fbond(3)

	alpha=0.5
	etotal = 0
	df(:,:)=0.0
	
	do idx = 1,natoms
		if(potential(1,1)%nparam/=7)then
			write(*,*)'wrong number of parameters for RGL potential. I need 7 parameters: alpha, p, q, A, dzeta, Ecoh, Z .'
			write(*,*)'you gave me only: ',potential%nparam,' parameters'
			write(*,*)'exiting...'
			stop
		elseif(potential(1,1)%type(1:3)/='RGL')then
			write(*,*)'--- BAD POTENTIAL!',potential(1,1)%type(1:7),potential(1,1)%pType,' ---'
			stop
		endif

		dEr(:) = (/0.0,0.0,0.0/)
		Eb(idx) = 0
		dEb = 0
		Er = 0
                do i = 1,neighn_classic(idx)
			j = neigh_classic(i,idx)
			dr1(:) = (ratom(:,idx)-atom(idx,i,ratom,natoms))
			r1 = sqrt(dr1(1)*dr1(1)+dr1(2)*dr1(2)+ dr1(3)*dr1(3))
			sr1(:) = dr1(:)/r1
                        
			if (j/=idx .or. sum(xl(:, neigh_b_classic(i, idx))*xl(:, neigh_b_classic(i, idx))) > 0.01)then ! &
                             ! .and. r1< 2*potential(imass(j),imass(idx))%cutoff) then
                                alpha = potential(imass(j),imass(idx))%params(1)
				p = potential(imass(j),imass(idx))%params(2)
				q = potential(imass(j),imass(idx))%params(3)
				A = potential(imass(j),imass(idx))%params(4)
				dzeta = potential(imass(j),imass(idx))%params(5)
				ro = potential(imass(j),imass(idx))%params(6)
				Z = potential(imass(j),imass(idx))%params(7)
! bonding term
				cEb = (dzeta**(1.0/alpha))*dexp(-q*(r1/ro-1)/alpha)
				Eb(idx) = Eb(idx) + cEb
! repulsive Born-Mayer term
				Er_ij = A*dexp(-p*(r1/ro-1))
				Er = Er + Er_ij
				dEr(:) = dEr(:) - Er_ij*p*sr1(:)/ro
			endif
		enddo
		Eb(idx) = dsqrt(Eb(idx))
		if(alpha == 0.5)then
			de = Er - Eb(idx)
			df(:,idx) = -2*dEr(:)
		else
			de = Er - (Eb(idx)**(2*alpha))
			df(:,idx) = -2*dEr(:)
		endif
		etotal = etotal + de
	enddo

! BONDING PART OF THE FORCEij DEPENDS ALSO ON Eb(j) OF NEIGBORS:	
	do idx = 1,natoms
		fbond=0
		do i = 1,neighn_classic(idx)
			j = neigh_classic(i,idx)
			dr1(:) = (ratom(:,idx)-atom(idx,i,ratom,natoms))
			r1 = sqrt(dr1(1)*dr1(1)+dr1(2)*dr1(2)+ dr1(3)*dr1(3))
			sr1(:) = dr1(:)/r1
			if (j/=idx .or. sum(xl(:, neigh_b_classic(i, idx))*xl(:, neigh_b_classic(i, idx))) > 0.01)then !&
                             ! .and. r1< 2*potential(imass(j),imass(idx))%cutoff) then
                                alpha = potential(imass(j),imass(idx))%params(1)
				p = potential(imass(j),imass(idx))%params(2)
				q = potential(imass(j),imass(idx))%params(3)
				A = potential(imass(j),imass(idx))%params(4)
				dzeta = potential(imass(j),imass(idx))%params(5)
				ro = potential(imass(j),imass(idx))%params(6)
				Z = potential(imass(j),imass(idx))%params(7)
! bonding term
				cEb = (dzeta**(1.0/alpha))*dexp(-q*(r1/ro-1)/alpha)
				dEb(:) = - cEb*q*sr1(:)/(ro*alpha)

				if(alpha == 0.5)then
					fbond(:) = fbond(:) + dEb(:)*(0.5/Eb(idx)+0.5/Eb(j))
				else
!for alpha!=0.5 it is not the optimal solution
					fbond(:) = fbond(:) + 2*alpha*(Eb(idx)**(2*(alpha-1)) + Eb(j)**(2*(alpha-1)))*dEb(:) 
				endif
			endif
		enddo

		df(:,idx) = df(:,idx) + fbond  !units - eV/A -> nN
	enddo
!	stop
end subroutine getRGLforce

