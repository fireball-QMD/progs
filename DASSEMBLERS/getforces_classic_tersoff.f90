! copyright info:
!
! @Copyright 2005
! Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio State University - Dave Drabold
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
!===============================================================
!Tersoff potential for the classical MD simulation - it is usable mainly for
! Si, Ge atoms, params are in files i.e. tersof_Si-Si.dat
!source: PRB37(1988)6991
!COMMENT TO ACCTUAL IMPLEMENTATION:
!the comlexity of this algorithm is N^4, but it can be N^3
!it will be better to rewrite it without cronecker deltas


subroutine getTersoffforce(ratom,ftot,de,distance)
	use classicMD, only: PotentialParam,potential,distanceNeigh, atom
	use neighbor_map, only: neigh_classic,neighn_classic,neigh_b_classic
	use interactions, only: imass
	use configuration, only: natoms,nspecies,xl
	implicit none
	interface
		function tersoff_fcutoff(ro,R,D)
			real, intent(in) :: ro,R,D
			real tersoff_fcutoff
		end function tersoff_fcutoff

		function tersoff_dfcutoff(ro,R,D)
			real, intent(in) :: ro,R,D
			real tersoff_dfcutoff
		end function tersoff_dfcutoff

		function delta(i,j)
			integer, intent(in) :: i,j
			integer :: delta
		end function
	end interface	

	real, dimension(3,natoms),intent(in) :: ratom
	real, intent(inout) :: distance(natoms,natoms)
	real, intent(out) :: ftot(3,natoms), de

	integer :: i,j,l,m,dk,nj,nl
	real :: cosTheta,A,B,lambda1,lambda2,lambda3,beta,n,c,d,h,R,rijd,rild,S,ksi
	real :: dksi_ij(3),rij(3),ril(3),dgijl(3),dexpTerm_ij(3),dbij(3),df(3),dcosTheta(3)
	real :: ksi_ij,gijl,expTerm_ijl,expTerm_ij,bij
	type(PotentialParam), pointer :: pot


	if(potential(1,1)%nparam/=12)then
		write(*,*)'wrong number of parameters for thersoff potential. '
		write(*,*)'I need 12 parameters: A,B,lambda1,lambda2,lambda3,beta,n,c,d,h,R,D,ksi.'
		write(*,*)'you gave me only: ',potential%nparam,' parameters'
		write(*,*)'exiting...'
		stop
	else if(potential(1,1)%type(1:7)/='Tersoff')then
		write(*,*)'--- BAD POTENTIAL!',potential(1,1)%type(1:7),potential(1,1)%pType,' ---'
		stop
	endif

	de = 0.0
	do dk = 1,natoms 
!we derive according dk
		df(:) = (/0.0,0.0,0.0/)
		do i = 1,natoms	
			do nj = 1,neighn_classic(i)
				j = neigh_classic(nj,i)
				if( j/=i .or. sum(xl(:, neigh_b_classic(nj, i))*xl(:, neigh_b_classic(nj, i))) > 0 )then
					pot => potential(imass(j),imass(i))
					A = 	pot%params(1)
					B = 	pot%params(2)
					lambda1 = pot%params(3)
					lambda2 = pot%params(4)
					lambda3 = pot%params(5)
					beta = 	pot%params(6)
					n = 	pot%params(7)
					c = 	pot%params(8)
					d = 	pot%params(9)
					h = 	pot%params(10)
					R = 	pot%params(11)
					ksi = 	pot%params(12)
					S = 	pot%cutoff

					ksi_ij = 0
					dksi_ij(:) = (/0.0,0.0,0.0/)
					expTerm_ij = 0

					dexpTerm_ij(:) = (/0.0,0.0,0.0/)
					rij(:) = (-ratom(:,i)+atom(i,nj,ratom,natoms))
					rijd = sqrt(rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3))

					do nl = 1,neighn_classic(i)
						l  = neigh_classic(nl,i)
						if((l /= i .or. sum(xl(:, neigh_b_classic(nl, i))*xl(:, neigh_b_classic(nl, i))) > 0) &
    	                      .and. (l /= j .or. sum(xl(:, neigh_b_classic(nl, i))*xl(:,neigh_b_classic(nl, i)))>0 ) )then
							rild = distanceNeigh(i,nl,ratom,natoms)
							ril(:) = (atom(i,nl,ratom,natoms)-ratom(:,i))

							cosTheta = sum(rij*ril)/(rijd*rild)

							dcosTheta = ((rij + ril)/(rijd*rild) - cosTheta*(rij/(rijd*rijd)+ril/(rild*rild)))*delta(i,dk) &
				   						- ( ril/(rijd*rild) - cosTheta*rij/(rijd*rijd))*delta(j,dk) &
				   						- ( rij/(rijd*rild) - cosTheta*ril/(rild*rild))*delta(l,dk)

							gijl = (1+c*c/(d*d)-(c*c)/(d*d+(h-cosTheta)**2))
							dgijl = -((c*c)/((d*d+(h-cosTheta)**2)**2))*2*(h-cosTheta)*dcosTheta
!this term is OK	
							expTerm_ijl = exp((lambda3**3)*(rijd-rild)**3)
			 				expTerm_ij = expTerm_ij + expTerm_ijl
							dexpTerm_ij = dexpTerm_ij + expTerm_ijl*((lambda3**3)*3*(rijd-rild)**2)*( &
								(rij/rijd - ril/rild)*delta(i,dk) & 
								- rij/rijd*delta(j,dk) + ril/rild*delta(l,dk) )

							ksi_ij = ksi_ij + tersoff_fcutoff(rijd,R,S)*expTerm_ij*gijl

							dksi_ij(:) = dksi_ij(:) + tersoff_dfcutoff(rijd,R,S)* &
								expTerm_ij*gijl*(delta(i,dk)-delta(j,dk))*rij/rijd + &
								tersoff_fcutoff(rijd,R,S)*expTerm_ij*dgijl + &
								tersoff_fcutoff(rijd,R,S)*dexpTerm_ij*gijl
						endif
					enddo

					bij = ksi*(1+(beta*ksi_ij)**n)**(-1/(2*n))
					if(ksi_ij /= 0)then
						dbij(:) = -0.5*((1+(beta*ksi_ij)**n)**(-1/(2*n)-1))*(beta**n)*(ksi_ij**(n-1))*dksi_ij(:)
					else if(ksi_ij == 0)then
    					dbij(:) = (/0.0,0.0,0.0/)
					else
						write(*,*)'ERROR in getforces_classic_tersoff.f90: ksi(i,j)=',ksi_ij,' is smaller than zero!'
						write(*,*)'stopping...'
						stop
					endif

					if(dk == 1)&
						de = de + tersoff_fcutoff(rijd,R,S)*(A*exp(-lambda1*rijd) - B*bij*exp(-lambda2*rijd))/2
						df(:) = df(:) + 0.5*( ( delta(i,dk)-delta(j,dk) )*rij(:)/rijd )*( &
						tersoff_dfcutoff(rijd,R,S)*( A*exp(-lambda1*rijd) - B*exp(-lambda2*rijd)*bij ) + &
						tersoff_fcutoff(rijd,R,S)*( -lambda1*A*exp(-lambda1*rijd)  + B*bij*lambda2*exp(-lambda2*rijd) ) &
						) - 0.5*tersoff_fcutoff(rijd,R,S)*B*exp(-lambda2*rijd)*dbij(:) ! 1.602 - units, eV/A -> nN
				endif
			enddo
		enddo
		ftot(:,dk) = df(:)
	enddo
end subroutine getTersoffforce

function delta(i,j)
	integer, intent(in) :: i,j
	integer :: delta

	if(i == j)then
		delta = 1
	else
		delta = 0
	endif
end function

function tersoff_fcutoff(ro,R,S)
	implicit none
	real, intent(in) :: ro
	real, intent(in) :: R,S
	real :: tersoff_fcutoff
	real, parameter :: pi=3.1415927
	if(ro<R)then
		tersoff_fcutoff=1
	elseif(ro>S)then
		tersoff_fcutoff=0
	else
		tersoff_fcutoff=0.5*(1+cos(Pi*(ro-R)/(S-R)))
	endif
end function tersoff_fcutoff

function tersoff_dfcutoff(ro,R,S)
	implicit none
	real, intent(in) :: ro
	real, intent(in) :: R,S
	real :: tersoff_dfcutoff
	real, parameter :: pi=3.1415927
	if(ro<R)then
		tersoff_dfcutoff=0	
	else if(ro>S)then
	    tersoff_dfcutoff=0
	else 
		tersoff_dfcutoff=-Pi*0.5*sin(Pi*(ro-R)/(S-R))/(S-R)
	endif
end function tersoff_dfcutoff

