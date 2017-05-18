! copyright info:
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
! reading of the parameters for the classical MD. readdata_classicMD() contains
! the main routine which read file Cdata/usePotentials.in and call the subroutine
! which read the required potential (in Cdata/POTENTIAL_NAME.dat)
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

subroutine readdata_classicMD()
	use configuration, only: nspecies,symbolA
	use classicMD, only: POTENTIALPARAM, potential, getLJforce, getRGLforce,freq_of_outputs !, getClassForcePt

	implicit none
	integer i,j
	character(5) :: elements
	character(2) :: element1,element2
	integer :: eof,pot=0
	character(256) :: line
	character(10) :: num
	logical exists, allIsTheSame,freqOfOutputsSet

	interface
		subroutine readPotential(element1,element2,potentialId)
			character(2), intent(in) :: element1, element2
			integer,intent(in) :: potentialId
		end subroutine

		subroutine checkPotsFilled()
		end subroutine
  		
		subroutine fillVdWoffcrosselements()
		end subroutine
	end interface


	allocate(potential(nspecies,nspecies))
	do i = 1, nspecies
		do j = 1, nspecies
			potential(i,j)%ptype = 0
		enddo
	enddo

	exists = .false.
	allIsTheSame = .false.
	freq_of_outputs = 1000  !implicite value
	inquire( file = 'Cdata/usePotential.in',exist = exists)
	write(*,*)'======================================================================='

	if ( exists )then
		open( 10, file = 'Cdata/usePotential.in', action = 'read')
		do
			freqOfOutputsSet=.false.
			read( 10 , '(a)' ,iostat = eof)line
			if ( eof < 0 ) exit
			i = index(line,':')+1
			elements = trim(line( i : index(line,' ') ))
			if( index(elements,'-') > 0 ) then      !the interaction between two elements
				element1 = elements(1:index(elements,'-')-1)
				element2 = elements(index(elements,'-')+1 :)
			else
				if( trim(line( 1 : index(line,' ') )) == 'freq_of_outputs' )then
					num = trim(line( index(line,' ')+1 :))
					read(num,'(i10)')freq_of_outputs
					freqOfOutputsSet=.true.
				elseif( trim(elements) == 'all' )then
					allIsTheSame = .true.
!					write(*,*)'all is the same ',line
				else
					element1 = elements( 1 :  )
					element2 = element1
				endif
			endif
			if(.not. freqOfOutputsSet)then
				do j = 1, nspecies
					if( allIsTheSame )element2 = trim(symbolA(j))
					do i = 1, nspecies
						if( allIsTheSame )element1 = trim(symbolA(i))
						if( trim(symbolA(i)) == trim(element1) .and.  trim(symbolA(j)) == trim(element2) )then
							select 	case(trim(line( index(line,' ')+1 :)))
								case('Lennard-Jones')
									pot = 1
								case('RGL')
									pot = 2
								case('Tersoff')
									pot = 3
								case('VdW')
									pot = 4
								case default
									write(*,*) 'can not identify the potential ',line( index(line,' ')+1 :),'in usePotential.in'
									stop
							end select
							call readPotential(element1,element2,pot)
						endif
					enddo
				enddo
				line=''
			endif
		enddo
!		if(pot > 0)then
!			select case(pot)
!			case (1)
!				getClassForcePt => getLJforce
!			case (2)
!				getClassForcePt => getRGLforce
!			case default
!				write(*,*) 'can not identify the potential ',line( index(line,' ')+1 :),'in usePotential.in'
!				stop
!			end select
!		endif
		close( 10 )
	else
		write(*,*)'=========================================================='
		write(*,*) 'File Cdata/usePotential.in does not exist. You must specify which potential'
		write(*,*)' you want to use in this file. Example for system using Lennard-Jones potential:'
		write(*,*) 'all Lennard-Jones'
		write(*,*)
		write(*,*) 'List of known potentials: Lennard-Jones,RGL'
		stop
	endif
  	if(pot==4) call fillVdWoffcrosselements()
	call checkPotsFilled()
	write(*,'(a,i5,a)')'frequency of outputs: each ',freq_of_outputs,' cycles'
end subroutine readdata_classicMD

subroutine readPotential(element1,element2,potentialId)
	use configuration, only: nspecies, symbolA
	implicit none
	character(2), intent(in) :: element1, element2
	integer,intent(in) :: potentialId
	logical :: exists
	character(64) :: filename
	character(20) :: msg
	integer :: npot
	external :: fillVdWparams
	logical, save :: vdwreaded = .false.

	interface
		subroutine readPotParam(element1,element2,filename,nparamIn,msg,potId)
			character(2), intent(in) :: element1,element2
			character(64), intent(in) :: filename
			integer, intent(in) :: nparamIn
			character(20), intent(in) :: msg
			integer, intent(in) :: potId
		end subroutine
	end interface
!	write(*,*)'====================',element1,element2,'========================='
	exists = .false.
	select 	case(potentialId)
		case(1)
			filename='Cdata/Lennard-Jones.dat'
			msg = 'LJ'
			npot=3
			inquire( file = filename ,exist = exists)
			if ( .not. exists )then
				write(*,*)'can not find file ',filename
				stop
			endif
		case(2)
			filename = 'Cdata/RGL.dat'
			msg='RGL'
			npot=8
			inquire( file = filename,exist = exists)
			if ( .not. exists )then
				write(*,*)'can not find file ',filename
				stop
			endif
		case(3)
			filename = 'Cdata/Tersoff.dat'
			msg='Tersoff'
			npot=13
			inquire( file = filename,exist = exists)
			if ( .not. exists )then
				write(*,*)'can not find file ',filename
				stop
			endif
 		case(4)
        		filename = 'Cdata/VdW.dat'
        		msg = 'VdW'
        		npot = 4
        		inquire( file = filename, exist = exists)
        		if (.not.exists)then
            			filename='vdw.optional'
            			inquire( file = filename, exist = exists)
            			if (.not.exists)then
                			write(*, *) 'can not find file ', filename
                			stop
            			endif
            			if(.not. vdwreaded)call readvdw(nspecies, symbolA, 1)
            			vdwreaded = .true.
        		endif
		case default
			write(*,*) 'can not identify the potential id=',potentialId
			stop
		end select
    	if(potentialId==4)then
        	call fillVdWparams(element1,element2,filename,npot, msg, potentialId)
    	else
		call readPotParam(element1,element2,filename,npot,msg,potentialId)
	endif
end subroutine

subroutine fillVdWparams(element1,element2,filename,nparam, msg, potId)
    use configuration, only: nspecies, symbolA
    use classicMD, only: POTENTIALPARAM, potential
    use interactions, only: C6, p_alpha, R0
    use neighbor_map, only: range_vdw

    implicit none
    interface
        subroutine readPotParam(element1, element2, filename, nparamIn, msg, potId)
            character(2),  intent(in) :: element1, element2
            character(64), intent(in) :: filename
            integer,       intent(in) :: nparamIn
            character(20), intent(in) :: msg
            integer,       intent(in) :: potId
        end subroutine
    end interface
    
    character(2), intent(in)  :: element1, element2
    character(64), intent(in) :: filename
    character(20), intent(in) :: msg
    integer, intent(in) :: potId, nparam
    integer :: i, j
    
    if(filename=='vdw.optional')then
        do i = 1, nspecies
		if(.not. allocated(potential(i, i) % params))then          
			potential(i,i) % cutoff = range_vdw
          		potential(i, i) % type = msg(1:10)
          		potential(i, i) % ptype = potId
          		potential(i, i) % nparam = nparam
          		allocate(potential(i, i) % params(potential(i, i) % nparam))
          		potential(i, i) % params(1) = C6(i)
          		potential(i, i) % params(2) = p_alpha(i)
          		potential(i, i) % params(3) = R0(i)
		endif
        enddo
    else if(element1==element2)then
        call readPotParam(element1, element2, filename, nparam, msg, potId)
    endif
end subroutine


subroutine fillVdWoffcrossElements()
    use configuration, only: nspecies
    use classicMD, only: POTENTIALPARAM, potential
    implicit none
    integer :: i,j
    do i = 1, nspecies
        do j = 1, nspecies
            if(i<j)then
!these potentials are used for calculating of neighbors only
                potential(i,j) % cutoff = (potential(i,i) % cutoff + potential(j,j) % cutoff)/2
                potential(j,i) % cutoff = potential(i,j) % cutoff
		potential(i, j) % ptype = 4
                potential(j, i) % ptype = 4
            endif
        enddo
    enddo
end subroutine

subroutine readPotParam(element1,element2,filename,nparamIn,msg,potId)
	use configuration, only: nspecies,symbolA
	use classicMD, only: POTENTIALPARAM, potential

	implicit none
	character(2), intent(in) :: element1,element2
	character(64), intent(in) :: filename
	integer, intent(in) :: nparamIn
	character(20), intent(in) :: msg
	integer, intent(in) :: potId
	integer ::  i,j
	integer :: nparam, eof
	real :: param(nparamIn)
	character(256) :: line
	character(10) :: forma,forma2

	interface
		subroutine fillTersofCrossPar()
		end subroutine
	end interface

	nparam = nparamIn-1
	open(unit = 1, file=filename)
!	write(*,*)'element1=',trim(element1)
	write(forma,'(a,a,a)')trim(element1),'-',trim(element2)
	write(forma2,'(a,a,a)')trim(element2),'-',trim(element1)
	do
		read( 1,'(a)', iostat = eof) line
		if((line(1:1) /= '#' .and. (index(line,trim(forma)) > 0 .or. index(line,trim(forma2)) > 0 )) .or. eof < 0 ) exit
	enddo
	close(1)
	if( eof < 0 )then
		write(*,*)'Can not find interaction for ',trim(forma),' in file ',filename
		stop
	endif
	forma=''
	line = line(index(line,'[')+1:index(line,']')-1)
	do i = 1, nparamIn
		j = index(line,',',.true.)
		if( j > 0 )then
			read(line(j+1:),*)param(nparam+2-i)
			line = line(1: j-1)
		else
			if( nparam+2-i /=1 )then
				write(*,*)
				write(*,'(a,i2,11a)') 'I need ',nparamIn,' parameters in format: "',trim(element1),'-',trim(element2), &
						' param=[sigma, ro, cutoff] . Please check file ',filename, &
					' if the potential for ',trim(element1),'-',trim(element2),' is setted correctly'
				write(*,'(a,i2,4a)')'You have only ',i,' parameters for ',trim(element1),'-',trim(element2)
				write(*,*)
				stop
			endif
			read(line(1:),*)param(1)
		endif
	enddo
	if( j > 0 ) then
		write(*,*)
		write(*,'(6a)')'Too much parameters in row with interaction for ',trim(element1),'-',trim(element2),' in file',filename
		write(*,'(a,i2,6a)')'I need ',nparamIn,' parameters exactly in format: "',trim(element1),&
							'-',trim(element2),' param=[params_sepatated_by_comas,...] . Please check file ',filename
		stop
	endif
!	write(*,*)'params=',param
	do i = 1, nspecies
		do j = 1, nspecies
			if ( trim(symbolA(i)) == trim(element1) .and. trim(symbolA(j)) == trim(element2) .and. potential(i,j)%ptype == 0 ) then
				potential(i,j)%type = msg
				potential(i,j)%ptype = potId
				potential(i,j)%nparam = nparam
				allocate(potential(i,j)%params( potential(i,j)%nparam ))
				potential(i,j)%params = param(1:nparam)
				potential(i,j)%cutoff = param(nparamIn)
				if( i /= j .and. potential(j,i)%ptype == 0 )then
					potential(j,i)%type = msg
					potential(j,i)%ptype = potId
					potential(j,i)%nparam = nparam
					allocate(potential(j,i)%params( potential(j,i)%nparam ))
					potential(j,i)%params = param(1:nparam)
					potential(j,i)%cutoff = param(nparamIn)
				endif
				exit
			endif
		enddo
	enddo

	if( potId == 3)then
		call fillTersofCrossPar()
	endif
end subroutine

subroutine fillTersofCrossPar()
	use configuration, only: nspecies,symbolA
	use classicMD, only: POTENTIALPARAM, potential

	do i = 1, nspecies
		do j = 1, nspecies
			if( i /= j .and. potential(j,i)%ptype == 0 )then
				potential(i,j)%params(1) = dsqrt( potential(i,i)%params(1)*potential(j,j)%params(1) )
				potential(i,j)%params(2) = dsqrt( potential(i,i)%params(2)*potential(j,j)%params(2) )
				potential(i,j)%params(11) = dsqrt( potential(i,i)%params(11)*potential(j,j)%params(11) )
				potential(i,j)%params(3) = ( potential(i,i)%params(3)+potential(j,j)%params(3) )/2.0
				potential(i,j)%params(4) = ( potential(i,i)%params(4)+potential(j,j)%params(4) )/2.0
				potential(i,j)%cutoff = dsqrt(potential(i,i)%cutoff*potential(j,j)%cutoff)

				potential(j,i)%params(1:nparam) = potential(i,j)%params(1:nparam)
				potential(j,i)%cutoff = potential(i,j)%cutoff
			endif
		enddo
	enddo
end subroutine

subroutine checkPotsFilled()
!this potential check that all possible interactions are defined - it check that potential->typ has defined type ("RGL,"LJ"...etc)
	use configuration, only: nspecies, symbolA
	use classicMD, only: potential

	implicit none
	integer i,j
	character(20) :: forma

	write(*,*)'Checking potentials'

	do i = 1, nspecies
		do j = 1, nspecies
			if( potential(i,j)%ptype == 0 ) then
				write(*,*)'The potential for ',symbolA(i),'-',symbolA(j),' is not defined in file Cdata/usePotentials.in'
				write(*,*)'Please determine the interaction for thise element(s) by row: element1-element2 paramtype'
				stop
			else
				if(potential(i,j)%ptype /= 4 .or. i==j)then
					write(forma,'(a,i2,a)')'(7a,f7.3,a,',potential(i,j)%nparam,'f7.3)'
					write(*,forma)'interaction for ',symbolA(i),'-',symbolA(j),' : ',trim(potential(i,j)%type), &
								', cutoff=',potential(i,j)%cutoff,', params=',potential(i,j)%params
				endif
			endif
		enddo
	enddo

	write(*,*)'Potentials checked'
end subroutine
