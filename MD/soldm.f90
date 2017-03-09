! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

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


! soldm.f90
! Program Description
! ===========================================================================
!       Computes one column of the force constant matrix 
! phi(ix,iatom,jx,jatom) by displacing atom jatom in the jx direction by a 
! displacement u0disp. The row is stored in phirow(ix,iatom) and is written to 
! disk file filephi.
!       Conversion of units: frequency of vibratinal modes is in [Hz] and [cm-1]
! Units of Hessian matrix: [(eV/Ang)/Ang]
! Units of Hessian matrix divided by mass: 1/(1.660*10^-27)[(eV/((Ang^2)*kg)]= 
! (1.602*10^-19)/(1.660*10^-27)[J/(Ang^2*kg)]=
! =(1.602*10^-19)/(1.660*10^-27*10^20)[((kg*m^2)/s^2)/(m^2*kg)]=
! =(1.602/1.66)*10^28[1/s^2]
! frequency of vibrational modes f [Hz], (lambda-eigenvalue of Hessian):
! f=(sqrt(lambda))/2*pi => sqrt((1.602/1.66)*10^28)[1/s] or [Hz]
! frequency of vibrational modes f [cm-1], (lambda-eigenvalue of Hessian):
! f=(sqrt(lambda))/(2*pi*c) => sqrt((1.602/1.66))*(1/299.792458)*10^6 [1/cm]
! from 1/cm to eV :  1.239 81 x 10-4	
! ===========================================================================
! Code written by:
! P. Jelinek & V. Zobac
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine soldm ( iephc )

   use dimensions
   use dynamo
   use configuration
   use interactions

   implicit none
 
! Argument Declaration and Description
! ===========================================================================

   integer,  intent (in) :: iephc

! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter :: mineig = 1.0d-6
!   real, parameter :: c = 29979245800 
!   real, parameter :: c2Hz = sqrt(1.602d0/1.66d0)*(1.0d14/(2.d0*3.1415927d0))
   real, parameter :: c2Hz = 15634980253554.3d0  
!   real, parameter :: c2cm1 = sqrt(1.602d0/1.66d0)*(1.0d0/29979245800d0)*(1.0d14/(2.d0*3.1415927d0))
   real, parameter :: c2cm1 = 521.526804171780d0  
   real, parameter :: cunit = 3.77127d-3
   real, parameter :: cm12eV = 1.23981d-4 
! Local Variable Declaration and Description
! ===========================================================================

   integer ndim
   integer info
   integer iatom, iiatom
   integer jatom, jjatom
   integer ix
   integer jx
   integer i
   integer j
   integer index
   integer lwork
   integer jj
   integer in1

   real norm

   real, dimension (:), allocatable ::  emode
   real, dimension (:), allocatable ::  rwork
   real, dimension (:,:), allocatable ::  evec
   real, dimension (:), allocatable :: sumrow

! e-ph coupling
   integer imode
   integer istate
   integer indexm
   
   real dem 
   real dep 
   real dEptot 
   real dEmtot 
   real E0
   real xmode
   real amp0
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
   write (*,*) ' Entering subroutine soldm '

   ndim = natoms_dm*ndx
   lwork = 100*ndim + 3*ndim*ndim  

! allocate eigenmode vector
   allocate ( emode(ndim) )
   allocate ( rwork(lwork) )
   allocate ( evec(ndim,ndim) )
   allocate (sumrow(ndim)) 

! write dynamical matrix to file
   open(unit=81,file='phidm.dat')
   do i=1,ndim
       write(81,400) (phidm(i,ix),ix=1,ndim)  
    enddo
   close (81)

! check if the sum of rows of Dynamical matrix are approximately zero
   sumrow=0
   do i=1,ndim
     do ix=1,ndim 
       sumrow(i)=sumrow(i) + phidm(i,ix) 
     enddo   
   enddo 
  
   do i=1,ndim
      if ( abs(sumrow(i)) .gt. 1e-6 ) then
         write(*,*)    '**WARNING**'            
         write(*,*) 'Absolute value of the sum of',i,'-th row dynamical matrix is' ,abs(sumrow(i))
         write(*,*) '  '  
      endif
   enddo

! calc matrix: 1/(M_i*M_j)*phi
   i = 0
   do iatom = 1,natoms_dm
      iiatom = jatoms_dm(iatom) 
      do ix = 1,ndx
         i = i + 1
         j = 0
         do jatom = 1,natoms_dm
            jjatom = jatoms_dm(jatom)
            do jx = 1,ndx
               j = j + 1
               phidm(i,j) = phidm(i,j)/(sqrt((xmass(jjatom)*xmass(iiatom))))
           enddo
         enddo
      end do
   end do 

!  The computed eigenvectors are normalized to have Euclidean norm equal
!  to 1 and largest component real.
!
!  VARIABLES:
!  EIGN   ... eigenvalues
!  BRA  ... u(j) left eigenvectors [u(j) = BRA(:,j),
!                  the j-th column of BRA]
!  KET  ... v(j) right eigenvectors [v(j) = KET(:,j),
!                  the j-th column of KET]


   write (*,*) '  '
   write (*,100) 
   write (*,*) '  '
   write (*,*) '   Diagonalize dynamical matrix ...'
   call dsyev ('V', 'U', ndim, phidm, ndim, emode, rwork, lwork, info )
   
   if (info .eq. 0) then
! ket |u>
      evec = Transpose (phidm)
! bra <u|
!      bra = Transpose(curm)
!      bra = Conjg(bra) 
    else
      write (*,*) '   ----  Error in dsyev subroutine  ----'
      write (*,*) '              info = ',info
      write (*,*) '  ---- the code will be terminated  ----   '
      stop
   endif

! write eigenmodes
   write (*,*) '  '
   write (*,*) '  Write eigenmodes of Dynamical matrix'
   index = 1
   do iatom = 1, natoms_dm
      do ix = 1, ndx
        if (emode(index) .lt. 0.0d0) then 
           write (*,*) ' ** WARNING **' 
           write (*,450) index, emode(index) 
           write (*,*) ' ' 
        endif
        if (emode(index) .lt. mineig) emode(index) = 0.0d0
! calc norm of eigenvector
        norm = 0.0d0
        do i = 1,ndim
           norm = norm + evec(index,i)**2.0d0
        enddo
        write (*,350) index, (sqrt(emode(index)))*c2Hz, norm
        write (*,360) index, (sqrt(emode(index)))*c2cm1
        write (*,300) ( evec(index,i), i=1,ndim ) 
        index = index + 1
      enddo
   enddo

  open (unit = 45, file = 'modes.xyz', status = 'unknown')
  do i=1,ndim 
        write(45,'(i5)') natoms      
        write(45,360) i, (sqrt(emode(i)))*c2cm1
    do iatom = 1, natoms 
        in1 = imass(iatom)
      if (ndx .eq. 1 ) then
         write (45,701) symbol(iatom), ratom0(:,iatom), (evec(i,jj),jj=iatom*ndx-(ndx-1),iatom*ndx), 0, 0     
      end if 
      if (ndx .eq. 2 ) then
         write (45,701) symbol(iatom), ratom0(:,iatom), (evec(i,jj),jj=iatom*ndx-(ndx-1),iatom*ndx),0
      end if
      if (ndx .eq. 3 ) then  
         write (45,701) symbol(iatom), ratom0(:,iatom), (evec(i,jj),jj=iatom*ndx-(ndx-1),iatom*ndx)   
      end if  
    end do
  end do 
  close (unit = 45)

! calculate e-ph coupling for selected modes
  if (iephc .eq. 1) then
    open (file = fileephc, unit=17, status='unknown')
    write (*,*) '  ====  calculationg e-ph coupling ===='
    do imode = 1, nvmodes
      indexm = idvmode(imode) 
      norm = 0.0d0
      do ix = 1, ndx*natoms
       norm = norm + evec(indexm,ix)*evec(indexm,ix)
      enddo 
      xmode = (sqrt(emode(indexm)))*c2cm1*cm12ev
      write (*,*) ' doing vib. mode = ',indexm, xmode
      amp0 = cunit*sqrt(temp_ephc/(mass_ephc*xmode))
      write (*,*) ' Amplitude = ', amp0
      write (17,600) indexm, xmode, amp0
      do istate = 1, nephc
        write (*,*) '   doing electronic state = ', eiglist(istate)	
        E0 = eigref(istate)
	dEptot = 0.0d0
	dEmtot = 0.0d0
        do ix = 1,ndx*natoms
         dem = (E0-deigen(istate,ix))/u0disp
!         write(*,*) 'dem:',E0,deigen(istate,ix), (E0-deigen(istate,ix)), u0disp,dem
!         write(*,*) 'vmode: ',evec(indexm,ix)
!         write(*,*) 'dE*u: ',dem*evec(indexm,ix)
         dEmtot = dEmtot + dem*evec(indexm,ix)
	end do ! ix
	write (*,460) abs(dEmtot*amp0)
! write to output file
      write (17,610) istate, E0, abs(dEmtot*amp0)   
      end do ! istate
    end do ! imode
    close (17)
  endif ! iephc

! Deallocate Arrays
! ===========================================================================

   deallocate ( emode )
   deallocate ( evec )
   deallocate ( rwork )
 
! Format Statements
! ===========================================================================
100     format ('===========================================================')
200     format ('  atom = ',i4,' ix = ', i2,'  eigenmode = ',f16.8 )
300     format ( 3f16.8 )
350     format ( '  Mode :',i4,' frequency  = ',E16.6,' [Hz] with norm :'f16.6 )
360     format ( '  Mode :',i4,' frequency  = ',f18.10,' [cm-1]')
450     format ( '  Mode :',i4,' has a negative value  = ',f16.6 )
460     format ( '    e-ph coupling = ',f16.6,' [eV] ' )
401     format (<ndim>f8.3)
400     format ( 9f8.3 )
600     format ( '  == Vib. mode no.',i6,' with frequency  = ',f16.6,' [eV] amplitude = ',f16.6,'[Ang]' )
610     format ( '     El. level no.',i6,' with energy  = ',f16.6,' [eV] e-ph couplong = ',f16.6,'[eV]' )
701     format (2x, a2,6(2x,f12.6))  
   return
 end subroutine soldm
