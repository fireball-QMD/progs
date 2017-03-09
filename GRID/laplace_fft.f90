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


! laplace_fft.f90
! Program Description
! ===========================================================================
! Subroutine that calculates the Hartree Potential via the poisson equation
! in reciprocal space.
! Uses FFTW library.
! d^2 Vh(r) = -4*Pi*n(r)   ; n(r)=sum_k(n(k)*exp(ikr))
! d^2 Vh(k)*exp(ikr) = i^2*k^2*Vh(k)*exp(ikr)= -1*k^2*Vh(k)*exp(ikr)
! final formulae for k-conponent of Vh:
! Vh(k) = 4*Pi*n(k)/k^2 
!
! ===========================================================================
! Code written by:
! ===========================================================================
! Jan Vlachy
! jan.vlachy@gmail.com
! and 
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
!
! Program Declaration
! ===========================================================================
 SUBROUTINE laplace_fft (icluster, ibias)   

!Necessary for the main calculation
   USE constants_fireball
   USE grid 
   USE configuration
   USE outputs
   USE bias

!Necessary for creating .xsf file
   USE charges
   USE interactions 
   

   IMPLICIT NONE
   INCLUDE 'fftw3.f'  
 
! Argument Declaration and Description
! ===========================================================================
! Input

!Just for .xsf
   INTEGER, INTENT (in) :: icluster 
   INTEGER, INTENT (in) :: ibias

!Output

 
! Local Parameters and Data Declaration
! ===========================================================================

   INTEGER, PARAMETER :: DPC = 8
   REAL  (kind = DPC),PARAMETER  :: EPSILON = 1E-15
   REAL  (kind = DPC),PARAMETER  :: ABOHR3 = ABOHR**3
 
   interface 
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  end interface   
! Local Variable Declaration and Description
! ===========================================================================

   !variable used by FFTW
   INTEGER(DPC) :: plan
   
   real(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: in3D
   COMPLEX(DPC), ALLOCATABLE, DIMENSION(:,:,:) :: four3D
   real, target, allocatable, dimension (:) :: resf
! FFTW is perhaps capable of transforming even 3D data in 1D
! arrays.
! Might be worth trying out later.

   INTEGER :: i,j,k
   INTEGER :: index
  
   
   INTEGER :: kx,ky,kz
   REAL(DPC), DIMENSION(3) :: kvector
   REAL(DPC)  :: kvector2 !kvector components squared
   REAL(DPC)  :: fourpi
   REAL(DPC)  :: rfact

   !For the Visualization part.
   INTEGER i0, j0, k0
   INTEGER iatom
   
   real dmax
   real, dimension (:), pointer   :: pmat
   real z
   real zdown
   real zup
   real vshift
   
   character (len=40) filename
   character (len=30) mssg

! Procedure
! ===========================================================================

!Initialization
! const 4*pi
   fourpi = 4.0d0*pi
!FFTW can save half of the space for the real data.
!That is why four3D uses just rm1 / 2.
   ALLOCATE (in3D(rm1,rm2,rm3))
   ALLOCATE (four3D(rm1/2 + 1,rm2,rm3))
   allocate (resf(0:nrm-1))   
   resf = vcaG
   
!Get data from 1D input to 3D
!drhoG, rm are from the 'grid' module
   DO k = 1, rm3
      DO j = 1, rm2
         DO i = 1, rm1
            index = i + rm1*(j-1) + rm1*rm2*(k-1) - 1
! (to remind the current units are 1/Ang^3)  
            in3D(i,j,k) = drhoG(index)   
         END DO
      END DO
   END DO
 

!Forward FFT  n(r) -> n(k)
!  
!If needed, there are more precise FFTW modes than FFTW_ESTIMATE.
!But so far, the 'ESTIMATE' has proved to be good enough
!for all purposes (even better).
   
   CALL dfftw_plan_dft_r2c_3d(plan,rm1,rm2,rm3,in3D,four3D,FFTW_ESTIMATE)!FFTW_MEASURE)
   CALL dfftw_execute(plan)
   CALL dfftw_destroy_plan(plan)


!Comparing coefficients in Fourier series of source term and
! solution
!
!FFT saves 'positive' frequencies in the first half of the arrays,
!the 'negative' ones in the second half.
    
  
   solve : DO i = 1,(rm1 / 2 + 1)

      IF (i <= (rm1 / 2 + 1)) THEN 
         kx = i-1
      ELSE
        kx = rm1-i+1
      END IF
          

      DO j = 1,rm2
         IF (j <= (rm2 / 2 + 1)) THEN 
           ky = j-1
        ELSE
           ky = rm2-j+1
        END IF
        
        DO k = 1,(rm3)
           IF (k <= (rm3 / 2 + 1)) THEN 
              kz = k-1
           ELSE
              kz = rm3-k+1
           ENDIF
           
            !Variable transform 
            kvector(1) = rlvec(1,1)*kx + rlvec(1,2)*ky + rlvec(1,3)*kz
            kvector(2) = rlvec(2,1)*kx + rlvec(2,2)*ky + rlvec(2,3)*kz
            kvector(3) = rlvec(3,1)*kx + rlvec(3,2)*ky + rlvec(3,3)*kz


            !EPSILON prevents the divide-by-zero blow-up
            !in the (1,1,1) term.
!            kvector = kvector
            kvector2 = kvector(1)*kvector(1) + kvector(2)*kvector(2) +&
                 & kvector(3)*kvector(3) + EPSILON
!            kvector2 = 1.0d0*kx*kx + 1.0d0*ky*ky + 1.0d0*kz*kz + EPSILON
!            write (*,*) kx,ky,kz
!            write (*,*) 'kvector =',kvector(:)
!            write (*,*) 'G2 =',kvector2 
! we divide both part of the complex number by the real number
            rfact = fourpi/kvector2
!            rfact = 1.0d0
!            write (*,*) four3D(i,j,k), rfact
            four3D(i,j,k) =  rfact*four3D(i,j,k)   
!            write (*,*) four3D(i,j,k)     
         END DO
      END DO
   END DO solve


   four3D(1,1,1) = 0.0d0 !zero Fourier term

   
!Backward FFT  Vh(k) -> Vh(r)
   CALL dfftw_plan_dft_c2r_3d(plan,rm1,rm2,rm3,four3D,in3D,FFTW_ESTIMATE)!MEASURE)
   CALL dfftw_execute(plan)
   CALL dfftw_destroy_plan(plan)
   
! convert units of the potential from a.u. to eV
! and divide by nrm as for FFT renormalization
   rfact = eq2/nrm
!   rfact = 1.0d0/nrm
      
!Transform back to 1D, returning potential
   DO k = 1,rm3
       DO j = 1,rm2
          DO i = 1,rm1
             index = i + (j-1)*rm1 + (k-1)*rm1*rm2 - 1
! renormalize 
             vcaG(index) = in3D(i,j,k)*rfact
          END DO
       END DO
    END DO

! Add bias voltage linear profile
    if (ibias .eq. 1) then
     zup = zb1 - g0(3)
     zdown = zb0 - g0(3)
     DO k = 1,rm3
      z = (k-1)*elvec(3,3)
      if (z .lt. zdown) then
       vshift =  VBias*0.5
      else if (z .gt. zup) then
      vshift = -1.0d0*VBias*0.5
      else  !  zb0 < z < zb1
       vshift = (0.5d0 - (z - zdown)/(zup-zdown))*VBias      
      end if !z
      DO j = 1,rm2
       DO i = 1,rm1
         index = (i-1) + (j-1)*rm1 + (k-1)*rm1*rm2 
         vcaG(index) = vcaG(index) + vshift
       END DO !i
      END DO  !j
     END DO   !k
    end if ! if ibias

    DEALLOCATE(in3D)
    DEALLOCATE(four3D)

! evaluate residua of vcaG
    dmax = 0.0d0
    do index = 0, nrm-1
     dmax = max(dmax, abs(vcaG(index)-resf(index)))
    enddo
    write (*,*) ' residua vcaG = ',dmax,dmax*dvol
    
   
! Data for Visualization in XCrysDen
! ===========================================================================
! write out xsf files
    if (iwrtxsf .eq. 1) then 
! write out hartree potential steaming from drho
     pmat => vcaG
     filename = 'fftpot.xsf'
     mssg = 'hartree_3d'
     call writeout_xsf (filename, mssg, pmat) 
    
! write out  hartree potential into vhartree.xsf file
! 9/12/08 added vna to dvnh
!  we used that for postprocessing to calculate workfunction
     resf = vcaG + vnaG     
     pmat => resf
     filename = 'vhartree.xsf'
     mssg = 'hartree_3d'
     call writeout_xsf (filename, mssg, pmat) 
    endif
  
! deallocate 
    deallocate (resf)  
 
! Format Statements
! ===========================================================================
100 format (2x, 'Hartree Energy = ',f14.7,' [eV]')
200  format('index = ',6i7)
300 format (e14.6)


END SUBROUTINE laplace_fft
