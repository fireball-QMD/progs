! copyright info:
!
!                             @Copyright 2009
!                FAST (Fireball Atomic Simulation Techniques)
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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


! delta_t_ks.f90
! Program Description
! ===========================================================================
! Calculate non-adiabatic coupling(d_{jk}  contribution V.d_{jk}
! using Kohn-Sham states at different
! time steps, and uses it for the MDET calculation
!
! ===========================================================================
! Code written by Jose Ortega Mateo & Enrique Abad Gonzalez
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine delta_t_ks (itime_step)

        use configuration
        use nonadiabatic
        use density
        use interactions
        use kpoints
        use MD

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer itime_step


! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: hbar = 0.65822d0

! Local Variable Declaration and Description
! ===========================================================================
        integer it
        integer imu, inu, jmu,jnu, iorbital
        integer in1, in2, in3
        integer iatom, jatom
        integer isorp, interaction
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
        real y, rcutoff_i, rcutoff_j, range
        real diff
        real delta
        real, dimension (3) :: r1, r2, r21, sighat
        real, dimension (3,3,3) :: deps
        real, dimension (3,3) :: eps
        real, dimension (numorb_max, numorb_max) :: sx
        real, dimension (3,numorb_max, numorb_max) :: spx
        real, dimension (norbitals, norbitals) :: s
        real, dimension (nele, nele) :: suma
        real, dimension (norbitals,nkpoints) :: deig
        real, dimension (norbitals,nkpoints) :: eig
        real, dimension (3,natoms) :: dvatom
        real, dimension (3,natoms) :: v
        real, dimension (3,natoms,nele,nele) :: dgks
        real, dimension (3,natoms,nele,nele) :: g
        complex, dimension (norbitals, norbitals, nkpoints) :: c_old
        complex, dimension (norbitals, norbitals, nkpoints) :: c_new
        complex, dimension (norbitals, norbitals, nkpoints) :: dc_na
        complex, dimension (norbitals, norbitals, nkpoints) :: caux
        complex aim
        complex a0
        complex a1

! ===========================================================================
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

! ===========================================================================
! Calculate non-adiabatic couplings using Kohn-Sham states at different
! time steps
! (1) Define overlap matrix between orbitals mu at time "t" and orbitals
! nu at time "t+dt"
        do iatom = 1, natoms
         r1(:) = ratom_old(:,iatom)
	 rcutoff_i = 0.0d0
         in1 = imass(iatom)
         do imu = 1, nssh(in1)
           if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
         end do
         do jatom = 1, natoms
          r2(:) = ratom(:,jatom)
          in2 = imass(jatom)
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
! check the cutoff radii
          rcutoff_j = 0.0d0
          in2 = imass(jatom)
          do imu = 1, nssh(in2)
              if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
          end do
          range = (rcutoff_i + rcutoff_j - 0.01d0)**2
	  range = sqrt(range)
!print *, 'llegamos al if'
          if (y .gt. range) then
!print *, 'y > range',y,range
           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
	     s(jmu,jnu) = 0.0d0
	    end do
	   end do
	  else
!print *, 'y < range',y,range
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
           else
            sighat(:) = r21(:)/y
           end if
           call epsilon (r2, sighat, eps)
           call deps2cent (r1, r2, eps, deps)
           isorp = 0
           interaction = 1
           in3 = in2
! iforce=0
           call doscentros (interaction, isorp, 0, in1, in2, in3, y, &
      &                     eps, deps, sx, spx)
           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
             s(jmu,jnu) = sx(imu,inu)
            end do
           end do
	  end if
         end do
        end do

! ===========================================================================
! Calculate overlap between Kohn-Sham states at different
! time steps
! Non-adiabatic term: dot pruduct sum
         sumb = 0.0d0
        do ikpoint = 1, nkpoints
         do ij = 1, nele
          do ik = 1, nele
           do imu = 1, norbitals
            do inu = 1, norbitals
            sumb(ik,ij) = sumb(ik,ij) +                                 & 
     &      bbnkre_old(imu,map_ks(ik),ikpoint)*bbnkre(inu,map_ks(ij),ikpoint)*s(imu,inu)
            end do
           end do
          end do
         end do
        end do
       
! Calculate non-adiabatic dot-product sum using Kohn-Sham states at different
! time steps
         do ij = 1, nele
          do ik = 1, nele
          dnac (ik,ij) = (sumb(ik,ij) - sumb(ij,ik))/(2.0d0*dt)
!	  write (552,*) ik, ij, dnac(ik,ij), sumb(ik,ij)/dt, sumb(ij,ik)/dt
!	  write(*,301) ij, ik, sumb(ij,ik)
          end do
         end do

! Calculate non-adiabatic dot-product sum using non-adiabatic couplings
! gks
! ===========================================================================
!       ddt = dt / nddt
!       deig = eigen_k - eigen_old 
!       dvatom = vatom - vatom_old
!       dgks = gks - gks_old
! ===========================================================================
!       delta = 0.5d0
! Interpolation 
!        v = vatom_old + dvatom*delta
! JOM-testv
         v = (ratom - ratom_old)/dt
!        g = gks_old + dgks*delta
! JOM-test
!        g = gks_old
!        v = vatom_old
         g = gks
!        v = vatom
!        eig = eigen_old + deig*delta
! Non-adiabatic term: dot pruduct sum
         suma = 0.0d0
         do ik = 1, nele
          do ij = 1, nele
           do iatom = 1, natoms
            do ix = 1, 3
              suma(ik,ij) = suma(ik,ij) + v(ix,iatom)*g(ix,iatom,ik,ij)	
            end do
           end do
          end do
         end do
! Compare both ways to calculate non-adaiabatic contribution
           do ik = 1, nele
            do ij = 1, nele
             diff = suma(ik,ij)-dnac(ik,ij)
             write(215,300)ik,ij,suma(ik,ij),dnac(ik,ij)
            end do
           end do
!        write(*,300)3,2,suma(3,2),dnac(3,2)
! Write c_na
!         do ia = 1, norbitals
!          do ik = 1, norbitals
!          write(*,100)ia,ik,cabs(c_na(ia,ik,1)),c_na(ia,ik,1)
!         end do
!        end do
! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================
100     format ('c_na',2i4,f8.4,2(f7.3))
!300     format ('NAC-SUMS',2i4,2f8.4)
300     format (2i4,2f14.8)
301     format ('S(t,tprime)',2i4,1f8.4)
400     format ('S',4f7.3)


        return
        end subroutine delta_t_ks

