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


! fewest_switches.f90
! Program Description
! ===========================================================================
!       This routine determines the hoppings between
!       Kohn-Sham states
!
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine fewest_switches (itime_step)

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


! Local Variable Declaration and Description
! ===========================================================================
        integer iele
        integer jele
        integer iband
        integer jband
        integer iatom
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
        integer iswitch
        real, dimension (nele, nele) :: d
        real, dimension (nele) :: prob
        real, dimension (3,natoms) :: v
!        real, dimension (3,natoms,nele,nele) :: g
        complex aim
        complex a0
        complex a1
        complex  akj
        real ajj
        real bkj
        real :: xrand
	real, dimension (nele, nele) :: suma



! Procedure
! ===========================================================================
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)
!----------------------------------------------------------
! Initialize seed for random_number
        call random_seed
!----------------------------------------------------------

!JOM-info : I assume that nkpoints = 1
!----------------------------------------------------------
        if (nkpoints .gt. 1) then
         write(*,*)'nkpoints=',nkpoints
         write(*,*)'nkpoints is greater then 1'
         write(*,*)'in subroutine fewest_switches'
         write(*,*)'not ready, must stop'
         stop
        end if
!----------------------------------------------------------
! JOM-info : calculate the dot product V.gks
! v is the velocity and g is the nonadiabatic coupling gks
         suma = 0.0d0
         do ik = 1, nele
          do ij = 1, nele
           do iatom = 1, natoms
            do ix = 1, 3
             suma(ik,ij) = suma(ik,ij) + vatom(ix,iatom)*gks(ix,iatom,ik,ij)
            end do
           end do
          end do
         end do
! but at what time step,  ?
!       v = vatom_old
!       g = gks_old
!       v = vatom
        v = (ratom - ratom_old)/dt
!        do iele = 1, nele
!         do jele = 1, nele
!           do ix = 1, 3
!             do iatom = 1, natoms
!               g(ix,iatom,iele,jele) = gks(ix,iatom,iele,jele)
!             end do
!           end do
!         end do
!        end do
!
!        d = 0.0d0
!        do ik = 1, nele
!         do ij = 1, nele
!          do iatom = 1, natoms
!           do ix = 1, 3
!           d(ik,ij) = d(ik,ij) + v(ix,iatom)*gks(ix,iatom,ik,ij)
!           end do
!!	  write (211,101) ik,ij,iatom,v(1,iatom),v(2,iatom),v(3,iatom),   &
!!      &  ((v(1,iatom)*gks(1,iatom,ik,ij))+(v(2,iatom)*gks(2,iatom,ik,ij))+(v(3,iatom)*gks(3,iatom,ik,ij))) ! ENRIQUE-JOM
!          end do
!         end do
!        end do

!	do iatom =1, natoms
!	  do ik = 1, nele
!	    do ij = 1, nele
!!	      write_suma = (v(1,iatom)*gks(1,iatom,ik,ij))+(v(2,iatom)*gks(2,iatom,ik,ij))+(v(3,iatom)*gks(3,iatom,ik,ij))
!           write(210,102)ik,ij,iatom,gks(:,iatom,ik,ij),   &
!        & sqrt((gks(1,iatom,ik,ij)*gks(1,iatom,ik,ij))+(gks(2,iatom,ik,ij)*gks(2,iatom,ik,ij))+(gks(3,iatom,ik,ij)*gks(3,iatom,ik,ij))), &
!        &  ((v(1,iatom)*gks(1,iatom,ik,ij))+(v(2,iatom)*gks(2,iatom,ik,ij))+(v(3,iatom)*gks(3,iatom,ik,ij))) ! ENRIQUE-JOM
!	    end do
!	  end do
!          write (211,101) iatom,v(1,iatom),v(2,iatom),v(3,iatom),sqrt((v(1,iatom)*v(1,iatom))+(v(2,iatom)*v(2,iatom))+(v(3,iatom)*v(3,iatom)))
!	end do
! ===========================================================================
! Calculate hopping probabilities for fewest switches
! map_ks(iele) gives back the corresponding adiabatic KS state
! we follow the possible transitions associated with states
! iele=1,nele
        do ikpoint = 1, nkpoints
         do ij = 1, nele
!----------------------------------------------------------
! Random numbers for Monte-Carlo
         call random_number(xrand)
!----------------------------------------------------------
! JOM-test
!        read (2121,*) xrand  
         write(*,*)'random',xrand

!----------------------------------------------------------
          ajj = real(conjg(c_na(ij,ij,ikpoint))*c_na(ij,ij,ikpoint))
          do ik = 1, nele
             akj = c_na(ij,ik,ikpoint)*conjg(c_na(ij,ij,ikpoint))
!            bkj = -2.0d0*real(conjg(akj)*dnac(ik,ij))
             bkj = -2.0d0*real(conjg(akj)*suma(ik,ij))
! JOM-warning: may be later we can "imporve" this by using eq(29) in JCP
! 101 4657 (1994)
!----------------------------------------------------------
!JOM-info : probability of the j ---> k transition
           prob(ik) = bkj*dt/ajj
           write(*,*)'prob',ij,ik,prob(ik)
           if (prob(ik) .lt. 0.0d0) then
            prob(ik) = 0.0d0
           end if
!----------------------------------------------------------
! JOM-test
!          if(prob(ik) .gt. 0.0001) then
!          write(*,*)'prob',ij,ik,prob(ik)
!          write(*,*)'akk', real(conjg(c_na(ia,ik,ikpoint))*c_na(ia,ik,ikpoint))
!          write(*,*)'ajj',ajj
!          write(*,*)'akj',akj
!          write(*,*)'bkj,dt',bkj,dt
!          end if
!----------------------------------------------------------
          end do ! do ik = 1, nele
!----------------------------------------------------------
! Monte-Carlo calculation for possible transitions
! JOM-warning : we should also allow transitions to states that are not
! fully occupied (ioccupy_na = 0, 1) [from states that are occupied
! ioccupy_na = 1, 2 ]. Use iocc for this (fix later)
!         iocc (ik) = ioccupy_na (ik, ikpoint)
!----------------------------------------------------------
          call mc_switch (xrand, nele, prob, ij, ikpoint, iswitch)
!----------------------------------------------------------
! JOM-test
!       write(*,*)'nele',nele,iele,ia,ij
!       write(*,*)'xrand',xrand
!       write(*,*)'iswitch',iswitch
!----------------------------------------------------------
          if (iswitch .ne. 0) then
             write(*,*)'SWITCH!!',ij, '--->',iswitch
!----------------------------------------------------------
! perform transition ij ---> iswitch
             call transition (itime_step, ij, iswitch, ikpoint)
!----------------------------------------------------------
             return  ! we can only have one switch
          end if

         end do !nele
        end do !kpoints
! ===========================================================================


  

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format ( 1i4, 4f8.4)
102     format ( 3i4, 5f8.4)


        return
        end subroutine fewest_switches

