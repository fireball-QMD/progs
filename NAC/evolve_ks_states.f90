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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! evolve_ks_states.f90
! Program Description
! ===========================================================================
!       This routine integrates the TD equations for the coefficients
!       c_na of the TD-wfs
!
! c_na (ia,ij,ikpoint): phi (ia) = \Sum_ij c_na(ia,ij)*psi(ij)
! use Runge-Kutta 4th order
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine evolve_ks_states (itime_step)

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
!       integer, parameter :: nddt = 1000
        real, parameter :: hbar = 0.6582119d0

! Local Variable Declaration and Description
! ===========================================================================
        integer iele
        integer jele
        integer iband
        integer jband
        integer it
        integer iatom
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
        integer imu
        real ddt
        real delta
        real, dimension (nele,nele) :: suma
        real, dimension (nele,nkpoints) :: deig
        real, dimension (nele,nkpoints) :: eig
        real, dimension (3,natoms) :: dvatom
        real, dimension (3,natoms) :: v
        real, dimension (3,natoms,nele,nele) :: dgks
        real, dimension (3,natoms,nele,nele) :: g
        real, dimension (nele,nele) :: ddnac
        real, dimension (nele,nele) :: nonac
!       complex, dimension (nele, nele, nkpoints) :: c_old
!       complex, dimension (nele, nele, nkpoints) :: c_new
        complex, dimension (nele, nele, nkpoints) :: dc_na
        complex, dimension (nele, nele, nkpoints) :: dc_aux
        complex, dimension (nele, nele, nkpoints) :: c_aux
        complex aim
        complex a0
        complex a1

! Procedure
! JOM-info : we will start with a simple Verlet-like algorithm 
! c (t + dt) = c(t - dt) + 2 dt d/dt c(t)
!
! Notice that d/dt c_{ak}(t) = eigen_{k}/ih c_{ak} (t) - \Sum_j V.d_{kj}
! ===========================================================================
! JOM-test
!	write(*,*) 'entramos en evolve states'
!       write(*,*)'c_na',c_na(2,2,1),cabs(c_na(2,2,1))
 
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

        if (itime_step .eq. 1) then
         allocate (vatom_old(3,natoms))
         allocate (gks_old(3,natoms,nele,nele))
         allocate (eigen_old(norbitals,nkpoints))
         allocate (eigen_1(nele,nkpoints))
         allocate (eigen_0(nele,nkpoints))
         allocate (dnac_old(nele,nele))

         vatom_old = vatom
         gks_old = gks
         eigen_old = eigen_k
         dnac_old = dnac
        end if
           write (212,*) ' ------ the energy eigenvalues ----'
           write (212,101) (eigen_k(imu,1), imu = 1, norbitals_new)

!          write (216, '(<nele>f16.10)') (eigen_k(imu,1), imu = map_ks(1),map_ks(1)+nele-1)
           write (216, '(<nele+2>f16.10)') (eigen_k(imu,1), imu = map_ks(1)-1,map_ks(1)+nele)
! ===========================================================================
! Calculate d/dt c_{ak} at different time steps in between t and t+dt
! tt(it) = t + dt/Nsteps * it . We need to interpolate the values for
! eigen_k, vatom, gks, using their values at t and t+dt. We start using
! a simple linear interpolation. 
! ===========================================================================

        do iele = 1, nele
          iband = map_ks(iele)
	  do ikpoint = 1, nkpoints
            eigen_1(iele,ikpoint) = eigen_k(iband,ikpoint)
            eigen_0(iele,ikpoint) = eigen_old(iband,ikpoint)
          end do
!          do jele = 1, nele
!            jband = map_ks(jele)
!            do ix = 1, 3
!              do iatom = 1, nele
!                gks_1(ix,iatom,iele,jele) = gks(ix,iatom,iband,jband)
!                gks_0(ix,iatom,iele,jele) = gks_old(ix,iatom,iband,jband)
!              end do
!            end do
!          end do
        end do

! Interpolation stuff
        ddt = dt / nddt
        deig = (eigen_1 - eigen_0)/nddt
        dvatom = (vatom - vatom_old)/nddt
        dgks = (gks - gks_old)/nddt
        ddnac = (dnac - dnac_old)/nddt
! ===========================================================================
!

! ===========================================================================
! Time Loop
        do it = 1, nddt
! step 1 !!!!!!!!!!!!!!
! Interpolation 
         v = vatom_old + dvatom*(it - 1)
         g = gks_old + dgks*(it - 1)
         eig = eigen_0 + deig*(it - 1)
         nonac = dnac + ddnac*(it - 1)
         c_aux = c_na
! Calculate derivative d/dt c_na(t)
         call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,natoms,nele,nkpoints,norbitals)
!        call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints) 
! JOM-test
!        write(*,*) 'primera llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
         dc_aux = dc_na/6.0d0
! step 2 !!!!!!!!!!!!!!
! Interpolation 
         v = vatom_old + dvatom*(it - 0.5d0)
         g = gks_old + dgks*(it - 0.5d0)
         eig = eigen_0 + deig*(it - 0.5d0)
         nonac = dnac + ddnac*(it - 0.5d0)
         c_aux = c_na + dc_na*ddt*0.5d0
! Calculate derivative d/dt c_na(t)
         call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,natoms,nele,nkpoints,norbitals)
!        call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints) 
! JOM-test
!        write(*,*) 'segunda llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
         dc_aux = dc_aux + dc_na/3.0d0
! step 3 !!!!!!!!!!!!!!
! Interpolation 
         v = vatom_old + dvatom*(it - 0.5d0)
         g = gks_old + dgks*(it - 0.5d0)
         eig = eigen_0 + deig*(it - 0.5d0)
         nonac = dnac + ddnac*(it - 0.5d0)
         c_aux = c_na + dc_na*ddt*0.5d0
! Calculate derivative d/dt c_na(t)
         call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,natoms,nele,nkpoints,norbitals)
!        call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints) 
! JOM-test
!        write(*,*) 'tercera llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
         dc_aux = dc_aux + dc_na/3.0d0
! step 4 !!!!!!!!!!!!!!
! Interpolation 
         v = vatom_old + dvatom*(it)
         g = gks_old + dgks*(it)
         eig = eigen_0 + deig*(it)
         nonac = dnac + ddnac*(it - 0.5d0)
         c_aux = c_na + dc_na*ddt
! Calculate derivative d/dt c_na(t)
         call dcdt_nac (v,g,nonac,eig,c_aux,dc_na,natoms,nele,nkpoints,norbitals)
!        call dcdt_nac (nonac,eig,c_aux,dc_na,natoms,nele,nkpoints) 
! JOM-test
!        write(*,*) 'cuarta llamada a dcdt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)
         dc_aux = dc_aux + dc_na/6.0d0
!
! Integrate coefficients c_na
         c_na = c_na + dc_aux*ddt
! JOM-test
!        write(*,*) 'hacemos c_na = c_na + dc_aux*dt, it =',it
!        write(*,*)'c_na',c_na(2,2,1),c_aux(2,2,1),dc_na(2,2,1)

! end Time loop
        end do
! JOM-test
!       write(*,*)'c_na',c_na(2,2,1),cabs(c_na(2,2,1)),dc_na(2,2,1)
!       do iele = 1, nele
!	 do jele = 1, nele
!         write(*,*)'c_na',c_na(iele,jele,1)
!	 end do
!	end do

	do iele = 1, nele
	 do jele = 1, nele
          write(214,102) iele, jele, real(c_na(iele,jele,1)),aimag(c_na(iele,jele,1))
	 end do
	end do
 

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, 4(2x,f11.5))
102     format ('c_na',2i4,f11.7,f11.7)
        return
        end subroutine evolve_ks_states

