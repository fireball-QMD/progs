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
        integer it
        integer iatom
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
        real ddt
        real delta
        real, dimension (norbitals, norbitals) :: suma
        real, dimension (norbitals,nkpoints) :: deig
        real, dimension (norbitals,nkpoints) :: eig
        real, dimension (3,natoms) :: dvatom
        real, dimension (3,natoms) :: v
        real, dimension (3,natoms,norbitals,norbitals) :: dgks
        real, dimension (3,natoms,norbitals,norbitals) :: g
        complex, dimension (norbitals, norbitals, nkpoints) :: c_old
        complex, dimension (norbitals, norbitals, nkpoints) :: c_new
        complex, dimension (norbitals, norbitals, nkpoints) :: dc_na
        complex, dimension (norbitals, norbitals, nkpoints) :: caux
        complex aim
        complex a0
        complex a1



! Procedure
! JOM-info : we will start with a simple Verlet-like algorithm 
! c (t + dt) = c(t - dt) + 2 dt d/dt c(t)
!
! Notice that d/dt c_{ak}(t) = eigen_{k}/ih c_{ak} (t) - \Sum_j V.d_{kj}
! ===========================================================================
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

        if (itime_step .eq. 1) then
         allocate (ratom_old(3,natoms))
         allocate (vatom_old(3,natoms))
         allocate (gks_old(3,natoms,norbitals,norbitals))
         allocate (eigen_old(norbitals,nkpoints))

         allocate (bbnkre_old(norbitals,norbitals,nkpoints))
         allocate (blowre_old(norbitals,norbitals,nkpoints))
         ratom_old = ratom
         vatom_old = vatom
         gks_old = gks
         eigen_old = eigen_k
         bbnkre_old = bbnkre
         blowre_old = blowre
        end if
! ===========================================================================
! Calculate d/dt c_{ak} at different time steps in between t and t+dt
! tt(it) = t + dt/Nsteps * it . We need to interpolate the values for
! eigen_k, vatom, gks, using their values at t and t+dt. We start using
! a simple linear interpolation. 
! ===========================================================================
! Interpolation stuff
        ddt = dt / nddt
        deig = eigen_k - eigen_old 
        dvatom = vatom - vatom_old
        dgks = gks - gks_old
! ===========================================================================
! Initialize : we need to carefully   define c_na and d/dt c_na for it=1
        c_old = c_na
         delta = 0.0d0
! Interpolation : values for it=0
         v = vatom_old
         g = gks_old
         eig = eigen_old
! Calculate derivative d/dt c_na at it=0
         call dcdt_nac (v,g,eig,c_na,dc_na,natoms,norbitals,nkpoints)
!----------------------------------------------------------
! Calculate derivative d/dt c_na at it=0.5
! first we need c_na for it=0.5
        c_na = c_na + ddt/2.0d0*dc_na
         delta = 0.50d0/real(nddt)
! Interpolation 
         v = vatom_old + dvatom*delta
         g = gks_old + dgks*delta
         eig = eigen_old + deig*delta
         call dcdt_nac (v,g,eig,c_na,dc_na,natoms,norbitals,nkpoints)
!----------------------------------------------------------
! Calculate c_na for it=1
        c_na = c_old + ddt*dc_na
!
! Now we have c_old = c_na(0) and c_na(1), and are ready to go to the
! main loop and calculate c_na(2), etc
!

! ===========================================================================
! Time Loop
        do it = 2, nddt
         delta = real(it)/real(nddt)
! Interpolation 
         v = vatom_old + dvatom*delta
         g = gks_old + dgks*delta
         eig = eigen_old + deig*delta
! call subroutine to
! Calculate derivative d/dt c_na(t)
         call dcdt_nac (v,g,eig,c_na,dc_na,natoms,norbitals,nkpoints)
! write
!        if (it .eq. 1) then
!         do ia = 1, norbitals
!          do ik = 1, norbitals
!        write(*,*)'dc_na',ia,ik,dc_na(ia,ik,1)
!          end do
!         end do
!        end if
! Integrate coefficients c_na
          c_new = c_old + 2.0d0*ddt*dc_na
! Update for next time step
          c_old = c_na
          c_na = c_new

! end Time loop
        end do





! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine evolve_ks_states

