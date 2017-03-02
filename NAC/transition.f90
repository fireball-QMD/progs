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

! transition.f90
! Program Description
! ===========================================================================
!       This routine performs the electronic transition ij --> is
!        and rescales the velocities after the transition between
!       Kohn-Sham states, in order to conserve total energy
! The velocities are re-scaled along the direction of the nonadiabatic
! coupling (see JCP 101, 4657 (1994) )
!
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine transition (itime_step, ij, is, ikp)

        use configuration
        use nonadiabatic
        use constants_fireball
        use density
        use options
        use energy
        use interactions

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: itime_step
        integer, intent(in) :: ij
        integer, intent(in) :: is
        integer, intent(in) :: ikp


! Local Parameters and Data Declaration
        real :: ejump   
! energy difference in eV, between "state is" and "state ij" (E_s - E_j)      
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix
        integer jband
        integer kband
        real :: aa
        real :: bb
        real :: cc
        real :: ener
        real :: alfa
        real :: tkinetic
        real :: etot_before
        real :: etot_after
        real, dimension (3,natoms) :: v
        real, dimension (3,natoms) :: g
        real, parameter :: tolaa = 1.0d-08

        real, dimension (norbitals, norbitals) :: bbnkre_before

! Procedure
! ===========================================================================
! switch from ij ---> is
         jband = map_ks(ij)
         kband = map_ks(is)
         ioccupy_na (jband, ikp) = ioccupy_na (jband,ikp) - 1
         foccupy_na (jband, ikp) = foccupy_na (jband,ikp) - 0.50d0
         ioccupy_na (kband, ikp) = ioccupy_na (kband,ikp) + 1
         foccupy_na (kband, ikp) = foccupy_na (kband,ikp) + 0.50d0
!----------------------------------------------------------
! Calculate energy jump
        if (itheory .eq. 0) then
                ejump = eigen_k (kband, ikp) - eigen_k (jband, ikp)
        else
                etot_before = etot
		bbnkre_before(:,:)=bbnkre(:,:,1) 
                call scf_loop(itime_step)
                call getenergy(itime_step)
                etot_after = etot
                ejump = etot_after - etot_before
                write(*,*)'ETOT-BEFORE=',etot_before
                write(*,*)'ETOT-AFTER=',etot_after
        end if
!----------------------------------------------------------
! transform energy shift from eV to "dynamical" units ( amu*(angs/fs)**2 )
        write(*,*)'ejump (eV)=', ejump
        ener = ejump*fovermp
        write(*,*)'ener (dynamical units)=', ener

! ===========================================================================
! Find out if transition ij --> is is accesible (i.e. if there is enough
! kinetic energy 
!----------------------------------------------------------
        v = vatom
        g (:,:) = gks (:,:,ij,is)
!----------------------------------------------------------
        do iatom = 1, natoms
         write(*,*)'g',  (g (ix,iatom), ix = 1,3)
        end do
!----------------------------------------------------------
!
        aa = 0.0d0
        bb = 0.0d0
        do iatom = 1, natoms
         do ix = 1, 3
          aa = aa + 0.50d0*g(ix,iatom)*g(ix,iatom)/xmass(iatom)
          bb = bb + v(ix,iatom)*g(ix,iatom)
         end do
        end do
!----------------------------------------------------------
        write(*,*)'aa, 4*ener*aa =', aa, 4.0d0*aa*ener
        write(*,*)'bb, bb**2=', bb, bb**2
        cc = bb**2 - 4.0d0*aa*ener
        write(*,*)'cc=', cc
! ===========================================================================
       if (aa .gt. tolaa) then
!----------------------------------------------------------
! the transition is accepted
        if (cc .ge. 0.0d0) then
         write(*,*)'transition accepted'
!----------------------------------------------------------
         if (bb .ge. 0.0d0) then
          alfa = (bb - sqrt(cc))/(2.0d0*aa)
         else
          alfa = (bb + sqrt(cc))/(2.0d0*aa)
         end if
!----------------------------------------------------------
        else !(cc .ge. 0.0d0)
! the transition is NOT accepted
         write(*,*)'transition NOT accepted'
!----------------------------------------------------------
! Revert transition
         ioccupy_na (jband, ikp) = ioccupy_na (jband,ikp) + 1
         foccupy_na (jband, ikp) = foccupy_na (jband,ikp) + 0.50d0
         ioccupy_na (kband, ikp) = ioccupy_na (kband,ikp) - 1
         foccupy_na (kband, ikp) = foccupy_na (kband,ikp) - 0.50d0
         bbnkre(:,:,1)=bbnkre_before(:,:) 
!----------------------------------------------------------
! Define alfa for re-scaling velocities
! JOM-info : two possibilities: a) Do nothing; b) Reflection
! a) Do nothing

! b) Reflection (see JCP 101, 4657 (1994), pag. 4664)
          alfa = bb/aa 
      
        end if !(cc .ge. 0.0d0)
       else  !(aa .gt. tolaa)
         alfa = 0.0d0
       end if  !(aa .gt. tolaa)
         write(*,*)'alfa=', alfa
! ===========================================================================
        do iatom = 1, natoms
         write(*,*)'vatom-B',  (vatom (ix,iatom), ix = 1,3)
        end do
!----------------------------------------------------------
        tkinetic = 0.0d0
        do iatom = 1, natoms
         tkinetic = tkinetic                                            &
     &       + (0.5d0/fovermp)*xmass(iatom)                             &
     &      *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
        end do
        write(*,*)'KINETIC=',tkinetic
!----------------------------------------------------------

!----------------------------------------------------------
! RESCALING VELOCITIES
        do iatom = 1, natoms
         do ix = 1, 3
          vatom (ix,iatom) = vatom (ix,iatom) - alfa*g(ix,iatom)/xmass(iatom)
         end do
        end do
!----------------------------------------------------------
!----------------------------------------------------------
        do iatom = 1, natoms
         write(*,*)'vatom-A',  (vatom (ix,iatom), ix = 1,3)
        end do
!----------------------------------------------------------
!----------------------------------------------------------
        tkinetic = 0.0d0
        do iatom = 1, natoms
         tkinetic = tkinetic                                            &
     &       + (0.5d0/fovermp)*xmass(iatom)                             &
     &      *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
        end do
        write(*,*)'KINETIC=',tkinetic
!----------------------------------------------------------




! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine transition

