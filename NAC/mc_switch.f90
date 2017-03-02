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

! mc_switch.f90
! Program Description
! ===========================================================================
!       This routine determines the hoppings switches between
!       Kohn-Sham states based on a Monte-Carlo approach (see J.C.
!       Tully, JCP 93, 1061 (1990).
!
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine mc_switch (xr, n, pr, ij, ikp, is)

        use density
        use nonadiabatic

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in) ::  xr
        integer, intent(in) :: n
        real, dimension (n), intent(in) :: pr
        integer, intent(in) :: ij
        integer, intent(in) :: ikp
! Output
        integer, intent(out) :: is


! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        integer :: j, k
        integer :: jband, kband
        real :: aux



! Procedure
! ===========================================================================
! Check that sum of probabilities is smaller than 1
        aux = 0.0d0
        do k = 1, n
! Consider only allowed transitions
         jband = map_ks(ij)
         kband = map_ks(k)
         if (ioccupy_na(jband,ikp) .gt. 0) then     
          if (ioccupy_na(kband,ikp) .lt. 2) then     
           aux = aux + pr(k)
          end if
         end if
        end do
        if (aux .gt. 1.0d0) then
         write(*,*)'sum of probabilities greater than 1'
         write(*,*)'in mc_switches.f90'
         write(*,*)'total probabilty', aux
         write(*,*)'for state', ij
         do k = 1, n
         write(*,*)'prob',k,pr(k)
         end do
!        write(*,*)'must stop'
!        stop
        end if
!----------------------------------------------------------
        is = 0
        aux = 0.0d0
        do k = 1, n
! Consider only allowed transitions
         jband = map_ks(ij)
         kband = map_ks(k)
         if (ioccupy_na(jband,ikp) .gt. 0) then     
          if (ioccupy_na(kband,ikp) .lt. 2) then    
           aux = aux + pr(k)
           if (aux .gt. xr) then
            is = k
         do j = 1, n
         write(*,*)'prob',j,pr(j)
         end do
            exit
           end if
          end if
         end if
        end do
!----------------------------------------------------------
! If is = 0, no switch (hopping) between states
! Otherwise, switch from current state ( "ij" = map_ks(iele) in
! fewest_switches subroutine) to state "is"
!----------------------------------------------------------

        return
        end subroutine mc_switch

