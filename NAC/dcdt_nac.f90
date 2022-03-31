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

! dcdt_nac.f90
! Program Description
! ===========================================================================
!       This routine gives the derivative wrt time for the coefficients
!       c_wf of the TD-wfs
!
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine dcdt_nac (nonac,eig,c_wf,dc_na,natoms,nele,nkpoints)


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer natoms,nele,nkpoints
        real, dimension (nele,nele) :: nonac
        real, dimension (nele,nkpoints) :: eig
        complex, dimension (nele, nele, nkpoints) :: c_wf


! Output
        complex, dimension (nele, nele, nkpoints) :: dc_na
! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: hbar = 0.6582119d0


! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
!        real, dimension (nele, nele) :: suma
        complex, dimension (nele, nele, nkpoints) :: caux
        complex aim
        complex a0
        complex a1



! Procedure
!
! Notice that d/dt c_{ak}(t) = eigen_{k}/ih c_{ak} (t) - \Sum_j V.d_{kj}
! ===========================================================================
        aim = cmplx(0.0d0, 1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

! ===========================================================================
! Non-adiabatic term: dot pruduct sum
!         suma = 0.0d0
!         do ik = 1, nele
!          do ij = 1, nele
!           do iatom = 1, natoms
!            do ix = 1, 3
!          suma(ik,ij) = suma(ik,ij) + v(ix,iatom)*g(ix,iatom,ik,ij)
!            end do
!           end do
!          end do
!         end do
! Non-adiabatic term: sum over j 
         caux = a0
         do ikpoint = 1, nkpoints
          do ia = 1, nele
           do ik = 1, nele
            do ij = 1, nele
          caux(ia,ik,ikpoint) = caux(ia,ik,ikpoint) + nonac(ik,ij)*c_wf(ia,ij,ikpoint)
            end do
           end do
          end do
         end do
! Calculate derivative d/dt c_na(t)
         dc_na = a0
         do ikpoint = 1, nkpoints
          do ia = 1, nele
           do ik = 1, nele
         dc_na(ia,ik,ikpoint) = - eig(ik,ikpoint)*aim/hbar*c_wf(ia,ik,ikpoint)   &
     &     - caux(ia,ik,ikpoint)  
! JOM-test
!        dc_na(ia,ik,ikpoint) = - eig(ik,ikpoint)*aim/hbar*c_na(ia,ik,ikpoint)
           end do
          end do
         end do
! write
!         do ia = 1, nele
!          do ik = 1, nele
!        write(*,*)'dc_na',ia,ik,dc_na(ia,ik,1), caux(ia,ik,1)
!          end do
!         end do





! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine dcdt_nac

