! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! build_gover1c.f90
! Program Description
! ===========================================================================
!       This routine takes the input data 
! and calculates the one center box for 
! the vector < G mu | nu > , where G stands for
! Gradient, in the variable gover1c(ix,imu,inu,iatom).
! (G wrt iatom-position)
! ===========================================================================
!
! JOM-warning : there is an overall minus (-) sign
! due to Grad wrt R (atomic position) = - Grad wrt r (electron position)
! and another overall minus (-) sign
! due to  < G mu | nu >  = - < mu | G nu > (this last one is the
! calculated in the notes).
! ===========================================================================
!
! JOM-info : the ordering of the orbitals in the different shells is
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                        1   2     3
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                  1     2   3     4      5

! while for the derivatives it is
! d/dx : 1 ; d/dy : 2 ; d/dz : 3.
!
! ===========================================================================
! Code written by:
! Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_gover1c (in1)
        use charges
        use density
        use interactions
        use nonadiabatic
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
!       integer, intent (in) :: iatom

! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ind1
        integer ind2
        integer inu
        integer issh
        integer jssh
        integer l1
        integer l2
        integer n1
        integer n2

!       real, parameter :: r3 = 1.732.........
        real ur3
        real ur5
        real ur15
        real f1
        real f2

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! constants_fireball
        ur3 = 1.0d0/dsqrt(3.0d0)
        ur5 = 1.0d0/dsqrt(5.0d0)
        ur15 = 1.0d0/dsqrt(15.0d0)
! Initialize
        gover1c = 0.0d0

! JOM-warning : there is an overall minus (-) sign
! due to Grad wrt R (atomic position) = - Grad wrt r (electron position)
! and another overall minus (-) sign
! due to  < G mu | nu >  = - < mu | G nu >


! Loop over shells i-atom        
        n1 = 0
        do issh = 1, nssh(in1)

! Number of orbitals per the shell in x-dim
         l1 = lssh(issh,in1)
         n2 = 0

! Loop over shells i y-dimension
         do jssh = 1, nssh(in1)

! Number of orbitals per the shell in y-dim
          l2 = lssh(jssh,in1)
! The F1 and F2 radial integrals
         f1 = f1nac1c(in1,issh,jssh)
         f2 = f2nac1c(in1,issh,jssh)
!        write(*,*) 'f1 =', in1, issh, jssh, f1
!        write(*,*) 'f2 =', in1, issh, jssh, f2
          if (l1.eq.0 .and. l2.eq.1) then
            gover1c(1, n1 + 1, n2 + 3) = ur3*(f1 + 2.0d0*f2)
            gover1c(2, n1 + 1, n2 + 1) = ur3*(f1 + 2.0d0*f2)
            gover1c(3, n1 + 1, n2 + 2) = ur3*(f1 + 2.0d0*f2)
          else if (l1.eq.1 .and. l2.eq.0) then
            gover1c(1, n1 + 3, n2 + 1) = ur3*f1
            gover1c(2, n1 + 1, n2 + 1) = ur3*f1
            gover1c(3, n1 + 2, n2 + 1) = ur3*f1
          else if (l1.eq.1 .and. l2.eq.2) then
            gover1c(1, n1 + 1, n2 + 1) = ur5*(f1 + 3.0d0*f2)
            gover1c(1, n1 + 2, n2 + 4) = ur5*(f1 + 3.0d0*f2)
            gover1c(2, n1 + 3, n2 + 1) = ur5*(f1 + 3.0d0*f2)
            gover1c(2, n1 + 2, n2 + 2) = ur5*(f1 + 3.0d0*f2)
            gover1c(3, n1 + 1, n2 + 2) = ur5*(f1 + 3.0d0*f2)
            gover1c(3, n1 + 3, n2 + 4) = ur5*(f1 + 3.0d0*f2)
!
            gover1c(3, n1 + 2, n2 + 3) = ur15*(2.0d0*f1 + 6.0d0*f2)
!
            gover1c(1, n1 + 3, n2 + 3) = ur15*(-1.0d0*f1 - 3.0d0*f2)
            gover1c(2, n1 + 1, n2 + 3) = ur15*(-1.0d0*f1 - 3.0d0*f2)
!
            gover1c(1, n1 + 3, n2 + 5) = ur5*(1.0d0*f1 + 3.0d0*f2)
            gover1c(2, n1 + 1, n2 + 5) = ur5*(-1.0d0*f1 - 3.0d0*f2)
            else if (l1.eq.2 .and. l2.eq.1) then
            gover1c(1, n1 + 1, n2 + 1) = ur5*(f1 - f2)
            gover1c(1, n1 + 4, n2 + 2) = ur5*(f1 - f2)
            gover1c(2, n1 + 1, n2 + 3) = ur5*(f1 - f2)
            gover1c(2, n1 + 2, n2 + 2) = ur5*(f1 - f2)
            gover1c(3, n1 + 2, n2 + 1) = ur5*(f1 - f2)
            gover1c(3, n1 + 4, n2 + 3) = ur5*(f1 - f2)
!
            gover1c(3, n1 + 3, n2 + 2) = 2.0d0*ur15*(f1 - f2)
!
            gover1c(1, n1 + 3, n2 + 3) = ur15*(f2 - f1)
            gover1c(2, n1 + 3, n2 + 1) = ur15*(f2 - f1)
!
            gover1c(1, n1 + 5, n2 + 3) = ur5*(f1 - f2)
            gover1c(2, n1 + 5, n2 + 1) = ur5*(f2 - f1)
          end if
          n2 = n2 + 2*l2 + 1
         end do !do jssh = 1, nssh(in1)
         n1 = n1 + 2*l1 + 1
        end do !do issh = 1, nssh(in1)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
        end subroutine build_gover1c
