! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

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

! find_neighPP_max.f90
! Program Description
! ===========================================================================
!       Finds all the maximum number of neighbors (Pseudopotential)to atoms 
! in the central cell.
!
!
!                   atom_VNL
!                   +     +
!                 +          +
!               +               + 
!             +                    +
!           +                         +
!        atom 1                          +
!                                       atom 2
! 
! Generaly overlap of psi_1 and psi_2 can be zero !!!. So we have different 
! number of pairs in case of PP interaction.
! Pseudopotential interaction are calculated only for distance rcutoff_i + rc_PP.
! This is different in common interactions, where we have rcutoff_i + rcutoff_j.
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine find_neighPP_max (nprocs, my_proc, iordern, icluster)
        use configuration
        use dimensions
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: iordern
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs

!$ volatile icluster, natoms, nprocs, my_proc, iordern
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in2
        integer jatom
        integer mbeta
        integer natomsp
        integer num_neigh
        integer num_neigh_vdw

        real distance2
        real range2
        real rcutoff_i
        real rcutoff_j

! Procedure
! ===========================================================================

        if (my_proc .eq. 0)                                                  &
     &   write (*,*) ' Determine maximum number of PP-neighbors. ' 

        if (icluster .eq. 1) mbeta_max = 0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if


! Note: We need two loops to determine max number of neighbors, namely 
! for <phi_i|Psi_j> and <Psi_i|phi_j>, where phi_i means atomic wavefunction 
! and Psi_i pseudowave function. We get larger one.

! Loop over all atoms.
        neighPP_max = -99
!$omp parallel do private (num_neigh, rcutoff_i, rcutoff_j)   &
!$omp&            private (in1, in2, distance2, range2)
! Loop for <phi_i|Psi_j>
        do iatom = iatomstart, iatomstart - 1 + natomsp
         num_neigh = 0
         rcutoff_i = 0.0d0
         in1 = imass(iatom)
         do imu = 1, nssh(in1)
          if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
         end do

! Loop over all possible neighbors (VNL atoms)
         do mbeta = 0, mbeta_max
          do jatom = 1, natoms

           in2 = imass(jatom)
           rcutoff_j = rc_PP(in2)

! Find the distance from (mbeta,jatom) to (0,iatom)
           distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &                + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &                + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2

! Add a small displacement to the sum of cutoffs. 
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2

           if (distance2 .le. range2) then
            num_neigh = num_neigh + 1
           end if
          end do
         end do

! Maximum number of neighbors thus far.
!$omp atomic
         neighPP_max = max(neighPP_max, num_neigh)
        end do

! Loop for <Psi_i|phi_j>
!$omp parallel do private (num_neigh, rcutoff_i, rcutoff_j)   &
!$omp&  private (in1, in2, imu, distance2, range2)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         num_neigh = 0
         in1 = imass(iatom)
         rcutoff_i = rc_PP(in1)

! Loop over all possible neighbors (VNL atoms)
         do mbeta = 0, mbeta_max
          do jatom = 1, natoms

           rcutoff_j = 0.0d0  
           in2 = imass(jatom)
           do imu = 1, nssh(in2)
            if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
           end do

! Find the distance from (mbeta,jatom) to (0,iatom)
           distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &                + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &                + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2

! Add a small displacement to the sum of cutoffs. 
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
           if (distance2 .le. range2) then
            num_neigh = num_neigh + 1
           end if
          end do
         end do

! Maximum number of neighbors thus far.
!$omp atomic
         neighPP_max = max(neighPP_max, num_neigh)
        end do
!       if (iordern .eq. 1) call find_neigh_max_ordern_final()

! Format Statements
! ===========================================================================
 
        return
      end subroutine find_neighPP_max
