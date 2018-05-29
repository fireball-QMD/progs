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


! common_neighborsPP.f90
! Program Description
! ===========================================================================
!       Given atom ialp, we find pairs of atoms that are common neighbors
! of atom ialp. Note that on atom ialp is placed VNL.
!
! Find all common neighbors (Pseudopotential) of atom alpha 
! (we call its i value ialp). Atom alpha is in the central cell, and 
! the common neighbors are at l-vector-sub-i,i and l-vector-sub-j,j
!
! Note:
! You should remember that in counting the common neighbors, we do not include
! any pairs in which one or both of the atoms is ontop atom ialp. Thus we are 
! keeping only third party common neighbors.
!
! ncomnPP(ialp) = num =total # of common neighbors to atom (0,ialp)
! ncomjPP(ialp,mneigh,1) = iatom value for mneigh'th pair of common neighbors
! ncomjPP(ialp,mneigh,2) = jatom value for mneigh'th pair of common neighbors
! ncombPP(ialp,mneigh,1) = mbeta-sub-iatom value for m'th pair of common neigh
! ncombPP(ialp,mneigh,2) = mbeta-sub-jatom value for m'th pair of common neigh
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
        subroutine common_neighborsPP (nprocs, my_proc, iordern, iwrtneigh, icluster )
        use configuration
        use dimensions
        use interactions
        use neighbor_map 
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: iordern
        integer, intent (in) :: iwrtneigh
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs
        integer, intent (in) :: icluster

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer ibeta
        integer imu
        integer in1, in2, indna
        integer ineigh
        integer jatom
        integer jbeta
        integer jneigh
        integer katom
        integer kbeta
        integer kneigh
        integer num_neigh
        integer iatomstart, natomsp
        integer mbeta

        real  distance
        real, dimension (3) :: dvec
        real, dimension (3) :: diff
        real, dimension (3) :: vec1, vec2, vec3, vec

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        if(icluster .eq. 1) mbeta_max = 0

!        if (my_proc .eq. 0) then
!         write (*,*) ' Welcome to common_neighborsPP - determine mapping of '
!         write (*,*) ' common neighbors for Pseudopotential. '
!        end if

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


! ****************************************************************************
!
!      P P        C O M M O N      N E I G H B O R S     L I S T
! ****************************************************************************
! SCHEMATIC VIEW:
!
!         <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!          ialp
!           * R2 (NL)            
!          / \                   nPPx  ->  <Psi_k|phi_i>
!         /   \                  nPPx  ->  <Psi_k|phi_j>
!   nPPx /     \ nPPx            
!       /       \ 
!      /         \
!     + <phi_i|   + |phi_j>
!   iatom       jatom
!    R1           R3
! ****************************************************************************


! Now we build the common neighbors list
!$omp parallel do private (num_neigh, ineigh, iatom, ibeta, in1, jneigh) &
!$omp&   private (jatom, jbeta, in2, kneigh, katom, kbeta) &
!$omp&   private (vec, vec1, vec2, vec3, diff)
        do ialp = iatomstart, iatomstart - 1 + natomsp
           
         vec2(:) = ratom(:,ialp)
! The variable num_neigh counts neighbors.
         num_neigh = 0

! 1.loop over <Psi_k|phi_i>  -> nPPx 
         do ineigh = 1, nPPxn(ialp)
          
          iatom = nPPx_j(ineigh,ialp)
          ibeta = nPPx_b(ineigh,ialp)
          in1 = imass(iatom)
! Keep only third party common neighbors. It means iatom /= ialp.
          if (.not. (iatom .eq. ialp .and. ibeta .eq. 0)) then
           vec1(:) = ratom(:,iatom) + xl(:,ibeta)
! Loop over all known neighbors of atom ialp again, and check to see if atoms 
! iatom and jatom are neighbors. If so, then we have found a pair of common 
! neighbors. Of course, if iatom is jatom in all respects, that does not count.
! 2.loop over <Psi_k|psi_j>  -> nPPx
           do jneigh = 1, nPPxn(ialp)

            jatom = nPPx_j(jneigh,ialp)
            jbeta = nPPx_b(jneigh,ialp)
            in2 = imass(jatom)

! Keep only third party common neighbors. It means jatom /= ialp .and.
! jatom /= iatom
            if (.not. (jatom .eq. ialp .and. jbeta .eq. 0)                  &
     &          .and. (ineigh .ne. jneigh)) then

             vec3(:) = ratom(:,jatom) + xl(:,jbeta)

! If we're here it means we've found third party common neighbors. Increase 
! number of instances.
             num_neigh = num_neigh + 1

! Check that # of common neighbor pairs is less than dimensions.
             if (num_neigh .gt. neighPP_max**2) then
               write (*,*) ' Oh no. In common_neighbors.f90, we have too many '
               write (*,*) ' common neighbors. In MODULES/dimensions.f90 '
               write (*,*) ' dimension neigh_max**2 = ', neighPP_max**2
               write (*,*) ' So far (*but still counting) we have '
               write (*,*) ' num_neigh = ', num_neigh, ' within neighbors.f90! '
! FIXME can't have a stop statement in an OpenMP loop
              stop
             end if ! if(num_neigh)

! Put iatom in first spot
             neighPP_comj(1,num_neigh,ialp) = iatom
             neighPP_comb(1,num_neigh,ialp) = ibeta
! Put jatom in second spot
             neighPP_comj(2,num_neigh,ialp) = jatom
             neighPP_comb(2,num_neigh,ialp) = jbeta

! We also need to know for a given ialp and (iatom,ibeta), what is the m value
! for (jatom,jbeta) with respect to iatom. That is, jatom is the m'th neighbor
! of iatom. What is m?

! Set to a crazy value, in case loop fails.
              diff = vec3 - vec1
              neighPP_comm(num_neigh,ialp) = -9999
              do kneigh = 1, neighPPn(iatom)
               katom = neighPP_j(kneigh,iatom)
               kbeta = neighPP_b(kneigh,iatom)
               vec(:) = xl(:,kbeta) + ratom(:,katom) - ratom(:,iatom)
               if ((abs(vec(1) - diff(1)) .lt. 1.0d-4) .and.                & 
     &             (abs(vec(2) - diff(2)) .lt. 1.0d-4) .and.                &
     &             (abs(vec(3) - diff(3)) .lt. 1.0d-4))                     &
     &          neighPP_comm(num_neigh,ialp) = kneigh
              end do ! do kneigh
! We check to see if it really was a neighbor. (could be very close to cutoff)
              if (neighPP_comm(num_neigh,ialp) .eq. -9999)                  &
     &         num_neigh = num_neigh - 1 
            end if ! if(.not.(jatom .eq. ialp)
           end do ! do jneigh 
          end if ! if(.not.(iatom .eq. ialp) 
         end do ! do ineigh 
         neighPP_comn(ialp) = num_neigh
        end do


        if (iordern .eq. 1)                                                  &
     &   call common_neighbors_ordern_final (natoms, nprocs, my_proc)

! Option for writing out the results for debugging.
        if (iwrtneigh .eq. 1 .and. my_proc.eq.0) then
         write (*,*) '  '
         write (*,*) ' Common neighbors (Pseudopotential) of each atom: '
         do ialp = 1, natoms
          num_neigh = neighPP_comn(ialp)
          write (*,*) '  '
          write (*,*) ' num_neigh = ', num_neigh
          write (*,100)
          write (*,101)
          do ineigh = 1, num_neigh
           iatom = neighPP_comj(1,ineigh,ialp)
           ibeta = neighPP_comb(1,ineigh,ialp)
           jatom = neighPP_comj(2,ineigh,ialp)
           jbeta = neighPP_comb(2,ineigh,ialp)
           kneigh = neighPP_comm(ineigh,ialp)

           distance =                                                        &
     &      sqrt(((ratom(1,iatom) + xl(1,ibeta))                             &
     &            - (ratom(1,jatom) + xl(1,jbeta)))**2                       &
     &           + ((ratom(2,iatom) + xl(2,ibeta))                           &
     &              - (ratom(2,jatom) + xl(2,jbeta)))**2                     &
     &           + ((ratom(3,iatom) + xl(3,ibeta))                           &
     &              - (ratom(3,jatom) + xl(3,jbeta)))**2)
 
           dvec(:) = xl(:,jbeta) + ratom(:,jatom) - xl(:,ibeta) - ratom(:,iatom)
           write (*,102) ialp, ineigh, iatom, ibeta, jatom, jbeta, kneigh,   &
     &                   distance, dvec
          end do
         end do
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, ' ialp  # iatom  ibeta jatom jbeta kneigh', 1x,       &
     &              ' distance           vector ')
101     format (2x, 75('='))
102     format (2x, 2i4, 2x, 4(i3, 3x), i4, 3x, f8.4, 2x, 3f8.4)


        return
      end subroutine common_neighborsPP
