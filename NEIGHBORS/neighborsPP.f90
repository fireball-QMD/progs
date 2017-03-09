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


! neighbors.f90
! Program Description
! ===========================================================================
!  Find all neighbors to atoms in the central cell in case of Kleiman Bylander 
! fully separable pseudopotential. The distance of neighbor PP pair must be 
! less then sum of wave function and PP radius cutoff   
! An atom is at lattice vector xl(mbeta), and basis ratom(iatom).  
! We refer to this as (mbeta,iatom). An atom in the central cell is 
! at (0,iatom).
! Find all PP-neighbors to atom (0,iatom).
!
! neighnPP(i)=# of neighbors of atom i
! neighjPP(i,m)= j-sub-m, the j value of the m'th neighbor.
! neighbPP(i,m)= beta-sub-m, the beta value for the m'th neighbor.
!
!       The important quantity here is neighj, which indicates which basis
! vector we have: This identifies the species a through the array 
! imass (neighjPP(iatom,ineigh)) is the type (1 or 2) of the ineigh'th 
! PP-neighbor u of atom iatom.  Furthermore, this atom is located at 
! xl(:,neighbPP(iatom,ineigh)) + ratom(:,neighjPP(iatom,ineigh)).
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
        subroutine neighborsPP (nprocs, my_proc, iordern, icluster,   &
     &                        iwrtneigh)
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
        integer, intent (in) :: iwrtneigh
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs


!$ volatile rcutoff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ineigh
        integer ibeta
        integer iatomstart
        integer imu
        integer in1
        integer in2
        integer jatom
        integer jneigh
        integer jbeta
        integer katom
        integer kneigh
        integer kbeta
        integer ialp
        integer mbeta
        integer natomsp
        integer num_neigh

        real distance
        real distance2
        real range2
        real rcutoff_i
        real rcutoff_j
        real, dimension (3) :: vec1, vec2, vec3, vec

        logical flag

! Procedure
! ===========================================================================
        if (my_proc .eq. 0)                                                  &
     &  write (*,*) ' Welcome to neighborsPP - determine mapping of PP neighbors. '
        if (icluster .eq. 1) mbeta_max = 0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
     &                + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   LIST nPP
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! First we'are going to build nPP list of neighbors. The list includes all
! neighbors of iatom with nonzero overlap <phi_i|Psi_j>, where
! phi_i .. is atomic wave function and Psi_j is pseudo wave function

! Loop over all atoms.
!$omp parallel do private (imu, in1, in2, jatom, mbeta, num_neigh) &
!$omp   private (rcutoff_i, rcutoff_j, distance2, distance, range2)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         num_neigh = 0
         rcutoff_i = 0.0d0
         in1 = imass(iatom)
! find max Rc of atomic wave function phi_i
         do imu = 1, nssh(in1)
          if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
         end do

! Loop over all possible neighbors (VNL atoms)
         do mbeta = 0, mbeta_max
          do jatom = 1, natoms

           in2 = imass(jatom) 
! Rc_PP of pseudowave functio Psi_j
           rcutoff_j = rc_PP(in2)

! Find the distance from (mbeta,jatom) to (0,iatom)
           distance2 = (ratom(1,iatom) - (xl(1,mbeta) + ratom(1,jatom)))**2  &
     &                + (ratom(2,iatom) - (xl(2,mbeta) + ratom(2,jatom)))**2 &
     &                + (ratom(3,iatom) - (xl(3,mbeta) + ratom(3,jatom)))**2
           distance = sqrt(distance2)

! Add a small displacement to the sum of cutoffs. 
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2

           if (distance2 .le. range2) then
            if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and.        &
     &          iatom .ne. jatom) then
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
            end if 
            num_neigh = num_neigh + 1
 
! The num_neigh'th neighbor to (0,iatom) at (mbeta,jatom)
            nPP_j(num_neigh,iatom) = jatom
            nPP_b(num_neigh,iatom) = mbeta
           end if ! if(distance2 .le. range2)

          end do ! do jatom 
         end do ! do mbeta
 
! The number of neighbors to atom iatom is num_neigh.
! Remember, that in the total count, the atom itself is a neighbor!
         nPPn(iatom) = num_neigh

        end do ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   LIST nPPx
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Now we'are going to build nPPx list of neighbors. The list includes all
! neighbors of iatom with nonzero overlap <Psi_i|phii_j>, where
! phi_j .. is atomic wave function centered on j-site 
! and Psi_i is pseudoatomic wave function centered now on i-site.
! Be carreful to distinguish between nPP (<phi_i|Psi_j>) 
! and nPPx (<Psi_i|phi_j>). These two lists are not identical if we have 
! more species in the system vwith different Rc and Rc_PP

! Loop over all atoms.
!$omp parallel do private (imu, in1, in2, jatom, mbeta, num_neigh) &
!$omp   private (rcutoff_i, rcutoff_j, distance2, distance, range2)
        do iatom = iatomstart, iatomstart - 1 + natomsp

         num_neigh = 0
         in1 = imass(iatom)
! Rc_PP of pseudowave functio Psi_i
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
           distance = sqrt(distance2)
! Add a small displacement to the sum of cutoffs. 
           range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2

           if (distance2 .le. range2) then
            if (distance2 .lt. 0.7d0 .and. distance .gt. 1.0d-4 .and.        &
     &          iatom .ne. jatom) then
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' WARNING - atoms dangerously close! '
             write (*,*) ' iatom, jatom, distance = ', iatom, jatom, distance
            end if 

            num_neigh = num_neigh + 1
! The num_neigh'th neighbor to (0,iatom) at (mbeta,jatom)
            nPPx_j(num_neigh,iatom) = jatom
            nPPx_b(num_neigh,iatom) = mbeta
           end if ! if(distance2 .le. range2)

          end do ! do jatom 
         end do ! do mbeta
 
! The number of neighbors to atom iatom is num_neigh.
! Remember, that in the total count, the atom itself is a neighbor!
         nPPxn(iatom) = num_neigh

        end do ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   LIST neighPP
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OK, we have two list of neighbors nPP <phi_i|Psi_j> and nPPx <Psi_i|phi_j>.
! But what we really need is a following list of pairs <phi_i|Vnl|phi_j> to
! be able map PP matrices to global Matrix <phi_i|V|phi_j>. 
! In this list, let's call it neighPP, we're going to include all possible 
! pairs 1center, 2center and 3center. Also we need to perform mapping from
! other lists to this central one.
!
! ****************************************************************************
! SCHEMATIC VIEW:
!
!           <Psi_k|        <phi_i|Vnl|phi_j> = <phi_i|Psi_k><Psi_k|phi_j>
!           ialp
!           * Rk (NL)            
!          / \                   nPP  ->  <phi_i|Psi_k>
!         /   \                  nPPx ->  <Psi_k|phi_j>
!    nPP /     \ nPPx            
!       /       \ 
!      /         \
!     + <phi_i|   + |phi_j>
!     iatom       jatom
!
! ****************************************************************************

! Loop over atoms
!$omp parallel do private (num_neigh, ialp, ibeta, jatom, jbeta, mbeta)     &
!$omp   private (katom, kbeta, kneigh, flag, distance, vec, vec1, vec2, vec3)
        do iatom = iatomstart, iatomstart - 1 + natomsp

! The variable num_neigh counts neighbors.
         num_neigh = 0
         vec1(:) = ratom(:,iatom)

! Loop over nPP ~ <phi_i|Psi_k>
         do ineigh = 1, nPPn(iatom)
           ialp = nPP_j(ineigh,iatom)
           ibeta = nPP_b(ineigh,iatom)
           vec2(:) = ratom(:,ialp) + xl(:,ibeta)

! Loop over nPPx ~ <Psi_k|phi_j>
           do jneigh = 1, nPPxn(ialp)
             jatom = nPPx_j(jneigh,ialp)
             jbeta = nPPx_b(jneigh,ialp)
! Keeping in mind ialp is not in base unit cell, so we have to add lattice 
! vector of ialp's unit cell to get real lattice vector with respect to 
! base unit cell
             vec3(:) = xl(:,jbeta) + xl(:,ibeta)
!             write (*,*) '----------------------------------------------'
!             write (*,*) 'distance3 =', sqrt(vec3(1)**2 + vec3(2)**2 + vec3(3)**2)
!             write (*,*) ' vec3 :', vec3(:)
!             write (*,*) '----------------------------------------------'

! now we have to find mbeta, which is the real lattice vector of jatom with
! respect to base unit cell.
! loop over all image unit cells
              do mbeta = 0 , mbeta_max
!                write (*,*) 'mbeta =', mbeta
                vec(:) = xl(:,mbeta)
!                vec(:) = ratom(:,jatom) + xl(:,mbeta)
!                write (*,*) 'vec_diff', vec3(1)-vec(1), vec3(2)-vec(2),vec3(3)-vec(3)
                
                distance = sqrt( ( vec(1) - vec3(1) )**2.0d0  +         &
     &                           ( vec(2) - vec3(2) )**2.0d0  +         &
     &                           ( vec(3) - vec3(3) )**2.0d0  ) 
                if(distance .lt. 0.0001) exit
              enddo

              if (mbeta .gt. mbeta_max) then 
               write (*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++'
               write (*,*) ' mbeta_max =', mbeta_max
               write (*,*) ' vec_look =', vec3(:)
               write (*,*) ' Fireball cannot find desired neighboring cell. '
               write (*,*) ' This propably means you need to increase number '
               write (*,*) ' of the periodic unit cell box and to recompile the code.'
               write (*,*) ' For details see INITIALIZERS/initboxes.f90 file.'
               write (*,*) ' Fireball is going to die!'
               write (*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++'
               stop
              endif

! This is not all. Now we have to check if the jatom is not already inlcuded 
! in the list, so let's serch in temporary list of pairs neighPP
              flag = .false.  
              do kneigh = 1, num_neigh
                 katom = neighPP_j(kneigh,iatom)
                 kbeta = neighPP_b(kneigh,iatom)
                 if( katom .eq. jatom .and. kbeta .eq. mbeta ) flag = .true.
              enddo ! do kneigh

! there isn't already included on the list 
              if(.not. flag ) then 
                 num_neigh = num_neigh + 1
! Check that # of common neighbor pairs is less than dimensions.
                 if (num_neigh .gt. neighPP_max**2) then
                   write (*,*) 'Oh no. In common_neighbors.f90, we have too many'
                   write (*,*) 'common neighbors. In MODULES/dimensions.f90 '
                   write (*,*) 'dimension neigh_max**2 = ', neighPP_max**2
                   write (*,*) 'So far (*but still counting) we have '
                   write (*,*) 'num_neigh = ',num_neigh,' within neighbors.f90!'
! FIXME can't have a stop statement in an OpenMP loop
                   stop
                 end if

! Put jatom in second spot
                 neighPP_j(num_neigh,iatom) = jatom
                 neighPP_b(num_neigh,iatom) = mbeta
              endif ! if(flag)
           end do ! do jneigh
         end do ! do ineigh
! Save number of pairs of iatom
         neighPPn(iatom) = num_neigh
        end do ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   MAP nPPx_point  (nPPx -> nPP)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now we need mapping between nPP and nPPx. In the case of ontop left case
! we have <phi_i|Psi_i><Psi_i|phi_j>, so we will calculate it 
! if <Psi_i|phi_j> /= 0. OK, but matrix sVNL gives you only <phi_j|Psi_i> 
! overlap. So we need contection between the nPPx and nPP lists.

! Loop over atoms
!$omp parallel do private (jatom, jbeta, ineigh, kneigh, katom, kbeta) &
!$omp  private (vec1, vec2, distance)
        do iatom = 1, natoms
           vec1(:) = ratom(:,iatom)
! Loop over nPPx-neighbors <Psi_i|phi_j> of iatom
           do ineigh = 1,nPPxn(iatom)
              jatom = nPPx_j(ineigh,iatom)
              jbeta = nPPx_b(ineigh,iatom)
! Loop over nPP-neighbors <phi_j|Psi_k> of jatom
              do kneigh = 1,nPPn(jatom)
                 katom = nPP_j(kneigh,jatom)
                 kbeta = nPP_b(kneigh,jatom)
                 vec2(:) = ratom(:,katom) + xl(:,kbeta) + xl(:,jbeta)
                 distance = sqrt( ( vec1(1) - vec2(1) )**2.0d0  +      &
     &                            ( vec1(2) - vec2(2) )**2.0d0  +      &
     &                            ( vec1(3) - vec2(3) )**2.0d0  )
                 if (distance .lt. 0.0001d0) then
                    nPPx_point(ineigh,iatom) = kneigh
                    exit
                 endif ! if (distance)
              enddo ! do jneigh
           enddo ! do ineigh
        enddo ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   MAP nPPx_map  (nPPx -> neighPP)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! We need also a knowledge how place a local matrix into the global matrix.
! So, we're going to map the list nPPx to the global list neighPP
! Loop over atoms
!$omp parallel do private (jatom, jbeta, ineigh, kneigh, katom, kbeta) &
!$omp  private (vec1, vec2, distance)
        do iatom = 1, natoms
! Loop over nPPx-neighbors <Psi_i|phi_j> of iatom
           do ineigh = 1,nPPxn(iatom)
              jatom = nPPx_j(ineigh,iatom)
              jbeta = nPPx_b(ineigh,iatom)
! Loop over neighPP-neighbors <phi_i|phi_j> of iatom
              do kneigh = 1,neighPPn(iatom)
                 katom = neighPP_j(kneigh,iatom)
                 kbeta = neighPP_b(kneigh,iatom)
! Test on iatom == katom
                 if (katom .eq. jatom .and. kbeta .eq. jbeta) then
                   nPPx_map(ineigh,iatom) = kneigh 
                   exit
                 endif ! if (katom)
              enddo ! do kneigh
           enddo ! do ineigh
        enddo ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   MAP nPP_map  (nPP -> neighPP)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The same like in the previuos case, but now we're going to map nPP list 
! into the slobal list neighPP
! Loop over atoms
!$omp parallel do private (jatom, jbeta, ineigh, kneigh, katom, kbeta) &
!$omp  private (vec1, vec2, distance)
        do iatom = 1, natoms
! Loop over nPPx-neighbors <Psi_i|phi_j> of iatom
           do ineigh = 1,nPPn(iatom)
              jatom = nPP_j(ineigh,iatom)
              jbeta = nPP_b(ineigh,iatom)
! Loop over neighPP-neighbors <phi_i|phi_j> of iatom
              do kneigh = 1,neighPPn(iatom)
                 katom = neighPP_j(kneigh,iatom)
                 kbeta = neighPP_b(kneigh,iatom)
! Test on iatom == katom
                 if(katom .eq. jatom .and. kbeta .eq. jbeta) then
                   nPP_map(ineigh,iatom) = kneigh 
                   exit
                 endif ! if (katom)
              enddo ! do kneigh
           enddo ! do ineigh
        enddo ! do iatom

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   SELF nPP_self  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! What's the neighbor of the atom itself?
        nPP_self = -999
!$omp parallel do private (ineigh, mbeta, jatom)
        do iatom = 1, natoms
          do ineigh = 1, nPPn(iatom)
            mbeta = nPP_b(ineigh,iatom)
            jatom = nPP_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0) nPP_self(iatom) = ineigh
          end do
        end do

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   SELF nPPx_self  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        nPPx_self = -999
        do iatom = 1, natoms
!$omp parallel do private (ineigh, mbeta, jatom)
          do ineigh = 1, nPPxn(iatom)
            mbeta = nPPx_b(ineigh,iatom)
            jatom = nPPx_j(ineigh,iatom)
            if (iatom .eq. jatom .and. mbeta .eq. 0) nPPx_self(iatom) = ineigh
          end do
        end do

!        if (iordern .eq. 1)                                                  &
!     &   call neighbors_ordern_final (natoms, nprocs, my_proc, ivdw)

        call writeout_neighborsPP (nprocs, iwrtneigh, basisfile)
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine neighborsPP
