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


! common_neighbors.f90
! Program Description
! ===========================================================================
!       Given atom ialp, we find pairs of atoms that are common neighbors
! of atom ialp.
!
! Find all common neighbors of atom alpha (we call its i value ialp).
! atom alpha is in the central cell, and the common neighbors are
! at l-vector-sub-i,i and l-vector-sub-j,j
!
! Note:
! You should remember that in counting the common neighbors, we do not include
! any pairs in which one or both of the atoms is ontop atom ialp. Thus we are 
! keeping only third party common neighbors.
!
! ncomn(ialp) = num =total # of common neighbors to atom (0,ialp)
! ncomj(ialp,mneigh,1) = iatom value for mneigh'th pair of common neighbors
! ncomj(ialp,mneigh,2) = jatom value for mneigh'th pair of common neighbors
! ncomb(ialp,mneigh,1) = mbeta-sub-iatom value for m'th pair of common neighbor
! ncomb(ialp,mneigh,2) = mbeta-sub-jatom value for m'th pair of common neighbor 
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine common_neighbors (nprocs, my_proc, iordern, iwrtneigh)
        use configuration
        use dimensions
        use interactions
        use neighbor_map
        use options, only : itheory_xc 
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: iordern
        integer, intent (in) :: iwrtneigh
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs

 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer ibeta
        integer imu
        integer in1, in2
        integer ineigh
        integer jatom
        integer jbeta
        integer jneigh
        integer katom
        integer kbeta
        integer kneigh
        integer num_neigh
        integer iatomstart, natomsp

        real distance
        real distance2
        real range2
        real rcutoff_i, rcutoff_j

        real, dimension (3) :: diff
        real, dimension (3) :: dvec
        real, dimension (3) :: vec, vec1, vec2

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
!     if (my_proc .eq. 0 .and. wrtout) then
!      write (*,*) ' Welcome to common_neighbors - determine mapping of '
!      write (*,*) ' common neighbors. '
!     end if

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

! Loop over all atoms
!$omp parallel do private (num_neigh, imu, iatom, jatom, katom, ibeta, jbeta) &
!$omp&   private (ineigh, jneigh, kbeta, in1, in2, rcutoff_i, rcutoff_j, vec) &
!$omp&   private (vec1, vec2, distance2, range2, diff)
      do ialp = iatomstart, iatomstart - 1 + natomsp

! The variable num_neigh counts neighbors.
         num_neigh = 0

! Loop over all known neighbors of ialp. Call these atoms ineigh.
         do ineigh = 1, neighn(ialp)
          iatom = neigh_j(ineigh,ialp)
          ibeta = neigh_b(ineigh,ialp)

! Keep only third party common neighbors.
          if (.not. (iatom .eq. ialp .and. ibeta .eq. 0)) then
           in1 = imass(iatom)
           rcutoff_i = 0.0d0
           do imu = 1, nssh(in1)
            if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
           end do
           vec1(:) = ratom(:,iatom) + xl(:,ibeta)

! Loop over all known neighbors of atom ialp again, and check to see if atoms 
! iatom and jatom are neighbors. If so, then we have found a pair of common 
! neighbors. Of course, if iatom is jatom in all respects, that does not count.
           do jneigh = 1, neighn(ialp)
            jatom = neigh_j(jneigh,ialp)
            jbeta = neigh_b(jneigh,ialp)
!############# BEGIN MODIFICATION: SYMMETRY: NEW (APRIL 2018)
! Keep only third party common neighbors.
            if (.not. (jatom .eq. ialp .and. jbeta .eq. 0)                   &
!           &          .and. (ineigh .ne. jneigh)) then
            &          .and. (ineigh .lt. jneigh)) then    
!Modification; instead of demanding that iatom and jatom are not the
!same atoms, we demandn that iatom's index is strictly less than jatom's. That way
!we cut in half the number of pairs of common neighbors
!############ END MODIFICATION: SYMMETRY: NEW (APRIL 2018)
             in2 = imass(jatom)
             rcutoff_j = 0.0d0
             do imu = 1, nssh(in2)
              if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
             end do
             vec2(:) = ratom(:,jatom) + xl(:,jbeta)

! Check distance from iatom to jatom.
             distance2 = (vec2(1) - vec1(1))**2 + (vec2(2) - vec1(2))**2     &
                        + (vec2(3) - vec1(3))**2
 
             range2 = (rcutoff_i + rcutoff_j - 0.01d0)**2
             if (distance2 .le. range2) then

! Found one!
              num_neigh = num_neigh + 1

! Check that # of common neighbor pairs is less than dimensions.
              if (num_neigh .gt. neigh_max**2) then
               write (*,*) ' Oh no. In common_neighbors.f90, we have too many '
               write (*,*) ' common neighbors. In MODULES/dimensions.f90 '
               write (*,*) ' dimension neigh_max**2 = ', neigh_max**2
               write (*,*) ' So far (*but still counting) we have '
               write (*,*) ' num_neigh = ', num_neigh, ' within neighbors.f90! '
! FIXME can't have a stop statement in an OpenMP loop
!              stop
              end if

! Put iatom in first spot
! Put jatom in second spot
              neigh_comj(1,num_neigh,ialp) = iatom
              neigh_comb(1,num_neigh,ialp) = ibeta
              neigh_comj(2,num_neigh,ialp) = jatom
              neigh_comb(2,num_neigh,ialp) = jbeta
              if (itheory_xc .eq. 4) then
              neigh_com_ng(1,num_neigh,ialp) = ineigh
              neigh_com_ng(2,num_neigh,ialp) = jneigh 
              end if

! We also need to know for a given ialp and (iatom,ibeta), what is the m value
! for (jatom,jbeta) with respect to iatom. That is, jatom is the m'th neighbor
! of iatom. What is m?

! Set to a crazy value, in case loop fails.
              diff = vec2 - vec1
              neigh_comm(num_neigh,ialp) = -9999
              do kneigh = 1, neighn(iatom)
               katom = neigh_j(kneigh,iatom)
               kbeta = neigh_b(kneigh,iatom)
               vec(:) = xl(:,kbeta) + ratom(:,katom) - ratom(:,iatom)
               if ((abs(vec(1) - diff(1)) .lt. 1.0d-4) .and.                 & 
     &             (abs(vec(2) - diff(2)) .lt. 1.0d-4) .and.                 &
     &             (abs(vec(3) - diff(3)) .lt. 1.0d-4))                      &
     &          neigh_comm(num_neigh,ialp) = kneigh
              end do

! We check to see if it really was a neighbor. (could be very close to cutoff)
              if (neigh_comm(num_neigh,ialp) .eq. -9999)                     &
     &         num_neigh = num_neigh - 1 
             end if
            end if
           end do
          end if
         end do
         neigh_comn(ialp) = num_neigh
        end do

        if (iordern .eq. 1)                                                  &
     &   call common_neighbors_ordern_final (natoms, nprocs, my_proc)

! Option for writing out the results for debugging.
        if (iwrtneigh .eq. 1 .and. my_proc.eq.0) then
         write (*,*) '  '
         write (*,*) ' Common neighbors of each atom: '
         do ialp = 1, natoms
          num_neigh = neigh_comn(ialp)
          write (*,*) '  '
          write (*,*) ' num_neigh = ', num_neigh
          write (*,100)
          write (*,101)
          do ineigh = 1, num_neigh
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           kneigh = neigh_comm(ineigh,ialp)

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
        end
