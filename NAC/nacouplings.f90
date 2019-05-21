! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad Autonoma de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold

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


! nacouplings.f90
! Program Description
! ===========================================================================
!       This routine calculates the NAC between selected Kohn-Sham eigenstates
! < Psi_i | d/dR Psi_j >
!
! ===========================================================================
! Code written by:
! Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!      
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine nacouplings ()
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use nonadiabatic
        use options
        !use qmmm_module, only : qmmm_struct
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
!       integer, intent (in) :: icluster


! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ix
        integer iatom
        integer iband
        integer jband
        integer ikpoint
        integer imu, inu
        integer ineigh
        integer in1, in2
        integer iorbital
        integer issh
        integer jatom
        integer jneigh
        integer mqn
        integer mbeta
        integer mmu
!       integer noccupy
        integer nnu
        integer katom

!       integer, dimension (norbitals) :: ioccupy
!       integer, dimension (norbitals, nkpoints) :: ioccupy_k

        real aux1, aux2, aux3
        real dot
        real gutr
        real pcharge
        real ztest
        real cmunu
        real tolnac
        real diff
! JOM-test
        real, dimension (3,natoms) :: ftest

!       real, dimension (norbitals, nkpoints) :: foccupy
!       real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec

        complex ai
        complex phase, phasex
        complex step1, step2

!       logical read_occupy

! Procedure
! ===========================================================================
! Initialize some things
        ai = cmplx(0.0d0,1.0d0)
        gks = 0.0d0
! tolerance for degeneracy in eigenstates for nonadiabatic coupling
        tolnac = 0.0001d0

! JOM-test : set some values to 0.0d0 for testinf purposes
!      write(*,*)'TESTING nacouplings, result not correct!!!!!'
!       gover = 0.0d0
!       gh_2c = 0.0d0
!       gh_3c = 0.0d0
!       gh_xc_3c = 0.0d0
!       gh_atm = 0.0d0
!       gh_pp_atm = 0.0d0
!       gh_pp_otl = 0.0d0
!       gh_pp_otr = 0.0d0
!       gh_pp_3c = 0.0d0
! JOM-test

! JOM I have separated gh_3c into gh_3c and gh_xc_3c for testing. Sum
! them here now
     !   gh_3c = gh_3c + gh_xc_3c
         gh_3c = gh_3c + gh_lrew_qmmm
! itheory = 1 : sum contributions
!        if (itheory .eq. 1) then
       !  gh_3c = gh_3c + gh_3c_ca + gh_lrew
       !  gh_2c = gh_2c + gh_2c_ca
       !  gh_atm = gh_atm + gh_atm_ca
!        end if

! JOM so far only icluster.eq.1 works
        if (icluster .ne. 1) then
         write(*,*)'icluster .ne. 1 in nacouplings; must stop'
         stop
        end if  

! JOM : loops on the selected KS-bands
        do iband = 1, nele
         do jband = 1, nele
! Loop over the special k points.
          do ikpoint = 1, nkpoints
           if (iband .ne. jband) then
           
! JOM check for degeneracy
             diff = abs(eigen_k(map_ks(iband),ikpoint) - eigen_k(map_ks(jband),ikpoint) )
             if (diff .lt. tolnac) then
             !  write(*,*)'TWO EIGENVALUES VERY CLOSE'
             !  write(*,*)'band', iband, eigen_k(map_ks(iband),ikpoint)
             !  write(*,*)'band', jband, eigen_k(map_ks(jband),ikpoint)
             !  write(*,*)'The nonadiabatic coupling is'   
             !  write(*,*)'NOT CALCULATED'   
             else
! LOOP-2C
! ****************************************************************************
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
             do iatom = 1, natoms
               in1 = imass(iatom)
               do ineigh = 1, neighn(iatom)
                mbeta = neigh_b(ineigh,iatom)
                jatom = neigh_j(ineigh,iatom)
                in2 = imass(jatom)
                vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

                dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
!         phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
                phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
! JOM : I guess we do not need now any phase (icluster.eq.1) but we may
! need it later; keep it just in case, buy without foccupy
!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                phase = phasex
                do imu = 1, num_orb(in1)
                 mmu = imu + degelec(iatom)
                 step1 = phase*bbnkre(mmu,map_ks(iband),ikpoint)
                 do inu = 1, num_orb(in2)
                   nnu = inu + degelec(jatom)
!               step2 = step1*bbnkre(nnu,map_ks(iband),ikpoint)
                   step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
! JOM : careful with this once we include periodicity
                   gutr = real(step2)
! Finally the expressions.........
                   cmunu = gutr
                   gks(:,iatom,iband,jband) = gks(:,iatom,iband,jband) +    &
     &   cmunu*( gover(:,imu,inu,ineigh,iatom)*eigen_k(map_ks(jband),ikpoint)    &
     &           - gh_2c(:,imu,inu,ineigh,iatom) )
  
                  gks(:,jatom,iband,jband) = gks(:,jatom,iband,jband) +    &
     &   cmunu*( - gover(:,imu,inu,ineigh,iatom)*eigen_k(map_ks(iband),ikpoint)  &
     &           + gh_2c(:,imu,inu,ineigh,iatom) )
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop to add 3-C contributions
                  do katom = 1, natoms
                    gks(:,katom,iband,jband) = gks(:,katom,iband,jband) -   &
     &             cmunu*gh_3c(:,katom,imu,inu,ineigh,iatom)            
                  end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 end do ! do inu
                end do ! do imu

! JOM : special case: atom-case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATOM case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                do imu = 1, num_orb(in1)
                 mmu = imu + degelec(iatom)
!              step1 = bbnkre(mmu,map_ks(iband),ikpoint)*spin
                 step1 = bbnkre(mmu,map_ks(iband),ikpoint)
                 do inu = 1, num_orb(in1)
                  nnu = inu + degelec(iatom)
                  step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
! JOM : careful with this once we include periodicity
                  gutr = real(step2)
! Finally the expressions.........
                  cmunu = gutr
                  gks(:,iatom,iband,jband) = gks(:,iatom,iband,jband) +    &
     &   cmunu*( - gh_atm(:,imu,inu,ineigh,iatom) )
!
                  gks(:,jatom,iband,jband) = gks(:,jatom,iband,jband) +    &
     &   cmunu*(  gh_atm(:,imu,inu,ineigh,iatom) )
                 end do ! do inu
                end do ! do imu

! JOM-test
1001    continue
! Finish loop over atoms and neighbors.
               end do ! do ineigh
              end do ! do iatoms

! JOM-test
!         go to 1002
! LOOPS-2C-PP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          PP-neighbors-2C
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over all atoms iatom in the unit cell
              do iatom = 1, natoms
               in1 = imass(iatom)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Ontop-Left case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over the neighbors of each iatom for ontop-Left case
               do ineigh = 1, nPPxn(iatom)
                mbeta = nPPx_b(ineigh,iatom)
                jatom = nPPx_j(ineigh,iatom)
                in2 = imass(jatom)

                vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
                dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
!         phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
                phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)

!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                phase = phasex
                do imu = 1, num_orb(in1)
                  mmu = imu + degelec(iatom)
                  step1 = phase*bbnkre(mmu,map_ks(iband),ikpoint)
                  do inu = 1, num_orb(in2)
                    nnu = inu + degelec(jatom)
                    step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
                    gutr = real(step2)
                    cmunu = gutr
! Finally the expressions.........
                    gks(:,iatom,iband,jband) = gks(:,iatom,iband,jband) +    &
     &   cmunu*( - gh_pp_otl(:,imu,inu,ineigh,iatom) )                   
!
                    gks(:,jatom,iband,jband) = gks(:,jatom,iband,jband) +    &
     &   cmunu*(  gh_pp_otl(:,imu,inu,ineigh,iatom) )                   
!
                  end do ! do inu
                end do  ! do imu
! Finish loop over neighbors.
               end do ! do ineigh
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Ontop-Right case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over the neighbors of each iatom for ontop-Left case
               do ineigh = 1, nPPn(iatom)
                mbeta = nPP_b(ineigh,iatom)
                jatom = nPP_j(ineigh,iatom)
                in2 = imass(jatom)

                vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
                dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
!         phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
                phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)

!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                phase = phasex
                do imu = 1, num_orb(in1)
                  mmu = imu + degelec(iatom)
                  step1 = phase*bbnkre(mmu,map_ks(iband),ikpoint)
                  do inu = 1, num_orb(in2)
                    nnu = inu + degelec(jatom)
                    step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
                    gutr = real(step2)
                    cmunu = gutr
! Finally the expressions.........
                    gks(:,iatom,iband,jband) = gks(:,iatom,iband,jband) +    &
     &   cmunu*( - gh_pp_otr(:,imu,inu,ineigh,iatom) )                   
!
                    gks(:,jatom,iband,jband) = gks(:,jatom,iband,jband) +    &
     &   cmunu*(  gh_pp_otr(:,imu,inu,ineigh,iatom) )                   
!
                  end do ! do inu
                end do  ! do imu
! Finish loop over neighbors.
               end do  ! do ineigh
! JOM : special case: atom-case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATOM case
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over the neighbors of each iatom for atom case
               do ineigh = 1, nPPn(iatom)
                 mbeta = nPP_b(ineigh,iatom)
                 jatom = nPP_j(ineigh,iatom)
                 in2 = imass(jatom)

                 vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
                 dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
!         phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
                 phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)

!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                 phase = phasex

                 do imu = 1, num_orb(in1)
                   mmu = imu + degelec(iatom)
!              step1 = bbnkre(mmu,map_ks(iband),ikpoint)*spin
                   step1 = bbnkre(mmu,map_ks(iband),ikpoint)
                   do inu = 1, num_orb(in1)
                     nnu = inu + degelec(iatom)
                     step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
! JOM : careful with this once we include periodicity
                     gutr = real(step2)
! Finally the expressions.........
                     cmunu = gutr
                     gks(:,iatom,iband,jband) = gks(:,iatom,iband,jband) +    &
     &   cmunu*( - gh_pp_atm(:,imu,inu,ineigh,iatom) )
!
                     gks(:,jatom,iband,jband) = gks(:,jatom,iband,jband) +    &
     &   cmunu*(  gh_pp_atm(:,imu,inu,ineigh,iatom) )
                
                   end do  ! do inu
                 end do  ! do imu

! Finish loop over atoms and neighbors.
               end do ! do ineigh
             end do  ! do iatom
! LOOPS-3C-PP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          PP-neighbors-3C
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over all atoms iatom in the unit cell
             do iatom = 1, natoms
              in1 = imass(iatom)
              do ineigh = 1, neighPPn(iatom)
                mbeta = neighPP_b(ineigh,iatom)
                jatom = neighPP_j(ineigh,iatom)
                in2 = imass(jatom)
                vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
                dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
!         phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
                phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
!              phase = phasex*foccupy(map_ks(iband),ikpoint)
                phase = phasex
                do imu = 1, num_orb(in1)
                  mmu = imu + degelec(iatom)
                  step1 = phase*bbnkre(mmu,map_ks(iband),ikpoint)
                  do inu = 1, num_orb(in2)
                    nnu = inu + degelec(jatom)
                    step2 = step1*bbnkre(nnu,map_ks(jband),ikpoint)
                    gutr = real(step2)
                    cmunu = gutr
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop to add 3-C-PP contributions
                    do katom = 1, natoms
                      gks(:,katom,iband,jband) = gks(:,katom,iband,jband) -   &
     &           cmunu*gh_pp_3c(:,katom,imu,inu,ineigh,iatom)            
                    end do ! do katom
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  end do  ! do inu
                end do  ! do imu
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              end do ! do ineigh
             end do ! do iatom
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             do katom = 1, natoms
               gks(:,katom,iband,jband) = gks(:,katom,iband,jband) /                 &
     &    (eigen_k(map_ks(iband),ikpoint) - eigen_k(map_ks(jband),ikpoint) )
             end do  ! do katom
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
! Finish if(diff .lt. tolnac))
            end if
           end if
! Finish loop over k-points.
          end do ! do kpoints
! Finish loop over iband and jband.
         end do ! do jband
        end do ! do iband
        
! dump dij 


     do iband = 1, nele
      do jband = 1, nele
  !      iband = 1
  !      jband = 2
        ikpoint = 1
        ! write (210,*)
         do katom = 1, natoms
            diff = eigen_k(map_ks(iband),ikpoint) - eigen_k(map_ks(jband),ikpoint)
      if (iband .gt. jband) then
            write(210,100) katom, iband, jband, diff, gks(:,katom,iband,jband),   &
         & sqrt(gks(1,katom,iband,jband)**2 + gks(2,katom,iband,jband)**2 + gks(3,katom,iband,jband)**2)  ! ENRIQUE-JOM
      end if 
         end do ! do katom
  !          write(210,*) "----------"
      end do
     end do  
! Provisionally, deallocate bbnkre, blowre, eigen_k here
!       deallocate (eigen_k)
!       deallocate (bbnkre)
!       deallocate (blowre)
! JOM-info icluster = 1
!       deallocate (bbnkim) 
!       deallocate (blowim)



! Format Statements
! ===========================================================================
100     format ( 3i5, 5f12.4)
101     format (2x, 4(2x,f11.5))
200     format ('gover',4i4, f8.4 )
201     format ('gh_2c',4i4, f8.4 )
202     format ('gh_atm',4i4, f8.4 )
300     format ('gh_2c',4i4, f8.4 )
301     format (2x, i4, f10.6)
800     format (2x,4i3,f12.6)
        return
      end subroutine nacouplings

