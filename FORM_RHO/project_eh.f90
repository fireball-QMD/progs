! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

 
! project_eh.f90
! Program Description
! ===========================================================================
!       This routine calculates the fermi energy.
!
! ===========================================================================
! Code written by:
! V. Zobac
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: zobac@fzu.cz
!
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
!        subroutine project_eh(ioccupy_k, foccupy, ikpoint)
        subroutine project_eh()


        use dimensions
        use kpoints
        use charges
        use scf  
        use configuration
        use MD
        use density    
        use interactions       
        use nonadiabatic ! vlada
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
     
! Output
 !    integer, intent(inout), dimension (norbitals, nkpoints) :: ioccupy_k  
 !    real, intent(inout), dimension (norbitals, nkpoints) :: foccupy
 !    integer, intent (in) :: ikpoint
 
! Local Parameters
! ===========================================================================

 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer jmu
        integer inu
        integer ielec
        integer jelec
        integer itime_step
        integer i 
        integer j 
        integer map_ks_fix(nele)
        real, dimension (norbitals) :: AAA
        integer nel
        integer save_map_ks
        integer loc_orb(nele)
        real pmax(nele)
        real norm_aux
        real, dimension (norbitals) :: map_proj_list
        real, dimension (norbitals) :: norm
        real, dimension (norbitals) :: save_wf

! Procedure
! ===========================================================================
!   itime_step = itime_step_g
!   flag_proj=1
    if (itime_step_g .ne. 0 ) then
!     write(33333,*) itime_step_g
!     write(11333,*) itime_step_g
!    write (171717,*) "ielec, loc_orb, map_ks(ielec),iocc(ielec)"
!    write(88888,'(i4,i4,<norbitals>f6.1)') itime_step, Kscf, (foccupy_na(j,1),j=1,norbitals)

! ==== w.f. projection ====
! loop over tracked electronic states
     do ielec=1,nele
! reset norm
      norm(:)=0.0d0

! loop over previous eigenstates
      do i=1,n_hist

        AAA(:)= wf_weight(i)*blowre_hist(:,ielec,i)
! loop over each eigenstate
        do imu=1,norbitals

! loop over dimension of the jmu-eigenstate
         norm_aux=0.0d0
          do inu=1,norbitals
            norm_aux = norm_aux + AAA(inu)*blowre(inu,imu,1)
          end do ! end do imu
          norm(imu) = norm(imu) + abs(norm_aux)
        end do  ! end do imu
      end do ! end do n_hist

! ===== estimate the best projection for each state ====
! which eigenstate has biggest overlap with the history?
      pmax(ielec) = 0.0d0
      do imu = 1, norbitals
        if (abs(norm(imu)) .gt. pmax(ielec)) then
          pmax(ielec) = abs(norm(imu))   
          loc_orb(ielec)=imu
        end if
      end do ! enddo imu
      
!      write(11333,*) ielec, map_ks(ielec), loc_orb(ielec),pmax(ielec)

     end do ! nele

! ==== check of degeneracy =====
     flag_proj = 0
! loop over electrons
     do ielec=1, nele

! check of multiple occupancies
       map_ks_fix(:) = 0
       do jelec= 1, nele
         if ((loc_orb(ielec) .eq. loc_orb(jelec)) .and. (ielec .ne. jelec)) then
            map_ks_fix(ielec) = map_ks_fix(ielec) + 1
            map_ks_fix(jelec) = 1
         endif
       end do ! do jelec
! no degeneracy:
       if (map_ks_fix(ielec) .eq. 0) then

! check the projection
!         if((pmax(ielec) .gt. 0.75) .and. (loc_orb(ielec) .eq. map_ks(ielec)) ) then
         if(loc_orb(ielec) .eq. map_ks(ielec)) then
           if(hist_fix(ielec) .eq. 0) then
             write(*,*) "save wf to hist: step, ielec", itime_step_g, ielec
! projection is good, so update the eigenstate history
             do i= 1,n_hist-1
               blowre_hist(:,ielec,i) = blowre_hist(:,ielec,i+1)
             end do !
             blowre_hist(:,ielec,n_hist) = blowre(:,loc_orb(ielec),1)
           else
             blowre_hist(:,ielec,n_hist)   = blowre(:,loc_orb(ielec),1)
             blowre_hist(:,ielec,n_hist-1) = blowre(:,loc_orb(ielec),1)
             blowre_hist(:,ielec,n_hist-2) = blowre(:,loc_orb(ielec),1)
             hist_fix(ielec) = 0  
           end if ! hist_fix(ielec)
         else
           hist_fix(ielec) = 1
         endif ! loc_orb(ielec) .eq. map_ks(ielec)

       else

! multiple degeneracy: stay with the previous time step and reset the reference projection
! reset degenerated states to previous time step
         loc_orb(ielec) = map_ks(ielec)
         flag_proj = 1
! set the ref. projection to given state
         do i= 1,n_hist
           blowre_hist(:,ielec,i) = blowre(:,loc_orb(ielec),1)
         end do !
         do jelec = 1, nele
! select the degenerate states
           if (map_ks_fix(jelec) .eq. 1) then
             loc_orb(jelec) = map_ks(jelec)
! set the ref. projection to the previous time step
             do i= 1,n_hist
              blowre_hist(:,jelec,i) = blowre(:,map_ks(jelec),1)
             end do !
!              blowre_hist(:,jelec,n_hist) = blowre(:,map_ks(jelec),1)
           endif !if map_ks
         end do ! do jelec
       end if ! if map_ks_fix

     end do ! do iele


! update map_ks from actual loc_orb (free of the degeneracy)
    map_ks = loc_orb
!write (1004,'(<norbitals>f10.4)') (foccupy_na(imu,1), imu = 1, norbitals)
    do ielec=1, nele
       if (iocc(ielec) .eq. 1.0) then
          foccupy_na(loc_orb(ielec),1) = iocc(ielec)*0.5d0
          ioccupy_na(loc_orb(ielec),1) = iocc(ielec)        
       end if
    end do ! do iele
!write (1005,'(<norbitals>f10.4)') (foccupy_na(imu,1), imu = 1, norbitals)

! ???
!     foccupy_na_o=foccupy_na 
     
! check degeneracy again
     do ielec=1, nele
       do jelec= 1, nele
         if ((loc_orb(ielec) .eq. loc_orb(jelec)) .and. (ielec .ne. jelec)) then
!         if (loc_orb(ielec) .eq. loc_orb(jelec)) then
          write(*,*) "WARNING_degeneracy", itime_step_g, loc_orb(ielec), loc_orb(jelec) 
         endif
       end do ! do jelec
!       write(33333,*) ielec, map_ks(ielec), loc_orb(ielec),pmax(ielec)
     end do ! ielec

     nel = 0
     do ielec=1, norbitals   
      if (foccupy_na(ielec,1) .eq. 0.5) then
         nel = nel +1
         loc_el(nel) = ielec 
      end if 
     end do

  end if ! if itime_step_g
  
200     format (' Band n = ', i4, ' k-points: ioccupy = ', i2)
201     format (' Band n = ', i4, ' foccupy = ', f6.3)

  end subroutine project_eh
