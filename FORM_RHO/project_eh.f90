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
        subroutine project_eh(ioccupy_k, foccupy, ikpoint)

        use dimensions
        use constants_fireball
        use kpoints
        use charges
        use scf  
        use configuration
        
        use density
      
        use interactions
        use neighbor_map
        

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
     
! Output
     integer, intent(inout), dimension (norbitals, nkpoints) :: ioccupy_k  
     real, intent(inout), dimension (norbitals, nkpoints) :: foccupy
     integer, intent (in) :: ikpoint
 
! Local Parameters
! ===========================================================================

 
! Local Variable Declaration and Description
! ===========================================================================
        integer iband
        integer imu
        integer inu
        integer ielec

        real pmax
        real dqorb
        real dqelec
        real dqhole
        real, dimension (norbitals) :: projel
        real, dimension (norbitals) :: projh

! Procedure
! ===========================================================================
! project electron and hole on actual wf
      projel(:) = 0.0d0
      projh(:) = 0.0d0
      do iband = 1, norbitals
        do inu = 1, norbitals
         projh(iband) = projh(iband) + (wf_hole(inu, ikpoint)*blowre(inu,iband,1))  
         projel(iband) = projel(iband) + (wf_elec(inu,ikpoint)*blowre(inu,iband,1))
        end do
      end do         
! inquire state occupations of electron  
      pmax = 0.0d0
      do imu = 1, norbitals
       if (abs(projel(imu)) .gt. pmax) then  
          pmax = abs(projel(imu)) 
          id_elec = imu
       end if   
      end do 

! inquire state occupations of hole
      pmax = 0.0d0
      do imu = 1, norbitals
       if (abs(projh(imu)) .gt. pmax) then  
          pmax = abs(projh(imu)) 
          id_hole = imu
       end if   
      end do 

! several checks need to be done before we switch electron
! 1. state where electron will be placed is higher than state from which is taken
      if (id_hole .gt. id_elec) then
        write (*,*) "hole higher than electron, no switch"
        return
      endif 

      dqhole = foccupy(id_hole,ikpoint) - occup_elec
      dqelec = foccupy(id_elec,ikpoint) + occup_elec  
! 2. hole and excited electron states are crossing; smearing effect
      if ((dqhole .lt. 0.0d0) .and. (dqelec .gt. 1.0d0)) then
       dqorb = dqelec - 1.0d0
       if (abs(dqhole) .gt. (dqelec-1.0d0)) dqorb = dqhole
       foccupy(id_hole,ikpoint) = foccupy(id_hole,ikpoint) - occup_elec + dqorb
       foccupy(id_elec,ikpoint) = foccupy(id_elec,ikpoint) + occup_elec - dqorb
       if (foccupy(id_hole,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_hole,ikpoint) = 1
       else
         ioccupy_k(id_hole,ikpoint) = 0
       end if
       if (foccupy(id_elec,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_elec,ikpoint) = 1
       else
         ioccupy_k(id_elec,ikpoint) = 0
       end if

! dump it
       write (*,300) id_elec, projel(id_elec), foccupy(id_elec,ikpoint)
       write (*,301) id_hole, projh(id_hole), foccupy(id_hole,ikpoint)
       write (*,301) id_hole-1, projh(id_hole-1), foccupy(id_hole-1,ikpoint)
       write (111,'(<norbitals>f6.1)')  (foccupy(iband,ikpoint), iband=1,norbitals)

       return
      endif 

! 3. hole state becomes empty; we use second lower state to be empty
     if (dqhole .lt. 0.0d0) then 
       foccupy(id_hole,ikpoint) = foccupy(id_hole,ikpoint) - occup_elec - dqhole
       foccupy(id_hole-1,ikpoint) = foccupy(id_hole-1,ikpoint) + dqhole
       foccupy(id_elec,ikpoint) = foccupy(id_elec,ikpoint) + occup_elec 
       if (foccupy(id_hole,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_hole,ikpoint) = 1
       else
         ioccupy_k(id_hole,ikpoint) = 0
       end if
       if (foccupy(id_hole-1,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_hole-1,ikpoint) = 1
       else
         ioccupy_k(id_hole-1,ikpoint) = 0
       end if
       if (foccupy(id_elec,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_elec,ikpoint) = 1
       else
         ioccupy_k(id_elec,ikpoint) = 0
       end if

! dump it
       write (*,300) id_elec, projel(id_elec), foccupy(id_elec,ikpoint)
       write (*,301) id_hole, projh(id_hole), foccupy(id_hole,ikpoint) 
       write (*,301) id_hole-1, projh(id_hole-1), foccupy(id_hole-1,ikpoint) 
       write (111,'(<norbitals>f6.1)')  (foccupy(iband,ikpoint), iband=1,norbitals)

       return
     endif

! 4. electron state becomes filled; we 
     if (dqelec .gt. 1.0d0) then 
       dqorb = dqelec - 1.0d0
       foccupy(id_elec,ikpoint) = foccupy(id_elec, ikpoint) + occup_elec - dqorb
       foccupy(id_hole,ikpoint) = foccupy(id_hole,ikpoint) - occup_elec + dqorb
       if (foccupy(id_hole,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_hole,ikpoint) = 1
       else
         ioccupy_k(id_hole,ikpoint) = 0
       end if
       foccupy(id_hole-1,ikpoint) = foccupy(id_hole,ikpoint) + dqorb
       if (foccupy(id_hole-1,ikpoint) .gt. 1.0d-5) then
         ioccupy_k(id_hole-1,ikpoint) = 1
       else
         ioccupy_k(id_hole-1,ikpoint) = 0
       end if

! dump it
       write (*,300) id_elec, projel(id_elec), foccupy(id_elec,ikpoint)
       write (*,301) id_hole, projh(id_hole), foccupy(id_hole,ikpoint) 
       write (*,301) id_hole-1, projh(id_hole-1), foccupy(id_hole-1,ikpoint) 
       write (111,'(<norbitals>f6.1)')  (foccupy(iband,ikpoint), iband=1,norbitals)

       return

     endif

 
! remove an electron charge; i.e. form a hole     
      foccupy(id_hole,ikpoint) = foccupy(id_hole,ikpoint) - occup_elec
      if (foccupy(id_hole,ikpoint) .lt. 0.0d0) foccupy(id_hole,ikpoint) = 0.0d0
      if (foccupy(id_hole,ikpoint) .gt. 1.0d-5) then
        ioccupy_k(id_hole,ikpoint) = 1
      else
        ioccupy_k(id_hole,ikpoint) = 0
      end if
! add an electron charge; i.e. form an excited state     
      foccupy(id_elec,ikpoint) = foccupy(id_elec, ikpoint) + occup_elec
      if (foccupy(id_elec,ikpoint) .gt. 1.0d-5) then
        ioccupy_k(id_elec,ikpoint) = 1
      else
        ioccupy_k(id_elec,ikpoint) = 0
      endif

      write (*,300) id_elec, projel(id_elec), foccupy(id_elec,ikpoint)
      write (*,301) id_hole, projh(id_hole), foccupy(id_hole,ikpoint) 
      write (111,'(<norbitals>f6.1)')  (foccupy(iband,ikpoint), iband=1,norbitals)


200     format (' Band n = ', i4, ' k-points: ioccupy = ', i2)
201     format (' Band n = ', i4, ' foccupy = ', f6.3)
300     format (' Excited e- :', i4, f12.6, f6.2)
301     format (' Left h+    :', i4, f12.6, f6.2)

  end subroutine project_eh
