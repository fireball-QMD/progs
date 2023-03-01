! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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

 
! writeout_charges.f90
! Program Description
! ===========================================================================
!       Write out the charges for restart capabilities.
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
        subroutine writeout_charges (natoms, ifixcharge, iqout, iwrtcharges, &
     &                               iwrtdensity, basisfile, symbol,ab)
        use charges
        use outputs, only : iwrtdipole
        use density
        use interactions
        use neighbor_map
        use scf, only : scf_achieved, Kscf
        use MD, only : itime_step_g
        use configuration, only : ratom, xl
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrtcharges
        integer, intent (in) :: iwrtdensity
        integer, intent (in) :: ab !ab = 0: before mixing. ab=1: after mixing

        character (len=40) basisfile
        character (len=2), dimension (natoms) :: symbol

 
! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter ::  Debye = 0.208194
        real, parameter ::  klambda = 1.0
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer in1, in2
        integer ineigh
        integer inu
        integer issh
        integer jatom
        integer mbeta
        real    Qtot, Qtot1, Qtot2
        real    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x, dip_res_x, dipQin_x
        real    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y, dip_res_y, dipQin_y
        real    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z, dip_res_z, dipQin_z
        real    dip_tot, dip_proy, dipQin_tot, dipTot_tot, dipIntra_tot, dipQout_tot, dip_res_tot
        real, dimension(3) :: r1,r2,Rbc,u21

        real, dimension(3) :: rmedio, raux
        real    w_suma
        real, dimension(3,3) :: bwrr, bwrr_inv, u_bwrr, ut_bwrr, v_bwrr, vt_bwrr, zero_bwrr
        real, dimension (natoms) :: c_k, dq_DP
        real, dimension (neigh_max) :: w_k
        real, dimension(3,natoms) :: intra_dip, res_dip
        real, dimension(3,1) :: intra_dip_aux, delta_ck
        integer n_bwrr, lda_bwrr, info, lda, i
        integer, dimension(3) :: ipiv
        real, dimension(3) :: lwork, s_bwrr, dummy
 


! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! ****************************************************************************
!
write(*,* ) 'MDC writeout_charges'
! W R I T E    O U T    D E N S I T Y  
! ****************************************************************************
! Write out the density for the first atom.
        if (iwrtdensity .eq. 1) then
         do iatom = 1, natoms
          write (*,*) '  '
          write (*,*) ' Write out the density matrix for iatom = ', iatom
          write (*,*) ' Number of neighbors = ', neighn(iatom)
          in1 = imass(iatom)
          write (*,*) ' Number of orbitals on iatom, num_orb(in1) = ',       &
     &     num_orb(in1)
          do ineigh = 1, neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           write (*,*) '  '
           write (*,*) ' iatom, jatom = ', iatom, jatom
           write (*,*) ' Number of orbitals on jatom, num_orb(in2) = ',      &
     &      num_orb(in2)
           write (*,*) ' ----------------------------------------------'
           do imu = 1, num_orb(in1)
            write (*,400) (rho(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          end do
         end do         
        end if ! iwrtdensity



        if (iwrtcharges .eq. 1) then
! ****************************************************************************
!
! W R I T E    O U T    L O W D I N    C H A R G E S
! ****************************************************************************
         if (iqout .eq. 1 .or. iqout .eq. 3) then
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' LOWDIN CHARGES (by shell): '
          write (*,500)
          write (*,501)
          write (*,500)
          do iatom = 1, natoms
           in1 = imass(iatom)
           write (*,502) iatom, symbol(iatom), nssh(in1),                   &
     &                   (Qin(issh,iatom), issh = 1, nssh(in1))
          end do

          write (*,500)
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' Total Lowdin charges for each atom: '
          write (*,500)
          do iatom = 1, natoms
           write (*,503) iatom, symbol(iatom), QLowdin_TOT(iatom)
          end do
          write (*,500)
          write (*,*) '  '
          write (*,*) '  '
         end if ! end if (iqout .ne. 2)


! ****************************************************************************
!
! W R I T E    O U T    M U L L I K E N    C H A R G E S
! ****************************************************************************
         if ((iqout .eq. 2) .or. (iqout .eq. 4)) then
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' MULLIKEN CHARGES (by shell): '
          write (*,500)
          write (*,501)
          write (*,500)
          do iatom = 1, natoms
           in1 = imass(iatom)
           write (*,502) iatom, symbol(iatom), nssh(in1),                   &
     &                   (Qin(issh,iatom), issh = 1, nssh(in1))
          end do

          write (*,500)
          write (*,*) '  '
          write (*,*) '  '
          write (*,*) ' Total Mulliken charges for each atom: '
          write (*,500)
          do iatom = 1, natoms
           write (*,503) iatom, symbol(iatom), QMulliken_TOT(iatom)
          end do
          write (*,500)
          write (*,*) '  '
          write (*,*) '  '
         end if !end if (iqout .eq. 2)

         write (*,500)
         write (*,*) '  '
        end if  !end if (iwrtcharges .eq. 1)

! ****************************************************************************
!
! W R I T E    O U T    T E M P O R A L    S E R I E S   O F    C H A R G E S
! ****************************************************************************
      if (iwrtcharges .eq. 2) then
      if ( scf_achieved ) then

      open(unit = 333, file = 'CHARGES.xyz', status = 'unknown', &
                                      & position = 'append')

       
      open(unit = 334, file = 'CHARGES_series', status = 'unknown', &
                                      & position = 'append')

       !  write (333,*) '  '
       !  write (333,*) '  '
       !  write (334,*) '  '
       !  write (334,*) '  '


      write(333,*) natoms
      write(333,*) 'Time step = ', itime_step_g
      do iatom = 1, natoms
       Qtot=0
       in1 = imass(iatom)
       do issh = 1,nssh(in1)
          Qtot = Qtot+Qin(issh,iatom)
       end do
     !      write (333,*) symbol(iatom),                                 &
     ! &                   ratom(1,iatom), ratom(2,iatom), ratom(3,iatom), &
          write(333,601)              (Qin(issh,iatom), issh = 1, nssh(in1)),         &
     &                   Qtot
           write(334,400,advance="no") Qtot
          
      
      end do !end do iatom = 1,natoms

        close(334)
         close(333)

       end if !end if ( scf_achieved )
      end if !end if (iwrtcharges .eq. 2)


!***************
                !IWRTCHARGES = 3
!***************
         if (iwrtcharges .eq. 3) then
    
      if (ab .eq. 1) then

      open(unit = 333, file = 'CHARGES.xyz', status = 'unknown', &
                                      & position = 'append')
 
      else ! ab .eq. 0
 
      open(unit = 333, file = 'CHARGES_no_mix.xyz', status = 'unknown', &
                                      & position = 'append')

      end if ! if ab .eq. 1

      !open(unit = 334, file = 'CHARGES_series', status = 'unknown', &
      !                                & position = 'append')


      write(333,*) natoms
      write(333,*) 'Time step = ', itime_step_g, 'Kscf = ', Kscf
      do iatom = 1, natoms
       Qtot=0
       in1 = imass(iatom)
       do issh = 1,nssh(in1)
          Qtot = Qtot+Qout(issh,iatom)
       end do
      !     write (333,*) symbol(iatom),&
     ! &                   ratom(1,iatom), ratom(2,iatom), ratom(3,iatom), &
      write(333,601)                  (Qout(issh,iatom), issh = 1, nssh(in1)),  &
     &                   Qtot
          ! write(334,400,advance="no") Qtot


      end do !end do iatom = 1,natoms

      !close(334)
      close(333)     
      
      end if !end if (iwrtcharges .eq. 3)

    

! ****************************************************************************
!
! W R I T E    O U T    C H A R G E S    F I L E
! ****************************************************************************
! Open the file CHARGES which contain the Lowdin charges for restart purposes.
        if (ifixcharge .eq. 0 .and. wrtout) then
         open (unit = 21, file = 'CHARGES', status = 'unknown')
         write (21,600) natoms, basisfile, iqout

! Write the output charges to the CHARGES file.
         do iatom = 1, natoms
          in1 = imass(iatom)
          write (21,601) (Qin(issh,iatom), issh = 1, nssh(in1))
         end do
         close (unit = 21)
        end if


! ****************************************************************************
! W R I T E   T H E   D I P O L E   (For testing purposes)
! ****************************************************************************
     if (iwrtdipole .gt. 0) then
      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0
      do iatom = 1, natoms
         Qtot=-Q0_TOT(iatom)
         in1 = imass(iatom)
         do issh = 1,nssh(in1)
             Qtot = Qtot+Qin(issh,iatom)
         end do
            
            !write(*,*) 'ratom is ',ratom(1,iatom),ratom(2,iatom),ratom(3,iatom)
            dip_x = dip_x+Qtot*ratom(1,iatom)
            dip_y = dip_y+Qtot*ratom(2,iatom)
            dip_z = dip_z+Qtot*ratom(3,iatom)

       enddo !end do iatom = 1,natoms
 
        dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )   

         dipQin_x=dip_x
         dipQin_y=dip_y
         dipQin_z=dip_z
         dipQin_tot=dip_tot

     open( unit = 1729, file = 'dipole_Qin', status = 'unknown',  &
                                   &     position = 'append')
       
       write(1729,444) dip_x/Debye, dip_y/Debye, dip_z/Debye, dip_tot/Debye

     close(1729)

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0
      
      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         do issh = 1,nssh(in1)
             Qtot = Qtot+Qin(issh,iatom)
         end do !end do issh = 1,nssh(in1)
         dip_x = dip_x-Q0_TOT(iatom)*r1(1)
         dip_y = dip_y-Q0_TOT(iatom)*r1(2)         
         dip_z = dip_z-Q0_TOT(iatom)*r1(3)
      end do !end do iatom = 1,natoms

       !dip_x = 0.0d0
       !dip_y = 0.0d0
       !dip_z = 0.0d0


      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)         
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)
                 
                !if (iatom .ne. jatom) then !.or. imu .eq. inu) then
                !else 
                  dip_x = dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y = dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z = dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))


                 !dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1)+ &
                 !           & dipc(2,imu,inu,ineigh,iatom)*u21(2)+ &
                 !           & dipc(3,imu,inu,ineigh,iatom)*u21(3)

                  !dip_x =  dip_x+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(1)+ &
                  !       & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  !dip_y =  dip_y+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(2)+ &
                  !       & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  !dip_z =  dip_z+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(3)+ &
                  !       & Rbc(3)*s_mat(imu,inu,ineigh,iatom))

                !end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      end do ! end do iatom = 1,natoms
     dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )

         dipTot_x=dip_x
         dipTot_y=dip_y
         dipTot_z=dip_z
         dipTot_tot=dip_tot

     open( unit = 666, file = 'dipole_Tot', status = 'unknown',  &
                                   &     position = 'append')

     write(666,444) dip_x/Debye, dip_y/Debye, dip_z/Debye, dip_tot/Debye
     close(666)

     !W R I T E     T H E     B I G     F I L E

               !write the intraatomic dipole


       dip_x = 0.0d0
       dip_y = 0.0d0
       dip_z = 0.0d0


      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)

                if ((iatom .eq. jatom) .and. (imu .ne. inu)) then
                !else 
                  dip_x = dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y = dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z = dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))


                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      end do ! end do iatom = 1,natoms
           dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )

         dipIntra_x=dip_x
         dipIntra_y=dip_y
         dipIntra_z=dip_z
         dipIntra_tot=dip_tot
     

!           if ( iwrtdipole .eq.2 .and. ab .eq. 0 ) then
           if ( iwrtdipole .eq.2) then
           
                       !DIP RES = DIP_TOT - PROY
      dip_res_x = 0.0d0
      dip_res_y = 0.0d0
      dip_res_z = 0.0d0
      
      do iatom = 1, natoms
         dip_res_x = 0.0d0
         dip_res_y = 0.0d0
         dip_res_z = 0.0d0
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)         
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)
                  if (iatom .ne. jatom) then !.or. imu .eq. inu) then
                   u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))

                    dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1)+ &
                            & dipc(2,imu,inu,ineigh,iatom)*u21(2)+ &
                            & dipc(3,imu,inu,ineigh,iatom)*u21(3)

                    dip_res_x = dip_res_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(1))

                    dip_res_y = dip_res_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(2))

                    dip_res_z = dip_res_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(3))

                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
         dip_res_tot = sqrt (dip_res_x**2 + dip_res_y**2 + dip_res_z**2 )
         !dipRes_x = dipRes_x + dip_res_x
         !dipRes_y = dipRes_y + dip_res_y
         !dipRes_z = dipRes_z + dip_res_z
         res_dip(1,iatom) = dip_res_x
         res_dip(2,iatom) = dip_res_y
         res_dip(3,iatom) = dip_res_z
      end do ! end do iatom = 1,natoms

          
         do iatom = 1, natoms
         dip_x = 0.0d0
         dip_y = 0.0d0
         dip_z = 0.0d0
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         Qtot=0.0d0
         Qtot1=0.0d0
         do issh = 1,nssh(in1)
             Qtot1 = Qtot1+Qin(issh,iatom)
         end do !end do issh = 1,nssh(in1)
         Qtot=0.0d0
         Qtot2=0.0d0
         do issh = 1,nssh(in1)
             Qtot2 = Qtot2+Qout(issh,iatom)
         end do !end do issh = 1,nssh(in1)
         jatom=iatom
         ineigh=neigh_self(iatom)
         in2=in1
         r2(:)=r1(:)
         Rbc(:)=(r1(:)+r2(:))/2.0d0
             do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)

                if ( (imu .ne. inu)) then
                !else 
                  dip_x = dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y = dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z = dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))
!Dintra escribir la matrix para de N*xyz
                end if !end if
               end do !end do inu
            end do !end do imu
            !write(*,*) 'dipc(1,imu,inu,ineigh,iatom)',dipc(1,imu,inu,ineigh,iatom)
            !write(*,*) 'dipc(2,imu,inu,ineigh,iatom)',dipc(2,imu,inu,ineigh,iatom)
            !write(*,*) 'dipc(3,imu,inu,ineigh,iatom)',dipc(3,imu,inu,ineigh,iatom)
            intra_dip(1,iatom) = dip_x + res_dip(1,iatom)!/Debye
            intra_dip(2,iatom) = dip_y + res_dip(2,iatom)!/Debye
            intra_dip(3,iatom) = dip_z + res_dip(3,iatom)!/Debye
            dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
           write(667,445) symbol(iatom), Qtot1, Qtot2, dip_x/Debye, dip_y/Debye, dip_z/Debye, dip_tot/Debye
      end do !end do iatom = 1,natoms

      dq_DP = 0.0d0

      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         rmedio = 0.0d0
         w_suma = 0.0d0

         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            w_k(ineigh)=exp(-klambda*((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            rmedio = rmedio + w_k(ineigh)*r2
            w_suma = w_suma + w_k(ineigh)

         end do !end ineigh = 1,natoms

         rmedio = rmedio / w_suma

         bwrr = 0.0d0

         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)
            do imu = 1,3 !xyz
              do inu = 1,3 !xyz
                bwrr(inu,imu) = bwrr(inu,imu) + w_k(ineigh)*r2(imu)*(r2(inu)-rmedio(inu))
              enddo !inu
            enddo !imu
         end do !end ineigh = 1,natoms

         !inversa de bwrr
      n_bwrr = 3
      lda_bwrr = n_bwrr
      lwork = n_bwrr
      !Allocate (a(lda,n), work(lwork), ipiv(n))
      bwrr_inv = bwrr

      !write(*,*) 'bwrr1', bwrr
      !call dgetrf(n_bwrr, n_bwrr, bwrr_inv, lda_bwrr, ipiv, info)
      !write(*,*) 'bwrr2', bwrr_inv

      !write(*,*) 'info', info
      !if (info==0) Then
      !Compute inverse of A
      !call dgetri(n_bwrr, bwrr_inv, lda_bwrr, ipiv, lwork, n_bwrr, info)
      !else
        !write (*, *) 'The factor U is singular, go to pseudoinverse'
        bwrr_inv = bwrr
        call dgesvd('A', 'S', n_bwrr, n_bwrr, bwrr_inv, lda_bwrr,s_bwrr, u_bwrr, lda_bwrr, vt_bwrr, lda_bwrr, dummy, lwork, info)
        write(*,*) 's_bwrr', s_bwrr
        !write(*,*) 'u_bwrr', u_bwrr
        !write(*,*) 'vt_bwrr', vt_bwrr
        zero_bwrr = 0.00
        do i = 1,3
          if (abs(s_bwrr(i)) .gt. 0.000001) then
            zero_bwrr(i,i) = 1/s_bwrr(i)
          endif
        enddo
        v_bwrr = transpose(vt_bwrr)
        ut_bwrr = transpose(u_bwrr)
        bwrr_inv=matmul(v_bwrr,zero_bwrr)
        bwrr_inv=matmul(bwrr_inv,ut_bwrr)
      !end if
      !write(*,*) 'bwrr3', bwrr_inv

      intra_dip_aux(1,1)=intra_dip(1,iatom)
      intra_dip_aux(2,1)=intra_dip(2,iatom)
      intra_dip_aux(3,1)=intra_dip(3,iatom)
      delta_ck = matmul(bwrr_inv,intra_dip_aux)
      !write(*,*) 'bwrr_inv', bwrr_inv
      !write(*,*) 'intra_dip_aux', intra_dip_aux
      !write(*,*) 'delta_ck', delta_ck

      do ineigh = 1,neighn(iatom)
        mbeta = neigh_b(ineigh,iatom)
        jatom = neigh_j(ineigh,iatom)
        r2(:) = ratom(:,jatom) + xl(:,mbeta)
        raux = r2(:) - rmedio(:)
        c_k(jatom) = w_k(ineigh) * (raux(1)*delta_ck(1,1)+raux(2)*delta_ck(2,1)+raux(3)*delta_ck(3,1))
        write(*,*) 'c_k',jatom,  c_k(jatom)
        dq_DP(jatom) = dq_DP(jatom) + c_k(jatom)
      enddo

      !write(*,*) 'Ha pasado por aqui'

      end do ! end do iatom = 1,natoms

      do i =1,natoms
        write(*,*) 'dq_DP (',i,') =',dq_DP(i)
      enddo
      
      do iatom = 1, natoms
        Q0_TOT(iatom) = 0
        in1 = imass(iatom)
        do issh = 1, nssh(in1)
          Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
        end do

      end do

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0

      open( unit = 800, file = 'CHARGES'//trim(estring), status ='unknown', position = 'append')
      do iatom = 1, natoms
         Qtot=-Q0_TOT(iatom)+dq_DP(iatom)
         in1 = imass(iatom)
         do imu = 1,num_orb(in1)
             Qtot = Qtot+Qout(imu,iatom)
         end do
         write(*,*) 'Q_new (',iatom,') =',-Qtot
             
         write(800,'(2x, i4,2x,f10.6)') iatom, Q0_TOT(iatom)+Qtot
         !write(*,*) 'ratom is ',ratom(1,iatom),ratom(2,iatom),ratom(3,iatom)
         dip_x = dip_x+Qtot*ratom(1,iatom)
         dip_y = dip_y+Qtot*ratom(2,iatom)
         dip_z = dip_z+Qtot*ratom(3,iatom)

      enddo !end do iatom = 1,natoms
      close(800)

        dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )   

         dipQout_x=dip_x
         dipQout_y=dip_y
         dipQout_z=dip_z
         dipQout_tot=dip_tot

               open( unit = 667, file = 'Charges_and_Dipoles', status = 'unknown',  &
                                   &     position = 'append')
             write(667,*)   '+++++++++++++++++++ N E W   S T E P +++++&
                               &+++++++++++++' 
             write(667,444) dipTot_x/Debye, dipTot_y/Debye, dipTot_z/Debye, dipTot_tot/Debye
             write(667,444) dipQin_x/Debye, dipQin_y/Debye, dipQin_z/Debye, dipQin_tot/Debye
             write(667,444) dipQout_x/Debye, dipQout_y/Debye, dipQout_z/Debye, dipQout_tot/Debye
             write(667,444) dipIntra_x/Debye, dipIntra_y/Debye, dipIntra_z/Debye, dipIntra_tot/Debye
    
     end if !end iwrtdipole .eq. 2
           
           
      
     end if !end if (iwrtdipole .gt. 0)

! END OF WRITING THE DIPOLE
!******************************************************************************
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
400     format (9f9.4)
500     format (70('='))
501     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Charges ')
502     format (3x, i5, 7x, a2, 5x, i2, 4x, 8(1x, f5.2))
503     format (3x, i5, 7x, a2, 4x, f10.4)
600     format (2x, i5, 2x, a40, 2x, i2)
601     format (2x, 10f14.8)
444     format (4f10.4)
445     format (a2,6f10.4)

        return
        end
