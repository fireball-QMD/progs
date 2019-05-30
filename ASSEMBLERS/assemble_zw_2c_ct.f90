! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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

 
! assemble_zw_2c_ct.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!This subroutine computes charge transfer contribution to the XC matrix
!elements using second order terms in the expansion for the XC energy
!around the neutral density
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_zw_2c_ct (nprocs, iforce, iordern)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use scf, only : Kscf
        use energy, only : uxcdcc_zw
        use integrals, only : xcnu1c
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer icount
        integer icount_sav
        integer ierror
        integer imu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer jcount
        integer jcount_sav
        integer jssh
        integer kforce
        integer matom
        integer matom2
        integer mbeta
        integer my_proc
        integer natomsp
        integer ix
        integer iy
        integer igamma
        integer :: count_l, count_l_ini, issh1, issh2
        integer l
 
        real dq1
        real dq2
        real dterm
        real dterm_1
        real dterm_2
        real dstn_temp
        real dxn
        real dqdc
        real rcutoff_j
        real rend
        real rend1
        real rend2
        real sterm_1
        real sterm_2
        real y
        real rcutoff_i
        real :: qmu, q0mu, dqmu
        real :: A,B
 
        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (3, nsh_max, nsh_max) :: bccapx
        real, dimension (nsh_max, nsh_max) :: bccax
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
! JOM-JIMM
!       real, dimension (numorb_max, numorb_max) :: stn1
!       real, dimension (numorb_max, numorb_max) :: stn2
        real stn1
        real stn2

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE


! Procedure
! ===========================================================================
! Initialize interactions to zero.
      dxcdcc_zw = 0.0d0
      bcca = 0.0d0
      bccax = 0.0d0
      bccapx = 0.0d0
      vxc_ca = 0.0d0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
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

! Loop over the atoms in the central cell.
!!$omp parallel do private (icount, icount_sav, in1, in2, in3, interaction)   &
!!$omp&            private (isorp, jatom, jcount, jcount_sav, kforce, matom)  &
!!$omp&            private (mbeta, dq1, dq2, dterm_1, dterm_2, dstn_temp)     &
!!$omp&            private (dxn, rcutoff_j, rend, rend1, rend2, sterm_1)      &
!!$omp&            private (sterm_2, stn_temp1, stn_temp2, y, bcca, bccapx)   &
!!$omp&            private (bccax, deps, eps, dipx, dippx, emnpl, r1, r2, r21)&
!!$omp&            private (sighat, stn1, stn2)



!+++++++++++ COMPUTE INTEGRALS G2NU  !IF KSCF = 1
       if (Kscf .eq. 1) then   !Only do this in the first step of the SCF loop! 
          !do iatom = iatomstart, iatomstart - 1 + natomsp
         g2nu = 0.0d0
         g2nup = 0.0d0
          do iatom = 1, natoms
           r1(:) = ratom(:,iatom)
           in1 = imass(iatom)
           matom = neigh_self(iatom)
           do ineigh = 1, neighn(iatom)        
             mbeta = neigh_b(ineigh,iatom)
             jatom = neigh_j(ineigh,iatom)
             r2(:) = ratom(:,jatom) + xl(:,mbeta)
             in2 = imass(jatom) 
             r21(:) = r2(:) - r1(:)
             y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
             if (y .lt. 1.0d-05) then
               sighat(1) = 0.0d0
               sighat(2) = 0.0d0
               sighat(3) = 1.0d0
             else
               sighat(:) = r21(:)/y
             end if
             call epsilon (r2, sighat, eps)
             call deps2cent (r1, r2, eps, deps)
             bcca = 0.0d0
             kforce = 1
             if (matom .ne. ineigh) then
             interaction =  14   !????
             in3 = in2
               isorp = 0
               !use doscentrosS !!! (for the time being...)
              call doscentrosS (interaction, isorp, kforce, in1, in2, in3,y, eps, bccax, bccapx)
               do issh1 = 1, nssh(in1)
                  do issh2 = 1,nssh(in2)
                    g2nu(issh1,issh2,ineigh,iatom)=bccax(issh1,issh2)  !store the integrals
                    g2nup(:,issh1,issh2,ineigh,iatom)=bccapx(:,issh1,issh2) !store the integrals
                  end do !end do issh2 
               end do !end do issh1
             else !if matom .eq. ineigh, case onecenter
                do issh1 = 1, nssh(in1)
                  do issh2 = 1, nssh(in2)
                   g2nu(issh1,issh2,matom,iatom) = xcnu1c(issh1,issh2,in1)
                  end do !end do issh2
                end do !end do issh1
             end if ! end if matom .ne. ineigh
           end do ! end do ineigh = 1, neighn(iatom)
          end do !end do iatom = iatomstart, iatomstart - 1 + natomsp
         end if ! end if Kscf .eq. 1
!++++++++++ END OF COMPUTING G2NU INTEGRALS


        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
! Find charge on iatom
         dq1 = 0.0d0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)         ! <==== loop 2 over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          matom2 = neigh_self(jatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
 
! ****************************************************************************
!
!                     VNA FOR ATOM CASE
! ****************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
! We need jatom because the interaction will include both neutral atom and
! isorp pieces. The isorp pieces will invole Qin. Here is a snippet from
! doscenatm:
!         scam(i,j)=scam(i,j)+temp(i,j)*(Qin(isorp,jk)-Qneutral(isorp,jk)
 
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
 
! Initialize bcca to zero for the charge atom interactions.
          bcca = 0.0d0
 
          do isorp = 1, nssh(in2)

           dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
 
! Now correct bccax by doing the stinky correction.
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
           in3=in1
           !do inu = 1,num_orb(in3)
             do imu = 1,num_orb(in1)
                issh1=orb2shell(imu,in1)
                !issh2=orb2shell(inu,in3)
                bcca(imu,imu) = bcca(imu,imu) + g2nu(issh1,isorp,ineigh,iatom)*dxn 
                !Here in the atom case we only have diagonal terms
             end do !end do imu = 1, mum_orb(in1)     
           !end do !end do inu = 1,num_orb(in3)

          !Double counting correction:
          do issh1 = 1,nssh(in1)
             dqdc = (Qin(issh1,iatom) - Qneutral(issh1,in1))
             !dqdc=0.0d0
             uxcdcc_zw = uxcdcc_zw - (Qin(issh1,iatom)-0.50*dqdc)*g2nu(issh1,isorp,ineigh,iatom)*dxn
             !change sign to make dxcdcc_zw force-like
             dxcdcc_zw(:,ineigh,iatom) = dxcdcc_zw(:,ineigh,iatom) &
             & + (Qin(issh1,iatom)-0.50*dqdc)*g2nup(:,issh1,isorp,ineigh,iatom)*dxn
          end do ! end do issh1 = 1,

!------------------------------------------------------------------------------------
!           count_l=1
!           do issh = 1, nssh(in1)
             !bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn 
             !For the time being, bccax(imu,inu) in the atom case is
             !zero whenever imu .neq. inu
!             l=lssh(issh,in1)
!             count_l_ini=count_l
!             qmu=Qin(issh,iatom)/(2.0*l+1)
!             q0mu=Qneutral(issh,in1)/(2.0*l+1)
!             
!             dqmu=qmu-q0mu
!             do imu = count_l_ini,count_l_ini+2*l
!                   bcca(imu,imu) = bcca(imu,imu) + g2nu(issh,isorp,ineigh,iatom)*dxn  
                   !xc charge-transfer part double counting (energy)
          !        uxcdcc_zw = uxcdcc_zw- &         
          !        & (qmu-0.5*dqmu)*g2nu(issh,isorp,ineigh,iatom)*dxn
                   !xc charge-transfer part double counting (forces)
          !        dxcdcc_zw(:,ineigh,iatom) = dxcdcc_zw(:,ineigh,iatom)- &
          !        & (qmu-0.5*dqmu)*g2nup(:,issh,isorp,ineigh,iatom)*dxn
!               count_l=imu+1
!             end do ! end do imu = count_l_ini, count_l_ini+l
!           end do !end do issh
!-----------------------------------------------------------------------------------------


          end do !end do isorp = 1, nssh(in2) (line 270)
          in3 = in1
          !bcca = 0.0d0
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
!$omp atomic
            vxc_ca(imu,inu,matom,iatom) = vxc_ca(imu,inu,matom,iatom) +             &
        &          bcca(imu,inu)
           end do !end do inu
          end do !end do imu

          !TEST DC
        !  do imu = 1,num_orb(in1)
        !     issh1 = orb2shell(imu,in1)
        !     l=lssh(issh1,in1)
        !     qmu=Qin(issh1,iatom)/(2.0*l+1)
        !     q0mu=Qneutral(issh1,in1)/(2.0*l+1)             
        !     dqmu=qmu-q0mu
!
        !            uxcdcc_zw = uxcdcc_zw -  &
      !           &       (qmu-0.50*dqmu)*vxc_ca(imu,imu,matom,iatom) 
        !         end do !end do imu
! ****************************************************************************
!
!                       GET VNA FOR ONTOP CASE 
! ****************************************************************************
! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
 
          else

             r21(:) = r2(:) - r1(:)
             y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

 
! Initialize bcca for charged atom interactions.
           bcca = 0.0d0


! .... The following integrals, ontopl and ontopr, can be computed, as
! in the 3-center case, by splitting in two with the coefficients from
! the Mulliken-Dipole projection.  
! For the vna_ontopl case, the potential is on the first atom (j):
! Charged atom piece
           in3 = in2
           do isorp = 1, nssh(in1)
!            call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
!     &                       eps, deps, bccax, bccapx)
            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
                !possibility: 
                   issh1=orb2shell(imu,in1)
                   issh2=orb2shell(inu,in3)
                   !Here is the split of the INTEGRAL gcab...---> 
                     A=0.5*s_mat(imu,inu,ineigh,iatom)-dip(imu,inu,ineigh,iatom)/y
                     B=0.5*s_mat(imu,inu,ineigh,iatom)+dip(imu,inu,ineigh,iatom)/y  
                      bcca(imu,inu) = bcca(imu,inu)+ &
                      & A*g2nu(issh1,isorp,matom,iatom)*dxn+B*g2nu(isorp,issh2,ineigh,iatom)*dxn
                 !end of possibility
             end do !end do inmu
            end do !end do inu
           end do !end do isorp
 
! For the vna_ontopr case, the potential is on the second atom (j):
! Charged atom piece
           in3 = in2
           do isorp = 1, nssh(in2)
!            call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
!     &                       eps, deps, bccax, bccapx)
 
            dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
             !possibility:
                  issh1=orb2shell(imu,in1)
                  issh2=orb2shell(inu,in3)
                 !Here is the split of the INTEGRAL gcab...--->    
                   A=0.5*s_mat(imu,inu,ineigh,iatom)-dip(imu,inu,ineigh,iatom)/y
                   B=0.5*s_mat(imu,inu,ineigh,iatom)+dip(imu,inu,ineigh,iatom)/y  
                   bcca(imu,inu) = bcca(imu,inu)+ &
                   & A*g2nu(issh1,isorp,ineigh,iatom)*dxn+B*g2nu(issh2,isorp,matom2,jatom)*dxn     
              !end of possibility
             end do !end do imu
            end do !end do inu
           end do !end do isorp
 
! Now put into vca.
           !bcca = 0.0d0
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
              vxc_ca(imu,inu,ineigh,iatom) =                                     &
      &        vxc_ca(imu,inu,ineigh,iatom) + bcca(imu,inu)
            end do
           end do
 
! End if for r1 .ne. r2 case
          end if
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do ! do ineigh
        end do ! do iatom

        !    vxc_ca = 0.0d0

! Format Statements
! ===========================================================================



        return
        end subroutine assemble_zw_2c_ct
