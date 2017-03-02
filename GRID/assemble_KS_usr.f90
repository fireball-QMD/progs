! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! assemble_KS_usr.f90
! Program Description
! ===========================================================================
!       This routine computes the two-center contribution to the total
! energy. This routine program computes the terms u0(iatom,ineigh) and
! uee00(iatom) for the short ranged potential, usr = u0 - uee. Here iatom
! is a basis atom in the central cell and ineigh is the ineigh'th neighbor
! to iatom, also the long rangeed contribution is calculated (the
! information comes from ewald.f). The results are converted to eV energy
! units.
 
! This routine computes derivatives only if iforce = 1. This routine also
! computes the force derivative with respect to ratom of the short-ranged
! energy per cell, thus dusr(3,iatom) = - d/d(ratom(3,iatom)) usr.
! Here ratom is the basis atom position in the central cell. The minus sign
! makes it force-like.
!
! The u0 interaction is:
! -1/2 * int d3r  (n(nuclear)*vion(r) + n0 * vh0),
! where n(nuclear) is the nuclear charge density, vion the local ion
! potential, n0 the neutral atom charge density, and vh0
! the hartree potential due to neutral atoms.
!
!
! ===========================================================================
! Original code from Otto F. Sankey with modifications by Alex A. Demkov
! and Jose Ortega (for charge transfer interactions).
 
! Code rewritten by:
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
        subroutine assemble_KS_usr (iforce, uiiuee)

        use charges
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
! Output
        real, intent (out) :: uiiuee
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ideriv
        integer in1, in2
        integer index
        integer index_coulomb
        integer ineigh
        integer interaction
        integer issh
        integer jatom
        integer jssh
        integer mbeta
        integer n1, n2
        integer imu, inu
        integer l1,l2
   
        real distance
        real dq1, dq2
        real dqi, dqj
        real eklr
        real qi, qj
        real QQ
        real u0tot
        real ue0tot
        real xforce
        real Zi, Zj
        real Qni, Qnj
 
        real, dimension (natoms, neigh_max) :: corksr
        real, dimension (nsh_max, nsh_max) :: coulomb
        real, dimension (nsh_max, nsh_max) :: coulombD
        real, dimension (3) :: dcorksr
        real, dimension (ME2c_max) :: dslist
        real, dimension (3) :: eta
        real, dimension (natoms) :: Q, Q0
        real, dimension (3) :: r1, r2
        real, dimension (ME2c_max) :: slist
        real, dimension (natoms, neigh_max) :: u0
        real, dimension (natoms) :: uee00
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' Welcome to assemble_usrG.f! '
        write (*,*) '  '
 
! Initialize arrays
        dxcv = 0.0d0
        dusr = 0.0d0
        u0 = 0.0d0
 
! Calculate delta charges (integer) into a real variable.
        do iatom = 1, natoms
         Q(iatom) = 0.0d0
         Q0(iatom) = 0.0d0
         in1 = imass(iatom)
         do issh = 1,nssh(in1)
          Q0(iatom) = Q0(iatom) + Qneutral(issh,in1)
         end do
         do imu = 1, num_orb(in1)
          Q(iatom) = Q(iatom) + Qin(imu,iatom)
         end do
        end do
   
! Loop over all atoms in the central cell.
        do iatom = 1, natoms
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
      
! Initialize the charge on iatom:
         qi = Q(iatom)
         Zi = Q0(iatom)
 
! Determine dqi and dq1:
         dqi = Q(iatom) - Q0(iatom)
         dq1 = dq(in1)
          
! Loop over all neighbors ineigh of iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = xl(:,mbeta) + ratom(:,jatom)
          in2 = imass(jatom)
 
! Initialize the charge on jatom:
          qj = Q(jatom)
          Zj = Q0(jatom)
          QQ = Zi*Zj - qi*qj
! Calculate the distance between the two atoms.
          distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2            &
     &                                       + (r2(3) - r1(3))**2)
 
! ****************************************************************************
!
! GET COULOMB INTERACTIONS 
! ****************************************************************************
! Now find the three coulomb integrals need to evaluate the neutral
! atom/neutral atom hartree interaction.
! Loop over all the non-zero integrals for this interaction:
          index_coulomb = nssh(in1)*nssh(in2)
          interaction = 12
          ideriv = 0
          do index = 1, index_coulomb
           call interpolate_1d (interaction, ideriv, in1, in2, index,        &
     &                          iforce, distance, slist(index), dslist(index))
          end do
 
! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
          n1 = nssh(in1)
          n2 = nssh(in2)
          call recoverC (n1, n2, slist, dslist, coulomb, coulombD)
         
! Actually, now we calculate not only the neutral atom contribution,
! but also the short-ranged contribution due to the transfer of charge
! between the atoms:
!
! (Eii - Eee)neut - SUM(short range)(n(i) + dn(i))*dn(j)*J(i,j),
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! SPECIAL CASE: SELF-INTERACTION
           uee00(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in1)
             uee00(iatom) = uee00(iatom)                                    &
     &         + Qneutral(issh,in1)*Qneutral(jssh,in2)*coulomb(issh,jssh)
            end do ! do jssh
           end do ! do issh

! put the half and the units in:
           uee00(iatom) = uee00(iatom)*(eq2/2.0d0)
           u0(iatom,ineigh) = 0.0d0

          else
! BONAFIDE TWO ATOM CASE
! Compute u0
           u0(iatom,ineigh) = 0.0d0
           imu = 1 
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             u0(iatom,ineigh) = u0(iatom,ineigh)                            &
     &         + Qneutral(issh,in1)*Qneutral(jssh,in2)*coulomb(issh,jssh)
            end do ! do jssh
           end do ! do issh

           u0(iatom,ineigh) = (eq2/2.0d0)*(Zi*Zj/distance - u0(iatom,ineigh))
! ***************************************************************************
!
!                                FORCES
! ***************************************************************************
           if (iforce .eq. 1) then
            eta(:) = (r2(:) - r1(:))/distance
 
! Put in the forces due to the charge transfer. This 'sumit' has the sign of
! d/rd1, and is NOT force-like. We put in force-like character later in dusr.
            xforce = 0.0d0

            imu = 1 
            do issh = 1, nssh(in1)
             do jssh = 1, nssh(in2)
               xforce = xforce                                               &
     &          + Qneutral(issh,in1)*Qneutral(jssh,in2)*coulombD(issh,jssh)
              end do ! do jssh
            end do ! do issh

            dusr(:,iatom) = dusr(:,iatom)                                    &
     &       - eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
            dusr(:,jatom) = dusr(:,jatom)                                    &
     &       + eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)

           end if                  ! end if (forces)
          end if                   ! end if (iatom .eq. jatom)
 
! End of loop over neighbors
         end do
 
! End of loop over iatom
        end do

 
! ***************************************************************************
!
! Compute the total cell value of uii-uee; uii-uee = sum u0(i,m) - sum uee00(i)
! Add long-range ewald interactions.
! ***************************************************************************
        u0tot = 0.0d0
        ue0tot = 0.0d0
        do iatom = 1, natoms
         ue0tot = ue0tot + uee00(iatom)
         do ineigh = 1, neighn(iatom)
          u0tot = u0tot + u0(iatom,ineigh)
         end do
        end do

        uiiuee = u0tot - ue0tot
        
! Format Statements
! ===========================================================================
 
        return
      end subroutine assemble_KS_usr
