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
 
! assemble_usr.f90
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
! -1/2 * int(slash) d3r  (n(nuclear)*vion(r) + n0 * vh0),
! where n(nuclear) is the nuclear charge density, vion the local ion
! potential, n0 the neutral atom charge density, and vh0
! the hartree potential due to neutral atoms.
!
! This routine also computes the xc double counting correction. First, we
! get the data for the atoms in1, in2 at r1,r2. This call will return five
! values: one for the neutral pair, and four with some predetermined charge
! transfer (dq of the charged shell for xc, set in CREATOR: e.g., Si.inc).
! We then interpolate for the current charge distribution of the pair.
! We calculate:
!   (n1+n2)*(exc(1+2)-muxc(1+2)) - n1*(exc(1)-xcmu(1))
!                                - n2*(exc(2)-xcmu(2))
!
! ===========================================================================
! Original code from Otto F. Sankey with modifications by Alex A. Demkov
! and Jose Ortega (for charge transfer interactions).
 
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_usr (itheory, itheory_xc, iforce,   &
     &                           uxcdcc, uiiuee)
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
        integer, intent (in) :: itheory
        integer, intent (in) :: itheory_xc
 
! Output
        real, intent (out) :: uiiuee
        real, intent (out) :: uxcdcc
 
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
        integer non2c
   
        real distance
        real dq1, dq2
        real dqi, dqj
        real dxc
        real dxc00, dxc0P, dxc0M, dxcP0, dxcM0
        real eklr
        real qi, qj
        real QQ
        real u0tot
        real ue0tot
        real xc
        real xc00, xc0P, xc0M, xcP0, xcM0
        real xforce
        real Zi, Zj

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
        write (*,*) ' Welcome to assemble_usr.f! '
        write (*,*) '  '
 
! Initialize arrays
        dxcv = 0.0d0
        dusr = 0.0d0
        u0 = 0.0d0
        uxcdcc = 0.0d0
 

! Calculate delta charges (integer) into a real variable.
        do iatom = 1, natoms
         Q(iatom) = 0.0d0
         Q0(iatom) = 0.0d0
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          Q(iatom) = Q(iatom) + Qin(issh,iatom)
          Q0(iatom) = Q0(iatom) + Qneutral(issh,in1)
         end do
        end do
   
! Loop over all atoms in the central cell.
!!$omp parallel do private (r1, in1, qi, Zi, dqi, dq1, mbeta, jatom, r2, in2) &
!!$omp&   private (qj, Zj, QQ, distance, index_coulomb, interaction, ideriv)  &
!!$omp&   private (slist, n1, n2, coulomb, coulombD)                          &
!!$omp&   private (non2c, dqj, dq2, xc, dxc, xc00, dxc00, xcM0, dxcM0)        &
!!$omp&   private (xcP0, dxcP0, xc0M, dxc0M, xc0P, dxc0P, eta, xforce)
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
          distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2       &
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
           call interpolate_1d (interaction, ideriv, in1, in2, index,   &
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
           if (itheory .eq. 1) then
            do issh = 1, nssh(in1)
             do jssh = 1, nssh(in1)
              uee00(iatom) = uee00(iatom)                               &
     &         + Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
             end do
            end do
           else if (itheory .eq. 0 .or. itheory .eq. 2) then
            do issh = 1, nssh(in1)
             do jssh = 1, nssh(in1)
              uee00(iatom) = uee00(iatom) +                             &
     &         Qneutral(issh,in1)*Qneutral(jssh,in2)*coulomb(issh,jssh)
             end do
            end do
           end if
 
! put the half and the units in:
           uee00(iatom) = uee00(iatom)*(eq2/2.0d0)
           u0(iatom,ineigh) = 0.0d0
           corksr(iatom,ineigh) = 0.0d0
          else
 
! BONAFIDE TWO ATOM CASE
! Compute u0
           u0(iatom,ineigh) = 0.0d0
           if (itheory .eq. 1) then
            do issh = 1, nssh(in1)
             do jssh = 1, nssh(in2)
              u0(iatom,ineigh) = u0(iatom,ineigh) +                     &
     &          Qin(issh,iatom)*Qin(jssh,jatom)*coulomb(issh,jssh)
             end do
            end do
           else if (itheory .eq. 0 .or. itheory .eq. 2) then
            do issh = 1, nssh(in1)
             do jssh = 1, nssh(in2)
              u0(iatom,ineigh) = u0(iatom,ineigh) +                     &
     &          Qneutral(issh,in1)*Qneutral(jssh,in2)*coulomb(issh,jssh)
             end do
            end do
           end if
           u0(iatom,ineigh)=(eq2/2.0d0)*(Zi*Zj/distance-u0(iatom,ineigh))
 
! This is a correction for the extention of the long-ranged sum to all the
! atoms in the system, which we do in order to use the Ewald summation
! technique. This energy is to be subtracted from the total, since it is
! correctly calculated by integrals, but overcounted by the ewald
! contribution. We include the minus HERE, which means we subtract it from
! etot.
           corksr(iatom,ineigh) = - (eq2/2.0d0)*QQ/distance
 
 
! ****************************************************************************
!
! XC DOUBLE COUNTING CORRECTION
! ****************************************************************************
           if (itheory_xc .eq. 0) then 

! First we get the data for the atoms in1, in2 at r1,r2. This call will return
! five values: one for the neutral pair, and four with some predetermined
! charge transfer (dq of the charged shell for xc, set in CREATOR: e.g. Si.inc).
! Then interpolate for the currnt charge distribution of the pair. Here is the
! key from creator:
!
! We calculate   (n1+n2)*(exc(1+2)-muxc(1+2)) - n1*(exc(1)-xcmu(1))
!                                             - n2*(exc(2)-xcmu(2))
! The catch comes in when we compute derivatives. We compute
! neutral,neutral for ideriv1. For other ideriv's we have the following KEY:
!
! (xy) means charge on (1,2). Case 1 (KEY=1), neutral neutral corresponds to
! (00) etc. KEY = 1,2,3,4,5 for ideriv=1,2,3,4,5
!
!                                   |
!                                   + (0+) KEY=5
!                                   |
!                                   |
!                                   |
!                       KEY=2       |  KEY=1
!                      (-0)         |(00)        (+0) KEY=3
!                     --+-----------O-----------+--
!                                   |
!                                   |
!                                   |
!                                   |
!                                   |
!                                   +(0-) KEY=4
!                                   |
            interaction = 8
            non2c = 1
            ideriv = 0
            call interpolate_1d (interaction, ideriv, in1, in2, non2c,  &
     &                           iforce, distance, xc00, dxc00)
 
            if (itheory .eq. 1) then
             ideriv = 1
             call interpolate_1d (interaction, ideriv, in1, in2, non2c, &
     &                            iforce, distance, xcM0, dxcM0)
             ideriv = 2
             call interpolate_1d (interaction, ideriv, in1, in2, non2c, &
     &                            iforce, distance, xcP0, dxcP0)
             ideriv = 3
             call interpolate_1d (interaction, ideriv, in1, in2, non2c, &
     &                            iforce, distance, xc0M, dxc0M)
             ideriv = 4
             call interpolate_1d (interaction, ideriv, in1, in2, non2c, &
     &                            iforce, distance, xc0P, dxc0P)
            end if
 
! Determine dqj and dq2:
            dqj = Q(jatom) - Q0(jatom)
            dq2 = dq(in2)
 
! The above is only to get this thing going.
! Now we should interpolate - the fast way is:
!
!   e(dqi,dqj) = exc(0,0) + dqi*(exc(1,0) - exc(0,0))/dQ
!                         + dqj*(exc(0,1) - exc(0,0))/dQ
!
! Here dQ is coming from creator.
! The good way is to use a three point Lagrange interpolation along the axis.
!
! Lagrange: f(x) = f(1)*L1(x) + f(2)*L2(x) + f(3)*L3(x)
!
! L1(x) = (x - x2)/(x1 - x2)*(x - x3)/(x1 - x3)
! L2(x) = (x - x3)/(x2 - x1)*(x - x3)/(x2 - x3)
! L3(x) = (x - x1)/(x3 - x1)*(x - x2)/(x3 - x2)
!
! in our case:
!
! L1(dq) = (1/2)*dq*(dq - 1)
! L2(dq) = -(dq + 1)*(dq - 1)
! L3(dq) = (1/2)*dq*(dq + 1)
!
! The interpolation does not depend on the (qi,qj) quadrant:
!
!  f(dqi,dqj)=f(0,dqj)+f(dqi,0)-f(0,0)
 
! Neutral case:
            xc = xc00
            dxc = dxc00
 
! Non-neutral case:
! Do this for DOGS only!
            if (itheory .eq. 1) then

! There are four cases, note the (ge,gt,le,lt) set:
! (+,+) case I
             if (dqi .gt. 0.0d0 .and. dqj .gt. 0.0d0) then
              xc = xc + ((xcP0 - xc00)/dq1)*dqi +                       &
     &                  ((xc0P - xc00)/dq2)*dqj
              dxc = dxc + ((dxcP0 - dxc00)/dq1)*dqi +                   &
     &                    ((dxc0P - dxc00)/dq2)*dqj
             end if
 
! (-,-) case II
             if (dqi .lt. 0.0d0 .and. dqj .lt. 0.0d0) then
              xc = xc + ((xc00 - xcM0)/dq1)*dqi +                       &
     &                  ((xc00 - xc0M)/dq2)*dqj
              dxc = dxc + ((dxc00 - dxcM0)/dq1)*dqi +                   &
     &                    ((dxc00 - dxc0M)/dq2)*dqj
             end if
 
! (+,-) case III
             if (dqi .gt. 0.0d0 .and. dqj .lt. 0.0d0) then
              xc = xc + ((xcP0 - xc00)/dq1)*dqi +                       &
     &                  ((xc00 - xc0M)/dq2)*dqj
              dxc = dxc + ((dxcP0 - dxc00)/dq1)*dqi +                   &
     &                    ((dxc00 - dxc0M)/dq2)*dqj
             end if
 
! (-,+) case IV
             if (dqi .lt. 0.0d0 .and. dqj .ge. 0.0d0) then
              xc = xc + ((xc00 - xcM0)/dq1)*dqi +                       &
     &                  ((xc0P - xc00)/dq2)*dqj
              dxc = dxc + ((dxc00 - dxcM0)/dq1)*dqi +                   &
     &                    ((dxc0P - dxc00)/dq2)*dqj
             end if
            end if
 
! Now we add the contribution to the total. Notice the one half factor, it is
! the sum over (iatom,jatom) with iatom not equal to jatom.
!!$omp atomic
            uxcdcc = uxcdcc + xc/2.0d0
           end if                          !end if(itheory_xc)
 
! ***************************************************************************
!
!                                FORCES
! ***************************************************************************
           if (iforce .eq. 1) then
            eta(:) = (r2(:) - r1(:))/distance
 
! Put in the forces due to the charge transfer. This 'sumit' has the sign of
! d/rd1, and is NOT force-like. We put in force-like character later in dusr.
            xforce = 0.0d0
            if (itheory .eq. 1) then
             do issh = 1, nssh(in1)
              do jssh = 1, nssh(in2)
               xforce = xforce +                                        &
     &         Qin(issh,iatom)*Qin(jssh,jatom)*coulombD(issh,jssh)
              end do
             end do
            else if (itheory .eq. 0 .or. itheory .eq. 2) then
             do issh = 1, nssh(in1)
              do jssh = 1, nssh(in2)
               xforce = xforce +                                        &
     &         Qneutral(issh,in1)*Qneutral(jssh,in2)*coulombD(issh,jssh)
              end do
             end do
            end if
              
            dusr(:,iatom) = dusr(:,iatom) -                             &
     &        eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)
            dusr(:,jatom) = dusr(:,jatom) +                             &
     &        eta(:)*(eq2/2.0d0)*(Zi*Zj/distance**2 + xforce)

! Now we add the corksr correction. Both of these are d/dr1
! derivatives and are NOT force-like.
            if (itheory .eq. 1) then
             dcorksr(:) = - eta(:)*(eq2/2.0d0)*QQ/distance**2
             dusr(:,iatom) = dusr(:,iatom) - dcorksr(:)
             dusr(:,jatom) = dusr(:,jatom) + dcorksr(:)
            end if
 
! XC-DOUBLE COUNTING FORCE:
            if (itheory_xc .eq. 0) then
             dxcv(:,iatom) = dxcv(:,iatom) + eta(:)*dxc/2.0d0
             dxcv(:,jatom) = dxcv(:,jatom) - eta(:)*dxc/2.0d0
            end if
           end if                  ! end if (forces)
          end if                   ! end if (iatom .eq. jatom)
 
! End of loop over neighbors
         end do
 
! End of loop over iatom
        end do

! Subtract the forces for the ewald interaction
! The variable fewald is already force-like.
        if (itheory .eq. 1) then 
         do iatom = 1, natoms
          dusr(:,iatom) = dusr(:,iatom) - (eq2/2.0d0)*fewald(:,iatom)
         end do
        end if
  
 
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
          if (itheory .eq. 1) u0tot = u0tot + corksr(iatom,ineigh)
         end do
        end do
 
! Notice that we add the corksr correction to etot.
! This is because of the subtraction minus discussed above.
        eklr = 0.0d0
        if (itheory .eq. 1) then
         do iatom = 1, natoms
          do jatom = iatom, natoms
           
! Calculate q(iatom)*q(jatom) - q0(iatom)*q0(jatom) = QQ
           QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
           if (iatom .eq. jatom) then
            eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ
           else
            eklr = eklr + (eq2/2.0d0)*ewald(iatom,jatom)*QQ +           &
     &                    (eq2/2.0d0)*ewald(jatom,iatom)*QQ
           end if
          end do
         end do
         u0tot = u0tot - eklr
        end if
        uiiuee = u0tot - ue0tot


! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_usr
