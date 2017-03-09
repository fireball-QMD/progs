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

 
! assemble_hxc_usr.f90
! Program Description
! ===========================================================================
! This routine also computes the xc double counting correction. First, we
! get the data for the atoms in1, in2 at r1,r2. This call will return five
! values: one for the neutral pair, and four with some predetermined charge
! transfer (dq of the charged shell for xc, set in CREATOR: e.g., Si.inc).
! We then interpolate for the current charge distribution of the pair.
! We calculate:
!   (n1+n2)*(exc(1+2)-muxc(1+2)) - n1*(exc(1)-xcmu(1))
!                                - n2*(exc(2)-xcmu(2))
 
! This routine computes derivatives only if iforce = 1. This routine also
! computes the force derivative with respect to ratom of the short-ranged
! energy per cell, thus dusr(3,iatom) = - d/d(ratom(3,iatom)) usr.
! Here ratom is the basis atom position in the central cell. The minus sign
! makes it force-like.
!
! ===========================================================================
! Code written by:
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
        subroutine assemble_hxc_usr (natoms, itheory, iforce, uxcdcc)
        use configuration
        use charges
        use constants_fireball
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
        integer, intent (in) :: itheory
        integer, intent (in) :: natoms
 
! Output
        real, intent (out) :: uxcdcc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ideriv
        integer in1
        integer in2
        integer ineigh
        integer interaction
        integer issh
        integer jatom
        integer mbeta
        integer non2c
 
        real distance
        real dq1
        real dq2
        real dqi
        real dqj
        real dxc
        real dxc00, dxc0P, dxc0M, dxcP0, dxcM0
        real xc
        real xc00, xc0P, xc0M, xcP0, xcM0
 
        real, dimension (3) :: eta
        real, dimension (natoms) :: Q
        real, dimension (natoms) :: Q0
        real, dimension (3) :: r1
        real, dimension (3) :: r2
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' Welcome to assemble_hxc_usr.f90! '
        write (*,*) '  '
 
! Initialize arrays
        dxcv = 0.0d0
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
        do iatom = 1, natoms
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Determine dqi and dq1:
         dqi = Q(iatom) - Q0(iatom)
         dq1 = dq(in1)

! Loop over all neighbors ineigh of iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = xl(:,mbeta) + ratom(:,jatom)
          in2 = imass(jatom)

! Determine dqj and dq2:
          dqj = Q(jatom) - Q0(jatom)
          dq2 = dq(in2)
 
! Calculate the distance between the two atoms.
          distance = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2            &
     &                                       + (r2(3) - r1(3))**2)
          if (iatom .eq. jatom .and. mbeta .eq. 0) then

! SPECIAL CASE: SELF-INTERACTION

          else        
 
! ****************************************************************************
!
! XC DOUBLE COUNTING CORRECTION
! ****************************************************************************
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
           call interpolate_1d (interaction, ideriv, in1, in2, non2c, iforce,&
     &                          distance, xc00, dxc00)
 
           if (itheory .eq. 1) then
            ideriv = 1
            call interpolate_1d (interaction, ideriv, in1, in2, non2c,       &
     &                           iforce, distance, xcM0, dxcM0)
            ideriv = 2
            call interpolate_1d (interaction, ideriv, in1, in2, non2c,       &
     &                           iforce, distance, xcP0, dxcP0)
            ideriv = 3
            call interpolate_1d (interaction, ideriv, in1, in2, non2c,       &
     &                           iforce, distance, xc0M, dxc0M)
            ideriv = 4
            call interpolate_1d (interaction, ideriv, in1, in2, non2c,       &
     &                           iforce, distance, xc0P, dxc0P)
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
             xc = xc + ((xcP0 - xc00)/dq1)*dqi + ((xc0P - xc00)/dq2)*dqj
             dxc = dxc + ((dxcP0 - dxc00)/dq1)*dqi + ((dxc0P - dxc00)/dq2)*dqj
            end if
 
! (-,-) case II
            if (dqi .lt. 0.0d0 .and. dqj .lt. 0.0d0) then
             xc = xc + ((xc00 - xcM0)/dq1)*dqi + ((xc00 - xc0M)/dq2)*dqj
             dxc = dxc + ((dxc00 - dxcM0)/dq1)*dqi + ((dxc00 - dxc0M)/dq2)*dqj
            end if
 
! (+,-) case III
            if (dqi .gt. 0.0d0 .and. dqj .lt. 0.0d0) then
             xc = xc + ((xcP0 - xc00)/dq1)*dqi + ((xc00 - xc0M)/dq2)*dqj
             dxc = dxc + ((dxcP0 - dxc00)/dq1)*dqi + ((dxc00 - dxc0M)/dq2)*dqj
            end if
 
! (-,+) case IV
            if (dqi .lt. 0.0d0 .and. dqj .ge. 0.0d0) then
             xc = xc + ((xc00 - xcM0)/dq1)*dqi + ((xc0P - xc00)/dq2)*dqj
             dxc = dxc + ((dxc00 - dxcM0)/dq1)*dqi + ((dxc0P - dxc00)/dq2)*dqj
            end if
           end if
 
! Now we add the contribution to the total. Notice the one half factor, it is
! the sum over (iatom,jatom) with iatom not equal to jatom.
           uxcdcc = uxcdcc + xc/2.0d0
 
! ***************************************************************************
!
!                                FORCES
! ***************************************************************************
           if (iforce .eq. 1) then
            eta = (r2 - r1)/distance
 
! XC-DOUBLE COUNTING FORCE:
            dxcv(:,iatom) = dxcv(:,iatom) + eta(:)*dxc/2.0d0
            dxcv(:,jatom) = dxcv(:,jatom) - eta(:)*dxc/2.0d0
           end if
          end if
 
! End of loop over neighbors
         end do
 
! End of loop over iatom
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
