! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
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


! doscentorsGS_overlap.f90
! Program Description
! ===========================================================================
! This subroutine calculates the 2c-matrix (sx) and its derivatives (spx).

! In this routine, gelementGS_overlap is called to get the overlap matrix element 
! - <psi|psi> using gaussian wavefactions
!
! ===========================================================================
! Code rewritten by:
! Hao Wang
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-9270
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine doscentrosG_overlap (in1,in2,y,eps,deps,sx,spx,rcutoff)
! ===========================================================================
!
        use dimensions
        use interactions
        use gaussG
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
 
        integer, intent (in) :: in1
        integer, intent (in) :: in2

        real, intent (in), dimension (3, 3) :: eps
        real, intent (in), dimension (3, 3, 3) :: deps
        real, intent (in), dimension (0:nspec_max-1, nsh_max) :: rcutoff
        real, intent (in) :: y    

! Output
! The variable sx is < phi(mu,r-r1) ! vna(r-ratm) ! phi(nu,r-r1)>
! which is a representation in molecular coordinates.

        real, intent (out), dimension (numorb_max, numorb_max) :: sx
        real, intent (out), dimension (3,numorb_max, numorb_max) :: spx
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer iindex
        integer jindex
        integer l1
        integer l2
        integer ish
        integer jsh
        integer im
        integer jm
        integer igausssignal
        integer iparam

        real rmatel
        real rmately
        real rc1
        real rc2
        real r12

        real, dimension (numorb_max, numorb_max) :: sm
        real, dimension (numorb_max, numorb_max) :: spm
        real, dimension (3, numorb_max, numorb_max) :: spmx

        real, dimension (3) :: eta

! Procedure
! ===========================================================================
! Initialize sm, scam and sx to zero.
          sm = 0.0d0
          sx = 0.0d0
         spm = 0.0d0
        spmx = 0.0d0
         spx = 0.0d0

! WANG gauss-XC
! Since we use three center integral subroutine (gelements.f90) 
! to calculate density of two-center case,
! so, we need to setup some initial values

        iindex=0
        do ish=1,nssh(in1)
        l1=lssh(ish,in1)
! JKT step rc1
        rc1=rcutoff(in1-1,ish)
        do im=-l1,l1

        iindex=iindex+1

        jindex=0
        do jsh=1,nssh(in2)
        l2=lssh(jsh,in2)
! JKT step rc2
        rc2=rcutoff(in2-1,jsh)
        do jm=-l2,l2

        jindex=jindex+1

        if (r12.gt.(rc1+rc2)) then
        rmatel=0.
        else

        call gelementsG_overlap (in1, in2, l1, l2, ish, jsh, im, jm, y, rmatel)

! endif of symmetry
         endif
!
         sm(iindex,jindex)=rmatel

         sx(iindex,jindex) = sm(iindex,jindex)
! im, jm
         end do
         end do
! ish, jsh
         end do
         end do
! ====  get sm (2c-matrix) ==================================================
!
! No rotation for spherical case

!! Rotate sm into crystal-coordinates: sm --> sx
        call rotate_fb (in1, in2, eps, sm, sx)
 
! ***************************************************************************
! Force, two-center case

        iindex=0
        do ish=1,nssh(in1)
        l1=lssh(ish,in1)
! JKT step
        rc1=rcutoff(in1-1,ish)
        do im=-l1,l1
        iindex=iindex+1
        jindex=0

        do jsh=1,nssh(in2)
        l2=lssh(jsh,in2)
! JKT step rc2
        rc2=rcutoff(in2-1,jsh)
        do jm=-l2,l2
        jindex=jindex+1

        if (r12.gt.(rc1+rc2)) then
        rmately=0.
        else

        call DgelementsG_overlap (l1, l2, in1, in2, ish, jsh, im, jm, y, rmately)

        end if

        spm(iindex,jindex)=rmately
! im, jm
        end do
        end do
! ish, jsh
        end do
        end do

! As long as epsilon1 is called with sighat in the second "spot" as
! call epsilon1(R1,sighat,spe), then eps(ix,3) = eta(ix).
         eta(:) = eps(:,3)

! First do the matrix.
         do inu = 1, nssh(in2)
          do imu = 1, nssh(in1)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
           if (y .gt. 1.0d-3) then
           spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
           spx(:,imu,inu) = spmx(:,imu,inu)
           end if
          end do
         end do

! No rotation for spherical case
         call rotated (in1, in2, eps, deps, sm, spmx, spx)

! Format Statements
! ===========================================================================
        return
        end
! ===========================================================================
