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


! dosgaussians.f90
! Program Description
! ===========================================================================
!      This subroutine calculates the 2c-matrix (sx) and its derivatives (spx).
!
! ===========================================================================
! Code rewritten by:
! Hao Wang
! Department of Physics and Astronomy
! Brigham Young University
! N311 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-9270
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine dosgaussians (in1, in2, in3, y, cost, eps, deps, sx, &
     &                           spx, rcutoff)

        use dimensions
        use interactions
        use gaussG
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
 
        integer, intent (in) :: in1, in2, in3
        real, intent (in), dimension (3, 3) :: eps
        real, intent (in), dimension (3, 3, 3) :: deps
        real, intent (in), dimension (0:nspec_max-1, nsh_max) :: rcutoff
        real, intent (in) :: y    
        real, intent (in) :: cost     

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
        integer iparam

        real x
        real rmatel
        real rmately
        real rc1
        real rc2
        real rc3
        real r12
        real r23
        real r31
        real sint

        real, dimension (numorb_max, numorb_max) :: sm
        real, dimension (numorb_max, numorb_max) :: spm
        real, dimension (3, numorb_max, numorb_max) :: spmx
        real, dimension (3) :: eta

! Procedure
! ===========================================================================
! For the atom case, in3 = in1, but for everything else in3 = in2.
! For the ontop case, in2 = in1 (left) or in3 (right).
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
!
        x = 0.0
!
! JKT step. why is rcutoff indexed this way?! (0:nspec-1)?!!!
! JKT step  get rc3...
        rc3=0.
        do ish=1,nssh(in3)
        if (rcutoff(in3-1,ish).gt.rc3) then
        rc3=rcutoff(in3-1,ish)
        end if
        end do
!
!
!           X  ^                                * (NA)
!              !                             +
! X-Z PLANE    !                          +
!              !                      +   THETA               -----> Z
!              O ---------------------------------------------O
!            1 END OF BC                                OTHER END OF BC
!               MU                                          NU
!
! JKT step. get r12,r23,r31
! JKT step. get r12,r23,r31
        sint=sqrt(abs(1.-cost**2))
        r12=y
        r23=sqrt((y/2.-x*cost)**2+(x*sint)**2)
        r31=sqrt((y/2.+x*cost)**2+(x*sint)**2)


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

! symmetry, if only one of im or jm is .lt. 0 then rmatel=0. by symmetry
        if ((im.lt.0 .and. jm.ge.0) .or. (jm.lt.0 .and. im.ge.0)) then
        rmatel=0.
! JKT step added else-if
        else if (r12.gt.(rc1+rc2) .or. r23.gt.(rc2+rc3) .or.            &
     &           r31.gt.(rc3+rc1)) then
        rmatel=0.
        else

! WANG gauss, XC
! igausssignal = -1 for XC part
! Here we only calculate the integral of two-center electron density n(mu,nu)
! OFS, PRB 40, 3979 (1989)

        call gelements_VXC (in1, in2, in3, l1, l2, ish, jsh, im, jm,    &
     &                      x, y, cost, rmatel)

! endif of symmetry
         endif
!
         sm(iindex,jindex)=rmatel

! im, jm
         end do
         end do
! ish, jsh
         end do
         end do
! ====  get sm (2c-matrix) ==================================================
! Rotate sm into crystal-coordinates: sm --> sx
        call rotate_fb (in1, in3, eps, sm, sx)
 
! ***************************************************************************
! Force, two-center case

! WANG gauss-XC, since we use three center integral to deal with two-center case,
! So, we need to setup some initial values

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

! if only one of im or jm is .lt. 0 then everything's 0. by symmetry
        if ((im.lt.0 .and. jm.ge.0) .or. (jm.lt.0 .and. im.ge.0)) then
        rmately=0.
! JKT step added else if
        else if (r12.gt.(rc1+rc2) .or. r23.gt.(rc2+rc3) .or.            &
     &           r31.gt.(rc3+rc1)) then
        rmately=0.
        else

        iparam=1

        call Dgelements_VXC (l1, l2, ish, jsh, im, jm, in1, in2, in3,   &
     &                       x, y, cost, rmately, iparam)

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
         do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)

! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".

! Note the minus sign. d/dr1 = - eta * d/dd.
           if (y .gt. 1.0d-3) then
           spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
           end if
          end do
         end do

         call rotated (in1, in3, eps, deps, sm, spmx, spx)

! Format Statements
! ===========================================================================
        return
        end
! ===========================================================================
