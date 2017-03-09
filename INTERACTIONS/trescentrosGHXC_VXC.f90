! copyright info:
!
!                             @Copyright 1998
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!            Jun Wang
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
!
! fireball-qmd is a free (GPLv3) open project.

!
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

! tresgaussians.f90
!
! Adapted from trescentros. This gets the same three-center matrix elements 
! gotten by trescentros but uses gaussian fits to wavefunctions and 
! na potential.
!
! Program Description for trescentros
! ===========================================================================
!       This subroutine calculates the (three-center) matrix elements (mu,nu).
! It calculates the matrix elements of a bondcharge with a NA interaction,
! using data which was stored in - threecint(nz,nx,kbc,jxx,ithet,ind1,ind2,ind3)
! This subroutine calculates the matrix elements molecular and crystal
! coordinates of a bondcharge with the neutral atom elsewhere (i.e. neither an
! "ontop" nor an atom-neutral atom).
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
! The KEY for isorp for x/c 3-center. Evaluate the spherically averaged wave
! functions for all orbitals of the two atoms.
! (x y z) = charge on (1 2 3).       1=Left  2=right  3=NA.
!
! isopr = 0    Neutral  (0 0 0)
!         1             (- 0 0)
!         2             (+ 0 0)
!         3             (0 - 0)
!         4             (0 + 0)
!         5             (0 0 -)
!         6             (0 0 +)
!
!
! ymin,ymax,numy: bond charge  distances grid
! xmin,xmax,numx: neutral atom distances grid
!
! Here are the data file characteristics:
!
! ===========================================================================
! Code originally written by Jose Ortega and Juergen Fritsch
 
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E. RM Dock
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
!
! Code re-rewritten for gaussians by:
! John Tomfohr
! john.tomfohr@asu.edu
! ===========================================================================
!
! Program Declaration
! ===========================================================================
! WANG gauss, XC
! Added two indexes ish00 and n_center
! ish00 = 0  for NA, the output bcnax is 3C-neutral atom potential matrix elements
! ish00 = -1 for XC, the output bcnax is 3C-electron density matrix elements
!
! n_center is control parameter for XC
! n_center = 2 for two center case
! n_center = 3 for three center case    
!
        subroutine trescentrosGHXC_VXC (in1, in2, indna, x, y, cost,    &
     &                                  eps, bcnax, rcutoff)
   
        use dimensions
        use interactions
        use gaussG
        use constants_fireball
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y
 
        real, intent (in), dimension (3, 3) :: eps

! JKT step. added rcutoff
        real, intent (in), dimension (0:nspec_max-1, nsh_max) :: rcutoff
 
! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcnax
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iindex
        integer jindex
        integer l1
        integer l2

        integer ish
        integer jsh
        integer im
        integer jm
        integer imu
        integer inu

        real rmatel
        real, dimension (numorb_max, numorb_max) :: bcnam

! JKT step
        real rc1
        real rc2
        real rc3

        real r12
        real r23
        real r31

        real sint

        real ptemp12
        real ptemp23
        real ptemp31

        real fermi12
        real fermi23
        real fermi31

! Procedure
! ===========================================================================
! Initialize bcnam, bcnax and smatel
        bcnam = 0.0
        bcnax = 0.0

! JKT step. why is rcutoff indexed this way?! (0:nspec-1)?!!!
! JKT step  get rc3...
        rc3=0.
        do ish=1,nssh(indna)
        if (rcutoff(indna-1,ish).gt.rc3) then
        rc3=rcutoff(indna-1,ish)
        end if
        end do        

! JKT step. get r12,r23,r31
        sint=sqrt(abs(1.-cost**2))
        r12=y
        r23=sqrt((y/2.-x*cost)**2+(x*sint)**2)
        r31=sqrt((y/2.+x*cost)**2+(x*sint)**2)

! JKT step. get ptemps. (pseudo-temperatures)
        ptemp12=ptemp*r12
        ptemp23=ptemp*r23
        ptemp31=ptemp*r31

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

        call gelements_VXC (in1, in2, indna, l1, l2, ish, jsh, im, jm,  &
     &                      x, y, cost, rmatel)

! endif of symmetry

         fermi12=2./(exp((r12-(rc1+rc2))/ptemp12)+1)-1.
         fermi23=2./(exp((r23-(rc2+rc3))/ptemp23)+1)-1.
         fermi31=2./(exp((r31-(rc3+rc1))/ptemp31)+1)-1.
         rmatel=rmatel*fermi12*fermi23*fermi31
         endif
!
         bcnam(iindex,jindex)=rmatel
!
! im, jm
         end do
         end do
! ish, jsh
         end do
         end do

! Rotate bcnam into crystal-coordinates: bcnam => bcnax

        call rotate_fb (in1, in2, eps, bcnam,  bcnax)

        return
        end

! ===========================================================================
