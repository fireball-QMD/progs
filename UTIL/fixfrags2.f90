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


! fixfrags2.f90
! Program Description
! ===========================================================================
!       This takes a force and removes internal movements within fragments
! that are fixed.
! 
! Eventually selected atom will be frozen in the initial position 
! (option ifrags = 0). ifrags is first variable read in FRAGMENTS file.  
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
!
! Modified by P.Jelinek (2003)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
         subroutine fixfrags2(ftot)
         use dimensions
         use configuration
         use fragments
         implicit none

! Argument Declaration and Description
! ===========================================================================
         real, intent(inout) :: ftot (3,natoms)      ! force

! Local Parameters and Data Declaration
! ===========================================================================
         real, parameter :: tiny = 0.0000001  ! cutoff

! Local Variable Declaration and Description
! ===========================================================================
         integer ifrag     ! loop over fragments
         integer imu       ! loop index
         integer num       ! total atoms in fragment
         integer iatom     ! atom of the fragment
         integer ix        ! loop index

         real rcmx
         real rcmy         ! center of mass
         real rcmz
         real temp         ! mass of a given atom
         real rx,ry,rz,rr  ! used in torque part of code
         real torquex
         real torquey      ! total torques
         real torquez
         real totalx
         real totaly       ! total translation
         real totalz
         real xmasstot     ! total mass of fragment

! Procedure
! ===========================================================================
         if(ifrags .ne. 0) then 
! loop over fragments
           do ifrag=1,numfrags
             num=fragsize(ifrag)
! calculate total translational force of fragments center of mass
             totalx=0
             totaly=0
             totalz=0
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                totalx=totalx+ftot(1,iatom)
                totaly=totaly+ftot(2,iatom)
                totalz=totalz+ftot(3,iatom)
             end do

! calculate rcm (center-of-mass)
             xmasstot=0
             rcmx=0
             rcmy=0
             rcmz=0
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                temp=xmass(iatom)
                xmasstot=xmasstot+temp
                rcmx=rcmx+temp*ratom(1,iatom)
                rcmy=rcmy+temp*ratom(2,iatom)
                rcmz=rcmz+temp*ratom(3,iatom)
             end do
             rcmx=rcmx/xmasstot
             rcmy=rcmy/xmasstot
             rcmz=rcmz/xmasstot
             totalx=totalx/xmasstot
             totaly=totaly/xmasstot
             totalz=totalz/xmasstot

! project out net translation
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                temp=xmass(iatom)
                ftot(1,iatom)=ftot(1,iatom)-totalx*temp
                ftot(2,iatom)=ftot(2,iatom)-totaly*temp
                ftot(3,iatom)=ftot(3,iatom)-totalz*temp
             end do

! calculate net torque
             torquex=0
             torquey=0
             torquez=0
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                temp=xmass(iatom)
                rx=(ratom(1,iatom)-rcmx)*temp
                ry=(ratom(2,iatom)-rcmy)*temp
                rz=(ratom(3,iatom)-rcmz)*temp
                rr=rx**2+ry**2+rz**2
                if(rr .gt. tiny)then
                   rr=1.0/rr
                   torquex=torquex+(ry*ftot(3,iatom)-rz*ftot(2,iatom))*rr
                   torquey=torquey+(rz*ftot(1,iatom)-rx*ftot(3,iatom))*rr
                   torquez=torquez+(rx*ftot(2,iatom)-ry*ftot(1,iatom))*rr
                end if
             end do
             torquex=torquex/num
             torquey=torquey/num
             torquez=torquez/num

! set ftots equal to translation (need larger force on heavier atoms)
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                temp=xmass(iatom)
                ftot(1,iatom)=totalx*temp
                ftot(2,iatom)=totaly*temp
                ftot(3,iatom)=totalz*temp
             end do

! add in rotation
             
             do imu=1,num
                iatom=fragatm(imu,ifrag)
                temp=xmass(iatom)
                rx=(ratom(1,iatom)-rcmx)*temp
                ry=(ratom(2,iatom)-rcmy)*temp
                rz=(ratom(3,iatom)-rcmz)*temp
                ftot(2,iatom)=ftot(2,iatom)-torquex*rz
                ftot(3,iatom)=ftot(3,iatom)+torquex*ry
                ftot(3,iatom)=ftot(3,iatom)-torquey*rx
                ftot(1,iatom)=ftot(1,iatom)+torquey*rz
                ftot(1,iatom)=ftot(1,iatom)-torquez*ry
                ftot(2,iatom)=ftot(2,iatom)+torquez*rx
             end do

          enddo !di ifrag 
       else
! remove forces on selected atoms
! loop over fragments
           do ifrag=1,numfrags
              num=fragsize(ifrag)
! set forces equal to zero
              do imu=1,num
                 iatom=fragatm(imu,ifrag)
                 do ix = 1,3
                   if(fragxyz(ix,iatom) .eq. 1) ftot(ix,iatom) = 0.0d0
                 end do
              end do
           enddo
        endif

! Format Statements
! ===========================================================================
         return
       end subroutine fixfrags2

