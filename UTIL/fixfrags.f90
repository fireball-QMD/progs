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


! fixfrags.f90
! Program Description
! ===========================================================================
!       This takes a offset or velocity and remove internal movements within
! fragments that are fixed. 
! jel-fr
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
! Modified by  P.Jelinek (2003)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
         subroutine fixfrags()
         use dimensions
         use fragments
         use configuration
         implicit none

! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================
         real, parameter :: tiny = 0.0000001  ! cutoff

! Local Variable Declaration and Description
! ===========================================================================
         integer ifrag     ! loop over fragments
         integer imu       ! loop index
         integer ix        ! loop index
         integer num       ! total atoms in fragment
         integer iatom     ! atom of the fragment

         real rcmx
         real rcmy         ! center
         real rcmz
         real rx,ry,rz,rr  ! used in torque part of code
         real temp
         real torquevx
         real torquevy     ! total velocity torques
         real torquevz
         real totalvx
         real totalvy      ! total velocity translation
         real totalvz
         real xmasstot

         real coster, sinister,phi

! Procedure
! ===========================================================================
         if(ifrags .ne. 0) then 
!        recover inital geometry (avoid slow numerical drift)
            ratom_frag = ratom_frag_save

! loop over fragments
            do ifrag = 1,numfrags
              num = fragsize(ifrag)
! calculate total translation of fragment center of mass
              rcmx = 0.0d0
              rcmy = 0.0d0
              rcmz = 0.0d0
              xmasstot = 0.0d0
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 xmasstot = xmasstot + temp
                 rcmx = rcmx + (xdot(0,1,iatom)-ratom_frag(1,iatom))*temp
                 rcmy = rcmy + (xdot(0,2,iatom)-ratom_frag(2,iatom))*temp
                 rcmz = rcmz + (xdot(0,3,iatom)-ratom_frag(3,iatom))*temp
              end do
              rcmx = rcmx / xmasstot
              rcmy = rcmy / xmasstot
              rcmz = rcmz / xmasstot
! update ratom_frag so that centers of old and new are the same
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 ratom_frag(1,iatom) = ratom_frag(1,iatom) + rcmx
                 ratom_frag(2,iatom) = ratom_frag(2,iatom) + rcmy
                 ratom_frag(3,iatom) = ratom_frag(3,iatom) + rcmz
              end do

! calculate new center of mass
              rcmx = 0.0d0
              rcmy = 0.0d0
              rcmz = 0.0d0
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 rcmx = rcmx + xdot(0,1,iatom)*temp
                 rcmy = rcmy + xdot(0,2,iatom)*temp
                 rcmz = rcmz + xdot(0,3,iatom)*temp
              end do
              rcmx = rcmx / xmasstot
              rcmy = rcmy / xmasstot
              rcmz = rcmz / xmasstot

! calculate net translational velocity of the center of mass
              totalvx = 0.0d0
              totalvy = 0.0d0
              totalvz = 0.0d0
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 totalvx = totalvx + xdot(1,1,iatom)*temp
                 totalvy = totalvy + xdot(1,2,iatom)*temp
                 totalvz = totalvz + xdot(1,3,iatom)*temp
              end do
              totalvx = totalvx / xmasstot
              totalvy = totalvy / xmasstot
              totalvz = totalvz / xmasstot

! project out
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 xdot(1,1,iatom) = xdot(1,1,iatom) - totalvx
                 xdot(1,2,iatom) = xdot(1,2,iatom) - totalvy
                 xdot(1,3,iatom) = xdot(1,3,iatom) - totalvz
              end do

! calculate net rotational velocity
              torquevx = 0.0d0
              torquevy = 0.0d0
              torquevz = 0.0d0
              do imu = 1,num
                 iatom=fragatm(imu,ifrag)
                 temp=xmass(iatom)
                 rx=(xdot(0,1,iatom)-rcmx)
                 ry=(xdot(0,2,iatom)-rcmy)
                 rz=(xdot(0,3,iatom)-rcmz)
                 rr=rx**2+ry**2+rz**2
                 if(rr .gt. tiny)then
                    torquevx = torquevx +                                   &
    &                  (ry*xdot(1,3,iatom) - rz*xdot(1,2,iatom))*temp/rr
                    torquevy = torquevy +                                   &
    &                  (rz*xdot(1,1,iatom) - rx*xdot(1,3,iatom))*temp/rr
                    torquevz = torquevz +                                   &
    &                  (rx*xdot(1,2,iatom) - ry*xdot(1,1,iatom))*temp/rr
                 end if
              end do
              torquevx = torquevx / xmasstot
              torquevy = torquevy / xmasstot
              torquevz = torquevz / xmasstot

! find net rotation about z axis
              phi = 0.0d0 
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 rx = xdot(0,1,iatom)-rcmx
                 ry = xdot(0,2,iatom)-rcmy
                 rr = rx**2+ry**2
                 if(rr .gt. tiny)then
                    phi = phi + sign(acos(rx/sqrt(rr)),ry)*temp
                 end if
                 rx = ratom_frag(1,iatom)-rcmx
                 ry = ratom_frag(2,iatom)-rcmy
                 rr=rx**2+ry**2
                 if(rr .gt. tiny)then
                    phi = phi - sign(acos(rx/sqrt(rr)),ry)*temp
                 end if
              end do
              phi = phi / xmasstot

! move ratom_frag by this
              coster = cos(phi)
              sinister = sin(phi)
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 rx = ratom_frag(1,iatom) - rcmx
                 ry = ratom_frag(2,iatom) - rcmy
                 ratom_frag(1,iatom) = (rx*coster-ry*sinister) + rcmx
                 ratom_frag(2,iatom) = (ry*coster+rx*sinister) + rcmy
              end do

! do rotation about y-axis
              phi = 0.0d0
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 rz = xdot(0,3,iatom) - rcmz
                 rx = xdot(0,1,iatom) - rcmx
                 rr = rz**2 + rx**2
                 if(rr .gt. tiny)then
                    phi = phi + sign(acos(rz/sqrt(rr)),rx)*temp
                 end if
                 rz = ratom_frag(3,iatom) - rcmz
                 rx = ratom_frag(1,iatom) - rcmx
                 rr = rz**2 + rx**2
                 if(rr .gt. tiny)then
                    phi = phi - sign(acos(rz/sqrt(rr)),rx)*temp
                 end if
              end do
              phi = phi / xmasstot

              coster = cos(phi)
              sinister = sin(phi)
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 rz = ratom_frag(3,iatom) - rcmz
                 rx = ratom_frag(1,iatom) - rcmx
                 ratom_frag(3,iatom) = (rz*coster - rx*sinister)+rcmz
                 ratom_frag(1,iatom) = (rx*coster + rz*sinister)+rcmx
              end do
              
! do rotation about x-axis
              phi = 0.0d0
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 ry = xdot(0,2,iatom) - rcmy
                 rz = xdot(0,3,iatom) - rcmz
                 rr = ry**2 + rz**2
                 if(rr .gt. tiny)then
                    phi = phi + sign(acos(ry/sqrt(rr)),rz)*temp
                 end if
                 ry = ratom_frag(2,iatom) - rcmy
                 rz = ratom_frag(3,iatom) - rcmz
                 rr = ry**2 + rz**2
                 if(rr .gt. tiny)then
                    phi = phi - sign(acos(ry/sqrt(rr)),rz)*temp
                 end if
              end do
              phi = phi / xmasstot
              
              coster = cos(phi)
              sinister = sin(phi)
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 ry = ratom_frag(2,iatom) - rcmy
                 rz = ratom_frag(3,iatom) - rcmz
                 ratom_frag(2,iatom) = (ry*coster-rz*sinister) + rcmy
                 ratom_frag(3,iatom) = (rz*coster+ry*sinister) + rcmz
              end do
              
! copy ratom_frags to xdot(0,*,*)
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 xdot(0,1,iatom) = ratom_frag(1,iatom)
                 xdot(0,2,iatom) = ratom_frag(2,iatom)
                 xdot(0,3,iatom) = ratom_frag(3,iatom)
              end do

! put net velocity into xdot(1,*,*)
! set velocities equal to net velocity translation
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 xdot(1,1,iatom) = totalvx
                 xdot(1,2,iatom) = totalvy
                 xdot(1,3,iatom) = totalvz
              end do
! add in rotation
              do imu = 1,num
                 iatom = fragatm(imu,ifrag)
                 temp = xmass(iatom)
                 rx = (ratom_frag(1,iatom) - rcmx)
                 ry = (ratom_frag(2,iatom) - rcmy)
                 rz = (ratom_frag(3,iatom) - rcmz)
                 xdot(1,2,iatom) = xdot(1,2,iatom) - torquevx*rz
                 xdot(1,3,iatom) = xdot(1,3,iatom) + torquevx*ry
                 xdot(1,3,iatom) = xdot(1,3,iatom) - torquevy*rx
                 xdot(1,1,iatom) = xdot(1,1,iatom) + torquevy*rz
                 xdot(1,1,iatom) = xdot(1,1,iatom) - torquevz*ry
                 xdot(1,2,iatom) = xdot(1,2,iatom) + torquevz*rx
              end do

           enddo
        else
! fix selected atoms
! loop over fragments
           do ifrag=1,numfrags
              num=fragsize(ifrag)
! set velocities equal to zero
              do imu=1,num
                 iatom=fragatm(imu,ifrag)
                 do ix = 1,3
                   if(fragxyz(ix,iatom) .eq. 1) xdot(1,ix,iatom) = 0.0d0
                 end do
              end do
           enddo
        endif
! Format Statements
! ===========================================================================
         return
       end subroutine fixfrags

