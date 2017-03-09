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


! epsilon.f90
! Program Description
! ===========================================================================
! input: r1,r2
! output: epsilon the metric tensor
!
! note: the third column of epsilon is eta(3)
!
! spe=epsilon backwards.
!
! R1vector points toward the point O while R2vector points
!   away from the point O.
!
!                         *O    (XP,YP,ZP)
!                      *    *
!       R1VECTOR    *        *  R2VECTOR
!                *            *
!             *                *
! (X,Y,Z)  *                    *
!       *
!
!            |  ^     ^       ^     ^        ^     ^   |
!            |  X-dot-XP      X-dot-YP       X-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!      spe = |  Y-dot-XP      Y-dot-YP       Y-dot-ZP  |
!            |                                         |
!            |  ^     ^       ^     ^        ^     ^   |
!            |  Z-dot-XP      Z-dot-YP       Z-dot-ZP  |
!            |                                         |
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
         subroutine epsilon(R1,R2,spe)
         implicit none
 
! Argument Declaration and Description
! ===========================================================================
         real, intent(in) :: r1(3)
         real, intent(in) :: r2(3)
         real, intent(out) :: spe(3,3)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
         integer i,j,ii,jj,kk,ix
         real r1mag,r2mag,ypmag,unit
         real XPHAT(3),YPHAT(3),ZPHAT(3),R1HAT(3)
 
! Procedure
! ===========================================================================
        r1mag=sqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
        r2mag=sqrt(r2(1)*r2(1)+r2(2)*r2(2)+r2(3)*r2(3))
!
        if (r2mag .lt. 1.0d-4) then
!        r2vector = (0,0,0) ----- set eps = unit matrix
         do i=1,3
          do j=1,3
           spe(i,j)=0.0e0
          end do
          spe(i,i)=1.0e0
         end do
         write(*,*) '  r2=0 : spe = unit matrix!'
         return
        end if

! zphat lies along r2vector
        zphat(1)=r2(1)/r2mag
        zphat(2)=r2(2)/r2mag
        zphat(3)=r2(3)/r2mag

! yphat = zphat-cross-r1hat
        if(r1mag.gt.1.0d-4)then
          r1hat(1)=r1(1)/r1mag
          r1hat(2)=r1(2)/r1mag
          r1hat(3)=r1(3)/r1mag
          yphat(1)=zphat(2)*r1hat(3)-zphat(3)*r1hat(2)
          yphat(2)=zphat(3)*r1hat(1)-zphat(1)*r1hat(3)
          yphat(3)=zphat(1)*r1hat(2)-zphat(2)*r1hat(1)
          ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
          if(ypmag.gt.0.000001)goto 3000
        end if
!
! zphat and r1hat are colinear or r1vector=(0,0,0)
!   find the first non-zero component of zphat
!
          if(abs(zphat(1)).gt.1.0d-4)then
! zphat(1) not equal to zero
            yphat(1)=-(zphat(2)+zphat(3))/zphat(1)
            yphat(2)=1.0e0
            yphat(3)=1.0e0
            ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
          else if(abs(zphat(2)).gt.1.0d-4)then
! zphat(2) not equal to zero
            yphat(1)=1.0e0
            yphat(2)=-(zphat(1)+zphat(3))/zphat(2)
            yphat(3)=1.0e0
            ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
          else
! zphat(3) not equal to zero
            yphat(1)=1.0e0
            yphat(2)=1.0e0
            yphat(3)=-(zphat(1)+zphat(2))/zphat(3)
            ypmag=sqrt(yphat(1)*yphat(1)+yphat(2)*yphat(2)+yphat(3)*yphat(3))
          end if

 3000   continue

        yphat(1)=yphat(1)/ypmag
        yphat(2)=yphat(2)/ypmag
        yphat(3)=yphat(3)/ypmag
! find pihat
        xphat(1)=yphat(2)*zphat(3)-yphat(3)*zphat(2)
        xphat(2)=yphat(3)*zphat(1)-yphat(1)*zphat(3)
        xphat(3)=yphat(1)*zphat(2)-yphat(2)*zphat(1)
!
! find epsilon matrix
        do ix=1,3
          spe(ix,1)=xphat(ix)
          spe(ix,2)=yphat(ix)
          spe(ix,3)=zphat(ix)
        end do
 
! test by computing spe*spe(dagger)
        do ii=1,3
          do jj=1,3
            unit=0.0e0
            do kk=1,3
              unit=unit+spe(ii,kk)*spe(jj,kk)
            end do
            if(ii.eq.jj.and.abs(unit-1.e0).gt.1.0d-4)then
              write(*,*) ' ******* error in  epsiln 1=spe*spedag=',unit
            end if
            if(ii.ne.jj.and.abs(unit).gt.1.0d-4)then
              write(*,*) ' ******* error in  epsiln 0=spe*spedag=',unit
            end if
          end do
        end do
        return
 
! format statements
! ===========================================================================
        return
        end
