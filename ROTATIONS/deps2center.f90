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


! deps2center.f90
! Program Description
! ===========================================================================
!       This subroutine sets up deps/dr1 in the two-center molecular system 
! defined by sighat=(r2-r1)/|r2-r1| and piprimehat=(r2-cross-r1)/|r2-cross-r1|.
! The eps matrix is obtained by calling epsiln(r1,rvec) where rvec=r2-r1.
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
        subroutine deps2cent(r1,r2,eps2,deps2)
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! input:
! r1    - position of atom 1
! r2    - posiiton of atom 2
! eps2  - 2 center epsilon with z=r2-r1
        real, intent(in) :: r1(3),r2(3),eps2(3,3)

! output:
! deps2 - deps/dr1
        real, intent(out) :: deps2(3,3,3)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i,ii,ix

        real r2mag2,r2mag,r1mag,denom
        real crossmag,dd,dot,term,ddinv,crossinv
        real crossa(3)
 
! Procedure
! ===========================================================================
        deps2=0.e0

! If we are doing an atom, r1=r2. Then set deps2 to zero.
        dd=sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2) 
        if(dd.lt.1.0d-4)return
        ddinv=1.0/dd 
 
        r2mag2= r2(1)**2 + r2(2)**2 + r2(3)**2
        r2mag=sqrt(r2mag2)
        r1mag=sqrt(r1(1)**2 + r1(2)**2 + r1(3)**2)
        crossa(1)=r2(2)*r1(3)-r2(3)*r1(2)
        crossa(2)=r2(3)*r1(1)-r2(1)*r1(3)
        crossa(3)=r2(1)*r1(2)-r2(2)*r1(1)
        crossmag=sqrt(crossa(1)**2+crossa(2)**2+crossa(3)**2)
        dot=0.e0
        do i=1,3
          dot=dot+r1(i)*r2(i)
        end do
        denom=r1mag*r2mag
 
! check to see if any atom is near origin
        if(denom.lt.1.0d-3)then
!          write(*,*) ' too close to origin'
!          write(*,*) '   r1*r2=',denom
!          write(*,*) ' setting deps2/dr1 to zero!'
          return
        endif
! check to see if atoms are colinear
        if(abs(crossmag).lt.1.0d-3)then
!          write(*,*) ' *********** warning in deps2cent'
!          write(*,*) '  r1vec and r2vec dangerously colinear'
!          write(*,*) ' r1=',r1
!          write(*,*) ' r2=',r2
!          write(*,*) '  r2vec-cross-r1vec =',crossmag
!          write(*,*) '  setting deps2/dr1 to zero'
          return
        endif
        crossinv=1.0/crossmag
!
! now calculate deps
        do ii=1,3
         do ix=1,3
          term=xlevi(ix,ii,1)*r2(1)+   &
     &         xlevi(ix,ii,2)*r2(2)+   &
     &         xlevi(ix,ii,3)*r2(3)

          deps2(ix,ii,1)=(eps2(ii,1)*eps2(ix,3)*ddinv)-(eps2(ii,1)*  &
     &     crossinv**2)*(r2mag2*r1(ix)-dot*r2(ix))+  &
     &     ddinv*crossinv*(delk(ii,ix)*(r2mag2-dot)-  &
     &     r1(ii)*r2(ix)-r2(ii)*r2(ix)+2.e0*r1(ix)*r2(ii))

          deps2(ix,ii,2)=crossinv*(term+(eps2(ii,2)*crossinv)*  &
     &        (dot*r2(ix)-r2mag2*r1(ix)))

          deps2(ix,ii,3)=-(delk(ii,ix)-eps2(ii,3)*eps2(ix,3))*ddinv
         end do
        end do
!
! Format Statements
! ===========================================================================
 
        return
        end
