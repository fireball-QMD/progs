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


! setterp_2d.f90
! Program Description
! ==========================================================================
!       This subroutine must be called in main before you interpolate 
! anything.
!
! ==========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ==========================================================================

! Program Declaration
! ===========================================================================
        subroutine setterp_2d (nspecies, itheory_xc, itheory)
        use dimensions
        use integrals
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: itheory

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer in2
        integer in3
        integer index
        integer isorp

        real xmin
        real ymin
 
! Procedure
! ===========================================================================
        write (*,*) ' Running setterp_2d. Set up two-dimensional interpolator. '
 
! Table some trivial things that otherwise get recomputed at every possible
! opportunity.

! For now, we assume ALWAYS that xmin = 0.0d0
        xmin = 0.0d0
        ymin = 0.0d0
        do in1 = 1, nspecies
         do in2 = 1, nspecies
          do in3 = 1, nspecies
           index = icon3c(in1,in2,in3)
           do isorp = 0, isorpmax
            hx_bcna(isorp,index) = (x3cmax_bcna(isorp,index) - xmin)         &
     &                            /(numx3c_bcna(isorp,index) - 1)
            hy_bcna(isorp,index) = (y3cmax_bcna(isorp,index) - ymin)         &
     &                            /(numy3c_bcna(isorp,index) - 1)
           end do
           if (itheory .ne. 3) then 
            if(itheory_xc .eq. 0) then
              do isorp = 0, ideriv_max
                 hx_xc3c(isorp,index) = (x3cmax_xc3c(isorp,index) - xmin)  &
     &                            /real(numx3c_xc3c(isorp,index) - 1)
                 hy_xc3c(isorp,index) = (y3cmax_xc3c(isorp,index) - ymin)  &
     &                            /real(numy3c_xc3c(isorp,index) - 1)
              end do
            else if(itheory_xc .ne. 3) then
              do isorp = 1, isorpmax_xc
                 hx_den3(isorp,index) = (x3cmax_den3(isorp,index) - xmin)   &
     &                            /real(numx3c_den3(isorp,index) - 1)
                 hy_den3(isorp,index) = (y3cmax_den3(isorp,index) - ymin)   &
     &                            /real(numy3c_den3(isorp,index) - 1)
              end do
            end if
           endif

          end do   ! end do in3
         end do  ! end do in2 
        end do ! end do in1
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine setterp_2d
