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


! trescentros.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the (three-center) matrix elements (mu,nu).
! It calculates the matrix elements of a bondcharge with a NA interaction,
! using data which was stored in - threecint(ithet,nz,nx,kbc,jxx,ind1,ind2,ind3)
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
        subroutine trescentrosS ( isorp, maxtype, in1, in2, indna, x, y,    &
     &                            cost, eps, bcnax, nspecies)
        use dimensions
        use interactions
        use integrals
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna
        integer, intent (in) :: isorp
        integer, intent (in) :: maxtype
        integer, intent (in) :: nspecies

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y
 
        real, intent (in), dimension (3, 3) :: eps

! Output
        real, intent (out), dimension (nsh_max, nsh_max) :: bcnax
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer iME
        integer index
        integer inu
        integer kforce
        integer nl
        integer nx
        integer ny

        real argument
        real cost2
        real dQ_Ldx
        real dQ_Ldy
        real Q_L
        real sint
        real xxmax
        real yymax
        real hx
        real hy
 
        real, dimension (0:ithetamax - 1, MES_max) :: bcnalist
        real, dimension (MES_max) :: hlist
        real, dimension (0:ithetamax - 1) :: p
 
! Procedure
! ===========================================================================
! Initialize bcnax
        do inu = 1, nssh(in2)
         do imu = 1, nssh(in1)
          bcnax(imu,inu) = 0.0d0
         end do
        end do
        kforce = 0

! Now interpolate.
! This subroutine calls the subroutine intrp1d as needed to find the value
! of the matrix elements for any given atomic configuration.
        if (ntheta .gt. ithetamax .or. ntheta .gt. 5) then
         write (*,*) ' too many thetas in trescentros.f90 '
         stop
        end if
 
        index = icon3c(in1,in2,indna)

! HAO Sep. 11, 2003, modified
         hx = hx_den3(isorp,index)
         hy = hy_den3(isorp,index)
         nx = numx3c_den3(isorp,index)
         ny = numy3c_den3(isorp,index)
         xxmax = x3cmax_den3(isorp,index)
         yymax = y3cmax_den3(isorp,index)
         if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
           write (*,*) ' What the heck is going on in trescentros!!! error!!! '
           write (*,*) ' x = ', x, ' Max of data = ', xxmax
           write (*,*) ' y = ', y, ' Max of data = ', yymax
           stop
         end if
        do iME = 1, index_maxS(in1,in2)
           call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                       den3S_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
           bcnalist(0,iME) = Q_L

           call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                       den3S_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
           bcnalist(1,iME) = Q_L

           call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                       den3S_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
           bcnalist(2,iME) = Q_L

           call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                       den3S_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
           bcnalist(3,iME) = Q_L

           call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                       den3S_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
           bcnalist(4,iME) = Q_L
        end do

! Now calculate the matrix elements in molecular coordinates.
! The variable cost is passed from the assembler.
        cost2 = cost**2
        argument = 1.0d0 - cost2 + 1.0d-05
        if (argument .ge. 0.0d0) then
         sint = sqrt(argument)
        else
         write (*,*) '  Trescentros: BAD SQRT ******* ERROR '
         write (*,*) '  *** STANDARD FIXUP, SET SINT = 0 '
         argument = 0.0d0
         sint = 0.0d0
        end if
 
! Set up Legendre polys in cos(theta).
        p(0) = 1.0d0
        p(1) = cost
        p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
        p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
        p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
 
! Multiply bcnalist pieces by the appropriate Legendre polynomials and
! form hlist.
        do iME = 1, index_maxS(in1,in2)
         hlist(iME) = 0.0d0
         do nl = 0, ntheta - 1
          hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
         end do
         if (mvalueS(iME,in1,in2) .eq. 1) then
          hlist(iME) = hlist(iME)*sint
         end if
        end do
 
! Now recover bcnax which is a two-dimensional array from hlist which
! is only one-dimensional.
        call recover_S (in1, in2, hlist, bcnax)

! Format Statements
! ===========================================================================
 
        return
      end subroutine trescentrosS
