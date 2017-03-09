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
        subroutine trescentros (interaction, isorp, maxtype, in1, in2, indna, &
     &                           x, y, cost, eps, bcnax, nspecies)
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
        integer, intent (in) :: interaction
        integer, intent (in) :: isorp
        integer, intent (in) :: maxtype
        integer, intent (in) :: nspecies

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y
 
        real, intent (in), dimension (3, 3) :: eps

! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcnax
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer iME
        integer index
        integer inu
        integer kforce
        integer nl

        real argument
        real cost2
        real dQ_Ldx
        real dQ_Ldy
        real Q_L
        real sint
        integer nx, ny
        real xxmax, yymax
        real hx, hy
 
        real, dimension (0:ntheta - 1, ME3c_max) :: bcnalist
        real, dimension (numorb_max, numorb_max) :: bcnam
        real, dimension (ME3c_max) :: hlist
        real, dimension (0:ntheta - 1) :: p
 
! Procedure
! ===========================================================================
        kforce = 0

! Now interpolate.
! This subroutine calls the subroutine intrp1d as needed to find the value
! of the matrix elements for any given atomic configuration.
        index = icon3c(in1,in2,indna)
        if (interaction .eq. 1) then

         hx = hx_bcna(isorp,index)
         hy = hy_bcna(isorp,index)
         nx = numx3c_bcna(isorp,index)
         ny = numy3c_bcna(isorp,index)
         xxmax = x3cmax_bcna(isorp,index)
         yymax = y3cmax_bcna(isorp,index)

         if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
           write (*,*) ' What the heck is going on in trescentros!!! error!!! '
           write (*,*) ' x = ', x, ' Max of data = ', xxmax
           write (*,*) ' y = ', y, ' Max of data = ', yymax
           stop
         end if

         do iME = 1, index_max3c(in1,in2)
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         bcna_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(0,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         bcna_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(1,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         bcna_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(2,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         bcna_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(3,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         bcna_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(4,iME) = Q_L
         end do

        else if (interaction .eq. 2) then

         hx = hx_xc3c(isorp,index)
         hy = hy_xc3c(isorp,index)
         nx = numx3c_xc3c(isorp,index)
         ny = numy3c_xc3c(isorp,index)
         xxmax = x3cmax_xc3c(isorp,index)
         yymax = y3cmax_xc3c(isorp,index)

         if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
           write (*,*) ' What the heck is going on in trescentros!!! error!!! ' 
           write (*,*) ' x = ', x, ' Max of data = ', xxmax
           write (*,*) ' y = ', y, ' Max of data = ', yymax
           stop
         end if

         do iME = 1, index_max3c(in1,in2)
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         xc3c_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(0,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(1,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         xc3c_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(2,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         xc3c_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(3,iME) = Q_L

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         xc3c_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(4,iME) = Q_L
         end do

        else if (interaction .eq. 3) then

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

         do iME = 1, index_max3c(in1,in2)
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(0,iME) = Q_L
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                         den3_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(1,iME) = Q_L
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(2,iME) = Q_L
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(3,iME) = Q_L
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(4,iME) = Q_L
         end do
        end if

! Now calculate the matrix elements in molecular coordinates.
! The variable cost is passed from the assembler.
        cost2 = cost**2
        argument = 1.0d0 - cost2
        if (argument .lt. 1.0d-5) argument = 1.0d-5
        sint = sqrt(argument) 

! Set up Legendre polys in cos(theta).
        if (ntheta .ne. 5) then
          write(*,*) ' ntheta must be 5, but it is ',ntheta
          stop     
        end if
        p(0) = 1.0d0
        p(1) = cost
        p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
        p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
        p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
 
! Multiply bcnalist pieces by the appropriate Legendre polynomials and
! form hlist.
        do iME = 1, index_max3c(in1,in2)
         hlist(iME) = 0.0d0
         do nl = 0, ntheta - 1
          hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
         end do
         if (mvalue(iME,in1,in2) .eq. 1) then
          hlist(iME) = hlist(iME)*sint
         end if
        end do
 
! Now recover bcnam which is a two-dimensional array from hlist which
! is only one-dimensional.
        call recover_3c (in1, in2, hlist, bcnam)

! Rotate bcnam into crystal-coordinates: bcnam => bcnax
        call rotate_fb (in1, in2, eps, bcnam, bcnax)
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine trescentros
