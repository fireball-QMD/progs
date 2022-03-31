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


! Dtrescentros.f90
! Program Description
! ===========================================================================
!       This routine is like trescentros.f, except that it calculates the
! forces.
!
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
! Code originally written by Jose Ortega and Juergen Fritsch.
 
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
        subroutine Dtrescentros (interaction, isorp, maxtype, in1, in2, indna,&
     &                           x, y, cost, eps, depsA, depsB, rhat, sighat, &
     &                           bcnax, f3naXa, f3naXb, f3naXc, nspecies)
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
        real, intent (in), dimension (3, 3, 3) :: depsA
        real, intent (in), dimension (3, 3, 3) :: depsB
        real, intent (in), dimension (3) :: rhat
        real, intent (in), dimension (3) :: sighat

! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcnax
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXa
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXb
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer iME
        integer index
        integer inu
        integer ix
        integer kforce
        integer nl
 
        real amt
        real argument
        real bmt
        real cost2
        real dQ_Ldx
        real dQ_Ldy
        real Q_L
        real sint
        real hx, hy
        integer nx, ny
        real xxmax, yymax, xinv
 
        real, dimension (0:ntheta - 1, ME3c_max) :: bcnalist
        real, dimension (numorb_max, numorb_max) :: bcnam
        real, dimension (0:ntheta - 1) :: dp
        real, dimension (ME3c_max) :: dphlist
        real, dimension (0:ntheta - 1, ME3c_max) :: dxbcnalist
        real, dimension (ME3c_max) :: dxhlist
        real, dimension (0:ntheta - 1, ME3c_max) :: dybcnalist
        real, dimension (ME3c_max) :: dyhlist
        real, dimension (3, numorb_max, numorb_max) :: f3naMa
        real, dimension (3, numorb_max, numorb_max) :: f3naMb
        real, dimension (ME3c_max) :: hlist
        real, dimension (0:ntheta - 1) :: p
        real, dimension (numorb_max, numorb_max) :: temp

! Procedure
! ===========================================================================
        kforce = 1
 
! Now interpolate.
! This subroutine calls the subroutine intrp1d as needed to find the value
! of the matrix elements for any given atomic configuration.
 
! For 5 angles - ntheta.

        index = icon3c(in1,in2,indna)

        if (interaction .eq. 1) then
         hx = hx_bcna(isorp,index)
         hy = hy_bcna(isorp,index)
         nx = numx3c_bcna(isorp,index)
         ny = numy3c_bcna(isorp,index)
         xxmax = x3cmax_bcna(isorp,index)
         yymax = y3cmax_bcna(isorp,index)

         if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
           write (*,*) ' What the heck is going on in Dtrescentros!!! error!!! '
           write (*,*) ' x = ', x, ' Max of data = ', xxmax
           write (*,*) ' y = ', y, ' Max of data = ', yymax
           stop
         end if

         do iME = 1, index_max3c(in1,in2)
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                         bcna_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(0,iME) = Q_L
          dxbcnalist(0,iME) = dQ_Ldx
          dybcnalist(0,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         bcna_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(1,iME) = Q_L
          dxbcnalist(1,iME) = dQ_Ldx
          dybcnalist(1,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         bcna_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(2,iME) = Q_L
          dxbcnalist(2,iME) = dQ_Ldx
          dybcnalist(2,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         bcna_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(3,iME) = Q_L
          dxbcnalist(3,iME) = dQ_Ldx
          dybcnalist(3,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, & 
     &                         bcna_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(4,iME) = Q_L
          dxbcnalist(4,iME) = dQ_Ldx
          dybcnalist(4,iME) = dQ_Ldy
         end do

        else if (interaction .eq. 2) then
         hx = hx_xc3c(isorp,index)
         hy = hy_xc3c(isorp,index)
         nx = numx3c_xc3c(isorp,index)
         ny = numy3c_xc3c(isorp,index)
         xxmax = x3cmax_xc3c(isorp,index)
         yymax = y3cmax_xc3c(isorp,index)

         if (x .gt. xxmax .or. y .gt. yymax .or. x .lt. 0 .or. y .lt. 0) then
           write (*,*) ' What the heck is going on in Dtrescentros!!! error!!! '
           write (*,*) ' x = ', x, ' Max of data = ', xxmax
           write (*,*) ' y = ', y, ' Max of data = ', yymax
           stop
         end if

         do iME = 1, index_max3c(in1,in2)
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_01(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(0,iME) = Q_L
          dxbcnalist(0,iME) = dQ_Ldx
          dybcnalist(0,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(1,iME) = Q_L
          dxbcnalist(1,iME) = dQ_Ldx
          dybcnalist(1,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(2,iME) = Q_L
          dxbcnalist(2,iME) = dQ_Ldx
          dybcnalist(2,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(3,iME) = Q_L
          dxbcnalist(3,iME) = dQ_Ldx
          dybcnalist(3,iME) = dQ_Ldy

          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  & 
     &                         xc3c_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy) 
          bcnalist(4,iME) = Q_L
          dxbcnalist(4,iME) = dQ_Ldx
          dybcnalist(4,iME) = dQ_Ldy
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
          dxbcnalist(0,iME) = dQ_Ldx
          dybcnalist(0,iME) = dQ_Ldy
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy,  &
     &                         den3_02(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(1,iME) = Q_L
          dxbcnalist(1,iME) = dQ_Ldx
          dybcnalist(1,iME) = dQ_Ldy
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_03(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(2,iME) = Q_L
          dxbcnalist(2,iME) = dQ_Ldx
          dybcnalist(2,iME) = dQ_Ldy
                                                                                                                             
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_04(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(3,iME) = Q_L
          dxbcnalist(3,iME) = dQ_Ldx
          dybcnalist(3,iME) = dQ_Ldy
                       
          call interpolate_2d (x, y, kforce, nx, ny, hx, hy, &
     &                         den3_05(1,1,iME,isorp,index), Q_L, dQ_Ldx, dQ_Ldy)
          bcnalist(4,iME) = Q_L
          dxbcnalist(4,iME) = dQ_Ldx
          dybcnalist(4,iME) = dQ_Ldy
         end do
        end if

! Now calculate the matrix elements in molecular coordinates.
! The variable cost is passed from the assembler.
        argument = 1.0d0 - cost**2
        if (argument .lt. 1.0d-5) argument = 1.0d-5
        sint = sqrt(argument) 

! Set up Legendre polys in cos(theta) and derivatives.
        if (ntheta .ne. 5) then
          write(*,*) ' ntheta must be 5, but it is ',ntheta
          stop
        end if
        p(0) = 1.0d0
        p(1) = cost
        cost2 = cost*cost
        p(2) = (3.0d0*cost2 - 1.0d0)/2.0d0
        p(3) = (5.0d0*cost2*cost - 3.0d0*cost)/2.0d0
        p(4) = (35.0d0*cost2*cost2 - 30.0d0*cost2 + 3.0d0)/8.0d0
 
        dp(0) = 0.0d0
        dp(1) = 1.0d0
        dp(2) = 3.0d0*cost
        dp(3) = (15.0d0*cost2 - 3.0d0)/2.0d0
        dp(4) = (35.0d0*cost*cost2 - 15.0d0*cost)/2.0d0
 
! ***************************************************************************
! Multiply bcnalist pieces by the appropriate Legendre polynomials and
! form hlist.
        do iME = 1, index_max3c(in1,in2)
         hlist(iME) = 0.0d0
         dphlist(iME) = 0.0d0
         dxhlist(iME) = 0.0d0
         dyhlist(iME) = 0.0d0
         do nl = 0, ntheta - 1
          hlist(iME) = hlist(iME) + p(nl)*bcnalist(nl,iME)
          dphlist(iME) = dphlist(iME) + dp(nl)*bcnalist(nl,iME)
          dxhlist(iME) = dxhlist(iME) + p(nl)*dxbcnalist(nl,iME)
          dyhlist(iME) = dyhlist(iME) + p(nl)*dybcnalist(nl,iME)
         end do
         if (mvalue(iME,in1,in2) .eq. 1) then
          dxhlist(iME) = dxhlist(iME)*sint
          dyhlist(iME) = dyhlist(iME)*sint
          if (sint .eq. 0.0d0) then
              write (*,*) ' Dividing by zero (sint = 0) in Dtrescentros.f '
              stop
          end if
          dphlist(iME) = dphlist(iME)*sint - cost*hlist(iME)/sint
          hlist(iME) = hlist(iME)*sint
         end if
        end do
 
! Now recover bcnam which is a two-dimensional array from hlist which
! is only one-dimensional.
        call recover_3c (in1, in2, hlist, bcnam)
 
! Rotate bcnam into crystal-coordinates: bcnam => bcnax
        call rotate_fb (in1, in2,  eps, bcnam, bcnax)
 
! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
        do ix = 1, 3
 
! The first piece will be the force with respect to atom 3.

         if (x .gt. 1.0d-03) then
           xinv = 1/x
         else
           xinv = x/(0.001**2) 
         end if
         amt = (sighat(ix) - cost*rhat(ix))*xinv
         do iME = 1, index_max3c(in1,in2)
          hlist(iME) = rhat(ix)*dxhlist(iME) + amt*dphlist(iME)
         end do
 
! Now recover f3naMa which is a two-dimensional array from hlist which
! is only one-dimensional.
         call recover_3c (in1, in2, hlist, temp)
         do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
           f3naMa(ix,imu,inu) = temp(imu,inu)
          end do
         end do
 
! The second piece will be the force with respect to atom 1.
         do iME = 1, index_max3c(in1,in2)
           hlist(iME) = - sighat(ix)*dyhlist(iME) - hlist(iME)/2.0d0
         end do
 
! Now recover f3naMb which is a two-dimensional array from hlist which
! is only one-dimensional.
         call recover_3c (in1, in2, hlist, temp)
         do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
           f3naMb(ix,imu,inu) = temp(imu,inu)
          end do
         end do
        end do
 
 
! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3naMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3naMa).
        call rotated (in1, in2, eps, depsA, bcnam, f3naMa, f3naXa)
 
! Force on the neutral atom with respect to atom 1 (f3naMb).
        call rotated (in1, in2, eps, depsB, bcnam, f3naMb, f3naXb)
 
! Make things force-like.
! Now determine the force f3naXc which is found from Newton's Laws.
        do inu = 1, num_orb(in2)
          do imu = 1, num_orb(in1)
           f3naXa(:,imu,inu) = - f3naXa(:,imu,inu)
           f3naXb(:,imu,inu) = - f3naXb(:,imu,inu)
           f3naXc(:,imu,inu) = - f3naXa(:,imu,inu) - f3naXb(:,imu,inu)
          end do
        end do

! Format Statements
! ===========================================================================
 
        return
      end subroutine Dtrescentros
