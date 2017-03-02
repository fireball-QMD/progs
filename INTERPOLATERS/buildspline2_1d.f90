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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! buildspline2_1d.f90
! Program Description
! ===========================================================================
!  1D cubic spline interpolation following Numerical Recipes in Fortran 77
!  Chapter 3.3
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine buildspline2_1d (x, y, n, spline )
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: n
        real, intent(in) :: x(n)
        real, intent(in):: y(n)
        real, intent(out):: spline(n)

! Local Variable Declaration and Description
! ===========================================================================
        integer i

        real p
        real sig
        real, dimension (n) :: u
        real yp1
        real ypn
        real h
        real qn
        real pn
        real un

! Procedure
! ===========================================================================
! Note : the points must be equally spaced!!!!
! calculate 1-th derivatives of extreme points
        h = x(1)-x(2)
        yp1 = (y(2) - y(1))/h 
        ypn = (y(n) - y(n-1))/h
        
! lower boundary conditions 'natural'
        if (yp1 .gt. .99e30) then 
           spline(1) = 0.0d0
           u(1) = 0.0d0
        else
           spline(1) = -0.5d0
           u(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif ! if (yp1)
        do i = 2,n-1
           sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
           p = sig*spline(i-1) + 2.0d0
           spline(i) = (sig - 1.0d0) / p
           u(i) = ( 6.0d0 * ((y(i+1)-y(i))/(x(i+1)-x(i))                 &
    &        - (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) -sig*u(i-1))/p
        enddo ! do i

! upper boundary conditions 'natural'
        if (yp1 .gt. .99e30) then 
           qn = 0.0d0
           un = 0.0d0
        else
           qn = 0.5d0
           un = (3.0d0/(x(n)-x(n-1)))*(ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
        endif ! if (yp1)
        spline(n) = (un-qn*u(n-1))/(qn*spline(n-1)+1.0d0)
        
! backsubstitution loop of the tridiagonal algorithm
        do i = n-1,1,-1
           spline(i) = spline(i)*spline(i+1) + u(i)
        enddo ! do i

! Format Statements
! ===========================================================================

        return
      end subroutine buildspline2_1d
