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

! twisterd.f90
! Program Description
! ===========================================================================
!       This routine prepares the derivatives of the D matrices for a given
! geometry of the matrix element.  This routine needs to know the variable
! epsilon.  This epsilon is either eps2 or eps3, depending on what matrix
! element we are dealing with (e.g. T is a 2C, while Vna is a 3C case).
!
! Input:
! The variable eps is a 3x3 output of the subroutine epsilon.
! The variable deps is a 3x3X3 output of either deps2cent, or deps3cent.
 
! Output:
! The variable ddmat is a 3x5x5 matrix  derivative rotating d-orbitals.
! The variable dpmat is a 3x3x3 matrix  derivative rotating p-orbitals.
!
! Here is the famous Ortega convention:
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! ===========================================================================
! Original code from Alex Demkov

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
        subroutine twisterd (eps, deps, ddmat, dpmat)
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in) :: eps (3, 3)
        real, intent(in) :: deps (3, 3, 3)
 
! Output
        real, intent(out) :: ddmat (3, 5, 5)
        real, intent(out) :: dpmat (3, 3, 3)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ix, jx, kx
 
        real aterm12, aterm32, aterm33, aterm13, aterm11
 
        real amat_term
 
! Procedure
! ===========================================================================
! Nothing for S orbitals.
! Set the dP/drk matrices: eps: x,y,z; pmat: y,z,x
        do ix = 1, 3
         dpmat(ix,1,1) = deps(ix,2,2)
         dpmat(ix,1,2) = deps(ix,2,3)
         dpmat(ix,1,3) = deps(ix,2,1)
 
         dpmat(ix,2,1) = deps(ix,3,2)
         dpmat(ix,2,2) = deps(ix,3,3)
         dpmat(ix,2,3) = deps(ix,3,1)
 
         dpmat(ix,3,1) = deps(ix,1,2)
         dpmat(ix,3,2) = deps(ix,1,3)
         dpmat(ix,3,3) = deps(ix,1,1)
        end do

        if (.not. haveDorbitals) return

! ***************************************************************************
! Set the dD/dr matrices according to the general formula
! (see p.3 notes 9/29/98):
!
! dD(mu|m)/dr_k ~ dLAMBDA(mu)_[a,b]/d_r = SUM_(i,j)
!
! a(mu)_[i,j]((deps(i,a)/dr_k)*eps(j,b) + eps(i,a)*(deps(j,b)/dr_k) )
!
! and m defines [a,b], e.g. m = -2 requires [1,2] (see p. 3 notes 9/29/98),
! Note the formula is "mu-independent" for a given m!
! ***************************************************************************
        do imu = 1, 5
         do kx = 1, 3
          aterm12 = 0.0d0
          aterm32 = 0.0d0 
          aterm33 = 0.0d0 
          aterm13 = 0.0d0 
          aterm11 = 0.0d0 
          do ix = 1, 3
           do jx = 1, 3
            if (amat(jx,ix,imu) .ne. 0.0d0) then
             amat_term = amat(jx,ix,imu)
             aterm12 = aterm12 + amat_term*(deps(kx,ix,1)*eps(jx,2) + eps(ix,1)*deps(kx,jx,2))
             aterm32 = aterm32 + amat_term*(deps(kx,ix,3)*eps(jx,2) + eps(ix,3)*deps(kx,jx,2))
             aterm33 = aterm33 + amat_term*(deps(kx,ix,3)*eps(jx,3) + eps(ix,3)*deps(kx,jx,3))
             aterm13 = aterm13 + amat_term*(deps(kx,ix,1)*eps(jx,3) + eps(ix,1)*deps(kx,jx,3))
             aterm11 = aterm11 + amat_term*(deps(kx,ix,1)*eps(jx,1) + eps(ix,1)*deps(kx,jx,1))
            end if
           end do
          end do
! Creating derivatives of D(mu|-2) ==> D(mu|1), LAMBDA_[12]
          ddmat(kx,imu,1) = 2.0d0*aterm12
! Creating derivatives of D(mu|-1) ==> D(mu|2)  , LAMBDA_[32]
          ddmat(kx,imu,2) = 2.0d0*aterm32
! Creating derivatives of D(mu|0)  ==> D(mu|3), LAMBDA_[33]
          ddmat(kx,imu,3) = sqrt(3.0d0)*aterm33
! Creating derivatives of D(mu|1)  ==> D(mu|4), LAMBDA_[13]
          ddmat(kx,imu,4) = 2.0d0*aterm13
! Creating derivatives of D(mu|2)  ==> D(mu|5), 2*LAMBDA_[11]+LAMBDA_[33]
          ddmat(kx,imu,5) = 2.0d0*aterm11 + aterm33
         end do
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
