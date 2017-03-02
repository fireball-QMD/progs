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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! assemble_sVNL.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center sVNL (separable  
! pseudopotential) interactions.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine assemble_sVNL (iforce)
        use dimensions
        use configuration
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer in1
        integer in2
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer jatom
        integer matom
        integer mbeta
 
        real y
 
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: sVNLx
        real, dimension (3, numorb_max, numorb_max) :: spVNLx
 
! Procedure
! ===========================================================================
! Loop over the atoms in the central cell.
!$omp parallel do private (matom, r1, in1, ineigh, mbeta, jatom, r2, in2)   &
!$omp   private (r21, y, sighat, eps, deps, isorp, interaction, sVNLx)      &
!$omp   private (spVNLx, imu, inu)
        do iatom = 1, natoms                           ! <==== loop over iatom
         matom = nPP_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPn(iatom)        ! <==== loop over i's neighbors
          mbeta = nPP_b(ineigh,iatom)
          jatom = nPP_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
! Find the unit vector in sigma direction.
          if (y .lt. 1.0d-05) then
           sighat(1) = 0.0d0
           sighat(2) = 0.0d0
           sighat(3) = 1.0d0
          else
           sighat(:) = r21(:)/y
          end if
 
          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)
 
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VNL TERMS - sVNLx and spVNLx IN CRYSTAL COORDINATES.
! SET-UP sVNL and spVNL THE NEIGHBOR STORING ARRAYS.
! ****************************************************************************
! Call doscentros and get 2-center "overlap" interactions of VNL.
! That is we get <phi_i | VNL_j>. We treat it just like an overlap.
          isorp = 0
          interaction = 5
          call doscentrosPP (interaction, isorp, y, eps, deps, iforce,       &
     &                       in1, in2, sVNLx, spVNLx)
 
! Now write arrays sVNL and spVNL (derivatives).
          if (ineigh .ne. matom) then
           do inu = 1, num_orbPP(in2)
            do imu = 1, num_orb(in1)
             sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
             if (iforce .eq. 1) then
              spVNL(:,imu,inu,ineigh,iatom) = spVNLx(:,imu,inu)
             end if
            end do
           end do
          else
           do inu = 1, num_orbPP(in2)
            do imu = 1, num_orb(in1)
             sVNL(imu,inu,ineigh,iatom) = sVNLx(imu,inu)
             spVNL(:,imu,inu,ineigh,iatom) = 0.0d0
            end do
           end do
          end if

! End loop over iatom and its neighbors - jatom.
         end do
        end do
 
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine assemble_sVNL
