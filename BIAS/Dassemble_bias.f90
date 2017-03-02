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

! Dassemble_bias.f90
! Program Description
! ===========================================================================
! 	This routine assembles the two-center forces
! related to the applied bias voltage.
! Scheme
!
!       eV/2       zb0
!  -----------------\
!                    \ eV*(1/2-(z-zb0)/(zb1-zb0))
!                     \
!                      \--------------------------
!  0                    zb1          -eV/2      Lz
!  ---> z-axis
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
! ===========================================================================
!
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_bias (nprocs, impi)
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use bias
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: impi
        integer, intent (in) :: nprocs

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer ierror
        integer imu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer ix
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp

        real sumS
        real y
        real z1
	real z2
	real za
        real Vb

        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize forces to zero.

! Determine which atoms are assigned to this processor.
        if (impi .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

        fbias = 0.0d0

! Loop over the atoms in the central cell.
!!$omp parallel do private (matom, r1, in1, ineigh, mbeta, jatom, r2, in2)    &
!!$omp&  private (r21, y, sighat, eps, deps, ix, sumT, sumS, inu, imu, muxc)  &
!!$omp&  private (isorp, kforce, interaction, in3, bcnax, bcnapx)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)     ! <==== loop 2 over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

! ****************************************************************************
!
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

! evaluate the bias potential
          z1 = ratom (3,iatom)
          z2 = ratom (3,jatom)
          za = (z1 + z2)/2.0d0
          if (za .gt. zb1) then
           Vb = -0.5d0*Vbias
          else if (za .gt. zb0) then
           Vb = (0.5d0 - (za - zb0)/(zb1-zb0))*Vbias
          else
           Vb = 0.5d0*Vbias
          endif

! ****************************************************************************
!
!  ASSEMBLE FORCES
! ****************************************************************************
! The derivatives are tpx and spx, where p means derivative and x means crytal
! coordinates. The derivatives are:  d/d(b1) <phi(r-b1) ! O ! phi(r-b2)>
! where O = 1 (for overlap), or T.  The derivative is a vector in crystal
! coordinates and is stored in tp(ix,imu,inu,iatom,ineigh). The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4) or spx(3,4,4).  The
! derivatives are computed only if iforce = 1 !
!
! Overlap and Kinetic.
! No need to call doscentros to get spx or tpx, because it was calculated in
! assemble_2c.f.
! Now do kinetic and repulsive overlap forces.
          do ix = 1, 3
           sumS = 0.0d0
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sumS = sumS + sp_mat(ix,imu,inu,ineigh,iatom)
            end do
           end do

! Now add sum to appropriate force term. see notes "the total band structure
! force", ofs 9/14/88. Also see notes "total force due to fr offsites"
!                                 and "total force due to fs offsites" 11/10/88
! The (-1.d0) makes it "force-like".
!!$omp critical
! Direct terms.
           fbias(ix,iatom) = fbias(ix,iatom) + (-1.0d0)*sumS*Vb

! Cross terms.
           fbias (ix,jatom) = fbias (ix,jatom) - (-1.0d0)*sumS*Vb
!!$omp end critical
          end do ! do ix

! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================
        return
      end subroutine Dassemble_bias
