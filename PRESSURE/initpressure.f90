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


! initpressure.f90
! Program Description
! ===========================================================================
!
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
        subroutine initpressure (natoms, nstepi, ipressdyn, ratom, &
     &                           vatom, a1vec, a2vec, a3vec, dt,  &
     &                           hdot, sdot, Wentz, unith)
        use dimensions
        use kpoints
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ipressdyn
        integer, intent (in) :: natoms
        integer, intent (in) :: nstepi
 
        real, intent (in) :: dt
 
        real, intent (in), dimension (3) :: a1vec, a2vec, a3vec
        real, intent (inout), dimension (3, natoms) :: ratom
        real, intent (in), dimension (3, natoms) :: vatom
 
! Output
        real, intent (inout), dimension (0:5, 3, 3) :: hdot
        real, intent (out), dimension (0:5, 3, natoms) :: sdot
        real, intent (out), dimension (3, 3) :: unith
        real, intent (out), dimension (3, 3) :: Wentz
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer idot
        integer ifactor
        integer ikpoint
        integer isum
        integer ix
        integer jx
 
        real factorial
        real xmag
 
        real, dimension (3, 3) :: hinv
        real, dimension (3, 3) :: hmet
        real, dimension (3, 3) :: hp
        real, dimension (3, 3) :: sigma
        real, dimension (3, 3) :: temp
        real, dimension (3) :: vector
 
        external factorial
 
        real volume
        real, dimension (3, 3) :: g
        real, dimension (3, 3) :: gdot
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' ************************************************ '
        write (*,*) ' Welcome to initpressure...... '
        write (*,*) '  '
 
! Initialize if and only if nstepi = 1
        if (nstepi .eq. 1) then
 
        sdot = 0.0d0
        hdot = 0.0d0
        scale_k = 0.0d0
 
! We need to set up the htensor stuff:
        hdot(0,:,1) = a1vec(:)
        hdot(0,:,2) = a2vec(:)
        hdot(0,:,3) = a3vec(:)
 
        write (*,*) '  '
        write (*,*) ' hdot: '
        write (*,100) hdot(0,1,:)
        write (*,101) hdot(0,2,:)
        write (*,102) hdot(0,3,:)
        write (*,*) '  '
 
! Now call hmetric to get g and gdot and also the volume.
! Initialize hmet and hp=d(hmet)/dt
        hmet = hdot(0,:,:)
        hp = hdot(1,:,:)
 
! Calculate scale_k's right now.
! The quantity scale_k is for a k-vector what sdot is for b-vector:
! scale_k = hmetric*kvec
        do ikpoint = 1, nkpoints
         do ix = 1, 3
          scale_k(:,ikpoint) =                                              &
     &     scale_k(:,ikpoint) + hmet(:,ix)*special_k(ix,ikpoint)
         end do
         write (*,*) '  '
         write (*,200) ikpoint, (scale_k(ix,ikpoint),ix=1,3)
        end do
 
        write (*,*) '  '
        write (*,*) ' call hmetric '
        call hmetric (hp, g, gdot, hmet, volume)
 
! Now get the inverse of hdot(0,i,k)
        call invert3x3 (hmet, hinv)
 
! We now have to put the b's into scale variables:
! s(iatom,ix) =sum (jx=1,3) hinv(ix,jx)*b(iatom,jx)
         do iatom = 1, natoms
          do ix = 1, 3
          sdot(0,:,iatom) = sdot(0,:,iatom) + hinv(:,ix)*ratom(ix,iatom)
          sdot(1,:,iatom) = sdot(1,:,iatom) + hinv(:,ix)*vatom(ix,iatom)
          end do
         end do
 
! Come here if nstepi .gt. 2!
        else
 
! ****************************************************************************
! Predict new positions,etc. before calculating the forces.
! Do the gear algorithm #5 explicitely for sdot, hdot.
         write (*,*) '  '
         write (*,*) ' GETTING GOING AFTER RESTART! '
 
! Predict new positions, etc. before calculating the forces.
! Here is the plan:
! a) calculate hinv
! b) put ratom's into s's: sdot(0,ix,i)=sum(jx)[hinv(ix,jx)*b(jx,i)
! c) put vatom's into sprime doing s'=hinv*v+hinv'*s
! d) then we predict positions and hdot and then get kpoints
         write (*,*) '  '
         write (*,*) ' Predict sdot for the first step: '
         write (*,*) '  '
         do iatom = 1, natoms
          do idot = 0, 4
           do isum = idot + 1, 5
            ifactor = isum - idot
 
! The variable factorial is a function in the MD subdirectory.
            sdot(idot,:,iatom) = sdot(idot,:,iatom)   &
     &       + (sdot(isum,:,iatom)*(dt**ifactor)/factorial(ifactor))
           end do
          end do
         end do
 
         write (*,*) '  '
         write (*,*) ' This is hdot(iorder=0,5,ix=1,3,jx=1,3): '
         write (*,*) '  '
         write (*,300) hdot(0,1,:)
         write (*,300) hdot(0,2,:)
         write (*,300) hdot(0,3,:)
         write (*,*) '  '
         write (*,300) hdot(1,1,:)
         write (*,300) hdot(1,2,:)
         write (*,300) hdot(1,3,:)
         write (*,*) '  '
         write (*,300) hdot(2,1,:)
         write (*,300) hdot(2,2,:)
         write (*,300) hdot(2,3,:)
         write (*,*) '  '
         write (*,300) hdot(3,1,:)
         write (*,300) hdot(3,2,:)
         write (*,300) hdot(3,3,:)
         write (*,*) '  '
         write (*,300) hdot(4,1,:)
         write (*,300) hdot(4,2,:)
         write (*,300) hdot(4,3,:)
         write (*,*) '  '
         write (*,300) hdot(5,1,:)
         write (*,300) hdot(5,2,:)
         write (*,300) hdot(5,3,:)
 
         write (*,*) '  '
         write (*,*) ' Predict hdot for the first step: '
         write (*,*) '  '
 
         do idot = 0, 4
          do isum = idot + 1, 5
           ifactor = isum - idot
           hdot(idot,:,1) = hdot(idot,:,1)                                   &
     &      + (hdot(isum,:,1)*(dt**ifactor)/factorial(ifactor))
           hdot(idot,:,2) = hdot(idot,:,2)                                   &
     &      + (hdot(isum,:,2)*(dt**ifactor)/factorial(ifactor))
           hdot(idot,:,3) = hdot(idot,:,3)                                   &
     &      + (hdot(isum,:,3)*(dt**ifactor)/factorial(ifactor))
          end do
         end do
 
! Call hmetric to get g and gdot.
! Initialize hmet and hp=d(hmet)/dt
         hmet = hdot(0,:,:)
         hp = hdot(1,:,:)
 
         call hmetric (hp, g, gdot, hmet, volume)
 
! Put sdot into b(ix,i)'s in order to calculate the real forces!
         do iatom = 1, natoms
          ratom(:,iatom) = 0.0d0
          do ix = 1, 3
           ratom(:,iatom) = ratom(:,iatom) + hdot(0,:,ix)*sdot(0,ix,iatom)
          end do
         end do
 
! Now we have to predict the k-points:
! How do we do that: we first take the predicted value of hdot(0,ix,iy),
! then we invert it and with the help of scale_k that we read in recalculate
! kspecial.
         hmet = hdot(0,:,:)
         call invert3x3 (hmet, hinv)

         do ikpoint = 1, nkpoints
          special_k(1,ikpoint) = 0.0d0
          do ix = 1, 3
           special_k(:,ikpoint) =                                            &
     &      special_k(:,ikpoint) + hinv(:,ix)*scale_k(ix,ikpoint)
          end do
         end do
 
! If nstepi .eq. 1 then we are done initializing.
! If nstepf .ge. 2 then we are done updating.
        end if
 
! Determine Wentz matrix if idyn = 3.
        if (ipressdyn .eq. 3) then
 
! We are doing a Wentzcovitch pressure run.
! We first construct
! sigma(alpha,beta) = d(volcel)/dh(alpha,beta)=[a2Xa3 a3Xa1 a1Xa2]
! for our REFERENCE configuration. The reference configuration is the
! a1, a2, a3 lattice vectors read in from the lattice vector file. We cannot
! call sigma to get it, since this will use the CURRENT lattice vectors, which
! change with time.
         write (*,*) '  '
         write (*,*) ' In initpressure.f ipressdyn = 3. Set up Wentz!'
 
         call cross (a2vec, a3vec, vector)
         sigma(:,1) = vector(:)
 
         call cross (a3vec, a1vec, vector)
         sigma(:,2) = vector(:)
 
         call cross (a1vec, a2vec, vector)
         sigma(:,3) = vector(:)
 
! Form sigma(dagger)*sigma.
         temp = 0.0d0
 
         do ix = 1, 3
          temp(1,:) = temp(1,:) + sigma(ix,1)*sigma(ix,:)
          temp(2,:) = temp(2,:) + sigma(ix,2)*sigma(ix,:)
          temp(3,:) = temp(3,:) + sigma(ix,3)*sigma(ix,:)
         end do
 
         write (*,*) ' The Original Wentz is:'
         do ix = 1, 3
          write (*,*) ' Wentz = ', (temp(ix,jx), jx = 1, 3)
         end do
 
! Wentzkovitch
         call invert3x3 (temp, Wentz)

         write (*,*) ' The New Wentz is:'
         do ix = 1, 3
          write (*,*) ' Wentz = ', (Wentz(ix,jx), jx = 1, 3)
         end do
 
! Projected force cell dynamics has been chosen.
        else if (ipressdyn .eq. 1) then
         xmag = sqrt(a1vec(1)**2 + a1vec(2)**2 + a1vec(3)**2)
         unith(:,1) = a1vec(:)/xmag
 
         xmag = sqrt(a2vec(1)**2 + a2vec(2)**2 + a2vec(3)**2)
         unith(:,2) = a2vec(:)/xmag
 
         xmag = sqrt(a3vec(1)**2 + a3vec(2)**2 + a3vec(3)**2)
         unith(:,3) = a3vec(:)/xmag
        end if
 
        write (*,*) ' ************************************************ '
 
! Format Statements
! ===========================================================================
100     format (2x, ' hdot(0,1,1-3) = ', 3f10.4)
101     format (2x, ' hdot(0,2,1-3) = ', 3f10.4)
102     format (2x, ' hdot(0,3,1-3) = ', 3f10.4)
200     format (2x, ' ikpoint = ', i4, ' scale_k = ', 3f12.6)
300     format (2x, 3f10.4)
 
        return
        end
