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

! writeout_cdcoeffs.f90
! Program Description
! ===========================================================================
!       This routine writes out the charge density coefficients.
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
        subroutine writeout_cd (icluster, iwrtcdcoefs, itime_step)
        use configuration
        use dimensions
        use density
        use interactions
        use kpoints
        use neighbor_map
        use charges
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: itime_step
        integer, intent (in) :: iwrtcdcoefs


! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: norbmax = 20

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ikpoint
        integer imu, inu
        integer in1
        integer inum
        integer iremainder
        integer issh
        integer lm

        integer iband
        integer jatom
        integer mbeta
        integer ineigh
        integer in2
        integer nnu,mmu
        integer mindex
        integer, dimension (:), allocatable :: ainx
        integer, dimension (:), allocatable :: oinx

        real gutr
        real pcharge
        real, dimension (:), allocatable :: pocharge
        real, dimension (3) :: vec
        real dot
        real mval

	complex ai
        complex phase
        complex step1, step2

        real zero
        character (len=3) :: zorb

! Procedure
! ===========================================================================
! Initialize
        zero = 0.0d0
	ai = cmplx(0.0d0,1.0d0)

! ****************************************************************************
!
!        W R I T E - O U T    C H A R G E    D E N S I T I E S
! ****************************************************************************
! This section will write to a file the coefficients for a charge density
! calculation. If iwrtcdcoefs = 0 skip over this section. This is rarely wanted.
        write (*,*) '  '
        write (*,*) ' Writing out the charge densities. '

        if (iwrtcdcoefs .ne. 2) then
         open (unit = 39, file = 'cdcoeffs.dat', status = 'unknown',         &
     &         position = 'append')
        else
         open (unit = 39, file = 'cdcoeffs.dat', status = 'unknown',         &
     &         form = 'unformatted', position = 'append')
        end if
        if (itime_step .eq. 1) then
         write (39,*) norbitals
         if (iwrtcdcoefs .eq. 1) then
          write (39,100) xl(:,1)
          write (39,100) xl(:,3)
          write (39,100) xl(:,5)
         else if (iwrtcdcoefs .eq. 2) then
          write (39,*) xl(:,1)
          write (39,*) xl(:,3)
          write (39,*) xl(:,5)
         else if (iwrtcdcoefs .eq. 3) then
          write (39,100) xl(:,1)
          write (39,100) xl(:,3)
          write (39,100) xl(:,5)
         end if
         write (39,*) nkpoints

! Loop over number of kpoints
         do ikpoint = 1, nkpoints
          if (iwrtcdcoefs .eq. 1) then
           write (39,101) special_k(:,ikpoint), weight_k(ikpoint)
          else if (iwrtcdcoefs .eq. 2) then
           write (39,*) special_k(:,ikpoint), weight_k(ikpoint)
          else if (iwrtcdcoefs .eq. 3) then
           write (39,101) special_k(:,ikpoint), weight_k(ikpoint)
          end if
         end do
        end if

! Loop over number of kpoints
        do ikpoint = 1, nkpoints
         do imu = 1, norbitals_new
          write (39,*) imu, eigen_k(imu,ikpoint)

! Write out the points
          inum = int(norbitals_new/4) ! We want int/int division
          iremainder = norbitals_new - (inum*4)
          do inu = 1, norbitals_new - iremainder, 4
           if (icluster .ne. 1) then
            if (iwrtcdcoefs .eq. 1) then
             write (39,102)                                                  &
     &       bbnkre(inu  ,imu,ikpoint), bbnkim(inu  ,imu,ikpoint),           &
     &       bbnkre(inu+1,imu,ikpoint), bbnkim(inu+1,imu,ikpoint),           &
     &       bbnkre(inu+2,imu,ikpoint), bbnkim(inu+2,imu,ikpoint),           &
     &       bbnkre(inu+3,imu,ikpoint), bbnkim(inu+3,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                    &
     &       bbnkre(inu  ,imu,ikpoint), bbnkim(inu  ,imu,ikpoint),           &
     &       bbnkre(inu+1,imu,ikpoint), bbnkim(inu+1,imu,ikpoint),           &
     &       bbnkre(inu+2,imu,ikpoint), bbnkim(inu+2,imu,ikpoint),           &
     &       bbnkre(inu+3,imu,ikpoint), bbnkim(inu+3,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 3) then
             write (39,102)                                                  &
     &       blowre(inu  ,imu,ikpoint), blowim(inu  ,imu,ikpoint),           &
     &       blowre(inu+1,imu,ikpoint), blowim(inu+1,imu,ikpoint),           &
     &       blowre(inu+2,imu,ikpoint), blowim(inu+2,imu,ikpoint),           &
     &       blowre(inu+3,imu,ikpoint), blowim(inu+3,imu,ikpoint)
            end if
           else
            if (iwrtcdcoefs .eq. 1) then
             write (39,102)                                                  &
     &       bbnkre(inu  ,imu,ikpoint), zero,                                &
     &       bbnkre(inu+1,imu,ikpoint), zero,                                &
     &       bbnkre(inu+2,imu,ikpoint), zero,                                &
     &       bbnkre(inu+3,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                    &
     &       bbnkre(inu  ,imu,ikpoint), zero,                                &
     &       bbnkre(inu+1,imu,ikpoint), zero,                                &
     &       bbnkre(inu+2,imu,ikpoint), zero,                                &
     &       bbnkre(inu+3,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 3) then
             write (39,102)                                                  &
     &       blowre(inu  ,imu,ikpoint), zero,                                &
     &       blowre(inu+1,imu,ikpoint), zero,                                &
     &       blowre(inu+2,imu,ikpoint), zero,                                &
     &       blowre(inu+3,imu,ikpoint), zero
            end if
           end if
          end do

          if (iremainder .eq. 1) then
           if (icluster .ne. 1) then
            if (iwrtcdcoefs .eq. 1) then
             write (39,103)                                                  &
     &        bbnkre(norbitals_new,imu,ikpoint),                             &
     &        bbnkim(norbitals_new,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                    &
     &        bbnkre(norbitals_new,imu,ikpoint),                             &
     &        bbnkim(norbitals_new,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 3) then
             write (39,103)                                                  &
     &        blowre(norbitals_new,imu,ikpoint),                             &
     &        blowim(norbitals_new,imu,ikpoint)
            end if
           else
            if (iwrtcdcoefs .eq. 1) then
             write (39,103) bbnkre(norbitals_new,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*) bbnkre(norbitals_new,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 3) then
             write (39,103) blowre(norbitals_new,imu,ikpoint), zero
            end if
           end if
          else if (iremainder .eq. 2) then
           if (icluster .ne. 1) then
            if (iwrtcdcoefs .eq. 1) then
             write (39,104)                                                  &
     &        bbnkre(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkre(norbitals_new  ,imu,ikpoint),                           &
     &        bbnkim(norbitals_new  ,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                   &
     &        bbnkre(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkre(norbitals_new  ,imu,ikpoint),                           &
     &        bbnkim(norbitals_new  ,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 3) then
             write (39,104)                                                  &
     &        blowre(norbitals_new-1,imu,ikpoint),                           &
     &        blowim(norbitals_new-1,imu,ikpoint),                           &
     &        blowre(norbitals_new  ,imu,ikpoint),                           &
     &        blowim(norbitals_new  ,imu,ikpoint)
            end if
           else
            if (iwrtcdcoefs .eq. 1) then
             write (39,104)                                                  &
     &        bbnkre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new  ,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                   &
     &        bbnkre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new  ,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 3) then
             write (39,104)                                                  &
     &        blowre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        blowre(norbitals_new  ,imu,ikpoint), zero
            end if
           end if
          else if (iremainder .eq. 3) then
           if (icluster .ne. 1) then
            if (iwrtcdcoefs .eq. 1) then
             write (39,105)                                                  &
     &        bbnkre(norbitals_new-2,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-2,imu,ikpoint),                           &
     &        bbnkre(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkre(norbitals_new  ,imu,ikpoint),                           &
     &        bbnkim(norbitals_new  ,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                   &
     &        bbnkre(norbitals_new-2,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-2,imu,ikpoint),                           &
     &        bbnkre(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkim(norbitals_new-1,imu,ikpoint),                           &
     &        bbnkre(norbitals_new  ,imu,ikpoint),                           &
     &        bbnkim(norbitals_new  ,imu,ikpoint)
            else if (iwrtcdcoefs .eq. 3) then
             write (39,105)                                                  &
     &        blowre(norbitals_new-2,imu,ikpoint),                           &
     &        blowim(norbitals_new-2,imu,ikpoint),                           &
     &        blowre(norbitals_new-1,imu,ikpoint),                           &
     &        blowim(norbitals_new-1,imu,ikpoint),                           &
     &        blowre(norbitals_new  ,imu,ikpoint),                           &
     &        blowim(norbitals_new  ,imu,ikpoint)
            end if
           else
            if (iwrtcdcoefs .eq. 1) then
             write (39,105)                                                  &
     &        bbnkre(norbitals_new-2,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new  ,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 2) then
             write (39,*)                                                   &
     &        bbnkre(norbitals_new-2,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        bbnkre(norbitals_new  ,imu,ikpoint), zero
            else if (iwrtcdcoefs .eq. 3) then
             write (39,105)                                                  &
     &        blowre(norbitals_new-2,imu,ikpoint), zero,                     &
     &        blowre(norbitals_new-1,imu,ikpoint), zero,                     &
     &        blowre(norbitals_new  ,imu,ikpoint), zero
            end if
           end if
          end if
         end do
        end do


! If you are monkeying with the occupation it is probable that you would like
! to know about charge localization on the snuffed sites, how much
! eigenvector there is in orbital mu, etc.
        write (*,*) '  '
! allocate aux array
        allocate ( pocharge(norbitals) )
        allocate ( ainx(norbitals) )
        allocate ( oinx(norbitals) )
! open files
        open (unit = 40, file = 'Q2atom.dat', status = 'unknown')
        open (unit = 41, file = 'Q2orbital.dat', status = 'unknown')

! set up aux arrays
        do iatom = 1, natoms
         in1 = imass(iatom)
         do imu = 1,num_orb(in1)
          mmu = imu + degelec(iatom)
          ainx(mmu) = iatom
          oinx(mmu) = imu
         end do
        end do

! Loop over bands
        do iband = 1, norbitals

         pocharge = 0.0d0
! Loop over the atoms and the neighbors
         do iatom = 1, natoms

          in1 = imass(iatom)
          pcharge = 0.0d0

          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

! Loop over the special k points.
           do ikpoint = 1, nkpoints
            dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) &
    &                                         + special_k(3,ikpoint)*vec(3)
            phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)
! Loop over all bands
            if (icluster .ne. 1) then
             do imu = 1, num_orb(in1)
              mmu = imu + degelec(iatom)
              step1 = phase*(bbnkre(mmu,iband,ikpoint)                      &
    &                        - ai*bbnkim(mmu,iband,ikpoint))
              do inu = 1, num_orb(in2)
               nnu = inu + degelec(jatom)
               step2 = step1*(bbnkre(nnu,iband,ikpoint)                     &
    &                         + ai*bbnkim(nnu,iband,ikpoint))
               gutr = real(step2)
! projected charges per atom
               pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom)
! projected charges per orbital
               pocharge(mmu) = pocharge(mmu) + gutr*s_mat(imu,inu,ineigh,iatom)
! Finished loop over bands
              end do
             end do
            else
             do imu = 1, num_orb(in1)
              mmu = imu + degelec(iatom)
              step1 = phase*bbnkre(mmu,iband,ikpoint)
              do inu = 1, num_orb(in2)
               nnu = inu + degelec(jatom)
               step2 = step1*bbnkre(nnu,iband,ikpoint)
               gutr = real(step2)
! projected charges per atom
               pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom)
! projected charges per orbital
               pocharge(mmu) = pocharge(mmu) + gutr*s_mat(imu,inu,ineigh,iatom)
! Finished loop over bands
              end do ! do inu
             end do ! do imu
            end if ! if (icluster)

! Finished loop over kpoints
           end do ! do ikpoint

! Finish remaining loops
          end do ! do ineigh
          write (40,200) iband, iatom, pcharge
         end do ! do iatom

! Writeout 3 orbitals which belongs to the band
         write (41,*) ''
         write (41,*) '=================================================='
         write (41,400) iband
         write (41,*) '=================================================='

! loop over x-maxim
         do in2 = 1,norbmax
! find maximum value
          mval = -1.0d0
          do inu = 1, norbitals
           if (mval .lt. pocharge(inu)) then
            mval = pocharge(inu)
            mindex = inu
           endif
          end do ! do inu

! find kind of orbital
          in1 = imass(ainx(mindex))
          inu = 1
          imu = oinx(mindex)
          iatom = ainx(mindex)
          do issh = 1, nssh(in1)
           lm = lssh(issh,in1)
           do inum = -lm,lm
            if (imu .eq. inu ) then
             if ( lm .eq. 0 ) then
              zorb = 's     '
             else if ( lm .eq. 1 ) then
              if (inum .eq. -1) zorb = 'px    '
              if (inum .eq. 0) zorb = 'pz    '
              if (inum .eq. 1) zorb = 'py    '
             else if ( lm .eq. 2) then
              if (inum .eq. -2) zorb = 'dxy   '
              if (inum .eq. -1) zorb = 'dyz   '
              if (inum .eq. 0) zorb = 'dz2   '
              if (inum .eq. 1) zorb = 'dxz   '
              if (inum .eq. 2) zorb = 'dx2-y2'
             endif
            end if
            inu = inu + 1
           end do !
          end do ! do issh
! writeout
          write (41,401) iatom, mindex, imu, zorb, pocharge(mindex)
          pocharge(mindex) = -10.0d0
         end do ! do imu

        end do ! do iband

! deallocate aux arrays
        deallocate (pocharge)
        deallocate (oinx)
        deallocate (ainx)

! close the open files
        close (unit = 39)
        close (unit = 40)
        close (unit = 41)

! Format Statements
! ===========================================================================
100     format (3f11.6)
101     format (4f11.5)
102     format (4(1x, '(', f14.11, ',', f14.11, ')'))
103     format (1(1x, '(', f14.11, ',', f14.11, ')'))
104     format (2(1x, '(', f14.11, ',', f14.11, ')'))
105     format (3(1x, '(', f14.11, ',', f14.11, ')'))
200     format (2x,' Band : ',i4,' Atom :',i4,' pcharge :',f8.4)
400     format (2x,' Band = ',i4)
401     format (4x,' Atom = ',i4,' orb_glob =',i4,' orb_loc =',i4,'',a6,' charge =',f8.4)


        return
        end
