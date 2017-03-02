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

! ewald.f90
! Program Description
! ===========================================================================
!       This routine calculates the Ewald sum for a crystal with a given
! basis. This is specially designed for molecules with a given dipole moment.
! See Ihm, Zunger, Cohen -- Momentum Space Formalism for Total Energy of Solids
! J. Phys C v.12 (79). The ewald sum calculated here is actually gamma(ewald)
! in paper. The terms are also found in M.T. Yin and M.l. Cohen, Phys Rev. B26,
! 3259 (1982).
!
! Input parameters:
!       natoms - basis atoms
!       ratom (3,natoms) - basis atom positions (Angstroms)
!       nz(natoms) - valence charge of the atom
!       a1vec(3) - x,y,z components of lattice vector a1 (Angstroms)
!       a2vec(3) - x,y,z components of lattice vector a2
!       a3vec(3) - x,y,z components of lattice vector a3
!       g1(3) - x,y,z components of reciprocal lattice vector g1 (1/Angstrom)
!       g2(3) - x,y,z components of reciprocal lattice vector g2
!       g3(3) - x,y,z components of reciprocal lattice vector g3
!       volcel - volume of the unit cell (Angstrom**3)
!       iforce - 0 (no forces) or 1 (forces)
!
! Output:
!       ewald(natoms,natom) = ewald gamma in eV for all the basis atoms.
!       ewald = gamma1 + gamma2 + gamma3 + gamma4 - vself
!
!       gamma1 = first term of eq. 21 (Yin, Cohen paper), sum over g term.
!       gamma2 = erf term of eq. 21, sum over l term
!       gamma3 = delta(s,s') term in eq. 21.
!       gamma4 = ztot1*ztot2 term in eq. 21.
!       dewald = ewald forces
!
! ===========================================================================
! Original code from Otto F. Sankey with modifications by Alex A. Demkov
! and Jose Ortega (for charge transfer interactions).
 
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
!
! subroutine rewritten by P. Jelinek (openMP optimization)
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_ewald (natoms, nprocs, my_proc, iforce, icluster,     &
     &                        itheory, iordern, a1vec, a2vec, a3vec)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use omp_lib
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: iforce
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: my_proc
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs
 
        real, intent (in), dimension (3) :: a1vec
        real, intent (in), dimension (3) :: a2vec
        real, intent (in), dimension (3) :: a3vec
!$ volatile a1vec, a2vec, a3vec
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer iiterstart
        integer ig1
        integer ig2
        integer ig3
        integer ig1mx
        integer ig2mx
        integer ig3mx
        integer il1
        integer il2
        integer il3
        integer il1mx
        integer il2mx
        integer il3mx
        integer in1
        integer issh
        integer iteration
        integer ix
        integer jatom
        integer natomsp
        integer niters
        integer nitersp

        real argument
        real derfcdr
        real distance
        real erfc
        real factor
        real factorf
        real g1mag2, g2mag2, g3mag2
        real gdotb
        real gmax
        real gmin2
        real gsq
        real kappa
        real QQ
        real r1mag2, r2mag2, r3mag2
        real rmax
        real rmin2
        real stuff
        real volcel
 
        real, dimension (3) :: cvec
        real, dimension (3) :: eta
!        real, dimension (3, natoms) :: fewald1, fewald2
! openMP local arrays
        real, allocatable, dimension (:,:,:) :: fewald1
        real, allocatable, dimension (:,:,:) :: fewald2
        real, allocatable, dimension (:,:,:,:) :: dewaldl
        real, allocatable, dimension (:,:,:) :: ewaldl
        real, dimension (3) :: g
        real, dimension (3) :: g1, g2, g3
        real, dimension (natoms) :: Q, Q0
        real, dimension (3) :: vecl
        integer nth
        integer ith
 
! Procedure
! ===========================================================================
! get number of threads
        nth = omp_get_max_threads()
! allocate matrices
        if (iforce .eq. 1) then 
         allocate ( fewald1 (3,natoms,nth) )
         allocate ( fewald2 (3,natoms,nth) )
         allocate ( dewaldl (3,natoms,natoms,nth) )
         dewaldl = 0.0d0
        endif
        allocate ( ewaldl (natoms,natoms,nth) )
 
! Initialize ewald, dewald to zero
        ewald = 0.0d0
        ewaldl = 0.0d0
        if (iforce .eq. 1) dewald = 0.0d0
        if (iforce .eq. 1) fewald = 0.0d0
 
! Calculate delta charges (integer) into a real variable.
        do iatom = 1, natoms
         Q(iatom) = 0.0d0
         Q0(iatom) = 0.0d0
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          Q(iatom) = Q(iatom) + Qin(issh,iatom)
          Q0(iatom) = Q0(iatom) + Qneutral(issh,in1)
         end do
        end do

        if (my_proc .eq. 0 .and. icluster .eq. 0) then
         write (*,*) '  '
         write (*,*) ' Computing the ewald energy and forces.'
         write (*,*) '  '
        end if
 
! Determine the reciprical lattice vectors. Do this by Ashcroft and Mermin
! physics.
! First get a2 X a3.
        call cross (a2vec, a3vec, cvec)
 
! Next find the volume of the cell.
! NOTE: volcel actually has a sign in the above. At this point the sign is
! important since we form g vectors by dividing by a1 dot (a2 X a3).
! Oh, you say. what difference does it make if we change the sign of g.
! it makes no difference in principle.
        volcel = a1vec(1)*cvec(1) + a1vec(2)*cvec(2) + a1vec(3)*cvec(3)
        g1(:) = 2.0d0*pi*cvec(:)/volcel
 
! Next we get a3 X a1, and g2.
        call cross (a3vec, a1vec, cvec)
        g2(:) = 2.0d0*pi*cvec(:)/volcel
 
! Finally we get a1 X a2, and g3.
        call cross (a1vec, a2vec, cvec)
        g3(:) = 2.0d0*pi*cvec(:)/volcel
        volcel = abs(volcel)
 
! Initialize gmax. This determines how far we sum g1, g2, and g3. See below
! why gmax = 5.0d0 is a reasonable criterion.
        gmax = 5.0d0
 
! Initialize rmax. This determines how far we sum a1, a2, and a3. See below
! why rmax = 5.0d0 is a reasonable criterion. The parameters a1, a2, a3 are
! the direct lattice vectors.
        rmax = 5.0d0
 
! Determine the magnitude of the vectors g1, g2, g3, a1, a2, a3.
        g1mag2 = g1(1)**2 + g1(2)**2 + g1(3)**2
        g2mag2 = g2(1)**2 + g2(2)**2 + g2(3)**2
        g3mag2 = g3(1)**2 + g3(2)**2 + g3(3)**2
 
        r1mag2 = a1vec(1)**2 + a1vec(2)**2 + a1vec(3)**2
        r2mag2 = a2vec(1)**2 + a2vec(2)**2 + a2vec(3)**2
        r3mag2 = a3vec(1)**2 + a3vec(2)**2 + a3vec(3)**2
 
! ****************************************************************************
! The parameter kappa is adjustable, chosen to make the sum's converge rapidly.
! The sum over g converges as exp (-g**2/(4*kappa*kappa)), while the
! sum over l converges as exp (-r**2*kappa**2). Lets set the arguments equal
! to determine a reasonable kappa value. We set them equal for the smallest
! g value and the smallest l value.
! First find the smallest rmag.
        rmin2 = r1mag2
        if (r2mag2 .lt. rmin2) rmin2 = r2mag2
        if (r3mag2 .lt. rmin2) rmin2 = r3mag2
 
! Next find the smallest gmag.
        gmin2 = g1mag2
        if (g2mag2 .lt. gmin2) gmin2 = g2mag2
        if (g3mag2 .lt. gmin2) gmin2 = g3mag2
 
! Now set rmin2*kappa**2 = gmin2/(4*kappa**2) and solve for kappa.
        kappa = sqrt(sqrt(gmin2/(4.0d0*rmin2)))
 
! ****************************************************************************
! In gamma1 we must sum over g vectors. The decay is exp(-g**2/4*kappa*kappa).
! We require the exponent for a given direction in g-space to be gmax**2.
! For instance gmax = 5.0, corresponding to an exponent of gmax**2 = 25.0 seems
! to be a reasonable choice. This gives us g = ig1mx*g1 where
! ig1mx**2 g1**2/(4*kappa**2) = gmax**2. Solve for ig1mx, and add 1.0 for
! good measure.
        ig1mx = int(gmax * sqrt(4.0d0*kappa**2/g1mag2) + 1.0d0)
 
! Now we do the same thing for g2 and g3
        ig2mx = int(gmax * sqrt(4.0d0*kappa**2/g2mag2) + 1.0d0)
        ig3mx = int(gmax * sqrt(4.0d0*kappa**2/g3mag2) + 1.0d0)
 
        if (ig1mx .le. 1) ig1mx = 2
        if (ig2mx .le. 1) ig2mx = 2
        if (ig3mx .le. 1) ig3mx = 2
 
! In gamma2 we must sum over l vectors. The asymptotic decay is
! exp(-kappa*kappa*r**2). We require the exponent for a given direction
! in r-space to be rmax**2. For instance rmax = 5.0, corresponding to an
! exponent of rmax**2 = 25.0 seems to be a reasonable choice. This gives us
! r = il1mx*a1 where il1mx**2 r**2 * kappa**2 = rmax**2. Solve for ir1mx, and
! add 1.0 for good measure.
        il1mx = int(rmax * sqrt(1.0d0/(kappa**2*r1mag2)) + 1.0d0)
 
! Now we do the same thing for r2 and r3
        il2mx = int(rmax * sqrt(1.0d0/(kappa**2*r2mag2)) + 1.0d0)
        il3mx = int(rmax * sqrt(1.0d0/(kappa**2*r3mag2)) + 1.0d0)
 
        if (il1mx .le. 1) il1mx = 2
        if (il2mx .le. 1) il2mx = 2
        if (il3mx .le. 1) il3mx = 2

! Compute the total number of loop iterations.
        niters = (natoms*(natoms + 1)) / 2

! Determine which iterations are assigned to this processor.
        if (iordern .eq. 1) then
         nitersp = niters/nprocs
         if (my_proc .lt. mod(niters,nprocs)) then
          nitersp = nitersp + 1
          iiterstart = nitersp*my_proc + 1
         else
          iiterstart = (nitersp + 1)*mod(niters,nprocs)                      &
     &                + nitersp*(my_proc - mod(niters,nprocs)) + 1
         end if
        else
         iiterstart = 1
         nitersp = niters
        end if
 
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
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

! The real answer: now compute gamma ewald.
! ***********************************************************************
! Compute gamma1:
! ***********************************************************************
! Initialize fewald1
        if (iforce .eq. 1) fewald1 = 0.0d0
 
! Sum over g vectors.  If we are doing only a cluster, then only the gamma
! point is considered in the sum.
        if (icluster .eq. 1) then
         ig1mx = 0
         ig2mx = 0
         ig3mx = 0
        end if
        do ig1 = -ig1mx, ig1mx
         do ig2 = -ig2mx, ig2mx
          do ig3 = -ig3mx, ig3mx
 
! skip the origin
           if (.not. (ig1 .eq. 0 .and. ig2 .eq. 0 .and. ig3 .eq. 0)) then
            g(:) = ig1*g1(:) + ig2*g2(:) + ig3*g3(:)
            gsq = g(1)*g(1) + g(2)*g(2) + g(3)*g(3)
            argument = gsq/(4.0d0*kappa**2)
 
! The variable stuff contains a number of factors, including the exponential
! which is expensive to compute. That is why its outside the iatom,jatom loop.
            stuff = 4.0d0*pi*exp(-argument)/(gsq*volcel)

! Sum over s and s', the basis indices.
!$omp parallel do private(ith, factor, factorf, gdotb, QQ, iatom, jatom)
            do iteration = iiterstart, iiterstart - 1 + nitersp

! get id threads
              ith = omp_get_thread_num () + 1

              call get_atom_indices (iteration, natoms, iatom, jatom)
              factor = 1.0d0*stuff
              factorf = 2.0d0*stuff
              if (jatom .eq. iatom) factor = 0.5d0*stuff
 
! g dot b:
              gdotb = g(1)*(ratom(1,iatom) - ratom(1,jatom))                 &
     &               + g(2)*(ratom(2,iatom) - ratom(2,jatom))                &
     &               + g(3)*(ratom(3,iatom) - ratom(3,jatom))

! Calculate q(iatom)*q(jatom) - q0(iatom)*q0(jatom) = QQ
              if (itheory .eq. 1) QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
              if (itheory .eq. 2)                                            &
     &         QQ = (Q(iatom) - Q0(iatom))*(Q(jatom) - Q0(jatom))


              ewaldl(iatom,jatom,ith) =                                      &
     &          ewaldl(iatom,jatom,ith) + factor*cos(gdotb)
              ewaldl(jatom,iatom,ith) =                                      &
     &          ewaldl(jatom,iatom,ith) + factor*cos(gdotb)
                
! d/dr (cos(gdotb)) = - sin(gdotb) * d/dr (gdotb) = - sin(gdotb) * g
! The variable fewald1 is a force-like derivative => multiply by -1.0d0 
              if (iforce .eq. 1) then
               do ix = 1, 3
                fewald1(ix,iatom,ith) =                                      &
     &           fewald1(ix,iatom,ith) + QQ*factorf*sin(gdotb)*g(ix)

                fewald1(ix,jatom,ith) =                                      &
     &           fewald1(ix,jatom,ith) - QQ*factorf*sin(gdotb)*g(ix)

! The variable dewald is not a force-like derivative
                dewaldl(ix,iatom,jatom,ith) =                                &
     &           dewaldl(ix,iatom,jatom,ith) - factor*sin(gdotb)*g(ix)

                dewaldl(ix,jatom,iatom,ith) =                                &
     &           dewaldl(ix,jatom,iatom,ith) + factor*sin(gdotb)*g(ix)
               end do ! do ix
              end if ! if(iforce)
            end do ! do iteration
           end if
          end do ! do ig3
         end do ! do ig2
        end do ! do ig1
 
! ***********************************************************************
! Compute gamma2:
! ***********************************************************************
! Initialize fewald2
        if (iforce .eq. 1) fewald2 = 0.0d0

! If we are doing only a cluster, then only the central cell is considered 
! in the sum.
        if (icluster .eq. 1) then
         il1mx = 0
         il2mx = 0
         il3mx = 0
         kappa = 0.0d0
        end if

! Now carry out the sum over the cells.
        do il1 = -il1mx, il1mx
         do il2 = -il2mx, il2mx
          do il3 = -il3mx, il3mx

! Sum over atoms iatom and atoms jatom. Note that we sum over jatom .ge. iatom
! which yields an extra factor of two for iatom .ne. jatom.
!$omp parallel do private( ith, factor, factorf, iatom, jatom, QQ, vecl, eta) &
!$omp  private ( distance, argument, derfcdr)
           do iteration = iiterstart, iiterstart - 1 + nitersp

! get id threads
             ith = omp_get_thread_num () + 1

             call get_atom_indices (iteration, natoms, iatom, jatom)
             factor = 1.0d0
             factorf = 2.0d0
             if (jatom .eq. iatom) factor = 0.5d0
 
! Calculate q(iatom)*q(jatom) - q0(iatom)*q0(jatom) = QQ
             if (itheory .eq. 1) QQ = Q(iatom)*Q(jatom) - Q0(iatom)*Q0(jatom)
             if (itheory .eq. 2)                                             &
     &        QQ = (Q(iatom) - Q0(iatom))*(Q(jatom) - Q0(jatom))
 
             vecl(:) = il1*a1vec(:) + il2*a2vec(:) + il3*a3vec(:)
             eta(:) = vecl(:) + ratom(:,iatom) - ratom(:,jatom)
             distance = sqrt(eta(1)**2 + eta(2)**2 + eta(3)**2)
 
! skip the infinite self term.
             if (distance .gt. 0.0001d0) then
              argument = kappa*distance
 
              ewaldl(iatom,jatom,ith) =                                      &
     &         ewaldl(iatom,jatom,ith) + factor*erfc(argument)/distance
              ewaldl(jatom,iatom,ith) =                                      &
     &         ewaldl(jatom,iatom,ith) + factor*erfc(argument)/distance
 
! The variable fewald2 is a force-like derivative => multiply by -1.0d0
              if (iforce .eq. 1) then
               derfcdr = (2.0d0*exp(-argument**2)*kappa/sqrt(pi)             &
     &                   + erfc(argument)/distance)/distance**2
               do ix = 1, 3
                fewald2(ix,iatom,ith) =                                      &
     &           fewald2(ix,iatom,ith) + QQ*eta(ix)*factorf*derfcdr
                fewald2(ix,jatom,ith) =                                      &
     &           fewald2(ix,jatom,ith) - QQ*eta(ix)*factorf*derfcdr

! The variable dewald is not a force-like derivative
                dewaldl(ix,iatom,jatom,ith) =                                &
     &           dewaldl(ix,iatom,jatom,ith) - eta(ix)*factor*derfcdr
                dewaldl(ix,jatom,iatom,ith) =                                &
     &           dewaldl(ix,jatom,iatom,ith) + eta(ix)*factor*derfcdr
               end do ! do ix 
              end if ! if (iforce)
             end if ! if (distance)
            end do ! do iteration
           end do ! do il3
          end do ! do il2
         end do ! do il1


! ***********************************************************************
! Compute gamma3:
! ***********************************************************************
! This term should remain constant always - unless the charges as a function
! of r are changing. Also if the parameter kappa changes corresponding to
! the lattice vectors.
! There are no forces.

! sum all threads
	ewald = 0.0d0
        do ith = 1, nth
         do iatom = 1, natoms
          do jatom = 1, natoms
           ewald(iatom,jatom) = ewald(iatom,jatom) + ewaldl(iatom,jatom,ith)
          end do ! do jatom
         end do ! do iatom
        end do ! do ith

! jel-FIXME: This loop over natomsp doesn't make sence, before we summed over
! all natoms and now we use natomsp ??!
!$omp parallel do
        do iatom = iatomstart, iatomstart - 1 + natomsp
         ewald(iatom,iatom) = ewald(iatom,iatom) - 2.0d0*kappa/sqrt(pi)
        end do
! ***********************************************************************
! Compute gamma4:
! ***********************************************************************
! gamma4 is zero!
 
! ***********************************************************************
! Combine ewald pieces
! ***********************************************************************
        if (iforce .eq. 1) then
! sum over threads 
          do ith = 1,nth 
!!$omp parallel do 
           do iatom = 1, natoms
            fewald(:,iatom) =                                                &
     &        fewald(:,iatom) + fewald1(:,iatom,ith) + fewald2(:,iatom,ith)
           end do ! do iatom

!!$omp parallel do private (jatom)
           do iatom = 1, natoms
            do jatom = 1, natoms
             dewald(:,jatom,iatom) =                                         &
     &         dewald(:,jatom,iatom) + dewaldl(:,jatom,iatom,ith)
            end do ! do jatom
           end do ! do iatom

          end do ! do ith 

        end if ! if (iforce)

        if (iordern .eq. 1)                                                  &
     &   call ewald_energy_ordern_final (natoms, iforce, icluster, itheory)

! deallocate local arrays
        deallocate ( ewaldl )
        if (iforce .eq. 1) then 
         deallocate ( fewald1 )
         deallocate ( fewald2 )
         deallocate ( dewaldl )
        endif 
 
! Format Statements
! ===========================================================================
 
        return
        end



! The function get_atom_indices uses a binary search to translate from an iteration of the loop
!    do iter = 1, (natoms * (natoms + 1)) / 2
! to the loop
!    do iatom = 1, natoms
!      do jatom = iatom, natoms
! If we call get_atom_indices with the successive indices iter = 1, ..., (natoms * (natoms + 1)) / 2,
! then we get back the successive values of iatom and jatom that we would have if we executed the
! iatom/jatom loop.
! This allows us to apply OpenMP parallelism to the iter loop instead of the iatom loop, so that
! the workload is evenly divided among the processors.

        subroutine get_atom_indices (iteration, natoms, iatom, jatom)
          integer, intent (in) :: iteration, natoms
          integer, intent (out) :: iatom, jatom

          integer niters, a, b

          niters = (natoms * (natoms+1)) / 2
          if (iteration .le. natoms) then
             iatom = natoms
          else
             a = 1
             b = natoms
             iatom = (a+b) / 2
             do while (a .lt. b)
                if (iteration-1 .lt. niters - iatom*(iatom+1) / 2) then
                   a = max(a+1,iatom)
                else if (iteration .gt. niters - iatom*(iatom-1) / 2) then
                   b = min(b-1,iatom)
                else
                   exit
                end if
                iatom = (a+b) / 2
             end do
          end if
          jatom = iteration - (niters - iatom*(iatom+1) / 2)
          iatom = natoms + 1 - iatom
          jatom = iatom + jatom - 1

          return
        end subroutine get_atom_indices
