! copyright info:
!
!                             @Copyright 2009
!                FAST (Fireball Atomic Simulation Techniques)
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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

! overlap_sign.f90
! Program Description
! ===========================================================================
! Calculate non-adiabatic coupling(d_{jk}  contribution V.d_{jk}
! using Kohn-Sham states at different
! time steps, and compares with the equivalent contribution obtained
! directly using the non-adiabtic couplings calculated in
! nacouplings.f90
!
! ===========================================================================
! Code written by Enrique Abad Gonzalez
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine overlap_numeric (itime_step)

        use configuration
        use nonadiabatic
        use interactions
        use density
        use kpoints
        implicit none

! Argument Declaration and Description
! ===========================================================================


! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        integer it
        integer imu, inu, jmu,jnu, iorbital
        integer in1, in2, in3
        integer iatom, jatom
        integer isorp, interaction
        integer ix
        integer ia
        integer ik, ij
        integer ikpoint
        integer itime_step 
        real y, rcutoff_i, rcutoff_j, range
        real diff
        real delta
        real, dimension (3) :: r1, r2, r21, sighat
        real, dimension (3,3,3) :: deps
        real, dimension (3,3) :: eps
        real, dimension (numorb_max, numorb_max) :: sx
        real, dimension (3,numorb_max, numorb_max) :: spx
        real, dimension (norbitals, norbitals) :: s
        real, dimension (nele, nele) :: suma
!        real, dimension (nele, nele) :: sumb



! ===========================================================================

        if ( (itime_step .eq. 1) ) then
         allocate (sumb(nele,nele))
         allocate (ratom_old(3,natoms))
         allocate (bbnkre_old(norbitals,norbitals,nkpoints))
         allocate (blowre_old(norbitals,norbitals,nkpoints))
         ratom_old = ratom
         bbnkre_old = bbnkre
         blowre_old = blowre
        end if



        do iatom = 1, natoms
         r1(:) = ratom_old(:,iatom)
	 rcutoff_i = 0.0d0
         in1 = imass(iatom)
         do imu = 1, nssh(in1)
           if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
         end do
         do jatom = 1, natoms
          r2(:) = ratom(:,jatom)
          in2 = imass(jatom)
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
! check the cutoff radii
          rcutoff_j = 0.0d0
          in2 = imass(jatom)
          do imu = 1, nssh(in2)
              if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
          end do
          range = (rcutoff_i + rcutoff_j - 0.01d0)**2
	  range = sqrt(range)
!print *, 'llegamos al if'
          if (y .gt. range) then
!print *, 'y > range',y,range
           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
	     s(jmu,jnu) = 0.0d0
	    end do
	   end do
	  else
!print *, 'y < range',y,range
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
           else
            sighat(:) = r21(:)/y
           end if
           call epsilon (r2, sighat, eps)
           call deps2cent (r1, r2, eps, deps)
           isorp = 0
           interaction = 1
           in3 = in2
! iforce=0
           call doscentros (interaction, isorp, 0, in1, in2, in3, y, &
      &                     eps, deps, sx, spx)
           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
             s(jmu,jnu) = sx(imu,inu)
            end do
           end do
	  end if
         end do
        end do


! ===========================================================================
! Calculate overlap between Kohn-Sham states at different
! time steps
! Non-adiabatic term: dot pruduct sum
         sumb = 0.0d0
     do ikpoint = 1, nkpoints
         do ij = 1, nele
          do ik = 1, nele
           do imu = 1, norbitals
            do inu = 1, norbitals
            sumb(ik,ij) = sumb(ik,ij) +                                 & 
     &      bbnkre_old(imu,map_ks(ik),ikpoint)*bbnkre(inu,map_ks(ij),ikpoint)*s(imu,inu)
            end do
           end do
          end do
         end do
       do imu = 1, nele
        write(*,400) (sumb(imu,inu),inu=1,nele)
       end do
   end do ! end loop kpoints





! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format ('c_na',2i4,f8.4,2(f7.3))
300     format ('NAC-SUMS',2i4,2f8.4)
301     format ('S(t,tprime)',2i4,1f8.4)
400     format ('S',4f7.3)


        return
        end subroutine overlap_numeric

