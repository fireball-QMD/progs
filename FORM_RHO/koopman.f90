! copyright info:
!
!                             @Copyright 2009
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

! koopman.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
! ===========================================================================
! Code written by:
! Enrique Abad Gonzalez & Yannick J Dappe
! Dpto. de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! Phone: +34-91-497-86-48
! ===========================================================================
!
! Program Declaration
! ===========================================================================

        subroutine koopman (natoms,nkpoints,ratom,ioccupy_k,foccupy)


        use neighbor_map
        use interactions
        use hartree_fock
        use density
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nkpoints
        real, intent (in), dimension (3,natoms) :: ratom
        integer, intent (in), dimension (norbitals, nkpoints) :: ioccupy_k
        real, intent(in), dimension (norbitals, nkpoints) :: foccupy
! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom
        integer ineigh
        integer in1, in2
        integer imu, inu
        integer jmu, jnu
        integer ieigen
        real    r21mod
        real,dimension(norbitals) :: deltae
        real,dimension(3) :: r21
        real,dimension(norbitals,norbitals) :: delta_n

! Procedure
! ===========================================================================

      open (unit = 142, file = 'koopman_shift', status = 'unknown')

! Test: we need to being working in a cluster calculation:

      if ( nkpoints .ne. 1 ) then
        write (*,*) 'Koopman gap correction need to be done with'
        write (*,*) 'gamma k-point ONLY. You are using many k-points'
        write (*,*) 'so this calculation is useless. Please, rerun the'
        write (*,*) 'code using gamma point only'
        STOP
      end if


! First of all we need to write the lowdin coefficients in order to use
! later the scissor operator, this way seems a little bit "ugly" but
! works preety well

!      do ikpoint = 1
!        open (unit=222+ikpoint, status='unknown')
!          do imu =1, 240    !norbitals
!            do inu =1, 240     !norbitals
!              write(222+ikpoint,*) blowre(inu,imu,ikpoint), blowim(inu,imu,ikpoint)
!            end do
!          end do
!        close (unit=222+ikpoint)
!      end do

! Now it's time to calculate the koopman correction

      deltae = 0.0

      do ieigen = 1, norbitals ! or should be norbitals_new??

      ! First of all we calculate delta_n for this eigenstate
        do jmu = 1, norbitals
          delta_n(jmu,ieigen) = blowre(jmu,ieigen,1)**2
        end do
        ! Now we calculate the deltae for this eigenstate

        do iatom = 1, natoms
          in1 = imass(iatom)
          do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            ! diagonal part
            deltae(ieigen) = deltae(ieigen) + (1-(2*ioccupy_k(ieigen,1)))* &
     &   0.5* Jialpha(imu,iatom)*(delta_n(jmu,ieigen)**2)
            ! (1-2*iocc_k(ieigen)) gives -1 for filled states and + 1
            ! for empty states (that's that we want)
!        print *, 'parte Jialfa'
!        print *, ieigen, deltae(ieigen)
            ! lo de la Uisigma igual habia que cambiarlo mas adelante
            ! y lo de efe tambien hay que hacerlo bien
!            deltae(ieigen) = deltae(ieigen) + efe(imu,iatom)*       &
            deltae(ieigen) = deltae(ieigen) + (1-(2*ioccupy_k(ieigen,1)))*0.5*0.5* &
     & (Uisigma(imu,imu,iatom)-Jialpha(imu,iatom))*(delta_n(jmu,ieigen)**2)
!        print *, 'parte f(U-J)'
!        print *, ieigen, deltae(ieigen)
            ! non-diagonal part
            do ineigh = 1, neighn(iatom)
              jatom = neigh_j(ineigh,iatom)
              in2 = imass(jatom)
              if (iatom.eq.jatom) then
                do inu = 1, num_orb(in2)
                  jnu = inu + degelec(jatom)
                  deltae(ieigen) = deltae(ieigen) + (1-(2*ioccupy_k(ieigen,1)))* &
     & 0.5*Uisigma(imu,inu,in1)*delta_n(jmu,ieigen)*delta_n(jnu,ieigen)
!        print *, 'parte no diagonal (U)'
!        print *, ieigen, deltae(ieigen)
                end do
              else
                do inu = 1, num_orb(in2)
                  jnu = inu + degelec(jatom)
                  r21(:) = ratom(:,jatom) - ratom(:,iatom)
                  r21mod = sqrt(r21(1)**2+r21(2)**2+r21(3)**2)
        ! VAMOS A PROBAR PRIMEROS VECIONS
!        if (r21mod .lt. 1.60 ) then
                  deltae(ieigen) = deltae(ieigen) + (1-(2*ioccupy_k(ieigen,1)))* &
     & 0.5*(Jijsigma(imu,in1,inu,in2)/r21mod)*delta_n(jmu,ieigen)*delta_n(jnu,ieigen)
!        print *, 'parte no diagonal (J)'
!        print *, ieigen, deltae(ieigen),iatom,jatom
!        end if ! fin TEST
                end do ! inu
              end if
            end do ! jatom
          end do ! imu
        end do !iatom

        write (142,*) deltae(ieigen)
!	print *, deltae(ieigen)


      end do  ! end loop on ieigen

      close (unit = 142)

! Format Statements
! ===========================================================================
101     format (f10.6)

        return
      end subroutine koopman

