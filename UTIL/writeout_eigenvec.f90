

        subroutine writeout_eigenvec ()
      !  use configuration
      !  use density
      !  use dimensions
      !  use interactions
      !  use kpoints
      !  use neighbor_map

       ! use charges
        use configuration
       ! use constants
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        implicit none

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom,ineigh,jatom
        integer iband
        integer imu,inu
        integer nnu,mmu
        integer in1,in2
        integer mbeta
        
        integer issh, l,lmu

        real dot
        real dens,densS,densP,densD,densPx,densPy,densPz,dd
        real densDz2, densDxz, densDyz, densDxy, densDx2y2

        real alpha
        complex phase,step1,step2,ai
        complex coef

        real,    dimension (3) :: k,R 
        integer :: ikpoint                   ! index of kpoint in a,b,c direction 

        complex, dimension (:), allocatable :: eigenv

! Procedure
! ===========================================================================

        ai = cmplx(0.0d0,1.0d0)
        allocate(eigenv(norbitals_new))

        write (*,*) ' ****************************************************** '
        write (*,*) '         Write Eigenvectors  (by Prokop Hapala)         '
        write (*,*) ' ****************************************************** '

      open (unit = 72, file = 'eigenvec_k.dat', status = 'unknown')
      open (unit = 73, file = 'eigenrho.dat', status = 'unknown')

      open (unit = 74, file = 'eigenrho_s.dat', status = 'unknown')
      open (unit = 75, file = 'eigenrho_p.dat', status = 'unknown')
      open (unit = 76, file = 'eigenrho_d.dat', status = 'unknown')

      open (unit = 77, file = 'eigenrho_px.dat', status = 'unknown')
      open (unit = 78, file = 'eigenrho_py.dat', status = 'unknown')
      open (unit = 79, file = 'eigenrho_pz.dat', status = 'unknown')

      open (unit = 80, file = 'eigenrho_dz2.dat', status = 'unknown')
      open (unit = 81, file = 'eigenrho_dxz.dat', status = 'unknown')
      open (unit = 82, file = 'eigenrho_dyz.dat', status = 'unknown')
      open (unit = 83, file = 'eigenrho_dxy.dat', status = 'unknown')
      open (unit = 84, file = 'eigenrho_dx2y2.dat', status = 'unknown')




! ==========================================================
!    eigenvec_k.dat        Write out band eigenvectors Uk
! ==========================================================


write (72,*) natoms, norbitals_new, nkpoints 
do iband = 1, norbitals_new
     write (72,*) " band ",iband
 do  ikpoint = 1, nkpoints
     k = special_k(:,ikpoint)
     eigenv(:) = (bbnkre(:,iband,ikpoint) + ai*bbnkim(:,iband,ikpoint)) 
     write (72,'(i10,f10.5,f10.5,f10.5,f10.5)') ikpoint,k(1),k(2),k(3),eigen_k(iband,ikpoint)
     do iatom = 1, natoms
          in1 = imass(iatom)
          R = ratom(:,iatom)       - ratom(:,natoms)
          dot = k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
          phase = cmplx(cos(dot),sin(-dot))   ! with -dot inverse fourier
          write (72,'(i5,f15.10,f15.10,f15.10,f15.10)',advance = 'no') iatom  ! , ratom(3,iatom), dot, real(phase), aimag(phase)
         do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               write (72,'(f15.10,f15.10)',advance = 'no') real(eigenv(mmu)), aimag(eigenv(mmu))
               eigenv(mmu) = phase*eigenv(mmu)
               alpha = ATAN2( real(eigenv(mmu)), aimag(eigenv(mmu)))
               if (alpha .lt. 0 ) alpha = alpha + 3.14159265 
               write (72,'(f15.10,f15.10,f15.10)',advance = 'no') real(eigenv(mmu)), aimag(eigenv(mmu)), alpha
         end do 
         do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
                ! coef = phase*(bbnkre(mmu,iband,ikpoint) + ai*bbnkim(mmu,iband,ikpoint) )
                ! coef = (bbnkre(mmu,iband,ikpoint) + ai*bbnkim(mmu,iband,ikpoint) )
                ! coef = bbnkre(mmu,iband,ikpoint)
                ! write (72,'(f15.10,f15.10,f15.10,f15.10)',advance = 'no') bbnkre(mmu,iband,ikpoint), bbnkim(mmu,iband,ikpoint), real(coef), aimag(coef)
          end do ! imu
          write (72,*)
      end do  ! iatom
  end do  ! ikpoint
end do  ! iband          
 close(72)

! =============================================================
!     eigenrho.dat    Write out Mulliken like PDOS of atom 
! =============================================================
write (73,*) natoms, norbitals_new, nkpoints
do iband = 1, norbitals_new
 do  ikpoint = 1, nkpoints
      k = special_k(:,ikpoint)
      write (73,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      do iatom = 1, natoms
          in1 = imass(iatom)
          dens = 0.0
          do imu = 1, num_orb(in1)
            mmu = imu + degelec(iatom)
           do ineigh = 1, neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            R = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
            dot = k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
            phase = cmplx(cos(dot),sin(dot))
            step1 = phase*(bbnkre(mmu,iband,ikpoint)- ai*bbnkim(mmu,iband,ikpoint))
              do inu = 1, num_orb(in2)
               nnu = inu + degelec(jatom)
               step2 = step1*(bbnkre(nnu,iband,ikpoint)+ai*bbnkim(nnu,iband,ikpoint))
               dens = dens + real(step2)*s_mat(imu,inu,ineigh,iatom)
              end do ! inu
             end do ! ineigh

          !   step1 = (bbnkre(mmu,iband,ikpoint)- ai*bbnkim(mmu,iband,ikpoint))
          !   dens = dens + step1*CONJG(step1)
            end do ! imu
          ! write (73,'(f15.10)',advance = 'no') dens*weight_k(ikpoint)
          write (73,'(f15.10)',advance = 'no') dens
      end do  ! iatom
      write (73,*)    
  end do  ! ikpoint
end do  ! iband   
  close(73)



! =============================================================
!     eigenrho_Lxyz.dat    Character projected PDOS
! =============================================================


write (74,*) natoms, norbitals_new, nkpoints
write (75,*) natoms, norbitals_new, nkpoints
write (76,*) natoms, norbitals_new, nkpoints
write (77,*) natoms, norbitals_new, nkpoints
write (78,*) natoms, norbitals_new, nkpoints
write (79,*) natoms, norbitals_new, nkpoints

write (80,*) natoms, norbitals_new, nkpoints
write (81,*) natoms, norbitals_new, nkpoints
write (82,*) natoms, norbitals_new, nkpoints
write (83,*) natoms, norbitals_new, nkpoints
write (84,*) natoms, norbitals_new, nkpoints

do iband = 1, norbitals_new
 do  ikpoint = 1, nkpoints
      k = special_k(:,ikpoint)
      write (74,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (75,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (76,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (77,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (78,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (79,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)

      write (80,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (81,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (82,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (83,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)
      write (84,'(f10.5,f10.5,f10.5,f10.5)',advance = 'no')  k(1),k(2),k(3),eigen_k(iband,ikpoint)

      do iatom = 1, natoms
          in1 = imass(iatom)

        imu = 0
        dens  = 0.0
        densS = 0.0
        densP = 0.0
        densD = 0.0
        densPx = 0.0
        densPy = 0.0
        densPz = 0.0


        densDz2 = 0.0
        densDxz = 0.0
        densDyz = 0.0 
        densDxy = 0.0
        densDx2y2 = 0.0

        do issh = 1,nssh(in1)
         l = lssh(issh,in1)
         do lmu = 1, (2*l+1)
          imu = imu + 1
           mmu = imu + degelec(iatom)

            dd = 0.0
            do ineigh = 1, neighn(iatom)
              mbeta = neigh_b(ineigh,iatom)
              jatom = neigh_j(ineigh,iatom)
              in2 = imass(jatom)
              R = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
              dot = k(1)*R(1) + k(2)*R(2) + k(3)*R(3)
              phase = cmplx(cos(dot),sin(dot))
              step1 = phase*(bbnkre(mmu,iband,ikpoint)- ai*bbnkim(mmu,iband,ikpoint))
              do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*(bbnkre(nnu,iband,ikpoint)+ai*bbnkim(nnu,iband,ikpoint))
                dd = dd + real(step2)*s_mat(imu,inu,ineigh,iatom)
              end do ! inu
            end do ! ineigh

              !dens = dens + dd
              if ( l .eq. 0) densS = densS + dd      ! s orbitals
              if ( l .eq. 1) then                   ! p orbitals 
                  densP = densP + dd
                  if ( lmu .eq. 1) densPx = densPx + dd
                  if ( lmu .eq. 2) densPx = densPx + dd
                  if ( lmu .eq. 3) densPy = densPy + dd
              end if
              if ( l .eq. 2) then                  ! d orbitals
                  densD = densD + dd
                  if ( lmu .eq. 1) densDz2   = densDz2   + dd
                  if ( lmu .eq. 2) densDxz   = densDxz   + dd
                  if ( lmu .eq. 3) densDyz   = densDyz   + dd
                  if ( lmu .eq. 4) densDxy   = densDxy   + dd
                  if ( lmu .eq. 5) densDx2y2 = densDx2y2 + dd
              end if
           enddo ! do lmu
        enddo ! do issh

        write (74,'(f15.10)',advance = 'no') densS
        write (75,'(f15.10)',advance = 'no') densP
        write (76,'(f15.10)',advance = 'no') densD
        write (77,'(f15.10)',advance = 'no') densPx
        write (78,'(f15.10)',advance = 'no') densPy
        write (79,'(f15.10)',advance = 'no') densPz

        write (80,'(f15.10)',advance = 'no') densDz2
        write (81,'(f15.10)',advance = 'no') densDxz
        write (82,'(f15.10)',advance = 'no') densDyz 
        write (83,'(f15.10)',advance = 'no') densDxy
        write (84,'(f15.10)',advance = 'no') densDx2y2


        write (1179,'(f15.10)',advance = 'no') dens
      end do  ! iatom
      write (74,*)    
      write (75,*)  
      write (76,*)  
      write (77,*)  
      write (78,*)  
      write (79,*)  

      write (80,*)    
      write (81,*)  
      write (82,*)  
      write (83,*)  
      write (84,*)  
  end do  ! ikpoint
end do  ! iband   

  close(74)
  close(75)
  close(76)
  close(77)
  close(78)
  close(79)

  close(80)
  close(81)
  close(82)
  close(83)
  close(84)



  deallocate(eigenv)

  return
end subroutine writeout_eigenvec

