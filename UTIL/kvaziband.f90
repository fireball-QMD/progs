



        subroutine kvaziband ()
!        use charges
        use configuration
!        use constants
        use density
        use dimensions
        use interactions
!        use neighbor_map
!        use kpoints
        implicit none

! Local Parameters and Data Declaration
! ===========================================================================

   COMPLEX, PARAMETER :: ai = (0.0,1.0) ! imaginary unit

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iband
        integer imu,mmu
        integer in1
        integer issh
        integer i

  complex projection,kwave
  real dot
  real normWave, normState


 real Emax, Emin
 integer nkvertexs, nkps
 integer ivtx, ikp

 real     ,dimension (:,:), allocatable :: kvertexs 
 real     ,dimension (:,:), allocatable :: kp            
 integer ,dimension (:),   allocatable :: lineNs   




! Procedure
! ===========================================================================
! Initialize some things


        write (*,*) ' ****************************************************** '
        write (*,*) '                   Computing Kvazibands                 '
        write (*,*) '                 = faking bandstructre :-)              '    
        write (*,*) ' ****************************************************** '



! ========= start READ INPUTS ============

  open (unit = 30, file = 'kvaziband.dat', status = 'unknown')
 
  open (unit = 27, file = 'kvaziband.optional', status = 'unknown')
! read energy range
   read (27,*) Emin,Emax
! read lattice
!   read (27,*) basecell(1,1),basecell(2,1),basecell(3,1)
!   read (27,*) basecell(1,1),basecell(2,1),basecell(3,1)
!   read (27,*) basecell(1,1),basecell(2,1),basecell(3,1)

! read nuber of samples and vertexes
   read (27,*) nkps, nkvertexs
   ALLOCATE ( kvertexs (3,nkvertexs) )
   ALLOCATE ( kp       (3,nkps)     )
   ALLOCATE ( lineNs   (nkvertexs-1) )
! read lines and vertexes
   do ivtx = 1,nkvertexs-1
     read (27,*) kvertexs(1,ivtx),kvertexs(2,ivtx),kvertexs(3,ivtx)
     read (27,*) lineNs(ivtx)
   enddo
   read (27,*) kvertexs(1,nkvertexs),kvertexs(2,nkvertexs),kvertexs(3,nkvertexs)
   
! generate kpoints from lines
   ikp = 1
   kp(:,1) = kvertexs(:,1)
   do ivtx = 1,nkvertexs-1
     Write (*,'(" vertex ",3f10.5)') kvertexs(1,ivtx),kvertexs(2,ivtx),kvertexs(3,ivtx)
     do i = 1,lineNs(ivtx)
     kp(:,ikp) = ((lineNs(ivtx)-i)*kvertexs(:,ivtx) + (i)*kvertexs(:,ivtx+1))/float(lineNs(ivtx)) 
      Write (*,'(3f10.5)') kp(1,ikp),kp(2,ikp),kp(3,ikp)
     ikp=ikp+1
     end do
   enddo

! ========= end READ INPUTS ============


do iband = 1, norbitals
 ! do ikpoint = 1, norbitals  ?????????
  Write (*,*) "DEBUG ", EMAX, EMIN, eigen_k(iband,1)
    if (  ( eigen_k(iband,1) .lt. Emax) .AND. ( eigen_k(iband,1) .gt. Emin) ) then
       Write (*,*) iband, eigen_k(iband,1) 


 write (30,'(2x,f20.10)',advance='no')   eigen_k(iband,1)
 do ikp = 1,nkps

   Projection  = 0.0d0
   NormState   = 0.0d0
   NormWave    = 0.0d0 
     do iatom = 1, natoms
          in1 = imass(iatom)
    dot = (kp(1,ikp)*ratom(1,iatom) + kp(2,ikp)*ratom(2,iatom) + kp(3,ikp)*ratom(3,iatom))   !  no 6.28318531!!!      exp(ikr), "Pi" is already inside "k"
          kwave = cos(dot) + ai*sin(dot)
          do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)

                Projection  = Projection +    kwave * bbnkre(mmu,iband,1)
                                        ! +   kwave * ( bbnkre(mmu,iband,1) + ai*bbnkim(mmu,iband,1)  )
                NormState   = NormState  +    abs(kwave * CONJG(kwave))
                NormWave    = NormWave   +    bbnkre(mmu,iband,1)*bbnkre(mmu,iband,1)
                                        ! +   bbnkim(mmu,iband,1)*bbnkim(mmu,iband,1)
          end do ! imu
      end do  ! iatom
   write (30,'(2x,f20.10)',advance='no')  abs(Projection)

end do ! ikp
write (30,*)

 end if ! EMAX > eigen_k(iband,1) > EMIN
end do ! iband



! ======== Finalizing subroutine ==========

   deALLOCATE ( kvertexs )
   deALLOCATE ( kp       )
   deALLOCATE ( lineNs   )

  close (27)
  close (30)


  return
end subroutine kvaziband

