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

! hampiece.f90
! Program Description
! ===========================================================================
!       This routine writes out pieces of the non-orthogonal hamiltonian.
! 
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine hamtrans ()

  use configuration
  use dimensions
  use interactions
  use neighbor_map
  use charges
  use module_dos
  
  implicit none
  
! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
  integer iatom
  integer jatom
  integer jatom0
  integer in1
  integer in2
  integer ineigh
  integer ineigh0
  integer imu
  integer inu
  integer issh
  integer mbeta
  integer mbeta0
  integer numorb
  integer matom

  real distance
  real, dimension (numorb_max, numorb_max) :: htemp
  real, dimension (numorb_max, numorb_max) :: stemp
  
  character(3) :: ind_aux
  
  real, dimension (3) :: dvec
 
! Procedure
! ===========================================================================
  write (*,*) ' '
  write (*,*) ' In hamtrans.f90 '
  write (*,*) ' Writing out pieces of the Hamiltonian. '
  write (*,*)
  write (*,*) ' The information needed to make a TB model or for doing '
  write (*,*) ' the transport calculation '
  write (*,*) ' Here we are storing non-orthogonal hamiltonian and '
  write (*,*) ' overlap matices between neighbors of each atom.'

! write intraatomic non-orthogonal interaction for transporte code
  do iatom = 1, natoms
     
     in1 = imass(iatom)
     
     htemp = 0.0d0           
     
     if (iatom .gt. 9 .and. iatom .lt. 100) write (ind_aux,'(i2)') iatom
     if (iatom .gt. 99) write (ind_aux,'(i3)') iatom
     if (iatom .le. 9)  write (ind_aux,'(i1)')iatom
     
     open (unit=10, file='Hamilt_'//ind_aux, status='unknown')
     write (10,*) 'Nonorthogonal intraatomic interactions of atom:', iatom
     write (10,*) num_orb(in1)
     write (10,*) efermi


! sum contributuion from regular neighbors
     matom = neigh_self(iatom)
     do imu = 1, num_orb(in1) 
        do inu = 1, num_orb(in1) 
           htemp(imu,inu) = htemp(imu,inu) + h_mat(imu,inu,matom,iatom) 
        enddo ! inu
     enddo ! imu

! sum contributuion from PP-neighbors
     matom = neighPP_self(iatom)
     do imu = 1, num_orb(in1) 
        do inu = 1, num_orb(in1) 
           htemp(imu,inu) = htemp(imu,inu) + vnl(imu,inu,matom,iatom) 
        enddo ! inu
     enddo ! imu

! write non-orthogonal intraatomic interactions
     do imu = 1,num_orb(in1)     
        write (10,'(15f16.8)') htemp(imu,1:num_orb(in1))
     end do
     
     close (10)

  enddo ! do iatom


! loop over atoms         
  do iatom = 1, natoms

     numorb = 0
     in1 = imass(iatom)
     do issh = 1, nssh(in1)
        numorb = numorb + 2*lssh(issh,in1) + 1
     end do

! write interaction for transport code
     if (iatom .gt. 9 .and. iatom .lt. 100) write (ind_aux,'(i2)') iatom
     if (iatom .gt. 99) write (ind_aux,'(i3)') iatom
     if (iatom .le. 9)  write (ind_aux,'(i1)')iatom

! open existing hamil_* file
     open (unit=10,file='Hamilt_'//ind_aux,position='append',status='old')
     write (10,*) 'Interatomics interactions'
     write (10,*) 'vecinos:', neighn_tot(iatom)-1

! open new overlap_* file
     open (unit=20, file='Overlap_'//ind_aux, status='unknown')
     write (20,*) 'Overlaps of atom:', iatom
     write (20,*) num_orb(in1)
     write (20,*) 'vecinos:', neighn_tot(iatom)-1

! loop over total list of neighbors 
     do ineigh = 1, neighn_tot(iatom)
        
        jatom = neighj_tot(ineigh,iatom) 
        mbeta = neighb_tot(ineigh,iatom) 
        in2 = imass(jatom) 
        dvec(:) = ratom(:,jatom) + xl(:,mbeta) - ratom(:,iatom) 

        htemp = 0.0d0
        stemp = 0.0d0

! loop over regular list of neighbors
        do ineigh0 = 1, neighn(iatom)
           jatom0 = neigh_j(ineigh0,iatom) 
           mbeta0 = neigh_b(ineigh0,iatom) 
           
! find identical neighbors
           if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then 
            do imu = 1, num_orb(in1) 
              do inu = 1, num_orb(in2) 
                htemp(imu,inu) = htemp(imu,inu) + h_mat(imu,inu,ineigh0,iatom) 
                stemp(imu,inu) = s_mat(imu,inu,ineigh0,iatom) 
              enddo ! inu
            enddo ! imu
           endif
        enddo ! do ineigh0


! loop over PP list of neighbors
        do ineigh0 = 1, neighPPn(iatom)
           jatom0 = neighPP_j(ineigh0,iatom) 
           mbeta0 = neighPP_b(ineigh0,iatom) 
           
! find identical neighbors
           if (jatom .eq. jatom0 .and. mbeta .eq. mbeta0) then 
            do imu = 1, num_orb(in1) 
             do inu = 1, num_orb(in2) 
              htemp(imu,inu) = htemp(imu,inu) + vnl(imu,inu,ineigh0,iatom) 
             enddo ! inu
            enddo ! imu
           endif
        enddo ! do ineigh0


! write in format of transporte code
! skip intraatomic interactions (yet written in header)
        if (jatom .eq. iatom .and. mbeta .eq. 0) then
        else

! write interatomic interaction
           write (10,700) 'vecino:', jatom, mbeta, (dvec(1:3)/lattice)
           write (10,*) num_orb(in2)            
           do imu = 1,num_orb(in1)  
              write (10,'(15f16.8)') htemp(imu,1:num_orb(in2))
           end do
! write overlap 
           write (20,700) 'vecino:', jatom, mbeta, (dvec(1:3)/lattice)
           write (20,*) num_orb(in2)            
           do imu = 1,num_orb(in1)             
              write (20,'(15f16.8)') stemp(imu,1:num_orb(in2))
           end do                                

         end if ! if(jatom.eq.iatom)

      end do  ! do ineigh

! close files
      close (10)
      close (20)
      
   end do ! do iatom


!!! write the file struc.inp (structural properties) for the STM calculation
   open (unit = 11, file = 'struc.inp', status = 'unknown')
! write number of atoms
   write (11,*) natoms
! write number of species
   write (11,*) nspecies
! write number of orbitals of each specie
   do in1 = 1, nspecies
      write (11,*) in1, num_orb(in1) 
   enddo
! write max. number of neighbors
   write (11,*) num_neig_maxtot
! write atomic coordinates and atomic type
   if (ishiftO .eq. 1) then 
      do iatom = 1, natoms
         write (11,200) (ratom(1:3,iatom)-shifter(1:3))/lattice,imass(iatom)
      end do
   else 
      do iatom = 1, natoms
         write (11,200) ratom(1:3,iatom)/lattice,imass(iatom)
      end do
   endif

! write lattice vector
   write (11,300) a1vec(1:3)/lattice 
   write (11,300) a2vec(1:3)/lattice 
   write (11,300) a3vec(1:3)/lattice 
   close (11)


        write (*,100)
 
! Format Statements
! ===========================================================================
100     format (75('*'))
200     format (3f18.8, 2x, i5)
300     format (3f18.8)
700     format (a8,2x,i4,2x,i4,6x,3f16.8)

        close (unit = 11)
        return
      end subroutine hamtrans
