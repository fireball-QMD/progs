! Copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! excitations.f90
!   Compute duipole transition 
!   form states of interval (frommin,frommax)      
!   to states   of interval (tomin,tomax)
!
! ===========================================================================
! Code written by:
! ===========================================================================
! P. ek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
 subroutine excitations

   use configuration
   use dimensions
   use interactions
   use neighbor_map
   use grid
   use density
   use charges
   use kpoints
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

!Output



! Local Parameters and Data Declaration
! ===========================================================================.

   real, parameter ::  Debye = 0.208194

   interface
    subroutine writeout_xsf (xsfname, message, aa)
     real, dimension (:), pointer, intent (in) :: aa
     character (len=40) xsfname
     character (len=30) message
    end
  end interface

   interface
    subroutine  project_orb(iband,ewfaux)
     integer iband
     real, dimension (:), pointer, intent (out) :: ewfaux
    end
  end interface

! Local Variable Declaration and Description
! ===========================================================================

! grid variables
   integer index
   integer i, j, k
   integer, dimension (3) :: nr

   integer ifrom,ito,ifromlast

   real, target, dimension (:), allocatable :: ewf2
   real, target, dimension (:), allocatable :: ewf1
   real, dimension (:), pointer   :: pewf1
   real, dimension (:), pointer   :: pewf2
   real, target, dimension (:), allocatable :: transmap

! output variables 
   character(40)   :: namewf
   character(4)    :: fromname
   character(4)    :: toname
   real, dimension (:), pointer   :: pmat
   character (len=30) mssg
   integer trannum 

! Dipole variables
   real x,y,z
   real dip_x
   real dip_y
   real dip_z
   real dip_tot
   real dqi

! input file varibales
   integer inputfile    
   real deltaEmin,deltaEmax   
   integer plottransmap             ! should be transition map plot to xsf?
!   integer fromorbmin,fromorbmax    ! range of orbitals from which excite
!   integer toorbmin,toorbmax        ! range of orbitals to which excite

! Debug Time
    real time_start,time_end, time_sum



! Procedure
! ===========================================================================


write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
write (*,*) " !!!!!!!  excitations   !!!!!!! "
write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
write (*,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "

  open ( unit = 10111, file = "Excitations.dat", status = 'unknown' )

  open (unit = 27, file = 'excitations.optional', status = 'unknown')

  read (27,*) plottransmap
  read (27,*) deltaEmin
  read (27,*) deltaEmax
  ! read (27,*) fromorbmin
  ! read (27,*) fromorbmax
  ! read (27,*) toorbmin
  ! read (27,*) toorbmax 

! allocate aux arrays
   allocate ( ewf1(0:nrm-1))
   allocate ( ewf2(0:nrm-1))
   allocate ( transmap(0:nrm-1))

   trannum = 0
! Loop over bands

 time_sum = 0


ifromlast = -1

   do ifrom = 1, norbitals
    do ito = 1, norbitals

!    Write (*,'(1f20.10,1f20.10,1f20.10,1f20.10,1f20.10)') eigen_k(ifrom,1),eigen_k(ito,1),efermi,(eigen_k(ito,1) - eigen_k(ifrom,1)),deltaEmax

    if (  ( eigen_k(ifrom,1) .lt. efermi) .AND.                        &
       &  ( eigen_k(ito,1)   .gt. efermi) .AND.                        &
       &  ( (eigen_k(ito,1) - eigen_k(ifrom,1)) .gt. deltaEmin) .AND.  &
       &  ( (eigen_k(ito,1) - eigen_k(ifrom,1)) .lt. deltaEmax)   ) then
       trannum = trannum + 1
 !   Write (*,*) " ================== DO TRANSITION:   "

 call cpu_time (time_start)

      if(ifromlast .ne. ifrom) then
      ! ewf1 = 0.0d0
      pewf1 => ewf1
      call project_orb(ifrom,pewf1)
      ifromlast = ifrom
      WRITE (*,*) "DEBUG: ifrom changed"
      end if    

      ! ewf2 = 0.0d0
      pewf2 => ewf2
      call project_orb(ito,pewf2)

    write (*,*) ' Dipole transition ',trannum,': ',ifrom,' --> ',ito,' dE= ',eigen_k(ito,1) - eigen_k(ifrom,1)

   transmap = 0.0d0
! calc dipole with the unit cell

   dip_x = 0.0d0
   dip_y = 0.0d0
   dip_z = 0.0d0

!--$OMP PARALLEL PRIVATE (x, y, z,i,j,k,dqi,index) &
!--$OMP & SHARED (pewf1,pewf2,elvec,dip_x,dip_y,dip_z)
  index = 0
!--$OMP DO
   do k = 0, rm3-1
    do j = 0, rm2-1
     do i = 0, rm1-1
      index = k*rm2*rm1+j*rm1+i
      ! dqi = ewf1(index)*ewf2(index) 
      dqi = pewf1(index)*pewf2(index) 
      transmap(index) = dqi 
      x = i*elvec(1,1) + j*elvec(2,1) + k*elvec(3,1)
      y = i*elvec(1,2) + j*elvec(2,2) + k*elvec(3,2)
      z = i*elvec(1,3) + j*elvec(2,3) + k*elvec(3,3)
      dip_x = dip_x + x*dqi
      dip_y = dip_y + y*dqi
      dip_z = dip_z + z*dqi
      ! index =  + index + 1
     enddo ! i
    enddo ! j
   enddo ! k
!--$OMP END DO
!--$OMP END PARALLEL

   dip_x = dip_x * dvol
   dip_y = dip_y * dvol
   dip_z = dip_z * dvol
   dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
  ! write (10111,100) ifrom,ito,dip_tot/Debye,dip_x/Debye,dip_y/Debye,dip_z/Debye
   write (10111,100) ifrom,ito,(eigen_k(ito,1) - eigen_k(ifrom,1)),       &
                  &   dip_tot/Debye,dip_x/Debye,dip_y/Debye,dip_z/Debye,   &
                  &   eigen_k(ito,1), eigen_k(ifrom,1)
! Write out Transmap
   if (plottransmap .eq. 1) then
     write (fromname,'(i4.4)') ifrom
     write (toname,'(i4.4)') ito
     namewf = 'T_'//fromname//'_'//toname//'.xsf'
     pmat => transmap
     mssg = 'wavefunc_3D'
     call writeout_xsf (namewf, mssg, pmat)
   end if

 call cpu_time (time_end)

write (*,*) " Time: ", (time_end - time_start), " [s]"
time_sum = time_sum + (time_end - time_start)

        end if ! deltaEmax

      enddo ! itoorb
   enddo ! ifromorb


write (*,*) " Total time to compute excitations: ", time_sum, " [s]"
write (*,*) " Average time per 1 excitation: ", (time_sum/float(trannum)), " [s]"


   close (10111)

   deallocate (ewf1)
   deallocate (ewf2)
   deallocate (transmap)

! Format Statements
! ===========================================================================
 ! 100 format (i4,1i4,1e20.10,1e20.10,1e20.10,1e20.10)
 ! 100 format (i4,1i4,1f20.10,1e20.10,1e20.10,1e20.10,1e20.10)
   100 format (i7,1i7,1f20.10,1f20.10,1f20.10,1f20.10,1f20.10,1f20.10,1f20.10)
   return
 end subroutine excitations

