!############################################
! velocity Verlest aldoritm for n time step
!############################################
 subroutine vverlet(ratoms, vatoms, fatoms, matoms, dt)

   implicit none

   real, dimension (:, :), intent(in)			:: ratoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: vatoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: fatoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: matoms	! positions of atoms
   integer,			intent(in)			:: dt		! time step
   
   integer :: iatoms	! atom's iteratot
   integer :: natoms	! number o atoms (size of ratoms)
   real :: dt2, dt22	! parameters fo verlet algorytm
   
   natoms = size(ratom,2)
  
   if( natoms .eq. size(vatoms,2) .and. natoms .eq. size(fatoms,2) .and. size(matoms,1)) then 
    
    dt2 = dt*0.5
    dt22 = dt*dt*0.5
    do iatoms = 1, natoms
     vatoms(:,iatoms) = vatoms(:,iatoms)  + dt2*fatoms(:,iatoms)/matoms(iatoms)
     ratoms(:,iatoms) = ratoms(:,iatoms) + dt*vatoms(:,iatoms) + dt22*fatoms(:,iatoms)/matoms(iatoms)
     vatoms(:,iatoms) = vatoms(:,iatoms) + dt2*fatoms(:,iatoms)/matoms(iatoms)
    end do
    
  else
  
  	write(*,*) "vverlet's error : sizes of input arrays are not correct"
  
  end if
      
 end subroutine vverlet
 
 
 
 
 !############################################
 ! velocity Verlest aldoritm for 0 time step
 !############################################
 subroutine vverlet0(ratoms, vatoms, fatoms, matoms, dt)
   implicit none
   real, dimension (:, :), intent(in)			:: ratoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: vatoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: fatoms	! positions of atoms
   real, dimension (:, :), intent(in)			:: matoms	! positions of atoms
   integer,			intent(in)			:: dt		! time step
   
   
   integer :: iatoms	! atom's iteratot
   integer :: natoms	! number o atoms (size of ratoms)
   real :: dt2, dt22	! parameters fo verlet algorytm
   
   natoms = size(ratom,2)
  
   if( natoms .eq. size(vatoms,2) .and. natoms .eq. size(fatoms,2) .and. size(matoms,1)) then 
    
    dt2 = dt*0.5
    dt22 = dt*dt*0.5
    do iatoms = 1, natoms
     ratoms(:,iatoms) = ratoms(:,iatoms) + dt*vatoms(:,iatoms) + dt22*fatoms(:,iatoms)/matoms(iatoms)
     vatoms(:,iatoms) = vatoms(:,iatoms) + dt2*fatoms(:,iatoms)/matoms(iatoms)
    end do
    
  else
  
  	write(*,*) "vverlet's error : sizes of input arrays are not correct"
  
  end if
      
 end subroutine vverlet0
 

 

 subroutine move_image (natoms, img)
   use configuration
   use neb
   use interactions

   implicit none
   integer,			intent(in)			:: natoms		! number of atoms in img'th image
   integer,			intent(in)			:: img			! 


   if(itime .le. nimg_neb)
   	! be carefull about velocitis. you should update it with ratom to the ratom_neb and vel_neb
   	call vverlet0(ratom,vatom,Fneb,imass,dt_neb)
   else
   	call vverlet(ratom,vatom,Fneb,imass,dt_neb)
   end if

 end subroutine move_image