module noseHoover

   integer  nnos                     ! number of chains (default is four)
  ! real xi(m),pxi(m),qh(m),sc(m),psc(m)
   real :: kT,gkT,gT

   real, allocatable, dimension(:) :: xi     ! friction coefficients
   real, allocatable, dimension(:) :: v_xi   ! derivative friction coefficients
   real, allocatable, dimension(:) :: Q_i    ! nose-hoover masses
   real, allocatable, dimension(:) :: G_i 
   real, allocatable, dimension(:) :: omega  ! characteristic frequencies
   logical debug 

end module 
