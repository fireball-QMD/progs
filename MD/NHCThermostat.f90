subroutine NHCThermostat(dt,natoms,xmass,vatom)
use noseHoover
implicit none

integer, intent(in):: natoms
real, intent(in), dimension(natoms) :: xmass
real, intent(inout), dimension(3,natoms) :: vatom
real, intent(in) :: dt
! this routine does the nose-hoover part of the integration from t=0 to t=dt/2
! get the total kinetic energy

! it is taken from (G.J. Martyna et.al. Mol Phys 1996 87(5) 1117-1157)

! glogs G_i of the text
! qmass Q_i of the text
! gkT = number of degrees of freedom times kT
! nyosh n_ys ""
! nresn n_c ""
! dt delta_t ""
! dt2 detla_t/2 
! dt22 delta_t^2/2 
! wdti w delta_t/n_c ""
! wdti2 w delta_t/2n_c ""  
! etc.

real akin,scale,www,aa
integer :: i,iresn,iyosh,inos,nresn=1,nyosh=5 ! can change this to other values
real, dimension(5) :: wdti2, wdti4, wdti8

     ! write(*,*)'nyosh',nyosh
     ! write(*,*)'wdti8',wdti8

if(nyosh == 1) then

 wdti2(1) = dt / ( 2.d0 * nresn )
 wdti4(1) = wdti2(1) * 0.5d0
 wdti8(1) = wdti4(1) * 0.5d0

else if(nyosh == 3) then

 www = 1.d0 / ( 2.d0 - 2.d0**(1.d0/3.d0) )

 wdti2(1) = www * dt / ( 2.d0 * nresn )
 wdti2(2) = ( 1.d0 - 2.d0 * www ) * dt / ( 2.d0 * nresn )
 wdti2(3) = wdti2(1)

 wdti4(1) = wdti2(1) * 0.5d0
 wdti4(2) = wdti2(2) * 0.5d0
 wdti4(3) = wdti4(1)

 wdti8(1) = wdti4(1) * 0.5d0
 wdti8(2) = wdti4(2) * 0.5d0
 wdti8(3) = wdti8(1)

else if(nyosh == 5) then

 www = 1.d0 / ( 4.d0 - 4.d0**(1.d0/3.d0) )

! the extra factor of 2 is because we are doing dt/2.0         
 wdti2(1) = www * dt / ( 2.0*nresn )
!wdti2(1) = www * dt / ( nresn )
 wdti2(2) = wdti2(1)
 wdti2(3) = ( 1.d0 - 4.d0 * www ) * dt / ( 2.0*nresn )
!wdti2(3) = ( 1.d0 - 4.d0 * www ) * dt / ( nresn )
 wdti2(4) = wdti2(1)
 wdti2(5) = wdti2(1)

 wdti4(1) = wdti2(1) * 0.5d0
 wdti4(2) = wdti4(1)
 wdti4(3) = wdti2(3) * 0.5d0
 wdti4(4) = wdti4(1)
 wdti4(5) = wdti4(1)

 wdti8(1) = wdti4(1) * 0.5d0
 wdti8(2) = wdti8(1)
 wdti8(3) = wdti4(3) * 0.5d0
 wdti8(4) = wdti8(1)
 wdti8(5) = wdti8(1)

end if

!write(*,*)'www',www
!write(*,*)'wdti2',wdti2
!write(*,*)'wdti4',wdti4
!write(*,*)'wdti8',wdti8
scale=1.0
call get_enk(natoms,xmass,vatom,akin)

! start the multiple time step procedure
do iresn=1,nresn
 do iyosh=1,nyosh
! update the thermostat velocities
! write(*,*)'v_xi(nnos) before',v_xi(nnos)
! write(*,*)'G_i(nnos) before',G_i(nnos)
! write(*,*)'wdti4(iyosh) before',wdti4(iyosh)

  G_i(nnos) = (Q_i(nnos-1)*v_xi(nnos-1)*v_xi(nnos-1) - kT)/Q_i(nnos)
  v_xi(nnos) = v_xi(nnos) + G_i(nnos)*wdti4(iyosh)

  do inos=nnos-1,2,-1
   aa=exp(-wdti8(iyosh)*v_xi(inos+1))
   v_xi(inos) = v_xi(inos)*aa
   G_i(inos) = (Q_i(inos-1)*v_xi(inos-1)*v_xi(inos-1) &
     - kT)/Q_i(inos)
   v_xi(inos) = v_xi(inos) + wdti4(iyosh)*G_i(inos)
   v_xi(inos) = v_xi(inos)*aa
  end do  

  aa=exp(-wdti8(iyosh)*v_xi(2))
  v_xi(1) = v_xi(1)*aa
  G_i(1)=(2*akin - gkT)/Q_i(1)
  v_xi(1) = v_xi(1) + wdti4(iyosh)*G_i(1)
  v_xi(1) = v_xi(1)*aa

! update the particle velocities
  aa=exp(-wdti2(iyosh)*v_xi(1))
  scale=scale*aa
  akin = akin*aa*aa

  if(debug) then
   write(*,*)'1. v_xi',v_xi
   write(*,*)'G_i',G_i
   write(*,*)'aa',aa
  end if

! update the thermostat positions
  do inos=1,nnos
   xi(inos)=xi(inos) + v_xi(inos)*wdti2(iyosh)
  end do

! now apply the other thermostat operator
  aa=exp(-wdti8(iyosh)*v_xi(2))
  v_xi(1) = v_xi(1)*aa
  G_i(1)=(2*akin - gkT)/Q_i(1)
  v_xi(1) = v_xi(1) + wdti4(iyosh)*G_i(1)
  v_xi(1) = v_xi(1)*aa

! update the thermostat velocities
  do inos=2,nnos-1
   aa = exp(-wdti8(iyosh)*v_xi(inos+1))
   v_xi(inos) = v_xi(inos)*aa
   G_i(inos) = (Q_i(inos-1)*v_xi(inos-1)*v_xi(inos-1) &
     - kT)/Q_i(inos)
   v_xi(inos) = v_xi(inos) + wdti4(iyosh)*G_i(inos)
   v_xi(inos) = v_xi(inos)*aa
  end do
  G_i(nnos) = (Q_i(nnos-1)*v_xi(nnos-1)*v_xi(nnos-1) &
     - kT)/Q_i(nnos)
  v_xi(nnos) = v_xi(nnos) + G_i(nnos)*wdti4(iyosh)
 end do
end do

!update the particle velocities
vatom = vatom*scale
if(debug) then
 write(*,*)'scale',scale
 write(*,*)'vatom',vatom
 write(*,*)'G_i',G_i
 write(*,*)'2. v_xi',v_xi
 write(*,*)'xi',xi
end if

end subroutine
