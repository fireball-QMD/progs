subroutine get_enk(natoms,xmass,vatom,akin)
use constants_fireball
implicit none
integer, intent(in):: natoms
real, intent(in), dimension(natoms) :: xmass
real, intent(in), dimension(3,natoms) :: vatom
real, intent(out) :: akin
integer iatom
akin=0.0
do iatom = 1, natoms
 akin = akin + (0.5d0/fovermp)*xmass(iatom)*(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
end do
!write(*,*)'akin',akin
return
end subroutine
