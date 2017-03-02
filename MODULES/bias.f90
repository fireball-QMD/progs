module bias
! Description
! module contains variables related to bias voltage
!

! value of bias voltage
  real :: Vbias
! matrix storing Hamiltonian of Vbias
  real, dimension (:, :, :, :), allocatable :: Vbias_mat
! and forces
  real, dimension (:, :, :, :, :), allocatable :: Vbiasp_mat
! z-axis limit defining the region where the bias voltage is applied
  real :: zb1
! z-axis limit defining the region where the bias voltage goes to  zero
  real :: zb0
! note: we supose linear shape of the bias voltage variation between zb1 and zb0
! so we have 3 regions:
! 1. z > zb1             V(z) = Vbias
! 2. zb1 > z > zb0       V(z) = (z-zb0)/(zb1-zb0)*Vbias
! 3. zb0 > z             V(z) = 0

end module bias
