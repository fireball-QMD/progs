 module energy 
! This module defines variables of energies   
! ===========================================================================

! Energy and force contributions and terms 
        real atomic_energy
        real etot
        real ebs
        real etotnew
        real etotold
        real etotper
        real etotxc_1c
        real getot
        real getot_initial
        real getotper
        real deltaE
        real deltaFmax
        real uiiuee
        real uxcdcc
        real uxcdcc_sn
        real uxcdcc_hf
        real uxcdcc_ols 
! Kohn-Sham
        real uxcdcc_ks 
        real uhdcc_ks
! Ext Hubbard
        real ehxcc
        real ehcoolc
        real Umuxc_1c
        real Uexc_1c
! QM/MM energy
        real eqmmm

 end module energy 
