 module outputs
! This module defines variables controlling output data
! ===========================================================================
! ---------------------------------------------------------------------------
! Toggles read from the output.input file.
! ---------------------------------------------------------------------------
        integer iwrtcdcoefs         ! write out wavefunction coefficients
        integer iwrtcharges         ! write out charges
        integer iwrtdensity         ! write out density matrix
        integer iwrtefermi          ! write out fermi levels
        integer iwrteigen           ! write out eigenvalues
        integer iwrtfpieces         ! write out force pieces
        integer iwrthampiece        ! write out Hamiltonian pieces
        integer iwrtcomponents      ! write out Hamiltonian components
        integer iwrtneigh           ! write out neighbor map
        integer iwrtneigh_com       ! write out common neighbor map
        integer iwrtxyz             ! write out xyz file
        integer iwrtdos             ! write out dos files CGP
        integer iwrtdosng           ! write out dosng (not green) file
        integer iwrthop             ! write out hoppings for STM CGP
        integer iwrtatom            ! write out hoppings for STM CGP
        integer iwrtpop             ! write out population
        integer iwrtHS              ! write out H & S
        integer iwrtvel             ! write out velocities
        integer iwrtNHC             ! write phase space coordinates/NHC Hamil
        integer iwrtden             ! write density projected on the grid
        integer iwrtewf             ! write projected eigenfunctions on the grid
        integer iwrtxsf             ! write xsf-format file (xcrysden)
        integer idensimport         ! importing density file 'rhoin' for projection
        integer iwrtdipole          ! write out dipolar moment
        integer iwrtVibcouplings    ! write out vibronic couplings
! TDSE output options
        integer iwrtpsit            ! write w.f. coeff of each electron in time
        integer iwrtqt              ! write charge occupancy in time
		integer iwrtkvaziband       ! bandstructure projection to different unitcel - prokop
		integer iwrtexcit			! computation of optical dipole transitions <psi_i|d|psi_j> by fermi golden rule on grid - prokop 

 end module outputs
