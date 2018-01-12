!                             @@Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! dftd3_corrections.f90
! Program Description
! ===========================================================================
!     Add DFTD3 corrections 
!
! ===========================================================================
! Code written by:
! Jesús I. Mendieta
! Departamento de físca teórica de la materia condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================

        subroutine dftd3_corrections

        use dftd3_api 
        use configuration
        use interactions
        use energy
        use forces
        use charges, only: nzx
        use options
       

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        type(dftd3_input) :: input
        type(dftd3_calc) :: dftd3
        integer :: iatom,in1,i
        real :: edisp
        integer :: atomic_mass(natoms)
        real :: grads(3, natoms)
        real :: stress(3,3)
        real :: lvs33(3,3)  

! Procedure
! ===========================================================================
        etot_dftd3 = 0.0d0
        ftot_dftd3 = 0.0d0
        
        write(*,*) 'DFTD3 corrections subroutine'
        call dftd3_init(dftd3, input)
        !write(*,*) 'finish dftd3_init'


        if (idftd3 .eq. 1) then
          call dftd3_set_functional(dftd3, dftd3_func, dftd3_version, dftd3_tz)
        endif

        if (idftd3 .eq. 2) then
          call dftd3_set_params(dftd3, dftd3_params, dftd3_version)
        endif
        !write(*,*) 'finish dftd3_set_functional'

        !dftd3_dispersion(this, coords, izp, disp, grads)
        !inputs   coords: coordinate, izp: Atomic number
        !outputs  disp: energy correction, grads: gradient correction

        do i = 1,3
          lvs33(1,i)=a1vec(i)
          lvs33(2,i)=a2vec(i)
          lvs33(3,i)=a3vec(i)
        enddo

        do iatom = 1, natoms
          in1 = imass(iatom)
          atomic_mass(iatom) = nzx(in1)
        enddo

        !atomic positions and lvs vectors in atomic units

        if (icluster .eq. 1) then
          call dftd3_dispersion(dftd3, ratom / 0.52917726, atomic_mass, edisp, grads)
        else 
          call dftd3_pbc_dispersion(dftd3,  ratom / 0.52917726,  atomic_mass,lvs33 / 0.52917726, edisp, grads, stress)
        endif

        !la distancias van en bohr para dftd3: ratom/0.52917726
        write(*,*) 'energy disp correction =', edisp
        !energia en a.u. para pasar a eV * 27.2107
        etot_dftd3 = edisp * 27.2107
        !1write(*,*) 'force disp correcion ='
        !write(*,*) grads
        !fuerzas en a.u. (Eh/ao) para pasar a eV/A * 51.421
        do iatom = 1, natoms
          ftot_dftd3(:,iatom) = - grads(:,iatom) * 51.421
        enddo

        endsubroutine dftd3_corrections
