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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! kspace_ordern_fk.f90
! Program Description
! ===========================================================================
!       Faking routines for kspace_ordern.f90 when exact diagonalization is
! used link in these dummies to satisfy the compiling.
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Qin_bcast ()
        write (*,*) '  '
        write (*,*) ' In Qin_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Qin_bcast

        subroutine Qneutral_bcast ()
        write (*,*) '  '
        write (*,*) ' In Qneutral_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Qneutral_bcast

        subroutine ewald_bcast ()
        write (*,*) '  '
        write (*,*) ' In ewald_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine ewald_bcast

        subroutine scf_bcast ()
        write (*,*) '  '
        write (*,*) ' In scf_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine scf_bcast

        subroutine ratom_bcast ()
        write (*,*) '  '
        write (*,*) ' In ratom_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine ratom_bcast

        subroutine ME_max_bcast ()
        write (*,*) '  '
        write (*,*) ' In ME_max_bcast '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine ME_max_bcast

! KSPACE_ORDERN
        subroutine kspace_ordern ()
        write (*,*) '  '
        write (*,*) ' In kspace_ordern! '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine kspace_ordern

        subroutine kspace_ordern_init ()
        write (*,*) '  '
        write (*,*) ' In kspace_ordern_init '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine kspace_ordern_init

        subroutine kspace_ordern_slave ()
        write (*,*) '  '
        write (*,*) ' In kspace_ordern_slave: '
        write (*,*) ' The iordern option is set to 1; however, we should not '
        write (*,*) ' not be in this routine.  The ordern routine was not '
        write (*,*) ' compiled into the executable. Change Makefile and '
        write (*,*) ' re-run or fix iordern in the options.input file. '
        stop
        return
        end subroutine kspace_ordern_slave

        subroutine ordern_init ()
        write (*,*) '  '
        write (*,*) ' In ordern_init '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine ordern_init

! ASSEMBLERS
        subroutine assemble_2c_ordern_final ()
        implicit none
        write (*,*) '  '
        write (*,*) ' In assemble_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_2c_ordern_final

        subroutine assemble_2c_ordern_init ()
        write (*,*) '  '
        write (*,*) ' In assemble_2c_ordern_init '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_2c_ordern_init

        subroutine assemble_3c_ordern_final ()
        implicit none
        write (*,*) '  '
        write (*,*) ' In assemble_3c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_3c_ordern_final

        subroutine assemble_ca_2c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In assemble_ca_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_ca_2c_ordern_final

        subroutine assemble_ca_3c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In assemble_ca_3c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_ca_3c_ordern_final

        subroutine assemble_eh_2c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In assemble_eh_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_eh_2c_ordern_final

        subroutine assemble_lr_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In assemble_lr_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_lr_ordern_final

! DASSEMBLERS
        subroutine Dassemble_ca_2c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_ca_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_ca_2c_ordern_final

        subroutine Dassemble_ca_3c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_ca_3c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_ca_3c_ordern_final

        subroutine Dassemble_eh_2c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_eh_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_eh_2c_ordern_final

        subroutine assemble_ordern_sub_ewald ()
        write (*,*) '  '
        write (*,*) ' In assemble_ordern_sub_ewald '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_ordern_sub_ewald

        subroutine assemble_ordern_sub_dewald ()
        write (*,*) '  '
        write (*,*) ' In assemble_ordern_sub_dewald '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine assemble_ordern_sub_dewald

        subroutine Dassemble_lr_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_lr_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_lr_ordern_final

        subroutine Dassemble_2c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_2c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_2c_ordern_final

        subroutine Dassemble_3c_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In Dassemble_3c_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine Dassemble_3c_ordern_final

        subroutine ewald_energy_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In ewald_energy_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine ewald_energy_ordern_final

        subroutine buildh_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In buildh_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine buildh_ordern_final

! NEIGHBORS
        subroutine common_neighbors_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In common_neighbors_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine common_neighbors_ordern_final

        subroutine find_neigh_max_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In find_neigh_max_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine find_neigh_max_ordern_final

        subroutine neighbors_ordern_final ()
        write (*,*) '  '
        write (*,*) ' In neighbors_ordern_final '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine neighbors_ordern_final

! READFILES
        subroutine readdata_ordern_init ()
        write (*,*) '  '
        write (*,*) ' In readdata_ordern_init '
        write (*,*) ' We should not be in this routine - the ordern '
        write (*,*) ' routine was not compiled into the executable. '
        write (*,*) ' Change Makefile and re-run. '
        stop
        return
        end subroutine readdata_ordern_init
