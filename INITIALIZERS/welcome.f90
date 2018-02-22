! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega
! Brigham Young University - Hao Wang

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


! welcome.f90
! Program Description
! ===========================================================================
!       This routine prints out the welcome and informs the user the
! dimensionality of everything.
!
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
        subroutine welcome
        implicit none

! Argument Declaration and Description
! ===========================================================================

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

! Procedure
! ===========================================================================
        write (*,*) '               Welcome to FIREBALL QMD !     '
        write (*,*) '          A fast local orbital QMD Package    '
        write (*,*) '              Version  1.2.0 17/5/2010    '
        write (*,*) '     brought to you by the Fireball Committee: '
        write (*,*) '  '
        write (*,*) '                 James P. Lewis '
        write (*,*) '        Department of Physics and Astronomy '
        write (*,*) '             Brigham Young University '
        write (*,*) '  '
        write (*,*) '                 Otto F. Sankey '
        write (*,*) '       Department of Physics and Astronomy '
        write (*,*) '           Arizona State University '
        write (*,*) '  '
        write (*,*) '                  Jose Ortega '
        write (*,*) '             Departmento de Fisica '
        write (*,*) '       Teorica de la Materia Condensada '
        write (*,*) '        Universidad Autonoma de Madrid '
        write (*,*) '  '
        write (*,*) '                Pavel Jelinek '
        write (*,*) '             Institute of Physics '
        write (*,*) '            Prague, Czech Republic '
        write (*,*) '  '

        write (*,*) '  '
        write (*,*) '             with contributions from: '
        write (*,*) '  '
        write (*,*) ' Alex A. Demkov (Motorola Physical Sciences Research Labs)'
        write (*,*) ' Gary B. Adams (Arizona State University) '
        write (*,*) ' Jian Jun Dong (Arizona State University) '
        write (*,*) ' David A. Drabold (Ohio University) '
        write (*,*) ' Peter A. Fedders (Washington University) '
        write (*,*) ' Kurt R. Glaesemann (University of Utah) '
        write (*,*) ' Kevin Schmidt (Arizona State University) '
        write (*,*) ' Spencer Shellman (University of Utah) '
        write (*,*) ' John Tomfohr (Arizona State University) '
        write (*,*) ' Hao Wang (Brigham Young University) '
        write (*,*) ' Daniel G. Trabada (Universidad Autonoma de Madrid) '
        write (*,*) ' Jesus I. Mendieta-Moreno (Universidad Autonoma de Madrid) '

        write (*,*) '  '
        write (*,*) '     Latest version Feb. 2018 '
        write (*,*) '  '
        write (*,*) '         See Copyright information: '
        write (*,*) '           !!!Proprietory Code!!! '
        write (*,*) '  '
        write (*,*) ' Usable only with permission from the Fireball executive '
        write (*,*) ' committee. This program is NOT, under any circumstances, '
        write (*,*) ' to be transfered to an unauthorized user. '
        write (*,*) '  '

        write (*,*) ' This source code is completely dynamical in memory. '

! Format Statements
! ===========================================================================

        return
        end
