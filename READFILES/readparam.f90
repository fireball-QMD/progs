! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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


! readparam.f90
! Program Description
! ===========================================================================
! A subroutine that reads in fireball.param file all directives for fireball
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readparam ()

        use outputs
        use options
        use configuration
        use kpoints
        use MD
        use scf
        use charges
        use barrier
        use nonadiabatic
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        character (len = 30) str
        logical isfile
        logical isstr

! Namelist for options control
  Namelist /output/ iwrtcdcoefs, iwrtcharges, iwrtdensity, iwrteigen,       &
  					 iwrtefermi, iwrtfpieces, iwrthampiece, iwrtcomponents, &
  					 iwrtneigh, iwrtneigh_com, iwrtxyz, iwrtdos, iwrthop,   &
  					 iwrtatom, iwrtpop, iwrtHS, iwrtvel, iwrtden, iwrtewf,  &
  					 iwrtxsf, idensimport, iwrtpsit, iwrtqt, iwrtkvaziband, iwrtexcit

  Namelist /option/ iharris, idogs, ihubbard, ihorsfield, imcweda, igsn, iquench, &
  					 iqout, qstate, icluster, iensemble, ifixcharge, ifixneigh,   &
  					 iumbrella, ivdw, ibarrier, iimage, idynmat, iharmonic, iks,  &
  					 iconstraints, iendtemp, ineb, itrans, basisfile, lvsfile,    &
  					 kptpreference, acfile, xvfile, nstepi, nstepf, dt,           &
  					 T_initial, T_final, max_scf_iterations, bmix, sigmatol,      &
  					 tempfe, itdse, ibias, rescal, xyz2line, imdet, iProjWF, nddt,&
					 igap, ialgmix, iclassicMD, icDFT, iqmmm, idipole, iephc,     &
                                         idftd3, dftd3_func, dftd3_version, dftd3_tz, dftd3_s6,       &
                                         dftd3_rs6, dftd3_s18, dftd3_rs18, dftd3_alp, mix_embedding,  &
                                         cut_embedding

! Procedure
! ===========================================================================

! default settings
! ------  DEFAULT OPTIONS  ------
        iharris = 0
        idogs = 1
        ihubbard = 0
        ihorsfield = 0
        imcweda = 1
        iks = 0
        igsn = 0
        iqout = 1
        qstate = 0.0d0
        iquench = -1
        icluster = 0
        iensemble = 0
        ifixcharge = 0
        ifixneigh = 0
  		iumbrella = 0
  		ibarrier = 0
  		ivdw = 0
  		iimage = 0
  		idynmat = 0
  		iharmonic = 0
  		iconstraints(1) = 0
        iconstraints(2) = 1
        iconstraints(3) = 1
        iconstraints(4) = 1
  		iendtemp = 0
  		ineb = 0
  		itrans = 0
		basisfile = 'input.bas'
		lvsfile = 'input.lvs'
		kptpreference = 'input.kpts'
		acfile = 'ac.dat'
		xvfile = 'xv.dat'
		nstepi = 1
  		nstepf = 100
  		dt = 0.5d0
  		T_initial = 0.0d0
  		T_final = 10.0d0
        max_scf_iterations = 200
        bmix = 0.04d0
        ialgmix = 1
        sigmatol = 1.0E-8
        tempfe = 100.0d0
        itdse = 0
        ibias = 0
        rescal = 1.0d0
        xyz2line = 1
!JOm-add
        imdet = 0
        iProjWF = 0 
        nddt = 1000
!CHROM	classic potential in interactions in MD or GO
		iclassicMD = 0
!END CHROM

! GAP ENRIQUE-FF
	igap = 0
! cDFT
        icDFT = 0

! these options are currently out of order
! to switch them on, you have to list them in the param_option Namelist
        ispin = 0
        ipathintegral = 0
  		ireducekpts = 0
  		iordern = 0
  		ithermoint = 0
  		igauss = 0

! QM/MM 
        iqmmm = 0
        mix_embedding = 0
        cut_embedding = 999
! e-ph coupling 
        iephc = 0
! DFTD3 corrections
        idftd3 = 0
        dftd3_func = 'b-lyp'
        dftd3_version = 4
        dftd3_tz = .false.
        dftd3_s6 = 1.0d0
        dftd3_rs6 = 0.4298
        dftd3_s18  = 2.6996
        dftd3_rs18  = 4.2359
        dftd3_alp = 14.0d0 
! Long range dipole
        idipole = 0
! ------  DEFAULT OUTPUTS  ------
        iwrtcdcoefs = 0
        iwrtcharges = 0
        iwrtdensity = 0
        iwrteigen = 0
        iwrtefermi = 0
        iwrtfpieces = 0
        iwrthampiece = 0
        iwrtcomponents = 0
        iwrtneigh = 0
        iwrtneigh_com = 0
        iwrtxyz = 0
        iwrtdos = 0
        iwrthop = 0
        iwrtatom = 0
        iwrtpop = 0
        iwrtHS = 0
        iwrtvel = 0
        iwrtden = 0
        iwrtewf = 0
        iwrtxsf = 0
        idensimport = 0
        iwrtpsit = 0
        iwrtqt = 0


        inquire (file = initfile, exist = isfile)
! file fireball.in exists so let's read it
        if (isfile) then

         write (*,*) '  '
         write (*,100)
         write (*,*) ' Now reading file fireball.in '
         write (*,*) '  '
         write (*,*) '  '

! open param file
         open (unit = 60, file = initfile, status = 'old')
! section OPTIONS
! check if the section OPTIONS is listed
         str = '&OPTION'
         call issection (str,isstr)
         if (isstr) read (60, NML=option)

! reset positon of the cursor in the file
         rewind 60

! section OUTPUT
! check if the section OUTPUT is listed
         str = '&OUTPUT'
         call issection (str,isstr)
         if (isstr) read (60, NML=output)

! close unit
         close (unit = 60)

        else

         write (*,*) '  '
         write (*,100)
         write (*,*) ' Fireball.in file does not exist. Default settings used. '
         write (*,*) '  '

        endif ! if (isfile)

! checkout the data consistency
        call checksum_options ()

! SECTION OPTIONS
        write (*,*) ' The name of the basisfile: '
        write (*,201) basisfile

        write (*,*) '  '
        write (*,100)
        write (*,*) ' The name of the lattice vector file: '
        write (*,201) lvsfile

        write (*,*) '  '
        write (*,100)
        write (*,*) ' The name of the k-points file: '
        write (*,201) kptpreference

        write (*,*) '  '
        write (*,100)
        write (*,*) ' The initial and final time steps for this '
        write (*,*) ' simulation (e.g. 101,200). N O T E: setting the initial '
        write (*,*) ' time step to a value less than or equal to zero causes '
        write (*,*) ' this program to NOT compute any any derivatives or '
        write (*,*) ' forces. Thus no simulation, and no calls to forces.'
        write (*,*) '  '
        write (*,*) ' nstepi, nstepf = ', nstepi, nstepf

        iforce = 1
        if (nstepi .le. 0) then
         write (*,*) '  '
         write (*,*) ' YOU HAVE CHOSEN NOT TO DO A SIMULATION. '
         iforce = 0
        end if

        write (*,100)
        write (*,*) '  '
        write (*,*) ' Input time step in fs (10**-15 sec): '
        write (*,*) ' dt = ', dt

        write (*,*) '  '
        write (*,*) ' The position-velocity filename: XXX.xv. '
        write (*,*) ' This file is input and output if nstepi .ne. 1.'
        write (*,201) xvfile

        write (*,*) '  '
        write (*,*) ' The file name for the acfile: XXX.ac. '
        write (*,*) ' (acfile records values of acceleration and other '
        write (*,*) '  time derivatives for last executed step.) '
        write (*,201) acfile

        write (*,*) ' You can have charge transfer in an SCF theory, or you '
        write (*,*) ' can run fixed charges as in the Harris functional. The '
        write (*,*) ' fixed charges are set by the info.dat file. These same '
        write (*,*) ' charges are the initial guess for the SCF charges. '
        write (*,*) '  '

        write (*,*) ' Input itheory - '
        write (*,*) ' This describes the level of theory that you want.'
        write (*,*) '  0 => Harris '
        write (*,*) '  1 => DOGS '
        write (*,*) '  2 => extended-Hubbard '
        write (*,*) '  3 => Kohn-Sham '
        write (*,*) ' itheory = ', itheory
        if (itheory .lt. 0 .or. itheory .gt. 3) then
         write (*,*) ' This selection for itheory is not valid! Change to an '
         write (*,*) ' appropriate value in the options.input file.'
         write (*,*) ' I am afraid that we are going to have to stop. '
         stop
        end if

! option 3: itheory_xc
        write (*,*)
        write (*,*) ' Input itheory_xc - '
        write (*,*) ' This describes the exchange-correlation theory used. '
        write (*,*) '  0 => Horsfield '
        write (*,*) '  1 => generalized Sankey-Niklewski '
        write (*,*) '  2 => McWEDA '
        write (*,*) ' itheory_xc = ', itheory_xc
        if (itheory_xc .lt. 0 .or. itheory_xc .gt. 3) then
         write (*,*) ' This selection for itheory_xc is not valid! Change to '
         write (*,*) ' an appropriate value in the options.input file. '
         write (*,*) ' I am afraid that we are going to have to stop. '
         stop
        end if

! option 5: iquot
        write (*,*) '  '
        write (*,*) ' Would you prefer to use Lowdin charges or Mulliken '
        write (*,*) ' charges?  Enter 1 - Lowdin or 2 - Mulliken: '
        write (*,*) ' 3 - Natural population analysis: '
        write (*,*) ' iqout = ', iqout
        if (itheory .eq. 2 .and. iqout .ne. 2) then
         write (*,*) ' You are using the extended-hubbard approximation. '
         write (*,*) ' As a result, you must use Mulliken charges! '
         write (*,*) ' Setting iqout = 2! '
         iqout = 2
        end if

! option 6: qstate
        write (*,*) '  '
        write (*,*) ' The option exists for changing the charge state of the '
        write (*,*) ' system so that it is non-neutral. Insert qstate, the '
        write (*,*) ' charge state that you want. '
        write (*,*) ' qstate = ', qstate

! option 7: iquench
        write (*,*) '  '
        write (*,*) ' This code allows for free dynamics or quenching. '
        write (*,*) '  '
        write (*,*) ' Insert iquench -'
        write (*,*) ' iquench =  0 ==> Free dynamics (i.e. pure Newton)'
        write (*,*) ' iquench = -1 ==> dynamical quenching '
        write (*,*) ' iquench = -2 ==> crude constant temperature MD '
        write (*,*) '                  (from velocity rescaling) '
        write (*,*) ' iquench = -3 ==> power quench '
        write (*,*) ' iquench = -4 ==> conjugate gradient minimization '
        write (*,*) ' iquench = -5 ==> Newton-CG minimization (l-bfgs-b)'
        write (*,*) ' iquench = +n ==> periodic quench every n steps '

        write (*,*) '  '
        write (*,*) ' In the event that you are searching for an energy '
        write (*,*) ' minimum, you can also set energy and force tolerances '
        write (*,*) ' (etol and ftot) which will stop execution after '
        write (*,*) ' tolerance is achieved. These are read from the '
        write (*,*) ' the quench.optional file. '
        write (*,*) ' iquench = ', iquench

        if (iquench .lt. -6) then
         write (*,*) ' bad iquench input'
         stop
        end if

! option 8: iensemble
! See Denis J. Evans et al., Phys Rev. A28, 1016 (1983) for velocity rescale
        write (*,*) '  '
        write (*,*) 'You can do various ensembles such as (0) NVE (also the'
        write (*,*) 'correct setting for doing quenching, conjugate '
        write (*,*) 'gradient, etc.), (1) NVT with velocity rescaling, (2) NVT'
        write (*,*) 'with a Nose-Hoover chain thermostat and velocity-verlet'
        write (*,*) 'integrator, (3) NVE with velocity-verlet integrator'
        write (*,*) 'Insert 0/1/2/3.'
        write (*,*) ' iensemble= ', iensemble
        if (iensemble .gt. 0 .and. iquench .ne. 0) then
         write (*,*) ' ******************** NOTE ********************* '
         write (*,*) ' You have chosen to do an NVT ensemble '
         write (*,*) ' but iquench = ', iquench
         write (*,*) ' STOPPING!'
         write (*,*) ' ******************** NOTE ********************* '
         STOP
        end if

! option 9: iBarrier
        write (*,*) '  '
        write (*,*) ' The option exists to calculate a crude energy barrier. '
        write (*,*) ' The additions described here are designed to "push" the '
        write (*,*) ' system of interest from a given initial configuration to '
        write (*,*) ' a given final configuration via a "path of least '
        write (*,*) ' resistance". '
        write (*,*) '  '
        write (*,*) ' At each time step we take the vector difference, for '
        write (*,*) ' each atom, between the current position (at time step #1 '
        write (*,*) ' this is the initial configuration) and the final, or '
        write (*,*) ' desired, configuration.  These difference vectors are '
        write (*,*) ' then transformed into unit vectors.  These unit vectors '
        write (*,*) ' describe the direction that each atom must move in order '
        write (*,*) ' to reach the final configuration. '
        write (*,*) '  '
        write (*,*) ' Next the calculated forces for each atom are examined. '
        write (*,*) ' The question is asked, "Does the force on iatom #i have '
        write (*,*) ' a nonnegative component in the desired direction?" '
        write (*,*) ' (Is force dot unit direction .ge. 0?).  If the answer is '
        write (*,*) ' yes, then we are satisfied that iatom is not being '
        write (*,*) ' pushed away from its desired final position and we do '
        write (*,*) ' nothing.  But, if the answer is no, we add enough force '
        write (*,*) ' IN THE DESIRED DIRECTION to make the component of the '
        write (*,*) ' force in the desired direction equal to a given value. '
        write (*,*) '  '
        write (*,*) ' This procedure insures that iatom always has some '
        write (*,*) ' component of acceleration towards the desired final '
        write (*,*) ' position of iatom.  Given enough time, iatom will arrive '
        write (*,*) ' at the desired final position.  (Two rare cases are '
        write (*,*) ' ignored here, the force on iatom may be exactly zero, '
        write (*,*) ' and the dot product may be exactly zero, but these two '
        write (*,*) ' cases could be handled if necessary.) '

        write (*,*) ' ibarrier = ', ibarrier

! option 10: The big 4
        write (*,*) '  '
        write (*,*) ' In basvec, we will check the big-4 constraints.'
        write (*,*) ' (1) rcm  = 0'
        write (*,*) ' (2) pcm  = 0'
        write (*,*) ' (3) (1/n) (1/2)*sum (m * v**2) = 3/2 kb * T_initial'
        write (*,*) ' (4) ltot = 0'
        write (*,*) '  '
        write (*,*) ' We allow the option of constraining (1)..(4), or '
        write (*,*) ' not constraining.'
        write (*,*) '  '
        write (*,*) ' If you want to constrain quantities, insert 1.'
        write (*,*) ' e.g. for no constriants,  insert 0, 0, 0, 0 '
        write (*,*) ' e.g. for all constraints, insert 1, 1, 1, 1 '
        write (*,*) ' Insert iconstraints(1-4): '

        write (*,101) iconstraints(:)
        write (*,100)

! option 11: ifixcharges
        write (*,*) '  '
        write (*,*) ' Do you want to fix the charge? '
        write (*,*) ' (This is necessary if you are doing a band-structure '
        write (*,*) ' calculation for example). '
        write (*,*) ' ifixcharge = ', ifixcharge

! option 12: ifixneigh
        write (*,*) '  '
        write (*,*) ' Do you want to fix the neighbors? '
        write (*,*) ' (This is necessary if you are doing a phase transition '
        write (*,*) ' type calculation for example). '
        write (*,*) ' ifixneigh = ', ifixneigh

! option 13: icluster
        write (*,*) '  '
        write (*,*) ' Are you doing a cluster (gas-phase molecule)? '
        write (*,*) ' In other words, do you want to set the coulomb '
        write (*,*) ' interaction between cells equal to zero? '
        write (*,*) ' icluster = ', icluster

! option 15: iumbrella
        write (*,*) '  '
        write (*,*) ' Are you using the umbrella sampling algorithm? '
        write (*,*) ' Insert 0/1 = N/Y '
        write (*,*) ' iumbrella = ', iumbrella

! option 16: ivdw
        write (*,*) '  '
        write (*,*) ' The option exists for adding in the long range van der '
        write (*,*) ' Waals interactions (dispersion terms).  This is done '
        write (*,*) ' by adding in an empirical terms proportional to 1/R**6, '
        write (*,*) ' according to a scheme that determines the C_6 parameters '
        write (*,*) ' by using a combination rule. '
        write (*,*) ' Are you including the van der Waals interactions? '
        write (*,*) ' Insert 0/1 = N/Y '
        write (*,*) ' ivdw = ', ivdw

! option 17: igauss
        write (*,*) '  '
        write (*,*) ' Do you want to use gaussian expansions to evaluate bcna '
        write (*,*) ' three-center integrals? '
        write (*,*) ' Note: you must first run gauss_create with the same '
        write (*,*) ' input files you used when you ran create!!! '
        write (*,*) ' igauss = ', igauss
        if (igauss .eq. 1 .and. itheory_xc .ne. 1) then
         write (*,*) ' ********************* WARNING *********************** '
         write (*,*) ' You have chosen to use the gaussian approach for the  '
         write (*,*) ' three-center interactions.  However, the three-center '
         write (*,*) ' exchange-correlation interactions must be evaluated   '
         write (*,*) ' in the old Sankey-Niklewski format (i.e. average      '
         write (*,*) ' densities).                                           '
         write (*,*) ' ********************* WARNING *********************** '
        end if

! option 19: iimage
        write (*,*) '  '
        write (*,*) ' Do you want to reimage the atoms to the central cell at'
        write (*,*) ' every N time steps? (zero means never).  This is needed'
        write (*,*) ' if the atoms move a long distance and you get a bad '
        write (*,*) ' neighbor map as a result.  Note: the coordinates '
        write (*,*) ' printed out will still not be the minimal image.  This '
        write (*,*) ' might be incompatible with a barrier or umbrella '
        write (*,*) ' calculation.  This option is also automatically set '
        write (*,*) ' to zero if you are doing a cluster calculation. '
        if (icluster .ne. 0) iimage=0
        write (*,*) ' iimage = ', iimage
        if (iimage .lt. 0) then
         write (*,*) ' I am afraid that we are going to have to stop. '
         stop
        else if (iimage .gt. 0) then
         write (*,*) ' Will image every ', iimage, ' time steps '
         if (ibarrier .eq.  1 .or.  iumbrella .eq.  1) then
          write (*,*) ' ********************** WARNING ************************'
          write (*,*) ' You have chosen to image the system, but you are doing'
          write (*,*) ' things that will be incompatable.  Check your results.'
          write (*,*) ' ********************** WARNING ************************'
         end if
        end if


! option 21: idynmat
        write (*,*) '  '
        write (*,*) ' You have the option of calculating the dynamical'
        write (*,*) ' matrix. You do this by performing ndim X natoms'
        write (*,*) ' calculations. We displace each atom u0 in'
        write (*,*) ' each of the ndim directions. '
        write (*,*) ' (Normally ndim=3, but may be less for 2d or '
        write (*,*) ' 1d molecules.) Each of these calculations'
        write (*,*) ' is performed in turn by the time step loop. So you'
        write (*,*) ' ultimately need to perform itime=1,ndim*natoms'
        write (*,*) ' time steps.'
        write (*,*) ' Do you wish to perform dynamical matrix calculation?'
        write (*,*) ' idynmat = ', idynmat

! e-ph coupling ?
        write (*,*) 'Do you want to calculate e-ph coupling?'
        write (*,*) ' (0 .. no; 1 .. yes)'
        write (*,*) ' iephc = ', iephc

! option 23: iharmonic
        write (*,*) '  '
        write (*,*) ' Do you want to apply an external field?  If so and'
        write (*,*) ' the field is a harmonic osciallator, put 1. Otherwise'
        write (*,*) ' put zero.'
        write (*,*) ' iharmonic = ', iharmonic

! option 26: iendtemp
        write (*,*) '  '
        write (*,*) ' In a MD simulation you have the option of '
        write (*,*) 'allowing the system to warm up or cool down during'
        write (*,*) 'the course of the simulation.  The change in '
        write (*,*) 'temperature will increment linearly over the '
        write (*,*) 'course of the simulation.  Thus in order to '
        write (*,*) 'have a smooth change in temperature a small '
        write (*,*) 'difference between the initial temp and the final '
        write (*,*) 'temperature should be used over a large number'
        write (*,*) 'of time steps.'
!               This is implemented in readquench.
        write (*,*) ' '
        write (*,*) 'Do you want to change the temperature over the'
        write (*,*) 'course of an MD simulation? (0=no,1=yes)'
        write (*,*) ' iendtemp = ', iendtemp


        if (iendtemp .eq. 1 .and. iensemble .eq. 0) then
            write(*,*) ''
            write(*,*) 'You have selected iendtemp =1 and iensemble=0.'
            write(*,*) 'In order to use iendtemp you must set iensemble=1.'
            write(*,*) ' STOPPING!'
            STOP
        end if

        write (*,*) '  '
        write (*,*) 'Do you want to use Nugged Elastic Band method'
        write (*,*) 'to find out minimum energy path? (0=no,1=yes)'
        write (*,*) ' ineb = ', ineb

        write (*,*) '  '
        write (*,*) 'Do you want to calculate transport?'
        write (*,*) '(0=no,1=yes)'
        write (*,*) ' itrans = ', itrans

        write (*,*) '  '
        write (*,*) 'Do you want to run time-dependent SE ?'
        write (*,*) '(0=no,1=yes)'
        write (*,*) ' itdse = ', itdse

        write (*,*) '  '
        write (*,*) 'Do you want to apply bias voltage ?'
        write (*,*) '(0=no,1=yes)'
        write (*,*) ' ibias = ', ibias

! JOM-add
        write (*,*) '  '
        write (*,*) 'Do you want to run nonadiabatic coupling calculation ?'
        write (*,*) '(0=no,1=yes)'
        write (*,*) ' imdet = ', imdet

        write (*,*) '  '
        write (*,*) 'The number of iteration steps for the time-dependent '
        write (*,*) 'integration for each dynamical time-step dt'
        write (*,*) ' nddt = ', nddt

        write (*,*) '  '
        write (*,*) 'Do you want recalate kpts, lvs and atom positions?'
        write (*,*) '(1.0d0=no)'
        write (*,*) ' rescalar = ', rescal

        write (*,*) '  '
        write (*,*) 'Do you want see the energy in the second line of the answer.xyz?'
        write (*,*) '(0=nothing, 1=energy,temperature 2=energy,temperature,time 3=energy)'
        write (*,*) ' xyz2line = ', xyz2line
! option 27: igap GAP ENRIQUE-FF
        write (*,*) ' You have the option of calculating the molecular'
        write (*,*) ' gap correction using Hartree-Fock or Koopman'
        write (*,*) ' aproximations. If igap = 0 no gap corrections'
        write (*,*) ' will be done. If igap = 1, Hartree-Fock correction'
        write (*,*) ' will be used. If igap = 2, Koopman displacement of'
        write (*,*) ' molecular levels will be calculated. If igap = 3'
        write (*,*) ' scissor operator (maybe using Koopman correction) will'
        write (*,*) ' be used.'
        write (*,*) ' igap = ', igap
! end GAP ENRIQUE-FF

! SECTION SCF
        if (iharris .ne. 1) then
        ! option 4: max_scf_iterations
         write (*,*) '  '
         write (*,*) ' The maximum number of iterations to self- '
         write (*,*) ' consistency.'
         write (*,*) ' max_scf_iterations = ', max_scf_iterations
         if (max_scf_iterations .lt. 1) then
          write (*,*) ' Variable max_scf_iterations set < 1. Set '
          write (*,*) ' max_scf_iterations > 1. Try again with larger value. '
          stop
         end if
         if (max_scf_iterations .gt. 400) then
          write (*,*) ' Variable max_scf_iterations set > 200. Set '
          write (*,*) ' max_scf_iterations < 200. Try again with smaller value. '
          stop
         end if
         write (*,*) '  '
         write (*,*) ' The mixing algorithm for electronic structure,'
         write (*,*) ' ialgmix = ', ialgmix
         write (*,*) ' 1 ... Anderson method ' 
         write (*,*) '   ( see V. Eyert, J. Comp. Phys. 124, 271 (1996))'
         write (*,*) ' 2 ... Broyden method '
         write (*,*) '   ( see D.Vanderbilt & S.G.Louie, Phys.Rev.B 30, 6118 (1984))'
         write (*,*) ' 3 ... Louie method '
         write (*,*) '   (see D.Vanderbilt & S.G.Louie, Phys.Rev.B 30, 6118 (1984))'
         write (*,*) ' 4 ... Pulay (RMM-DIIS) method ' 
         write (*,*) '   (see F. Eckert et al, J. Comp. Chem. 18, 1473 (1997))'
         write (*,100)
         write (*,*) '  '
         write (*,*) ' The convergence criteria for electronic structure,'
         write (*,*) ' sigmatol = ', sigmatol
         write (*,100)
         write (*,*) '  '
         write (*,*) ' the mixing factor in SCF iterations, '
         write (*,*) ' bmix = ',bmix
         write (*,100)
         write (*,*) ' The Fermi temperature. This variable determine'
         write (*,*) ' smearing of (un)occupied states (1 eV = 11604 K)'
         write (*,*) ' tempfe = ', tempfe
         write (*,100)
        endif

! QM/MM
        if (iqmmm .eq. 1) then
         write (*,100)
         write (*,*) ' iqmmm = 1'
         write (*,*) 'QM/MM with electrostatic embedding'
         write (*,100)
        end if
!JIMM
! XYZ Dipole
        if (idipole .eq. 1) then
         write (*,100)
         write (*,*) ' idipole = 1'
         write (*,*) 'Long range interactions with XYZ dipole'
         write (*,100)
        end if
	
        if (idipole .eq. 1 .and. icluster .eq. 0) then
         write (*,*) ' ******************** NOTE ********************* '
         write (*,*) ' idipole = 1 theory is not compatible with periodic systems '
         write (*,*) ' STOPPING!'
         write (*,*) ' ******************** NOTE ********************* '
         STOP
        end if
!JIMM
!DFTD3
        if (idftd3 .eq. 1) then
         write (*,100)
         write (*,*) ' idftd3 = 1', 'dftd3_func = ', dftd3_func, 'dftd3_version = ', dftd3_version, 'tz =', dftd3_tz
         write (*,*) ' DFTD3 corrections, J. Chem. Phys. 132, 154104 (2010)'
         write (*,100)
        end if

        if (idftd3 .eq. 2) then
         write (*,100)
	 write (*,*) ' choose your own parameters for dftd3 (van der Waals) '
         write (*,*) ' idftd3 = 2 ', 'dftd3_s6 =', dftd3_s6,'dftd3_rs6 =', dftd3_rs6,'dftd3_s18 =', dftd3_s18, &
                     'dftd3_rs18 =', dftd3_rs18,'dftd3_alp =',dftd3_alp
         write (*,*) ' DFTD3 corrections, J. Chem. Phys. 132, 154104 (2010)'
         dftd3_params(1)=dftd3_s6
         dftd3_params(2)=dftd3_rs6
         dftd3_params(3)=dftd3_s18
         dftd3_params(4)=dftd3_rs18
         dftd3_params(5)=dftd3_alp
        end if


! SECTION OUTPUTS
        if (iwrtcdcoefs .gt. 0)                                              &
     &   write (*,*) ' Writing out wavefunction coefficients. '
        if (iwrtcharges .gt. 0)                                              &
     &   write (*,*) ' Writing out charges (Lowdin or Mulliken). '
        if (iwrtdensity .gt. 0) write (*,*) ' Writing out density matrix. '
        if (iwrteigen .gt. 0) write (*,*) ' Writing out eigenvalues. '
        if (iwrtefermi .gt. 0) write (*,*) ' Writing out fermi occupations. '
        if (iwrtfpieces .gt. 0) write (*,*) ' Writing out force components. '
        if (iwrthampiece .gt. 0)                                             &
     &   write (*,*) ' Writing out Hamiltonian matrix elements. '
        if (iwrtcomponents .gt. 0)                                           &
     &   write (*,*) ' Writing out components of band-structure energy. '
        if (iwrtneigh .gt. 0) write (*,*) ' Writing out neighbor map. '
        if (iwrtneigh_com .gt. 0)                                            &
     &   write (*,*) ' Writing out common neighbor map. '
        if (iwrtxyz .gt. 0) write (*,*) ' Writing out xyz file - answer.xyz '
        if (iwrtdos .gt. 0) write (*,*) ' Writing out dos files '
        if (iwrthop .gt. 0) write (*,*) ' Writing out hopping values for STM '
        if (iwrtatom.gt. 0) write (*,*) ' Writing out the Atomo_i files '
        if (iwrtpop .gt. 0) write (*,*) ' Writing out population file '
        if (iwrtHS .gt. 0) write (*,*) ' Writing out H & S file '
        if (iwrtvel .gt. 0) write (*,*) ' Writing out VELOCITY.dat file '
        if (iwrtden .gt. 0) write (*,*) ' Writing out density projected on the grid '
        if (iwrtewf .gt. 0) write (*,*) ' Writing out eigenfunctions projected on the grid '
        if (iwrtxsf .gt. 0) write (*,*) ' Writing out xsf-format file  '
        if (idensimport .gt. 0) write (*,*) ' Importing density file for projection  '


! writeout resume of the input variables into param.dat file
        open (unit = 50, file = 'param.dat', status = 'unknown')
        write (50, *) ' FILES:'
        write (50, *) '  basisfile         : ',basisfile
        write (50, *) '  lvsfile           : ',lvsfile
        write (50, *) '  kptsfile          : ',kptpreference
        write (50, *) '  acfile            : ',acfile
        write (50, *) '  xvfile            : ',xvfile
        write (50,100)
        write (50, *) ''
        write (50, *) ' TIME STEP:'
        write (50, *) '  nstepi            : ',nstepi
        write (50, *) '  nstepf            : ',nstepf
        write (50, *) '  dtime             : ',dt
        write (50,100)
        write (50, *) ''
        write (50, *) ' OPTIONS:'
        write (50, *) '  iharris           : ',iharris
        write (50, *) '  idogs             : ',idogs
        write (50, *) '  ihubbard          : ',ihubbard
        write (50, *) '  ihorsfield        : ',ihorsfield
        write (50, *) '  imcweda           : ',imcweda
        write (50, *) '  igsn              : ',igsn
        write (50, *) '  iks               : ',iks
        write (50, *) '  iqout             : ',iqout
        write (50, *) '  qstate            : ',qstate
        write (50, *) '  icluster          : ',icluster
        write (50, *) '  iensemble         : ',iensemble
        write (50, *) '  ifixcharge        : ',ifixcharge
        write (50, *) '  ifixneigh         : ',ifixneigh
        write (50, *) '  iumbrella         : ',iumbrella
        write (50, *) '  ibarrier          : ',ibarrier
        write (50,300) iconstraints(:)
        write (50, *) '  iharmonic         : ',iharmonic
        write (50, *) '  iimage            : ',iimage
        write (50, *) '  idynmat           : ',idynmat
        write (50, *) '  iephc             : ',iephc
        write (50, *) '  ivdw              : ',ivdw
        write (50, *) '  ineb              : ',ineb
        write (50, *) '  itrans            : ',itrans
        write (50, *) '  itdse             : ',itdse
        write (50, *) '  imdet             : ',imdet
        write (50, *) '  nddt              : ',nddt
        write (50, *) '  ibias             : ',ibias
        write (50, *) '  rescalar          : ',rescal
        write (50, *) '  icDFT             : ',icdft
        write (50, *) '  iqmmm             : ',iqmmm
        write (50, *) '  mix_embedding     : ',mix_embedding
        write (50, *) '  cut_embedding     : ',cut_embedding
        write (50, *) '  idipole           ; ',idipole
        write (50, *) '  idftd3            : ',idftd3
        write (50, *) '  dftd3_func        : ',dftd3_func
        write (50, *) '  dftd3_version     : ',dftd3_version
        write (50, *) '  dftd3_tz          : ',dftd3_tz
        write (50,100)
        write (50, *) ''
        write (50, *) ' SCF'
        write (50, *) '  max_scf_iterations  : ',max_scf_iterations
        write (50, *) '  ialgmix             : ',ialgmix
        write (50, *) '  bmix                : ',bmix
        write (50, *) '  tempfe              : ',tempfe
        write (50, *) '  sigmatol            : ',sigmatol
        write (50, *) ''
        write (50, *) ' TEMPERATURE'
        write (50, *) '  T_initial           : ',T_initial
        write (50, *) '  T_final             : ',T_final

        close (50)

        write (*,100)

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, 4i4)
200     format (a30)
201     format (2x, ' file = ', a30)
300     format ('   iconstraints      : ',4i4)

        return
      end subroutine readparam
