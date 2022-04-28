# copyright info:
#
#                             @Copyright 2006
#                           Fireball Committee
# West Virginia University - James P. Lewis, Chair
# Arizona State University - Otto F. Sankey
# Universidad de Madrid - Jose Ortega
# Academy of Sciences of the Czech Republic - Pavel Jelinek

# Other contributors, past and present:
# Auburn University - Jian Jun Dong
# Arizona State University - Gary B. Adams
# Arizona State University - Kevin Schmidt
# Arizona State University - John Tomfohr
# Lawrence Livermore National Laboratory - Kurt Glaesemann
# Motorola, Physical Sciences Research Labs - Alex Demkov
# Motorola, Physical Sciences Research Labs - Jun Wang
# Ohio University - Dave Drabold
# University of Regensburg - Juergen Fritsch

#
# RESTRICTED RIGHTS LEGEND
# Use, duplication, or disclosure of this software and its documentation
# by the Government is subject to restrictions as set forth in subdivision
# { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
# clause at 52.227-7013.

# Do you optimize or do you debug? (OPT/DEBUG)
MODE = OPT 
#MODE = DEBUG
# Do you use vendor optimized blas and lapack libraries? (YES/NO)
USEBLAS = YES 
#USEBLAS = NO 
# What machine? (AIX/AIX_LLNL/BLUEHORIZON/LINUX/LAPTOP/MARYLOU/MARYLOU10/
# SGI_LANL/SIMPLE/ICEBOX/ALPHA/TRU64/ILINUX/GFORTRAN)
#MACHINE = ILINUX_MAC
#MACHINE = ILINUX_fast05
#MACHINE = ILINUX_beem
#MACHINE = ILINUX_fast05
MACHINE = ILINUX_fast05.static
MACHINE = ILINUX_gnu
#MACHINE = SGI-krejci
#MACHINE = SGI
#MACHINE = GFORTRAN

# Pick your parallel type: <blank for none>, MPI, MPICH (MUST set method below)
# or OPENMP 
#PARALLEL = MPI
#PARALLEL = OPENMP
PARALLEL = 
# Note: the linear-scaling option is in the testing stages only right now, so
# setting to ORDERN here may mess things up.
# Pick the method that you are using: <blank for nothing special>, SCALAPACK
# (parallel diagonalizer) or ORDERN (linear-scaling BTN minimizer). 
# Your version of BLACS must use the same version of MPI!  You can get BLACS
# from http://www.netlib.org/blacs/archives/ precompiled or from 
# http://www.netlib.org/blacs/BLACS/Papers.html as source code.  You can get 
# SCALAPACK from http://www.netlib.org/scalapack/ precompiled or as source code.
# other options of METHOD
# DOUBLE ... double precision
# LAPACK95 ... lapack_95 
# GAMMA ... use real artithemtics for gamma k-point calc (large systems)  
METHOD = DOUBLE 
#METHOD = GAMMA 
#METHOD = SCALAPACK
# Do you want to do thermodynamic integration (NO,YES)?  If so 
# you'll need a compatible C compiler.
CCOMPILE = NO
# Do you use TclMD option? (YES/NO)
# Note: you must have TclMD library installed and Tcl/Tk to use this.
# You can get Tcl/Tk from http://dev.scriptics.com/software/tcltk/
TCLMD = NO

# Choose a profiler - none, VAMPIR
PROFILER = -pg

# Set some defaults for non-MPI
ORDERN = 
FORM_RHO = denmat.o denmat_es.o denmat_KS.o fermie.o denmata_ordern_fk.o \
	denmatb_ordern_fk.o denmatc_ordern_fk.o ss12_fk.o build_rho.o \
	build_rho_KS.o build_Ji.o koopman.o build_nij.o project_eh.o project_wfmdet.o
KSPACE = kspace.o kspace_KS.o kspace_ordern_fk.o diag_error.o diag_k.o \
	diag_k_KS.o
INITMPI = init_noMPI.o
UTIL_SPARSE = 
UTIL = anderson.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o  writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o

	#bcast_k-MPI.o # this last one has to be added in order to at least call the dummy routine

DASSEMBLERS = Dassemble_2c.o Dassemble_3c.o Dassemble_ca_2c.o \
	Dassemble_ca_3c.o Dassemble_eh_2c.o Dassemble_hxc_2c.o \
	Dassemble_hxc_3c.o Dassemble_lr.o  Dassemble_snxc_on.o \
	Dassemble_olsxc_on.o Dassemble_olsxc_2c.o Dassemble_olsxc_3c.o \
	Dassemble_snxc_2c.o Dassemble_snxc_3c.o Dassemble_2c_PP.o \
	Dassemble_3c_PP.o Dassemble_ca_olsxc_on.o Dassemble_ca_snxc_on.o \
	Dassemble_ca_snxc_3c.o Dassemble_ca_olsxc_3c.o \
	Dassemble_ca_snxc_2c.o Dassemble_ca_olsxc_2c.o Dassemble_qmmm.o Dassemble_qmmm_dip.o \
	Dassemble_ca_2c_dip.o Dassemble_ca_3c_dip.o Dassemble_lr_dip.o \
	getforces_mcweda.o getforces_eh.o getforces_hxc.o \
	getforces_KS.o getforces_classic.o getforces_classic_RGL.o getforces_classic_vdw.o \
	getforces.o getforces_socket.o getforces_classic_tersoff.o \
	getforces_zw.o Dassemble_zw_2c_ct.o \
	Dassemble_zw_3c_ct.o Dassemble_zw_2c_na.o Dassemble_zw_3c_na.o \
        Dassemble_zw_on_na.o

INTERACTIONS = cl_value.o Dtrescentros.o doscentros.o doscentrosPP.o  \
	doscentrosDipY.o doscentrosDipX.o  \
	get_ewald.o get_vdw.o getHarmonic.o trescentros.o unocentros.o \
	smoother.o  doscentrosS.o trescentrosS.o DtrescentrosS.o \
	dosgaussians.o \
	gelements_VXC.o Dgelements_VXC.o \
	gelementsG_VXC.o  DgelementsG_VXC.o \
	gelementsG_VNA.o  DgelementsG_VNA.o \
	gelementsGS_VXC.o DgelementsGS_VXC.o \
	gelementsG_VNA_SH.o DgelementsG_VNA_SH.o \
	gelementsG_VNA_SH_RNA.o DgelementsG_VNA_SH_RNA.o \
	trescentrosGHXC_VXC.o  DtrescentrosGHXC_VXC.o \
	trescentrosG_VNA.o  DtrescentrosG_VNA.o \
	trescentrosG_VXC.o  DtrescentrosG_VXC.o \
	trescentrosGS_VXC.o DtrescentrosGS_VXC.o \
	trescentrosG_VNA_SH.o DtrescentrosG_VNA_SH.o \
	doscentrosGS_overlap.o gelementsGS_overlap.o DgelementsGS_overlap.o \
	doscentrosG_overlap.o gelementsG_overlap.o DgelementsG_overlap.o \
	internalLambda.o tester2c.o

# Set some defaults for non-TclMD stuff 
VISUALIZATION = noTclMD.o nograce.o noxmgr.o
# SCALAPACK 
ifneq (,$(findstring SCALAPACK, ${METHOD}))
	KSPACE = blacsaba.o diag_error.o kspace_ordern_fk.o kspace_MPI.o \
		kspace_MPI_slave.o pclagetter.o pclaputter.o diag_k.o \
		diag_k_KS.o kspace_KS.o
	INITMPI = mpi_declarations.o init_MPI.o

	UTIL = anderson2.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o  compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o
endif

# OPENMP 
#ifneq (,$(findstring OPENMP, ${PARALLEL}))

#DASSEMBLERS = Dassemble_2c.o Dassemble_3c.o Dassemble_ca_2c.o \
	Dassemble_ca_3c.o Dassemble_eh_2c.o Dassemble_hxc_2c.o \
	Dassemble_hxc_3c.o Dassemble_lr_OMP.o  Dassemble_snxc_on.o \
	Dassemble_olsxc_on.o Dassemble_olsxc_2c.o Dassemble_olsxc_3c.o \
	Dassemble_snxc_2c.o Dassemble_snxc_3c.o Dassemble_2c_PP.o \
	Dassemble_3c_PP.o Dassemble_ca_olsxc_on.o Dassemble_ca_snxc_on.o \
	Dassemble_ca_snxc_3c.o Dassemble_ca_olsxc_3c.o \
	Dassemble_ca_snxc_2c.o Dassemble_ca_olsxc_2c.o Dassemble_qmmm.o Dassemble_qmmm_dip.o\
	getforces_mcweda.o getforces_eh.o getforces_hxc.o \
	getforces_KS.o getforces_classic.o getforces_classic_vdw.o getforces_classic_RGL.o getforces_classic_tersoff.o

#INTERACTIONS = cl_value.o Dtrescentros.o doscentros.o doscentrosPP.o  \
	get_ewald_OMP.o get_vdw.o getHarmonic.o trescentros.o unocentros.o \
	smoother.o  doscentrosS.o trescentrosS.o DtrescentrosS.o \
	dosgaussians.o \
	gelements_VXC.o Dgelements_VXC.o \
	gelementsG_VXC.o  DgelementsG_VXC.o \
	gelementsG_VNA.o  DgelementsG_VNA.o \
	gelementsGS_VXC.o DgelementsGS_VXC.o \
	gelementsG_VNA_SH.o DgelementsG_VNA_SH.o \
	gelementsG_VNA_SH_RNA.o DgelementsG_VNA_SH_RNA.o \
	trescentrosGHXC_VXC.o  DtrescentrosGHXC_VXC.o \
	trescentrosG_VNA.o  DtrescentrosG_VNA.o \
	trescentrosG_VXC.o  DtrescentrosG_VXC.o \
	trescentrosGS_VXC.o DtrescentrosGS_VXC.o \
	trescentrosG_VNA_SH.o DtrescentrosG_VNA_SH.o \
	doscentrosGS_overlap.o gelementsGS_overlap.o DgelementsGS_overlap.o \
	doscentrosG_overlap.o gelementsG_overlap.o DgelementsG_overlap.o \
	internalLambda.o tester2c.o

#endif
# ORDERN
ifneq (,$(findstring ORDERN, ${METHOD}))
	ORDERN = ordern.o allocate_ordern.o assemble_2c_ordern_final.o \
		assemble_2c_ordern_init.o assemble_3c_ordern_final.o \
		assemble_ca_2c_ordern_final.o assemble_ca_3c_ordern_final.o \
		assemble_eh_2c_ordern_final.o \
		Dassemble_2c_ordern_final.o \
		Dassemble_3c_ordern_final.o \
		Dassemble_ca_2c_ordern_final.o Dassemble_ca_3c_ordern_final.o \
	FORM_RHO = chebft.o denmata_ordern.o denmatb_ordern.o denmatc_ordern.o \
		denmat_fk.o formrho_sparse.o ss12.o build_rho.o build_rho_KS.o
	KSPACE = eandg.o xeandg.o formc_compact.o formsh_compact.o initguess.o \
		kspace_fk.o kspace_ordern.o kspace_ordern_init.o \
		kspace_ordern_slave.o ordern_init.o qralg.o set_dimensions.o \
		set_maxdimension.o getsendrecv.o getstepsize.o
	UTIL_SPARSE = build_transpose.o lanc.o sparse_add.o sparse_copy.o \
		sparse_getdimension.o sparse_getpacksize.o sparse_mask.o \
		sparse_mult.o sparse_norm2.o sparse_pack.o \
		sparse_pack_elements.o sparse_pack_indices.o sparse_unpack.o \
		sparse_unpack_elements.o sparse_unpack_indices.o \
		sparse_vecprod.o 
	INITMPI = mpi_declarations.o init_MPI.o
endif
# LAPACK_95
ifneq (,$(findstring LAPACK95, ${METHOD}))
	KSPACE = kspace_l95.o kspace_ordern_fk.o diag_error.o diag_k.o \
	diag_k_KS.o
	UTIL = anderson_l95.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o  compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o
endif

# DOUBLE
ifneq (,$(findstring DOUBLE, ${METHOD}))
	KSPACE = kspace2.o kspace_KS.o kspace_ordern_fk.o diag_error.o diag_k.o \
	diag_k_KS.o
	UTIL = anderson2.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o
endif

# REAL ARITHMETICS (only gamma k-point)
ifneq (,$(findstring GAMMA, ${METHOD}))
	KSPACE = kspace_withnok.o kspace_KS.o kspace_ordern_fk.o diag_error.o diag_k.o \
	diag_k_KS.o
	UTIL = anderson2.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o
endif

# MPI-k
ifneq (,$(findstring MPI-k, ${PARALLEL}))
	KSPACE = kspace_MPI-k.o kspace_ordern_fk.o diag_error.o diag_k-MPI.o
	UTIL = anderson2.o fixfrags.o fixfrags2.o hampiece.o push_atoms.o \
	writeout_ac.o writeout_cd.o writeout_charges.o writeout_dipole.o writeout_comph.o \
	writeout_neighbors.o writeout_xv.o writeout_neighborsPP.o \
	hamtrans.o mixer.o diag_k-MPI_slave.o bcast_k-MPI.o bcast_k-MPI_slave.o \
	den2mesh.o ew2mesh.o ew2mesh_gamma.o postscf.o den2mesh_import.o \
	broyden.o louie.o pulay.o compute_neutral_rho_grid.o \
	ew2mesh_fourier.o ew2mesh_kscan.o ew2mesh_ARPES.o \
	project_orb.o project_orb_complex.o excitations.o kvaziband.o \
	band_project.o ARPES_LCAO.o writeout_eigenvec.o writeout_eigenvec_G.o \
	writeCoefsLCAO.o compute_all_charges.o
	INITMPI = mpi_declarations.o init_MPI.o
endif

# Note: These are the machines files used by the Lewis Research Group. You need 
# to create a MACHINES file specific to your architecture.  You can use one of
# the machines file below as a starting point.  
include MACHINES/${MACHINE}
#include MACHINES/ERROR
#include MACHINES/AIX
#include MACHINES/SGI
#include MACHINES/BIGBEN
#include MACHINES/GFORTRAN
#include MACHINES/ILINUX
#include MACHINES/ILINUX_CLS1
#include MACHINES/ILINUX_MAC
#include MACHINES/ILINUX2
#include MACHINES/MARYLOU
#include MACHINES/MARYLOU4
#include MACHINES/MARYLOUX
#include MACHINES/PLINUX
#include MACHINES/RACHEL
#include MACHINES/SIMPLE
#include MACHINES/ILINUX_PARALLEL_CLS1

ALLOCATIONS = allocate_f.o allocate_h.o allocate_neigh.o allocate_rho.o \
	allocate_umb.o allocate_steered.o reallocate_f.o reallocate_h.o reallocate_neigh.o \
	reallocate_rho.o allocate_dos.o allocate_grid.o allocate_trans.o 

ASSEMBLERS = assemble_olsxc_1c.o assemble_hxc_1c.o assemble_2c.o assemble_3c.o \
	assemble_ca_2c.o assemble_3c_PP.o assemble_2c_PP.o\
	assemble_ca_3c.o assemble_eh_2c.o assemble_eh_usr.o assemble_F.o \
	assemble_hxc_2c.o assemble_hxc_3c.o assemble_lr.o assemble_sVNL.o \
	assemble_usr.o buildh.o assemble_olsxc_on.o assemble_olsxc_off.o \
	build_olsxc_on.o build_olsxc_off.o average_rho.o average_ca_rho.o\
	build_snxc_on.o build_snxc_off.o assemble_snxc_on.o \
	assemble_snxc_off.o build_ca_snxc_on.o build_ca_olsxc_on.o \
	assemble_h.o assemble_mcweda.o assemble_hxc.o assemble_eh.o \
	getenergy.o getenergy_hxc.o getenergy_mcweda.o getenergy_eh.o \
	assemble_h_ks.o getenergy_KS.o assemble_S.o assemble_2c_S.o \
	assemble_hartree.o assemble_scissor.o assemble_qmmm.o assemble_qmmm_dip.o\
	assemble_ca_2c_dip.o assemble_ca_3c_dip.o assemble_lr_dip.o \
        assemble_zw_1c_na.o assemble_zw_2c_ct.o assemble_zw_3c_ct.o \
        assemble_xczw.o assemble_zw_off_na.o assemble_zw_on_na.o \
        build_zw_off_na.o build_zw_on_na.o assemble_1c_vdip.o \
        getenergy_zw.o


GRID = assemble_KS_den0.o assemble_KS_den.o assemble_KS_usr.o laplace_fft.o \
	assemble_KS_dcc.o assemble_KS_mat.o mixer_KS.o writeout_charges_KS.o \
	writeout_xsf.o assemble_KS_vna.o get_Hort.o
#	writeout_xsf.o psi2mesh.o psi22mesh.o assemble_KS_vna.o

INITIALIZERS = diagnostics.o initatomicE.o initconstraints.o initcharges.o \
	initconstants.o initboxes.o initkpoints.o initmasses.o initneighbors.o \
	welcome.o make_mu2shell.o make_munu.o make_munuPP.o restart.o \
        make_munuDipY.o make_munuDipX.o \
	zero_ang_mom.o initamat.o make_munuS.o initNH.o getkpoints.o \
	greatKAuto.o greatKsubsAuto.o initgrid.o initdenmat.o \
	get_info_orbital.o initbasics.o initcharges_KS.o initcDFT.o 

INTERPOLATERS = buildspline_1d.o interpolate_1d.o interpolate_2d.o \
	recover_2c.o recover_3c.o recover_PP.o recoverC.o setterp_2d.o \
	recover_2cDipY.o recover_2cDipX.o \
	recover_S.o buildspline2_1d.o getpsi.o getYlm.o getvna.o

LOOPS = main_loop.o main_loop_MD.o main_loop_CG.o scf_loop.o scf_loop_harris.o \
	main_loop_NEB.o main_loop_DM.o scf_loop_ks.o main_loop_importrho.o \
	main_loop_TDSE.o main_loop_MDET.o main_loop_MIN.o main_loop_NAC.o \
	main_loop_FIRE.o main_loop_socket.o

MAIN = fireball.o 

MAIN_SERVER = fireball_server.o

MAIN_SERVER_AMBER = fireball_server_amber.o

MD = cross.o factorial.o predictor.o gaussT.o corrector.o imaged.o setgear.o \
	phimat.o bvec.o soldm.o NHCThermostat.o get_enk.o writeHNose.o resetNHC.o \
	move_ions.o

MODULES = barrier.o charges.o configuration.o constants_fireball.o density.o \
	dimensions.o forces.o fragments.o gaussG.o integrals.o interactions.o \
	kpoints.o neighbor_map.o umbrella.o  steered.o optimization.o module_dos.o \
	dynamo.o cproc.o noseHoover.o scf.o grid.o wavefunction.o neb_module.o \
	vnneutral.o transport.o matmultmod.o outputs.o options.o energy.o \
	MD.o  mpi_main.o tdse.o bias.o nonadiabatic.o hartree.o sockets.o fsockets.o fb_socket.o

MODULES_C =  $(MODULES) classicMD.o

NEIGHBORS = backnay.o common_neighbors.o find_neigh_max.o find_neigh_max_class.o \
	mpairnay.o neighbors.o neighbors_pairs.o find_neighPP_max.o neighborsPP.o \
	common_neighborsPP.o num_neigh_tot.o

PRESSURE = hmetric.o initpressure.o invert3x3.o

READFILES = append_string.o read_1c.o read_2c.o read_3c.o readbasis.o \
	readbarrier.o readdata_2c.o readdata_3c.o readfragments.o \
	readheader_2c.o readheader_3c.o readinfo.o readlvs.o \
	readparam.o readphi.o readpressure.o readquench.o \
	readsa.o readvdw.o readcgo.o readdos.o \
	readgaussG.o findFdata.o readgrid.o read_wf.o read_vna.o \
	readtrans.o readbind.o readhop.o readdata.o readdata_hxc.o \
	readdata_mcweda.o readdata_xczw.o readdata_eh.o checksum_options.o readdata_KS.o \
	getsections.o readdata_classicMD.o readhartree.o readephc.o

ROTATIONS = chooser.o chooserd.o deps2center.o deps3center.o makeDmat.o \
	makeDmatPP.o rotate.o rotated.o rotatedPP.o twister.o twisterd.o \
	rotatePP.o epsilon.o

ifneq (,$(findstring YES, ${CCOMPILE}))
THERMOINT = cclient.o 
else
THERMOINT =
endif

SOCKETS = get_geometry.o create_socket.o send_geometry.o sendrecv.o soc_init.o

SOLVESH_DIAG = $(KSPACE) $(BLAS)

UMBRELLA = assemble_umbrella.o Dassemble_umbrella.o get_umbrella.o \
	readumbrella.o assemble_steered.o Dassemble_steered.o \
        get_steered.o readsteered.o

XC = ceperley_alder.o cepal.o

CG = cgo.o bfgs.o l-bfgs-b.o FIRE.o

DOS = dos.o invierte.o writeout_dos.o writeout_dosng.o hoppings.o writeout_atom.o \
	hamilt_atom.o 

NEB = initneb.o neb.o

TDSE = ete_loop.o psi2es.o eigenHS.o tddenmat.o tddiag_k.o \
	diag_Sk.o diag_Hk.o allocate_tdse.o tdbc.o propTpsi.o \
	readtdse.o initpsi.o ortho_H.o get_QLow.o get_QMul.o \
	postete.o wrtout_psiT.o

TRANS = assemble_t12_fit.o assemble_t12_bare.o  calcG.o  assemble_Hsam.o \
	assemble_Gsam.o assemble_Dxx.o sqrt_mat.o interpolate_hop.o gethop.o
	
BIAS = assemble_bias.o Dassemble_bias.o allocate_bias.o reallocate_bias.o \
	readbias.o

NAC = allocate_nac.o assemble_G_S.o nacouplings.o build_gover1c.o init_mdet.o \
      mdetdenmat.o getforces_mdet.o save_mdetstuff.o evolve_ks_states.o \
      deallocate_nac.o delta_t_ks.o dcdt_nac.o Dassemble_2c_mdet.o \
      Dassemble_2c_PP_mdet.o Dassemble_olsxc_on_mdet.o Dassemble_olsxc_2c_mdet.o \
      Dassemble_3c_mdet.o Dassemble_3c_PP_mdet.o Dassemble_olsxc_3c_mdet.o \
      fewest_switches.o mc_switch.o transition.o \
      Dassemble_ca_2c_mdet.o Dassemble_ca_3c_mdet.o \
      Dassemble_lr_mdet.o \
      Dassemble_ca_olsxc_on_mdet.o Dassemble_ca_olsxc_2c_mdet.o \
      Dassemble_ca_olsxc_3c_mdet.o move_correc.o move_predic.o overlap_sign.o \
      check_swap.o overlap_numeric.o getnac.o MCsolar.o Dassemble_qmmm_mdet.o \
      Dassemble_ca_2c_mdet_dip.o Dassemble_qmmm_mdet_dip.o \
      Dassemble_ca_3c_mdet_dip.o Dassemble_lr_mdet_dip.o vibcouplings.o

QMMM =  main_loop_MDET_qmmm.o main_loop_MD_qmmm.o fireball_qmmm_loop.o


DFTD3 = common.o sizes.o pars.o core.o api.o dftd3_corrections.o

#NAC = allocate_nac.o assemble_G_S.o nacouplings.o build_gover1c.o init_mdet.o \
#      mdetdenmat.o getforces_mdet.o save_mdetstuff.o evolve_ks_states.o \
#      deallocate_nac.o dcdt_nac.o Dassemble_2c_mdet.o \
#      Dassemble_2c_PP_mdet.o Dassemble_olsxc_on_mdet.o Dassemble_olsxc_2c_mdet.o \
#      Dassemble_3c_mdet.o Dassemble_3c_PP_mdet.o Dassemble_olsxc_3c_mdet.o \
#      fewest_switches.o mc_switch.o transition.o \
#      Dassemble_ca_2c_mdet.o Dassemble_ca_3c_mdet.o Dassemble_lr_mdet.o \
#      Dassemble_ca_olsxc_on_mdet.o Dassemble_ca_olsxc_2c_mdet.o \
#      Dassemble_ca_olsxc_3c_mdet.o move_correc.o move_predic.o


OBJECTS_COM = $(DFTD3) $(INITMPI) $(ORDERN) $(ALLOCATIONS) $(ASSEMBLERS) \
        $(DASSEMBLERS) $(INITIALIZERS) $(INTERACTIONS) $(INTERPOLATERS) \
        $(LOOPS) $(MD) $(NEIGHBORS) $(PRESSURE) $(READFILES) \
        $(ROTATIONS) $(SOLVESH_DIAG) $(FORM_RHO) $(UMBRELLA) $(UTIL) \
        $(UTIL_SPARSE) $(VISUALIZATION) $(XC) $(CG) $(DOS) $(THERMOINT) \
        $(NEB) $(TRANS) $(GRID) $(TDSE) $(BIAS) $(NAC)


OBJECTS = $(MODULES_C) $(OBJECTS_COM) $(MAIN)

OBJECTS_QMMM = $(MODULES_C) $(OBJECTS_COM) $(QMMM)

OBJECTS_SERVER = $(MODULES_C)  \
	$(OBJECTS_COM) $(MAIN_SERVER)

OBJECTS_SERVER_AMBER = $(MODULES_C) $(OBJECTS_COM) $(MAIN_SERVER_AMBER)

fireball.x: $(OBJECTS)
	$(F90) -o  fireball.x $(FFLAGS) $(OBJECTS) $(VISFLAGS) $(PARLFLAGS) \
	$(LFLAGS) 

.PHONY: clean veryclean extraclean

clean:
	rm -f -r core *.o .nfs* rii_files fireball.x.ip*  *.mod ldtmp* \
			 *.vo *~ *.il

veryclean: clean
	rm -f fireball.x libfireball.a

extraclean: veryclean

all:
	make fireball.x

libfireball: $(OBJECTS_QMMM)
	ar rv libfireball.a $(OBJECTS_QMMM)
	ranlib libfireball.a

server: $(OBJECTS_SERVER)
	$(F90)  -o  fireball_server.x $(FFLAGS) $(OBJECTS_SERVER) $(VISFLAGS) $(PARLFLAGS) \
	$(LFLAGS) 
	
server_amber: $(OBJECTS_SERVER_AMBER)
	$(F90)  -o  fireball_server $(FFLAGS) $(OBJECTS_SERVER_AMBER) $(VISFLAGS) $(PARLFLAGS) \
	$(LFLAGS) 

#*****************************************************************************
# modules
# *****************************************************************************
dimensions.o : MODULES/dimensions.f90
	$(F90) $(FFLAGS) -c MODULES/dimensions.f90
barrier.o : MODULES/barrier.f90
	$(F90) $(FFLAGS) -c MODULES/barrier.f90
charges.o : MODULES/charges.f90 dimensions.o
	$(F90) $(FFLAGS) -c MODULES/charges.f90
configuration.o : MODULES/configuration.f90
	$(F90) $(FFLAGS) -c MODULES/configuration.f90
constants_fireball.o : MODULES/constants_fireball.f90
	$(F90) $(FFLAGS) -c MODULES/constants_fireball.f90
density.o : MODULES/density.f90
	$(F90) $(FFLAGS) -c MODULES/density.f90
forces.o : MODULES/forces.f90
	$(F90) $(FFLAGS) -c MODULES/forces.f90
fragments.o : MODULES/fragments.f90
	$(F90) $(FFLAGS) -c MODULES/fragments.f90
gaussG.o : MODULES/gaussG.f90 dimensions.o
	$(F90) $(FFLAGS) -c MODULES/gaussG.f90
integrals.o : MODULES/integrals.f90 dimensions.o
	$(F90) $(FFLAGS) -c MODULES/integrals.f90
interactions.o : MODULES/interactions.f90
	$(F90) $(FFLAGS) -c MODULES/interactions.f90
kpoints.o : MODULES/kpoints.f90
	$(F90) $(FFLAGS) -c MODULES/kpoints.f90
mpi_declarations.o : MODULES/mpi_declarations.f90 
	$(F90) $(FFLAGS) -c MODULES/mpi_declarations.f90
neighbor_map.o : MODULES/neighbor_map.f90 dimensions.o
	$(F90) $(FFLAGS) -c MODULES/neighbor_map.f90
ordern.o : MODULES/ordern.f90 dimensions.o
	$(F90) $(FFLAGS) -c MODULES/ordern.f90
umbrella.o : MODULES/umbrella.f90
	$(F90) $(FFLAGS) -c MODULES/umbrella.f90
steered.o : MODULES/steered.f90
	$(F90) $(FFLAGS) -c MODULES/steered.f90
optimization.o : MODULES/optimization.f90
	$(F90) $(FFLAGS) -c MODULES/optimization.f90
module_dos.o : MODULES/module_dos.f90
	$(F90) $(FFLAGS) -c MODULES/module_dos.f90
dynamo.o : MODULES/dynamo.f90
	$(F90) $(FFLAGS) -c MODULES/dynamo.f90
cproc.o : MODULES/cproc.f90
	$(F90) $(FFLAGS) -c MODULES/cproc.f90
noseHoover.o : MODULES/noseHoover.f90
	$(F90) $(FFLAGS) -c MODULES/noseHoover.f90
scf.o : MODULES/scf.f90
	$(F90) $(FFLAGS) -c MODULES/scf.f90
grid.o : MODULES/grid.f90
	$(F90) $(FFLAGS) -c MODULES/grid.f90
wavefunction.o : MODULES/wavefunction.f90
	$(F90) $(FFLAGS) -c MODULES/wavefunction.f90
neb_module.o : MODULES/neb_module.f90
	$(F90) $(FFLAGS) -c MODULES/neb_module.f90
vnneutral.o : MODULES/vnneutral.f90
	$(F90) $(FFLAGS) -c MODULES/vnneutral.f90
transport.o : MODULES/transport.f90
	$(F90) $(FFLAGS) -c MODULES/transport.f90
matmultmod.o : MODULES/matmultmod.f90
	$(F90) $(FFLAGS) -c MODULES/matmultmod.f90
options.o : MODULES/options.f90
	$(F90) $(FFLAGS) -c MODULES/options.f90
outputs.o : MODULES/outputs.f90
	$(F90) $(FFLAGS) -c MODULES/outputs.f90
MD.o : MODULES/MD.f90
	$(F90) $(FFLAGS) -c MODULES/MD.f90
energy.o : MODULES/energy.f90
	$(F90) $(FFLAGS) -c MODULES/energy.f90
mpi_main.o : MODULES/mpi_main.f90
	$(F90) $(FFLAGS) -c MODULES/mpi_main.f90
tdse.o : MODULES/tdse.f90
	$(F90) $(FFLAGS) -c MODULES/tdse.f90
bias.o : MODULES/bias.f90
	$(F90) $(FFLAGS) -c MODULES/bias.f90
nonadiabatic.o : MODULES/nonadiabatic.f90
	$(F90) $(FFLAGS) -c MODULES/nonadiabatic.f90
hartree.o : MODULES/hartree.f90
	$(F90) $(FFLAGS) -c MODULES/hartree.f90
fsockets.o : MODULES/fsockets.f90
	$(F90) $(FFLAGS) -c MODULES/fsockets.f90
fb_socket.o : MODULES/fb_socket.f90
	$(F90) $(FFLAGS) -c MODULES/fb_socket.f90
sockets.o : MODULES/sockets.c
	$(CC) $(CFLAGS) -c MODULES/sockets.c
# *****************************************************************************
# modules_c
# *****************************************************************************
classicMD.o : MODULES/classicMD.f90
	$(F90) $(FFLAGS) -c MODULES/classicMD.f90


# *****************************************************************************
# main
# *****************************************************************************
fireball.o : fireball.f90 $(MODULES)
	$(F90) $(FFLAGS) -I/usr/local/include -c fireball.f90
ifneq (,$(findstring SCALAPACK, ${METHOD}))
fireball_server_amber.o : server/fireball_server_amber.f90 $(MODULES) 
	$(F90) -c server/fireball_server_amber.f90
endif
ifneq (,$(findstring SCALAPACK, ${METHOD}))
fireball_server.o : server/fireball_server.f90 $(MODULES) 
	$(F90) -c server/fireball_server.f90
endif
ifneq (,$(findstring OPENMP, ${PARALLEL}))
fireball_server.o : server/fireball_server_OPENMP.f90 $(MODULES) 
	$(F90) -c server/fireball_server_OPENMP.f90 -o fireball_server.o
endif
# *****************************************************************************
# allocate objects
# *****************************************************************************
allocate_f.o : ALLOCATIONS/allocate_f.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_f.f90
allocate_h.o : ALLOCATIONS/allocate_h.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_h.f90
allocate_neigh.o : ALLOCATIONS/allocate_neigh.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_neigh.f90
allocate_ordern.o : ALLOCATIONS/allocate_ordern.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_ordern.f90
allocate_rho.o : ALLOCATIONS/allocate_rho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_rho.f90
allocate_umb.o : ALLOCATIONS/allocate_umb.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_umb.f90
allocate_steered.o : ALLOCATIONS/allocate_steered.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_steered.f90
reallocate_f.o : ALLOCATIONS/reallocate_f.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/reallocate_f.f90
reallocate_h.o : ALLOCATIONS/reallocate_h.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/reallocate_h.f90
reallocate_neigh.o : ALLOCATIONS/reallocate_neigh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/reallocate_neigh.f90
reallocate_rho.o : ALLOCATIONS/reallocate_rho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/reallocate_rho.f90
allocate_dos.o : ALLOCATIONS/allocate_dos.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_dos.f90
allocate_grid.o : ALLOCATIONS/allocate_grid.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_grid.f90
allocate_trans.o : ALLOCATIONS/allocate_trans.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ALLOCATIONS/allocate_trans.f90

# *****************************************************************************
# assembler objects
# *****************************************************************************
average_rho.o : ASSEMBLERS/average_rho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/average_rho.f90
average_ca_rho.o : ASSEMBLERS/average_ca_rho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/average_ca_rho.f90
assemble_hxc_1c.o : ASSEMBLERS/assemble_hxc_1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_hxc_1c.f90
assemble_olsxc_1c.o : ASSEMBLERS/assemble_olsxc_1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_olsxc_1c.f90
assemble_2c.o : ASSEMBLERS/assemble_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_2c.f90
assemble_2c_PP.o : ASSEMBLERS/assemble_2c_PP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_2c_PP.f90
assemble_2c_ordern_final.o : ASSEMBLERS/assemble_2c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_2c_ordern_final.f90
assemble_2c_ordern_init.o : ASSEMBLERS/assemble_2c_ordern_init.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_2c_ordern_init.f90
assemble_3c.o : ASSEMBLERS/assemble_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_3c.f90
assemble_3c_PP.o : ASSEMBLERS/assemble_3c_PP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_3c_PP.f90
assemble_3c_ordern_final.o : ASSEMBLERS/assemble_3c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_3c_ordern_final.f90
assemble_ca_2c.o : ASSEMBLERS/assemble_ca_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_2c.f90
assemble_ca_2c_ordern_final.o : ASSEMBLERS/assemble_ca_2c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_2c_ordern_final.f90
assemble_ca_3c.o : ASSEMBLERS/assemble_ca_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_3c.f90
assemble_ca_3c_ordern_final.o : ASSEMBLERS/assemble_ca_3c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_3c_ordern_final.f90
assemble_eh_2c.o : ASSEMBLERS/assemble_eh_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_eh_2c.f90
assemble_eh_2c_ordern_final.o : ASSEMBLERS/assemble_eh_2c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_eh_2c_ordern_final.f90
assemble_eh_usr.o : ASSEMBLERS/assemble_eh_usr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_eh_usr.f90
assemble_F.o : ASSEMBLERS/assemble_F.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_F.f90
assemble_hxc_2c.o : ASSEMBLERS/assemble_hxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_hxc_2c.f90
assemble_hxc_3c.o : ASSEMBLERS/assemble_hxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_hxc_3c.f90
assemble_lr.o : ASSEMBLERS/assemble_lr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_lr.f90
assemble_sVNL.o : ASSEMBLERS/assemble_sVNL.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_sVNL.f90
assemble_usr.o : ASSEMBLERS/assemble_usr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_usr.f90
buildh.o : ASSEMBLERS/buildh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/buildh.f90
assemble_olsxc_on.o : ASSEMBLERS/assemble_olsxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_olsxc_on.f90
assemble_olsxc_off.o : ASSEMBLERS/assemble_olsxc_off.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_olsxc_off.f90
build_olsxc_on.o : ASSEMBLERS/build_olsxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_olsxc_on.f90
build_ca_olsxc_on.o : ASSEMBLERS/build_ca_olsxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_ca_olsxc_on.f90
build_olsxc_off.o : ASSEMBLERS/build_olsxc_off.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_olsxc_off.f90
assemble_snxc_on.o : ASSEMBLERS/assemble_snxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_snxc_on.f90
assemble_snxc_off.o : ASSEMBLERS/assemble_snxc_off.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_snxc_off.f90
build_snxc_on.o : ASSEMBLERS/build_snxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_snxc_on.f90
build_ca_snxc_on.o : ASSEMBLERS/build_ca_snxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_ca_snxc_on.f90
build_snxc_off.o : ASSEMBLERS/build_snxc_off.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_snxc_off.f90
assemble_h.o : ASSEMBLERS/assemble_h.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_h.f90
assemble_hxc.o : ASSEMBLERS/assemble_hxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_hxc.f90
assemble_eh.o : ASSEMBLERS/assemble_eh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_eh.f90
assemble_mcweda.o : ASSEMBLERS/assemble_mcweda.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_mcweda.f90
getenergy.o : ASSEMBLERS/getenergy.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy.f90
getenergy_mcweda.o : ASSEMBLERS/getenergy_mcweda.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy_mcweda.f90
getenergy_hxc.o : ASSEMBLERS/getenergy_hxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy_hxc.f90
getenergy_eh.o : ASSEMBLERS/getenergy_eh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy_eh.f90
assemble_h_ks.o : ASSEMBLERS/assemble_h_ks.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_h_ks.f90
getenergy_KS.o : ASSEMBLERS/getenergy_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy_KS.f90
assemble_S.o : ASSEMBLERS/assemble_S.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_S.f90
assemble_2c_S.o : ASSEMBLERS/assemble_2c_S.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_2c_S.f90
assemble_hartree.o : ASSEMBLERS/assemble_hartree.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_hartree.f90
assemble_scissor.o : ASSEMBLERS/assemble_scissor.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_scissor.f90
assemble_ca_2c_dip.o : ASSEMBLERS/assemble_ca_2c_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_2c_dip.f90
assemble_ca_3c_dip.o : ASSEMBLERS/assemble_ca_3c_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_ca_3c_dip.f90
assemble_lr_dip.o : ASSEMBLERS/assemble_lr_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_lr_dip.f90
assemble_xczw.o : ASSEMBLERS/assemble_xczw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_xczw.f90
assemble_zw_1c_na.o : ASSEMBLERS/assemble_zw_1c_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_zw_1c_na.f90
assemble_zw_2c_ct.o : ASSEMBLERS/assemble_zw_2c_ct.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_zw_2c_ct.f90
assemble_zw_3c_ct.o : ASSEMBLERS/assemble_zw_3c_ct.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_zw_3c_ct.f90
assemble_zw_off_na.o : ASSEMBLERS/assemble_zw_off_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_zw_off_na.f90
assemble_zw_on_na.o : ASSEMBLERS/assemble_zw_on_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_zw_on_na.f90
build_zw_on_na.o : ASSEMBLERS/build_zw_on_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_zw_on_na.f90
build_zw_off_na.o : ASSEMBLERS/build_zw_off_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/build_zw_off_na.f90
getenergy_zw.o : ASSEMBLERS/getenergy_zw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/getenergy_zw.f90
assemble_1c_vdip.o : ASSEMBLERS/assemble_1c_vdip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ASSEMBLERS/assemble_1c_vdip.f90
# *****************************************************************************
# Grid objects
# *****************************************************************************
assemble_KS_den0.o : GRID/assemble_KS_den0.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_den0.f90
assemble_KS_mat.o : GRID/assemble_KS_mat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_mat.f90
assemble_KS_den.o : GRID/assemble_KS_den.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_den.f90
assemble_KS_usr.o : GRID/assemble_KS_usr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_usr.f90
assemble_KS_dcc.o : GRID/assemble_KS_dcc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_dcc.f90
laplace_fft.o : GRID/laplace_fft.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/laplace_fft.f90
mixer_KS.o : GRID/mixer_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/mixer_KS.f90
writeout_charges_KS.o : GRID/writeout_charges_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/writeout_charges_KS.f90
writeout_xsf.o : GRID/writeout_xsf.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/writeout_xsf.f90
#psi2mesh.o : GRID/psi2mesh.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c GRID/psi2mesh.f90
#psi22mesh.o : GRID/psi22mesh.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c GRID/psi22mesh.f90
assemble_KS_vna.o : GRID/assemble_KS_vna.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/assemble_KS_vna.f90
get_Hort.o : GRID/get_Hort.f90 $(MODULES)
	$(F90) $(FFLAGS) -c GRID/get_Hort.f90

# *****************************************************************************
# Dassembler objects
# *****************************************************************************
Dassemble_2c.o : DASSEMBLERS/Dassemble_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_2c.f90
Dassemble_2c_ordern_final.o : DASSEMBLERS/Dassemble_2c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_2c_ordern_final.f90
Dassemble_3c.o : DASSEMBLERS/Dassemble_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_3c.f90
Dassemble_3c_ordern_final.o : DASSEMBLERS/Dassemble_3c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_3c_ordern_final.f90
Dassemble_ca_2c.o : DASSEMBLERS/Dassemble_ca_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_2c.f90
Dassemble_ca_2c_ordern_final.o : DASSEMBLERS/Dassemble_ca_2c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_2c_ordern_final.f90
Dassemble_ca_3c.o : DASSEMBLERS/Dassemble_ca_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_3c.f90
Dassemble_ca_3c_ordern_final.o : DASSEMBLERS/Dassemble_ca_3c_ordern_final.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_3c_ordern_final.f90
Dassemble_eh_2c.o : DASSEMBLERS/Dassemble_eh_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_eh_2c.f90
Dassemble_hxc_2c.o : DASSEMBLERS/Dassemble_hxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_hxc_2c.f90
Dassemble_hxc_3c.o : DASSEMBLERS/Dassemble_hxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_hxc_3c.f90
Dassemble_lr.o : DASSEMBLERS/Dassemble_lr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_lr.f90
Dassemble_lr_OMP.o : DASSEMBLERS/Dassemble_lr_OMP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_lr_OMP.f90
Dassemble_snxc_on.o : DASSEMBLERS/Dassemble_snxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_snxc_on.f90
Dassemble_olsxc_on.o : DASSEMBLERS/Dassemble_olsxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_olsxc_on.f90
Dassemble_olsxc_2c.o : DASSEMBLERS/Dassemble_olsxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_olsxc_2c.f90
Dassemble_olsxc_3c.o : DASSEMBLERS/Dassemble_olsxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_olsxc_3c.f90
Dassemble_snxc_2c.o : DASSEMBLERS/Dassemble_snxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_snxc_2c.f90
Dassemble_snxc_3c.o : DASSEMBLERS/Dassemble_snxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_snxc_3c.f90
Dassemble_2c_PP.o : DASSEMBLERS/Dassemble_2c_PP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_2c_PP.f90
Dassemble_3c_PP.o : DASSEMBLERS/Dassemble_3c_PP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_3c_PP.f90
Dassemble_ca_snxc_on.o : DASSEMBLERS/Dassemble_ca_snxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_snxc_on.f90
Dassemble_ca_olsxc_on.o : DASSEMBLERS/Dassemble_ca_olsxc_on.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_olsxc_on.f90
Dassemble_ca_snxc_3c.o : DASSEMBLERS/Dassemble_ca_snxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_snxc_3c.f90
Dassemble_ca_olsxc_3c.o : DASSEMBLERS/Dassemble_ca_olsxc_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_olsxc_3c.f90
Dassemble_ca_snxc_2c.o : DASSEMBLERS/Dassemble_ca_snxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_snxc_2c.f90
Dassemble_ca_olsxc_2c.o : DASSEMBLERS/Dassemble_ca_olsxc_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_olsxc_2c.f90
Dassemble_ca_2c_dip.o : DASSEMBLERS/Dassemble_ca_2c_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_2c_dip.f90
Dassemble_ca_3c_dip.o : DASSEMBLERS/Dassemble_ca_3c_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_ca_3c_dip.f90
Dassemble_lr_dip.o : DASSEMBLERS/Dassemble_lr_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_lr_dip.f90
Dassemble_zw_2c_ct.o : DASSEMBLERS/Dassemble_zw_2c_ct.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_zw_2c_ct.f90
Dassemble_zw_2c_na.o : DASSEMBLERS/Dassemble_zw_2c_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_zw_2c_na.f90
Dassemble_zw_3c_ct.o : DASSEMBLERS/Dassemble_zw_3c_ct.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_zw_3c_ct.f90
Dassemble_zw_3c_na.o : DASSEMBLERS/Dassemble_zw_3c_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_zw_3c_na.f90
Dassemble_zw_on_na.o : DASSEMBLERS/Dassemble_zw_on_na.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/Dassemble_zw_on_na.f90
getforces_zw.o : DASSEMBLERS/getforces_zw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_zw.f90
getforces_mcweda.o : DASSEMBLERS/getforces_mcweda.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_mcweda.f90
getforces_mdet.o : NAC/getforces_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/getforces_mdet.f90
getforces_eh.o : DASSEMBLERS/getforces_eh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_eh.f90
getforces_hxc.o : DASSEMBLERS/getforces_hxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_hxc.f90
getforces_KS.o : DASSEMBLERS/getforces_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_KS.f90
getforces_classic.o : DASSEMBLERS/getforces_classic.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_classic.f90
getforces_classic_RGL.o : DASSEMBLERS/getforces_classic_RGL.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_classic_RGL.f90
getforces_classic_tersoff.o : DASSEMBLERS/getforces_classic_tersoff.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_classic_tersoff.f90
getforces_classic_vdw.o : DASSEMBLERS/getforces_classic_vdw.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_classic_vdw.f90
getforces.o : DASSEMBLERS/getforces.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces.f90
getforces_socket.o : DASSEMBLERS/getforces_socket.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DASSEMBLERS/getforces_socket.f90
#******************************************************************************
# Dos objects
#******************************************************************************
dos.o : DOS/dos.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DOS/dos.f90
invierte.o : DOS/invierte.f $(MODULES)
	$(F90) $(FFLAGS) -c DOS/invierte.f
writeout_dos.o : DOS/writeout_dos.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DOS/writeout_dos.f90
writeout_dosng.o : DOS/writeout_dosng.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DOS/writeout_dosng.f90
writeout_atom.o : DOS/writeout_atom.f90 $(MODULES)
	$(F90) $(FFLAGS) -c DOS/writeout_atom.f90
hoppings.o : DOS/hoppings.f90
	$(F90) $(FFLAGS) -c DOS/hoppings.f90
hamilt_atom.o : DOS/hamilt_atom.f90
	$(F90) $(FFLAGS) -c DOS/hamilt_atom.f90

# *****************************************************************************
# form_rho objects
# *****************************************************************************
chebft.o : FORM_RHO/chebft.f90 $(MODULES)
	$(F90) $(DFLAGS) -c FORM_RHO/chebft.f90
denmat.o : FORM_RHO/denmat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmat.f90
denmat_fk.o : FORM_RHO/denmat_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmat_fk.f90
denmat_es.o : FORM_RHO/denmat_es.f90
	$(F90) $(FFLAGS) -c FORM_RHO/denmat_es.f90
denmata_ordern.o : FORM_RHO/denmata_ordern.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmata_ordern.f90
denmatb_ordern.o : FORM_RHO/denmatb_ordern.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmatb_ordern.f90
denmatc_ordern.o : FORM_RHO/denmatc_ordern.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmatc_ordern.f90
denmata_ordern_fk.o : FORM_RHO/denmata_ordern_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmata_ordern_fk.f90
denmatb_ordern_fk.o : FORM_RHO/denmatb_ordern_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmatb_ordern_fk.f90
denmatc_ordern_fk.o : FORM_RHO/denmatc_ordern_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmatc_ordern_fk.f90
fermie.o : FORM_RHO/fermie.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/fermie.f90
formrho_sparse.o : FORM_RHO/formrho_sparse.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/formrho_sparse.f90
ss12.o : FORM_RHO/ss12.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/ss12.f90
ss12_fk.o : FORM_RHO/ss12_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/ss12_fk.f90
build_rho.o : FORM_RHO/build_rho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/build_rho.f90
denmat_KS.o : FORM_RHO/denmat_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/denmat_KS.f90
build_rho_KS.o : FORM_RHO/build_rho_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/build_rho_KS.f90
build_Ji.o : FORM_RHO/build_Ji.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/build_Ji.f90
koopman.o : FORM_RHO/koopman.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/koopman.f90
build_nij.o : FORM_RHO/build_nij.f90 $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/build_nij.f90
project_eh.o : FORM_RHO/project_eh.f90  $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/project_eh.f90
project_wfmdet.o : FORM_RHO/project_wfmdet.f90  $(MODULES)
	$(F90) $(FFLAGS) -c FORM_RHO/project_wfmdet.f90
# *****************************************************************************
# initializers objects
# *****************************************************************************
initconstraints.o : INITIALIZERS/initconstraints.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initconstraints.f90
diagnostics.o : INITIALIZERS/diagnostics.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/diagnostics.f90
initamat.o : INITIALIZERS/initamat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initamat.f90
initcharges.o : INITIALIZERS/initcharges.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initcharges.f90
initconstants.o : INITIALIZERS/initconstants.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initconstants.f90
initkpoints.o : INITIALIZERS/initkpoints.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initkpoints.f90
initmasses.o : INITIALIZERS/initmasses.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initmasses.f90
initneighbors.o : INITIALIZERS/initneighbors.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initneighbors.f90
welcome.o : INITIALIZERS/welcome.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/welcome.f90
init_MPI.o : INITIALIZERS/init_MPI.f90 $(MODULES) 
	$(F90) $(FFLAGS) $(MODMPI) -c INITIALIZERS/init_MPI.f90
init_noMPI.o : INITIALIZERS/init_noMPI.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/init_noMPI.f90
initatomicE.o : INITIALIZERS/initatomicE.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initatomicE.f90
initboxes.o : INITIALIZERS/initboxes.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initboxes.f90
make_mu2shell.o : INITIALIZERS/make_mu2shell.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/make_mu2shell.f90
make_munu.o : INITIALIZERS/make_munu.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/make_munu.f90
make_munuPP.o : INITIALIZERS/make_munuPP.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/make_munuPP.f90
make_munuS.o : INITIALIZERS/make_munuS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/make_munuS.f90
make_munuDipY.o : INITIALIZERS/make_munuDipY.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/make_munuDipY.f90
make_munuDipX.o : INITIALIZERS/make_munuDipX.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/make_munuDipX.f90
restart.o : INITIALIZERS/restart.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/restart.f90
zero_ang_mom.o : INITIALIZERS/zero_ang_mom.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/zero_ang_mom.f90
initNH.o : INITIALIZERS/initNH.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/initNH.f90
getkpoints.o : INITIALIZERS/getkpoints.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INITIALIZERS/getkpoints.f90
greatKAuto.o : INITIALIZERS/greatKAuto.f90 
	$(F90) $(FFLAGS) -c INITIALIZERS/greatKAuto.f90
greatKsubsAuto.o : INITIALIZERS/greatKsubsAuto.f90  
	$(F90) $(FFLAGS) -c INITIALIZERS/greatKsubsAuto.f90
initgrid.o : INITIALIZERS/initgrid.f90
	$(F90) $(FFLAGS) -c INITIALIZERS/initgrid.f90
initdenmat.o : INITIALIZERS/initdenmat.f90
	$(F90) $(FFLAGS) -c INITIALIZERS/initdenmat.f90
get_info_orbital.o : INITIALIZERS/get_info_orbital.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/get_info_orbital.f90
initbasics.o : INITIALIZERS/initbasics.f90
	$(F90) $(FFLAGS) -c INITIALIZERS/initbasics.f90
initcharges_KS.o : INITIALIZERS/initcharges_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initcharges_KS.f90
initcDFT.o : INITIALIZERS/initcDFT.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INITIALIZERS/initcDFT.f90

# *****************************************************************************
# interactions objects
# *****************************************************************************
cl_value.o : INTERACTIONS/cl_value.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INTERACTIONS/cl_value.f90
Dtrescentros.o : INTERACTIONS/Dtrescentros.f90 $(MODULES) 
	$(F90) $(FFLAGS) -c INTERACTIONS/Dtrescentros.f90
DtrescentrosS.o : INTERACTIONS/DtrescentrosS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosS.f90
doscentros.o : INTERACTIONS/doscentros.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentros.f90
doscentrosPP.o : INTERACTIONS/doscentrosPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosPP.f90
doscentrosS.o : INTERACTIONS/doscentrosS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosS.f90
doscentrosDipY.o : INTERACTIONS/doscentrosDipY.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosDipY.f90
doscentrosDipX.o : INTERACTIONS/doscentrosDipX.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosDipX.f90
get_ewald.o : INTERACTIONS/get_ewald.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/get_ewald.f90
get_ewald_OMP.o : INTERACTIONS/get_ewald_OMP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/get_ewald_OMP.f90
get_vdw.o : INTERACTIONS/get_vdw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/get_vdw.f90
getHarmonic.o : INTERACTIONS/getHarmonic.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/getHarmonic.f90
internalLambda.o : INTERACTIONS/internalLambda.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/internalLambda.f90
trescentros.o : INTERACTIONS/trescentros.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentros.f90
trescentrosS.o : INTERACTIONS/trescentrosS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosS.f90
unocentros.o : INTERACTIONS/unocentros.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/unocentros.f90
smoother.o : INTERACTIONS/smoother.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/smoother.f90

dosgaussians.o : INTERACTIONS/dosgaussians.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/dosgaussians.f90
gelements_VXC.o : INTERACTIONS/gelements_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelements_VXC.f90
Dgelements_VXC.o : INTERACTIONS/Dgelements_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/Dgelements_VXC.f90
trescentrosGHXC_VXC.o : INTERACTIONS/trescentrosGHXC_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosGHXC_VXC.f90
DtrescentrosGHXC_VXC.o : INTERACTIONS/DtrescentrosGHXC_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosGHXC_VXC.f90

gelementsG_VNA.o : INTERACTIONS/gelementsG_VNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsG_VNA.f90
gelementsG_VXC.o : INTERACTIONS/gelementsG_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsG_VXC.f90
gelementsGS_VXC.o : INTERACTIONS/gelementsGS_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsGS_VXC.f90
DgelementsG_VNA.o : INTERACTIONS/DgelementsG_VNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsG_VNA.f90
DgelementsG_VXC.o : INTERACTIONS/DgelementsG_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsG_VXC.f90
DgelementsGS_VXC.o : INTERACTIONS/DgelementsGS_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsGS_VXC.f90
gelementsGS_overlap.o : INTERACTIONS/gelementsGS_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsGS_overlap.f90
DgelementsGS_overlap.o : INTERACTIONS/DgelementsGS_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsGS_overlap.f90
doscentrosGS_overlap.o : INTERACTIONS/doscentrosGS_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosGS_overlap.f90

gelementsG_overlap.o : INTERACTIONS/gelementsG_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsG_overlap.f90
DgelementsG_overlap.o : INTERACTIONS/DgelementsG_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsG_overlap.f90
doscentrosG_overlap.o : INTERACTIONS/doscentrosG_overlap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/doscentrosG_overlap.f90

trescentrosG_VNA.o : INTERACTIONS/trescentrosG_VNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosG_VNA.f90
trescentrosG_VXC.o : INTERACTIONS/trescentrosG_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosG_VXC.f90
trescentrosGS_VXC.o : INTERACTIONS/trescentrosGS_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosGS_VXC.f90
DtrescentrosG_VNA.o : INTERACTIONS/DtrescentrosG_VNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosG_VNA.f90
DtrescentrosG_VXC.o : INTERACTIONS/DtrescentrosG_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosG_VXC.f90
DtrescentrosGS_VXC.o : INTERACTIONS/DtrescentrosGS_VXC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosGS_VXC.f90

trescentrosG_VNA_SH.o : INTERACTIONS/trescentrosG_VNA_SH.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/trescentrosG_VNA_SH.f90
DtrescentrosG_VNA_SH.o : INTERACTIONS/DtrescentrosG_VNA_SH.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DtrescentrosG_VNA_SH.f90
gelementsG_VNA_SH.o : INTERACTIONS/gelementsG_VNA_SH.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsG_VNA_SH.f90
DgelementsG_VNA_SH.o : INTERACTIONS/DgelementsG_VNA_SH.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsG_VNA_SH.f90
gelementsG_VNA_SH_RNA.o : INTERACTIONS/gelementsG_VNA_SH_RNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/gelementsG_VNA_SH_RNA.f90
DgelementsG_VNA_SH_RNA.o : INTERACTIONS/DgelementsG_VNA_SH_RNA.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/DgelementsG_VNA_SH_RNA.f90
tester2c.o : INTERACTIONS/tester2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERACTIONS/tester2c.f90

# *****************************************************************************
# interpolation objects
# *****************************************************************************
buildspline_1d.o : INTERPOLATERS/buildspline_1d.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/buildspline_1d.f90
interpolate_1d.o : INTERPOLATERS/interpolate_1d.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/interpolate_1d.f90
interpolate_2d.o : INTERPOLATERS/interpolate_2d.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/interpolate_2d.f90
setterp_2d.o : INTERPOLATERS/setterp_2d.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/setterp_2d.f90
recover_2c.o : INTERPOLATERS/recover_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_2c.f90
recoverC.o : INTERPOLATERS/recoverC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recoverC.f90
recover_PP.o : INTERPOLATERS/recover_PP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_PP.f90
recover_S.o : INTERPOLATERS/recover_S.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_S.f90
recover_3c.o : INTERPOLATERS/recover_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_3c.f90
recover_2cDipY.o : INTERPOLATERS/recover_2cDipY.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_2cDipY.f90
recover_2cDipX.o : INTERPOLATERS/recover_2cDipX.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/recover_2cDipX.f90
buildspline2_1d.o : INTERPOLATERS/buildspline2_1d.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/buildspline2_1d.f90
getpsi.o : INTERPOLATERS/getpsi.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/getpsi.f90
getYlm.o : INTERPOLATERS/getYlm.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/getYlm.f90
getvna.o : INTERPOLATERS/getvna.f90 $(MODULES)
	$(F90) $(FFLAGS) -c INTERPOLATERS/getvna.f90

# *****************************************************************************
# loops objects
# *****************************************************************************
scf_loop.o : LOOPS/scf_loop.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/scf_loop.f90
scf_loop_harris.o : LOOPS/scf_loop_harris.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/scf_loop_harris.f90
main_loop.o : LOOPS/main_loop.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop.f90
main_loop_MD.o : LOOPS/main_loop_MD.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_MD.f90
main_loop_NEB.o : LOOPS/main_loop_NEB.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_NEB.f90
main_loop_CG.o : LOOPS/main_loop_CG.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_CG.f90
main_loop_importrho.o : LOOPS/main_loop_importrho.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_importrho.f90	
main_loop_DM.o : LOOPS/main_loop_DM.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_DM.f90
scf_loop_ks.o : LOOPS/scf_loop_ks.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/scf_loop_ks.f90
main_loop_TDSE.o : LOOPS/main_loop_TDSE.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_TDSE.f90
main_loop_MDET.o : NAC/main_loop_MDET.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/main_loop_MDET.f90
main_loop_MIN.o : LOOPS/main_loop_MIN.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_MIN.f90
main_loop_NAC.o : LOOPS/main_loop_NAC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_NAC.f90
main_loop_FIRE.o : LOOPS/main_loop_FIRE.f90 $(MODULES)
	$(F90) $(FFLAGS) -c LOOPS/main_loop_FIRE.f90
main_loop_socket.o : LOOPS/main_loop_socket.f90
	$(F90) $(FFLAGS) -c LOOPS/main_loop_socket.f90

# *****************************************************************************
# molecular dynamics objects
# *****************************************************************************
cross.o : MD/cross.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/cross.f90
factorial.o : MD/factorial.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/factorial.f90
gaussT.o : MD/gaussT.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/gaussT.f90
imaged.o : MD/imaged.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/imaged.f90
corrector.o :MD/corrector.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/corrector.f90
predictor.o :MD/predictor.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/predictor.f90
setgear.o :MD/setgear.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/setgear.f90
phimat.o :MD/phimat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/phimat.f90
bvec.o :MD/bvec.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/bvec.f90
soldm.o : MD/soldm.f90  $(MODULES)
	$(F90) $(FFLAGS) -c MD/soldm.f90
NHCThermostat.o: MD/NHCThermostat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/NHCThermostat.f90
get_enk.o: MD/get_enk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/get_enk.f90
writeHNose.o: MD/writeHNose.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/writeHNose.f90
resetNHC.o: MD/resetNHC.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/resetNHC.f90
move_ions.o: MD/move_ions.f90 $(MODULES)
	$(F90) $(FFLAGS) -c MD/move_ions.f90

# *****************************************************************************
# neighbors objects
# *****************************************************************************
backnay.o : NEIGHBORS/backnay.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/backnay.f90
common_neighbors.o : NEIGHBORS/common_neighbors.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/common_neighbors.f90
find_neigh_max.o : NEIGHBORS/find_neigh_max.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/find_neigh_max.f90
find_neigh_max_class.o : NEIGHBORS/find_neigh_max_class.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c NEIGHBORS/find_neigh_max_class.f90
mpairnay.o : NEIGHBORS/mpairnay.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/mpairnay.f90
neighbors.o : NEIGHBORS/neighbors.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/neighbors.f90
neighbors_pairs.o : NEIGHBORS/neighbors_pairs.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/neighbors_pairs.f90
neighborsPP.o : NEIGHBORS/neighborsPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/neighborsPP.f90
find_neighPP_max.o : NEIGHBORS/find_neighPP_max.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/find_neighPP_max.f90
common_neighborsPP.o : NEIGHBORS/common_neighborsPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/common_neighborsPP.f90
num_neigh_tot.o : NEIGHBORS/num_neigh_tot.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NEIGHBORS/num_neigh_tot.f90

# *****************************************************************************
# pressure objects
# *****************************************************************************
hmetric.o : PRESSURE/hmetric.f90 $(MODULES)
	$(F90) $(FFLAGS) -c PRESSURE/hmetric.f90
initpressure.o : PRESSURE/initpressure.f90 $(MODULES)
	$(F90) $(FFLAGS) -c PRESSURE/initpressure.f90
invert3x3.o : PRESSURE/invert3x3.f90 $(MODULES)
	$(F90) $(FFLAGS) -c PRESSURE/invert3x3.f90


# *****************************************************************************
# readfiles objects
# *****************************************************************************
append_string.o : READFILES/append_string.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/append_string.f90
#read_1c.o : READFILES/read_1c.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c READFILES/read_1c.f90
read_1c.o : NAC/read_1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/read_1c.f90
read_2c.o : READFILES/read_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/read_2c.f90
read_3c.o : READFILES/read_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/read_3c.f90
readbarrier.o : READFILES/readbarrier.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readbarrier.f90
readbasis.o : READFILES/readbasis.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readbasis.f90
readdata_2c.o : READFILES/readdata_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata_2c.f90
readdata_3c.o : READFILES/readdata_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata_3c.f90
readfragments.o : READFILES/readfragments.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readfragments.f90
readgaussG.o : READFILES/readgaussG.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readgaussG.f90
readheader_2c.o : READFILES/readheader_2c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readheader_2c.f90
readheader_3c.o : READFILES/readheader_3c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readheader_3c.f90
readinfo.o : READFILES/readinfo.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readinfo.f90
readlvs.o : READFILES/readlvs.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readlvs.f90
readparam.o : READFILES/readparam.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readparam.f90
#readparam.o : NAC/readparam.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c NAC/readparam.f90
readphi.o : READFILES/readphi.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readphi.f90
readpressure.o : READFILES/readpressure.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readpressure.f90
readquench.o : READFILES/readquench.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readquench.f90
readsa.o : READFILES/readsa.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readsa.f90
readvdw.o : READFILES/readvdw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readvdw.f90
readcgo.o : READFILES/readcgo.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readcgo.f90
readdos.o : READFILES/readdos.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdos.f90
findFdata.o : READFILES/findFdata.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/findFdata.f90
readgrid.o : READFILES/readgrid.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readgrid.f90
read_wf.o : READFILES/read_wf.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/read_wf.f90
read_vna.o : READFILES/read_vna.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/read_vna.f90
readtrans.o : READFILES/readtrans.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readtrans.f90
readbind.o : READFILES/readbind.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readbind.f90
readhop.o : READFILES/readhop.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readhop.f90
readdata.o : READFILES/readdata.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata.f90
readdata_hxc.o : READFILES/readdata_hxc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata_hxc.f90
#readdata_mcweda.o : READFILES/readdata_mcweda.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c READFILES/readdata_mcweda.f90
#readdata_xczw.o : READFILES/readdata_xczw.f90 $(MODULES)
#	$(F90) $(FFLAGS) -c READFILES/readdata_xczw.f90
readdata_mcweda.o : NAC/readdata_mcweda.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/readdata_mcweda.f90
readdata_xczw.o : NAC/readdata_xczw.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/readdata_xczw.f90
readdata_eh.o : READFILES/readdata_eh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata_eh.f90
checksum_options.o : READFILES/checksum_options.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/checksum_options.f90
readdata_KS.o : READFILES/readdata_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readdata_KS.f90
getsections.o : READFILES/getsections.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/getsections.f90
readdata_classicMD.o : READFILES/readdata_classicMD.f90 $(MODULES_C)
	$(F90) $(FFLAGS) -c READFILES/readdata_classicMD.f90
readhartree.o : READFILES/readhartree.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readhartree.f90
readephc.o : READFILES/readephc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c READFILES/readephc.f90

# *****************************************************************************
# rotations objects
# *****************************************************************************
chooser.o : ROTATIONS/chooser.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/chooser.f90
chooserd.o : ROTATIONS/chooserd.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/chooserd.f90
deps2center.o : ROTATIONS/deps2center.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/deps2center.f90
deps3center.o : ROTATIONS/deps3center.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/deps3center.f90
makeDmat.o : ROTATIONS/makeDmat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/makeDmat.f90
makeDmatPP.o : ROTATIONS/makeDmatPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/makeDmatPP.f90
rotate.o : ROTATIONS/rotate.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/rotate.f90
rotatePP.o : ROTATIONS/rotatePP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/rotatePP.f90
rotated.o : ROTATIONS/rotated.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/rotated.f90
rotatedPP.o : ROTATIONS/rotatedPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/rotatedPP.f90
twister.o : ROTATIONS/twister.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/twister.f90
twisterd.o : ROTATIONS/twisterd.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/twisterd.f90
epsilon.o : ROTATIONS/epsilon.f90 $(MODULES)
	$(F90) $(FFLAGS) -c ROTATIONS/epsilon.f90

# *****************************************************************************
# socket objects
# *****************************************************************************
get_geometry.o : SOCKETS/get_geometry.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOCKETS/get_geometry.f90
create_socket.o : SOCKETS/create_socket.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOCKETS/create_socket.f90
send_geometry.o : SOCKETS/send_geometry.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOCKETS/send_geometry.f90
sendrecv.o : SOCKETS/sendrecv.c SOCKETS/mysocks.h
	$(CC) $(CFLAGS) -c SOCKETS/sendrecv.c
soc_init.o : SOCKETS/soc_init.c SOCKETS/mysocks.h
	$(CC) $(CFLAGS) -c SOCKETS/soc_init.c

# *****************************************************************************
# thermodynamic integration objects
# *****************************************************************************
cclient.o : THERMOINT/cclient.c
	$(CC) $(CFLAGS) -c THERMOINT/cclient.c
	
# *****************************************************************************
# solvesh_diagonalization objects
# *****************************************************************************
blacsaba.o : SOLVESH_DIAG/blacsaba.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/blacsaba.f90
blas.o : SOLVESH_DIAG/blas.f # No MODULES
	$(F77) $(DFLAGS) -c SOLVESH_DIAG/blas.f
kspace.o : SOLVESH_DIAG/kspace.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/kspace.f90
kspace_withnok.o : SOLVESH_DIAG/kspace_withnok.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/kspace_withnok.f90
kspace_MPI.o : SOLVESH_DIAG/kspace_MPI.f90 $(MODULES)  
	$(F90) $(FFLAGS) $(MODMPI) -c SOLVESH_DIAG/kspace_MPI.f90
kspace_MPI_slave.o : SOLVESH_DIAG/kspace_MPI_slave.f90 $(MODULES) 
	$(F90) $(FFLAGS) $(MODMPI) -c SOLVESH_DIAG/kspace_MPI_slave.f90
kspace_ordern_fk.o : SOLVESH_DIAG/kspace_ordern_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/kspace_ordern_fk.f90
minilapack.o : SOLVESH_DIAG/minilapack.f # No MODULES
	$(F77) $(DFLAGS) -c SOLVESH_DIAG/minilapack.f
diag_error.o : SOLVESH_DIAG/diag_error.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/diag_error.f90
pclagetter.o : SOLVESH_DIAG/pclagetter.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/pclagetter.f90
pclaputter.o : SOLVESH_DIAG/pclaputter.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/pclaputter.f90
diag_k.o : SOLVESH_DIAG/diag_k.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/diag_k.f90
diag_k-MPI.o : SOLVESH_DIAG/diag_k-MPI.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/diag_k-MPI.f90
diag_k-MPI_slave.o : SOLVESH_DIAG/diag_k-MPI_slave.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/diag_k-MPI_slave.f90
bcast_k-MPI.o : SOLVESH_DIAG/bcast_k-MPI.f90 $(MODULES)
	$(F90) $(FFLAGS) $(MODMPI) -c SOLVESH_DIAG/bcast_k-MPI.f90
bcast_k-MPI_slave.o : SOLVESH_DIAG/bcast_k-MPI_slave.f90 $(MODULES)
	$(F90) $(FFLAGS) $(MODMPI) -c SOLVESH_DIAG/bcast_k-MPI_slave.f90
kspace_MPI-k.o : SOLVESH_DIAG/kspace_MPI-k.f90 $(MODULES)  
	$(F90) $(FFLAGS) $(MODMPI) -c SOLVESH_DIAG/kspace_MPI-k.f90
kspace_KS.o : SOLVESH_DIAG/kspace_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/kspace_KS.f90
diag_k_KS.o : SOLVESH_DIAG/diag_k_KS.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_DIAG/diag_k_KS.f90

# *****************************************************************************
# solvesh_ordern objects
# *****************************************************************************
eandg.o : SOLVESH_ORDERN/eandg.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/eandg.f90
xeandg.o : SOLVESH_ORDERN/xeandg.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/xeandg.f90
getstepsize.o : SOLVESH_ORDERN/getstepsize.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/getstepsize.f90
formc_compact.o : SOLVESH_ORDERN/formc_compact.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/formc_compact.f90
formsh_compact.o : SOLVESH_ORDERN/formsh_compact.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/formsh_compact.f90
initguess.o : SOLVESH_ORDERN/initguess.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/initguess.f90
kspace_fk.o : SOLVESH_ORDERN/kspace_fk.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/kspace_fk.f90
kspace_ordern.o : SOLVESH_ORDERN/kspace_ordern.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/kspace_ordern.f90
kspace_ordern_init.o : SOLVESH_ORDERN/kspace_ordern_init.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/kspace_ordern_init.f90
kspace_ordern_slave.o : SOLVESH_ORDERN/kspace_ordern_slave.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/kspace_ordern_slave.f90
getsendrecv.o: SOLVESH_ORDERN/getsendrecv.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/getsendrecv.f90
ordern_init.o : SOLVESH_ORDERN/ordern_init.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/ordern_init.f90
qralg.o : SOLVESH_ORDERN/qralg.f # No MODULES
	$(F77) $(DFLAGS) -c SOLVESH_ORDERN/qralg.f
set_dimensions.o : SOLVESH_ORDERN/set_dimensions.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/set_dimensions.f90
set_maxdimension.o : SOLVESH_ORDERN/set_maxdimension.f90 $(MODULES)
	$(F90) $(FFLAGS) -c SOLVESH_ORDERN/set_maxdimension.f90


# *****************************************************************************
# lapack_95 utilities objects
# *****************************************************************************
kspace_l95.o : SOLVESH_DIAG/kspace_l95.f90 $(MODULES)
	$(F90) $(FFLAGS) $(FLAG95) $(MODLAP95) -c SOLVESH_DIAG/kspace_l95.f90
anderson_l95.o : UTIL/anderson_l95.f90 $(MODULES)
	$(F90) $(FFLAGS) $(FLAG95) $(MODLAP95) -c UTIL/anderson_l95.f90

# *****************************************************************************
# double precision utilities objects
# *****************************************************************************
kspace2.o : SOLVESH_DIAG/kspace2.f90 $(MODULES)
	$(F90) $(FFLAGS) $(FLAG95) $(MODLAP95) -c SOLVESH_DIAG/kspace2.f90
anderson2.o : UTIL/anderson2.f90 $(MODULES)
	$(F90) $(FFLAGS) $(FLAG95) $(MODLAP95) -c UTIL/anderson2.f90
kspaceX.o : SOLVESH_DIAG/kspaceX.f90 $(MODULES)
	$(F90) $(FFLAGS) $(FLAG95) $(MODLAP95) -c SOLVESH_DIAG/kspaceX.f90

# *****************************************************************************
# sparse utilities objects
# *****************************************************************************
build_transpose.o : UTIL_SPARSE/build_transpose.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/build_transpose.f90
lanc.o : UTIL_SPARSE/lanc.f $(MODULES)
	$(F77) $(DFLAGS) -c UTIL_SPARSE/lanc.f
sparse_add.o : UTIL_SPARSE/sparse_add.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_add.f90
sparse_copy.o : UTIL_SPARSE/sparse_copy.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_copy.f90
sparse_getdimension.o : UTIL_SPARSE/sparse_getdimension.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_getdimension.f90
sparse_getpacksize.o : UTIL_SPARSE/sparse_getpacksize.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_getpacksize.f90
sparse_mult.o : UTIL_SPARSE/sparse_mult.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_mult.f90
sparse_mask.o : UTIL_SPARSE/sparse_mask.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_mask.f90
sparse_norm2.o : UTIL_SPARSE/sparse_norm2.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_norm2.f90
sparse_pack.o : UTIL_SPARSE/sparse_pack.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_pack.f90
sparse_pack_elements.o : UTIL_SPARSE/sparse_pack_elements.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_pack_elements.f90
sparse_pack_indices.o : UTIL_SPARSE/sparse_pack_indices.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_pack_indices.f90
sparse_unpack.o : UTIL_SPARSE/sparse_unpack.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_unpack.f90
sparse_unpack_elements.o : UTIL_SPARSE/sparse_unpack_elements.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_unpack_elements.f90
sparse_unpack_indices.o : UTIL_SPARSE/sparse_unpack_indices.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_unpack_indices.f90
sparse_vecprod.o : UTIL_SPARSE/sparse_vecprod.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL_SPARSE/sparse_vecprod.f90


# *****************************************************************************
# umbrella sampling objects
# *****************************************************************************
assemble_umbrella.o : UMBRELLA/assemble_umbrella.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/assemble_umbrella.f90
Dassemble_umbrella.o : UMBRELLA/Dassemble_umbrella.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/Dassemble_umbrella.f90
get_umbrella.o : UMBRELLA/get_umbrella.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/get_umbrella.f90
readumbrella.o : UMBRELLA/readumbrella.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/readumbrella.f90
assemble_steered.o : UMBRELLA/assemble_steered.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/assemble_steered.f90
Dassemble_steered.o : UMBRELLA/Dassemble_steered.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/Dassemble_steered.f90
get_steered.o : UMBRELLA/get_steered.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/get_steered.f90
readsteered.o : UMBRELLA/readsteered.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UMBRELLA/readsteered.f90

# *****************************************************************************
# utilities objects
# *****************************************************************************
anderson.o : UTIL/anderson.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/anderson.f90
fixfrags.o : UTIL/fixfrags.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/fixfrags.f90
fixfrags2.o : UTIL/fixfrags2.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/fixfrags2.f90
hampiece.o : UTIL/hampiece.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/hampiece.f90
push_atoms.o : UTIL/push_atoms.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/push_atoms.f90
writeout_ac.o : UTIL/writeout_ac.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_ac.f90
writeout_cd.o : UTIL/writeout_cd.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_cd.f90
writeout_charges.o : UTIL/writeout_charges.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_charges.f90
writeout_dipole.o : UTIL/writeout_dipole.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_dipole.f90
compute_neutral_rho_grid.o : UTIL/compute_neutral_rho_grid.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/compute_neutral_rho_grid.f90
writeout_comph.o : UTIL/writeout_comph.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_comph.f90
writeout_neighbors.o : UTIL/writeout_neighbors.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_neighbors.f90
writeout_xv.o : UTIL/writeout_xv.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_xv.f90
writeout_neighborsPP.o : UTIL/writeout_neighborsPP.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_neighborsPP.f90
hamtrans.o : UTIL/hamtrans.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/hamtrans.f90
mixer.o : UTIL/mixer.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/mixer.f90
den2mesh.o : UTIL/den2mesh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/den2mesh.f90
ew2mesh.o : UTIL/ew2mesh.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ew2mesh.f90
ew2mesh_gamma.o : UTIL/ew2mesh_gamma.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ew2mesh_gamma.f90
postscf.o : UTIL/postscf.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/postscf.f90
den2mesh_import.o : UTIL/den2mesh_import.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/den2mesh_import.f90	
pulay.o : UTIL/pulay.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/pulay.f90	
louie.o : UTIL/louie.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/louie.f90	
broyden.o : UTIL/broyden.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/broyden.f90	
# added by prokop
ew2mesh_ARPES.o : UTIL/ew2mesh_ARPES.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ew2mesh_ARPES.f90
ew2mesh_fourier.o : UTIL/ew2mesh_fourier.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ew2mesh_fourier.f90
ew2mesh_kscan.o : UTIL/ew2mesh_kscan.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ew2mesh_kscan.f90	
project_orb.o : UTIL/project_orb.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/project_orb.f90	
project_orb_complex.o : UTIL/project_orb_complex.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/project_orb_complex.f90
excitations.o : UTIL/excitations.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/excitations.f90	
kvaziband.o : UTIL/kvaziband.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/kvaziband.f90
ARPES_LCAO.o : UTIL/ARPES_LCAO.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/ARPES_LCAO.f90
band_project.o : UTIL/band_project.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/band_project.f90
writeCoefsLCAO.o : UTIL/writeCoefsLCAO.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeCoefsLCAO.f90
writeout_eigenvec.o : UTIL/writeout_eigenvec.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_eigenvec.f90
writeout_eigenvec_G.o : UTIL/writeout_eigenvec_G.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/writeout_eigenvec_G.f90
compute_all_charges.o : UTIL/compute_all_charges.f90 $(MODULES)
	$(F90) $(FFLAGS) -c UTIL/compute_all_charges.f90
# *****************************************************************************
# visualization objects
# *****************************************************************************
noTclMD.o : VISUALIZATION/noTclMD.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/noTclMD.f90
tclmdtransfer.o : VISUALIZATION/tclmdtransfer.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/tclmdtransfer.f90
graceinit.o : VISUALIZATION/graceinit.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/graceinit.f90
graceupdate.o : VISUALIZATION/graceupdate.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/graceupdate.f90
nograce.o : VISUALIZATION/nograce.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/nograce.f90
noxmgr.o : VISUALIZATION/noxmgr.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/noxmgr.f90
xmgrinit.o : VISUALIZATION/xmgrinit.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/xmgrinit.f90
xmgrupdate.o : VISUALIZATION/xmgrupdate.f90 $(MODULES)
	$(F90) $(FFLAGS) -c VISUALIZATION/xmgrupdate.f90


# *****************************************************************************
# xc objects
# *****************************************************************************
ceperley_alder.o : XC/ceperley_alder.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/ceperley_alder.f90
cepal.o : XC/cepal.f90 $(MODULES)
	$(F90) $(FFLAGS) -c XC/cepal.f90

# *****************************************************************************
# cg objects
# *****************************************************************************
cgo.o : CG/cgo.f90 $(MODULES)
	$(F90) $(FFLAGS) -c CG/cgo.f90
bfgs.o : CG/bfgs.f90 $(MODULES)
	$(F90) $(FFLAGS) -c CG/bfgs.f90
l-bfgs-b.o : CG/l-bfgs-b.f90 $(MODULES)
	$(F90) $(FFLAGS) -c CG/l-bfgs-b.f90
FIRE.o : CG/FIRE.f90  $(MODULES)
	$(F90) $(FFLAGS) -c CG/FIRE.f90

# *****************************************************************************
# NEB objects
# *****************************************************************************
initneb.o : NEB/initneb.f90
	$(F90) $(FFLAGS) -c NEB/initneb.f90
neb.o : NEB/neb.f90
	$(F90) $(FFLAGS) -c NEB/neb.f90

# *****************************************************************************
# TDSE objects
# *****************************************************************************
ete_loop.o : TDSE/ete_loop.f90
	$(F90) $(FFLAGS) -c TDSE/ete_loop.f90
diag_Sk.o : TDSE/diag_Sk.f90
	$(F90) $(FFLAGS) -c TDSE/diag_Sk.f90
diag_Hk.o : TDSE/diag_Hk.f90
	$(F90) $(FFLAGS) -c TDSE/diag_Hk.f90
tddiag_k.o : TDSE/tddiag_k.f90
	$(F90) $(FFLAGS) -c TDSE/tddiag_k.f90
eigenHS.o : TDSE/eigenHS.f90
	$(F90) $(FFLAGS) -c TDSE/eigenHS.f90
psi2es.o : TDSE/psi2es.f90
	$(F90) $(FFLAGS) -c TDSE/psi2es.f90
tddenmat.o : TDSE/tddenmat.f90
	$(F90) $(FFLAGS) -c TDSE/tddenmat.f90
allocate_tdse.o : TDSE/allocate_tdse.f90
	$(F90) $(FFLAGS) -c TDSE/allocate_tdse.f90
tdbc.o : TDSE/tdbc.f90
	$(F90) $(FFLAGS) -c TDSE/tdbc.f90
propTpsi.o : TDSE/propTpsi.f90
	$(F90) $(FFLAGS) -c TDSE/propTpsi.f90
readtdse.o : TDSE/readtdse.f90
	$(F90) $(FFLAGS) -c TDSE/readtdse.f90
initpsi.o : TDSE/initpsi.f90
	$(F90) $(FFLAGS) -c TDSE/initpsi.f90
ortho_H.o : TDSE/ortho_H.f90
	$(F90) $(FFLAGS) -c TDSE/ortho_H.f90
get_QLow.o : TDSE/get_QLow.f90
	$(F90) $(FFLAGS) -c TDSE/get_QLow.f90
get_QMul.o : TDSE/get_QMul.f90
	$(F90) $(FFLAGS) -c TDSE/get_QMul.f90
postete.o : TDSE/postete.f90
	$(F90) $(FFLAGS) -c TDSE/postete.f90
wrtout_psiT.o : TDSE/wrtout_psiT.f90
	$(F90) $(FFLAGS) -c TDSE/wrtout_psiT.f90

# *****************************************************************************
# TRANS objects
# *****************************************************************************
assemble_t12_fit.o : TRANS/assemble_t12_fit.f90
	$(F90) $(FFLAGS) -c TRANS/assemble_t12_fit.f90
assemble_t12_bare.o : TRANS/assemble_t12_bare.f90
	$(F90) $(FFLAGS) -c TRANS/assemble_t12_bare.f90
calcG.o : TRANS/calcG.f90
	$(F90) $(FFLAGS) -c TRANS/calcG.f90
assemble_Hsam.o : TRANS/assemble_Hsam.f90
	$(F90) $(FFLAGS) -c TRANS/assemble_Hsam.f90
assemble_Gsam.o : TRANS/assemble_Gsam.f90
	$(F90) $(FFLAGS) -c TRANS/assemble_Gsam.f90
assemble_Dxx.o : TRANS/assemble_Dxx.f90
	$(F90) $(FFLAGS) -c TRANS/assemble_Dxx.f90
sqrt_mat.o : TRANS/sqrt_mat.f90
	$(F90) $(FFLAGS) -c TRANS/sqrt_mat.f90
interpolate_hop.o : TRANS/interpolate_hop.f90
	$(F90) $(FFLAGS) -c TRANS/interpolate_hop.f90
gethop.o : TRANS/gethop.f90
	$(F90) $(FFLAGS) -c TRANS/gethop.f90
	
# *****************************************************************************
# BIAS objects
# *****************************************************************************
assemble_bias.o : BIAS/assemble_bias.f90
	$(F90) $(FFLAGS) -c BIAS/assemble_bias.f90
Dassemble_bias.o : BIAS/Dassemble_bias.f90
	$(F90) $(FFLAGS) -c BIAS/Dassemble_bias.f90
allocate_bias.o : BIAS/allocate_bias.f90
	$(F90) $(FFLAGS) -c BIAS/allocate_bias.f90
reallocate_bias.o : BIAS/reallocate_bias.f90
	$(F90) $(FFLAGS) -c BIAS/reallocate_bias.f90
readbias.o : BIAS/readbias.f90
	$(F90) $(FFLAGS) -c BIAS/readbias.f90
# *****************************************************************************
# NAC objects
# *****************************************************************************
allocate_nac.o : NAC/allocate_nac.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/allocate_nac.f90
deallocate_nac.o : NAC/deallocate_nac.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/deallocate_nac.f90
assemble_G_S.o : NAC/assemble_G_S.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/assemble_G_S.f90
nacouplings.o : NAC/nacouplings.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/nacouplings.f90
build_gover1c.o : NAC/build_gover1c.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/build_gover1c.f90
init_mdet.o : NAC/init_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/init_mdet.f90
mdetdenmat.o : NAC/mdetdenmat.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/mdetdenmat.f90
save_mdetstuff.o : NAC/save_mdetstuff.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/save_mdetstuff.f90
evolve_ks_states.o : NAC/evolve_ks_states.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/evolve_ks_states.f90
delta_t_ks.o : NAC/delta_t_ks.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/delta_t_ks.f90
dcdt_nac.o : NAC/dcdt_nac.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/dcdt_nac.f90
Dassemble_2c_mdet.o : NAC/Dassemble_2c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_2c_mdet.f90
Dassemble_ca_2c_mdet.o : NAC/Dassemble_ca_2c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_2c_mdet.f90
Dassemble_ca_2c_mdet_dip.o : NAC/Dassemble_ca_2c_mdet_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_2c_mdet_dip.f90
Dassemble_2c_PP_mdet.o : NAC/Dassemble_2c_PP_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_2c_PP_mdet.f90
Dassemble_olsxc_on_mdet.o : NAC/Dassemble_olsxc_on_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_olsxc_on_mdet.f90
Dassemble_ca_olsxc_on_mdet.o : NAC/Dassemble_ca_olsxc_on_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_olsxc_on_mdet.f90
Dassemble_olsxc_2c_mdet.o : NAC/Dassemble_olsxc_2c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_olsxc_2c_mdet.f90
Dassemble_ca_olsxc_2c_mdet.o : NAC/Dassemble_ca_olsxc_2c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_olsxc_2c_mdet.f90
Dassemble_3c_mdet.o : NAC/Dassemble_3c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_3c_mdet.f90
Dassemble_ca_3c_mdet.o : NAC/Dassemble_ca_3c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_3c_mdet.f90
Dassemble_ca_3c_mdet_dip.o : NAC/Dassemble_ca_3c_mdet_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_3c_mdet_dip.f90
Dassemble_lr_mdet.o : NAC/Dassemble_lr_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_lr_mdet.f90
Dassemble_lr_mdet_dip.o : NAC/Dassemble_lr_mdet_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_lr_mdet_dip.f90
Dassemble_3c_PP_mdet.o : NAC/Dassemble_3c_PP_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_3c_PP_mdet.f90
Dassemble_olsxc_3c_mdet.o : NAC/Dassemble_olsxc_3c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_olsxc_3c_mdet.f90
Dassemble_ca_olsxc_3c_mdet.o : NAC/Dassemble_ca_olsxc_3c_mdet.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/Dassemble_ca_olsxc_3c_mdet.f90
fewest_switches.o : NAC/fewest_switches.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/fewest_switches.f90
mc_switch.o : NAC/mc_switch.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/mc_switch.f90
transition.o : NAC/transition.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/transition.f90
move_correc.o : NAC/move_correc.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/move_correc.f90
move_predic.o : NAC/move_predic.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/move_predic.f90
overlap_sign.o : NAC/overlap_sign.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/overlap_sign.f90
check_swap.o : NAC/check_swap.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/check_swap.f90
overlap_numeric.o : NAC/overlap_numeric.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/overlap_numeric.f90
getnac.o : NAC/getnac.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/getnac.f90
MCsolar.o : NAC/MCsolar.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/MCsolar.f90
vibcouplings.o : NAC/vibcouplings.f90 $(MODULES)
	$(F90) $(FFLAGS) -c NAC/vibcouplings.f90
# *****************************************************************************
# qmmm
# *****************************************************************************
# qmm_module_null do nothing, only if you want libfireball.a use qmmm_module
qmmm_module_null.o : QMMM/qmmm_module_null.f90
	$(F90) $(FFLAGS) -c QMMM/qmmm_module_null.f90
qmmm_module.o: ../sqm/qmmm_module.F90
	cd ../sqm && $(MAKE) qmmm_module.o
	cp ../sqm/qmmm_module.mod .
	cp ../sqm/qmmm_module.o .
fireball_qmmm_loop.o : QMMM/fireball_qmmm_loop.f90
	$(F90) $(FFLAGS) -c QMMM/fireball_qmmm_loop.f90
assemble_qmmm.o : QMMM/assemble_qmmm.f90
	$(F90) $(FFLAGS) -c QMMM/assemble_qmmm.f90
assemble_qmmm_dip.o : QMMM/assemble_qmmm_dip.f90
	$(F90) $(FFLAGS) -c QMMM/assemble_qmmm_dip.f90
Dassemble_qmmm.o : QMMM/Dassemble_qmmm.f90
	$(F90) $(FFLAGS) -c QMMM/Dassemble_qmmm.f90
Dassemble_qmmm_dip.o : QMMM/Dassemble_qmmm_dip.f90
	$(F90) $(FFLAGS) -c QMMM/Dassemble_qmmm_dip.f90
Dassemble_qmmm_mdet.o : QMMM/Dassemble_qmmm_mdet.f90
	$(F90) $(FFLAGS) -c QMMM/Dassemble_qmmm_mdet.f90
Dassemble_qmmm_mdet_dip.o : QMMM/Dassemble_qmmm_mdet_dip.f90 $(MODULES)
	$(F90) $(FFLAGS) -c QMMM/Dassemble_qmmm_mdet_dip.f90
main_loop_MDET_qmmm.o : QMMM/main_loop_MDET_qmmm.f90
	$(F90) $(FFLAGS) -c QMMM/main_loop_MDET_qmmm.f90
main_loop_MD_qmmm.o : QMMM/main_loop_MD_qmmm.f90
	$(F90) $(FFLAGS) -c QMMM/main_loop_MD_qmmm.f90
# *****************************************************************************
# dftd3
# *****************************************************************************
api.o: dftd3/api.f90
	$(F90) $(FFLAGS) -c dftd3/api.f90
core.o: dftd3/core.f90
	$(F90) $(FFLAGS) -c dftd3/core.f90
sizes.o: dftd3/sizes.f90
	$(F90) $(FFLAGS) -c dftd3/sizes.f90
pars.o: dftd3/pars.f90
	$(F90) $(FFLAGS) -c dftd3/pars.f90
common.o: dftd3/common.f90
	$(F90) $(FFLAGS) -c dftd3/common.f90
dftd3_corrections.o : dftd3/dftd3_corrections.f90
	$(F90) $(FFLAGS) -c dftd3/dftd3_corrections.f90
