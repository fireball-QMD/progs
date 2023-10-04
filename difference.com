echo 'fireball.f90'
diff fireball.f90 ../progs.older/fireball.f90
echo 'Makefile'
diff Makefile ../progs.older/Makefile
echo 'ALLOCATIONS/allocate_f.f90'
diff ALLOCATIONS/allocate_f.f90 ../progs.older/ALLOCATIONS/allocate_f.f90
echo 'ALLOCATIONS/allocate_h.f90'
diff ALLOCATIONS/allocate_h.f90 ../progs.older/ALLOCATIONS/allocate_h.f90
echo 'ALLOCATIONS/allocate_neigh.f90'
diff ALLOCATIONS/allocate_neigh.f90 ../progs.older/ALLOCATIONS/allocate_neigh.f90
echo 'ALLOCATIONS/allocate_ordern.f90'
diff ALLOCATIONS/allocate_ordern.f90 ../progs.older/ALLOCATIONS/allocate_ordern.f90
echo 'ALLOCATIONS/allocate_rho.f90'
diff ALLOCATIONS/allocate_rho.f90 ../progs.older/ALLOCATIONS/allocate_rho.f90
echo 'ALLOCATIONS/allocate_umb.f90'
diff ALLOCATIONS/allocate_umb.f90 ../progs.older/ALLOCATIONS/allocate_umb.f90
echo 'ALLOCATIONS/reallocate_f.f90'
diff ALLOCATIONS/reallocate_f.f90 ../progs.older/ALLOCATIONS/reallocate_f.f90
echo 'ALLOCATIONS/reallocate_h.f90'
diff ALLOCATIONS/reallocate_h.f90 ../progs.older/ALLOCATIONS/reallocate_h.f90
echo 'ALLOCATIONS/reallocate_neigh.f90'
diff ALLOCATIONS/reallocate_neigh.f90 ../progs.older/ALLOCATIONS/reallocate_neigh.f90
echo 'ALLOCATIONS/reallocate_rho.f90'
diff ALLOCATIONS/reallocate_rho.f90 ../progs.older/ALLOCATIONS/reallocate_rho.f90
echo 'ASSEMBLERS/assemble_2c.f90'
diff ASSEMBLERS/assemble_2c.f90 ../progs.older/ASSEMBLERS/assemble_2c.f90
echo 'ASSEMBLERS/assemble_2c_PP.f90'
diff ASSEMBLERS/assemble_2c_PP.f90 ../progs.older/ASSEMBLERS/assemble_2c_PP.f90
echo 'ASSEMBLERS/assemble_2c_ordern_final.f90'
diff ASSEMBLERS/assemble_2c_ordern_final.f90 ../progs.older/ASSEMBLERS/assemble_2c_ordern_final.f90
echo 'ASSEMBLERS/assemble_2c_ordern_init.f90'
diff ASSEMBLERS/assemble_2c_ordern_init.f90 ../progs.older/ASSEMBLERS/assemble_2c_ordern_init.f90
echo 'ASSEMBLERS/assemble_3c.f90'
diff ASSEMBLERS/assemble_3c.f90 ../progs.older/ASSEMBLERS/assemble_3c.f90
echo 'ASSEMBLERS/assemble_3c_PP.f90'
diff ASSEMBLERS/assemble_3c_PP.f90 ../progs.older/ASSEMBLERS/assemble_3c_PP.f90
echo 'ASSEMBLERS/assemble_3c_ordern_final.f90'
diff ASSEMBLERS/assemble_3c_ordern_final.f90 ../progs.older/ASSEMBLERS/assemble_3c_ordern_final.f90
echo 'ASSEMBLERS/assemble_F.f90'
diff ASSEMBLERS/assemble_F.f90 ../progs.older/ASSEMBLERS/assemble_F.f90
echo 'ASSEMBLERS/assemble_ca_2c.f90'
diff ASSEMBLERS/assemble_ca_2c.f90 ../progs.older/ASSEMBLERS/assemble_ca_2c.f90
echo 'ASSEMBLERS/assemble_ca_2c_ordern_final.f90'
diff ASSEMBLERS/assemble_ca_2c_ordern_final.f90 ../progs.older/ASSEMBLERS/assemble_ca_2c_ordern_final.f90
echo 'ASSEMBLERS/assemble_ca_3c.f90'
diff ASSEMBLERS/assemble_ca_3c.f90 ../progs.older/ASSEMBLERS/assemble_ca_3c.f90
echo 'ASSEMBLERS/assemble_ca_3c_ordern_final.f90'
diff ASSEMBLERS/assemble_ca_3c_ordern_final.f90 ../progs.older/ASSEMBLERS/assemble_ca_3c_ordern_final.f90
echo 'ASSEMBLERS/assemble_eh_2c.f90'
diff ASSEMBLERS/assemble_eh_2c.f90 ../progs.older/ASSEMBLERS/assemble_eh_2c.f90
echo 'ASSEMBLERS/assemble_eh_usr.f90'
diff ASSEMBLERS/assemble_eh_usr.f90 ../progs.older/ASSEMBLERS/assemble_eh_usr.f90
echo 'ASSEMBLERS/assemble_hxc_2c.f90'
diff ASSEMBLERS/assemble_hxc_2c.f90 ../progs.older/ASSEMBLERS/assemble_hxc_2c.f90
echo 'ASSEMBLERS/assemble_hxc_3c.f90'
diff ASSEMBLERS/assemble_hxc_3c.f90 ../progs.older/ASSEMBLERS/assemble_hxc_3c.f90
echo 'ASSEMBLERS/assemble_hxc_usr.f90'
diff ASSEMBLERS/assemble_hxc_usr.f90 ../progs.older/ASSEMBLERS/assemble_hxc_usr.f90
echo 'ASSEMBLERS/assemble_lr.f90'
diff ASSEMBLERS/assemble_lr.f90 ../progs.older/ASSEMBLERS/assemble_lr.f90
echo 'ASSEMBLERS/assemble_lr_ordern_final.f90'
diff ASSEMBLERS/assemble_lr_ordern_final.f90 ../progs.older/ASSEMBLERS/assemble_lr_ordern_final.f90
echo 'ASSEMBLERS/assemble_olsxc_1c.f90'
diff ASSEMBLERS/assemble_olsxc_1c.f90 ../progs.older/ASSEMBLERS/assemble_olsxc_1c.f90
echo 'ASSEMBLERS/assemble_olsxc_on.f90'
diff ASSEMBLERS/assemble_olsxc_on.f90 ../progs.older/ASSEMBLERS/assemble_olsxc_on.f90
echo 'ASSEMBLERS/assemble_olsxc_off.f90'
diff ASSEMBLERS/assemble_olsxc_off.f90 ../progs.older/ASSEMBLERS/assemble_olsxc_off.f90
echo 'ASSEMBLERS/assemble_snxc_on.f90'
diff ASSEMBLERS/assemble_snxc_on.f90 ../progs.older/ASSEMBLERS/assemble_snxc_on.f90
echo 'ASSEMBLERS/assemble_snxc_off.f90'
diff ASSEMBLERS/assemble_snxc_off.f90 ../progs.older/ASSEMBLERS/assemble_snxc_off.f90
echo 'ASSEMBLERS/build_snxc_off.f90'
diff ASSEMBLERS/build_snxc_off.f90 ../progs.older/ASSEMBLERS/build_snxc_off.f90
echo 'ASSEMBLERS/build_snxc_on.f90'
diff ASSEMBLERS/build_snxc_on.f90 ../progs.older/ASSEMBLERS/build_snxc_on.f90
echo 'ASSEMBLERS/build_ca_snxc_on.f90'
diff ASSEMBLERS/build_ca_snxc_on.f90 ../progs.older/ASSEMBLERS/build_ca_snxc_on.f90
echo 'ASSEMBLERS/build_olsxc_off.f90'
diff ASSEMBLERS/build_olsxc_off.f90 ../progs.older/ASSEMBLERS/build_olsxc_off.f90
echo 'ASSEMBLERS/build_olsxc_on.f90'
diff ASSEMBLERS/build_olsxc_on.f90 ../progs.older/ASSEMBLERS/build_olsxc_on.f90
echo 'ASSEMBLERS/build_ca_olsxc_on.f90'
diff ASSEMBLERS/build_ca_olsxc_on.f90 ../progs.older/ASSEMBLERS/build_ca_olsxc_on.f90
echo 'ASSEMBLERS/average_rho.f90'
diff ASSEMBLERS/average_rho.f90 ../progs.older/ASSEMBLERS/average_rho.f90
echo 'ASSEMBLERS/average_ca_rho.f90'
diff ASSEMBLERS/average_ca_rho.f90 ../progs.older/ASSEMBLERS/average_ca_rho.f90

echo 'ASSEMBLERS/assemble_sVNL.f90'
diff ASSEMBLERS/assemble_sVNL.f90 ../progs.older/ASSEMBLERS/assemble_sVNL.f90
echo 'ASSEMBLERS/assemble_usr.f90'
diff ASSEMBLERS/assemble_usr.f90 ../progs.older/ASSEMBLERS/assemble_usr.f90
echo 'ASSEMBLERS/buildh.f90'
diff ASSEMBLERS/buildh.f90 ../progs.older/ASSEMBLERS/buildh.f90
echo 'CG/cgo.f90'
diff CG/cgo.f90 ../progs.older/CG/cgo.f90
echo 'DASSEMBLERS/Dassemble_2c.f90'
diff DASSEMBLERS/Dassemble_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_2c.f90
echo 'DASSEMBLERS/Dassemble_2c_ordern_final.f90'
diff DASSEMBLERS/Dassemble_2c_ordern_final.f90 ../progs.older/DASSEMBLERS/Dassemble_2c_ordern_final.f90
echo 'DASSEMBLERS/Dassemble_3c.f90'
diff DASSEMBLERS/Dassemble_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_3c.f90
echo 'DASSEMBLERS/Dassemble_ca_2c.f90'
diff DASSEMBLERS/Dassemble_ca_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_2c.f90
echo 'DASSEMBLERS/Dassemble_ca_2cordern_final.f90'
diff DASSEMBLERS/Dassemble_ca_2c_ordern_final.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_2c_ordern_final.f90
echo 'DASSEMBLERS/Dassemble_ca_3c.f90'
diff DASSEMBLERS/Dassemble_ca_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_3c.f90
echo 'DASSEMBLERS/Dassemble_ca_3c_ordern_final.f90'
diff DASSEMBLERS/Dassemble_ca_3c_ordern_final.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_3c_ordern_final.f90
echo 'DASSEMBLERS/Dassemble_eh_2c.f90'
diff DASSEMBLERS/Dassemble_eh_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_eh_2c.f90
echo 'DASSEMBLERS/Dassemble_hxc_2c.f90'
diff DASSEMBLERS/Dassemble_hxc_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_hxc_2c.f90
echo 'DASSEMBLERS/Dassemble_hxc_3c.f90'
diff DASSEMBLERS/Dassemble_hxc_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_hxc_3c.f90
echo 'DASSEMBLERS/Dassemble_lr.f90'
diff DASSEMBLERS/Dassemble_lr.f90 ../progs.older/DASSEMBLERS/Dassemble_lr.f90
echo 'DASSEMBLERS/Dassemble_olsxc_2c.f90'
diff DASSEMBLERS/Dassemble_olsxc_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_olsxc_2c.f90
echo 'DASSEMBLERS/Dassemble_olsxc_3c.f90'
diff DASSEMBLERS/Dassemble_olsxc_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_olsxc_3c.f90
echo 'DASSEMBLERS/Dassemble_olsxc_on.f90'
diff DASSEMBLERS/Dassemble_olsxc_on.f90 ../progs.older/DASSEMBLERS/Dassemble_olsxc_on.f90
echo 'DASSEMBLERS/Dassemble_snxc_2c.f90'
diff DASSEMBLERS/Dassemble_snxc_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_snxc_2c.f90
echo 'DASSEMBLERS/Dassemble_snxc_3c.f90'
diff DASSEMBLERS/Dassemble_snxc_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_snxc_3c.f90
echo 'DASSEMBLERS/Dassemble_snxc_on.f90'
diff DASSEMBLERS/Dassemble_snxc_on.f90 ../progs.older/DASSEMBLERS/Dassemble_snxc_on.f90
echo 'DASSEMBLERS/Dassemble_ca_olsxc_2c.f90'
diff DASSEMBLERS/Dassemble_ca_olsxc_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_olsxc_2c.f90
echo 'DASSEMBLERS/Dassemble_ca_olsxc_3c.f90'
diff DASSEMBLERS/Dassemble_ca_olsxc_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_olsxc_3c.f90
echo 'DASSEMBLERS/Dassemble_ca_olsxc_on.f90'
diff DASSEMBLERS/Dassemble_ca_olsxc_on.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_olsxc_on.f90
echo 'DASSEMBLERS/Dassemble_ca_snxc_2c.f90'
diff DASSEMBLERS/Dassemble_ca_snxc_2c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_snxc_2c.f90
echo 'DASSEMBLERS/Dassemble_ca_snxc_3c.f90'
diff DASSEMBLERS/Dassemble_ca_snxc_3c.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_snxc_3c.f90
echo 'DASSEMBLERS/Dassemble_ca_snxc_on.f90'
diff DASSEMBLERS/Dassemble_ca_snxc_on.f90 ../progs.older/DASSEMBLERS/Dassemble_ca_snxc_on.f90
echo 'DOS/dos.f90'
diff DOS/dos.f90 ../progs.older/DOS/dos.f90
echo 'DOS/hamilt_atom.f90'
diff DOS/hamilt_atom.f90 ../progs.older/DOS/hamilt_atom.f90
echo 'DOS/hoppings.f90'
diff DOS/hoppings.f90 ../progs.older/DOS/hoppings.f90
echo 'DOS/invierte.f'
diff DOS/invierte.f ../progs.older/DOS/invierte.f
echo 'DOS/writeout_atom.f90'
diff DOS/writeout_atom.f90 ../progs.older/DOS/writeout_atom.f90
echo 'DOS/writeout_dos.f90'
diff DOS/writeout_dos.f90 ../progs.older/DOS/writeout_dos.f90
echo 'FORM_RHO/chebft.f90'
diff FORM_RHO/chebft.f90 ../progs.older/FORM_RHO/chebft.f90
echo 'FORM_RHO/denmat.f90'
diff FORM_RHO/denmat.f90 ../progs.older/FORM_RHO/denmat.f90
echo 'FORM_RHO/denmat_fk.f90'
diff FORM_RHO/denmat_fk.f90 ../progs.older/FORM_RHO/denmat_fk.f90
echo 'FORM_RHO/denmata_ordern.f90'
diff FORM_RHO/denmata_ordern.f90 ../progs.older/FORM_RHO/denmata_ordern.f90
echo 'FORM_RHO/denmata_ordern_fk.f90'
diff FORM_RHO/denmata_ordern_fk.f90 ../progs.older/FORM_RHO/denmata_ordern_fk.f90
echo 'FORM_RHO/denmatb_ordern.f90'
diff FORM_RHO/denmatb_ordern.f90 ../progs.older/FORM_RHO/denmatb_ordern.f90
echo 'FORM_RHO/denmatb_ordern_fk.f90'
diff FORM_RHO/denmatb_ordern_fk.f90 ../progs.older/FORM_RHO/denmatb_ordern_fk.f90
echo 'FORM_RHO/denmatc_ordern.f90'
diff FORM_RHO/denmatc_ordern.f90 ../progs.older/FORM_RHO/denmatc_ordern.f90
echo 'FORM_RHO/denmatc_ordern_fk.f90'
diff FORM_RHO/denmatc_ordern_fk.f90 ../progs.older/FORM_RHO/denmatc_ordern_fk.f90
echo 'FORM_RHO/fermie.f90'
diff FORM_RHO/fermie.f90 ../progs.older/FORM_RHO/fermie.f90
echo 'FORM_RHO/formrho_sparse.f90'
diff FORM_RHO/formrho_sparse.f90 ../progs.older/FORM_RHO/formrho_sparse.f90
echo 'FORM_RHO/ss12.f90'
diff FORM_RHO/ss12.f90 ../progs.older/FORM_RHO/ss12.f90
echo 'FORM_RHO/ss12_fk.f90'
diff FORM_RHO/ss12_fk.f90 ../progs.older/FORM_RHO/ss12_fk.f90
diff FUTURE/experiments.f90 ../progs.older/FUTURE/experiments.f90
diff FUTURE/sa.f90 ../progs.older/FUTURE/sa.f90
diff FUTURE/template.f90 ../progs.older/FUTURE/template.f90
echo 'INITIALIZERS/diagnostics.f90'
diff INITIALIZERS/diagnostics.f90 ../progs.older/INITIALIZERS/diagnostics.f90
echo 'INITIALIZERS/init_MPI.f90'
diff INITIALIZERS/init_MPI.f90 ../progs.older/INITIALIZERS/init_MPI.f90
echo 'INITIALIZERS/init_noMPI.f90'
diff INITIALIZERS/init_noMPI.f90 ../progs.older/INITIALIZERS/init_noMPI.f90
echo 'INITIALIZERS/initatomicE.f90'
diff INITIALIZERS/initatomicE.f90 ../progs.older/INITIALIZERS/initatomicE.f90
echo 'INITIALIZERS/initboxes.f90'
diff INITIALIZERS/initboxes.f90 ../progs.older/INITIALIZERS/initboxes.f90
echo 'INITIALIZERS/initcharges.f90'
diff INITIALIZERS/initcharges.f90 ../progs.older/INITIALIZERS/initcharges.f90
echo 'INITIALIZERS/initconstants.f90'
diff INITIALIZERS/initconstants.f90 ../progs.older/INITIALIZERS/initconstants.f90
echo 'INITIALIZERS/initconstraints.f90'
diff INITIALIZERS/initconstraints.f90 ../progs.older/INITIALIZERS/initconstraints.f90
echo 'INITIALIZERS/initkpoints.f90'
diff INITIALIZERS/initkpoints.f90 ../progs.older/INITIALIZERS/initkpoints.f90
echo 'INITIALIZERS/initmasses.f90'
diff INITIALIZERS/initmasses.f90 ../progs.older/INITIALIZERS/initmasses.f90
echo 'INITIALIZERS/initneighbors.f90'
diff INITIALIZERS/initneighbors.f90 ../progs.older/INITIALIZERS/initneighbors.f90
echo 'INITIALIZERS/make_mu2shell.f90'
diff INITIALIZERS/make_mu2shell.f90 ../progs.older/INITIALIZERS/make_mu2shell.f90
echo 'INITIALIZERS/make_munu.f90'
diff INITIALIZERS/make_munu.f90 ../progs.older/INITIALIZERS/make_munu.f90
echo 'INITIALIZERS/make_munuPP.f90'
diff INITIALIZERS/make_munuPP.f90 ../progs.older/INITIALIZERS/make_munuPP.f90
echo 'INITIALIZERS/restart.f90'
diff INITIALIZERS/restart.f90 ../progs.older/INITIALIZERS/restart.f90
echo 'INITIALIZERS/welcome.f90'
diff INITIALIZERS/welcome.f90 ../progs.older/INITIALIZERS/welcome.f90
echo 'INITIALIZERS/zero_ang_mom.f90'
diff INITIALIZERS/zero_ang_mom.f90 ../progs.older/INITIALIZERS/zero_ang_mom.f90
echo 'INTERACTIONS/Dtrescentros.f90'
diff INTERACTIONS/Dtrescentros.f90 ../progs.older/INTERACTIONS/Dtrescentros.f90
echo 'INTERACTIONS/DtrescentrosS.f90'
diff INTERACTIONS/DtrescentrosS.f90 ../progs.older/INTERACTIONS/DtrescentrosS.f90
echo 'INTERACTIONS/DtrescentrosGHXC_VXC.f90'
diff INTERACTIONS/DtrescentrosGHXC_VXC.f90 ../progs.older/INTERACTIONS/DtrescentrosGHXC_VXC.f90
echo 'INTERACTIONS/DtrescentrosGS_VXC.f90'
diff INTERACTIONS/DtrescentrosGS_VXC.f90 ../progs.older/INTERACTIONS/DtrescentrosGS_VXC.f90
echo 'INTERACTIONS/DtrescentrosG_VXC.f90'
diff INTERACTIONS/DtrescentrosG_VXC.f90 ../progs.older/INTERACTIONS/DtrescentrosG_VXC.f90
echo 'INTERACTIONS/DtrescentrosG_VNA.f90'
diff INTERACTIONS/DtrescentrosG_VNA.f90 ../progs.older/INTERACTIONS/DtrescentrosG_VNA.f90
echo 'INTERACTIONS/DtrescentrosG_VNA_SH.f90'
diff INTERACTIONS/DtrescentrosG_VNA_SH.f90 ../progs.older/INTERACTIONS/DtrescentrosG_VNA_SH.f90
echo 'INTERACTIONS/DgelementsG_overlap.f90'
diff INTERACTIONS/DgelementsG_overlap.f90 ../progs.older/INTERACTIONS/DgelementsG_overlap.f90
echo 'INTERACTIONS/DgelementsGS_overlap.f90'
diff INTERACTIONS/DgelementsGS_overlap.f90 ../progs.older/INTERACTIONS/DgelementsGS_overlap.f90
echo 'INTERACTIONS/Dgelements_VXC.f90'
diff INTERACTIONS/Dgelements_VXC.f90 ../progs.older/INTERACTIONS/Dgelements_VXC.f90
echo 'INTERACTIONS/DgelementsG_VXC.f90'
diff INTERACTIONS/DgelementsG_VXC.f90 ../progs.older/INTERACTIONS/DgelementsG_VXC.f90
echo 'INTERACTIONS/DgelementsGS_VXC.f90'
diff INTERACTIONS/DgelementsGS_VXC.f90 ../progs.older/INTERACTIONS/DgelementsGS_VXC.f90
echo 'INTERACTIONS/DgelementsG_VNA.f90'
diff INTERACTIONS/DgelementsG_VNA.f90 ../progs.older/INTERACTIONS/DgelementsG_VNA.f90
echo 'INTERACTIONS/DgelementsG_VNA_SH.f90'
diff INTERACTIONS/DgelementsG_VNA_SH.f90 ../progs.older/INTERACTIONS/DgelementsG_VNA_SH.f90
echo 'INTERACTIONS/DgelementsG_VNA_SH_RNA.f90'
diff INTERACTIONS/DgelementsG_VNA_SH_RNA.f90 ../progs.older/INTERACTIONS/DgelementsG_VNA_SH_RNA.f90
echo 'INTERACTIONS/cl_value.f90'
diff INTERACTIONS/cl_value.f90 ../progs.older/INTERACTIONS/cl_value.f90
echo 'INTERACTIONS/dosgaussians.f90'
diff INTERACTIONS/dosgaussians.f90 ../progs.older/INTERACTIONS/dosgaussians.f90
echo 'INTERACTIONS/doscentros.f90'
diff INTERACTIONS/doscentros.f90 ../progs.older/INTERACTIONS/doscentros.f90
echo 'INTERACTIONS/doscentrosPP.f90'
diff INTERACTIONS/doscentrosPP.f90 ../progs.older/INTERACTIONS/doscentrosPP.f90
echo 'INTERACTIONS/gelementsG_overlap.f90'
diff INTERACTIONS/gelementsG_overlap.f90 ../progs.older/INTERACTIONS/gelementsG_overlap.f90
echo 'INTERACTIONS/gelementsGS_overlap.f90'
diff INTERACTIONS/gelementsGS_overlap.f90 ../progs.older/INTERACTIONS/gelementsGS_overlap.f90
echo 'INTERACTIONS/gelements_VXC.f90'
diff INTERACTIONS/gelements_VXC.f90 ../progs.older/INTERACTIONS/gelements_VXC.f90
echo 'INTERACTIONS/gelementsG_VXC.f90'
diff INTERACTIONS/gelementsG_VXC.f90 ../progs.older/INTERACTIONS/gelementsG_VXC.f90
echo 'INTERACTIONS/gelementsGS_VXC.f90'
diff INTERACTIONS/gelementsGS_VXC.f90 ../progs.older/INTERACTIONS/gelementsGS_VXC.f90
echo 'INTERACTIONS/gelementsG_VNA.f90'
diff INTERACTIONS/gelementsG_VNA.f90 ../progs.older/INTERACTIONS/gelementsG_VNA.f90
echo 'INTERACTIONS/gelementsG_VNA_SH.f90'
diff INTERACTIONS/gelementsG_VNA_SH.f90 ../progs.older/INTERACTIONS/gelementsG_VNA_SH.f90
echo 'INTERACTIONS/gelementsG_VNA_SH_RNA.f90'
diff INTERACTIONS/gelementsG_VNA_SH_RNA.f90 ../progs.older/INTERACTIONS/gelementsG_VNA_SH_RNA.f90
echo 'INTERACTIONS/get_ewald.f90'
diff INTERACTIONS/get_ewald.f90 ../progs.older/INTERACTIONS/get_ewald.f90
echo 'INTERACTIONS/get_vdw.f90'
diff INTERACTIONS/get_vdw.f90 ../progs.older/INTERACTIONS/get_vdw.f90
echo 'INTERACTIONS/smoother.f90'
diff INTERACTIONS/smoother.f90 ../progs.older/INTERACTIONS/smoother.f90
echo 'INTERACTIONS/trescentros.f90'
diff INTERACTIONS/trescentros.f90 ../progs.older/INTERACTIONS/trescentros.f90
echo 'INTERACTIONS/trescentrosS.f90'
diff INTERACTIONS/trescentrosS.f90 ../progs.older/INTERACTIONS/trescentrosS.f90
echo 'INTERACTIONS/trescentrosGHXC_VXC.f90'
diff INTERACTIONS/trescentrosGHXC_VXC.f90 ../progs.older/INTERACTIONS/trescentrosGHXC_VXC.f90
echo 'INTERACTIONS/trescentrosGS_VXC.f90'
diff INTERACTIONS/trescentrosGS_VXC.f90 ../progs.older/INTERACTIONS/trescentrosGS_VXC.f90
echo 'INTERACTIONS/trescentrosG_VXC.f90'
diff INTERACTIONS/trescentrosG_VXC.f90 ../progs.older/INTERACTIONS/trescentrosG_VXC.f90
echo 'INTERACTIONS/trescentrosG_VNA.f90'
diff INTERACTIONS/trescentrosG_VNA.f90 ../progs.older/INTERACTIONS/trescentrosG_VNA.f90
echo 'INTERACTIONS/trescentrosG_VNA_SH.f90'
diff INTERACTIONS/trescentrosG_VNA_SH.f90 ../progs.older/INTERACTIONS/trescentrosG_VNA_SH.f90
echo 'INTERACTIONS/unocentros.f90'
diff INTERACTIONS/unocentros.f90 ../progs.older/INTERACTIONS/unocentros.f90
echo 'INTERPOLATERS/buildspline_1d.f90'
diff INTERPOLATERS/buildspline_1d.f90 ../progs.older/INTERPOLATERS/buildspline_1d.f90
echo 'INTERPOLATERS/interpolate_1d.f90'
diff INTERPOLATERS/interpolate_1d.f90 ../progs.older/INTERPOLATERS/interpolate_1d.f90
echo 'INTERPOLATERS/interpolate_2d.f90'
diff INTERPOLATERS/interpolate_2d.f90 ../progs.older/INTERPOLATERS/interpolate_2d.f90
echo 'INTERPOLATERS/recoverC.f90'
diff INTERPOLATERS/recoverC.f90 ../progs.older/INTERPOLATERS/recoverC.f90
echo 'INTERPOLATERS/recover_2c.f90'
diff INTERPOLATERS/recover_2c.f90 ../progs.older/INTERPOLATERS/recover_2c.f90
echo 'INTERPOLATERS/recover_3c.f90'
diff INTERPOLATERS/recover_3c.f90 ../progs.older/INTERPOLATERS/recover_3c.f90
echo 'INTERPOLATERS/recover_PP.f90'
diff INTERPOLATERS/recover_PP.f90 ../progs.older/INTERPOLATERS/recover_PP.f90
echo 'INTERPOLATERS/setterp_2d.f90'
diff INTERPOLATERS/setterp_2d.f90 ../progs.older/INTERPOLATERS/setterp_2d.f90
echo 'MD/bvec.f90'
diff MD/bvec.f90 ../progs.older/MD/bvec.f90
echo 'MD/corrector.f90'
diff MD/corrector.f90 ../progs.older/MD/corrector.f90
echo 'MD/cross.f90'
diff MD/cross.f90 ../progs.older/MD/cross.f90
echo 'MD/factorial.f90'
diff MD/factorial.f90 ../progs.older/MD/factorial.f90
echo 'MD/gaussT.f90'
diff MD/gaussT.f90 ../progs.older/MD/gaussT.f90
echo 'MD/phimat.f90'
diff MD/phimat.f90 ../progs.older/MD/phimat.f90
echo 'MD/predictor.f90'
diff MD/predictor.f90 ../progs.older/MD/predictor.f90
echo 'MD/setgear.f90'
diff MD/setgear.f90 ../progs.older/MD/setgear.f90
echo 'MD/soldm.f90'
diff MD/soldm.f90 ../progs.older/MD/soldm.f90
echo 'MODULES/barrier.f90'
diff MODULES/barrier.f90 ../progs.older/MODULES/barrier.f90
echo 'MODULES/cgoptim.f90'
diff MODULES/cgoptim.f90 ../progs.older/MODULES/cgoptim.f90
echo 'MODULES/charges.f90'
diff MODULES/charges.f90 ../progs.older/MODULES/charges.f90
echo 'MODULES/configuration.f90'
diff MODULES/configuration.f90 ../progs.older/MODULES/configuration.f90
echo 'MODULES/constants_fireball.f90'
diff MODULES/constants_fireball.f90 ../progs.older/MODULES/constants_fireball.f90
echo 'MODULES/density.f90'
diff MODULES/density.f90 ../progs.older/MODULES/density.f90
echo 'MODULES/dimensions.f90'
diff MODULES/dimensions.f90 ../progs.older/MODULES/dimensions.f90
echo 'MODULES/dynamo.f90'
diff MODULES/dynamo.f90 ../progs.older/MODULES/dynamo.f90
echo 'MODULES/forces.f90'
diff MODULES/forces.f90 ../progs.older/MODULES/forces.f90
echo 'MODULES/fragments.f90'
diff MODULES/fragments.f90 ../progs.older/MODULES/fragments.f90
echo 'MODULES/gaussG.f90'
diff MODULES/gaussG.f90 ../progs.older/MODULES/gaussG.f90
echo 'MODULES/integrals.f90'
diff MODULES/integrals.f90 ../progs.older/MODULES/integrals.f90
echo 'MODULES/interactions.f90'
diff MODULES/interactions.f90 ../progs.older/MODULES/interactions.f90
echo 'MODULES/kpoints.f90'
diff MODULES/kpoints.f90 ../progs.older/MODULES/kpoints.f90
echo 'MODULES/mpi_declarations.f90'
diff MODULES/mpi_declarations.f90 ../progs.older/MODULES/mpi_declarations.f90
echo 'MODULES/neighbor_map.f90'
diff MODULES/neighbor_map.f90 ../progs.older/MODULES/neighbor_map.f90
echo 'MODULES/ordern.f90'
diff MODULES/ordern.f90 ../progs.older/MODULES/ordern.f90
echo 'MODULES/scf.f90'
diff MODULES/scf.f90 ../progs.older/MODULES/scf.f90
echo 'MODULES/umbrella.f90'
diff MODULES/umbrella.f90 ../progs.older/MODULES/umbrella.f90
echo 'NEIGHBORS/backnay.f90'
diff NEIGHBORS/backnay.f90 ../progs.older/NEIGHBORS/backnay.f90
echo 'NEIGHBORS/common_neighbors.f90'
diff NEIGHBORS/common_neighbors.f90 ../progs.older/NEIGHBORS/common_neighbors.f90
echo 'NEIGHBORS/common_neighborsPP.f90'
diff NEIGHBORS/common_neighborsPP.f90 ../progs.older/NEIGHBORS/common_neighborsPP.f90
echo 'NEIGHBORS/find_neigh_max.f90'
diff NEIGHBORS/find_neigh_max.f90 ../progs.older/NEIGHBORS/find_neigh_max.f90
echo 'NEIGHBORS/find_neighPP_max.f90'
diff NEIGHBORS/find_neighPP_max.f90 ../progs.older/NEIGHBORS/find_neighPP_max.f90
echo 'NEIGHBORS/mpairnay.f90'
diff NEIGHBORS/mpairnay.f90 ../progs.older/NEIGHBORS/mpairnay.f90
echo 'NEIGHBORS/neighbors.f90'
diff NEIGHBORS/neighbors.f90 ../progs.older/NEIGHBORS/neighbors.f90
echo 'NEIGHBORS/neighborsPP.f90'
diff NEIGHBORS/neighborsPP.f90 ../progs.older/NEIGHBORS/neighborsPP.f90
echo 'NEIGHBORS/num_neigh_tot.f90'
diff NEIGHBORS/num_neigh_tot.f90 ../progs.older/NEIGHBORS/num_neigh_tot.f90
diff PRESSURE/hmetric.f90 ../progs.older/PRESSURE/hmetric.f90
diff PRESSURE/initpressure.f90 ../progs.older/PRESSURE/initpressure.f90
diff PRESSURE/invert3x3.f90 ../progs.older/PRESSURE/invert3x3.f90
echo 'READFILES/append_string.f90'
diff READFILES/append_string.f90 ../progs.older/READFILES/append_string.f90
echo 'READFILES/read_1c.f90'
diff READFILES/read_1c.f90 ../progs.older/READFILES/read_1c.f90
echo 'READFILES/read_2c.f90'
diff READFILES/read_2c.f90 ../progs.older/READFILES/read_2c.f90
echo 'READFILES/read_3c.f90'
diff READFILES/read_3c.f90 ../progs.older/READFILES/read_3c.f90
echo 'READFILES/readbarrier.f90'
diff READFILES/readbarrier.f90 ../progs.older/READFILES/readbarrier.f90
echo 'READFILES/readbasis.f90'
diff READFILES/readbasis.f90 ../progs.older/READFILES/readbasis.f90
echo 'READFILES/readdata_2c.f90'
diff READFILES/readdata_2c.f90 ../progs.older/READFILES/readdata_2c.f90
echo 'READFILES/readdata_3c.f90'
diff READFILES/readdata_3c.f90 ../progs.older/READFILES/readdata_3c.f90
echo 'READFILES/readgaussG.f90'
diff READFILES/readgaussG.f90 ../progs.older/READFILES/readgaussG.f90
echo 'READFILES/readheader_2c.f90'
diff READFILES/readheader_2c.f90 ../progs.older/READFILES/readheader_2c.f90
echo 'READFILES/readheader_3c.f90'
diff READFILES/readheader_3c.f90 ../progs.older/READFILES/readheader_3c.f90
echo 'READFILES/readinfo.f90'
diff READFILES/readinfo.f90 ../progs.older/READFILES/readinfo.f90
echo 'READFILES/readkpoints.f90'
diff READFILES/readkpoints.f90 ../progs.older/READFILES/readkpoints.f90
echo 'READFILES/readlvs.f90'
diff READFILES/readlvs.f90 ../progs.older/READFILES/readlvs.f90
echo 'READFILES/readoptions.f90'
diff READFILES/readoptions.f90 ../progs.older/READFILES/readoptions.f90
echo 'READFILES/readoutput.f90'
diff READFILES/readoutput.f90 ../progs.older/READFILES/readoutput.f90
echo 'READFILES/readphi.f90'
diff READFILES/readphi.f90 ../progs.older/READFILES/readphi.f90
echo 'READFILES/readpressure.f90'
diff READFILES/readpressure.f90 ../progs.older/READFILES/readpressure.f90
echo 'READFILES/readquench.f90'
diff READFILES/readquench.f90 ../progs.older/READFILES/readquench.f90
echo 'READFILES/readsa.f90'
diff READFILES/readsa.f90 ../progs.older/READFILES/readsa.f90
echo 'READFILES/readscf.f90'
diff READFILES/readscf.f90 ../progs.older/READFILES/readscf.f90
echo 'READFILES/readscript.f90'
diff READFILES/readscript.f90 ../progs.older/READFILES/readscript.f90
echo 'READFILES/readvdw.f90'
diff READFILES/readvdw.f90 ../progs.older/READFILES/readvdw.f90
echo 'READFILES/readdos.f90'
diff READFILES/readcgo.f90 ../progs.older/READFILES/readcgo.f90
echo 'READFILES/readcgo.f90'
diff READFILES/readdos.f90 ../progs.older/READFILES/readdos.f90
diff ROTATIONS/chooser.f90 ../progs.older/ROTATIONS/chooser.f90
diff ROTATIONS/chooserd.f90 ../progs.older/ROTATIONS/chooserd.f90
diff ROTATIONS/deps2center.f90 ../progs.older/ROTATIONS/deps2center.f90
diff ROTATIONS/deps3center.f90 ../progs.older/ROTATIONS/deps3center.f90
diff ROTATIONS/epsilon.f90 ../progs.older/ROTATIONS/epsilon.f90
diff ROTATIONS/makeDmat.f90 ../progs.older/ROTATIONS/makeDmat.f90
diff ROTATIONS/makeDmatPP.f90 ../progs.older/ROTATIONS/makeDmatPP.f90
diff ROTATIONS/rotate.f90 ../progs.older/ROTATIONS/rotate.f90
diff ROTATIONS/rotatePP.f90 ../progs.older/ROTATIONS/rotatePP.f90
diff ROTATIONS/rotated.f90 ../progs.older/ROTATIONS/rotated.f90
diff ROTATIONS/rotatedPP.f90 ../progs.older/ROTATIONS/rotatedPP.f90
diff ROTATIONS/twister.f90 ../progs.older/ROTATIONS/twister.f90
diff ROTATIONS/twisterd.f90 ../progs.older/ROTATIONS/twisterd.f90
echo 'SOLVESH_DIAG/blacsaba.f90'
diff SOLVESH_DIAG/blacsaba.f90 ../progs.older/SOLVESH_DIAG/blacsaba.f90
echo 'SOLVESH_DIAG/diag_error.f90'
diff SOLVESH_DIAG/diag_error.f90 ../progs.older/SOLVESH_DIAG/diag_error.f90
echo 'SOLVESH_DIAG/kspace.f90'
diff SOLVESH_DIAG/kspace.f90 ../progs.older/SOLVESH_DIAG/kspace.f90
echo 'SOLVESH_DIAG/kspace2.f90'
diff SOLVESH_DIAG/kspace2.f90 ../progs.older/SOLVESH_DIAG/kspace2.f90
echo 'SOLVESH_DIAG/kspace_l95.f90'
diff SOLVESH_DIAG/kspace_l95.f90 ../progs.older/SOLVESH_DIAG/kspace_l95.f90
echo 'SOLVESH_DIAG/kspace_withnok.f90'
diff SOLVESH_DIAG/kspace_withnok.f90 ../progs.older/SOLVESH_DIAG/kspace_withnok.f90
echo 'SOLVESH_DIAG/kspace_MPI.f90'
diff SOLVESH_DIAG/kspace_MPI.f90 ../progs.older/SOLVESH_DIAG/kspace_MPI.f90
echo 'SOLVESH_DIAG/kspace_MPI_slave.f90'
diff SOLVESH_DIAG/kspace_MPI_slave.f90 ../progs.older/SOLVESH_DIAG/kspace_MPI_slave.f90
echo 'SOLVESH_DIAG/kspace_ordern_fk.f90'
diff SOLVESH_DIAG/kspace_ordern_fk.f90 ../progs.older/SOLVESH_DIAG/kspace_ordern_fk.f90
echo 'SOLVESH_DIAG/pclagetter.f90'
diff SOLVESH_DIAG/pclagetter.f90 ../progs.older/SOLVESH_DIAG/pclagetter.f90
echo 'SOLVESH_DIAG/pclaputter.f90'
diff SOLVESH_DIAG/pclaputter.f90 ../progs.older/SOLVESH_DIAG/pclaputter.f90
echo 'SOLVESH_ORDERN/eandg.f90'
diff SOLVESH_ORDERN/eandg.f90 ../progs.older/SOLVESH_ORDERN/eandg.f90
echo 'SOLVESH_ORDERN/formc_compact.f90'
diff SOLVESH_ORDERN/formc_compact.f90 ../progs.older/SOLVESH_ORDERN/formc_compact.f90
echo 'SOLVESH_ORDERN/formsh_compact.f90'
diff SOLVESH_ORDERN/formsh_compact.f90 ../progs.older/SOLVESH_ORDERN/formsh_compact.f90
echo 'SOLVESH_ORDERN/getsendrecv.f90'
diff SOLVESH_ORDERN/getsendrecv.f90 ../progs.older/SOLVESH_ORDERN/getsendrecv.f90
echo 'SOLVESH_ORDERN/getstepsize.f90'
diff SOLVESH_ORDERN/getstepsize.f90 ../progs.older/SOLVESH_ORDERN/getstepsize.f90
echo 'SOLVESH_ORDERN/initguess.f90'
diff SOLVESH_ORDERN/initguess.f90 ../progs.older/SOLVESH_ORDERN/initguess.f90
echo 'SOLVESH_ORDERN/kspace_fk.f90'
diff SOLVESH_ORDERN/kspace_fk.f90 ../progs.older/SOLVESH_ORDERN/kspace_fk.f90
echo 'SOLVESH_ORDERN/kspace_ordern.f90'
diff SOLVESH_ORDERN/kspace_ordern.f90 ../progs.older/SOLVESH_ORDERN/kspace_ordern.f90
echo 'SOLVESH_ORDERN/kspace_ordern_init.f90'
diff SOLVESH_ORDERN/kspace_ordern_init.f90 ../progs.older/SOLVESH_ORDERN/kspace_ordern_init.f90
echo 'SOLVESH_ORDERN/kspace_ordern_slave.f90'
diff SOLVESH_ORDERN/kspace_ordern_slave.f90 ../progs.older/SOLVESH_ORDERN/kspace_ordern_slave.f90
echo 'SOLVESH_ORDERN/ordern_init.f90'
diff SOLVESH_ORDERN/ordern_init.f90 ../progs.older/SOLVESH_ORDERN/ordern_init.f90
echo 'SOLVESH_ORDERN/set_dimensions.f90'
diff SOLVESH_ORDERN/set_dimensions.f90 ../progs.older/SOLVESH_ORDERN/set_dimensions.f90
echo 'SOLVESH_ORDERN/set_maxdimension.f90'
diff SOLVESH_ORDERN/set_maxdimension.f90 ../progs.older/SOLVESH_ORDERN/set_maxdimension.f90
echo 'SOLVESH_ORDERN/xeandg.f90'
diff SOLVESH_ORDERN/xeandg.f90 ../progs.older/SOLVESH_ORDERN/xeandg.f90
echo 'UMBRELLA/Dassemble_umbrella.f90'
diff UMBRELLA/Dassemble_umbrella.f90 ../progs.older/UMBRELLA/Dassemble_umbrella.f90
echo 'UMBRELLA/assemble_umbrella.f90'
diff UMBRELLA/assemble_umbrella.f90 ../progs.older/UMBRELLA/assemble_umbrella.f90
echo 'UMBRELLA/get_umbrella.f90'
diff UMBRELLA/get_umbrella.f90 ../progs.older/UMBRELLA/get_umbrella.f90
echo 'UMBRELLA/readumbrella.f90'
diff UMBRELLA/readumbrella.f90 ../progs.older/UMBRELLA/readumbrella.f90
echo 'UTIL/anderson2.f90'
diff UTIL/anderson2.f90 ../progs.older/UTIL/anderson2.f90
echo 'UTIL/anderson77.f90'
diff UTIL/anderson77.f90 ../progs.older/UTIL/anderson77.f90
echo 'UTIL/anderson.f90'
diff UTIL/anderson.f90 ../progs.older/UTIL/anderson.f90
echo 'UTIL/anderson_l95.f90'
diff UTIL/anderson_l95.f90 ../progs.older/UTIL/anderson_l95.f90
echo 'UTIL/fixfrags.f90'
diff UTIL/fixfrags.f90 ../progs.older/UTIL/fixfrags.f90
echo 'UTIL/fixfrags2.f90'
diff UTIL/fixfrags2.f90 ../progs.older/UTIL/fixfrags2.f90
echo 'UTIL/hampiece.f90'
diff UTIL/hampiece.f90 ../progs.older/UTIL/hampiece.f90
echo 'UTIL/hamtrans.f90'
diff UTIL/hamtrans.f90 ../progs.older/UTIL/hamtrans.f90
echo 'UTIL/mixer.f90'
diff UTIL/mixer.f90 ../progs.older/UTIL/mixer.f90
echo 'UTIL/push_atoms.f90'
diff UTIL/push_atoms.f90 ../progs.older/UTIL/push_atoms.f90
echo 'UTIL/writeout_ac.f90'
diff UTIL/writeout_ac.f90 ../progs.older/UTIL/writeout_ac.f90
echo 'UTIL/writeout_cd.f90'
diff UTIL/writeout_cd.f90 ../progs.older/UTIL/writeout_cd.f90
echo 'UTIL/writeout_charges.f90'
diff UTIL/writeout_charges.f90 ../progs.older/UTIL/writeout_charges.f90
echo 'UTIL/writeout_comph.f90'
diff UTIL/writeout_comph.f90 ../progs.older/UTIL/writeout_comph.f90
echo 'UTIL/writeout_neighbors.f90'
diff UTIL/writeout_neighbors.f90 ../progs.older/UTIL/writeout_neighbors.f90
echo 'UTIL/writeout_neighborsPP.f90'
diff UTIL/writeout_neighborsPP.f90 ../progs.older/UTIL/writeout_neighborsPP.f90
echo 'UTIL/writeout_xv.f90'
diff UTIL/writeout_xv.f90 ../progs.older/UTIL/writeout_xv.f90
echo 'UTIL_SPARSE/build_transpose.f90'
diff UTIL_SPARSE/build_transpose.f90 ../progs.older/UTIL_SPARSE/build_transpose.f90
echo 'UTIL_SPARSE/isnan.f90'
diff UTIL_SPARSE/isnan.f90 ../progs.older/UTIL_SPARSE/isnan.f90
echo 'UTIL_SPARSE/sparse_add.f90'
diff UTIL_SPARSE/sparse_add.f90 ../progs.older/UTIL_SPARSE/sparse_add.f90
echo 'UTIL_SPARSE/sparse_copy.f90'
diff UTIL_SPARSE/sparse_copy.f90 ../progs.older/UTIL_SPARSE/sparse_copy.f90
echo 'UTIL_SPARSE/sparse_getdimension.f90'
diff UTIL_SPARSE/sparse_getdimension.f90 ../progs.older/UTIL_SPARSE/sparse_getdimension.f90
echo 'UTIL_SPARSE/sparse_getpacksize.f90'
diff UTIL_SPARSE/sparse_getpacksize.f90 ../progs.older/UTIL_SPARSE/sparse_getpacksize.f90
echo 'UTIL_SPARSE/sparse_mask.f90'
diff UTIL_SPARSE/sparse_mask.f90 ../progs.older/UTIL_SPARSE/sparse_mask.f90
echo 'UTIL_SPARSE/sparse_mult.f90'
diff UTIL_SPARSE/sparse_mult.f90 ../progs.older/UTIL_SPARSE/sparse_mult.f90
echo 'UTIL_SPARSE/sparse_norm2.f90'
diff UTIL_SPARSE/sparse_norm2.f90 ../progs.older/UTIL_SPARSE/sparse_norm2.f90
echo 'UTIL_SPARSE/sparse_pack.f90'
diff UTIL_SPARSE/sparse_pack.f90 ../progs.older/UTIL_SPARSE/sparse_pack.f90
echo 'UTIL_SPARSE/sparse_pack_elements.f90'
diff UTIL_SPARSE/sparse_pack_elements.f90 ../progs.older/UTIL_SPARSE/sparse_pack_elements.f90
echo 'UTIL_SPARSE/sparse_pack_indices.f90'
diff UTIL_SPARSE/sparse_pack_indices.f90 ../progs.older/UTIL_SPARSE/sparse_pack_indices.f90
echo 'UTIL_SPARSE/sparse_unpack.f90'
diff UTIL_SPARSE/sparse_unpack.f90 ../progs.older/UTIL_SPARSE/sparse_unpack.f90
echo 'UTIL_SPARSE/sparse_unpack_elements.f90'
diff UTIL_SPARSE/sparse_unpack_elements.f90 ../progs.older/UTIL_SPARSE/sparse_unpack_elements.f90
echo 'UTIL_SPARSE/sparse_unpack_indices.f90'
diff UTIL_SPARSE/sparse_unpack_indices.f90 ../progs.older/UTIL_SPARSE/sparse_unpack_indices.f90
echo 'UTIL_SPARSE/sparse_vecprod.f90'
diff UTIL_SPARSE/sparse_vecprod.f90 ../progs.older/UTIL_SPARSE/sparse_vecprod.f90
diff VISUALIZATION/graceinit.f90 ../progs.older/VISUALIZATION/graceinit.f90
diff VISUALIZATION/graceupdate.f90 ../progs.older/VISUALIZATION/graceupdate.f90
diff VISUALIZATION/noTclMD.f90 ../progs.older/VISUALIZATION/noTclMD.f90
diff VISUALIZATION/nograce.f90 ../progs.older/VISUALIZATION/nograce.f90
diff VISUALIZATION/noxmgr.f90 ../progs.older/VISUALIZATION/noxmgr.f90
diff VISUALIZATION/tclmdtransfer.f90 ../progs.older/VISUALIZATION/tclmdtransfer.f90
diff VISUALIZATION/xmgrinit.f90 ../progs.older/VISUALIZATION/xmgrinit.f90
diff VISUALIZATION/xmgrupdate.f90 ../progs.older/VISUALIZATION/xmgrupdate.f90
diff XC/ceperley_alder.f90 ../progs.older/XC/ceperley_alder.f90

