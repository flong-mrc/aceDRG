global_
_lib_name         ?
_lib_version      ?
_lib_update       ?
# ------------------------------------------------
#
# ---   LIST OF MONOMERS ---
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
NH4	NH4	'.		'	non-polymer	5	1	.
# ------------------------------------------------------
# ------------------------------------------------------
#
# --- DESCRIPTION OF MONOMERS ---
#
data_comp_NH4
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
NH4           N     N    NT       1      -0.581       0.302      -0.782
NH4         HN1     H     H       0      -0.581       0.302       0.262
NH4         HN2     H     H       0      -0.581       1.287      -1.130
NH4         HN3     H     H       0       0.272      -0.190      -1.130
NH4         HN4     H     H       0      -1.434      -0.190      -1.130
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
NH4           N         HN1      SINGLE       n     0.912  0.0200
NH4           N         HN2      SINGLE       n     0.912  0.0200
NH4           N         HN3      SINGLE       n     0.912  0.0200
NH4           N         HN4      SINGLE       n     0.912  0.0200
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
NH4         HN1           N         HN2     109.380    3.00
NH4         HN1           N         HN3     109.380    3.00
NH4         HN1           N         HN4     109.380    3.00
NH4         HN2           N         HN3     109.380    3.00
NH4         HN2           N         HN4     109.380    3.00
NH4         HN3           N         HN4     109.380    3.00
