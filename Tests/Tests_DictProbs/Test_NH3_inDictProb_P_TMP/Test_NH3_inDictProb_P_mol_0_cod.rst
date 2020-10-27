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
NH3	NH3	'.		'	non-polymer	4	1	.
# ------------------------------------------------------
# ------------------------------------------------------
#
# --- DESCRIPTION OF MONOMERS ---
#
data_comp_NH3
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
NH3           N     N   NT3       0      -0.211       0.868       0.434
NH3         HN1     H     H       0      -0.169       0.892       1.478
NH3         HN2     H     H       0      -0.169       1.860       0.110
NH3         HN3     H     H       0       0.669       0.409       0.110
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
NH3           N         HN1      SINGLE       n     0.928  0.0200
NH3           N         HN2      SINGLE       n     0.928  0.0200
NH3           N         HN3      SINGLE       n     0.928  0.0200
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
NH3         HN1           N         HN2     106.383    3.00
NH3         HN1           N         HN3     106.383    3.00
NH3         HN2           N         HN3     106.383    3.00
