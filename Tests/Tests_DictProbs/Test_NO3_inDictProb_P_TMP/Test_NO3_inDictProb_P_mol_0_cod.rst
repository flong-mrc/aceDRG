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
NO3	NO3	'.		'	non-polymer	4	4	.
# ------------------------------------------------------
# ------------------------------------------------------
#
# --- DESCRIPTION OF MONOMERS ---
#
data_comp_NO3
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
NO3          O3     O    OC      -1       0.112      -0.088      -0.003
NO3           N     N     N       1      -1.111       0.537      -0.000
NO3          O1     O     O       0      -2.158      -0.141       0.000
NO3          O2     O    OC      -1      -1.182       1.910       0.002
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.aromatic
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
NO3           N          O1      DOUBLE       n     1.226  0.0162
NO3           N          O2      SINGLE       n     1.226  0.0162
NO3          O3           N      SINGLE       n     1.226  0.0162
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
NO3          O1           N          O2     120.000    1.59
NO3          O1           N          O3     119.996    1.59
NO3          O2           N          O3     119.996    1.59
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
NO3    plan-1           N   0.020
NO3    plan-1          O1   0.020
NO3    plan-1          O2   0.020
NO3    plan-1          O3   0.020
