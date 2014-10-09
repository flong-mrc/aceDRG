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
UNL	UNL	'.		'	non-polymer	20	20	.
# ------------------------------------------------------
# ------------------------------------------------------
#
# --- DESCRIPTION OF MONOMERS ---
#
data_comp_UNL
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
UNL          O1     O   OH1       0       0.000       0.000       0.000
UNL         H10     H     H   0.000       0.685      -0.652       0.198
UNL          C1     C   CH1   0.000      -1.283      -0.626       0.058
UNL          H1     H     H   0.000      -1.330      -1.478      -0.635
UNL          C2     C   CH1   0.000      -2.404       0.398      -0.257
UNL          H2     H     H   0.000      -2.167       1.382       0.172
UNL          O2     O   OH1   0.000      -2.626       0.495      -1.665
UNL         HO2     H     H   0.000      -1.851       0.891      -2.086
UNL          C3     C   CH1   0.000      -3.627      -0.233       0.452
UNL          H3     H     H   0.000      -4.240      -0.796      -0.266
UNL          O3     O   OH1   0.000      -4.409       0.776       1.094
UNL         HO3     H     H   0.000      -4.788       1.362       0.425
UNL          O4     O    O2   0.000      -1.582      -1.050       1.409
UNL          C4     C   CH1   0.000      -3.010      -1.184       1.496
UNL          H4     H     H   0.000      -3.351      -0.905       2.503
UNL          C5     C   CH2   0.000      -3.417      -2.627       1.191
UNL          H5     H     H   0.000      -2.991      -2.930       0.232
UNL         H5A     H     H   0.000      -4.506      -2.695       1.141
UNL          O5     O   OH1   0.000      -2.932      -3.488       2.224
UNL         HO5     H     H   0.000      -3.189      -4.399       2.030
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
UNL          C5          O5      single     1.421   0.012
UNL          O5         HO5      single     0.850   0.020
UNL          C5          H5      single     0.979   0.020
UNL          C5         H5A      single     0.979   0.020
UNL          C4          C5      single     1.511   0.011
UNL          C4          H4      single     0.985   0.020
UNL          C3          C4      single     1.525   0.014
UNL          C3          O3      single     1.418   0.010
UNL          C3          H3      single     0.986   0.020
UNL          O3         HO3      single     0.847   0.020
UNL          C2          C3      single     1.528   0.010
UNL          C2          H2      single     1.020   0.020
UNL          C2          O2      single     1.415   0.010
UNL          O2         HO2      single     0.847   0.020
UNL          C1          C2      single     1.523   0.011
UNL          C1          O4      single     1.417   0.013
UNL          O1          C1      single     1.394   0.010
UNL          O4          C4      single     1.442   0.010
UNL          C1          H1      single     0.995   0.020
UNL          O1         H10      single     0.850   0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
UNL          C1          O1         H10     107.959    0.55
UNL          C2          C1          O4     104.548    1.49
UNL          C2          C1          O1     110.864    3.00
UNL          C2          C1          H1     112.228    2.27
UNL          O4          C1          O1     110.422    1.87
UNL          O4          C1          H1     109.421    2.78
UNL          O1          C1          H1     110.059    1.02
UNL          C3          C2          H2     109.137    1.14
UNL          C3          C2          O2     113.788    1.06
UNL          C3          C2          C1     101.308    0.88
UNL          H2          C2          O2     109.236    1.88
UNL          H2          C2          C1     110.546    3.00
UNL          O2          C2          C1     114.605    2.71
UNL          C2          O2         HO2     108.559    3.00
UNL          C4          C3          O3     110.899    2.25
UNL          C4          C3          H3     110.929    1.65
UNL          C4          C3          C2     103.251    0.72
UNL          O3          C3          H3     110.556    1.72
UNL          O3          C3          C2     114.175    0.79
UNL          H3          C3          C2     108.404    1.32
UNL          C3          O3         HO3     108.965    3.00
UNL          C1          O4          C4     109.248    1.57
UNL          C5          C4          H4     108.898    1.24
UNL          C5          C4          C3     114.835    1.46
UNL          C5          C4          O4     110.314    1.88
UNL          H4          C4          C3     109.305    1.32
UNL          H4          C4          O4     109.145    2.35
UNL          C3          C4          O4     106.153    0.78
UNL          O5          C5          H5     109.300    1.07
UNL          O5          C5         H5A     109.300    1.07
UNL          O5          C5          C4     111.386    2.14
UNL          H5          C5         H5A     108.237    1.64
UNL          H5          C5          C4     109.309    1.09
UNL         H5A          C5          C4     109.309    1.09
UNL          C5          O5         HO5     109.099    2.67
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
UNL       sp3_sp3_1          O4          C1          C2          C3      60.000   10.00     3
UNL       sp3_sp3_2          O4          C1          C2          H2     180.000   10.00     3
UNL       sp3_sp3_3          O4          C1          C2          O2     -60.000   10.00     3
UNL       sp3_sp3_4          H1          C1          C2          C3     -60.000   10.00     3
UNL       sp3_sp3_5          H1          C1          C2          H2      60.000   10.00     3
UNL       sp3_sp3_6          H1          C1          C2          O2     180.000   10.00     3
UNL       sp3_sp3_7          O1          C1          C2          C3     180.000   10.00     3
UNL       sp3_sp3_8          O1          C1          C2          H2     -60.000   10.00     3
UNL       sp3_sp3_9          O1          C1          C2          O2      60.000   10.00     3
UNL      sp3_sp3_10          C1          C2          C3          C4     -60.000   10.00     3
UNL      sp3_sp3_11          C1          C2          C3          O3      60.000   10.00     3
UNL      sp3_sp3_12          C1          C2          C3          H3     180.000   10.00     3
UNL      sp3_sp3_13          O2          C2          C3          C4     180.000   10.00     3
UNL      sp3_sp3_14          O2          C2          C3          O3     -60.000   10.00     3
UNL      sp3_sp3_15          O2          C2          C3          H3      60.000   10.00     3
UNL      sp3_sp3_16          H2          C2          C3          C4      60.000   10.00     3
UNL      sp3_sp3_17          H2          C2          C3          O3     180.000   10.00     3
UNL      sp3_sp3_18          H2          C2          C3          H3     -60.000   10.00     3
UNL      sp3_sp3_19          C2          C3          C4          O4      60.000   10.00     3
UNL      sp3_sp3_20          C2          C3          C4          H4     180.000   10.00     3
UNL      sp3_sp3_21          C2          C3          C4          C5     -60.000   10.00     3
UNL      sp3_sp3_22          H3          C3          C4          O4     -60.000   10.00     3
UNL      sp3_sp3_23          H3          C3          C4          H4      60.000   10.00     3
UNL      sp3_sp3_24          H3          C3          C4          C5     180.000   10.00     3
UNL      sp3_sp3_25          O3          C3          C4          O4     180.000   10.00     3
UNL      sp3_sp3_26          O3          C3          C4          H4     -60.000   10.00     3
UNL      sp3_sp3_27          O3          C3          C4          C5      60.000   10.00     3
UNL      sp3_sp3_28          C3          C4          O4          C1     -60.000   10.00     3
UNL      sp3_sp3_29          C5          C4          O4          C1      60.000   10.00     3
UNL      sp3_sp3_30          H4          C4          O4          C1     180.000   10.00     3
UNL      sp3_sp3_31          C4          C5          O5         HO5     180.000   10.00     3
UNL      sp3_sp3_32          H5          C5          O5         HO5     -60.000   10.00     3
UNL      sp3_sp3_33         H5A          C5          O5         HO5      60.000   10.00     3
UNL      sp3_sp3_34          O4          C4          C5          O5     180.000   10.00     3
UNL      sp3_sp3_35          O4          C4          C5          H5     -60.000   10.00     3
UNL      sp3_sp3_36          O4          C4          C5         H5A      60.000   10.00     3
UNL      sp3_sp3_37          C3          C4          C5          O5      60.000   10.00     3
UNL      sp3_sp3_38          C3          C4          C5          H5     180.000   10.00     3
UNL      sp3_sp3_39          C3          C4          C5         H5A     -60.000   10.00     3
UNL      sp3_sp3_40          H4          C4          C5          O5     -60.000   10.00     3
UNL      sp3_sp3_41          H4          C4          C5          H5      60.000   10.00     3
UNL      sp3_sp3_42          H4          C4          C5         H5A     180.000   10.00     3
UNL      sp3_sp3_43          C4          C3          O3         HO3     180.000   10.00     3
UNL      sp3_sp3_44          H3          C3          O3         HO3     -60.000   10.00     3
UNL      sp3_sp3_45          C2          C3          O3         HO3      60.000   10.00     3
UNL      sp3_sp3_46          C3          C2          O2         HO2     180.000   10.00     3
UNL      sp3_sp3_47          C1          C2          O2         HO2     -60.000   10.00     3
UNL      sp3_sp3_48          H2          C2          O2         HO2      60.000   10.00     3
UNL      sp3_sp3_49          C2          C1          O4          C4      60.000   10.00     3
UNL      sp3_sp3_50          O1          C1          O4          C4     180.000   10.00     3
UNL      sp3_sp3_51          H1          C1          O4          C4     -60.000   10.00     3
UNL      sp3_sp3_52          C2          C1          O1         H10     180.000   10.00     3
UNL      sp3_sp3_53          H1          C1          O1         H10     -60.000   10.00     3
UNL      sp3_sp3_54          O4          C1          O1         H10      60.000   10.00     3
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
UNL    chir_01    C4    C5    C3    O4    positiv
UNL    chir_02    C3    C4    O3    C2    negativ
UNL    chir_03    C2    C3    O2    C1    positiv
UNL    chir_04    C1    C2    O4    O1    positiv
