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
UNL	UNL	'.		'	non-polymer	32	32	.
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
UNL          O1     O     O       0       0.430      -3.214      -2.210
UNL          O2     O     O   0.000       1.264      -2.839       2.032
UNL          O3     O     O   0.000      -1.500      -6.196       0.706
UNL          N1     N  NR16   0.000      -0.977      -4.460      -0.827
UNL          N2     N  NR16   0.000      -0.001      -4.589       1.375
UNL          C1     C    CT   0.000      -0.111      -2.661       0.168
UNL          C2     C   CH2   0.000      -0.601      -1.208       0.505
UNL          C3     C   CH1   0.000       0.142       0.096       0.095
UNL          C4     C   CH2   0.000      -1.495      -3.127       0.376
UNL          C5     C   CR6   0.000      -0.169      -3.374      -1.141
UNL          C6     C   CR6   0.000       0.544      -3.350       1.231
UNL          C7     C   CH3   0.000      -0.597       1.397       0.540
UNL          C8     C   CH3   0.000       1.546       0.111       0.599
UNL          C9     C    C1   0.000      -1.560      -4.669       0.131
UNL         C10     C   CR6   0.000      -0.844      -5.179       0.440
UNL         C11     C    C2   0.000      -1.865      -5.519       1.006
UNL          H1     H     H   0.000      -1.598      -4.703      -1.334
UNL          H2     H     H   0.000       0.197      -5.025       2.072
UNL          H3     H     H   0.000      -0.639      -1.113       1.548
UNL          H4     H     H   0.000      -1.546      -1.132       0.247
UNL          H5     H     H   0.000       0.094       0.138      -0.876
UNL          H6     H     H   0.000      -1.771      -2.967       1.266
UNL          H7     H     H   0.000      -2.100      -2.695      -0.238
UNL          H8     H     H   0.000      -0.626       1.438       1.572
UNL          H9     H     H   0.000      -1.564       1.350       0.226
UNL         H10     H     H   0.000      -0.228       2.174       0.282
UNL         H11     H     H   0.000       1.593       0.127       1.521
UNL         H12     H     H   0.000       2.004      -0.709       0.224
UNL         H13     H     H   0.000       2.078       0.829       0.233
UNL         H14     H     H   0.000      -1.333      -4.963      -0.712
UNL         H15     H     H   0.000      -1.925      -6.461       0.802
UNL         H16     H     H   0.000      -2.130      -5.184       1.960
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
UNL          O1          C5      double     1.218   0.011
UNL          O2          C6      double     1.218   0.011
UNL          O3         C10      double     1.225   0.013
UNL          N1          C5      single     1.370   0.010
UNL          N1         C10      single     1.372   0.011
UNL          N2          C6      single     1.370   0.010
UNL          N2         C10      single     1.372   0.011
UNL          C1          C2      single     1.548   0.010
UNL          C1          C4      single     1.551   0.010
UNL          C1          C5      single     1.520   0.010
UNL          C1          C6      single     1.520   0.010
UNL          C2          C3      single     1.517   0.018
UNL          C3          C7      single     1.516   0.020
UNL          C3          C8      single     1.516   0.020
UNL          C4          C9      single     1.499   0.011
UNL          C9         C11      double     1.291   0.020
UNL          N1          H1      single     0.877   0.020
UNL          N2          H2      single     0.877   0.020
UNL          C2          H3      single     0.980   0.019
UNL          C2          H4      single     0.980   0.019
UNL          C3          H5      single     0.987   0.017
UNL          C4          H6      single     0.982   0.014
UNL          C4          H7      single     0.982   0.014
UNL          C7          H8      single     0.974   0.020
UNL          C7          H9      single     0.974   0.020
UNL          C7         H10      single     0.974   0.020
UNL          C8         H11      single     0.974   0.020
UNL          C8         H12      single     0.974   0.020
UNL          C8         H13      single     0.974   0.020
UNL          C9         H14      single     0.948   0.020
UNL         C11         H15      single     0.949   0.020
UNL         C11         H16      single     0.949   0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
UNL          C5          N1         C10     125.880    1.92
UNL          C5          N1          H1     117.283    2.19
UNL         C10          N1          H1     116.838    2.16
UNL          C6          N2         C10     125.880    1.92
UNL          C6          N2          H2     117.283    2.19
UNL         C10          N2          H2     116.838    2.16
UNL          C2          C1          C4     109.313    1.94
UNL          C2          C1          C5     108.254    1.00
UNL          C2          C1          C6     108.254    1.00
UNL          C4          C1          C5     108.346    0.82
UNL          C4          C1          C6     108.346    0.82
UNL          C5          C1          C6     114.092    0.38
UNL          C1          C2          C3     112.921    2.66
UNL          C1          C2          H3     108.589    0.39
UNL          C1          C2          H4     108.589    0.39
UNL          C3          C2          H3     107.989    0.99
UNL          C3          C2          H4     107.989    0.99
UNL          H3          C2          H4     107.577    1.40
UNL          C2          C3          C7     110.888    1.67
UNL          C2          C3          C8     110.888    1.67
UNL          C2          C3          H5     107.533    0.87
UNL          C7          C3          C8     110.499    2.07
UNL          C7          C3          H5     107.911    2.01
UNL          C8          C3          H5     107.911    2.01
UNL          C1          C4          C9     113.426    1.02
UNL          C1          C4          H6     108.862    0.22
UNL          C1          C4          H7     108.862    0.22
UNL          C9          C4          H6     108.956    0.79
UNL          C9          C4          H7     108.956    0.79
UNL          H6          C4          H7     107.646    0.56
UNL          O1          C5          N1     120.198    1.22
UNL          O1          C5          C1     121.682    0.95
UNL          N1          C5          C1     118.120    0.73
UNL          O2          C6          N2     120.198    1.22
UNL          O2          C6          C1     121.682    0.95
UNL          N2          C6          C1     118.120    0.73
UNL          C3          C7          H8     109.485    0.79
UNL          C3          C7          H9     109.485    0.79
UNL          C3          C7         H10     109.485    0.79
UNL          H8          C7          H9     109.300    1.49
UNL          H8          C7         H10     109.300    1.49
UNL          H9          C7         H10     109.300    1.49
UNL          C3          C8         H11     109.485    0.79
UNL          C3          C8         H12     109.485    0.79
UNL          C3          C8         H13     109.485    0.79
UNL         H11          C8         H12     109.300    1.49
UNL         H11          C8         H13     109.300    1.49
UNL         H12          C8         H13     109.300    1.49
UNL          C4          C9         C11     125.278    2.02
UNL          C4          C9         H14     117.905    2.92
UNL         C11          C9         H14     116.817    2.73
UNL          O3         C10          N1     121.882    1.20
UNL          O3         C10          N2     121.882    1.20
UNL          N1         C10          N2     116.236    1.09
UNL          C9         C11         H15     120.031    2.51
UNL          C9         C11         H16     120.031    2.51
UNL         H15         C11         H16     119.938    2.37
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
UNL       sp2_sp2_1          C1          C5          N1         C10       0.000   10.00     2
UNL       sp2_sp2_2          C1          C5          N1          H1     180.000   10.00     2
UNL       sp2_sp2_3          O1          C5          N1         C10     180.000   10.00     2
UNL       sp2_sp2_4          O1          C5          N1          H1       0.000   10.00     2
UNL       sp2_sp3_1          N1          C5          C1          C6       0.000   10.00     6
UNL       sp2_sp3_2          N1          C5          C1          C2     120.000   10.00     6
UNL       sp2_sp3_3          N1          C5          C1          C4    -120.000   10.00     6
UNL       sp2_sp3_4          O1          C5          C1          C6     180.000   10.00     6
UNL       sp2_sp3_5          O1          C5          C1          C2     -60.000   10.00     6
UNL       sp2_sp3_6          O1          C5          C1          C4      60.000   10.00     6
UNL       sp2_sp3_7          N2          C6          C1          C5       0.000   10.00     6
UNL       sp2_sp3_8          N2          C6          C1          C2     120.000   10.00     6
UNL       sp2_sp3_9          N2          C6          C1          C4    -120.000   10.00     6
UNL      sp2_sp3_10          O2          C6          C1          C5     180.000   10.00     6
UNL      sp2_sp3_11          O2          C6          C1          C2     -60.000   10.00     6
UNL      sp2_sp3_12          O2          C6          C1          C4      60.000   10.00     6
UNL       sp2_sp2_5          C1          C6          N2         C10       0.000   10.00     2
UNL       sp2_sp2_6          C1          C6          N2          H2     180.000   10.00     2
UNL       sp2_sp2_7          O2          C6          N2         C10     180.000   10.00     2
UNL       sp2_sp2_8          O2          C6          N2          H2       0.000   10.00     2
UNL       sp2_sp2_9          N1         C10          N2          C6       0.000   10.00     2
UNL      sp2_sp2_10          N1         C10          N2          H2     180.000   10.00     2
UNL      sp2_sp2_11          O3         C10          N2          C6     180.000   10.00     2
UNL      sp2_sp2_12          O3         C10          N2          H2       0.000   10.00     2
UNL      sp2_sp2_13          N2         C10          N1          C5       0.000   10.00     2
UNL      sp2_sp2_14          N2         C10          N1          H1     180.000   10.00     2
UNL      sp2_sp2_15          O3         C10          N1          C5     180.000   10.00     2
UNL      sp2_sp2_16          O3         C10          N1          H1       0.000   10.00     2
UNL       sp3_sp3_1          C4          C1          C2          C3     180.000   10.00     3
UNL       sp3_sp3_2          C4          C1          C2          H3     -60.000   10.00     3
UNL       sp3_sp3_3          C4          C1          C2          H4      60.000   10.00     3
UNL       sp3_sp3_4          C5          C1          C2          C3      60.000   10.00     3
UNL       sp3_sp3_5          C5          C1          C2          H3     180.000   10.00     3
UNL       sp3_sp3_6          C5          C1          C2          H4     -60.000   10.00     3
UNL       sp3_sp3_7          C6          C1          C2          C3     -60.000   10.00     3
UNL       sp3_sp3_8          C6          C1          C2          H3      60.000   10.00     3
UNL       sp3_sp3_9          C6          C1          C2          H4     180.000   10.00     3
UNL      sp3_sp3_10          C2          C1          C4          C9     180.000   10.00     3
UNL      sp3_sp3_11          C2          C1          C4          H6     -60.000   10.00     3
UNL      sp3_sp3_12          C2          C1          C4          H7      60.000   10.00     3
UNL      sp3_sp3_13          C5          C1          C4          C9      60.000   10.00     3
UNL      sp3_sp3_14          C5          C1          C4          H6     180.000   10.00     3
UNL      sp3_sp3_15          C5          C1          C4          H7     -60.000   10.00     3
UNL      sp3_sp3_16          C6          C1          C4          C9     -60.000   10.00     3
UNL      sp3_sp3_17          C6          C1          C4          H6      60.000   10.00     3
UNL      sp3_sp3_18          C6          C1          C4          H7     180.000   10.00     3
UNL      sp3_sp3_19          C1          C2          C3          C7     180.000   10.00     3
UNL      sp3_sp3_20          C1          C2          C3          C8     -60.000   10.00     3
UNL      sp3_sp3_21          C1          C2          C3          H5      60.000   10.00     3
UNL      sp3_sp3_22          H3          C2          C3          C7      60.000   10.00     3
UNL      sp3_sp3_23          H3          C2          C3          C8     180.000   10.00     3
UNL      sp3_sp3_24          H3          C2          C3          H5     -60.000   10.00     3
UNL      sp3_sp3_25          H4          C2          C3          C7     -60.000   10.00     3
UNL      sp3_sp3_26          H4          C2          C3          C8      60.000   10.00     3
UNL      sp3_sp3_27          H4          C2          C3          H5     180.000   10.00     3
UNL      sp3_sp3_28          C2          C3          C7          H8     180.000   10.00     3
UNL      sp3_sp3_29          C2          C3          C7          H9     -60.000   10.00     3
UNL      sp3_sp3_30          C2          C3          C7         H10      60.000   10.00     3
UNL      sp3_sp3_31          C8          C3          C7          H8      60.000   10.00     3
UNL      sp3_sp3_32          C8          C3          C7          H9     180.000   10.00     3
UNL      sp3_sp3_33          C8          C3          C7         H10     -60.000   10.00     3
UNL      sp3_sp3_34          H5          C3          C7          H8     -60.000   10.00     3
UNL      sp3_sp3_35          H5          C3          C7          H9      60.000   10.00     3
UNL      sp3_sp3_36          H5          C3          C7         H10     180.000   10.00     3
UNL      sp3_sp3_37          C2          C3          C8         H11     180.000   10.00     3
UNL      sp3_sp3_38          C2          C3          C8         H12     -60.000   10.00     3
UNL      sp3_sp3_39          C2          C3          C8         H13      60.000   10.00     3
UNL      sp3_sp3_40          C7          C3          C8         H11      60.000   10.00     3
UNL      sp3_sp3_41          C7          C3          C8         H12     180.000   10.00     3
UNL      sp3_sp3_42          C7          C3          C8         H13     -60.000   10.00     3
UNL      sp3_sp3_43          H5          C3          C8         H11     -60.000   10.00     3
UNL      sp3_sp3_44          H5          C3          C8         H12      60.000   10.00     3
UNL      sp3_sp3_45          H5          C3          C8         H13     180.000   10.00     3
UNL      sp2_sp3_13         C11          C9          C4          H6       0.000   10.00     6
UNL      sp2_sp3_14         C11          C9          C4          C1     120.000   10.00     6
UNL      sp2_sp3_15         C11          C9          C4          H7    -120.000   10.00     6
UNL      sp2_sp3_16         H14          C9          C4          H6     180.000   10.00     6
UNL      sp2_sp3_17         H14          C9          C4          C1     -60.000   10.00     6
UNL      sp2_sp3_18         H14          C9          C4          H7      60.000   10.00     6
UNL      sp2_sp2_17         H15         C11          C9          C4     180.000   10.00     2
UNL      sp2_sp2_18         H15         C11          C9         H14       0.000   10.00     2
UNL      sp2_sp2_19         H16         C11          C9          C4       0.000   10.00     2
UNL      sp2_sp2_20         H16         C11          C9         H14     180.000   10.00     2
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
UNL    chi1    C1    C2    C4    C5    both
UNL    chi2    C3    C2    C7    C8    both
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
UNL    plan-1         C10   0.020
UNL    plan-1          C5   0.020
UNL    plan-1          H1   0.020
UNL    plan-1          N1   0.020
UNL    plan-2         C10   0.020
UNL    plan-2          C6   0.020
UNL    plan-2          H2   0.020
UNL    plan-2          N2   0.020
UNL    plan-3          C1   0.020
UNL    plan-3          C5   0.020
UNL    plan-3          N1   0.020
UNL    plan-3          O1   0.020
UNL    plan-4          C1   0.020
UNL    plan-4          C6   0.020
UNL    plan-4          N2   0.020
UNL    plan-4          O2   0.020
UNL    plan-5         C11   0.020
UNL    plan-5          C4   0.020
UNL    plan-5          C9   0.020
UNL    plan-5         H14   0.020
UNL    plan-6         C10   0.020
UNL    plan-6          N1   0.020
UNL    plan-6          N2   0.020
UNL    plan-6          O3   0.020
UNL    plan-7         C11   0.020
UNL    plan-7          C9   0.020
UNL    plan-7         H15   0.020
UNL    plan-7         H16   0.020
