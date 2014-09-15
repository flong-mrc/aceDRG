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
UNL	UNL	'.		'	non-polymer	24	24	.
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
UNL          O1     O     O       0       4.321       1.591       0.115
UNL          O2     O     O   0.000       2.016       5.476       0.217
UNL          N1     N     N   0.000       0.744       3.620       0.253
UNL          N2     N     N   0.000       1.475       0.196       0.227
UNL          N3     N     N   0.000       3.080       3.525       0.182
UNL          N4     N     N   0.000      -0.333       1.514       0.274
UNL          C1     C  CR56   0.000       1.919       1.507       0.215
UNL          C2     C  CR56   0.000       0.750       2.266       0.245
UNL          C3     C   CR6   0.000       3.188       2.138       0.172
UNL          C4     C   CR6   0.000       1.952       4.280       0.230
UNL          C5     C  CR15   0.000       0.117       0.279       0.259
UNL          C6     C   CH3   0.000      -0.536       4.359       0.275
UNL          C7     C   CH3   0.000       2.309      -1.046       0.175
UNL          C8     C   CH3   0.000       4.339       4.223       0.193
UNL          H1     H     H   0.000      -0.475      -0.459       0.272
UNL          H2     H     H   0.000      -0.436       5.202       0.739
UNL          H3     H     H   0.000      -0.788       4.608      -0.688
UNL          H4     H     H   0.000      -1.255       3.798       0.694
UNL          H5     H     H   0.000       3.146      -0.852       0.666
UNL          H6     H     H   0.000       2.528      -1.320      -0.837
UNL          H7     H     H   0.000       1.864      -1.798       0.593
UNL          H8     H     H   0.000       4.208       5.042      -0.102
UNL          H9     H     H   0.000       4.650       4.318       1.167
UNL         H10     H     H   0.000       4.922       3.795      -0.375
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
UNL          O1          C3      double     1.224   0.019
UNL          O2          C4      double     1.218   0.012
UNL          N1          C2      single     1.376   0.011
UNL          N1          C4      single     1.376   0.011
UNL          N1          C6      single     1.465   0.010
UNL          N2          C1      single     1.384   0.010
UNL          N2          C5      single     1.343   0.013
UNL          N2          C7      single     1.460   0.014
UNL          N3          C3      single     1.405   0.012
UNL          N3          C4      single     1.390   0.013
UNL          N3          C8      single     1.469   0.010
UNL          N4          C2      single     1.354   0.011
UNL          N4          C5      double     1.338   0.020
UNL          C1          C2      double     1.369   0.016
UNL          C1          C3      single     1.422   0.010
UNL          C5          H1      single     0.943   0.020
UNL          C6          H2      single     0.969   0.020
UNL          C6          H3      single     0.969   0.020
UNL          C6          H4      single     0.969   0.020
UNL          C7          H5      single     0.969   0.020
UNL          C7          H6      single     0.969   0.020
UNL          C7          H7      single     0.969   0.020
UNL          C8          H8      single     0.969   0.020
UNL          C8          H9      single     0.969   0.020
UNL          C8         H10      single     0.969   0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
UNL          C2          N1          C4     119.603    0.59
UNL          C2          N1          C6     121.111    0.73
UNL          C4          N1          C6     119.286    1.20
UNL          C1          N2          C5     105.782    0.93
UNL          C1          N2          C7     127.300    1.10
UNL          C5          N2          C7     126.918    0.90
UNL          C3          N3          C4     126.504    0.90
UNL          C3          N3          C8     117.108    1.11
UNL          C4          N3          C8     116.389    0.81
UNL          C2          N4          C5     103.393    0.94
UNL          N2          C1          C2     105.605    0.90
UNL          N2          C1          C3     131.594    0.88
UNL          C2          C1          C3     122.801    1.13
UNL          N1          C2          N4     126.343    0.73
UNL          N1          C2          C1     121.841    0.76
UNL          N4          C2          C1     111.816    0.88
UNL          O1          C3          N3     121.297    1.02
UNL          O1          C3          C1     126.725    1.16
UNL          N3          C3          C1     111.977    1.50
UNL          O2          C4          N1     121.462    0.79
UNL          O2          C4          N3     121.264    0.81
UNL          N1          C4          N3     117.274    0.65
UNL          N2          C5          N4     113.403    0.89
UNL          N2          C5          H1     123.065    1.13
UNL          N4          C5          H1     123.531    1.23
UNL          N1          C6          H2     109.463    1.35
UNL          N1          C6          H3     109.463    1.35
UNL          N1          C6          H4     109.463    1.35
UNL          H2          C6          H3     109.417    2.07
UNL          H2          C6          H4     109.417    2.07
UNL          H3          C6          H4     109.417    2.07
UNL          N2          C7          H5     109.514    1.29
UNL          N2          C7          H6     109.514    1.29
UNL          N2          C7          H7     109.514    1.29
UNL          H5          C7          H6     109.415    1.67
UNL          H5          C7          H7     109.415    1.67
UNL          H6          C7          H7     109.415    1.67
UNL          N3          C8          H8     109.479    1.25
UNL          N3          C8          H9     109.479    1.25
UNL          N3          C8         H10     109.479    1.25
UNL          H8          C8          H9     109.417    2.07
UNL          H8          C8         H10     109.417    2.07
UNL          H9          C8         H10     109.417    2.07
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
UNL     P_sp2_sp2_1          C1          C2          N1          C4       0.000   10.00     2
UNL     P_sp2_sp2_2          C1          C2          N1          C6     180.000   10.00     2
UNL     P_sp2_sp2_3          N4          C2          N1          C4     180.000   10.00     2
UNL     P_sp2_sp2_4          N4          C2          N1          C6       0.000   10.00     2
UNL     P_sp2_sp2_5          N2          C1          C2          N4       0.000   10.00     2
UNL     P_sp2_sp2_6          N2          C1          C2          N1     180.000   10.00     2
UNL     P_sp2_sp2_7          C3          C1          C2          N4     180.000   10.00     2
UNL     P_sp2_sp2_8          C3          C1          C2          N1       0.000   10.00     2
UNL     P_sp2_sp2_9          C2          C1          C3          N3       0.000   10.00     2
UNL    P_sp2_sp2_10          C2          C1          C3          O1     180.000   10.00     2
UNL    P_sp2_sp2_11          N2          C1          C3          N3     180.000   10.00     2
UNL    P_sp2_sp2_12          N2          C1          C3          O1       0.000   10.00     2
UNL    P_sp2_sp2_13          C1          C3          N3          C4       0.000   10.00     2
UNL    P_sp2_sp2_14          C1          C3          N3          C8     180.000   10.00     2
UNL    P_sp2_sp2_15          O1          C3          N3          C4     180.000   10.00     2
UNL    P_sp2_sp2_16          O1          C3          N3          C8       0.000   10.00     2
UNL    P_sp2_sp2_17          N1          C4          N3          C3       0.000   10.00     2
UNL    P_sp2_sp2_18          N1          C4          N3          C8     180.000   10.00     2
UNL    P_sp2_sp2_19          O2          C4          N3          C3     180.000   10.00     2
UNL    P_sp2_sp2_20          O2          C4          N3          C8       0.000   10.00     2
UNL    P_sp2_sp2_21          C2          C1          N2          C5       0.000   10.00     2
UNL    P_sp2_sp2_22          C2          C1          N2          C7     180.000   10.00     2
UNL    P_sp2_sp2_23          C3          C1          N2          C5     180.000   10.00     2
UNL    P_sp2_sp2_24          C3          C1          N2          C7       0.000   10.00     2
UNL    P_sp2_sp2_25          N2          C1          C2          N4       0.000   10.00     2
UNL    P_sp2_sp2_26          N2          C1          C2          N1     180.000   10.00     2
UNL    P_sp2_sp2_27          C3          C1          C2          N4     180.000   10.00     2
UNL    P_sp2_sp2_28          C3          C1          C2          N1       0.000   10.00     2
UNL    P_sp2_sp2_29          C1          C2          N4          C5       0.000   10.00     2
UNL    P_sp2_sp2_30          N1          C2          N4          C5     180.000   10.00     2
UNL    P_sp2_sp2_31          N2          C5          N4          C2       0.000   10.00     2
UNL    P_sp2_sp2_32          H1          C5          N4          C2     180.000   10.00     2
UNL    P_sp2_sp2_33          N3          C4          N1          C2       0.000   10.00     2
UNL    P_sp2_sp2_34          N3          C4          N1          C6     180.000   10.00     2
UNL    P_sp2_sp2_35          O2          C4          N1          C2     180.000   10.00     2
UNL    P_sp2_sp2_36          O2          C4          N1          C6       0.000   10.00     2
UNL       sp2_sp3_1          C2          N1          C6          H2     150.000   10.00     6
UNL       sp2_sp3_2          C2          N1          C6          H3     -90.000   10.00     6
UNL       sp2_sp3_3          C2          N1          C6          H4      30.000   10.00     6
UNL       sp2_sp3_4          C4          N1          C6          H2     -30.000   10.00     6
UNL       sp2_sp3_5          C4          N1          C6          H3      90.000   10.00     6
UNL       sp2_sp3_6          C4          N1          C6          H4    -150.000   10.00     6
UNL    P_sp2_sp2_37          N4          C5          N2          C1       0.000   10.00     2
UNL    P_sp2_sp2_38          N4          C5          N2          C7     180.000   10.00     2
UNL    P_sp2_sp2_39          H1          C5          N2          C1     180.000   10.00     2
UNL    P_sp2_sp2_40          H1          C5          N2          C7       0.000   10.00     2
UNL       sp2_sp3_7          C1          N2          C7          H5     150.000   10.00     6
UNL       sp2_sp3_8          C1          N2          C7          H6     -90.000   10.00     6
UNL       sp2_sp3_9          C1          N2          C7          H7      30.000   10.00     6
UNL      sp2_sp3_10          C5          N2          C7          H5     -30.000   10.00     6
UNL      sp2_sp3_11          C5          N2          C7          H6      90.000   10.00     6
UNL      sp2_sp3_12          C5          N2          C7          H7    -150.000   10.00     6
UNL      sp2_sp3_13          C3          N3          C8          H8     150.000   10.00     6
UNL      sp2_sp3_14          C3          N3          C8          H9     -90.000   10.00     6
UNL      sp2_sp3_15          C3          N3          C8         H10      30.000   10.00     6
UNL      sp2_sp3_16          C4          N3          C8          H8     -30.000   10.00     6
UNL      sp2_sp3_17          C4          N3          C8          H9      90.000   10.00     6
UNL      sp2_sp3_18          C4          N3          C8         H10    -150.000   10.00     6
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
UNL    plan-1          C1   0.020
UNL    plan-1          C2   0.020
UNL    plan-1          C3   0.020
UNL    plan-1          C4   0.020
UNL    plan-1          C5   0.020
UNL    plan-1          C6   0.020
UNL    plan-1          C7   0.020
UNL    plan-1          C8   0.020
UNL    plan-1          H1   0.020
UNL    plan-1          N1   0.020
UNL    plan-1          N2   0.020
UNL    plan-1          N3   0.020
UNL    plan-1          N4   0.020
UNL    plan-1          O1   0.020
UNL    plan-1          O2   0.020
