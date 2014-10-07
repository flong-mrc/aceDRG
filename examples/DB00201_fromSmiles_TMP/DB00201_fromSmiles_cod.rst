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
UNL          O1     O     O       0       0.000       0.000       0.000
UNL          C6     C   CR6   0.000       0.795       0.453       1.037
UNL          N3     N   NR6   0.000       0.412      -0.409       2.119
UNL          C7     C   CH3   0.000       1.660      -0.942       2.808
UNL          H9     H     H   0.000       2.254      -1.492       2.118
UNL          H8     H     H   0.000       2.239      -0.137       3.192
UNL          H7     H     H   0.000       1.389      -1.583       3.613
UNL          C2     C  CR56   0.000      -0.961      -0.412       1.834
UNL          N1     N  NRD5   0.000      -1.596      -1.469       1.427
UNL          C1     C  CR15   0.000      -2.327      -1.396       0.283
UNL          H3     H     H   0.000      -2.832      -2.235      -0.178
UNL          N2     N   NR6   0.000       0.445       1.890       1.011
UNL          C5     C   CH3   0.000       1.784       2.508       1.064
UNL          H6     H     H   0.000       2.278       2.210       1.953
UNL          H5     H     H   0.000       2.353       2.197       0.226
UNL          H4     H     H   0.000       1.689       3.564       1.052
UNL          C4     C   CR6   0.000      -0.915       2.045       0.493
UNL           O     O     O   0.000      -0.973       2.811      -0.484
UNL          C3     C  CR56   0.000      -1.437       0.498       0.670
UNL           N     N   NR5   0.000      -2.338      -0.217      -0.168
UNL           C     C   CH3   0.000      -3.461       0.450      -0.832
UNL          H2     H     H   0.000      -4.095       0.875      -0.101
UNL          H1     H     H   0.000      -3.089       1.209      -1.467
UNL           H     H     H   0.000      -4.000      -0.259      -1.402
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
UNL           N           C      single     1.460   0.014
UNL          C1           N      single     1.343   0.013
UNL          N1          C1      double     1.338   0.020
UNL          C2          N1      single     1.354   0.011
UNL          C2          C3      double     1.369   0.016
UNL          C3           N      single     1.384   0.010
UNL          C4          C3      single     1.422   0.010
UNL          C4           O      double     1.224   0.019
UNL          N2          C4      single     1.405   0.012
UNL          N2          C5      single     1.469   0.010
UNL          C6          N2      single     1.390   0.013
UNL          O1          C6      double     1.218   0.012
UNL          C6          N3      single     1.376   0.011
UNL          N3          C2      single     1.376   0.011
UNL          N3          C7      single     1.465   0.010
UNL           C           H      single     0.969   0.020
UNL           C          H1      single     0.969   0.020
UNL           C          H2      single     0.969   0.020
UNL          C1          H3      single     0.943   0.020
UNL          C5          H4      single     0.969   0.020
UNL          C5          H5      single     0.969   0.020
UNL          C5          H6      single     0.969   0.020
UNL          C7          H7      single     0.969   0.020
UNL          C7          H8      single     0.969   0.020
UNL          C7          H9      single     0.969   0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
UNL          N2          C6          O1     121.264    0.81
UNL          N2          C6          N3     117.274    0.65
UNL          O1          C6          N3     121.462    0.79
UNL          C6          N3          C2     119.603    0.59
UNL          C6          N3          C7     119.286    1.20
UNL          C2          N3          C7     121.111    0.73
UNL          N3          C7          H7     109.463    1.35
UNL          N3          C7          H8     109.463    1.35
UNL          N3          C7          H9     109.463    1.35
UNL          H7          C7          H8     109.417    2.07
UNL          H7          C7          H9     109.417    2.07
UNL          H8          C7          H9     109.417    2.07
UNL          N1          C2          C3     111.816    0.88
UNL          N1          C2          N3     126.343    0.73
UNL          C3          C2          N3     121.841    0.76
UNL          C1          N1          C2     103.393    0.94
UNL           N          C1          N1     113.403    0.89
UNL           N          C1          H3     123.065    1.13
UNL          N1          C1          H3     123.531    1.23
UNL          C4          N2          C5     117.108    1.11
UNL          C4          N2          C6     126.504    0.90
UNL          C5          N2          C6     116.389    0.81
UNL          N2          C5          H4     109.479    1.25
UNL          N2          C5          H5     109.479    1.25
UNL          N2          C5          H6     109.479    1.25
UNL          H4          C5          H5     109.417    2.07
UNL          H4          C5          H6     109.417    2.07
UNL          H5          C5          H6     109.417    2.07
UNL          C3          C4           O     126.725    1.16
UNL          C3          C4          N2     111.977    1.50
UNL           O          C4          N2     121.297    1.02
UNL          C2          C3           N     105.605    0.90
UNL          C2          C3          C4     122.801    1.13
UNL           N          C3          C4     131.594    0.88
UNL           C           N          C1     126.918    0.90
UNL           C           N          C3     127.300    1.10
UNL          C1           N          C3     105.782    0.93
UNL           N           C           H     109.514    1.29
UNL           N           C          H1     109.514    1.29
UNL           N           C          H2     109.514    1.29
UNL           H           C          H1     109.415    1.67
UNL           H           C          H2     109.415    1.67
UNL          H1           C          H2     109.415    1.67
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
UNL     P_sp2_sp2_1          C3          C2          N1          C1       0.000   10.00     2
UNL     P_sp2_sp2_2          N3          C2          N1          C1     180.000   10.00     2
UNL     P_sp2_sp2_3           N          C1          N1          C2       0.000   10.00     2
UNL     P_sp2_sp2_4          H3          C1          N1          C2     180.000   10.00     2
UNL     P_sp2_sp2_5          N1          C1           N          C3       0.000   10.00     2
UNL     P_sp2_sp2_6          N1          C1           N           C     180.000   10.00     2
UNL     P_sp2_sp2_7          H3          C1           N          C3     180.000   10.00     2
UNL     P_sp2_sp2_8          H3          C1           N           C       0.000   10.00     2
UNL     P_sp2_sp2_9          C2          C3           N          C1       0.000   10.00     2
UNL    P_sp2_sp2_10          C2          C3           N           C     180.000   10.00     2
UNL    P_sp2_sp2_11          C4          C3           N          C1     180.000   10.00     2
UNL    P_sp2_sp2_12          C4          C3           N           C       0.000   10.00     2
UNL    P_sp2_sp2_13          N3          C6          N2          C4       0.000   10.00     2
UNL    P_sp2_sp2_14          N3          C6          N2          C5     180.000   10.00     2
UNL    P_sp2_sp2_15          O1          C6          N2          C4     180.000   10.00     2
UNL    P_sp2_sp2_16          O1          C6          N2          C5       0.000   10.00     2
UNL    P_sp2_sp2_17          C3          C4          N2          C6       0.000   10.00     2
UNL    P_sp2_sp2_18          C3          C4          N2          C5     180.000   10.00     2
UNL    P_sp2_sp2_19           O          C4          N2          C6     180.000   10.00     2
UNL    P_sp2_sp2_20           O          C4          N2          C5       0.000   10.00     2
UNL    P_sp2_sp2_21          C2          C3          C4          N2       0.000   10.00     2
UNL    P_sp2_sp2_22          C2          C3          C4           O     180.000   10.00     2
UNL    P_sp2_sp2_23           N          C3          C4          N2     180.000   10.00     2
UNL    P_sp2_sp2_24           N          C3          C4           O       0.000   10.00     2
UNL    P_sp2_sp2_25          N1          C2          C3           N       0.000   10.00     2
UNL    P_sp2_sp2_26          N1          C2          C3          C4     180.000   10.00     2
UNL    P_sp2_sp2_27          N3          C2          C3           N     180.000   10.00     2
UNL    P_sp2_sp2_28          N3          C2          C3          C4       0.000   10.00     2
UNL    P_sp2_sp2_29          C3          C2          N3          C6       0.000   10.00     2
UNL    P_sp2_sp2_30          C3          C2          N3          C7     180.000   10.00     2
UNL    P_sp2_sp2_31          N1          C2          N3          C6     180.000   10.00     2
UNL    P_sp2_sp2_32          N1          C2          N3          C7       0.000   10.00     2
UNL       sp2_sp3_1          C1           N           C           H     150.000   10.00     6
UNL       sp2_sp3_2          C1           N           C          H1     -90.000   10.00     6
UNL       sp2_sp3_3          C1           N           C          H2      30.000   10.00     6
UNL       sp2_sp3_4          C3           N           C           H     -30.000   10.00     6
UNL       sp2_sp3_5          C3           N           C          H1      90.000   10.00     6
UNL       sp2_sp3_6          C3           N           C          H2    -150.000   10.00     6
UNL       sp2_sp3_7          C4          N2          C5          H4     150.000   10.00     6
UNL       sp2_sp3_8          C4          N2          C5          H5     -90.000   10.00     6
UNL       sp2_sp3_9          C4          N2          C5          H6      30.000   10.00     6
UNL      sp2_sp3_10          C6          N2          C5          H4     -30.000   10.00     6
UNL      sp2_sp3_11          C6          N2          C5          H5      90.000   10.00     6
UNL      sp2_sp3_12          C6          N2          C5          H6    -150.000   10.00     6
UNL    P_sp2_sp2_33          N2          C6          N3          C2       0.000   10.00     2
UNL    P_sp2_sp2_34          N2          C6          N3          C7     180.000   10.00     2
UNL    P_sp2_sp2_35          O1          C6          N3          C2     180.000   10.00     2
UNL    P_sp2_sp2_36          O1          C6          N3          C7       0.000   10.00     2
UNL      sp2_sp3_13          C6          N3          C7          H7     150.000   10.00     6
UNL      sp2_sp3_14          C6          N3          C7          H8     -90.000   10.00     6
UNL      sp2_sp3_15          C6          N3          C7          H9      30.000   10.00     6
UNL      sp2_sp3_16          C2          N3          C7          H7     -30.000   10.00     6
UNL      sp2_sp3_17          C2          N3          C7          H8      90.000   10.00     6
UNL      sp2_sp3_18          C2          N3          C7          H9    -150.000   10.00     6
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
UNL    plan-1           C   0.020
UNL    plan-1          C1   0.020
UNL    plan-1          C2   0.020
UNL    plan-1          C3   0.020
UNL    plan-1          C4   0.020
UNL    plan-1          C5   0.020
UNL    plan-1          C6   0.020
UNL    plan-1          C7   0.020
UNL    plan-1          H3   0.020
UNL    plan-1           N   0.020
UNL    plan-1          N1   0.020
UNL    plan-1          N2   0.020
UNL    plan-1          N3   0.020
UNL    plan-1           O   0.020
UNL    plan-1          O1   0.020
