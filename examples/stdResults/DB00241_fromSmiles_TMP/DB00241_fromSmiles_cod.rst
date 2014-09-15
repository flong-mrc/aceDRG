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
UNL          O2     O     O       0       0.000       0.000       0.000
UNL         C10     C   CR6   0.000       0.380       0.146       1.167
UNL          C4     C    CT   0.000      -0.418      -0.885       1.879
UNL          C5     C   CH2   0.000       0.523      -1.896       1.212
UNL          H9     H     H   0.000       1.396      -1.283       0.980
UNL         H10     H     H   0.000      -0.009      -2.134       0.288
UNL          C6     C    C1   0.000       0.959      -3.153       1.873
UNL         H11     H     H   0.000       0.250      -3.717       2.455
UNL          C7     C    C2   0.000       2.207      -3.597       1.760
UNL         H13     H     H   0.000       2.929      -3.044       1.179
UNL         H12     H     H   0.000       2.503      -4.513       2.247
UNL          C3     C   CH2   0.000      -0.249      -0.497       3.311
UNL          H7     H     H   0.000      -1.212      -0.107       3.647
UNL          H8     H     H   0.000       0.494       0.304       3.334
UNL          C1     C   CH1   0.000       0.184      -1.605       4.214
UNL          H3     H     H   0.000      -0.499      -2.459       4.104
UNL          C2     C   CH3   0.000       1.605      -2.092       4.213
UNL          H6     H     H   0.000       1.850      -2.446       3.248
UNL          H5     H     H   0.000       2.246      -1.293       4.474
UNL          H4     H     H   0.000       1.705      -2.874       4.917
UNL           C     C   CH3   0.000       0.112      -1.075       5.636
UNL          H2     H     H   0.000       0.756      -0.240       5.733
UNL          H1     H     H   0.000      -0.882      -0.782       5.852
UNL           H     H     H   0.000       0.413      -1.834       6.310
UNL          N1     N  NR16   0.000       0.306       1.488       1.348
UNL         H15     H     H   0.000       0.980       2.060       1.897
UNL          C9     C   CR6   0.000      -0.770       1.961       0.709
UNL          O1     O     O   0.000      -0.613       2.630      -0.337
UNL           N     N  NR16   0.000      -1.636       0.973       0.550
UNL         H14     H     H   0.000      -2.421       1.150      -0.109
UNL          C8     C   CR6   0.000      -1.618      -0.272       1.162
UNL           O     O     O   0.000      -2.678      -0.639       1.669
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
UNL          C1           C      SINGLE     1.516   0.020
UNL          C1          C2      SINGLE     1.516   0.020
UNL          C3          C1      SINGLE     1.517   0.018
UNL          C4          C3      SINGLE     1.548   0.010
UNL          C4          C5      SINGLE     1.551   0.010
UNL          C5          C6      SINGLE     1.499   0.011
UNL          C6          C7      DOUBLE     1.291   0.020
UNL          C4          C8      SINGLE     1.520   0.010
UNL          C8           O      DOUBLE     1.218   0.011
UNL           N          C8      SINGLE     1.370   0.010
UNL          C9           N      SINGLE     1.372   0.011
UNL          C9          O1      DOUBLE     1.225   0.013
UNL          N1          C9      SINGLE     1.372   0.011
UNL         C10          N1      SINGLE     1.370   0.010
UNL         C10          C4      SINGLE     1.520   0.010
UNL          O2         C10      DOUBLE     1.218   0.011
UNL           C           H      SINGLE     0.974   0.020
UNL           C          H1      SINGLE     0.974   0.020
UNL           C          H2      SINGLE     0.974   0.020
UNL          C1          H3      SINGLE     0.987   0.017
UNL          C2          H4      SINGLE     0.974   0.020
UNL          C2          H5      SINGLE     0.974   0.020
UNL          C2          H6      SINGLE     0.974   0.020
UNL          C3          H7      SINGLE     0.980   0.019
UNL          C3          H8      SINGLE     0.980   0.019
UNL          C5          H9      SINGLE     0.982   0.014
UNL          C5         H10      SINGLE     0.982   0.014
UNL          C6         H11      SINGLE     0.948   0.020
UNL          C7         H12      SINGLE     0.949   0.020
UNL          C7         H13      SINGLE     0.949   0.020
UNL           N         H14      SINGLE     0.877   0.020
UNL          N1         H15      SINGLE     0.877   0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
UNL          N1         C10          C4     118.120    0.73
UNL          N1         C10          O2     120.198    1.22
UNL          C4         C10          O2     121.682    0.95
UNL          C3          C4          C5     109.313    1.94
UNL          C3          C4          C8     108.254    1.00
UNL          C3          C4         C10     108.254    1.00
UNL          C5          C4          C8     108.346    0.82
UNL          C5          C4         C10     108.346    0.82
UNL          C8          C4         C10     114.092    0.38
UNL          C4          C5          C6     113.426    1.02
UNL          C4          C5          H9     108.862    0.22
UNL          C4          C5         H10     108.862    0.22
UNL          C6          C5          H9     108.956    0.79
UNL          C6          C5         H10     108.956    0.79
UNL          H9          C5         H10     107.646    0.56
UNL          C5          C6          C7     125.278    2.02
UNL          C5          C6         H11     117.905    2.92
UNL          C7          C6         H11     116.817    2.73
UNL          C6          C7         H12     120.031    2.51
UNL          C6          C7         H13     120.031    2.51
UNL         H12          C7         H13     119.938    2.37
UNL          C1          C3          C4     112.921    2.66
UNL          C1          C3          H7     107.989    0.99
UNL          C1          C3          H8     107.989    0.99
UNL          C4          C3          H7     108.589    0.39
UNL          C4          C3          H8     108.589    0.39
UNL          H7          C3          H8     107.577    1.40
UNL           C          C1          C2     110.499    2.07
UNL           C          C1          C3     110.888    1.67
UNL           C          C1          H3     107.911    2.01
UNL          C2          C1          C3     110.888    1.67
UNL          C2          C1          H3     107.911    2.01
UNL          C3          C1          H3     107.533    0.87
UNL          C1          C2          H4     109.485    0.79
UNL          C1          C2          H5     109.485    0.79
UNL          C1          C2          H6     109.485    0.79
UNL          H4          C2          H5     109.300    1.49
UNL          H4          C2          H6     109.300    1.49
UNL          H5          C2          H6     109.300    1.49
UNL          C1           C           H     109.485    0.79
UNL          C1           C          H1     109.485    0.79
UNL          C1           C          H2     109.485    0.79
UNL           H           C          H1     109.300    1.49
UNL           H           C          H2     109.300    1.49
UNL          H1           C          H2     109.300    1.49
UNL          C9          N1         C10     125.880    1.92
UNL          C9          N1         H15     116.838    2.16
UNL         C10          N1         H15     117.283    2.19
UNL           N          C9          O1     121.882    1.20
UNL           N          C9          N1     116.236    1.09
UNL          O1          C9          N1     121.882    1.20
UNL          C8           N          C9     125.880    1.92
UNL          C8           N         H14     117.283    2.19
UNL          C9           N         H14     116.838    2.16
UNL          C4          C8           O     121.682    0.95
UNL          C4          C8           N     118.120    0.73
UNL           O          C8           N     120.198    1.22
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
UNL       sp2_sp2_1          C4         C10          N1          C9       0.000   10.00     2
UNL       sp2_sp2_2          C4         C10          N1         H15     180.000   10.00     2
UNL       sp2_sp2_3          O2         C10          N1          C9     180.000   10.00     2
UNL       sp2_sp2_4          O2         C10          N1         H15       0.000   10.00     2
UNL       sp2_sp2_5           N          C9          N1         C10       0.000   10.00     2
UNL       sp2_sp2_6           N          C9          N1         H15     180.000   10.00     2
UNL       sp2_sp2_7          O1          C9          N1         C10     180.000   10.00     2
UNL       sp2_sp2_8          O1          C9          N1         H15       0.000   10.00     2
UNL       sp2_sp2_9          N1          C9           N          C8       0.000   10.00     2
UNL      sp2_sp2_10          N1          C9           N         H14     180.000   10.00     2
UNL      sp2_sp2_11          O1          C9           N          C8     180.000   10.00     2
UNL      sp2_sp2_12          O1          C9           N         H14       0.000   10.00     2
UNL      sp2_sp2_13          C4          C8           N          C9       0.000   10.00     2
UNL      sp2_sp2_14          C4          C8           N         H14     180.000   10.00     2
UNL      sp2_sp2_15           O          C8           N          C9     180.000   10.00     2
UNL      sp2_sp2_16           O          C8           N         H14       0.000   10.00     2
UNL       sp2_sp3_1           N          C8          C4         C10       0.000   10.00     6
UNL       sp2_sp3_2           N          C8          C4          C5     120.000   10.00     6
UNL       sp2_sp3_3           N          C8          C4          C3    -120.000   10.00     6
UNL       sp2_sp3_4           O          C8          C4         C10     180.000   10.00     6
UNL       sp2_sp3_5           O          C8          C4          C5     -60.000   10.00     6
UNL       sp2_sp3_6           O          C8          C4          C3      60.000   10.00     6
UNL       sp3_sp3_1           H           C          C1          C3     180.000   10.00     3
UNL       sp3_sp3_2           H           C          C1          C2     -60.000   10.00     3
UNL       sp3_sp3_3           H           C          C1          H3      60.000   10.00     3
UNL       sp3_sp3_4          H1           C          C1          C3      60.000   10.00     3
UNL       sp3_sp3_5          H1           C          C1          C2     180.000   10.00     3
UNL       sp3_sp3_6          H1           C          C1          H3     -60.000   10.00     3
UNL       sp3_sp3_7          H2           C          C1          C3     -60.000   10.00     3
UNL       sp3_sp3_8          H2           C          C1          C2      60.000   10.00     3
UNL       sp3_sp3_9          H2           C          C1          H3     180.000   10.00     3
UNL      sp3_sp3_10           C          C1          C2          H4     180.000   10.00     3
UNL      sp3_sp3_11           C          C1          C2          H5     -60.000   10.00     3
UNL      sp3_sp3_12           C          C1          C2          H6      60.000   10.00     3
UNL      sp3_sp3_13          C3          C1          C2          H4      60.000   10.00     3
UNL      sp3_sp3_14          C3          C1          C2          H5     180.000   10.00     3
UNL      sp3_sp3_15          C3          C1          C2          H6     -60.000   10.00     3
UNL      sp3_sp3_16          H3          C1          C2          H4     -60.000   10.00     3
UNL      sp3_sp3_17          H3          C1          C2          H5      60.000   10.00     3
UNL      sp3_sp3_18          H3          C1          C2          H6     180.000   10.00     3
UNL      sp3_sp3_19          C2          C1          C3          C4     180.000   10.00     3
UNL      sp3_sp3_20          C2          C1          C3          H7     -60.000   10.00     3
UNL      sp3_sp3_21          C2          C1          C3          H8      60.000   10.00     3
UNL      sp3_sp3_22           C          C1          C3          C4      60.000   10.00     3
UNL      sp3_sp3_23           C          C1          C3          H7     180.000   10.00     3
UNL      sp3_sp3_24           C          C1          C3          H8     -60.000   10.00     3
UNL      sp3_sp3_25          H3          C1          C3          C4     -60.000   10.00     3
UNL      sp3_sp3_26          H3          C1          C3          H7      60.000   10.00     3
UNL      sp3_sp3_27          H3          C1          C3          H8     180.000   10.00     3
UNL      sp3_sp3_28          C1          C3          C4          C8     180.000   10.00     3
UNL      sp3_sp3_29          C1          C3          C4          C5     -60.000   10.00     3
UNL      sp3_sp3_30          C1          C3          C4         C10      60.000   10.00     3
UNL      sp3_sp3_31          H7          C3          C4          C8      60.000   10.00     3
UNL      sp3_sp3_32          H7          C3          C4          C5     180.000   10.00     3
UNL      sp3_sp3_33          H7          C3          C4         C10     -60.000   10.00     3
UNL      sp3_sp3_34          H8          C3          C4          C8     -60.000   10.00     3
UNL      sp3_sp3_35          H8          C3          C4          C5      60.000   10.00     3
UNL      sp3_sp3_36          H8          C3          C4         C10     180.000   10.00     3
UNL      sp3_sp3_37          C3          C4          C5          C6     180.000   10.00     3
UNL      sp3_sp3_38          C3          C4          C5          H9     -60.000   10.00     3
UNL      sp3_sp3_39          C3          C4          C5         H10      60.000   10.00     3
UNL      sp3_sp3_40          C8          C4          C5          C6      60.000   10.00     3
UNL      sp3_sp3_41          C8          C4          C5          H9     180.000   10.00     3
UNL      sp3_sp3_42          C8          C4          C5         H10     -60.000   10.00     3
UNL      sp3_sp3_43         C10          C4          C5          C6     -60.000   10.00     3
UNL      sp3_sp3_44         C10          C4          C5          H9      60.000   10.00     3
UNL      sp3_sp3_45         C10          C4          C5         H10     180.000   10.00     3
UNL       sp2_sp3_7          C7          C6          C5          H9       0.000   10.00     6
UNL       sp2_sp3_8          C7          C6          C5          C4     120.000   10.00     6
UNL       sp2_sp3_9          C7          C6          C5         H10    -120.000   10.00     6
UNL      sp2_sp3_10         H11          C6          C5          H9     180.000   10.00     6
UNL      sp2_sp3_11         H11          C6          C5          C4     -60.000   10.00     6
UNL      sp2_sp3_12         H11          C6          C5         H10      60.000   10.00     6
UNL      sp2_sp2_17          C5          C6          C7         H12     180.000   10.00     2
UNL      sp2_sp2_18          C5          C6          C7         H13       0.000   10.00     2
UNL      sp2_sp2_19         H11          C6          C7         H12       0.000   10.00     2
UNL      sp2_sp2_20         H11          C6          C7         H13     180.000   10.00     2
UNL      sp2_sp3_13          N1         C10          C4          C8       0.000   10.00     6
UNL      sp2_sp3_14          N1         C10          C4          C3     120.000   10.00     6
UNL      sp2_sp3_15          N1         C10          C4          C5    -120.000   10.00     6
UNL      sp2_sp3_16          O2         C10          C4          C8     180.000   10.00     6
UNL      sp2_sp3_17          O2         C10          C4          C3     -60.000   10.00     6
UNL      sp2_sp3_18          O2         C10          C4          C5      60.000   10.00     6
loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
UNL    chir_01    C1    C    C2    C3    both
UNL    chir_02    C4    C3    C5    C8    both
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
UNL    plan-1         C10   0.020
UNL    plan-1          C4   0.020
UNL    plan-1          N1   0.020
UNL    plan-1          O2   0.020
UNL    plan-2          C5   0.020
UNL    plan-2          C6   0.020
UNL    plan-2          C7   0.020
UNL    plan-2         H11   0.020
UNL    plan-3          C6   0.020
UNL    plan-3          C7   0.020
UNL    plan-3         H12   0.020
UNL    plan-3         H13   0.020
UNL    plan-4         C10   0.020
UNL    plan-4          C9   0.020
UNL    plan-4         H15   0.020
UNL    plan-4          N1   0.020
UNL    plan-5          C9   0.020
UNL    plan-5           N   0.020
UNL    plan-5          N1   0.020
UNL    plan-5          O1   0.020
UNL    plan-6          C8   0.020
UNL    plan-6          C9   0.020
UNL    plan-6         H14   0.020
UNL    plan-6           N   0.020
UNL    plan-7          C4   0.020
UNL    plan-7          C8   0.020
UNL    plan-7           N   0.020
UNL    plan-7           O   0.020
