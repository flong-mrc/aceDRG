# ======================================================================
# acedrg.tcl --
#
# CCP4Interface
#
# ======================================================================= 
#------------------------------------------------------------------------------
proc acedrg_prereq { } {
#------------------------------------------------------------------------------
  if { ![file exists [FindExecutable acedrg]] } {
    return 0
  }
  return 1
}

#----------------------------------------------------------------------
proc acedrg_setup {typedefVar arrayname } {
#----------------------------------------------------------------------
  upvar #0  $typedefVar typedef
  upvar #0  $arrayname  array

  set typedef(_acedrg_mode) { menu { "Generate Ligands/Monomers" }
                              { LIG_GEN } }

  set typedef(_liggen_mode) { menu {"MMCIF" 
                                    "SMILES"
                                    "MDL(MOLFILE)/SDF"
                                    "(SYBYL)MOL2"}
                              { MMC SMI  MOL MOL2} }

  DefineMenu _input_file_type [list "SMILES file " "MMCIF file " "MOL file " "(small mol) CIF file " ] \
                              [list SMILEIN MMCIIN MOLIN CIFIN  ] 

  return 1
}

#---------------------------------------------------------------------
proc acedrg_run { arrayname } {
#--------------------------------------------------------------------- 
  upvar #0 $arrayname array
 
  switch [GetValue $arrayname ACEDRG_MODE] \
      LIG_GEN {  
               # check if input files exist
               if { [file exists [GetFullFileName $array(SMILESIN) $array(DIR_SMILESIN) ] ] } {
                    set array(INPUT_FILES) " SMILESIN "
               }

               if { [file exists [GetFullFileName $array(MMCIFIN) $array(DIR_MMCIFIN) ] ] } {
                    set array(INPUT_FILES) " MMCIFIN "
               }

               if { [file exists [GetFullFileName $array(MOLIN) $array(DIR_MOLIN) ] ] } {
                    set array(INPUT_FILES) " MOLIN "
               }

               if { [file exists [GetFullFileName $array(CIFIN) $array(DIR_CIFIN) ] ] } {
                    set array(INPUT_FILES) " CIFIN "
               }

               if { [file exists [GetFullFileName $array(MOL2IN) $array(DIR_MOL2IN) ] ] } {
                    set array(INPUT_FILES) " MOL2IN "
               }

               if { $array(XYZOUT) != "" } {
                    set array(OUTPUT_FILES) " XYZOUT "
               } 
 
               if { $array(MMCIFOUT) != "" } {
                   append array(OUTPUT_FILES) " MMCIFOUT "
               }

               set array(OUT_ROOTDIR) [GetDefaultDirPath]
           }

  return 1

}

#-----------------------------------------------------------------------
proc acedrg_task_window {arrayname} {
#-----------------------------------------------------------------------
  upvar #0 $arrayname array

  if { [CreateTaskWindow  $arrayname \
        "Run ACEDRG " "ACEDRG" \
        ] == 0 } return

#=PROTOCOL==============================================================

  OpenFolder protocol

  CreateTitleLine line TITLE

  CreateLine line \
    message "Do you want to generate ligands/monomers or derive molecules/bonds/angles" \
    help acedrg_mode \
    label "Select a work mode for acedrg " \
    widget ACEDRG_MODE 

#=FILES================================================================

  # LIGAND/MONOMER generator mode
  OpenSubFrame frame -toggle_display ACEDRG_MODE open [list LIG_GEN ]
 
  CreateLine line \
    message "What is the format of your input file, SMILES, mmCif or Mol " \
    help liggen_mode \
    label "Select a format for your input file " \
    widget LIGGEN_MODE  \
    toggle_display ACEDRG_MODE open [list LIG_GEN ]

  CloseSubFrame

  OpenFolder file

  OpenSubFrame frame -toggle_display LIGGEN_MODE open [list SMI ]
 
  # input file format: SMILES 

  CreateInputFileLine line \
        "Enter input SMILES file name " \
        "Input SMILES file" \
         SMILESIN DIR_SMILESIN \
         -fileout MMCIFOUT DIR_MMCIFOUT "_acedrg" \
         -fileout XYZOUT DIR_XYZOUT "_acedrg" \
         toggle_display LIGGEN_MODE open [list SMI ] 
  
  CloseSubFrame

  OpenSubFrame frame -toggle_display LIGGEN_MODE open [list  MMC ]
 
  # input file format: mmCIF

  CreateInputFileLine line \
        "Enter input mmCif file name " \
        "Input mmCif file" \
         MMCIFIN DIR_MMCIFIN \
         -fileout MMCIFOUT DIR_MMCIFOUT "_acedrg" \
         -fileout XYZOUT DIR_XYZOUT "_acedrg" \
         toggle_display LIGGEN_MODE open [list  MMC ] 
    
  CloseSubFrame

  OpenSubFrame frame -toggle_display LIGGEN_MODE open [list MOL ]
 
  # input file format: MOL

  CreateInputFileLine line \
        "Enter input Mol file name " \
        "Input Mol file" \
         MOLIN DIR_MOLIN \
         -fileout MMCIFOUT DIR_MMCIFOUT "_acedrg" \
         -fileout XYZOUT DIR_XYZOUT "_acedrg" \
         toggle_display LIGGEN_MODE open [list MOL ]    
  
  CloseSubFrame

  #OpenSubFrame frame -toggle_display LIGGEN_MODE open [list MOL2 ]
 
  # input file format: SYBYL_MOL2

  CreateInputFileLine line \
        "Enter input Mol2 file name " \
        "Input Mol2 file" \
         MOL2IN DIR_MOL2IN \
         -fileout MMCIFOUT DIR_MMCIFOUT "_acedrg" \
         -fileout XYZOUT DIR_XYZOUT "_acedrg" \
         toggle_display LIGGEN_MODE open [list MOL2 ]    
  
  CloseSubFrame

  OpenSubFrame frame -toggle_display ACEDRG_MODE open [list LIG_GEN ]

  CreateLine line \
    message "Enter the ligand/monomer three-letter ID " \
    label "Ligand/monomer three-letter ID " \
    widget MONOM_ID -oblig
 
  #CreateOutputFileLine line \
  #  "Enter root name for output files and tables "  \
  #  "OUTPUT ROOT NAME " \
  #  OUT_ROOT DIR_OUT_ROOT 
  CreateOutputFileLine line \
    "Enter file name for the output dictionary file (mmcif) "  \
    "Solution mmCif " \
    MMCIFOUT DIR_MMCIFOUT \
    toggle_display ACEDRG_MODE open [list LIG_GEN ]

  CreateOutputFileLine line \
    "Enter file name for the output coordinate file (PDB) "  \
    "Solution PDB " \
    XYZOUT DIR_XYZOUT \
    toggle_display ACEDRG_MODE open [list LIG_GEN ]

  CloseSubFrame

  #-------------------------------- Required parameters

}







