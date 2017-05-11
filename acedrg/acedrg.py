#!/usr/bin/env  ccp4-python
# Python script
#
#
#     Copyright (C) 2014 --- 2019 Fei Long,  G. Murshudov
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Library.
#
#====================================================================
## The date of last modification: 21/07/2016
#

import os,os.path,sys
import platform
import glob,shutil
import re,string
from optparse import OptionParser 
import time
import math
import select
import random

from rdkit      import rdBase

from rdkit      import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops
from rdkit.Chem import Pharm3D 
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Geometry import rdGeometry 

from exebase    import CExeCode

from covLink    import CovLink
from covLink    import CovLinkGenerator

from utility    import listComp
from utility    import listComp2
from utility    import listCompDes
from utility    import listCompAcd
from utility    import setBoolDict
from utility    import splitLineSpa

if os.name != 'nt':
    import fcntl
import signal

class Acedrg(CExeCode ):

    def __init__( self, t_argvs):

        if len(sys.argv)==1:
            print "Look for help: %s -h "%(os.path.basename(sys.argv[0]))

        self.versionInfo       = {}

        self.errMessage       = []
        self.errLevel         = 0

        self.scrDir           = ""
        self.outRoot          = ""
        self.baseRoot         = ""

        self.geneInFileName   = ""
        self.geneInFileType   = ""
        self.inStdCifName     = ""
        self.inStdCifDir      = ""
        self.inMmCifName      = ""
        self.inSmiName        = ""
        self.inMdlName        = ""
        self.inPdbName        = ""
        self.inLigandPdbName  = ""
        self.iniLigandPdbName = ""
        self.inSDFName        = ""
        self.inMol2Name       = ""
        self.outRstCifName    = ""
        self.outRstPdb        = ""
        self.outAtmTypeName   = ""
        self.monomRoot        = ""


        self.acedrgTables     = ""
        self.libmol           = ""
        self.libmolLogName    = ""
        self.libmolAT1        = ""
        self.libmolAT2        = ""
        self.libmolMatched    = ""
        self.libcheck         = ""
        self.libcheckLogName  = ""
        self.refmac           = ""
        self.refmacLogName    = ""

        self.refmacXYSList    = {}
        self.refmacMinFValueList = []
        #self.refmacMinFValueList["value"] =100000.00
        #self.refmacMinFValueList["fileName"] =""
        

        self.linkInstructions = ""
        self.linkInsMap       = {}
        
        self.funcGroupTable   = ""

        self.numConformers    = 1
        self.useCifCoords     = False

        self.inputPara        = {}
        self.inputPara["PH"]  = [False, 0.0]

        self.elemSet          = []

        self.workMode         = 0
        self.useExistCoords   = False

        self.molGen           = False
        self.repCrds          = False

        self.HMO              = False
        
        self.testMode         = False

        self.outCifGlobSect    = []

        self.allBondsAndAngles = {}

        self.allBondsAndAngles["atomClasses"] = {}
        self.allBondsAndAngles["bonds"]       = {}
        self.allBondsAndAngles["angles"]      = {}


        self.acedrg = os.path.abspath(sys.argv[0])
        self.acedrgDir = os.path.dirname(os.path.dirname(self.acedrg))
        #print self.acedrgDir
        inputOptionsP         = self.InputParser(t_argvs) 

        #if inputOptionsP.geneInFileName:
        #    self.checInputFormat()          

        self.runExitCode      = 0 

        self.setWorkMode(inputOptionsP)
        self.setInputProcPara(inputOptionsP)

        self.checkDependency()
        self.checkVersionInfo()
        self.setOutCifGlobSec()

        if os.path.isfile(self.funcGroupTable):
            self.rdKit = AcedrgRDKit(self.funcGroupTable)
        else:
            self.rdKit = AcedrgRDKit()
        self.rdKit.setProcPara(self.inputPara)

        #print "input RDKit: userExistCoords ? ", self.rdKit.useExistCoords
        #print "input RDKit: number of optimization iters ", self.rdKit.numRDKitOptmSteps
        #print "input RDKit: number of initial conformers ", self.rdKit.numInitConformers
        #print "input RDKit: number of output conformers ", self.rdKit.numConformers

        self.fileConv         = FileTransformer()   
            
        self.chemCheck        = ChemCheck() 
        
        self.initMmcifMolMap  = {}

        self.cLinkMap         = {}

        #self.execute()  
        self.executeWithRDKit()  
        
    def InputParser(self, t_argvs):

        usage = "\n\
        (1) Input a mmCif file to get a ligand/Dictionary/restraint files (mmCif fomat) \n\
            acedrg -c your_mmCif_file -r your_monomer_name(optional) -o output_name_root(optional) \n\n\
        (2) Input a SMILES string in a command line (the SMILES string is a pair of quotation marks) or \n \
            a file containing a SMILES string  to get a ligand/Dictionary/restraint files (mmCif fomat) \n\
            acedrg -i your_SMILES_file -r your_monomer_name(optional) -o output_name_root(optional) \n\n\
        (3) Input a MOL file to get a ligand/Dictionary/restraint files (mmCif fomat) \n\
            acedrg -m your_mol_file -r your_monomer_name(optional) -o output_name_root(optional) \n\n\
        (4) Input a MOL2 file to get a ligand/Dictionary/restraint files (mmCif fomat) \n\
            acedrg -g your_mol_file -r your_monomer_name(optional) -o output_name_root(optional) \n\n"
    
        """
        (4) Input a small molecule CIF file to get two files describing your molecules, atom types, \n\
            bond lengths and angles (text format)\n\
            acedrg -e -b your_CIF_file -r your_monomer_name(optional) -o output_name_root(optional) \n\n\
         5) Input a directory containing a number of  molecule CIF files to get files describing your\n\
            molecules, atom types, bond lengths and angles, and their statistics \n\
            acedrg -e -d your_dir_of_CIF_files -o output_name_root(optional) \n\n"
         """

        self.inputParser = OptionParser(usage=usage)

        # Options 
        # input file format 
        self.inputParser.add_option("-b",  "--stdcif", dest="inStdCifName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input small molecule CIF File containing coordinates and crystal information")

        self.inputParser.add_option("-c",  "--mmcif", dest="inMmCifName", metavar="FILE", 
                                    action="store", type="string", help="Input MMCIF File containing coordinates and bonds")

        self.inputParser.add_option("-d",  "--cifdir", dest="inStdCifDir", metavar="FILE", 
                                    action="store", type="string", 
                                    help="An input directory to store a group of small molecule CIF Files containing coordinates and crystal information")

        self.inputParser.add_option("-e",  "--molgen", dest="molGen",  
                                    action="store_true",  default=False,  
                                    help="The option when the user want to generate molecules and values of the associated bonds and angles in a cif file")

        self.inputParser.add_option("-f",  "--infile", dest="geneInFileName",  
                                    action="store",  type="string",   
                                    help="user's input file, regardless what format it is")

        self.inputParser.add_option("-g",  "--mol2", dest="inMol2Name", metavar="FILE", 
                                    action="store", type="string", help="Input SYBL_MOL2 File containing coordinates and bonds")

        self.inputParser.add_option("-i",  "--smi", dest="inSmiName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input  File containing a SMILE string")

        self.inputParser.add_option("-j",  "--numInitConf", dest="numInitConformers", 
                                    action="store", type="int", 
                                    help="Number of initial conformers to try")

        self.inputParser.add_option("-k",  "--multiconf", dest="numConformers", 
                                    action="store", type="int", 
                                    help="Number of final output conformers")

        self.inputParser.add_option("-l",  "--numOptmStep", dest="numRDKitOptmSteps", 
                                    action="store", type="int", 
                                    help="Number of RDKit optimization steps requested")

        self.inputParser.add_option("-m",  "--mol", dest="inMdlName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File of MOL format containing coordinates and bonds")

        self.inputParser.add_option("-n",  "--typeOut", dest="typeOut",  
                                    action="store_true",  default=False,  
                                    help="The option when the user want to output two kinds of atom types (CCP4 and Acedrg) only")

        self.inputParser.add_option("-o",  "--out", dest="outRoot", 
                                    action="store", type="string", 
                                    help="A name root that users want their output files called(without extension)")

        self.inputParser.add_option("-x",  "--pdb", dest="inLigandPdbName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File of PDB format containing coordinates of the ligand")

        self.inputParser.add_option("-p",  "--coords", dest="useExistCoords",  
                                    action="store_true",  default=False,
                                    help="Using existing coordinates in the input file")

        self.inputParser.add_option("-r",  "--res", dest="monomRoot",  
                                    action="store", type="string", 
                                    help="The name of the chemical components users want to put into output files(e.g. PDB or MMCIF)")

        self.inputParser.add_option("-s",  "--sdf", dest="inSdfName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File of SDF format containing coordinates and bonds")

        self.inputParser.add_option("-t",  "--tab", dest="acedrgTables", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input path that stores all bond and angle tables (if no input, default CCP4 location will be used)")

        self.inputParser.add_option("-u",  "--hmo", dest="HMO",  
                                    action="store_true",  default=False,
                                    help="Calculate bond-orders and charges using HMO ")

        self.inputParser.add_option("-v",  "--version", dest="versionInfo",  
                                    action="store_true",  default=False,
                                    help="The option for checking version information of acedrg")

        self.inputParser.add_option("-y", "--repcrd",
                  action="store_true", dest="repCrd", default=False,
                  help="Use this keyword if you want to replace the atomic coordinates in the input mmCif with those in the input PDB")

        self.inputParser.add_option("-z",  "--noGeoOpt", dest="noGeoOpt",  
                                    action="store_true",  default=False,
                                    help="The option for not doing geometry optimization on coordinates")

        self.inputParser.add_option("-L",  "--linkInstruction", dest="linkInstructions", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File that gives the instructions to build a link")

        self.inputParser.add_option("-T",  "--Test", dest="testMode",  
                                    action="store_true",  default=False,  
                                    help="The option for internal tests")

        (inputOptionsP, inputOptionsU) = self.inputParser.parse_args(t_argvs)

        if inputOptionsU:    
            print "Line arguments for those keywords are missing "
            for a_opt in inputOptionsU:
                print "Keyword %s ?"%a_opt
            sys.exit()

        return inputOptionsP 

    def checkDependency(self):
    
        if self.acedrgDir !="" and os.path.isdir(self.acedrgDir):
            tLibmol = os.path.join(self.acedrgDir, "libexec","libmol")
            if platform.system()=="Windows": tLibmol += ".exe"
            if glob.glob(tLibmol):
                self.libmol = tLibmol

        if not self.libmol  and os.environ.has_key("CCP4"):
            tLibmol = os.path.join(os.environ['CCP4'], "libexec","libmol")
            if platform.system()=="Windows": tLibmol += ".exe"
            if glob.glob(tLibmol):
                self.libmol = tLibmol

        if not self.libmol :
            if os.environ.has_key("LIBMOL_ROOT"):
                tLibmol = os.path.join(os.environ['LIBMOL_ROOT'], "libexec","libmol")
                if platform.system()=="Windows": tLibmol += ".exe"
                if glob.glob(tLibmol):
                    self.libmol = tLibmol
                else:
                    print "libmol could not be found at %s"%tLibmol
                    sys.exit()

        if not self.libmol: 
            print "can not find libmol at libexec/"
            sys.exit()
                                
        if os.environ.has_key("CCP4"):
            tRefmac = os.path.join(os.environ['CBIN'], "refmac5")
            if platform.system()=="Windows": tRefmac += ".exe"
            if not glob.glob(tRefmac):
                print "refmac5 could not be found"
                sys.exit()
            else:
                self.refmac = tRefmac

            tLibcheck = os.path.join(os.environ['CBIN'], "libcheck")
            if platform.system()=="Windows": tLibcheck += ".exe"
            if not glob.glob(tLibcheck):
                print "libcheck could not be found"
                sys.exit()
            else:
                self.libcheck = tLibcheck
            
            if not self.libmol:
                tLibmol = os.path.join(os.environ['CCP4'], "libexec", "libmol")
                if platform.system()=="Windows": tLibmol += ".exe"
                if glob.glob(tLibmol):
                    self.libmol = tLibmol
            if not self.acedrgTables:
                tAcedrgTables = os.path.join(os.environ['CCP4'], "share","acedrg","tables")
                if glob.glob(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
        else :
            print "You need to install CCP4 suite"
            print "or activate ccp4.setup"
            sys.exit()

        if not self.acedrgTables:
            if os.environ.has_key("LIBMOL_ROOT"):
                tAcedrgTables = os.path.join(os.environ['LIBMOL_ROOT'], "share","acedrg","tables")
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
                else:
                    print "Tables required could not be found at %s"%tAcedrgTables
                    sys.exit()
                                
        if not self.acedrgTables: 
            tAcedrgTables = os.path.join(self.acedrgDir, "share","acedrg","tables")
            # print tAcedrgTables
            if os.path.isdir(tAcedrgTables):
                self.acedrgTables = tAcedrgTables
            else:
                print "Bond and angle tables required could not be found "
                sys.exit()

        if self.acedrgTables:
            tFuncGroupTable = os.path.join(self.acedrgTables, "funSmi.table")
            if os.path.isfile(tFuncGroupTable):
                self.funcGroupTable = tFuncGroupTable
            
        #print "The path to Acedrg tables is at ", self.acedrgTables
        #print "Libmol used is at ", self.libmol
        
    def checkVersionInfo(self):
  
        # Acedrg version info 
        self.versionInfo["man"] = os.path.join(self.acedrgTables, "manifest.txt")
        if not os.path.isfile(self.versionInfo["man"]):
            print "Version infomation is not available."
        else:
            # print self.versionInfo["man"]
            try:
                vInfo = open(self.versionInfo["man"], "r")
            except IOError:
                print self.versionInfo["man"], " Could not be opened for reading"
                print "Version infomation is not available."
            else:
                for aL in vInfo.readlines():
                    if aL.find(":") !=-1:
                        strs = aL.strip().split(":")
                        if len(strs)==2:
                            self.versionInfo[strs[0]] = strs[1]
       
        self.versionInfo["RDKit_VERSION"] = rdBase.rdkitVersion 
        # Refmac version info 
        self._log_name    = os.path.join(self.scrDir, "refmac_version.log")
        self.runRefmacVersionInfo()
        if os.path.isfile(self._log_name):
            fRV = open(self._log_name, "r")
            allLs = fRV.readlines()
            fRV.close()
            for aL in allLs:
                strs = aL.strip().split()
                if len(strs)== 4 and aL.find("Program") !=-1:
                    self.versionInfo["REFMAC_NAME"] = strs[1].strip()[0:-1]
                    #print "REFMAC NAME : ", self.versionInfo["REFMAC_NAME"] 
                    self.versionInfo["REFMAC_VERSION"] = strs[3].strip()
                    #print "REFMAC VERSION : ", self.versionInfo["REFMAC_VERSION"] 
                    
        else: 
            print "Refmac version info is not available"

        if not self.versionInfo.has_key("REFMAC_VERSION"):
            print "Refmac version info is not available"
           
             
    def setOutCifGlobSec(self):
         
        self.outCifGlobSect.append("#loop_\n")
        self.outCifGlobSect.append("#_software\n")
        self.outCifGlobSect.append("#_version\n")
        self.outCifGlobSect.append("#_purpose\n")
        if self.versionInfo.has_key("ACEDRG_VERSION"):
            self.outCifGlobSect.append("#%s%s%s\n"%("acedrg".ljust(30), self.versionInfo["ACEDRG_VERSION"].ljust(20), '\"dictionary generator\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("#%s%s%s\n"%("acedrg".ljust(30), "?".ljust(20), '\"dictionary generator\"'.ljust(40)))
        
        if self.versionInfo.has_key("DATABASE_VERSION"):
            self.outCifGlobSect.append("#%s%s%s\n"%("acedrg_database".ljust(30), self.versionInfo["DATABASE_VERSION"].ljust(20), '\"data source\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("#%s%s%s\n"%("acedrg_database".ljust(30), "?".ljust(20), '\"data source\"'.ljust(40)))

        self.outCifGlobSect.append("#%s%s%s\n"%("rdkit".ljust(30), rdBase.rdkitVersion.ljust(20), '\"chemistry perception\"' )) 
  
        if self.versionInfo.has_key("REFMAC_NAME"):
            self.outCifGlobSect.append("#%s%s%s\n"%(self.versionInfo["REFMAC_NAME"].ljust(30), self.versionInfo["REFMAC_VERSION"].ljust(20),\
                                       '\"optimization tool\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("#%s%s%s\n"%("refmac".ljust(30), "?".ljust(20), '\"optimization tool\"'.ljust(40)))
        
         
    def checInputFormat(self):
        
        try:
            tF =open(self.geneInFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFName
            sys.exit()
        else:
            
            tLines = tF.readlines()
            tF.close()
            for aL in tLines:
                strs = aL.strip().split()
                if aL.find("_chem_comp_atom") != -1: 
                    self.inFileType = "MMCIF"
                    break     
                elif aL.find("_atom_site_") != -1:
                    self.geneInName = "CIF" 
                    break
                elif len(strs)==8:
                     if strs[7].find("V2000") != -1 or strs[7].find("V3000") != -1 :
                         self.geneInName = "MOL" 
                         break
                
    def setMonoRoot(self, tDataDesc=None):
   
        if len(self.monomRoot) !=0 and self.monomRoot.find("UNL")==-1:
            return
 
        if tDataDesc:
            for aIdx in tDataDesc.keys():
                if tDataDesc[aIdx][0].find("_chem_comp.id") !=-1:
                    self.monomRoot = tDataDesc[aIdx][1].strip()
                    return
        
        if self.inMdlName != "":
            try:
                tF =open(self.inMdlName, "r")
            except IOError:
                print "%s can not be open for reading "%self.inMdlName
                sys.exit()
            else:
                tL = tF.readline().strip()
                tF.close()
                #print tL
                if tL.find("#") !=-1:
                    self.monomRoot="UNL"
                else:
                    if len(tL) ==3:
                        self.monomRoot= tL
                    elif len(tL) >3:
                        self.monomRoot= tL[:4]
                
        elif self.inMmCifName != "":
            try:
                tF =open(self.inMmCifName, "r")
            except IOError:
                print "%s can not be open for reading "%self.inMmCifName
                sys.exit()
            else:
                tLs = tF.readlines()
                tF.close()
                for aL in tLs:
                    if aL.find("data_comp_") != -1:
                        strs = aL.strip().split("_") 
                        if len(strs) ==3:
                            if len(strs[2]) ==3:
                                self.monomRoot = strs[2]
                                break
        if len(self.monomRoot) < 3:
            self.monomRoot = "UNL"


    def setWorkMode(self, t_inputOptionsP = None):

        # Sequnence for check the locations of acedrg tables 
        # (1) Check if the user provides the location
        # (2) If not, check CCP4 suite default location.
        # (3) If not, check if the environment variable LIBMOL_ROOT is defined in the user's machine.
        # (4) If none of them, program exits.   
        if not t_inputOptionsP.acedrgTables: 
            if self.acedrgDir !="" and os.path.isdir(self.acedrgDir):
                tAcedrgTables = os.path.join(self.acedrgDir, "share","acedrg","tables")
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
            elif  not self.acedrgTables and os.environ.has_key("CCP4"):
                tAcedrgTables = os.path.join(os.environ['CCP4'], "share","acedrg","tables")
                if glob.glob(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
            else:
                print "You need to install CCP4 suite"
                sys.exit()
        else:
            if os.path.isdir(t_inputOptionsP.acedrgTables):
                self.acedrgTables = t_inputOptionsP.acedrgTables
 
        if not self.acedrgTables or not glob.glob(self.acedrgTables):
            if os.environ.has_key("LIBMOL_ROOT"):
                tAcedrgTables = os.path.join(os.environ['LIBMOL_ROOT'], "share","acedrg","tables")
                print tAcedrgTables
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables

        if not t_inputOptionsP.molGen and not t_inputOptionsP.repCrd and not t_inputOptionsP.typeOut\
           and not t_inputOptionsP.HMO and not t_inputOptionsP.linkInstructions: 
            if not t_inputOptionsP.noGeoOpt:
                if t_inputOptionsP.inMmCifName:
                    self.inMmCifName = t_inputOptionsP.inMmCifName
                    self.workMode    = 11            
                elif t_inputOptionsP.inSmiName: 
                    self.inSmiName = t_inputOptionsP.inSmiName
                    self.workMode    = 12            
                elif t_inputOptionsP.inMdlName: 
                    self.inMdlName = t_inputOptionsP.inMdlName
                    self.workMode    = 13
                elif t_inputOptionsP.inSdfName: 
                    self.inSdfName = t_inputOptionsP.inSdfName
                    self.workMode    = 14
                elif t_inputOptionsP.inMol2Name: 
                    self.inMol2Name = t_inputOptionsP.inMol2Name
                    self.workMode    = 15
                elif t_inputOptionsP.inLigandPdbName: 
                    self.inLigandPdbName = t_inputOptionsP.inLigandPdbName
                    self.workMode    = 16
                    self.useExistCoords = True
            else:
                if t_inputOptionsP.inMmCifName:
                    self.inMmCifName = t_inputOptionsP.inMmCifName
                    self.workMode    = 111            
                elif t_inputOptionsP.inSmiName: 
                    self.inSmiName = t_inputOptionsP.inSmiName
                    self.workMode    = 121            
                elif t_inputOptionsP.inMdlName: 
                    self.inMdlName = t_inputOptionsP.inMdlName
                    self.workMode    = 131
                elif t_inputOptionsP.inSdfName: 
                    self.inSdfName = t_inputOptionsP.inSdfName
                    self.workMode    = 141
                elif t_inputOptionsP.inMol2Name: 
                    self.inMol2Name = t_inputOptionsP.inMol2Name
                    self.workMode    = 151
                elif t_inputOptionsP.inLigandPdbName: 
                    self.inLigandPdbName = t_inputOptionsP.inPdbName
                    self.workMode    = 161
                    self.useExistCoords = True
                  

        elif t_inputOptionsP.molGen:
            if t_inputOptionsP.inStdCifName: 
                self.inStdCifName = t_inputOptionsP.inStdCifName
                self.workMode     = 21
            elif t_inputOptionsP.inStdCifDir: 
                self.inStdCifDir = t_inputOptionsP.inStdCifDir
                self.workMode    = 22
        elif t_inputOptionsP.typeOut:
            if t_inputOptionsP.inMmCifName:
                self.inMmCifName = t_inputOptionsP.inMmCifName
                self.workMode    = 31
            elif t_inputOptionsP.inSmiName:
                self.inSmiName = t_inputOptionsP.inSmiName
                self.workMode    = 32
            elif t_inputOptionsP.inMdlName:
                self.inMdlName = t_inputOptionsP.inMdlName
                self.workMode    = 33
            elif t_inputOptionsP.inSdfName:
                self.inSdfName = t_inputOptionsP.inSdfName
                self.workMode    = 34
            if t_inputOptionsP.inStdCifName:
                self.inStdCifName = t_inputOptionsP.inStdCifName
                self.workMode     = 35
        elif t_inputOptionsP.repCrd :
            if t_inputOptionsP.inPdbName and t_inputOptionsP.inMmCifName: 
                self.inMmCifName = t_inputOptionsP.inMmCifName
                self.inPdbName = t_inputOptionsP.inPdbName
                self.workMode     = 41
        elif t_inputOptionsP.HMO :
            if t_inputOptionsP.inMmCifName:
                self.inMmCifName = t_inputOptionsP.inMmCifName
                self.workMode    = 51            
            elif t_inputOptionsP.inSmiName: 
                self.inSmiName = t_inputOptionsP.inSmiName
                self.workMode    = 52            
            elif t_inputOptionsP.inMdlName: 
                self.inMdlName = t_inputOptionsP.inMdlName
                self.workMode    = 53
            elif t_inputOptionsP.inSdfName: 
                self.inSdfName = t_inputOptionsP.inSdfName
                self.workMode    = 54
            elif t_inputOptionsP.inMol2Name: 
                self.inMol2Name = t_inputOptionsP.inMol2Name
                self.workMode    = 55
        elif t_inputOptionsP.linkInstructions :
            # Need to add more options here
            self.linkInstructions = t_inputOptionsP.linkInstructions.strip()
            self.workMode = 61
 

        if t_inputOptionsP.testMode :
            self.testMode = True

        if t_inputOptionsP.monomRoot:
            self.monomRoot   = t_inputOptionsP.monomRoot
        else:
            #print self.monomRoot
            self.setMonoRoot()
            # self.monomRoot   = "UNL"
        #print "monomRoot: ", self.monomRoot
        if t_inputOptionsP.outRoot:
            self.outRoot   = t_inputOptionsP.outRoot
            self.baseRoot  = os.path.basename(t_inputOptionsP.outRoot)
            tDir = os.path.dirname(t_inputOptionsP.outRoot)
            if tDir:
                if platform.system()=="Windows":
                    tStrGrp = os.path.abspath(tDir).strip().split("\\")
                    tCurDir = "\\"
                else:
                    tStrGrp = os.path.abspath(tDir).strip().split("/")
                    tCurDir = "/"
                for tSub in tStrGrp:
                    tCurDir = os.path.join(tCurDir, tSub)
                    if not os.path.isdir(tCurDir):
                        os.mkdir(tCurDir)
        else:
            self.outRoot   = "AcedrgOut"
            self.baseRoot   = "AcedrgOut"
        

        if t_inputOptionsP.numConformers:
            self.numConformers = t_inputOptionsP.numConformers
            
      
        if t_inputOptionsP.useExistCoords:
            self.useExistCoords = t_inputOptionsP.useExistCoords

        self.scrDir = self.outRoot + "_TMP"
        if not os.path.isdir(self.scrDir):
            os.mkdir(self.scrDir)

        if self.inSmiName !="":
            if not os.path.isfile(self.inSmiName): 
                tName = os.path.join(self.scrDir, "inputSmiles.smi")
                try :
                    tF =open(tName, "w")
                except IOError:
                    print "Could not put the input smiles in a file "
                    sys.exit(1)
                else:
                    tF.write(self.inSmiName+ "\n")
                    tF.close()
                    self.inSmiName = tName
                    print tName

    def setInputProcPara(self, t_inputOptionsP = None):
        
        if t_inputOptionsP.useExistCoords:
            self.inputPara["useExistCoords"]    = t_inputOptionsP.useExistCoords

        if t_inputOptionsP.numInitConformers:
            self.inputPara["numInitConformers"] = t_inputOptionsP.numInitConformers

        if t_inputOptionsP.numRDKitOptmSteps:
            self.inputPara["numRDKitOptmSteps"] = t_inputOptionsP.numRDKitOptmSteps

        if t_inputOptionsP.numConformers:
            self.inputPara["numConformers"]     = t_inputOptionsP.numConformers
  
        if t_inputOptionsP.numInitConformers and t_inputOptionsP.numConformers:
            if t_inputOptionsP.numConformers > t_inputOptionsP.numInitConformers:
                self.inputPara["numInitConformers"] = t_inputOptionsP.numConformers


    def printJobs(self):

        print "=====================================================================" 
        if self.versionInfo.has_key("ACEDRG_VERSION"):
            print "| ACEDRG version:  %s|"%self.versionInfo["ACEDRG_VERSION"].ljust(49)
        else: 
            print "=====================================================================" 
            print "| ACEDRG version is not available                                   |"    
        if self.versionInfo.has_key("DATABASE_VERSION"):
            print "| ACEDRG database: %s|"%self.versionInfo["DATABASE_VERSION"].ljust(49)
        else:
            print "| ACEDRG Database version is not available                          |"
        print "| RDKit version:  %s|"%rdBase.rdkitVersion.ljust(50)
        if self.versionInfo.has_key("REFMAC_NAME") and self.versionInfo.has_key("REFMAC_VERSION"):
            print "| %s  %s|"%((self.versionInfo["REFMAC_NAME"] + ":").ljust(15), self.versionInfo["REFMAC_VERSION"].ljust(49))
        print "=====================================================================" 
        if self.workMode in [11, 12, 13, 14, 15, 16, 111, 121, 131, 141, 151, 161] :
        #if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 or self.workMode==15 \
        #   or self.workMode == 111 or self.workMode==121 or self.workMode ==131 or self.workMode==141 or self.workMode==151 :
            print "=====================================================================" 
            print "| Your job is  generating the dictionary (cif) and coord(pdb) files |"
            print "| for your ligand and/or monomer                                    |"
            print "=====================================================================" 
            if self.workMode==11 or self.workMode==111:
                print "Input file: %s"%os.path.basename(self.inMmCifName)
            if self.workMode==12 or self.workMode==121:
                print "Input file: %s"%os.path.basename(self.inSmiName)
            if self.workMode==13 or self.workMode==131:
                print "Input file: %s"%os.path.basename(self.inMdlName)
            if self.workMode==14 or self.workMode==141:
                print "Input file: %s"%os.path.basename(self.inSdfName)
            if self.workMode==15 or self.workMode==151:
                print "Input file: %s"%os.path.basename(self.inMol2Name)
            if self.workMode==16 or self.workMode==161:
                print "Input file: %s"%os.path.basename(self.inLigandPdbName)
            print "Output dictionary file: %s"%self.outRoot + ".cif"
            if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 or self.workMode==15:
                print "Output coordinate file: %s"%self.outRoot + ".pdb"

        if self.workMode == 21 or self.workMode==22:
            print "=====================================================================" 
            print "| Your job is to generate molecules (sets of connected atoms), to   |"
            print "| get unique bonds and angles within the molecules and cluster them |" 
            print "| by specificted desigend atom types.                               |" 
            print "=====================================================================" 
            if self.workMode == 21:
                print "Input file: %s"%os.path.basename(self.inStdCifName)
                print "Output molecules: %s"%self.outRoot + "_all_mols.txt"
                print "Output bonds and angles : %s"%self.outRoot + "_unique_bond_and_angles.txt"
            elif self.workMode==22:
                print "Input directory where cif files are: %s"%self.inStdCifDir
                print "Output bond and angle file : %s "%(self.outRoot + "_all_bonds_and_angles.table")

        if self.workMode in [51, 52, 53, 54, 55] :
            print "=====================================================================" 
            print "| Your job is to calculate bond-orders and charges using HMO.      |"
            print "=====================================================================" 

        if self.workMode == 61:
            print "=====================================================================" 
            print "| Your job is to generate full descriptions for two monomers, their |"
            print "| modifications, and a link between them.                           |" 
            print "=====================================================================" 

    def runLibmol(self, tIn=None, tIdxMol=-1):
        print "workMode ", self.workMode
        self._cmdline = self.libmol
        if tIdxMol !=-1: 
            self._log_name       = os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(tIdxMol) + "_cod.log")
        else:
            self._log_name       = os.path.join(self.scrDir, self.baseRoot +  "_cod.log")

        if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 \
           or self.workMode==15 or self.workMode ==16\
           or self.workMode == 111 or self.workMode == 121 or self.workMode == 131 \
           or self.workMode==141 or self.workMode==151 or self.workMode==161:
            print "===================================================================" 
            print "| Generate the dictionary file using the internal database        |"
            print "===================================================================" 
            if tIdxMol !=-1: 
                self.outRstCifName   = os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(tIdxMol) + "_cod.rst")   
                self.outRstPdbName   = os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(tIdxMol) + "_cod.pdb")
            else:
                self.outRstCifName   = os.path.join(self.scrDir, self.baseRoot + "_cod.rst")   
                self.outRstPdbName   = os.path.join(self.scrDir, self.baseRoot + "_cod.pdb")
        if self.workMode == 51:
            self.outRstCifName   =  self.baseRoot + "_bondOrder.list"
        if self.workMode == 11 or self.workMode == 12 or self.workMode==16\
           or self.workMode == 111 or self.workMode == 121 or self.workMode==161:
            if tIn:
                self.inMmCifName    = tIn
            self._cmdline +=" -c %s -D %s "%(self.inMmCifName, self.acedrgTables)
            self._cmdline += " -r %s -o %s "%(self.monomRoot, self.outRstCifName)
            #print self._cmdline
            #os.system(self._cmdline)
            self.runExitCode = self.subExecute()
        if self.workMode == 13 or self.workMode == 14 or self.workMode == 131 or self.workMode == 141:
            if tIn:
                self.inMmCifName = tIn
            self._cmdline += " -s %s  -D %s "%(self.inMmCifName, self.acedrgTables)
            self._cmdline += " -r %s -o %s "%(self.monomRoot, self.outRstCifName)
            #print self._cmdline
            #os.system(self._cmdline)
            self.runExitCode = self.subExecute()
        if self.workMode == 15 or self.workMode == 151 :
            if tIn:
                self.inMmCifName = tIn
            self._cmdline += " -k %s  -D %s "%(self.inMmCifName, self.acedrgTables)
            self._cmdline += " -r %s -o %s "%(self.monomRoot, self.outRstCifName)
            # print self._cmdline  
            self.runExitCode = self.subExecute()
           

        if self.workMode == 21 :
            if tIn:
                self.inStdCifName = tIn
            self.outRstCifName  = self.outRoot + ".cif"
	    self.outMolsName    = self.monomRoot + "_all_mols.txt"
            self.outBondsAndAnglesName  = self.monomRoot + "_unique_bond_and_angles.txt"

            self._cmdline += " -b %s  "%self.inStdCifName
            self._cmdline += " -m yes -r %s -o %s "%(self.monomRoot, self.outRoot)
            #print self._cmdline
            #print "self.outRoot ", self.outRoot 
            self.subExecute()

        if self.workMode == 22 :
            if os.path.isdir(self.inStdCifDir):
                bTable = self.outRoot + "_all_atoms_bonds_angles.table"
                self.workMode  == 21
                tAllBondInMols = []
                tExcludedCif   = []
                tempCifStr = os.path.join(self.inStdCifDir,"*.cif")
                for aCif in glob.glob(tempCifStr):
                    tMonomRoot = os.path.basename(aCif).split(".")[0]
                    print "Generate molecules, bonds and angles from ", aCif
                    self._cmdline = self.libmol
                    self.outMolsName            = os.path.join(self.scrDir,tMonomRoot + "_all_mols.txt")
                    self.outBondsAndAnglesName  = os.path.join(self.scrDir, tMonomRoot + "_unique_bond_and_angles.txt")
                    self._log_name = os.path.join(self.scrDir,tMonomRoot + "_cod.log")
                    self._cmdline += " -b %s  "%aCif
                    tempMonRt = os.path.join(self.scrDir, tMonomRoot)
                    self._cmdline += " -m yes -r %s -o %s.cif "%(tMonomRoot, tempMonRt)
                    #print self._cmdline
                    self.subExecute()
                    if os.path.isfile(self.outBondsAndAnglesName):
                        tAllBondInMols.append(self.outBondsAndAnglesName)
                    else:
                        print "No bonds and angles are generated from ", aCif
                        tExcludedCif.append(aCif)
                if len(tAllBondInMols) !=0:
                    self.getBondsAndAngles(bTable, tAllBondInMols)
                if len(tExcludedCif) !=0:
                    print "The following cif files have been excluded from calculations because their large R factors"
                    for aExCif in tExcludedCif:
                        print aExCif

        if self.workMode == 31 or self.workMode == 32 or self.workMode == 33 or self.workMode == 34 or self.workMode == 35 : 
            #print "===================================================================" 
            #print "| Generate atom types of Acedrg style                             |"
            #print "===================================================================" 
            if self.workMode == 31 or self.workMode == 32 or self.workMode==33:
                if tIn:
                    inFileName    = tIn
                else:
                    inFileName    = self.inMmCifName
                if self.workMode == 31 or self.workMode == 32:
                    if tIdxMol !=-1:  
                        aStr = "1"
                    else:
                        aStr = "2"
                    self.outAtmTypeName = os.path.join(self.scrDir, "atomTypes_"+aStr+ ".txt")
                elif self.workMode == 33:
                    print "===================================================================" 
                    print "| Generate atom types of Acedrg style                             |"
                    print "===================================================================" 
                    self.outAtmTypeName = self.outRoot + "_atomTypes.txt"
                self._cmdline +=" -A yes -D %s -c %s  -o %s "%(self.acedrgTables, inFileName, self.outAtmTypeName)
                #print self._cmdline
                self.subExecute()
                
        if self.workMode == 41 :

            self._cmdline += " -p %s -c %s "%(self.inPdbName, self.inMmCifName)
            self._cmdline += " -y y -o %s "%(self.outRstCifName)
            #print self._cmdline
            self.subExecute()

        if self.workMode == 51 :

            self._cmdline += " -Z yes -c %s -o %s "%(self.inMmCifName, self.outRstCifName)
            #print self._cmdline
            self.subExecute()
        if self.workMode ==900:
            self._cmdline += " -T yes -X %s -Y %s"%(self.libmolAT1, self.libmolAT2)
            self._cmdline += " -o %s "%self.libmolMatched
            self.subExecute()
        
    def getBondsAndAngles(self, tFName, tMolTabs):
        
        for aMonTab in tMolTabs:
            try:
                f1 = open(aMonTab, "r")
            except IOError:
                print "%s has not been found for reading"%aMonTab
            else:
                f1_lines = f1.readlines()
                f1.close()
                aSetStrs = os.path.basename(aMonTab).strip().split("_")
                aFileIdx = aSetStrs[0]
                lBo = False
                lAn = False
                for aL in f1_lines:
                    strGrp = aL.strip().split()
                    if lBo and len(strGrp)==6:
                         # get sorted bonds
                         aElem = self.getElemFromAtomClass(strGrp[0])
                         bElem = self.getElemFromAtomClass(strGrp[1])
                         aSet  = [strGrp[0], aElem, strGrp[2]]                  
                         bSet  = [strGrp[1], bElem, strGrp[3]]   
                         tSets = [aSet, bSet]
                         tSets.sort(listComp)
                         if not self.allBondsAndAngles["bonds"].has_key(tSets[0][1]):         # class 1
                             self.allBondsAndAngles["bonds"][tSets[0][1]] = {}
                         if not self.allBondsAndAngles["bonds"][tSets[0][1]].has_key(tSets[1][1]): # class 2 
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]] = {}
                         if not self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]].has_key(tSets[0][0]): # id 1
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]] ={}
                         if not self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]].has_key(tSets[1][0]): # id 2
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]] = {}
                         if not self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]].has_key(strGrp[5]):
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]] = {}
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["observations"] = []
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["stats"]        = {}

                         aTup = [tSets[0][2], tSets[1][2], float(strGrp[4]), aFileIdx]              # atom1_id, atom2_id, bond length and the file index
                         self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["observations"].append(aTup) 
                

                         # Get sorted atom types
                         if not self.allBondsAndAngles["atomClasses"].has_key(tSets[0][1]):            # atom 1 element type
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]] = {}
                         if not self.allBondsAndAngles["atomClasses"][tSets[0][1]].has_key(tSets[0][0]): # atom 1 class type
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]][tSets[0][0]] = 1
                         else:
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]][tSets[0][0]] +=1

                         if not self.allBondsAndAngles["atomClasses"].has_key(tSets[1][1]):            # atom 2 element type
                              self.allBondsAndAngles["atomClasses"][tSets[1][1]] = {}
                         if not self.allBondsAndAngles["atomClasses"][tSets[1][1]].has_key(tSets[1][0]): # atom 2 class type
                              self.allBondsAndAngles["atomClasses"][tSets[1][1]][tSets[1][0]] = 1
                         else:
                              self.allBondsAndAngles["atomClasses"][tSets[1][1]][tSets[1][0]] +=1
                    elif lAn and len(strGrp)==7:
                         # get sorted bond-angles
                         cenElem = self.getElemFromAtomClass(strGrp[0])
                         aElem   = self.getElemFromAtomClass(strGrp[1])
                         bElem   = self.getElemFromAtomClass(strGrp[2])
                         aSet  = [strGrp[0], cenElem, strGrp[3]]                  
                         bSet  = [strGrp[1], aElem, strGrp[4]]   
                         cSet  = [strGrp[2], bElem, strGrp[5]]   
                         tSets = [bSet, cSet]
                         tSets.sort(listComp)
                         if not self.allBondsAndAngles["angles"].has_key(aSet[1]):         # center atom element  
                             self.allBondsAndAngles["angles"][aSet[1]] = {}
                         if not self.allBondsAndAngles["angles"][aSet[1]].has_key(tSets[0][1]): # atom1 element 
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]] = {}
                         if not self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]].has_key(tSets[1][1]) : # atom2 elment 
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]] = {}
                         if not self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]].has_key(aSet[0]): # center atom class
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]] = {}
                         if not self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]].has_key(tSets[0][0]): #  atom1 class
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]] = {}
                         if not self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]].has_key(tSets[1][0]): #  atom1 class
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]][tSets[1][0]] = {}
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]][tSets[1][0]]["observations"] = []
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]][tSets[1][0]]["stats"]        = {}
                         aTup = [aSet[2], tSets[0][2], tSets[1][2], float(strGrp[6]), aFileIdx]              # cen_atom_id, atom1_id, atom2_id, bond length and the file index
                         self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]][tSets[1][0]]["observations"].append(aTup)
                    elif aL.find("Bond_length") != -1:
                        lBo = True
                    elif aL.find("Angle_") !=-1:
                        lAn = True
                        lBo = False

        # output atom types, bonds and angles
        try:
            f2 = open(tFName, "w")
        except IOError:
            print "%s has not been found for writing"%tFName
        else:
            if len(self.allBondsAndAngles["atomClasses"]):
                f2.write("All unique atom types: \n")
                for aElem in sorted(self.allBondsAndAngles["atomClasses"].iterkeys()):
                    for aCl in sorted(self.allBondsAndAngles["atomClasses"][aElem].iterkeys()):
                        f2.write("%s        %s     \n"%(aElem,aCl))   
            if len(self.allBondsAndAngles["bonds"]):
                f2.write("All bond lengths: \n")  
                for aElem in sorted(self.allBondsAndAngles["bonds"].iterkeys()):
                    for bElem in sorted(self.allBondsAndAngles["bonds"][aElem].iterkeys()):
                        for aCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem].iterkeys()):
                            for bCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl].iterkeys()):
                                for rP in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl].iterkeys()):
                                    self.getStatsForOneBondClass(self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["observations"], \
                                                            self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["stats"])                                  
                                    for aSet in self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["observations"]:
                                        if self.testMode:
                                            f2.write("%s       %s        %s      %s      %s      %s    %s     %7.5f    %s\n" \
                                                     %(aElem, bElem, aCl, bCl, aSet[0],  aSet[1], rP, aSet[2], aSet[3]))   
                                        else:
                                            f2.write("%s       %s        %s      %s      %s      %s    %s     %7.5f\n" \
                                                     %(aElem, bElem, aCl, bCl, aSet[0],  aSet[1], rP, aSet[2]))   

                f2.write("\nAll bond lengths stats: \n")
                if self.testMode: 
                    bCaseIdx =0
                    casesDir = os.path.join(self.scrDir, "CasesNSTATS")
                    if not glob.glob(casesDir):
                        os.mkdir(casesDir)
                     
                
                for aElem in sorted(self.allBondsAndAngles["bonds"].iterkeys()):
                    for bElem in sorted(self.allBondsAndAngles["bonds"][aElem].iterkeys()):
                        for aCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem].iterkeys()):
                            for bCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl].iterkeys()):
                                for rP in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl].iterkeys()):
                                    ave = self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["stats"]["mean"]
                                    sig = self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["stats"]["sig"]
                                    num = self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["stats"]["nObs"]
                                    f2.write("%s      %s    %7.5f      %7.5f      %d\n" %(aCl, bCl, ave, sig, num)) 
                                    if self.testMode and num > 20 and aElem !="H" and bElem != "H" :
                                        bCaseIdx +=1
                                        tCasesNameR = aElem + "_" + bElem + str(bCaseIdx) + ".txt"
                                        tCasesName  = os.path.join(casesDir, tCasesNameR)
                                        tCases      = open(tCasesName, "w")
                                        for aSet in self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl][rP]["observations"]:
                                            tCases.write("%s      %s      %s      %s    %s     %7.5f    %s\n" \
                                                     %(aCl, bCl, aSet[0],  aSet[1], rP, aSet[2], aSet[3]))        
                                        tCases.close()
                                        
  
            if len(self.allBondsAndAngles["angles"]):
                f2.write("\n\nAll bond angles: \n")
                for cenElem in sorted(self.allBondsAndAngles["angles"].iterkeys()):
                    for aElem in sorted(self.allBondsAndAngles["angles"][cenElem].iterkeys()):
                        for bElem in sorted(self.allBondsAndAngles["angles"][cenElem][aElem].iterkeys()):
                            for cenC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem].iterkeys()):
                                for aC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC].iterkeys()):
                                    for bC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC].iterkeys()):
                                        self.getStatsForOneAngleClass(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["observations"], \
                                                                      self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["stats"])
                                        for aSet in self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["observations"]:
                                            if self.testMode:
                                                f2.write("%s      %s       %s        %s      %s      %s      %s        %s      %s     %7.5f    %s\n" \
                                                          %(cenElem, aElem, bElem, cenC, aC, bC, aSet[0],  aSet[1], aSet[2], aSet[3], aSet[4]))   
                                            else:
                                                f2.write("%s      %s       %s        %s      %s      %s      %s        %s      %s     %7.5f\n" \
                                                          %(cenElem, aElem, bElem, cenC, aC, bC, aSet[0],  aSet[1], aSet[2], aSet[3]))   
                f2.write("\nAll bond angle stats : \n")
                for cenElem in sorted(self.allBondsAndAngles["angles"].iterkeys()):
                    for aElem in sorted(self.allBondsAndAngles["angles"][cenElem].iterkeys()):
                        for bElem in sorted(self.allBondsAndAngles["angles"][cenElem][aElem].iterkeys()):
                            for cenC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem].iterkeys()):
                                for aC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC].iterkeys()):
                                    for bC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC].iterkeys()):
                                        ave  = self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["stats"]["mean"]
                                        sig  = self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["stats"]["sig"]
                                        nObs = self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC][bC]["stats"]["nObs"]
                                        f2.write("%s      %s      %s   %7.5f    %7.5f   %d\n" \
                                                  %(cenC, aC, bC, ave, sig, nObs))   
            f2.close()
   
    def getElemFromAtomClass(self, tAtm):
        
        tElem = ""
        if len(tAtm) != 0:
            tStrs = tAtm.strip().split("(")
            if len(tStrs) != 0:
                if tStrs[0].find("[")==-1:
                    tElem = tStrs[0]    
                else :
                    tElem = tStrs[0].split("[")[0]
        return tElem         


    def getStatsForOneBondClass(self, tBondObjList, tBondStatsDict):

        
        numObjs = len(tBondObjList)

        if numObjs:
     
            tBondStatsDict["nObs"] = numObjs
            tBondStatsDict["max"]  = 0.0
            tBondStatsDict["min"]  = 100.0

            sum = 0.0
            for aTup in tBondObjList:
                if len(aTup)==4:
                    sum = sum + aTup[2]
                    if aTup[2] > tBondStatsDict["max"]:
                        tBondStatsDict["max"] = aTup[2]
                    if aTup[2] < tBondStatsDict["min"]:
                        tBondStatsDict["min"] = aTup[2]
                else: 
                    print "Bug. error in ", aTup
                    sys.exit(1)

            tBondStatsDict["mean"] = sum/numObjs

            sum_diff_sq = 0.0

            for aTup in tBondObjList:
                sum_diff_sq +=((aTup[2]- tBondStatsDict["mean"])*(aTup[2]- tBondStatsDict["mean"]))
                              
            if numObjs > 1:
                tBondStatsDict["sig"] = math.sqrt(sum_diff_sq/(numObjs-1))
            else:
                tBondStatsDict["sig"] = 0.0  
            
            
    def getStatsForOneAngleClass(self, tAngleObjList, tAngleStatsDict):
        
        numObjs = len(tAngleObjList)

        if numObjs:
     
            tAngleStatsDict["nObs"] = numObjs
            tAngleStatsDict["max"]  = 0.0
            tAngleStatsDict["min"]  = 360.0

            sum = 0.0
            for aTup in tAngleObjList:
                if len(aTup)==5:
                    if (aTup[3] > 180.0):
                        aTup[3] = 360.0-aTup[3]
                    elif (aTup[3] < 0.0):
                        aTup[3] = math.fabs(aTup[3])

                    sum = sum + aTup[3]
                    if aTup[3] > tAngleStatsDict["max"]:
                        tAngleStatsDict["max"] = aTup[3]
                    if aTup[3] < tAngleStatsDict["min"]:
                        tAngleStatsDict["min"] = aTup[3]
                else: 
                    print "Bug. error in ", aTup
                    sys.exit(1)

            tAngleStatsDict["mean"] = sum/numObjs

            sum_diff_sq = 0.0

            for aTup in tAngleObjList:
                sum_diff_sq +=((aTup[3]- tAngleStatsDict["mean"])*(aTup[3]- tAngleStatsDict["mean"]))

            if numObjs > 1:
                tAngleStatsDict["sig"] = math.sqrt(sum_diff_sq/(numObjs-1))
            else:
                tAngleStatsDict["sig"] = 0.0  
            
            
    def setLibcheckBat(self, tSmiName, tOutRoot):

        self._cmdline = self.libcheck 

        libcheckBatName = os.path.join(self.scrDir, self.baseRoot + "_libcheck.bat")
        try:
            libcheckBat = open(libcheckBatName,"w")
        except IOError:
            print libcheckBatName, " could not be opened for write "
            sys.exit()
        else:
            libcheckBat.write(" N \n")
            libcheckBat.write("FILE_SMILE %s\n"%tSmiName) 
            libcheckBat.write("FILE_O  %s \n\n"%tOutRoot)
            libcheckBat.close()

            self._cmdline +=" < %s "%(libcheckBatName)
            #print self.cmdline
    
    def runLibcheck(self, tSmiName=None):

        if tSmiName:
            self.inSmiName = tSmiName
        self._log_name        = os.path.join(self.scrDir, self.baseRoot + "_libcheck.log")
        self.libcheckOutRoot  = os.path.join(self.scrDir, self.baseRoot + "_libcheck")
        self.libcheckLibName  = self.libcheckOutRoot +".lib"
        self.setLibcheckBat(self.inSmiName,  self.libcheckOutRoot) 
        self.subExecute()

    def setRefmacCom1(self, tPdbIn, tLibIn, tPdbOut, tStage=2):
     
        self._cmdline = self.refmac 
        self._cmdline += "  xyzin %s libin %s xyzout %s"%(tPdbIn, tLibIn, tPdbOut)
        if platform.system()!="Windows": self._cmdline += " <<eof\n"
        if tStage==0:
            tLog = os.path.join(self.scrDir, "init_cif.log")
            if platform.system()=="Windows":
                self._cmdline = "(ECHO ligand && ECHO ncyc 40 ECHO make hout yes ECHO make hydr full ECHO end) | " + self._cmdline
            else:
                self._cmdline += "ligand    \n"
                self._cmdline += "ncyc 40   \n"
                self._cmdline += "make hout yes    \n"
                self._cmdline += "make hydr full \n"
        elif tStage==1:
            if platform.system()=="Windows":
                self._cmdline = "(ECHO make hout yes && ECHO make hydr full && ECHO ncyc 40 && ECHO vdwr 2.0 && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "make hout yes    \n"
                self._cmdline += "make hydr full \n"
                self._cmdline += "ncyc 40    \n"
                self._cmdline += "vdwr 2.0    \n"
        elif tStage==2:
            if platform.system()=="Windows":
                self._cmdline += "(ECHO mode newe && ECHO make hydr no && ECHO make news && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "mode newe \n"
                self._cmdline += "make hydr no"
                self._cmdline += "make news"
        elif tStage==3:
            if platform.system()=="Windows":
                self._cmdline += "(ECHO make hout yes && ECHO make hydr full && ECHO ncyc 40 && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "make hout yes    \n"
                self._cmdline += "make hydr full \n"
                self._cmdline += "ncyc 40    \n"
        else: 
            if platform.system()=="Windows":
                self._cmdline += "(ECHO mode newe && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "mode newe \n"

        if platform.system()!="Windows":
            self._cmdline += "end       \n"
            self._cmdline += "eof       \n"

        #print self._cmdline

    def setRefmacCom2(self, tPdbIn, tLibIn, tPdbOut, tStage=2):
     
        self._cmdline = self.refmac
        self._cmdline += "  xyzin %s libin %s xyzout %s"%(tPdbIn, tLibIn, tPdbOut)
        if platform.system()!="Windows": self._cmdline += " <<eof\n"
        if tStage==1:
            if platform.system()=="Windows":
                self._cmdline += "(ECHO mode newe && ECHO make hout yes && ECHO make hydr full && ECHO ncyc 40 && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "mode newe \n"
                self._cmdline += "make hout yes    \n"
                self._cmdline += "make hydr full \n"
                self._cmdline += "ncyc 40    \n"
        elif tStage==2: 
            if platform.system()=="Windows":
                self._cmdline += "(ECHO mode newe && ECHO make hydr no && ECHO make news && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "mode newe \n"
                self._cmdline += "make hydr no"
                self._cmdline += "make news"
        elif tStage==3: 
            if platform.system()=="Windows":
                self._cmdline += "(ECHO make hout yes && ECHO make hydr full && ECHO ncyc 40 && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "make hout yes    \n"
                self._cmdline += "make hydr full \n"
                self._cmdline += "ncyc 40    \n"
        else:
            if platform.system()=="Windows":
                self._cmdline += "(ECHO mode newe && ECHO end) | " + self._cmdline
            else:
                self._cmdline += "mode newe \n"

        if platform.system()!="Windows":
            self._cmdline += "end       \n"
            self._cmdline += "eof       \n"
        #print self._cmdline

    def runRefmac(self, tPdbIn, tLibIn, tStage=2):
      
        print "| Stage  %d        |"%tStage
        self._log_name    = os.path.join(self.scrDir, self.baseRoot + "_refmac_stage_"+ str(tStage) + ".log")
        self.refmacXYZOUTName = os.path.join(self.scrDir, self.baseRoot + "_refmac_stage_"+ str(tStage) + ".pdb")
        self.setRefmacCom1(tPdbIn, tLibIn, self.refmacXYZOUTName, tStage)
        self.subExecute()
    
    def runRefmac(self, tPdbIn, tLibIn, tRoot, tStage=2):
      
        #print "| Stage  %d        |"%tStage
        #self._log_name    = os.path.join(self.scrDir, tRoot + "_refmac_stage_"+ str(tStage) + ".log")
        #if tStage ==0:
        #    self._log_name    = os.path.join(self.scrDir, "init_cif.log")
        #    self.refmacXYZOUTName = os.path.join(self.scrDir, tRoot+ "_init" + ".cif")
        #else:
        self._log_name    = os.path.join(self.scrDir, tRoot + "_refmac_stage_"+ str(tStage) + ".log")
        self.refmacXYZOUTName = os.path.join(self.scrDir, tRoot+ "_refmac_stage_"+ str(tStage) + ".pdb")
        self.setRefmacCom1(tPdbIn, tLibIn, self.refmacXYZOUTName, tStage)
        self.subExecute()

    def runRefmacVersionInfo(self):

        if self._log_name == "":
            self._log_name    = os.path.join(self.scrDir, "refmac_version.log")

        self._cmdline = self.refmac  + " -i "
        self.subExecute()

    def runGeoOpt(self):
        
        # Geometrical optimization
        if os.path.isfile(self.outRstCifName) and os.path.isfile(self.outRstPdbName):
            if os.path.getsize(self.outRstPdbName) > 100 and os.path.getsize(self.outRstCifName) > 50:
                stageNow = 1
                self.runRefmac(self.outRstPdbName, self.outRstCifName, stageNow)
                if os.path.isfile(self.refmacXYZOUTName):
                    tPdb1 = self.refmacXYZOUTName
                    stageNow = 2
                    self.runRefmac(tPdb1, self.outRstCifName, stageNow)
                    if os.path.isfile(self.refmacXYZOUTName):
                        tPdb1 = self.refmacXYZOUTName
                        stageNow = 3
                        self.runRefmac(tPdb1, self.outRstCifName, stageNow)
                        if os.path.isfile(self.refmacXYZOUTName):
                            finPdb = self.outRoot + ".pdb"
                            finRst = self.outRoot + ".cif"
                            #os.system("cp %s   %s "%(self.refmacXYZOUTName, finPdb))
                            shutil.copy(self.refmacXYZOUTName, finPdb)
                            #os.system("cp %s   %s "%(self.outRstCifName, finRst))
                            if os.path.isfile(finPdb):
                                self.inPdbName        = finPdb 
                                self.inMmCifName      = self.outRstCifName
                                self.outRstCifName    = finRst
                                self.transCoordsPdbToCif(self.inPdbName, self.inMmCifName, self.outRstCifName)
                                print "===================================================================" 
                            else:
                                print "Failed to produce %s after final geometrical optimization"%finPdb
                        else:
                            print "Failed to produce the  coordinates at stage 3 optimization" 
                    else:
                        print "Failed to produce the coordinates at stage 2 optimization" 
                else:
                    print "Failed to produce the coordinates at stage 1 optimization" 
            else:
                print "No dictionary file produced ! " 
        else:
            print "No dictionary file produced ! " 
                        
    def runGeoOpt(self, tRoot, tPdbIn, tCifLibIn):
       
        # Geometrical optimization
        if os.path.isfile(tPdbIn) and os.path.isfile(tCifLibIn):
            if os.path.getsize(tPdbIn) > 100 and os.path.getsize(tCifLibIn) > 50:
                stageNow = 0
                #self.runRefmac(tPdbIn, tCifLibIn, tRoot, stageNow)
                #tPdbIn1 = self.refmacXYZOUTName 
                tPdbIn1 = tPdbIn 
                #stageNow = 1
                self.runRefmac(tPdbIn1, tCifLibIn, tRoot, stageNow)
                #self.refmacXYZOUTName
                if not os.path.isfile(self.refmacXYZOUTName):
                    print "Failed to produce the coordinates for input file %s in optimization"%tPdbIn 
                else:
                    tFValue = -1.0 
                    tFValue = self.getRefmacFValue(self._log_name)
                    if tFValue >= -0.000001:
                        if len(tRoot) !=0:
                            if not self.refmacXYSList.has_key(tRoot):
                                 self.refmacXYSList[tRoot] = {}
                            self.refmacXYSList[tRoot]["log"]    = self._log_name 
                            self.refmacXYSList[tRoot]["xyz"]    = self.refmacXYZOUTName 
                            self.refmacXYSList[tRoot]["fvalue"] = tFValue 
                            # print "| %s|"%("Done: " + tRoot).ljust(64) 
                        aRefPair = [tFValue, self.refmacXYZOUTName]
                        self.refmacMinFValueList.append(aRefPair)
            else:
                print "No dictionary file for optimization ! " 
        else:
            print "No input pdb and/or dictionary files for optimization ! " 

    def runGeoOptOneMolFull(self, tIdxMol):

        print "Number of final output conformers for molecule %d is %d "%(tIdxMol+1, self.numConformers)
        tmpStr = ""
        if self.numConformers == 1:
            tmpStr = "_tmp"
        #nConf = self.rdKit.molecules[tIdxMol].GetNumConformers()
        #print "Number of intial conformers for the molecule  ", nConf
        print "Number of intial conformers for refmac geo-opt  ", len(self.rdKit.selecConformerIds)

        inPdbNamesRoot =[]
        #for idxConf in range(nConf): 
        idxC = 1
        #print "Number of selct conformers ", self.rdKit.selecConformerIds
        for idxConf in self.rdKit.selecConformerIds : 
            tPdbRoot = "mol_" + str(tIdxMol+1) + "_conf_" + str(idxC) + tmpStr
            aConfPdb = os.path.join(self.scrDir, tPdbRoot + "_init.pdb")
            #print "PDB root ", tPdbRoot
            #print aConfPdb
            self.fileConv.MolToPDBFile(aConfPdb, tIdxMol, self.rdKit.molecules[tIdxMol], self.fileConv.dataDescriptor,self.monomRoot, idxConf,  self.rdKit.repSign)
            if os.path.isfile(aConfPdb):
                inPdbNamesRoot.append(tPdbRoot)
            idxC+=1
        aLibCifIn = self.outRstCifName 
        for aFRoot in inPdbNamesRoot: 
            aPdbIn    = os.path.join(self.scrDir, aFRoot + "_init.pdb")
            #print "|%s%s|"%("Input XYZ : ".ljust(12), aPdbIn.ljust(53))
            #print "|%s%s|"%("Input LIB : ".ljust(12), aLibCifIn.ljust(53))
            self.runGeoOpt(aFRoot, aPdbIn, aLibCifIn) 
            if  self.runExitCode :
                print "Geometrical optimization fails to produce the final coordinates for %s after geometrical optimization"%aPdbIn
        if len(self.refmacMinFValueList) > 0 :
            #self.refmacMinFValueList.sort(listComp2)
            #for aPair in self.refmacMinFValueList:
            #    print "======"
            #    print "FValue: ", aPair[0], "  File name ", aPair[1]  
            if self.numConformers==1: 
                #print "Come to output final info"
                self.getFinalOutputFiles("", self.rdKit.molecules[tIdxMol], aLibCifIn, self.refmacMinFValueList[0][1], self.fileConv.ccp4DataDes,self.fileConv.strDescriptors,self.fileConv.delocBondList)
            else:
                for i in range(self.numConformers):
                    aRoot = "Mol" + str(tIdxMol) + "_conformer" + str(i+1) 
                    self.getFinalOutputFiles(aRoot, self.rdKit.molecules[tIdxMol], aLibCifIn, self.refmacMinFValueList[i][1], self.fileConv.ccp4DataDes,self.fileConv.strDescriptors,self.fileConv.delocBondList)
            
    def getRefmacFValue(self, tLogName):

        aFva = -1.0

        try:
            aLogIn = open(tLogName, "r")
        except IOError:
            print "Error %s can not be opened for reading"%tLogName
        else:
            allLs = aLogIn.readlines()
            aLogIn.close()

            for aL in allLs:
                if aL.find("fvalues") !=-1:
                    allStrs= aL.strip().split()
                    if len(allStrs)==5:
                        aFva = float(allStrs[2])
                        #print aL
                        #print "fvalue taken : ", aFva

        
        #print "Fvalue taken from %s is %10.5f"%(tLogName, aFva)     

        return aFva
        
    def getFinalOutputFiles(self, tRoot, tMol, tInCif, tInPdb, tDataDescriptor=None, tStrDescriptors=None, tDelocList=None):
       
        iIter = 0
        if os.path.isfile(tInPdb):
            if tRoot !="":
                finPdb = self.outRoot + "_"+ tRoot + ".pdb"
                finRst = self.outRoot +  "_"+ tRoot +".cif"
                tStr = tRoot.strip().split("_")[-1]
            else:
                finPdb = self.outRoot +  ".pdb"
                finRst = self.outRoot +  ".cif"
            #print "tInPdb ", tInPdb 
            #print "finPdb ", finPdb

            shutil.copy(tInPdb, finPdb)

            if os.path.isfile(finPdb):
                self.inPdbName        = finPdb 
                self.inMmCifName      = tInCif
                self.outRstCifName    = finRst
                self.transCoordsPdbToCif(self.inPdbName, self.inMmCifName, self.outRstCifName, tMol, tDataDescriptor, tStrDescriptors, tDelocList)
        else:
            print "Failed to produce %s after final geometrical optimization"%tInPdb

    def transCoordsPdbToCif(self, tPdbInName, tCifInName, tCifOutName, tMol=-1, tDataDescriptor=None, tStrDescriptors=None,tDelocList=None):

        cifCont = {}
        cifCont['head'] = []
        monoId = ""
        if self.monomRoot.find("UNL") ==-1:
            monoId = self.monomRoot
        elif tDataDescriptor:
            monoId = tDataDescriptor[-1].strip().split()[0]
        else:
            monoId = "UNL"
      
        if tDataDescriptor:
            cifCont['head']   = ["#\n", "data_comp_list\n", "loop_\n"]
            for aL in tDataDescriptor:
                cifCont['head'].append(aL+"\n")   
            cifCont['head'].append("#\n")
            # monoId = tDataDescriptor[-1].strip().split()[0]
            cifCont['head'].append("data_comp_%s\n"%monoId)
            cifCont['head'].append("#\n")
            cifCont['head'].append("loop_\n")
        
        cifCont['atoms']  = []
        cifCont['others'] = []
        cifCont['others1'] = []
        cifCont['others2'] = []
        if tDelocList:
            if len(tDelocList) !=0:
                cifCont['bonds'] = []
        pdbAtoms = {}
        try:
            tPdbIn = open(tPdbInName, "r")
        except IOError:
            print "%s can not be opened for reading"%tPdbInName
            sys.exit()
        else:
            try: 
                tCifIn = open(tCifInName, "r")
            except IOError:
                print "%s can not be opened for reading"%tCifInName
                sys.exit()
            else:
                allPdbLines  = tPdbIn.readlines()
                tPdbIn.close()
                for aLine in allPdbLines:
                    aLine = aLine.strip()
                    if aLine.find("ATOM") !=-1 or aLine.find("HETATM") !=-1 :
                        tName = aLine[12:16].strip()
                        tX    = aLine[30:38].strip()
                        tY    = aLine[38:46].strip()
                        tZ    = aLine[46:54].strip()
                        pdbAtoms[tName] = []
                        pdbAtoms[tName].append(tX)
                        pdbAtoms[tName].append(tY)
                        pdbAtoms[tName].append(tZ)
    
                allCifLines = tCifIn.readlines()
                tCifIn.close()
                lAtom   = False
                lBond   = False
                lOther  = False
                lOther1 = False
                lOther2 = False
                
                if tDataDescriptor:
                    cifCont['atoms'].append("_chem_comp_atom.comp_id\n")
                    for aLine in allCifLines:
                        if not lAtom and aLine.find("_chem_comp_atom.comp_id") != -1:
                            lAtom  = True
                        elif lAtom and aLine.find("loop") != -1:
                            lAtom  = False
                            lBond  = False
                            lOther1 = True
                            lOther2 = False
                            cifCont['others1'].append(aLine)
                        elif lBond and aLine.find("loop") != -1:
                            lAtom  = False
                            lBond  = False
                            lOther1 = False 
                            lOther2 = True
                            cifCont['others2'].append(aLine)
                        elif tDelocList and len(tDelocList)!=0 and lOther1 and aLine.find("_chem_comp_bond.comp_id") != -1:
                            lBond  = True
                            lOther1 = False
                            cifCont['bonds'].append(aLine)
                        elif lAtom:
                            strGrp = aLine.split()
                            if len(strGrp) ==1:
                                cifCont['atoms'].append(aLine)
                            else:
                                tName  = ""
                                tID  = strGrp[1].strip()
                                if tID.find("\"") !=-1:
                                    for aC in tID:
                                        if aC !="\"":
                                            tName +=aC
                                else:
                                    tName = tID
                            
                                if pdbAtoms.has_key(tName):
                                    bLine = "%s%s%s%s%s%s%s%s\n"%(monoId.ljust(8), tID.ljust(8), strGrp[2].ljust(8), \
                                                                strGrp[3].ljust(8), strGrp[4].ljust(8), \
                                                                pdbAtoms[tName][0].ljust(12), pdbAtoms[tName][1].ljust(12), \
                                                                pdbAtoms[tName][2].ljust(12)) 
                                    cifCont['atoms'].append(bLine)
                                else: 
                                    print "Bug. can not find atom %s in Pdb file %s "%(tName, tPdbInName) 
                                    sys.exit()
                        elif lBond:
                            cifCont['bonds'].append(aLine)
                        elif lOther1:
                            cifCont['others1'].append(aLine)
                        elif lOther2:
                            cifCont['others2'].append(aLine)
                else :
                    lStart = True
                    for aLine in allCifLines:
                        if aLine.find("_chem_comp_atom.z") != -1:
                            cifCont['head'].append(aLine)
                            lStart = False
                            lAtom  = True
                        elif lAtom and aLine.find("loop") != -1:
                            lAtom  = False
                            lOther = True
                            cifCont['others'].append(aLine)
                        elif lStart:  
                            cifCont['head'].append(aLine)
                        elif lAtom:
                            strGrp = aLine.split()
                            tName  = ""
                            tID  = strGrp[1].strip()
                            if tID.find("\"") !=-1:
                                for aC in tID:
                                    if aC !="\"":
                                        tName +=aC
                            else:
                                tName = tID
                            
                            if pdbAtoms.has_key(tName):
                                bLine = "%s%s%s%s%s%s%s%s\n"%(strGrp[0].ljust(8), tID.ljust(8), strGrp[2].ljust(8), \
                                                              strGrp[3].ljust(8), strGrp[4].ljust(8), \
                                                              pdbAtoms[tName][0].ljust(12), pdbAtoms[tName][1].ljust(12), \
                                                              pdbAtoms[tName][2].ljust(12)) 
                                cifCont['atoms'].append(bLine)
                            else: 
                                print "Bug. can not find atom %s in Pdb file %s "%(tName, tPdbInName) 
                                sys.exit()
                        elif lOther:
                            cifCont['others'].append(aLine)
 
                if cifCont.has_key("bonds") and len(cifCont['bonds']) !=0:
                    idxMap = {}
                    tHead  = []
                    tBLs   = []
                    iB=0
                    for aL in cifCont['bonds']:
                        strGrp = aL.strip().split()
                        if len(strGrp)==1:
                            if aL.find("_chem_comp_bond.atom_id_1") !=-1:
                                idxMap["atom1"]= iB
                                #print "atom 1 ", iB
                            elif aL.find("_chem_comp_bond.atom_id_2") !=-1:
                                idxMap["atom2"]= iB
                                #print "atom 2 ", iB
                            elif aL.find("_chem_comp_bond.type") !=-1:
                                idxMap["bondType"]= iB
                            tHead.append(aL)
                            iB+=1
                        else:
                            tBLs.append(aL)

    
                    cifCont["bonds"] = []
                    for aL in tHead:
                        cifCont["bonds"].append(aL)

                    for aB in tBLs:
                        for aDe in tDelocList:
                            aDe[0] = aDe[0].strip()
                            aDe[1] = aDe[1].strip()
                            #print strGrp[idxMap["atom1"]]
                            #print strGrp[idxMap["atom2"]]
                            strGrp = aB.strip().split()
                            if (strGrp[idxMap["atom1"]].strip()==aDe[0] and strGrp[idxMap["atom2"]].strip()==aDe[1]) \
                                or (strGrp[idxMap["atom2"]].strip()==aDe[0] and strGrp[idxMap["atom1"]].strip()==aDe[1]):
                                 strGrp[idxMap["bondType"]] = "deloc"
                                 aB = ""                                             
                                 for aS in strGrp:
                                     aB+=("%s"%aS.ljust(len(aS)+8))
                                 aB+="\n"
                                 break
                        cifCont["bonds"].append(aB)
                                  
                try:
                    tOutCif = open(tCifOutName, "w")
                except IOError:
                    print "%s can not be opened for reading"%tCifOutName
                    sys.exit()
                else:
                    if len(self.outCifGlobSect):
                        for aL in self.outCifGlobSect:
                            tOutCif.write(aL)
                            
                    for aL in cifCont['head']:
                        tOutCif.write(aL)
                    for aL in cifCont['atoms']:
                        tOutCif.write(aL)
                    if tDataDescriptor:
                        for aL in cifCont['others1']:
                            strGrp = aL.strip().split()
                            if len(strGrp)==1:
                                tOutCif.write(aL)
                            else:
                                aL1 = monoId + aL[3:]
                                tOutCif.write(aL1)
                        if cifCont.has_key('bonds'):
                            for aL in  cifCont['bonds']:
                                strGrp = aL.strip().split()
                                if len(strGrp)==1:
                                    tOutCif.write(aL)
                                else:
                                    aL1 = monoId + aL[3:]
                                    tOutCif.write(aL1)
                        for aL in cifCont['others2']:
                            strGrp = aL.strip().split()
                            if len(strGrp)==1:
                                tOutCif.write(aL)
                            else:
                                aL1 = monoId + aL[3:]
                                tOutCif.write(aL1)
                    else:
                        for aL in cifCont["others"]:
                            strGrp = aL.strip().split()
                            if len(strGrp)==1:
                                tOutCif.write(aL)
                            else:
                                aL1 = monoId + aL[3:]
                                tOutCif.write(aL1)
                    
                    if tStrDescriptors and tStrDescriptors.has_key("props") and len(tStrDescriptors["props"]) !=0:
                        tOutCif.write("loop_"+"\n")
                        for aL in tStrDescriptors["props"]:
                            tOutCif.write(aL+"\n")
                        for aL in tStrDescriptors["entries"]:
                            tOutCif.write(aL+"\n")
                    """
                    elif tStrDescriptors.has_key("defProps") and tStrDescriptors.has_key("defSmiles"):
                        for aProp in tStrDescriptors["defProps"]:
                            tOutCif.write(aProp+"\n")
                        for aL in tStrDescriptors["defSmiles"]:
                            tOutCif.write(aL+"\n")
                    else:
                        if tMol !=-1:
                            aList = ["loop_", "_pdbx_chem_comp_descriptor.comp_id", "_pdbx_chem_comp_descriptor.type", \
                                     "_pdbx_chem_comp_descriptor.program", "_pdbx_chem_comp_descriptor.program_version",\
                                     "_pdbx_chem_comp_descriptor.descriptor"]
                      
                            for aProp in aList:
                                tOutCif.write(aProp+"\n")
                                aSmi             = tMol.GetProp("SmilesOut")
                                aSmilen          = len(aSmi)
                            tOutCif.write("%s%s%s%s \"%s\"\n"%(monoId.ljust(10), "SMILES".ljust(10), "RDKit".ljust(12), "1.00".ljust(6), \
                                    aSmi.ljust(aSmilen)))           
                    """
                    tOutCif.close()
                        
    def outEnergyGeoMap(self, tIdxMol):

        aListFName = os.path.join(self.scrDir, self.baseRoot + "_energy_vs_Conformers.list")
        try: 
            aListF = open(aListFName, "w")
        except IOError:
            print "%s can not be opened for reading"%aListFName
        else:
            tmpStr = ""
            if self.numConformers == 1:
                tmpStr = "_tmp"
            aListF.write("RDKit-Energy\tOriginal-Conformer-Id\tPDB-File-Name\n")
            nId = 1
            for aEng in sorted(self.rdKit.conformerEngMap.iterkeys()):
                for aCId in self.rdKit.conformerEngMap[aEng]:
                    aPdbName = "mol_" + str(tIdxMol+1) + "_conf_" + str(nId) + tmpStr + "_init.pdb"
                    nId +=1   
                    aListF.write("%8.4f\t%d\t%s\n"%(aEng, aCId, aPdbName))

            aListF.write("REFMAC-Energy\tPDB-File-Name-With-Original-Conformer-Id\n")
            for aPair in self.refmacMinFValueList:
                aPair[1]= aPair[1].strip().split("/")[-1].strip()
                aListF.write("%8.4f\t%s\n"%(aPair[0], aPair[1]))
            aListF.close()
       
            
            
            

    def execute(self):
        
        self.printJobs()
        if self.workMode == 11 or self.workMode == 111:

            # Stage 1: dictionary generation using a mmcif file 
            if os.path.isfile(self.inMmCifName):
                if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                    print "The input system contains metal or other heavier element"
                    print "The current version deals only with the atoms in the set of 'organic' elements" 
                    sys.exit()
                self.runLibmol(self.inMmCifName)
            else:
                print "The input %s does not exist"%self.inMmCifName
                sys.exit()
            
            if self.workMode == 11 and not self.runExitCode :
                # Stage 2: optimization
                self.runGeoOpt() 
            elif not self.workMode == 111:
                self.printExitInfo() 
  
        if self.workMode == 12 or self.workMode == 121:
          
            # Stage 1: Transfer the SMILE input into a mmcif file
            if os.path.isfile(self.inSmiName):
                self.runLibcheck(self.inSmiName)
            else:
                print "%s does not exist"%self.inSmiName
                sys.exit()
            # Stage 2: dictionary generation using a mmcif file
            if os.path.isfile(self.libcheckLibName):
                if not self.chemCheck.isOrganic(self.libcheckLibName, 11):
                    print "The input system contains metal or other heavier element"
                    print "The current version deals only with the atoms in the set of 'organic' elements" 
                    sys.exit()
                self.runLibmol(self.libcheckLibName)
            else:
                print "The input smiles contain metal elements "
                sys.exit()

            if self.workMode == 12 and not self.runExitCode:
                # Stage 3: optimization
                self.runGeoOpt()            
            elif not self.workMode == 121:
                self.printExitInfo() 

        if self.workMode == 13 or self.workMode == 131 :
        
            """    
            # Stage 1: dictionary generation using  a mdl file
            if os.path.isfile(self.inMdlName):
                self.runLibmol(self.inMdlName)
            else:
                print "can not find %s to read "%self.inMdlName
                sys.exit() 
            
            if self.workMode == 13 and not self.runExitCode: 
                # Stage 2: optimization
                self.runGeoOpt()            
            elif not self.workMode == 131:
                self.printExitInfo() 
            """
  
    
        if self.workMode == 14 or self.workMode == 141 :
            
            # Stage 1: dictionary generation using  a sdf file
            if os.path.isfile(self.inSdfName):
                self.runLibmol(self.inSdfName)
            else:
                print "%s does not exist"%self.inSdfName
                sys.exit()
           
            if self.workMode == 14 and not self.runExitCode: 
                # Stage 2: optimization
                self.runGeoOpt()            
            elif not self.workMode == 141:
                self.printExitInfo() 

        if self.workMode == 15 or self.workMode == 151 :
            
            # Stage 1: dictionary generation using  a mol2 file
            if os.path.isfile(self.inMol2Name):
                self.runLibmol(self.inMol2Name)
            else:
                print "%s does not exist"%self.inMol2Name
                sys.exit()
           
            if self.workMode == 15 and not self.runExitCode: 
                # Stage 2: optimization
                self.runGeoOpt()            
            else:
                self.printExitInfo() 

        if self.workMode == 21:
            
            # Stage 1: generate molecules and the associated bond and bond-angle values 
            # using a small molecule cif file
            if os.path.isfile(self.inStdCifName):
                self.runLibmol(self.inStdCifName)
            else:
                print "Can not find the input file ", self.inStdCifName 
            
        if self.workMode == 22:
            
            # 1. Generate molecules using the small molecule cif files at the input directory. 
            # 2. Generate atom classes for atoms in the molecules.
            # 3. Obtain unique bond lengths and angles and cluster them according to their 
            #    component atoms in tables.
            if os.path.isdir(self.inStdCifDir):
                self.runLibmol()
            else:
                print "Can not find the input directory ", self.inStdCifDir
                
        if self.workMode == 31:
            print "work mode ", self.workMode
            if os.path.isfile(self.inMmCifName):
                self.runLibmol()    
        
        if self.workMode == 32:
            print "work mode ", self.workMode
            if os.path.isfile(self.inSmiName):
                self.runLibcheck(self.inSmiName)
            else:
                print "%s does not exist"%self.inSmiName
                sys.exit()

            if os.path.isfile(self.libcheckLibName):
                self.runLibmol(self.libcheckLibName)
            else:
                print "%s does not exist"%self.inMmCifName
                sys.exit()
            
        if self.workMode == 33:
            print "work mode ", self.workMode
            if os.path.isfile(self.inMdlName):
                self.runLibmol()    
        

        if self.workMode == 34:
            print "work mode ", self.workMode
            if os.path.isfile(self.inSdfName):
                self.runLibmol()    

        if self.workMode == 35:
            print "work mode ", self.workMode
            if os.path.isfile(self.inStdCifName):
                self.runLibmol()    

        if self.workMode ==111 or self.workMode ==121 or self.workMode ==131 or self.workMode ==141:
            if os.path.isfile(self.outRstCifName):
                tCif = self.outRoot + ".cif"
                #os.system("cp %s %s"%(self.outRstCifName, tCif)) 
                shutil.copy(self.outRstCifName, tCif) 
            else:
                print "acedrg failed to generate a dictionary file"       

    def executeWithRDKit(self):
        
        self.printJobs()
        self.rdKit.useExistCoords  = self.useExistCoords 
        if self.useExistCoords or self.workMode==16 or self.workMode==161:
            print "One of output conformers will using input coordinates as initial ones"
        elif self.workMode !=0 and self.workMode != 61 :
            print "Input coordinates will be ignored"

        # Stage 1: initiate a mol file for RDKit obj
       
        if self.workMode == 11 or self.workMode == 111:
            # The input file is an mmcif file     
            if os.path.isfile(self.inMmCifName) and self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                self.fileConv.mmCifReader(self.inMmCifName)
                if len(self.fileConv.dataDescriptor):
                    self.setMonoRoot(self.fileConv.dataDescriptor) 
                if len(self.fileConv.atoms) !=0 and len(self.fileConv.bonds) !=0 :
                    # Option A: 
                    if self.useExistCoords :
                        aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                        self.fileConv.MmCifToMolFile(self.inMmCifName, aIniMolName, 2)
                        if os.path.isfile(aIniMolName) :
                            if len(self.fileConv.chiralPre) !=0:
                            # Chiral centers defined in the original cif file
                                self.rdKit.chiralPre =[]
                                for aChi in self.fileConv.chiralPre:
                                    self.rdKit.chiralPre.append(aChi) 
                            self.rdKit.initMols("mol", aIniMolName, self.monomRoot, \
                                                self.chemCheck, self.inputPara["PH"], self.numConformers, 0,\
                                                self.fileConv.nameMapingCifMol, self.fileConv.inputCharge) 
                    elif self.fileConv.strDescriptors.has_key("props") \
                       and self.fileConv.strDescriptors.has_key("entries"):
                        iProp = 0
                        iType = -1
                        iDes  = -1
                        iProg = -1
                        lSmi  = False
                        aSmi  = ""
                        for aProp in self.fileConv.strDescriptors["props"]:
                            if aProp.find("_pdbx_chem_comp_descriptor.type") !=-1:
                                iType = iProp
                            elif aProp.find("_pdbx_chem_comp_descriptor.descriptor") !=-1:
                                iDes = iProp
                            elif aProp.find("_pdbx_chem_comp_descriptor.program") !=-1\
                                 and aProp.find("_version")==-1:
                                iProg = iProp
                            iProp +=1
                        #print "prog col ", iProg
                        if iType > 0  and iDes > 0 :
                            # Get Canonical one 
                            for aEnt in self.fileConv.strDescriptors["entries"]:
                                strGrp = []
                                if aEnt.find("\"") !=-1:
                                   strGrpT = []
                                   lStart = False
                                   lEnd   = False
                                   strT = ""
                                   for aC in aEnt:
                                       if lStart == False and aC=="\"":
                                           lStart = True
                                           if len(strT.strip()) !=0:
                                               strGrpT.append(strT)
                                               strT = ""
                                           strT+=(aC)
                                       elif lStart == True and aC=="\"":
                                           lStart = False
                                           strT+=(aC)
                                           strGrpT.append(strT)
                                           strT = ""
                                       else:
                                           strT +=(aC) 
                                   if len(strT) !=0:
                                       strGrpT.append(strT)
                                       strT = ""
                                   if len(strGrpT) !=0:
                                       for aElem in strGrpT:
                                           if aElem[0].find("\"") ==-1:
                                               strGrpTT=aElem.strip().split()
                                               for aSubE in strGrpTT:
                                                   strGrp.append(aSubE)
                                           else:
                                               strGrp.append(aElem)
                                else:
                                    strGrp = aEnt.strip().split()
                                #print strGrp
                                #print strGrp[iType].upper()
                                if len(strGrp)==len(self.fileConv.strDescriptors["props"]):
                                    if strGrp[iType].upper().find("CANONICAL") !=-1\
                                        and strGrp[iProg].upper().find("OPENEYE") !=-1:
                                        if strGrp[iDes][0].find('\"') !=-1 or strGrp[iDes][0].find("\'") !=-1: 
                                            aSmi = strGrp[iDes][1:-1]
                                        else:
                                            aSmi = strGrp[iDes][0:]
                                        lSmi = True    
                                        break
                            if not lSmi:   # Get non-openEye and "CANONICAL" SMILES
                                for aEnt in self.fileConv.strDescriptors["entries"]:
                                    strGrp = aEnt.strip().split()
                                    if len(strGrp)==len(self.fileConv.strDescriptors["props"]):
                                        if strGrp[iType].upper().find("CANONICAL") !=-1 :
                                            if strGrp[iDes][0].find('\"') !=-1 or strGrp[iDes][0].find("\'") !=-1: 
                                                aSmi = strGrp[iDes][1:-1]
                                            else:
                                                aSmi = strGrp[iDes][0:]
                                            lSmi = True    
                                            break
                            if not lSmi:   # Get any smiles
                                for aEnt in self.fileConv.strDescriptors["entries"]:
                                    strGrp = aEnt.strip().split()
                                    if len(strGrp)==len(self.fileConv.strDescriptors["props"]):
                                        if strGrp[iType].upper().find("SMILES") !=-1:
                                            lSmi = True    
                                            aSmi = strGrp[iDes][1:-1]
                                            break
                       
                        if lSmi :
                            print "Smiles str  ", aSmi
                            aIniSmiName = os.path.join(self.scrDir, self.baseRoot + "_init.smi")
                            #print "Smiles file ", aIniSmiName
                            #print "Smiles str  ", aSmi
                            fSmi = open(aIniSmiName, "w")
                            fSmi.write(aSmi+"\n")
                            fSmi.close()
                            self.rdKit.reSetSmi = False
                            self.rdKit.initMols("smi", aIniSmiName, self.monomRoot,\
                                                self.chemCheck, self.inputPara["PH"], self.numConformers, 0)
                            if len(self.rdKit.molecules) !=0:
                                self.workMode = 31
                                tSmiMmcifName = os.path.join(self.scrDir, "tSmi.cif")
                                self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[0], tSmiMmcifName, self.chemCheck, self.monomRoot)
                                self.runLibmol(tSmiMmcifName, 1)
                                # Mode 2           
                                aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                                if os.path.isfile(self.inMmCifName):
                                    aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                                    self.fileConv.MmCifToMolFile(self.inMmCifName, aIniMolName, 2)
                                    if os.path.isfile(aIniMolName) :
                                        self.rdKit.initMols("mol", aIniMolName, self.monomRoot, \
                                             self.chemCheck, self.inputPara["PH"], self.numConformers, 2,\
                                             self.fileConv.nameMapingCifMol, self.fileConv.inputCharge) 
                                        tMolMmcifName = os.path.join(self.scrDir, "tSmi2.cif")
                                        self.rdKit.MolToSimplifiedMmcif(self.rdKit.moleculesB[0], tMolMmcifName, self.chemCheck, self.monomRoot)
                                        self.runLibmol(tMolMmcifName)
                            self.libmolAT1     = os.path.join(self.scrDir, "atomTypes_1.txt")
                            self.libmolAT2     = os.path.join(self.scrDir, "atomTypes_2.txt")
                            self.libmolMatched = os.path.join(self.scrDir, "matchedType.txt")

                            if os.path.isfile(self.libmolAT1) and os.path.isfile(self.libmolAT2):
                                 self.workMode = 900
                                 self.runLibmol()
                                 self.rdKit.reSetSmi = True
                                 self.rdKit.molecules = []
                                 self.rdKit.initMols("smi", aIniSmiName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)  
                                 self.fileConv.mergeAtomNames(self.libmolMatched, self.rdKit.molecules[0])
                                 self.fileConv.addAtomOrigChiralSign(self.rdKit.molecules[0])
                                 self.useCifCoords = True
                            self.workMode = 11
                    else: 
                        aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                        if os.path.isfile(self.inMmCifName):
                            aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                            self.fileConv.MmCifToMolFile(self.inMmCifName, aIniMolName, 1)
                            if os.path.isfile(aIniMolName) :
                                if len(self.fileConv.chiralPre) !=0:
                                    # Chiral centers defined in the original cif file
                                    self.rdKit.chiralPre =[]
                                    for aChi in self.fileConv.chiralPre:
                                        self.rdKit.chiralPre.append(aChi) 
                                    #self.rdKit.reSetChirals = True
                                self.rdKit.initMols("mol", aIniMolName, self.monomRoot, self.chemCheck, self.inputPara["PH"],\
                                                self.numConformers, 0, self.fileConv.nameMapingCifMol,\
                                                self.fileConv.inputCharge) 
                    
        if self.workMode==51 : 
            aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
            if os.path.isfile(self.inMmCifName):
                aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                self.fileConv.MmCifToMolFile(self.inMmCifName, aIniMolName)
                if os.path.isfile(aIniMolName) :
                    self.rdKit.initMols("mol", aIniMolName, self.monomRoot, self.chemCheck, self.inputPara["PH"],\
                                         self.numConformers, 0, self.fileConv.nameMapingCifMol,\
                                         self.fileConv.inputCharge) 
                    
        if self.workMode == 12 or self.workMode == 121 or self.workMode==52 :
          
            # The input file is  a SMILES file
            if os.path.isfile(self.inSmiName) and self.chemCheck.isOrganic(self.inSmiName, self.workMode):
                self.rdKit.reSetSmi = True
                self.rdKit.initMols("smi", self.inSmiName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)
                if len(self.rdKit.monoName) !=0:
                    self.monomRoot = self.rdKit.monoName  
         
        if self.workMode == 13 or self.workMode == 131 or self.workMode==53:
                
            # The input file is  a mol file
            if os.path.isfile(self.inMdlName) and self.chemCheck.isOrganic(self.inMdlName, self.workMode):
                aTmpMoleFile = os.path.join(self.scrDir, self.baseRoot + "_edited.mol")
                self.fileConv.CheckElemSymbolsInMolFile(self.inMdlName, aTmpMoleFile)
                self.rdKit.initMols("mol", aTmpMoleFile, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)           
                        
        if self.workMode == 14 or self.workMode == 141 or self.workMode==54 :
            
            # The input file is a sdf file
            if os.path.isfile(self.inSdfName) and self.chemCheck.isOrganic(self.inSdfName, self.workMode):
                self.rdKit.initMols("sdf", self.inSdfName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)

        if self.workMode == 15 or self.workMode == 151 or self.workMode==55 :
            
            # The input file is a mol2 file
            if os.path.isfile(self.inMol2Name):
                if self.chemCheck.isOrganic(self.inMol2Name, self.workMode):
                    self.runLibmol(self.inMol2Name)
                    print self.outRstCifName
                    if not self.runExitCode and os.path.isfile(self.outRstCifName):
                        aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                        if os.path.isfile(self.outRstCifName):
                            aIniMolName = os.path.join(self.scrDir, self.baseRoot + "_initTransMol.mol")
                            #print aIniMolName
                            self.fileConv.MmCifToMolFile(self.outRstCifName, aIniMolName)
                            if len(self.fileConv.atoms) !=0 and len(self.fileConv.bonds) !=0 \
                               and os.path.isfile(aIniMolName) :
                                self.rdKit.initMols("mol", aIniMolName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers, 0, self.fileConv.nameMapingCifMol)    
            else: 
                print self.inMol2Name, " can not be found for reading "
                sys.exit()

        if self.workMode == 16 or self.workMode == 161 :
            
            # The input file is a pdb file
            # if os.path.isfile(self.inLigandPdbName) and self.chemCheck.isOrganic(self.inLigandPdbName, self.workMode):
            if os.path.isfile(self.inLigandPdbName):

                # RDKit way of generating molecules
                #self.fileConv.getAtomNamesInPDB(self.inLigandPdbName)
                #self.rdKit.initMols("pdb", self.inLigandPdbName, self.monomRoot, self.chemCheck,\
                #                    self.inputPara["PH"], self.numConformers, 0, self.fileConv.nameMapingPDBMol)

                # Refmac + RDKit way of generating molecules
                bRoot = os.path.basename(self.inLigandPdbName)
                pRoot = bRoot.strip().split(".")[0]
                self.iniLigandPdbName = os.path.join(self.scrDir, pRoot + ".pdb") 
                self.fileConv.checkAndAddCryst1InPDB(self.inLigandPdbName, self.iniLigandPdbName)
                self.monomRoot = self.fileConv.getResNameFromPDB(self.inLigandPdbName)
                if os.path.isfile(self.iniLigandPdbName):
                    curStage = 0
                    aLibIn   = ""
                    self.runRefmac(self.iniLigandPdbName, aLibIn, self.monomRoot, curStage)
                    #print "initial input cif is  ", self.refmacXYZOUTName
                else:
                    print "Can not add line with 'CRYST1' to the temp PDB file ", tPDBName
                    sys.exit() 

                if os.path.isfile(self.refmacXYZOUTName):
                    self.inMmCifName = self.refmacXYZOUTName
                    self.useExistCoords = True
                    """
                    aOutLibCif = self.refmacXYZOUTName
                    self.refmacXYZOUTName = ""
                    self.fileConv.mmCifReader(aOutLibCif)
                    if len(self.fileConv.atoms) !=0 and len(self.fileConv.bonds) !=0 :
                        self.useExistCoords = True
                        #self.fileConv.MmCifToMolFile(self.inMmCifName, aIniMolName, 2)
                        if os.path.isfile(aIniMolName) :
                            if len(self.fileConv.chiralPre) !=0:
                            # Chiral centers defined in the original cif file
                                self.rdKit.chiralPre =[]
                                for aChi in self.fileConv.chiralPre:
                                    self.rdKit.chiralPre.append(aChi) 
                            self.rdKit.initMols("mol", aIniMolName, self.monomRoot, \
                                                self.chemCheck, self.inputPara["PH"], self.numConformers, 0,\
                                                self.fileConv.nameMapingCifMol, self.fileConv.inputCharge)    
                    """
                else: 
                    print "Failed to generate initial dictionary file ", self.refmacXYZOUTName
                    sys.exit()

        if self.workMode in [11, 12, 13, 14, 15]:
            self.workMode = 11
        elif self.workMode in [111, 121, 131, 141, 151]:
            self.workMode = 111
        if self.workMode in [51, 52, 53, 54, 55]:
            self.workMode = 51
        if self.workMode in [11,  111, 51]:
            #print len(self.rdKit.molecules)
            if len(self.rdKit.molecules):
                print "Ligand ID ", self.monomRoot
                self.fileConv.getCCP4DataDescritor(self.rdKit.molecules[0],  self.chemCheck, self.monomRoot)
                self.rdKit.hasCCP4Type = self.fileConv.hasCCP4Type
            for iMol in range(len(self.rdKit.molecules)):
                self.inMmCifName =  os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(iMol) + ".cif")
                self.initMmcifMolMap[iMol] = self.inMmCifName
                #if self.monomRoot.upper() in self.chemCheck.aminoAcids:
                #    self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[iMol], self.inMmCifName, self.chemCheck, self.monomRoot, "L-peptide")
                #else:
                if self.workMode in [11,  111]:
                    print "Using coords ", self.rdKit.useExistCoords
                self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[iMol], self.inMmCifName, self.chemCheck, self.monomRoot, self.fileConv.chiralPre)
                if os.path.isfile(self.inMmCifName):
                    if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                        print "The input system contains metal or other heavier element"
                        print "The current version deals only with the atoms in the set of 'organic' elements" 
                        sys.exit()
                    self.runLibmol(self.inMmCifName, iMol)
                else:
                    print "The input %s does not exist"%self.inMmCifName
                    sys.exit()

                if self.workMode == 11:
                    if  not self.runExitCode :
                        # Stage 2: optimization
                        print "===================================================================" 
                        print "| Geometrical Optimization                                        |"
                        print "===================================================================" 
                        
                        if len(self.rdKit.molecules) != 0 and os.path.isfile(self.outRstCifName):
                            inPdbNamesRoot = {} 
                            for idxMol in range(len(self.rdKit.molecules)): 
                                if not inPdbNamesRoot.has_key(idxMol):
                                    inPdbNamesRoot[idxMol] = []
                                print "Number of atoms in molecule %d is %d "%(idxMol+1, self.rdKit.molecules[idxMol].GetNumAtoms())
                                #nConf = self.rdKit.molecules[idxMol].GetNumConformers()
                                self.runGeoOptOneMolFull(idxMol)
                                if not self.useExistCoords:
                                    self.outEnergyGeoMap(idxMol)
                    else:
                        print "Error: No dictionary cif file is generated by Acedrg "

        elif self.workMode == 16 or self.workMode == 161:
            if os.path.isfile(self.inMmCifName):
                if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                    print "The input system contains metal or other heavier element"
                    print "The current version deals only with the atoms in the set of 'organic' elements" 
                    sys.exit()
                self.runLibmol()
            else:
                print "The input %s does not exist"%self.inMmCifName
                sys.exit()
            if self.workMode == 16: 
                if not self.runExitCode :
                    # Stage 2: optimization
                    print "===================================================================" 
                    print "| Geometrical Optimization                                        |"
                    print "===================================================================" 
                        
                    if os.path.isfile(self.outRstCifName):
                        self.refmacXYZOUTName = ""
                        self.runRefmac(self.iniLigandPdbName, self.outRstCifName, self.baseRoot, 1)
                        if not self.runExitCode and os.path.isfile(self.refmacXYZOUTName):
                            finPdb = self.outRoot + ".pdb"
                            finRst = self.outRoot + ".cif"
                            shutil.copy(self.refmacXYZOUTName, finPdb)
                            if os.path.isfile(finPdb):
                                self.inPdbName        = finPdb
                                self.inMmCifName      = self.outRstCifName
                                self.outRstCifName    = finRst
                                self.transCoordsPdbToCif(self.inPdbName, self.inMmCifName, self.outRstCifName)
                                #print "==================================================================="
                            else:
                                print "Failed to produce %s after final geometrical optimization"%finPdb

        if self.workMode == 21:
            
            # Stage 1: generate molecules and the associated bond and bond-angle values 
            # using a small molecule cif file
            if os.path.isfile(self.inStdCifName):
                self.runLibmol(self.inStdCifName)
            else:
                print "Can not find the input file ", self.inStdCifName 
            
        if self.workMode == 22:
            
            # 1. Generate molecules using the small molecule cif files at the input directory. 
            # 2. Generate atom classes for atoms in the molecules.
            # 3. Obtain unique bond lengths and angles and cluster them according to their 
            #    component atoms in tables.
            if os.path.isdir(self.inStdCifDir):
                self.runLibmol()
            else:
                print "Can not find the input directory ", self.inStdCifDir
                
        if self.workMode == 31:
            print "work mode ", self.workMode
            if os.path.isfile(self.inMmCifName):
                self.runLibmol()    
        
        if self.workMode == 32:
            print "work mode ", self.workMode
            if os.path.isfile(self.inSmiName):
                self.runLibcheck(self.inSmiName)
            else:
                print "%s does not exist"%self.inSmiName
                sys.exit()

            if os.path.isfile(self.libcheckLibName):
                self.runLibmol(self.libcheckLibName)
            else:
                print "%s does not exist"%self.inMmCifName
                sys.exit()
            
        if self.workMode == 33:
            if os.path.isfile(self.inMdlName) and self.chemCheck.isOrganic(self.inMdlName, self.workMode):
                self.rdKit.initMols("mol", self.inMdlName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)
                tMolMmcifName = os.path.join(self.scrDir, "tmpMol.cif")
                if len(self.rdKit.molecules)==1:
                    self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[0], tMolMmcifName, self.chemCheck, self.monomRoot)
                    print "Number of molecules ", len(self.rdKit.molecules)
                    self.runLibmol(tMolMmcifName)
                else:
                    print "Do not process atom-types of multiple molecules at the same time, do nothing"

        if self.workMode == 34:
            if os.path.isfile(self.inSdfName):
                self.runLibmol()    

        if self.workMode == 35:
            if os.path.isfile(self.inStdCifName):
                self.runLibmol()    
        
        if self.workMode == 61:
            aCLinkGenerator = CovLinkGenerator(self.linkInstructions, self.scrDir, self.outRoot, self.versionInfo)

        if self.workMode ==111 or self.workMode ==121 or self.workMode ==131 or self.workMode ==141:
            if os.path.isfile(self.outRstCifName):
                tCif = self.outRoot + ".cif"
                #os.system("cp %s %s"%(self.outRstCifName, tCif)) 
                shutil.copy(self.outRstCifName, tCif) 
            else:
                print "acedrg failed to generate a dictionary file"     

    def printExitInfo(self):

        print "Error: check log file at %s"%self._log_name
    
class AcedrgRDKit():

    def __init__(self, tGFFName = "", tProcessParaSet=None):

        self.molecules  = []
        self.moleculesA = []
        self.moleculesB = []

        self.funcGroupTab     = {}
        self.stdFuncGroupMols = []
        self.funcGroups       = {}       
        # key1: index of molecule, key2: stdFSmi, value: tuple of atom idx tuples 
        if tGFFName and os.path.isfile(tGFFName):
            self.setupFuncGroupTab(tGFFName)

        self.nMaxIters               = 20
        self.numRDKitOptmSteps       = 5000
        self.numInitConformers       = 500
        self.maxNumInitConformers    = 10000
        self.numConformers           = 1
        self.useExistCoords          = False
        self.conformerEngMap         = {}
        self.numSelectForRefConfs    = 25
        self.selecConformerIds       = []

        self.hasCCP4Type     = False
   
        self.chemCheck = ChemCheck()

        self.torsins   = []

        self.defaultComfId = -10

        # TMP variables
        self.smiOrig          = ""
        self.smiMod           = ""
        self.repSign          = ""
        self.reSetSmi         = False
        self.reSetChirals     = False

        self.monoName         = ""
        self.longName         = ""

    
    def setProcPara(self, tProcessParaSet):
        
        if tProcessParaSet.has_key("useExistCoords"):
            self.useExistCoords = tProcessParaSet["useExistCoords"]
        else:
            self.useExistCoords = False
  
        if tProcessParaSet.has_key("numRDKitOptmSteps"):
            self.numRDKitOptmSteps = tProcessParaSet["numRDKitOptmSteps"]
        else:
            self.numRDKitOptmSteps = 5000

        if tProcessParaSet.has_key("numInitConformers"):
            self.numInitConformers = tProcessParaSet["numInitConformers"]
        else:
            if self.useExistCoords:
                self.numInitConformers = 1
            else:
                self.numInitConformers = 500

        if tProcessParaSet.has_key("numConformers"):
            self.numConformers = tProcessParaSet["numConformers"]
        else:
            self.numConformers = 1

        if self.numConformers > self.numInitConformers:
            self.numInitConformers = self.numConformers
  
        if self.numInitConformers > self.maxNumInitConformers:
            self.numInitConformers = self.maxNumInitConformers
 
        self.numSelectForRefConfs = self.numInitConformers/20
        if self.numSelectForRefConfs < 25:
            self.numSelectForRefConfs = 25
        

    def setRepSign(self):

        # TEMP function, select the deleted chemical eleement among "F, CL, BR, I, AT", 
        # which will be used to solve temporally the chiral problem for 3 valence atoms,
        # where RDKit fails

        tList = [ "F", "Cl", "Br", "I", "At"] 
        atomElems = []
        aMol     = Chem.MolFromSmiles(self.smiOrig)
        if aMol:
           allAtoms = aMol.GetAtoms()
           for aAt in allAtoms:
               aSymb = aAt.GetSymbol().strip()
               if not aSymb in atomElems:
                   atomElems.append(aSymb)

           for aE in tList:
               if not aE in atomElems:
                   self.repSign = aE
                   break
           #print "Elements in the molecule : ", atomElems
           #print "ReplacedSymbol           : ", self.repSign 

  
    def modifySmiTmp(self, tMode=0):

        self.setRepSign()
  
        if len(self.smiOrig) !=0:
            unitTmp = ""
            lInUT = False
            lBr   = False  
   
            for aStr in self.smiOrig:
                if aStr.find("[") !=-1:
                    lInUT=True
                    self.smiMod +=aStr
                elif lInUT:
                    if aStr.find("N") !=-1 :
                        unitTmp +=aStr
                        lBr  =True
                    elif lBr:
                        if aStr.find("@") !=-1:
                            unitTmp +=aStr
                        elif aStr.find("H") !=-1:
                            unitTmp +=aStr
                        elif aStr.isdigit() :
                            unitTmp +=aStr
                        elif aStr.find("+") !=-1 or aStr.find("+") !=-1:    # stop and put everything back
                            if len(unitTmp) !=0:
                                self.smiMod +=unitTmp
                            unitTmp = ""
                            self.smiMod += aStr
                            lInUT = False
                            lBr   = False
                        elif aStr.find("]") !=-1:
                            unitTmp +=aStr
                            lBr  = False
                    elif not lBr and len(unitTmp):
                        if aStr.isdigit():
                            unitTmp +=aStr
                            self.smiChiRep(unitTmp)
                            unitTmp=""
                        else:       
                            self.smiChiRep(unitTmp)
                            self.smiMod +=aStr
                        unitTmp=""
                        lInUT=False
                    else :
                        self.smiMod +=aStr
                        lInUT = False
                else :
                    self.smiMod +=aStr
        #print "Orignal smiles is ", self.smiOrig
        #print "New smiles is ", self.smiMod
        #if self.smiOrig != self.smiMod:
        #     print "Input SMILES modified"
        #else:
        #     print "No changes for Input SMILES"
 

    def smiChiRep(self, tUniTmp):

        #print "Initial subUnit ", tUniTmp
        tNewStr =""
        tStrs =[]
        if tUniTmp.find("N@@") !=-1:
           tStrs =tUniTmp.strip().split('N@@')
           tRepStr  = "N@+1" 
        elif tUniTmp.find("N@") !=-1 and tUniTmp.find("N@@")==-1:
           tStrs =tUniTmp.strip().split('N@')
           tRepStr  = "N@@+1" 
           
        tNewStrs = tStrs[0] + tRepStr + tStrs[1] + "(%s)"%self.repSign
        self.smiMod +=tNewStrs
        #print "New repUnit ", tRepStr
        #print "New Unit ", tNewStrs
           
    def initMols(self, tFileType, tFileName, tMonoRoot, tChemCheck, tPH, tNumConf, tMode=0, tNameMap=None, tChargeList=None):

        aMolT = None
        aMolName = ""
        if tFileType == "mol" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                    aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                    aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                #print "Molecule name:  ", aMolName
                # self.reSetChirals     =    True
                aMolT                   =    Chem.MolFromMolFile(tFileName)
                #aMolT = Chem.AddHs(aMolT1)
            else:
                print "File %s does not exist "%tFileName
                sys.exit()

        elif tFileType == "mol2" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                     aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                     aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                #print "Molecule name:  ", aMolName
                aMolT     = Chem.MolFromMol2File(tFileName)
                self.chemCheck.addHs(aMolT)
            else:
                print "File %s does not exist "%tFileName
                sys.exit()

        elif tFileType =="smi" :
            if os.path.isfile(tFileName):
                # SMILES string in a file
                try:
                    fSmi = open(tFileName, "r")
                except IOError:
                    print  "% can not be open for reading "%tFileName
                    sys.exit()
                else:
                    aSmiStrT = fSmi.read()
                    strGrp = aSmiStrT.strip().split()
                    if len(strGrp) > 0:
                        aSmiStr = strGrp[0].strip()
                        if len(strGrp) > 1:
                            strGrp[1] = strGrp[1].strip()
                            if len(strGrp[1]) == 3:
                                self.monoName = strGrp[1]
                            else:
                                self.longName = strGrp[1]
     
                    else:
                        print "SMILES string formation error "
                        sys.exit()
                        
            else:
                # SMILES string in from a commandline  
                aSmiStr = tFileName.strip()
            if len(aSmiStr):
                self.smiOrig = aSmiStr
                print self.smiOrig
                if self.reSetSmi:
                    self.modifySmiTmp()
                    self.smiMod = self.smiMod.strip()
                    aMolT     = Chem.MolFromSmiles(self.smiMod.strip())
                else:
                    aMolT     = Chem.MolFromSmiles(self.smiOrig)
                if aMolT:
                    aMolT.SetProp("SmilesIn", self.smiOrig)
                else:
                    print "No molecule is generated using SMILES str ", self.smiOrig 
                    sys.exit()
        elif tFileType == "sdf" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                    aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                    aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                #print "Molecule name:  ", aMolName
                aMolT     = Chem.MolFromMolFile(tFileName)
                #aMolT = Chem.AddHs(aMolT1)
            else:
                print "File %s does not exist "%tFileName

        elif tFileType == "pdb" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                    aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                    aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                print "Molecule name:  ", aMolName
                # Test. H atoms removed as usual. Maybe change that later
                aMolT     = Chem.MolFromPDBFile(tFileName)
                if aMolT:
                    self.removeWater(aMolT)
        if not aMolT:
            print "Molecules can not generated  from file %s ! "%tFileName
            print "Check your file format "
            sys.exit()
        else:
            aMolT.SetProp("ResidueName", tMonoRoot)
            self.setOtherMolInfo(aMolT, tNumConf, tChemCheck, tPH, tNameMap, tMode, tChargeList)
            """
            lenA = len(self.moleculesA)
            lenB = len(self.moleculesB)
            lenG = len(self.molecules)
            if lenA !=0 and lenB !=0 and lenA==lenB :
                self.mergeAtomNames()
                self.molecules = []
            """

    def setNamesForAtomsInMol(self, tMol, tChemCheck, tNameMap, tStage=0):
   
        dictAtomTypes = {}
        dictAtomNames = {}

        if not tNameMap:
            for aA in tMol.GetAtoms():
                aElem = aA.GetSymbol()
                if not dictAtomTypes.has_key(aElem):
                    dictAtomTypes[aElem] = []
                dictAtomTypes[aElem].append(aA.GetIdx())
            
            for aElem in dictAtomTypes.keys():
                i = 1
                for aIdx in dictAtomTypes[aElem]:
                    aName = aElem + str(i)
                    dictAtomNames[aIdx] = aName
                    i = i+1

            for aAtom in tMol.GetAtoms():
                aAtom.SetProp("Name", dictAtomNames[aAtom.GetIdx()]) 
  
        else:
            if tStage ==0:
                # for non-H atoms
                for aAtom in tMol.GetAtoms():
                    if aAtom.GetSymbol() != "H":
                        aAtom.SetProp("Name", tNameMap["nonH"][aAtom.GetIdx()]) 
            elif tStage ==1:
                if tMol.GetProp("ResidueName") in tChemCheck.aminoAcids:
                    self.setNamesForHAtomsInMol_PP(tMol, tChemCheck)
                else:
                    self.setNamesForHAtomsInMol(tMol, tNameMap, tChemCheck)
         

    def setNamesForHAtomsInMol(self, tMol, tNameMap, tChemCheck):
      
        #tIdxHs = {}
        #tIdxHs = []
        HConns  = {}
        allAtoms = tMol.GetAtoms()
        #for aAtom in allAtoms:
        #    if aAtom.GetSymbol() !="H":
        #        print "Atom ", aAtom.GetProp("Name"), " of serial number ", aAtom.GetIdx()
        nh = 0
        for aAtom in allAtoms:
            if aAtom.GetSymbol() =="H":
                nh +=1
                idxH = aAtom.GetIdx()
                #print "H Atom  of serial number ", aAtom.GetIdx()
                idxB  = -1
                idxE  = -1
                idxC  = -1
                aSetBonds = aAtom.GetBonds()
                if len(aSetBonds)==1: 
                    idxB = aSetBonds[0].GetBeginAtomIdx()
                    idxE = aSetBonds[0].GetEndAtomIdx()
                    if idxH == idxB:
                        idxC = idxE
                    elif idxH==idxE: 
                        idxC = idxB
                    #print "Bonding  to ", " atom ", tMol.GetAtomWithIdx(idxC).GetProp("Name")
                    if idxC !=-1 :
                        nonH_Id = tNameMap["nonH"][idxC] 
                        if not HConns.has_key(nonH_Id):
                            HConns[nonH_Id] = []
                        HConns[nonH_Id].append(idxH)  
                    else:
                        print "Can not find the non-H atom that connects to H atom of serial number ", idxH   
                else:
                    print "One H atom connect more than one atoms, check"
                    sys.exit()
        #print "Total number of H atoms is ", nh
        #for aKey in HConns.keys():
        #    print "Atom ", aKey, " bonds to ", len(HConns[aKey]), " H atoms "

        if len(HConns.keys()) !=0:  
            # Check if total numbers of H atoms are different between the original mol and current mol  
            # Names for H atoms in the original file
            origHNames =[]
            for aNonH in HConns.keys():
                if  tNameMap["H"].has_key(aNonH):
                    c3 = len(tNameMap["H"][aNonH])
                    for i in range(c3):
                        origHNames.append(tNameMap["H"][aNonH][i])
            numOrigH = len(origHNames)
            nExtra = numOrigH + 1
            for aK in HConns.keys():
                c1 =  len(HConns[aK])
                if tNameMap["H"].has_key(aK):
                    c2 = len(tNameMap["H"][aK])
                    if c1 == c2:
                        for i in range(c1):
                            tMol.GetAtomWithIdx(HConns[aK][i]).SetProp("Name",tNameMap["H"][aK][i])
                    elif c1 < c2:
                        for i in range(c1):
                            tMol.GetAtomWithIdx(HConns[aK][i]).SetProp("Name",tNameMap["H"][aK][i])
                    elif c1 > c2:
                        #print "Number of H atoms bond to atom ", aK, " has changed "
                        for i in range(c2):
                            tMol.GetAtomWithIdx(HConns[aK][i]).SetProp("Name",tNameMap["H"][aK][i])
                            #print "A H atom has been set to existing name ", tMol.GetAtomWithIdx(HConns[aK][i]).GetProp("Name")   
                        # Decide the root section of additional H atoms
                        tRootId =""
                        hRootId =""
                        for aC in aK:
                            if not aC.isdigit():
                                tRootId = tRootId + aC
                        if len(tRootId) !=2:
                            hRootId ="H"
                        else:
                            hRootId = "H" + tRootId[1:]
                        print "H rootName is ", hRootId         
                        # Find largest digit number of existing H atoms attached to 
                        # the same non-H atom.
                        idxMax = 0 
                        for aEH in tNameMap["H"][aK]:
                            iNum = 0
                            for aChar in aEH:
                                if aChar.isdigit():
                                    break
                                else:
                                    iNum+=1
                            if iNum < len(aEH) and aEH[iNum:].isdigit():
                                aInt = int(aEH[iNum:])
                                if aInt > idxMax:
                                    idxMax = aInt
                        #print "idxMax ", idxMax        
                        cDiff = c1 - c2
                        #print "cDiff is ", cDiff
                        for j in range(cDiff):
                            tHName = hRootId + str(idxMax+1)
                            if not tHName in origHNames: 
                                tMol.GetAtomWithIdx(HConns[aK][c2+j]).SetProp("Name",tHName)
                                aNewHName = tMol.GetAtomWithIdx(HConns[aK][c2+j]).GetProp("Name")
                                print "An added H is named as ", aNewHName 
                                origHNames.append(aNewHName)
                            else:
                                nExtra = 2
                                tHName = hRootId + str(idxMax+nExtra)
                                while True:
                                    if not tHName in origHNames:
                                        tMol.GetAtomWithIdx(HConns[aK][c2+j]).SetProp("Name",tHName)
                                        origHNames.append(tHName)
                                        print "An added H is named as ", tHName 
                                        nExtra +=1
                                        break
                                    else:
                                        nExtra +=1
                                        tHName = hRootId + str(idxMax+nExtra)
        
                else:
                    hRootId = "H" 
                    nExtra = 2
                    for i in range(c1):
                        tHName = hRootId + str(nExtra)
                        while True:
                            if not tHName in origHNames:
                                tMol.GetAtomWithIdx(HConns[aK][i]).SetProp("Name",tHName)
                                origHNames.append(tHName)
                                print "An added H is named as ", tHName 
                                nExtra +=1
                                break
                            else:
                                nExtra +=1
                                tHName = hRootId + str(nExtra)
        
        for aAtom in allAtoms:
            if aAtom.GetSymbol() =="H":
                aSetBonds = aAtom.GetBonds()
                if aAtom.HasProp("Name"):
                    print "\nH atom idx ", aAtom.GetIdx(), "  its matched Name ", aAtom.GetProp("Name")
                    print "It is in the following bonds: "
                    for aB in aSetBonds:
                        print "Bond ", aB.GetIdx()
                        print "Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name"))
                        print "Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name"))
                else:
                    print "\nH atom without name, its idx is ", aAtom.GetIdx()

        #if len(tIdxHs) !=0:
            
        #    hNameDone = [] 
        #    nH = 1
        #    for aIdxH in tIdxHs:
        #        aHName = "H" + str(nH)
        #        if aHName in hNameDone:
        #            nH += 1
        #            aHName = "H" + str(nH)
        #        tMol.GetAtomWithIdx(aIdxH).SetProp("Name", aHName)
        #        print "not matched name ", tMol.GetAtomWithIdx(aIdxH).GetProp("Name")
        #        print "should be name ", tNameMap[tMol.GetAtomWithIdx(aIdxH).GetIdx()]
        #       
        #        nH += 1
            
    def setNamesForHAtomsInMol_PP(self, tMol, tChemCheck):
      
        tIdxHs = {}
        tIdxHs1 = []
        
        for aAtom in tMol.GetAtoms():
             if not aAtom.HasProp("Name") or aAtom.GetProp("Name") =="":
                 if aAtom.GetSymbol() =="H":
                     tIdxHs1.append(aAtom.GetIdx())

        if len(tIdxHs1) !=0:
            
            for aIdxH in tIdxHs1:
                #print "For atom ", aIdxH
                aB = tMol.GetAtomWithIdx(aIdxH).GetBonds()
                if len(aB) ==0:
                    print "Bug: a H atom of index %d does not bond to any atom "%aIdxH
                    sys.exit()
                elif len(aB) >1 : 
                    print "Bug: a H atom of index %d bond to more than one atom "%aIdxH
                    sys.exit()
                else:
                    tIdxBH = -1
                    tBH1 = aB[0].GetBeginAtomIdx()
                    tBH2 = aB[0].GetEndAtomIdx()
                    if aIdxH == tBH1:
                        tIdxBH = tBH2
                    elif aIdxH == tBH2:
                        tIdxBH = tBH1
                    else:
                        print "Bug: a H atom of index %d is not in the bond of index %d "%(aB.GetIdx())
                        sys.exit()

                    # print "It bonds atom ", tIdxBH
   
                    if not tIdxHs.has_key(tIdxBH):
                        tIdxHs[tIdxBH] =[]
                    tIdxHs[tIdxBH].append(aIdxH)
         
            hNameDone = []           
            #elemList = ["C", "N", "S", "O"]
            elemList = ["C"]
            charList = ["A", "B", "C", "D", "E", "F"]
            numList  = [ "3", "2", "1"]
            numListCN  = [ "", "2", "3"]
            numListCN2  = [ "1", "2", "3", "4", "5", "6"]
            lH       = False
            lNumOnly = False
            rName = tMol.GetProp("ResidueName")
         
            for aIdxBH in tIdxHs.keys():
                twoP = self.getTwoParts(tMol.GetAtomWithIdx(aIdxBH).GetProp("Name"))
                #print "twoP[0] : ", twoP[0]
                #print "twoP[1] : ", twoP[1]
                #print "LH ", lH
                aHName = ""
                mainElem = tMol.GetAtomWithIdx(aIdxBH).GetSymbol()
                if mainElem in elemList :
                    if len(twoP[0])==1 and len(twoP[1]) ==0:
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            if not lH:
                                aHName = "H"
                                tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                                lH = True
                            else:     
                                aHName = "H1"
                                tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                            hNameDone.append(aHName)
                        else :
                            i=0
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                aHName = "H" + numListCN[i]
                                if not aHName in hNameDone:
                                   tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                                else:
                                   i +=1 
                                   tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                                hNameDone.append(aHName)
                                if i==0:
                                    lH = True
                                i +=1 
                    elif len(twoP[0])==1 and len(twoP[1]) !=0:
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            aHName = "H" + twoP[1]
                            if not aHName in hNameDone:
                                tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                            else:
                                aHName =aHName+"G"
                                tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName+"G")
                            hNameDone.append(aHName)
                        else :
                            i=0
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                aHName = "H" + twoP[1]+ charList[i]
                                if not aHName in hNameDone: 
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                                else:
                                    aHName = "H" + twoP[1]+ charList[i+1]
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", aHName)
                                hNameDone.append(aHName)
                                i +=1
                    elif len(twoP[0]) > 1 and len(twoP[1]) ==0:
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0][1:])
                        else :
                            i=0
                            aLen = len(tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()])
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                if aLen < 4:
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0][1:]+ numList[i])
                                else:
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0][1:]+ str(aLen-i))
                                i +=1
                    elif len(twoP[0]) > 1 and len(twoP[1]) !=0:
                        name1 = tMol.GetAtomWithIdx(aIdxBH).GetProp("Name")[1:] 
                        #print "Name ", tMol.GetAtomWithIdx(aIdxBH).GetProp("Name")
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1)
                        else :
                            i=0
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                # print i
                                tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1 + numListCN2[i])
                                i +=1
                else: 
                    name1 = tMol.GetAtomWithIdx(aIdxBH).GetProp("Name")
                    name2 = ""
                    if len(name1) > 1:
                        name2 = name1[1:]
                    #print "name1 ", name1
                    #print "name2 ", name2
                    #print "twoP ", twoP 
                    if len(twoP[0]) ==1 and len(twoP[1]) ==0:
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            if lH:
                                 tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0])
                            else:
                                 tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H")
                                 lH=True
                        else :
                            i=0
                            aLen = len(tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()])
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                if mainElem.find("N") !=-1 :
                                    if not lH:
                                        
                                        tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H")
                                        lH = True
                                    else:
                                        if not lNumOnly:
                                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + numListCN2[i]) 
                                        else:
                                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0] + str(aLen-i))
                                    if i==(aLen-1) and not lNumOnly:
                                        lNumOnly = True
                                else:
                                    if aLen < 4:
                                        tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0] + numList[i])
                                    else:
                                        tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + twoP[0] + str(aLen-i))
                                i +=1 
                    elif len(twoP[0])==1 and len(twoP[1]) !=0 :
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1)
                        else :
                            i=0
                            for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                if len(name1) >2:
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name2 + charList[i])
                                else:
                                    tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1 + charList[i])
                                i +=1
                    elif len(twoP[0]) >= 2 :
                        if len(tIdxHs[aIdxBH]) ==1:
                            curIdxH= tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()][0]
                            tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1[1:])
                        else :
                             i=0
                             for curIdxH in tIdxHs[tMol.GetAtomWithIdx(aIdxBH).GetIdx()]:
                                 tMol.GetAtomWithIdx(curIdxH).SetProp("Name", "H" + name1[1:]+ numListCN2[i])
                                 i +=1
                    else:
                       print "Atom name %s can not be handled "%name1 
                       sys.exit()

    def getTwoParts(self, tStr):

        t1   = ""
        t2   = ""
        lD   = False
        for aChar in tStr:
            if aChar.isdigit():
                lD = True
            if lD:
                t2 +=aChar
            else:
                t1 +=aChar

                         
        return [t1, t2]
                    
   
    def setInitGeomOneMol(self, tMol, tConfId=-1, tMaxIters=1000, tMaxMult=20):

        # Set one initial conformers

        AllChem.EmbedMolecule(tMol)
        tFailure =0
        tFailure = self.optOneConformer(tMol, tMaxIters, tMaxMult, tConfId)
       
        return tFailure

    def setInitGeomOneMol2(self, tMol, tMaxIters=1000, tMaxMult=20):
        AllChem.EmbedMolecule(tMol)
        tFailure = 0
        tFailure=AllChem.UFFOptimizeMolecule(tMol, maxIters=tMaxIters)
        if tFailure ==1:
            nFac = 2
            while tFailure ==1 and nFac < tMaxMult:
                tFailure=AllChem.UFFOptimizeMolecule(tMol, maxIters=nFac*tMaxIters)
                nFac += 1
        return tFailure

    def optOneConformer(self, tMol, tMaxIters=1000, tMaxMult=20, tConfId=-1):
        tFailure = 0
        
        tFailure=AllChem.UFFOptimizeMolecule(tMol, maxIters=tMaxIters, vdwThresh=10.0, confId=tConfId)
        if tFailure ==1:
            nFac = 2
            while tFailure ==1 and nFac < tMaxMult:
                maxIters=nFac*tMaxIters
                tFailure=AllChem.UFFOptimizeMolecule(tMol, maxIters)
                nFac += 1
        return tFailure

    def setInitConformersAllMols(self, tNConfId=1, tConfId=-1, tMaxIters=1000, tMaxMult=20):

        tMols = []
        for aMolT in self.molecules:
            Chem.SanitizeMol(aMolT)
            aMol = Chem.AddHs(aMolT)
            
            #  Setup conformers for the molecules
            if tNConfId == 1:
                #  One initail conformer 
                lFailure =self.optOneConformer(aMol)
                # rdmolfiles.MolToPDBFile(aMol, "Test_a.pdb")
                lFailure = self.setInitGeomOneMol(aMol)
                if not lFailure:
                    # Setup chiral information associated with this conformer 
                    rdmolops.AssignAtomChiralTagsFromStructure(aMol)
                else:
                    print "Geometry of molecule %d can not be optimized within %d circles "%(len(self.allMols), self.nMaxMult*self.nMaxIters)
            elif tNConfId > 1:
                #  multiple conformers
                confIds =AllChem.EmbedMultipleConfs(aMol, tNConfId)
                for aId in confIds:
                    aFailure = self.optOneConformer(aMol, aId)
                    if not aFailure:
                        rdmolops.AssignAtomChiralTagsFromStructure(aMol, aId)
                    else:
                        print "Conformer: ", aId, " is not optimized"
            tMols.append(aMol)

        if len(tMols) != len(self.molecules):
            print "Bug in setInitConformersAllMols "
            sys.exit()
        else:
            self.molecules = []
            for aM in tMols:
                self.molecules.append(aM) 
           
    def setInitConformersOneMol(self, tMol):
            
        #  Setup conformers for the molecules
        if self.useExistCoords:
            confIds =AllChem.EmbedMultipleConfs(tMol, self.numInitConformers, maxAttempts=0, randomSeed=-1, clearConfs=False)
            tReq    = self.numInitConformers
            nNewCon = len(confIds)
            while nNewCon==0 and tReq <=5 :
                tReq +=1
                confIds =AllChem.EmbedMultipleConfs(tMol, tReq, maxAttempts=0, randomSeed=-1, clearConfs=False)
                nNewCon = len(confIds)             
                print tReq
            #for aConf in tMol.GetConformers():
            #    print "!!!!!!!!!! Conf ", aConf.GetId()
            #    for aAtom in tMol.GetAtoms():
            #        aIdx = aAtom.GetIdx()
            #        name  = aAtom.GetProp("Name") 
            #        aPos = aConf.GetAtomPosition(aIdx)
            #        print "For atom ", name
            #        print "Its coordinates are : "
            #        print "x:  ", aPos.x           
            #        print "y:  ", aPos.y           
            #        print "z:  ", aPos.z           
        else:
            confIds =AllChem.EmbedMultipleConfs(tMol, self.numInitConformers, maxAttempts=0, randomSeed=-1, clearConfs=True)
        print "Number of initial conformers requested", self.numInitConformers
        print "Number of number of opt step requested for each conformer ", self.numRDKitOptmSteps
        print "Number of new conformers ", len(confIds)
        nConf = tMol.GetNumConformers()
        print "Number of initial conformers obtained", nConf
        allConfs = tMol.GetConformers()
        allAtoms = tMol.GetAtoms()
        print Chem.MolToMolBlock(tMol)
        for aConf in allConfs:
            aCIdx =  aConf.GetId()
            aWCoordList = []
            lNorm = self.checkH_Abnormal(aConf, allAtoms, aWCoordList)
            if not lNorm:
                print "Conformer ", aCIdx, " has abormal coordinates."
                self.forceH_coords(aConf, allAtoms, aWCoordList)
            aFailure = self.optOneConformer(tMol, self.numRDKitOptmSteps, self.nMaxIters, aCIdx)
            if not aFailure:
                #print "opted id ", aId
                aForceField = AllChem.UFFGetMoleculeForceField(tMol, confId=aCIdx)
                aEng        = aForceField.CalcEnergy()
                if not self.conformerEngMap.has_key(aEng):
                    self.conformerEngMap[aEng] = []
                self.conformerEngMap[aEng].append(aCIdx)
                print "Conf : ", aCIdx, " Engergy : ", aEng
                rdmolops.AssignAtomChiralTagsFromStructure(tMol, aCIdx)
        if len(self.conformerEngMap):
            print "Current conformers have %d energy levels from UFF force field "%len(self.conformerEngMap)
            print "They are : "
            for aEng in sorted(self.conformerEngMap.iterkeys()):
                print "Energy ", aEng
                for aCid in self.conformerEngMap[aEng]:
                    print aCid

        print "The following conformers are selected for refinement: "
        nSelect = 2
        if not self.useExistCoords: 
            nSelect = self.numSelectForRefConfs
        else:
            if tMol.GetNumConformers() <2:
                nSelect = self.numSelectForRefConfs

        nID= 0
        for aEng in sorted(self.conformerEngMap.iterkeys()):
            for aCId in self.conformerEngMap[aEng]:
                if nID < self.numSelectForRefConfs:
                    self.selecConformerIds.append(aCId)
                    print "Conformer ID: ", aCId, " UFF energy : ", aEng
                    nID +=1
                else:
                    break
        print "Number of conformers selected for refinement is ",  len(self.selecConformerIds)

    def checkH_Abnormal(self, tConf, tAtoms, tWrongCoordsMap):

        lN = True
        for aAtom in tAtoms:
            aIdx = aAtom.GetIdx()
            aPos = tConf.GetAtomPosition(aIdx)
            if aAtom.GetSymbol().strip().upper()=="H" and\
               (math.isnan(aPos.x) or math.isinf(aPos.x)or\
                math.isnan(aPos.y) or math.isinf(aPos.y)or\
                math.isnan(aPos.z) or math.isinf(aPos.z)):
                if not aIdx in tWrongCoordsMap:
                    tWrongCoordsMap.append(aIdx)
                lN = False
                print "for atom ", aIdx
                print "Its position is :"
                print "x = ", aPos.x
                print "Is it nan ? ", math.isnan(aPos.x)
                print "Is it inf ", math.isinf(aPos.x)
                print "y = ", aPos.y
                print "Is it nan ? ", math.isnan(aPos.y)
                print "Is it inf ", math.isinf(aPos.y)
                print "z = ", aPos.z
                print "Is it nan ? ", math.isnan(aPos.z)
                print "Is it inf ", math.isinf(aPos.z)
               
        return lN
              
    def forceH_coords(self, tConf, tAtoms, tWrongCoordList, tReadInCoordMap=None):
        
        if tReadInCoordMap: 
            for aIdx in tWrongCoordList:
                aPos   = rdGeometry.Point3D()
                if tReadInCoordMap.has_key(aIdx):
                    aPos.x = tReadInCoordMap[aIdx][0]
                    aPos.y = tReadInCoordMap[aIdx][1]
                    aPos.z = tReadInCoordMap[aIdx][2]
                else:
                    aPos.x = random.random()
                    aPos.y = random.random()
                    aPos.z = random.random()
                tConf.SetAtomPosition(aIdx, aPos)
        else:
            # Better methods required
            for aIdx in tWrongCoordList:
                aPos = rdGeometry.Point3D()
                # Need the better methods
                aPos.x = random.random()
                aPos.y = random.random()
                aPos.z = random.random()
                tConf.SetAtomPosition(aIdx, aPos)

        # Now output all coordinates in the conformer
        for aAtom in tAtoms:
            aIdx = aAtom.GetIdx()
            name  = aAtom.GetProp("Name") 
            aPos = tConf.GetAtomPosition(aIdx)
            print "For atom ", name
            print "Its coordinates are : "
            print "x:  ", aPos.x           
            print "y:  ", aPos.y           
            print "z:  ", aPos.z           
            
    def setOtherMolInfo(self, tMol, tNumConf, tChemCheck, tPH, tNameMap, tMode=0, tChargeList=None):
        
        print "A molecule with residue name %s is generated"%tMol.GetProp("ResidueName")
        nAtoms      =  tMol.GetNumAtoms()
        #print "Number of atoms in the molecule is ", nAtoms
        initAtoms   =  tMol.GetAtoms()
        # check Atom elements
        elemList = []
        for aAtom in initAtoms:
            #print "Atom's idx ", aAtom.GetIdx()
            #print "Atom's element type ", aAtom.GetSymbol()
            elemList.append(aAtom.GetSymbol())

        #print "props in the mol ", tMol.GetPropNames()[0] 

        if not len(elemList) :
            print "No atoms in from your file, check your file format"
        elif not tChemCheck.isOrganic1(elemList):
            print "Your molecule contains METAL or other NON-ORGANIC elements "
            print "The molecule contains atoms of the following elements "
            aLine = ""
            for aElem in elemList:
                aLine.append(aElem+ "   ")
            print aLine
            sys.exit()
            
        #print tChemCheck.isOrganic1(elemList)
        #print "Number of atoms in this molecule is initially ", nAtoms
        self.setNamesForAtomsInMol(tMol, ChemCheck,  tNameMap, 0)
      
        #print "Before Sanitize the molecule"
        #print "Number of all atoms in this molecule is ", tMol.GetNumAtoms()
        for aAtom in tMol.GetAtoms():
            idxA  = aAtom.GetIdx()
            elemA = aAtom.GetSymbol()
            name  = aAtom.GetProp("Name") 
            charge = aAtom.GetFormalCharge()
            val    = aAtom.GetTotalValence()
            #print "For atom of index ", idxA
            #print "Its element symbol  ", elemA
            #print "Its name ", name
            #print "Its formal charge ", charge
            #print "Its total valence ", val

        Chem.SanitizeMol(tMol)
        Chem.Kekulize(tMol)
        # Make sure an atom in the molecule has the same charge in the input file. 
        # RDKit sometimes change it when initiating the molecule
        if tChargeList:
            for aAtom in tMol.GetAtoms():
                name  = aAtom.GetProp("Name")
                if tChargeList.has_key(name): 
                    aAtom.SetFormalCharge(int(tChargeList[name]))
                else:
                    aAtom.SetFormalCharge(0)
        """
        # Check
        print "==========================================="
        print "\nBefore setup formal charges in the molecule"
        allAtoms = tMol.GetAtoms()
        for aAtom in tMol.GetAtoms():
            idxA  = aAtom.GetIdx()
            elemA = aAtom.GetSymbol()
            name  = aAtom.GetProp("Name") 
            exHs    = aAtom.GetNumExplicitHs()
            imHs    = aAtom.GetNumImplicitHs()
            charge = aAtom.GetFormalCharge()
            val    = aAtom.GetTotalValence()
            print "\nFor atom of index ", idxA
            print "Its element symbol  ", elemA
            print "Its name ", name
            print "Its explicit Hs ", exHs
            print "Its mplicit Hs ", imHs
            print "Its formal charge ", charge
            print "Its total valence ", val
            #print "It is in the following bonds: "
            #aSetBonds = aAtom.GetBonds()
            #for aB in aSetBonds:
            #    print "Bond ", aB.GetIdx()
            #    print "Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name"))
            #    print "Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name"))
        print "==========================================="
        """

        if tPH[0] :
            self.setAllFormalChargeFuncGroupAtoms(tMol, tPH[1])
        else:
            self.setAllFormalChargeFuncGroupAtoms(tMol)
        
        """
        print "\nAfter setup formal charges in the molecule"
        for aAtom in tMol.GetAtoms():
            idxA   = aAtom.GetIdx()
            elemA  = aAtom.GetSymbol()
            name   = aAtom.GetProp("Name") 
            Hs     = aAtom.GetNumExplicitHs()
            imHs   = aAtom.GetNumImplicitHs()
            charge = aAtom.GetFormalCharge()
            val    = aAtom.GetTotalValence()
            print "\nFor atom of index ", idxA
            print "Its element symbol  ", elemA
            print "Its name ", name
            print "Its explicit Hs ", Hs
            print "Its mplicit Hs ", imHs
            print "Its formal charge ", charge
            print "Its total valence ", val
        print "==========================================="
        """
        tMol.UpdatePropertyCache()
        if self.useExistCoords:
            aMol = Chem.AddHs(tMol, explicitOnly=False, addCoords=True)
        else:
            aMol = Chem.AddHs(tMol)

        # Make SMILES before Hs are added 
        if tMol.HasProp('SmilesIn'):
            #print "Input SMILES : ", tMol.GetProp("SmilesIn")
            tMol.SetProp("SmilesOut", tMol.GetProp("SmilesIn"))
        #else:
            #isomericSmiles   = True
            #tKekuleSmiles    = True
            #rootedAtAtom     = -1
            #canonical        = True
            #allBondsExplicit = False
            #allHsExplicit    = False
            #tMol.SetProp("SmilesOut", Chem.MolToSmiles(aMol, isomericSmiles=False, kekuleSmiles=True, rootedAtAtom=-1, canonical=False, allBondsExplicit=False, allHsExplicit=False))
            #tMol.SetProp("SmilesOut", Chem.MolToSmiles(tMol, isomericSmiles=True, kekuleSmiles=False, canonical=False))
            #print "Output SMILES ", tMol.GetProp("SmilesOut") 
 
        allAtoms = aMol.GetAtoms()
        #print "After sanitize the molecule"
        #print "Number of atom is ", len(allAtoms)

        # Further: give names to those newly added H atoms
        self.setNamesForAtomsInMol(aMol, tChemCheck, tNameMap, 1)
        if self.reSetChirals:
            self.reAssignChirals(aMol) 
            """
            for aAtom in allAtoms:
                print "******************************"
                idxA  = aAtom.GetIdx()
                elemA = aAtom.GetSymbol()
                name  = aAtom.GetProp("Name") 
                charge = aAtom.GetFormalCharge()
                print "For atom of index ", idxA
                print "Its element symbol  ", elemA
                print "Its name ", name
                print "Its formal charge ", charge
        
                # Assign chiral centers 
                #if name == "C1":
                #    aAtom.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                #else:
                #    aAtom.SetChiralTag(rdchem.ChiralType.CHI_UNSPECIFIED)
                #aCT = aAtom.GetChiralTag()
                #print "Chiral center assigned ", aCT
            """
    
        for aAtom in allAtoms:
            aIdx = aAtom.GetIdx()
            tChemCheck.checkChiralCenters(aMol,aIdx)
            #print "Atom ", aAtom.GetProp("Name")
            #print "Is it a temporal chiral center ", aAtom.HasProp("TmpChiral") 

        self.setInitConformersOneMol(aMol)
       
        #rdmolfiles.MolToPDBFile(aMol, "Test_2.pdb")
        aSetTorsions = []
        self.assignTorsions(aMol, aSetTorsions)

        # print "Number of torsions in this molecule is ", len(aSetTorsions)
        tCurMols = [] 
        if tMode ==0: 
            self.molecules.append(aMol)
            tCurMols = self.molecules 
        elif tMode ==1:
            self.moleculesA.append(aMol) 
            tCurMols = self.moleculesA 
        elif tMode ==2:
            self.moleculesB.append(aMol) 
            tCurMols = self.moleculesB 
        # Check
        for aMol in tCurMols:
            conformers =  aMol.GetConformers()
            print "Number of all atoms in this molecule is ", aMol.GetNumAtoms()
            allAtoms = aMol.GetAtoms()
            #print "Number of atom is ", len(allAtoms)
            rdmolops.AssignStereochemistry(aMol, cleanIt=False, force=False, flagPossibleStereoCenters=True)
            for i in range (len(conformers)):
                aConf = conformers[i]
                if i==0:
                    print "----------------------------------------------------"
                    print "| In conformer %d                                  |"%(i+1)
                    print "----------------------------------------------------"
                    for aAtom in allAtoms:
                        print "******************************"
                        idxA  = aAtom.GetIdx()
                        elemA = aAtom.GetSymbol()
                        name  = aAtom.GetProp("Name") 
                        charge = aAtom.GetFormalCharge()
                        print "For atom of index ", idxA
                        print "Its element symbol  ", elemA
                        print "Its name ", name
                        print "Its formal charge ", charge
                        pos = aConf.GetAtomPosition(idxA)
                        if elemA.find("H") ==-1:
                            print "the coordinates are:"
                            print "X : %7.4f."%pos.x
                            print "Y : %7.4f"%pos.y
                            print "Z : %7.4f"%pos.z
                        aSetBonds = aAtom.GetBonds()
                        print "It is in the following bonds: "
                        for aB in aSetBonds:
                            print "Bond ", aB.GetIdx()
                            print "Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name"))
                            print "Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name"))
                    #print "******************************"

    def reAssignChirals(self, tMol):

        cleanIt=False
        force=False
        flagPossibleStereoCenters=True
        rdmolops.AssignStereochemistry(tMol, cleanIt, force, flagPossibleStereoCenters)
        allAtoms = tMol.GetAtoms()
        for aAtom in allAtoms:
            elem  = aAtom.GetSymbol()
            name  = aAtom.GetProp("Name") 
            if aAtom.HasProp('_CIPCode'):
                print "Chiral center ", name
                print "CIP rank %s : Stero code %s"%(aATom.GetProp("_CIPRank"), aAtom.GetProp("_CIPCode")) 
        sys.exit()   

    def mergeAtomNames(self):

        if len(self.moleculesA) == len(self.moleculesB):
            for iMolA in range(len(self.moleculesA)):
                allAtomsA =  self.moleculesA[iMolA].GetAtoms()
                allAtomsB =  self.moleculesB[iMolA].GetAtoms()
                aList =[]
                bList =[]
                aLen = len(allAtomsA)
                bLen = len(allAtomsB)

                for iAtmsA in range(aLen):
                    aCR = allAtomsA[iAtmsA].GetProp("_CIPRank")
                    aList.append([iAtmsA, int(aCR), allAtomsA[iAtmsA]])
                    if iAtmsA < bLen:
                        bCR = allAtomsB[iAtmsA].GetProp("_CIPRank")
                        bList.append([iAtmsA, int(bCR), allAtomsB[iAtmsA]])

                aList.sort(listCompDes)
                bList.sort(listCompDes)
                """
                for iA in range(aLen):
                    print "A: atom ",  aList[iA][2].GetProp("Name")
                    print "CIPRank ",  aList[iA][1]
                print "\n"
                for iB in range(bLen):
                    print "B: atom ",  bList[iB][2].GetProp("Name")
                    print "CIPRank ",  bList[iB][1]
                """
                for iP in range(aLen):
                    if iP < bLen:
                        # print "AList : atom ", aList[iP][2].GetProp("Name"), " is with CIPRank ", aList[iP][1]
                        aList[iP][2].SetProp("Name", bList[iP][2].GetProp("Name"))
                        #print "AList: new name ", aList[iP][2].GetProp("Name")
                        #print "BList: atom ", bList[iP][2].GetProp("Name"), " is with CIPRank ", bList[iP][1]

    def modifyMol(self, tMol, tAllAtoms, tAllBonds, tAllChirals, tDelAtomIdxs, tAtomsBondedToDel):

        allAtoms1 = tMol.GetAtoms()
        allBonds = tMol.GetBonds()       
        print "Element to be deleted in the mol is ", self.repSign
        for aBond in tMol.GetBonds():
            lBond = True
            atom1 = aBond.GetBeginAtom()
            idx1  = atom1.GetIdx()
            symb1 = atom1.GetSymbol()
            atom2 = aBond.GetEndAtom()
            idx2  = atom2.GetIdx()
            symb2 = atom2.GetSymbol()
            if symb1.find(self.repSign) !=-1:
                tAtomsBondedToDel.append(idx2)
                lBond = False
            elif symb2.find(self.repSign) !=-1:
                tAtomsBondedToDel.append(idx1)
                lBond = False
            if lBond:
                tAllBonds.append(aBond) 

        if len(tAtomsBondedToDel) !=0:
            print "Those atoms are linked the del atoms : "
            for aIdx in tAtomsBondedToDel:
                print tMol.GetAtomWithIdx(aIdx).GetProp("Name")
            
        for aAt in tMol.GetAtoms():
            if aAt.GetSymbol().find(self.repSign) ==-1:
                if aAt.GetIdx() in tAtomsBondedToDel:
                    tC_pre = aAt.GetFormalCharge()
                    tC     = tC_pre -1
                    aAt.SetFormalCharge(tC)
                tAllAtoms.append(aAt)
            else:
                tDelAtomIdxs.append(aAt.GetIdx())

    def MolToSimplifiedMmcif(self, tMol, tMmcifName, tChemCheck, tMonoName="UNL", tChiDes=None, tGroupName="non-polymer", tIdxConform=0):
        
        # A simplified mmcif file contains:
        # (1) Header section      
        # (2) Description of atoms in the molecule
        # (3) Description of bonds in the molecule
        # (4) Description of torsion angles in the molecules
        # (5) Description of chiral centers in the molecules 
    
        # This file is mainly used as an input file for Acedrg
        #print "Ligand ID ", tMonoName 
        #print "Group Name ", tGroupName 

        allAtoms         = []
        allAtoms1        = tMol.GetAtoms()
        delAtomIdxs      = []
        atomsBondedToDel = []

        allBonds   = []
        allChirals = []

        nAt  = 0
        nHAt = 0

        if self.reSetSmi:
            self.modifyMol(tMol, allAtoms, allBonds, allChirals, delAtomIdxs, atomsBondedToDel)
            nAt  = len(allAtoms)
            nHAt = tMol.GetNumHeavyAtoms() - len(delAtomIdxs)
        else:
            allAtoms = tMol.GetAtoms()
            allBonds = tMol.GetBonds()
            nAt  = len(allAtoms)
            nHAt = tMol.GetNumHeavyAtoms()

        #print "number of atoms with pseudo-atoms is ", tMol.GetNumAtoms()
        #print "number of atoms in the initial smiles is  ", len(allAtoms)
  
        try:
            aMmCif = open(tMmcifName, "w")
        except IOError:
            print tMmcifName, " Could not be opened for reading"
        else:

            # Header section 
           
            aMmCif.write("global_\n")
            aMmCif.write("_lib_name         ?\n")
            aMmCif.write("_lib_version      ?\n")
            aMmCif.write("_lib_update       ?\n")
            aMmCif.write("# ------------------------------------------------\n")
            aMmCif.write("#\n")
            
            # Monomer description 
            aMmCif.write("# ---   LIST OF MONOMERS ---\n")
            aMmCif.write("#\n")
            aMmCif.write("data_comp_list\n")
            aMmCif.write("loop_\n")
            aMmCif.write("_chem_comp.id\n")
            aMmCif.write("_chem_comp.three_letter_code\n")
            aMmCif.write("_chem_comp.name\n")
            aMmCif.write("_chem_comp.group\n")
            aMmCif.write("_chem_comp.number_atoms_all\n")
            aMmCif.write("_chem_comp.number_atoms_nh\n")
            aMmCif.write("_chem_comp.desc_level\n")
            aMmCif.write("%s       %s        %s        %s       %d     %d   %s\n" \
                         %(tMonoName, tMonoName, "\'.             \'",  tGroupName, nAt, nHAt, "."))
            aMmCif.write("# ------------------------------------------------------\n")
            aMmCif.write("# ------------------------------------------------------\n")
            aMmCif.write("#\n")
            aMmCif.write("# --- DESCRIPTION OF MONOMERS ---\n")
            aMmCif.write("#\n")
            aMmCif.write("data_comp_%s\n"%tMonoName)
            aMmCif.write("#\n")
        
            aConformer =  tMol.GetConformer(tIdxConform)
            rdmolops.AssignAtomChiralTagsFromStructure(tMol, confId=tIdxConform)

            # Atom section
            aMmCif.write("loop_\n")
            aMmCif.write("_chem_comp_atom.comp_id\n")
            aMmCif.write("_chem_comp_atom.atom_id\n")
            aMmCif.write("_chem_comp_atom.type_symbol\n")
            aMmCif.write("_chem_comp_atom.type_energy\n")
            aMmCif.write("_chem_comp_atom.charge\n")
            aMmCif.write("_chem_comp_atom.x\n")
            aMmCif.write("_chem_comp_atom.y\n")
            aMmCif.write("_chem_comp_atom.z\n")
            nTetraChi = 0 
            for aAtom in allAtoms:
                pos = aConformer.GetAtomPosition(aAtom.GetIdx())
                aChi = aAtom.GetChiralTag()
                if aChi != rdchem.ChiralType.CHI_UNSPECIFIED:
                    #print "Atom ", aAtom.GetProp("Name")
                    #print "Chiral center ? ", aChi
                    nTetraChi +=1
                elif aAtom.HasProp("TmpChiral") !=0:
                    nTetraChi +=1
                aMmCif.write("%s         %s      %s    %s     %3.2f   %5.4f    %5.4f     %5.4f\n" \
                             %(tMonoName, aAtom.GetProp("Name"), aAtom.GetSymbol(),  \
                               aAtom.GetSymbol(), float(aAtom.GetFormalCharge()), pos.x, pos.y, pos.z))
            # Bond section
            aMmCif.write("#\n")
            aMmCif.write("_chem_comp_bond.comp_id\n")
            aMmCif.write("_chem_comp_bond.atom_id_1\n")
            aMmCif.write("_chem_comp_bond.atom_id_2\n")
            aMmCif.write("_chem_comp_bond.type\n")
            aMmCif.write("_chem_comp_bond.aromatic\n")
            aMmCif.write("_chem_comp_bond.value_dist\n")
            aMmCif.write("_chem_comp_bond.value_dist_esd\n")
            #print "N Bonds ", len(allBonds)
            for aBond in allBonds:
                atom1 = aBond.GetBeginAtom()
                name1 = atom1.GetProp("Name")
                idx1  = atom1.GetIdx()
                symb1 = atom1.GetSymbol()
                atom2 = aBond.GetEndAtom()
                name2 = atom2.GetProp("Name")
                idx2  = atom2.GetIdx()
                symb2 = atom2.GetSymbol()

                bType = ""
                #print "Bond between atoms %s and %s is %s "%(name1, name2, aBond.GetBondType())
                if aBond.HasProp("SpecialBond"):
                    bType = aBond.GetProp("SpecialBond")
                else:
                    bType = aBond.GetBondType()

                isAro = "n"
                if aBond.GetIsAromatic():
                    isAr  = "y" 
                if self.reSetSmi:
                    if symb1.find(self.repSign)==-1 and symb2.find(self.repSign)==-1:
                        bLen  = rdMolTransforms.GetBondLength(aConformer, idx1, idx2)
                        dBlen = 0.20
                        aMmCif.write("%s       %s       %s       %s      %s     %5.4f     %5.4f\n" \
                                     %(tMonoName, name1, name2,  bType, \
                                       isAro, bLen, dBlen))
                else:
                    bLen  = rdMolTransforms.GetBondLength(aConformer, idx1, idx2)
                    dBlen = 0.20
                    aMmCif.write("%s       %s       %s       %s      %s     %5.4f     %5.4f\n" \
                                 %(tMonoName, name1, name2,  bType, \
                                   isAro, bLen, dBlen))

            # chiral center section
            lChiPre = False
            print " Number of chiral centers get from the comformer ", nTetraChi
            if tChiDes:
                print "number of chiral center predefined ", len(tChiDes)
                if len(tChiDes):
                    # The molecule contain chiral centers
                    aMmCif.write("#\n")
                    aMmCif.write("_chem_comp_chir.comp_id\n")
                    aMmCif.write("_chem_comp_chir.id\n")
                    aMmCif.write("_chem_comp_chir.atom_id_centre\n")
                    aMmCif.write("_chem_comp_chir.atom_id_1\n")
                    aMmCif.write("_chem_comp_chir.atom_id_2\n")
                    aMmCif.write("_chem_comp_chir.atom_id_3\n")
                    aMmCif.write("_chem_comp_chir.volume_sign\n")
                    for aChiral in tChiDes:
                        aMmCif.write(aChiral+"\n")
                    lChiPre=True
            elif nTetraChi !=0 and not lChiPre: 
                # The molecule contain chiral centers
                aMmCif.write("#\n")
                aMmCif.write("_chem_comp_chir.comp_id\n")
                aMmCif.write("_chem_comp_chir.id\n")
                aMmCif.write("_chem_comp_chir.atom_id_centre\n")
                aMmCif.write("_chem_comp_chir.atom_id_1\n")
                aMmCif.write("_chem_comp_chir.atom_id_2\n")
                aMmCif.write("_chem_comp_chir.atom_id_3\n")
                aMmCif.write("_chem_comp_chir.volume_sign\n")

                # Set all chiral atom according to format of mmCif 
                rdmolops.AssignStereochemistry(tMol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
                chiralIdx = 1
                for aAtom in allAtoms:      
                    aCT = aAtom.GetChiralTag()
                    aTmpCT = aAtom.HasProp("TmpChiral")
                    aSetAtomNames = []
                    if aCT != rdchem.ChiralType.CHI_UNSPECIFIED or aTmpCT !=0 :
                        aCTName = "chir_" + str(chiralIdx)
                        print "Atom %s "%aAtom.GetProp("Name")
                        print "Chiraltype ",  aCT
                        chiralIdx+=1
                        aIdx = aAtom.GetIdx()
                        outSign = ""
                        if aAtom.HasProp('_CIPCode'):
                            print "Chiraltype ",  aCT
                            print "Stereo type(CIPCode) ", aAtom.GetProp("_CIPCode") 
                            aCVol = EmbedLib.ComputeChiralVolume(tMol, aAtom.GetIdx(), tIdxConform)
                            aCVol = float(aCVol)
                            #if str(aAtom.GetProp("_CIPCode")).strip()=="S":
                                #outSign = "negative"
                            #    outSign = "positive"
                            #elif str(aAtom.GetProp("_CIPCode")).strip()=="R":
                            #    #outSign = "positive"
                            #    outSign = "negative"
                            """
                            if aCVol > 0.000001:
                                if aAtom.HasProp("pdb_stereo") and aAtom.GetProp("pdb_stereo")=="N":
                                    outSign = "both"
                                else: 
                                    outSign = "positive"
                            elif aCVol < -0.000001:
                                if aAtom.HasProp("pdb_stereo") and aAtom.GetProp("pdb_stereo")=="N":
                                    outSign = "both"
                                else: 
                                    outSign = "negative"
                            else:
                                outSign = "both"
                            """
                        #else:
                        #    outSign = "both"

                        #print "combined with original stereo sign, the chiral sign is ", outSign

                        aSetBonds = aAtom.GetBonds()
                        #print "It is in the following bonds: "
                        #for aB in aSetBonds:
                        #    print "Bond ", aB.GetIdx()
                        #    print "Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name"))
                        #    print "Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name"))
                        nNB = len(aAtom.GetNeighbors())
                        #print "Number of NB ", nNB
                        aSetName = []
       
                        for bIdx in range(nNB):
                            atmIdx1 = aSetBonds[bIdx].GetBeginAtomIdx()
                            atmIdx2 = aSetBonds[bIdx].GetEndAtomIdx() 
                            name1   = allAtoms[atmIdx1].GetProp("Name")  
                            name2   = allAtoms[atmIdx2].GetProp("Name")
                            symb1   = allAtoms[atmIdx1].GetSymbol()  
                            symb2   = allAtoms[atmIdx2].GetSymbol()
                            
                            if self.reSetSmi and len(self.repSign) !=0:
                                if not atmIdx1 in delAtomIdxs and not atmIdx2 in delAtomIdxs :
                                    if symb1.find(self.repSign)==-1 and symb2.find(self.repSign)==-1:
                                        if atmIdx1 == aIdx:                      
                                            aSetName.append([allAtoms1[atmIdx2], allAtoms1[atmIdx2].GetProp("_CIPRank")])
                                        elif atmIdx2 == aIdx:
                                            aSetName.append([allAtoms1[atmIdx1], allAtoms1[atmIdx1].GetProp("_CIPRank")])
                                        else:
                                            print "Bug! atom is not in bonds obtained by aAtom.GetBonds()"%(aAtom.GetProp("Name"))
                            else:
                                if atmIdx1 == aIdx:                      
                                    aSetName.append([allAtoms[atmIdx2], allAtoms[atmIdx2].GetProp("_CIPRank")])
                                elif atmIdx2 == aIdx:
                                    aSetName.append([allAtoms[atmIdx1], allAtoms[atmIdx1].GetProp("_CIPRank")])
                                else:
                                    print "Bug! atom is not in bonds obtained by aAtom.GetBonds()"%(aAtom.GetProp("Name"))

                        """
                        print "Atoms linked with the chiral center %s  are : "%aAtom.GetProp("Name")
                        print "Before sorting "
                        print "----------------------------------------------"
                        for aPair in aSetName:
                            print "Atom %s : CIPRank %s "%(aPair[0].GetProp("Name"), aPair[1])
                        print "----------------------------------------------"
                                         
                        print "After sorting "
                        """
                        aSetName.sort(listCompDes)
                        print "----------------------------------------------"
                        for aPair in aSetName:
                            print "Atom %s : CIPRank %s "%(aPair[0].GetProp("Name"), aPair[1])
                        print "----------------------------------------------"

                        aPass = self.checkUncertainChirals(aSetName)
                        if not aPass:
                            aMmCif.write("%s       %s       %s       %s      %s      %s      %s\n" \
                                         %(tMonoName, aCTName, aAtom.GetProp("Name"),  aSetName[0][0].GetProp("Name"), \
                                           aSetName[1][0].GetProp("Name"), aSetName[2][0].GetProp("Name"), "both"))
                            print "Chiral sign : both "
                        else:
                            posCen = aConformer.GetAtomPosition(aAtom.GetIdx())
                            coordsCen = [float(posCen.x), float(posCen.y), float(posCen.z)] 
                            pos1 = aConformer.GetAtomPosition(aSetName[0][0].GetIdx())
                            coords1 = [float(pos1.x), float(pos1.y), float(pos1.z)] 
                            pos2 = aConformer.GetAtomPosition(aSetName[1][0].GetIdx())
                            coords2 = [float(pos2.x), float(pos2.y), float(pos2.z)] 
                            pos3 = aConformer.GetAtomPosition(aSetName[2][0].GetIdx())
                            coords3 = [float(pos3.x), float(pos3.y), float(pos3.z)] 
                       
                            #print "Chiral volume by RDKit", aCVol
                            cenName = aAtom.GetProp("Name") 
                            tChemCheck.outChiSign[cenName] = tChemCheck.getChiralVolumeSign(posCen, pos1, pos2, pos3)
                            print "chiral sign : ", tChemCheck.outChiSign[cenName]
                            if len(tChemCheck.outChiSign[cenName]) !=0:                 
                                aMmCif.write("%s       %s       %s       %s      %s      %s      %s\n" \
                                             %(tMonoName, aCTName, aAtom.GetProp("Name"),  aSetName[0][0].GetProp("Name"), \
                                               aSetName[1][0].GetProp("Name"), aSetName[2][0].GetProp("Name"), tChemCheck.outChiSign[cenName]))
                            else:
                                print "Failed to calculate chiral volume for atom %s "%aAtom.GetProp("Name")
            aMmCif.close()
  
    def checkUncertainChirals(self, tNBAtoms):
        
        tPass = True
        tOneBondOs = 0
        tNBRanks   = {}
        for aPair in tNBAtoms:
            if aPair[0].GetSymbol().find("O") !=-1 and len(aPair[0].GetBonds())==1:
                tOneBondOs +=1 
            if not tNBRanks.has_key(aPair[1]):
                tNBRanks[aPair[1]] = []
            tNBRanks[aPair[1]].append(1)

        if tOneBondOs > 1 :
            tPass = False
        else :
            for aKey in tNBRanks.keys():
                if len(tNBRanks[aKey]) > 1:
                    tPass = False

        return tPass
        
    def setupFuncGroupTab(self, tFuncFileName):

        try:
            tFuncF = open(tFuncFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFuncFileName
            sys.exit()
        else:
            allFLs = tFuncF.readlines()
            tFuncF.close()

            for aL in allFLs:
                if aL.find("#")==-1 :
                    strGrp = aL.strip().split(":") 
                    if len(strGrp) ==3:
                        self.funcGroupTab[strGrp[0]] = [strGrp[1], float(strGrp[2])]

        # Generate molecules and  check
        #print "Number of functional groups is ", len(self.funcGroupTab)
        for aKey in sorted(self.funcGroupTab.iterkeys()):
            #print "SMARTS string for functional group %s is %s "%(aKey, self.funcGroupTab[aKey][0])                  
            #print "Its pKa value is  %5.3f "%(self.funcGroupTab[aKey][1]) 
            #print "Number of atoms in Smarts ",  Chem.MolFromSmarts(self.funcGroupTab[aKey][0]).GetNumAtoms()                
            #print "Number of atoms in Smiles ",  Chem.MolFromSmiles(self.funcGroupTab[aKey][0]).GetNumAtoms()                
            self.stdFuncGroupMols.append([aKey, Chem.MolFromSmarts(self.funcGroupTab[aKey][0]),  Chem.MolFromSmiles(self.funcGroupTab[aKey][0])])

    def getAllFuncGroupAtoms(self, tMol):
        self.funcGroups = {}
        for aTri in self.stdFuncGroupMols:
            atmGrpIdxs = tMol.GetSubstructMatches(aTri[1])
            if len(atmGrpIdxs)==0:
                atmGrpIdxs = tMol.GetSubstructMatches(aTri[2])
            if len(atmGrpIdxs):
                #print "search for the functional group %s "%aTri[0]
                #print "number of substruct ", atmGrpIdxs
                #print "Those atoms in this molecule are found in functional group %s: "%aTri[0]
                allAs = tMol.GetAtoms()
                #for i in range(len(atmGrpIdxs)):
                #    for aIdx in atmGrpIdxs[i]:
                #        print "atom ", allAs[aIdx].GetProp("Name")
                if not self.funcGroups.has_key(aTri[0]):
                    self.funcGroups[aTri[0]] = []
                for oneAtmIdxGrp in atmGrpIdxs:
                    #print aTri[0], "add func atoms"
                    self.funcGroups[aTri[0]].append(oneAtmIdxGrp)
                #print ""

    # The following methods are migrated from Acedrg c++ section

    def setAllFormalChargeFuncGroupAtoms(self, tMol, tPH=7.0):
      
        self.getAllFuncGroupAtoms(tMol)
        if len (self.funcGroups):
            for aFuncG in self.funcGroups.keys():
                # print "check ", aFuncG 
                if aFuncG.find("CARBOXY-AMINO-TERS") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeC_A_T(tMol, aFuncG,  self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-TER")!=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeC_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-ASP")!=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeC_ASP(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("AMINO-TER") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeA_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-AMINO-GLY") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeGLY(tMol, aFuncG,  self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("VAR-A-TER") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeV_A_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-LYS") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeNH_LYS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-HIS") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeNH_HIS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-PRO") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeNH_PRO(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-ARG") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeNH_ARG(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SH-CYS") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeSH_CYS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-TYR") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeBZ_TYR(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-NH") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeBZ_NH(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-NHA") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeBZ_NHA(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SO3") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeSO3(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SO4") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargeSO4(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO4") !=-1:
                    print "Doing ", aFuncG
                    self.setFormalChargePO4(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO3R") !=-1:
                    print "Doing ", aFuncG
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargePO3R(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO3N") !=-1:
                    print "Doing ", aFuncG
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargePO3N(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NITROMETHANE0") !=-1:
                    print "Doing ", aFuncG
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargeNITROMETHANE(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NITROMETHANE1") !=-1:
                    print "Doing ", aFuncG
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargeC_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)

    def setFormalChargeC_T(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                            and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0:     
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
     

    def setFormalChargeC_A_T(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0: 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
         
    def setFormalChargeC_ASP(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                #print "a set of atoms: "
                #print aSetIdxs
                for aIdx in aSetIdxs:
                    #print "Atom : ", tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                    #print "exH  : ", tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()
                    #print "imH  : ", tMol.GetAtomWithIdx(aIdx).GetNumImplicitHs()
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                            and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0:  
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
     
    def setFormalChargeA_T(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:    # implicit H included in the count
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            #print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                            #       tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                            # print "Its total valence is ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
         
    def setFormalChargeV_A_T(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:    # implicit H included in the count
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            #print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                            #       tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                            #print "Its total valence is ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()

    def setFormalChargeGLY(self, tMol, tFunG, tAtomIdxs, tPH):
   
        # This is for GLY only       
        nNonHAtoms = tMol.GetNumHeavyAtoms()
        for aSetIdxs in tAtomIdxs:
            if len(aSetIdxs) == nNonHAtoms:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                            and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0: 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
         
    def setFormalChargeNH_LYS(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N" \
                        and len(tMol.GetAtomWithIdx(aIdx).GetBonds())==1:
                        tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                        tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                        #print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                        #       tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                        #print "Its total valence is ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()

    def setFormalChargeNH_HIS(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH <= tPka  :
            lNSet = False
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        #print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        #print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        #print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        #print "num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() 
                        #print "num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors())
                        if (tMol.GetAtomWithIdx(aIdx).GetTotalValence() > 3 \
                            and len(tMol.GetAtomWithIdx(aIdx).GetBonds())==2) \
                            or \
                           (tMol.GetAtomWithIdx(aIdx).GetTotalValence() == 3 and \
                            len(tMol.GetAtomWithIdx(aIdx).GetBonds())==2 and \
                            tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()==0) :
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                        else:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(0)      # only one N has positive charge
                            #print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                            #   tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                        #print "Its total valence is ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
        elif tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        print "num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() 
                        print "num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors())
                        if tMol.GetAtomWithIdx(aIdx).GetFormalCharge() > 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0:    
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(0)       # de-protonation
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
         
    def setFormalChargeNH_PRO(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH <= tPka  :
            lNSet = False
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        print "num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() 
                        print "num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors())
                        if tMol.GetAtomWithIdx(aIdx).GetTotalValence() == 3 \
                            and tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print "atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                  tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :      
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()


    def setFormalChargeNH_ARG(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "Pka ", tPka
        #print "PH  ", tPH
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                nC =0
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        print "num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() 
                        print "num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors())
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() == 1     \
                           and tMol.GetAtomWithIdx(aIdx).GetTotalValence()==3 \
                           and len(tMol.GetAtomWithIdx(aIdx).GetBonds())==1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache() 
                            print "Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge()
                            break   
                             
    def setFormalChargeSH_CYS(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="S":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0:      # implicit H included in the count
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)

    def setFormalChargeSO3(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH > tPka :
            # Check if a -1 charge is already put on one of Os 
            nH =0
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0:      # implicit H included in the count
                            nH +=1
                if nH==2 :  # need to deprotonation for One of Os  
                    for aIdx in aSetIdxs:
                        if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                            if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                               and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :      
                                tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                                tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                                break

    def setFormalChargeSO4(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH > tPka :
            # Check if a -1 charge is already put on one of Os 
            nH =0
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0:      # implicit H included in the count
                            nH +=1
                if nH==2 :  # need to deprotonation for One of Os  
                    for aIdx in aSetIdxs:
                        if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                            if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                               and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :
                                tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                                tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)

    def setFormalChargeBZ_TYR(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        #print "PH ", tPH
        #print "tPka ", tPka
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)

    def setFormalChargeBZ_NH(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if len(tMol.GetAtomWithIdx(aIdx).GetBonds())==3: 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()

    def setFormalChargeBZ_NHA(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH < tPka :
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if len(tMol.GetAtomWithIdx(aIdx).GetBonds())==2 \
                           and tMol.GetAtomWithIdx(aIdx).GetTotalValence()==3: 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()

    def setFormalChargePO4(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        if tPH > tPka :
            # Put charge -1 to all 3 singly bonded Os 
            nH =0
            for aSetIdxs in tAtomIdxs:
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O" and\
                       tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                       and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :
                        # need to deprotonation for the O atom 
                        tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                        tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)

    def setFormalChargePO3R(self, tMol, tFunG, tAtomIdxs, tPH):
        tPka = self.funcGroupTab[tFunG][1]
        tIdxOs = []
        # Should exclude PO4
        for aSetIdxs in tAtomIdxs:
            nO=0
            #print ""
            for aIdx in aSetIdxs:
                #print "atom %s : "%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                    tIdxOs.append(aIdx)
                    nO+=1
            if nO==3:
                if tPH > tPka :
                    # Put charge -1 to 2 singly  bonded Os 
                    for aIdx in tIdxOs:
                        print "\nHere are details of O: "           
                        print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        print "num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()
                        print "num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors()) 
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :
                            # need to deprotonation for the O atom 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print "Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge()
                        #self.setOxyBondDeloc(tMol, aIdx)  
            tIdxOs = []
                                
    def setFormalChargePO3N(self, tMol, tFunG, tAtomIdxs, tPH):
        tPka = self.funcGroupTab[tFunG][1]
        tIdxOs = []
        # Should exclude PO4
        for aSetIdxs in tAtomIdxs:
            nO=0
            #print ""
            for aIdx in aSetIdxs:
                #print "atom %s : "%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                    tIdxOs.append(aIdx)
                    nO+=1
            if nO==3:
                if tPH > tPka :
                    # Put charge -1 to 3 singly  bonded Os 
                    for aIdx in tIdxOs:
                        # check and de-protonate O if needed
                        self.deProtonationO(tMol, aIdx)
                        
            tIdxOs = []
                                
    def setFormalChargeNITROMETHANE(self, tMol, tFunG, tAtomIdxs, tPH):
   
        tPka = self.funcGroupTab[tFunG][1]
        if tPH < tPka :
            lDeP = False
            for aSetIdxs in tAtomIdxs:
                idxN = -1
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        idxN = aIdx 
                for aIdx in aSetIdxs:
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        print "\nHere are details of O: "
                        print "atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                        print "total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence()
                        print "number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds())
                        if tMol.GetAtomWithIdx(aIdx).GetFormalCharge() !=1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                        tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                        print "Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge()  
                    elif tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if not lDeP:
                            # check and de-protonate O if needed
                            self.deProtonationO(tMol, aIdx)
                            lDeP = True
                        else: 
                            aBond = tMol.GetBondBetweenAtoms(idxN, aIdx)           
                            aBond.SetBondType(rdchem.BondType.DOUBLE)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            #self.setOxyBondDeloc(tMol, aIdx)
                           
    def deProtonationO(self, tMol, tIdx):
        print "\nHere are details of O: "           
        print "atom %s :"%tMol.GetAtomWithIdx(tIdx).GetProp("Name")
        print "total Val ", tMol.GetAtomWithIdx(tIdx).GetTotalValence()
        print "number of bonds ", len(tMol.GetAtomWithIdx(tIdx).GetBonds())
        print "num of H ", tMol.GetAtomWithIdx(tIdx).GetTotalNumHs()
        print "num of neighbor ", len(tMol.GetAtomWithIdx(tIdx).GetNeighbors()) 
        if tMol.GetAtomWithIdx(tIdx).GetTotalNumHs() != 0\
            and tMol.GetAtomWithIdx(tIdx).GetNumExplicitHs()==0 :
            tMol.GetAtomWithIdx(tIdx).SetFormalCharge(-1)
            tMol.GetAtomWithIdx(tIdx).SetNoImplicit(True)
            print "Its formal charge ", tMol.GetAtomWithIdx(tIdx).GetFormalCharge()
            #self.setOxyBondDeloc(tMol, tIdx)  
     
    def setOxyBondDeloc(self, tMol, tIdxO):
  
        oBonds = tMol.GetAtomWithIdx(tIdxO).GetBonds()
        for aBond in oBonds:
            tSet = []
            tSet.append(aBond.GetBeginAtom().GetSymbol())
            tSet.append(aBond.GetEndAtom().GetSymbol())
            if not ("O" in tSet and "H" in tSet):
                aBond.SetProp("SpecialBond", "DELOC")
         
    def assignTorsions(self, tMol, tAllTorsions):

        allAtoms    = tMol.GetAtoms()
        for aBond in tMol.GetBonds():
            atom1Idx   = aBond.GetBeginAtomIdx()
            atom1      = aBond.GetBeginAtom()
            atom2Idx   = aBond.GetEndAtomIdx()
            atom2      = aBond.GetEndAtom()
            atom1NBs   = atom1.GetNeighbors()
            atom1NBIdx = []
            atom2NBs   = atom2.GetNeighbors()
            atom2NBIdx = []
            for aNB in atom1NBs:
                idx1 = aNB.GetIdx()
                if idx1 != atom2Idx:
                    atom1NBIdx.append(idx1)
            for aNB in atom2NBs:
                idx2 = aNB.GetIdx()
                if idx2 != atom1Idx:
                    atom2NBIdx.append(idx2)
             
            if len(atom1NBIdx) > 0 and len(atom2NBIdx) > 0:
                for aEndOneAtomIdx in atom1NBIdx:
                    for aEndTwoAtomIdx in atom2NBIdx:
                        self.assignOneTorsion(aBond, allAtoms[aEndOneAtomIdx],\
                                               allAtoms[atom1Idx], allAtoms[atom2Idx],\
                                               allAtoms[aEndTwoAtomIdx], tAllTorsions)

    def assignOneTorsion(self, tBond, tAtom1, tAtom2, tAtom3, tAtom4, tAllTorsions):
        
        # Do not do actually torsion value calculations here. leave them 
        # to the corresponding "libmol" section
        # Here we just set-up other properties of a torsion
        aTor = {}
        aTor["compAtoms"] =[tAtom1.GetProp("Name"), tAtom2.GetProp("Name"), tAtom3.GetProp("Name"), tAtom4.GetProp("Name")]        
        if tBond.GetStereo()==Chem.rdchem.BondStereo.STEREOE\
           or tBond.GetStereo()==Chem.rdchem.BondStereo.STEREOZ:
            aTor["periodic"] = 1
            aTor["hybrid"]   = "const"
        else:
            if tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP2\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                aTor["periodic"] = 2
                aTor["hybrid"]   = "sp2_sp2"
            elif tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP3\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                aTor["periodic"] = 3
                aTor["hybrid"]   = "sp3_sp3"
            elif (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP3)\
               or (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP3\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP):
                aTor["periodic"] = 3
                aTor["hybrid"]   = "sp_sp3"
            elif (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP2\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP3)\
               or (tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP2\
               and tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP3):
                aTor["periodic"] = 6
                aTor["hybrid"]   = "sp2_sp3"
            elif (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP)\
               or (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP):
                aTor["periodic"] = 1
                aTor["hybrid"]   = "sp_sp"
            elif (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP2)\
               or (tAtom2.GetHybridization() == Chem.rdchem.HybridizationType.SP2\
               and tAtom3.GetHybridization() == Chem.rdchem.HybridizationType.SP):
                aTor["periodic"] = 1
                aTor["hybrid"]   = "sp_sp2"
            else:
                aTor["periodic"] = 1
                aTor["hybrid"]   = "others"

            tAllTorsions.append(aTor) 
    
    def removeWater(self, tMol):
       
        print "Number of atoms ", tMol.GetNumAtoms() 
        # Remove water molecules read from the input files (usually pdb format)
        aEdMol = rdchem.EditableMol(tMol)
        for aAt in tMol.GetAtoms():
            if aAt.GetSymbol() =="O":
                numBs = aAt.GetTotalDegree()
                numHs = aAt.GetTotalNumHs()
                if numBs == numHs:
                    print "Atom's element type ", aAt.GetSymbol()
                    print "Number of bonds it has ", numBs
                    print "Number of H  it connects ", numHs
                    aEdMol.RemoveAtom(aAt.GetIdx())
                    print "atom ", aAt.GetIdx(), " is removed "
        tMol = aEdMol.GetMol()
        print "Number of atoms after removing water", tMol.GetNumAtoms() 

class FileTransformer :

    def __init__(self):

        self.dataDescriptor  = {} 
        self.strDescriptors  = {} 
        self.strDescriptors["defProps"]    = ["loop_", "_pdbx_chem_comp_descriptor.comp_id", "_pdbx_chem_comp_descriptor.type", \
                               "_pdbx_chem_comp_descriptor.program", "_pdbx_chem_comp_descriptor.program_version",\
                               "_pdbx_chem_comp_descriptor.descriptor"]
        self.strDescriptors["defSmiles"]   = []

        self.cifGlobalLines                = []

        self.atoms           = []
        self.bonds           = []
        self.chirals         = []
        self.chiralPre       = []

        self.group           =""

        self.atomTypeMaps    = {}
        self.atomTypeMaps["atomTypes"]    = {}
        self.atomTypeMaps["connections"] = {}

        self.allBlockLs      = []

        self.bondTypeMmcifToMol = {}
        self.bondTypeMmcifToMol["SING"] = "1"
        self.bondTypeMmcifToMol["1"]    = "1"
        self.bondTypeMmcifToMol["DOUB"] = "2"
        self.bondTypeMmcifToMol["2"]    = "2"
        self.bondTypeMmcifToMol["TRIP"] = "3"
        self.bondTypeMmcifToMol["3"]    = "3"
        self.bondTypeMmcifToMol["AROM"] = "4"
        self.bondTypeMmcifToMol["AR"]   = "4"
        self.bondTypeMmcifToMol["AM"]   = "1"
        self.bondTypeMmcifToMol["DELO"] = "1"
        self.bondTypeMmcifToMol["ANY"]  = "8"
     
        self.bondTypeMolToMmcif = {}
        self.bondTypeMolToMmcif["1"] = "SING"
        self.bondTypeMolToMmcif["2"] = "DOUB"
        self.bondTypeMolToMmcif["3"] = "TRIP"
        self.bondTypeMolToMmcif["4"] = "AROM"
        self.bondTypeMolToMmcif["5"] = "DELO"
        self.bondTypeMolToMmcif["8"] = "ANY"

        self.cryst1      = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00                         " 

        self.ccp4DataDes = ["_chem_comp.id", "_chem_comp.three_letter_code", \
                            "_chem_comp.name", "_chem_comp.group",    \
                            "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh", \
                            "_chem_comp.desc_level"]
 
        self.rdkitSmiles  = {}

        self.ccp4MmCifDataMap = {}

        self.nameMapingCifMol = {}             # e.g. self.nameMapingCifMol[1]   = name 
                                               # where 1 : atom serial number   name : atom name
        self.nameMapingCifMol = {}            
        self.nameMapingCifMol["nonH"] = {}            
        self.nameMapingCifMol["H"]    = {}            
                                             

        self.nameMapingPDBMol = {}             # e.g. self.nameMapingPDBMol[1]   = name 
                                               # where 1 : atom serial number   name : atom name
        self.connInPDB        = {}

        self.delocBondList    = []             # recording bonds of the deloc type to put back late on 

        self.inputCharge      = {}             # input charges in a mmCif file

        self.PdbForMols       = {}

        self.hasCCP4Type      = False

    def mmCifReader(self, tFileName):
        """Read a detailed mmicif file to get basic information"""

        try:
            tFile = open(tFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFileName
            sys.exit()
        else:
            tAllLs = tFile.readlines()
            tFile.close()

            self.allBolcLs = []
            aBlockLs     = []
            all2ColLines = []
            filtedAllLines = []
            for aL in tAllLs:
                aL = aL.strip()
                if len(aL):
                    if aL[0].find("_") !=-1:
                        a2ColList = []
                        self.aLineToAlist(aL,a2ColList)
                        if len(a2ColList)==2:
                            all2ColLines.append(a2ColList)
                        else:
                            filtedAllLines.append(aL)
                    else:
                        filtedAllLines.append(aL)

            for aL in filtedAllLines :
                #if aL.find("data_") != -1 and aL.find("#") ==-1:
                #    if len(aBlockLs) !=0:
                #        self.allBlockLs.append(aBlockLs)
                #    aBlockLs     = []
                if aL.find("loop_") !=-1 and aL.find("#") ==-1:
                    if len(aBlockLs) !=0:
                        self.allBlockLs.append(aBlockLs)
                    aBlockLs     = []
                else:
                    if len(aL.strip()) >0:                   
                        if aL.strip()[0].find("#") == -1:
                            aBlockLs.append(aL)  

            # Last block
            if len(aBlockLs):
                self.allBlockLs.append(aBlockLs)
       
            print "Number of Blocks ", len(self.allBlockLs)
            i = 1
            for aBlc in self.allBlockLs:
                print "-------------"
                print "Block ", i
                i = i + 1
                for aL in aBlc:
                    print aL  
    

            if len(self.allBlockLs):
                for aBlk in self.allBlockLs:
                    self.parseOneMmCifBlk(aBlk)

            if len(all2ColLines) !=0:
                self.parserAll2Cols(all2ColLines)

            self.TmpChemCheck()
           
            # check
            """
            idKey = "_chem_comp_atom.atom_id"
            for aAtom in self.atoms:
                if aAtom.has_key(idKey):
                    print "==============================="
                    print "For atom ", aAtom[idKey], " : "
                    print "-------------------------------"  
                    for aKey in aAtom.keys():
                        print "label : ", aKey, " Value : ", aAtom[aKey]
                    print "==============================="
            """

    def TmpChemCheck(self):

        """ 
        Temporal function
        Exclude some of cases where RDKit does not accept certain valence values, e.g. those in ligands with B
        Cancel this function when RDKit make the corresponding changes 
        """ 
 
        if len(self.bonds) > 0 :
            atomVals = {}
            
            for aBond in self.bonds:
                if aBond.has_key("_chem_comp_bond.atom_id_1") and\
                   aBond.has_key("_chem_comp_bond.atom_id_2") and\
                   aBond.has_key("_chem_comp_bond.value_order"):
                    aOrd = 0
                    if aBond["_chem_comp_bond.value_order"].upper().find("SING") !=-1:
                        aOrd = 1
                    elif aBond["_chem_comp_bond.value_order"].upper().find("DOUB") !=-1:
                        aOrd = 2
                    elif aBond["_chem_comp_bond.value_order"].upper().find("TRIP") !=-1:
                        aOrd = 3
                    if not atomVals.has_key(aBond["_chem_comp_bond.atom_id_1"]):
                        atomVals[aBond["_chem_comp_bond.atom_id_1"]] = 0
                    if not atomVals.has_key(aBond["_chem_comp_bond.atom_id_2"]):
                        atomVals[aBond["_chem_comp_bond.atom_id_2"]] = 0
                    atomVals[aBond["_chem_comp_bond.atom_id_1"]] += aOrd
                    #print aBond["_chem_comp_bond.atom_id_1"], " : ", aOrd
                    atomVals[aBond["_chem_comp_bond.atom_id_2"]] += aOrd
                    #print aBond["_chem_comp_bond.atom_id_2"], " : ", aOrd
             
            # find IDs of B, Br etc
            speAtomIds = []
            for aAtom in self.atoms:
                if aAtom.has_key("_chem_comp_atom.atom_id") and\
                   aAtom.has_key("_chem_comp_atom.type_symbol"):
                    if aAtom["_chem_comp_atom.type_symbol"].strip() =="B"\
                       or aAtom["_chem_comp_atom.type_symbol"].strip()=="BR":
                        speAtomIds.append(aAtom["_chem_comp_atom.atom_id"])
                        if aAtom.has_key("_chem_comp_atom.charge"):
                            if aAtom["_chem_comp_atom.charge"].find(".") !=-1:
                                aAtom["_chem_comp_atom.charge"] = aAtom["_chem_comp_atom.charge"].strip().split(".")[0]
                            elif aAtom["_chem_comp_atom.charge"].find("?") !=-1:
                                aAtom["_chem_comp_atom.charge"] =0
                            aCh = int(aAtom["_chem_comp_atom.charge"])
                            if aCh !=0 and atomVals.has_key(aAtom["_chem_comp_atom.atom_id"]):
                                atomVals[aAtom["_chem_comp_atom.atom_id"]] +=aCh
                                
            if len(speAtomIds):
                for aId in speAtomIds:
                    if atomVals.has_key(aId):
                        if atomVals[aId] > 3: 
                            print "Acedrg stops because RDKit does not accept the following: "
                            print "%s  has total valence of %d "%(aId, atomVals[aId])
                            sys.exit()        
        
    def parserAll2Cols(self, t2ColLines):

        iC =0
      
        atomProp = {}
        bondProp = {}
        for aPair in t2ColLines:
            if aPair[0].find("_chem_comp.") !=-1:
               self.dataDescriptor[iC]=[aPair[0], aPair[1]]    
               iC+=1
            elif aPair[0].find("_chem_comp_atom.") !=-1:
                atomProp[aPair[0]] = aPair[1]
            elif aPair[0].find("_chem_comp_bond.") !=-1:
                bondProp[aPair[0]] = aPair[1]
  
        if len(atomProp) !=0:
            self.atoms.append(atomProp)
        if len(bondProp) !=0:
            self.bonds.append(bondProp)
           
    def parseOneMmCifBlk(self, tBlk):

        if len(tBlk):
            for aL in tBlk:
                if aL.find("_chem_comp.") !=-1:
                    self.getDataDescriptor(tBlk)
                    break
                elif aL.find("_chem_comp_atom.") !=-1:
                    self.getProp(tBlk, "atom")
                    break
                elif aL.find("_chem_comp_bond.") !=-1:
                    self.getProp(tBlk, "bond")
                    break
                elif aL.find("_chem_comp_chir.") !=-1:
                    self.getProp(tBlk, "chiral")
                    break
                elif aL.find("_pdbx_chem_comp_descriptor") !=-1:
                    self.getProp(tBlk, "strDescriptor")
                    break


    def getDataDescriptor(self, tBlk):

        nAll = 0;
        nC   = 0
        l2  = False

        for aL in tBlk:
            if aL.find("_chem_comp.") !=-1:
                nAll+=1
                strGrp = aL.strip().split()
                if len(strGrp) == 2:
                    l2 = True
                    break
                elif len(strGrp) == 1:
                    nC += 1
            else:
                break
 
        if nC==nAll:
           l2 = False

        if l2: # 2 col format
            iC =0
            for aL in tBlk:
                if aL.find("\"") !=-1:
                    strGrp = aL.strip().split("\"")
                    if len(strGrp) >= 2: 
                        self.dataDescriptor[iC]=[strGrp[0], "\""+ strGrp[1] + "\""]
                        iC +=1
                elif aL.find("\'") !=-1 :
                    strGrp = aL.strip().split("\'")
                    if len(strGrp) >= 2: 
                        self.dataDescriptor[iC]=[strGrp[0], "\'"+ strGrp[1] + "\'"]
                        iC +=1
                else:
                    strGrp = aL.strip().split()
                    if len(strGrp) == 2:
                        self.dataDescriptor[iC]=[strGrp[0], strGrp[1]]
                        iC +=1
                                           
        else:   # multiple col format          
            colIdx = []
            for aL in tBlk:
                if aL.find("\'") !=-1:
                    strGrp1 = aL.strip().split("\'")
                    strGrp  = []
                    if len(strGrp1) == 3: 
			strGrp10 = strGrp1[0].strip().split()
                        strGrp12 = strGrp1[2].strip().split()
                        for aS in strGrp10:
                            strGrp.append(aS)
                        strGrp.append("\'" + strGrp1[1] + "\'")
                        for aS in strGrp12:
                            strGrp.append(aS)
                        if len(strGrp)==len(colIdx):
                            for i in range(len(strGrp)):
                                self.dataDescriptor[i]=[colIdx[i], strGrp[i]]
                elif aL.find("\"") !=-1:
                    strGrp1 = aL.strip().split("\"")
                    strGrp  = []
                    if len(strGrp1) == 3: 
			strGrp10 = strGrp1[0].strip().split()
                        strGrp12 = strGrp1[2].strip().split()
                        for aS in strGrp10:
                            strGrp.append(aS)
                        strGrp.append("\"" + strGrp1[1] + "\"")
                        for aS in strGrp12:
                            strGrp.append(aS)
                        if len(strGrp)==len(colIdx):
                            for i in range(len(strGrp)):
                                self.dataDescriptor[i]=[colIdx[i], strGrp[i]]

                else:
                    strGrp = aL.strip().split()
                    if len(strGrp) == 1:
                        colIdx.append(strGrp[0])
        """
        # Check
        print "Two  colum format :"
        for i in sorted(self.dataDescriptor.iterkeys()):
            print "%s%s"%(self.dataDescriptor[i][0].ljust(60), self.dataDescriptor[i][1].ljust(40))
        print "\n"

        print "Multple  colum format :"
        for i in sorted(self.dataDescriptor.iterkeys()):
            print self.dataDescriptor[i][0]
        aSt = ""
        for i in sorted(self.dataDescriptor.iterkeys()):
            aSt+=(self.dataDescriptor[i][1].strip() + "\t")
        print aSt
        print "\n" 
        """

    def getCCP4DataDescritor(self, tMol, tChemCheck, tMonomRoot="UNL"):

        # Get a CCP4 monomer lib data descriptor
        s1 =tMonomRoot
        s2 =tMonomRoot
        s3 =".         "
        s4 ="non-polymer"
        s5 = str(tMol.GetNumAtoms())
        s6 = str(tMol.GetNumHeavyAtoms())
        s7 ="."

        for aKey in sorted(self.dataDescriptor.iterkeys()):
            if self.dataDescriptor[aKey][0].find("_chem_comp.id") !=-1:
                s1 = self.dataDescriptor[aKey][1]  
                s2 = self.dataDescriptor[aKey][1]  
            elif self.dataDescriptor[aKey][0].find("_chem_comp.name") !=-1:
                s3 = self.dataDescriptor[aKey][1]
            elif self.dataDescriptor[aKey][0].find("_chem_comp.group") !=-1:
                s4 = self.dataDescriptor[aKey][1]  
            elif self.dataDescriptor[aKey][0].find("_chem_comp.type") !=-1:
                if tMol.GetProp("ResidueName") in tChemCheck.aminoAcids:
                    s4 = "L-PEPTIDE"
                elif self.dataDescriptor[aKey][1].upper().find("DNA ") !=-1 :
                    s4 = "DNA"
                elif self.dataDescriptor[aKey][1].upper().find("RNA ") !=-1 :
                    s4 = "RNA"
                else:
                    s4 = self.dataDescriptor[aKey][1]  
            
        aLine = "%s%s%s%s%s%s%s"%(s1.ljust(len(s1)+5), s2.ljust(len(s2)+5), \
                                    s3.ljust(len(s3)+5), s4.ljust(len(s4)+5), \
                                    s5.ljust(len(s5)+5), s6.ljust(len(s6)+5), \
                                    s7.ljust(len(s7)+5))
        self.ccp4DataDes.append(aLine)

        # set a SMILES string here from the molecule
        if tMol.HasProp("SmilesOut"):
            aSmi             = tMol.GetProp("SmilesOut")
            aLine = "%s%s%s%s\"%s\"\n"%(tMonomRoot.ljust(10), "SMILES".ljust(10), "RDKit".ljust(12), "1.00".ljust(6), \
                                        aSmi)                                    
            self.strDescriptors["defSmiles"].append(aLine)

        # Check
        #for aL in self.ccp4DataDes:
        #    print aL

    def getProp(self, tBlk, tProp):

        # multiple col format          
        colIdx = []
        iStrD = 0
        for aL in tBlk:
            aL = aL.strip()
            strGrp = []
            self.aLineToAlist(aL, strGrp)
            #strGrp = aL.strip().split()
            if len(strGrp) == 1 and aL[0].find("_") !=-1 and strGrp[0].find(";")==-1:
                colIdx.append(strGrp[0])
            else:
                #if aL.find("\"") !=-1:
                #    if tProp =="atom" or tProp =="bond" or tProp =="chiral" or tProp =="strDescriptor" : 
                #        strGrp1 = aL.strip().split("\"")
                #        strGrp  = []
                #        for aStrC in strGrp1:
                #            strGrp2 = aStrC.strip().split()
                #            if len(strGrp2) !=0:
                #                for aS in strGrp2:
                #                    strGrp.append(aS)

                if len(strGrp)==len(colIdx):
                    if tProp =="atom" or tProp =="bond" or tProp =="chiral" : 
                        aProp = {}
                        for i in range(len(strGrp)):
                            aProp[colIdx[i]] = strGrp[i]
                        if tProp =="atom":
                            self.atoms.append(aProp)
                        elif tProp =="bond":
                            self.bonds.append(aProp)
                        elif tProp =="chiral":
                            self.chirals.append(aProp)
                    elif tProp =="strDescriptor":
                        if not self.strDescriptors.has_key("props"):
                            self.strDescriptors["props"] = []
                            for iProp in colIdx: 
                                self.strDescriptors["props"].append(iProp)
                        if not self.strDescriptors.has_key("entries"):
                            self.strDescriptors["entries"] = []
                        self.strDescriptors["entries"].append(aL)
                else :
                     if not self.strDescriptors.has_key("entries"):
                         self.strDescriptors["entries"] = []
                     self.strDescriptors["entries"].append(aL)

    
        if tProp =="atom":
            tAtoms = []
            tHAtoms = []
            for aAtom in self.atoms:
                if aAtom["_chem_comp_atom.type_symbol"] !="H":
                    tAtoms.append(aAtom)
                else:
                    tHAtoms.append(aAtom)

                if aAtom.has_key("_chem_comp_atom.atom_id"):
                    tCharge =0.0
                    if aAtom.has_key("_chem_comp_atom.charge"):
                        print "aAtom['_chem_comp_atom.charge'] ", aAtom["_chem_comp_atom.charge"] 
                        if aAtom["_chem_comp_atom.charge"].find("?") ==-1:
                            tCharge = float(aAtom["_chem_comp_atom.charge"])
                    elif aAtom.has_key("_chem_comp_atom.partial_charge"):
                        tCharge = float(aAtom["_chem_comp_atom.partial_charge"])
                        aAtom["_chem_comp_atom.charge"] = aAtom["_chem_comp_atom.partial_charge"]
                    if math.fabs(float(tCharge)) >=0.001:
                        self.inputCharge[aAtom["_chem_comp_atom.atom_id"]] = tCharge
                if aAtom.has_key("_chem_comp_atom.type_energy"):
                    self.hasCCP4Type = True
            self.atoms = []
            for aAtom in tAtoms:
                self.atoms.append(aAtom)
            for aAtom in tHAtoms:
                self.atoms.append(aAtom)
        elif tProp =="chiral":
            for aChiral in self.chirals:
                aStr = ""
                for i in range(len(colIdx)):
                    aStr +="%s"%aChiral[colIdx[i]].ljust(8)
                self.chiralPre.append(aStr)

        if len(colIdx):
            if tProp =="atom":
                for aAtom in self.atoms:
                    for i in range(len(colIdx)):
                        if colIdx[i].find("chem_comp_atom.atom_id") != -1:
                            print "Prop %s is %s "%(colIdx[i], aAtom[colIdx[i]].ljust(8))
        # Check
        if len(self.inputCharge) !=0:
            print "The following atoms have charges "
            for aName in self.inputCharge.keys():
                print "Name : ", aName, " charge : ", self.inputCharge[aName]

        if len(colIdx):
            for i in range(len(colIdx)):
                print colIdx[i]
            if tProp =="atom":
                for aAtom in self.atoms:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aAtom[colIdx[i]].ljust(8)
                    print aStr
            elif tProp =="bond":
                for aBond in self.bonds:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aBond[colIdx[i]].ljust(8)
                    print aStr
            elif tProp =="chiral":
                for aChiral in self.chirals:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aChiral[colIdx[i]].ljust(8)
                    print aStr

        if tProp =="strDescriptor":
            if len(self.strDescriptors["props"]):
                for aProp in self.strDescriptors["props"]:
                    print aProp
                for aEn in self.strDescriptors["entries"]:
                    print aEn.strip()

    def getCCP4MmCifMap(self):

        self.ccp4DataDis = ["_chem_comp.id", "_chem_comp.three_letter_code", \
                            "_chem_comp.name", "_chem_comp.group",    \
                            "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh", \
                            "_chem_comp.desc_level"]
        self.ccp4MmCifDataMap = {}

        if len(self.dataDescriptor):
            nFind =0
            for aIdx in self.dataDescriptor.keys():
                if self.dataDescriptor[aIdx][0] in self.ccp4DataDis:
                    self.ccp4MmCifDataMap[self.dataDescriptor[aIdx][0]]=aIdx
                    nFind +=1    
                                
    def DelocBondConvertor(self):

        doneList = []
        delocMap = {}

        aBarg = ""

        for i in range(len(self.bonds)):
            aB2 = ""
            if self.bonds[i].has_key("_chem_comp_bond.value_order"):
                aBarg = "_chem_comp_bond.value_order"
                aB2 = self.bonds[i]["_chem_comp_bond.value_order"]
            elif self.bonds[i].has_key("_chem_comp_bond.type"):             
                aBarg = "_chem_comp_bond.type"
                aB2 = self.bonds[i]["_chem_comp_bond.type"]
            if len(aB2) !=0:
                if aB2.lower().find("delo") !=-1:
                    a1 = self.bonds[i]["_chem_comp_bond.atom_id_1"]
                    a2 = self.bonds[i]["_chem_comp_bond.atom_id_2"]
                    if not delocMap.has_key(a1) :
                        delocMap[a1] = []
                    if not delocMap.has_key(a2) :
                        delocMap[a2] = []
                    delocMap[a1].append(i)
                    delocMap[a2].append(i)
                    self.delocBondList.append([a1.strip(), a2.strip()])

        if len(delocMap) !=0 and aBarg !="":
            for aKey in sorted(delocMap.iterkeys()):
                if len(delocMap[aKey])==2 or len(delocMap[aKey])==3:
                    if not delocMap[aKey][0] in doneList and not delocMap[aKey][1] in doneList:
                        self.bonds[delocMap[aKey][0]][aBarg] = "1"
                        self.bonds[delocMap[aKey][1]][aBarg] = "2"
                        doneList.append(delocMap[aKey][0])
                        doneList.append(delocMap[aKey][1])
        """
        if len(self.delocBondList) !=0:
            print "Deloc bonds exist between the following atom pairs "
            for aPair in self.delocBondList:
                print "Atom %s and %s "%(aPair[0], aPair[1]) 
        """

    def mergeAtomNames(self, tFileName, tMol):

        try:
            tFile = open(tFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFileName
            sys.exit()
        else:
            nameMap = {}
            for aL in tFile.readlines():
                aL = aL.strip()
                if len(aL) !=0:
                    if aL.find("#")==-1:
                        strGrp= aL.split()
                        if len(strGrp)==2:
                            nameMap[strGrp[0].strip()] = strGrp[1].strip()
            tFile.close()

            for aAtom in tMol.GetAtoms():
                oldName = aAtom.GetProp("Name") 
                if nameMap.has_key(oldName):
                    aAtom.SetProp("Name", nameMap[oldName]) 
                    print "Atom: old name %s : new name %s"%(oldName, aAtom.GetProp("Name"))

    def addAtomOrigChiralSign(self, tMol):

        nameMap = {}
   
        for i in range(len(self.atoms)):
            if self.atoms[i].has_key("_chem_comp_atom.atom_id"):
                nameMap[self.atoms[i]["_chem_comp_atom.atom_id"].strip()] = i
 
        for aAtom in tMol.GetAtoms():
            atomId = aAtom.GetProp("Name").strip()
            if nameMap.has_key(atomId):
                if self.atoms[nameMap[atomId]].has_key("_chem_comp_atom.pdbx_stereo_config"):
                    tSign = self.atoms[nameMap[atomId]]["_chem_comp_atom.pdbx_stereo_config"].strip() 
                    aAtom.SetProp("pdb_stereo", tSign)
                    #print "atom %s now has pdb_stereo %s"%(aAtom.GetProp("Name"), aAtom.GetProp("pdb_stereo"))

    def setAtomDicts(self, tFileName, tAtomTypeMaps):

        try:
            tFile = open(tFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFileName
            sys.exit()
        else:
            lA = False
            lC = False
            for aL in tFile.readlines():
                aL = aL.strip()
                if len(aL) !=0:
                    if aL.find("ATOMS:") !=-1:
                        lA = True
                        lC = False
                    elif aL.find("CONNECTIONS:") !=-1:
                        lC = True
                        lA = False
                    else:
                        strGrp = aL.strip().split()
                        if len(strGrp)==4: 
                            if lA :
                                idxA = int(strGrp[0])
                                if not tAtomTypeMaps["atomTypes"].has_key(strGrp[3]):
                                    tAtomTypeMaps["atomTypes"][strGrp[3]]=[]
                                tAtomTypeMaps["atomTypes"][strGrp[3]].append(idxA)
                                if not tAtomTypeMaps["atoms"].has_key(idxA):
                                    tAtomTypeMaps["atoms"][idxA]={}
                                tAtomTypeMaps["atoms"][idxA]["elem"]  = strGrp[1]
                                tAtomTypeMaps["atoms"][idxA]["name"]  = strGrp[2]
                                tAtomTypeMaps["atoms"][idxA]["class"] = strGrp[3]
                        if len(strGrp)==2: 
                            if lC :
                                idx0 = int(strGrp[0])
                                idx1 = int(strGrp[1])
                                if not tAtomTypeMaps["connections"].has_key(idx0):
                                    tAtomTypeMaps["connections"][idx0]=[]
                                if not tAtomTypeMaps["connections"].has_key(idx1):
                                    tAtomTypeMaps["connections"][idx1]=[]
                                tAtomTypeMaps["connections"][idx0].append(idx1)
                                tAtomTypeMaps["connections"][idx1].append(idx0)
          
            tFile.close() 

            # Check 
            print "Here are the details for atoms: "
            nAtms =0
            for aT in  sorted(tAtomTypeMaps["atomTypes"].iterkeys()):
                tAtomTypeMaps["atomTypes"][aT].sort()
                print " %d atoms has atom-type %s "%(len(tAtomTypeMaps["atomTypes"][aT]), aT)
                nAtms += len(tAtomTypeMaps["atomTypes"][aT])
                for aA in tAtomTypeMaps["atomTypes"][aT]:
                    print "atom %s  of serial number %d "%(tAtomTypeMaps["atoms"][aA]["name"], aA)                                               
                    print "which bonds the following atoms :"
                    print tAtomTypeMaps["connections"][aA]
                    for aNA in tAtomTypeMaps["connections"][aA]:
                        print "atom %s "%tAtomTypeMaps["atoms"][aNA]["name"]
            print "Total number of atoms is ", nAtms
    
    def AtomDictMapping(self, atomSet1, atomSet2, tMol):
 
        allIdxs  = []
        nonHIdxs = []
        for aIdx in sorted(atomSet1["atoms"].iterkeys()):
            allIdxs.append(aIdx)
            if atomSet1["atoms"][aIdx]["elem"].find("H")==-1:
                nonHIdxs.append(aIdx)
                
        print "total number of atoms ", len(allIdxs)
        print "total number of nonH atoms ", len(nonHIdxs)

        doneList        = []
        matchedAtoms    = {}
        remainClasses   = []
       
        # Match unique classses
        for aClass in atomSet1["atomTypes"].keys():
            if len(atomSet1["atomTypes"][aClass]) == 1 and aClass in atomSet2["atomTypes"]:
                if len(atomSet2["atomTypes"][aClass]) == 1 :
                    aIdx1 = atomSet1["atomTypes"][aClass][0]
                    aIdx2 = atomSet2["atomTypes"][aClass][0]
                    matchedAtoms[aIdx1] = aIdx2
                    doneList.append(aIdx2)
                else:
                    print "atom type %s appears in sys 1 and sys 2 different times "%aClass 
                    sys.exit()
            else:
                remainClasses.append(aClass)
  
        print "number of matched atoms ", len(doneList)
        print "Matched atoms:  "
        for aIdx in sorted(matchedAtoms.iterkeys()):
            bIdx = matchedAtoms[aIdx]
            print atomSet1["atoms"][aIdx]["class"]
            print "Atom %s to Atom %s "%(atomSet1["atoms"][aIdx]["name"],\
                                         atomSet2["atoms"][bIdx]["name"]) 

        print "Number of atoms to be matched ", len(atomSet1["atoms"])-len(doneList)
        print "Un-matched atoms:  "
        for aIdx in allIdxs:
            if not aIdx in doneList:
                print "Atom : %s of class %s "%(atomSet1["atoms"][aIdx]["name"], atomSet1["atoms"][aIdx]["class"])

        # Match symmetrical ones (Non-H) 
        # 1. singly bonded classes
        for aC in remainClasses:
            nAC1 = len(atomSet1["atomTypes"][aC])
            if aC in atomSet2["atomTypes"].keys():
                nAC2 = len(atomSet2["atomTypes"][aC])
            else:
                print "Check %s in sys1 but not in sys2 "%aC
                sys.exit()

            if nAC1 > 1 and nAC1==nAC2:
                aIdx0 = atomSet1["atomTypes"][aC][0]
                if len(atomSet1["connections"][aIdx0])==1 and atomSet1["atoms"][aIdx0]["elem"].find("H")==-1:
                    for i in range(nAC1):
                        aIdx1 = atomSet1["atomTypes"][aC][i]
                        aIdx2 = atomSet2["atomTypes"][aC][i]
                        matchedAtoms[aIdx1] = aIdx2
                        print "Atom %s to Atom %s "%(atomSet1["atoms"][aIdx1]["name"],\
                                                     atomSet2["atoms"][aIdx2]["name"])
                        doneList.append(aIdx2)

        # 2. Linked symmetrical classes
        # Match the rest using connections 
        tmpDoneList1 = matchedAtoms.keys()
        tmpDoneList2 = []
        for aIdx in tmpDoneList1:
            for aNA in atomSet1["connection"][aIdx]:
                if not aNA1 in doneList and not aNA1 in tmpDoneList \
                   and atomSet1["atoms"][aNA1]["elem"].find("H")==-1 :
                    group1 = []
                    group2 = []
                    aClassNA = atomSet1["atoms"][aNA1]["class"]
                    if not atomSet2["atomTypes"].has_key(aClassNA):
                        print "Discrepany in sys 1 and 2 at %s"%aClassNA
                    else:      
                        for aNNA1 in atomSet1["connection"][aNA1]:
                            if aNNA1 != aIdx and atomSet1["atoms"][aNNA1]["elem"].find("H")==-1 :
                                group1.append([aNA1, aNNA1])
                            #if atomSet2["atomTypes"].has_key(aClassNA) and atomSet2["atomTypes"].has_key(aClassNA):
                      
                         
                                             
       
    def MmCifToMolFile(self, tInFileName, tOutMolName, tMode=0):

        if tMode==0: 
            self.mmCifReader(tInFileName)

        #print "Num of atoms ", len(self.atoms)
        #print "Num of bonds ", len(self.bonds)
        if not len(self.atoms) or not len(self.bonds):
            print "No atoms and/or bonds from the input file, check !"
            sys.exit()
            
        try:
            tOutFile = open(tOutMolName, "w")
        except IOError:
            print "%s can not be open for reading "%tOutMolName
            sys.exit()
        else:
            nId   = -1 
            nName = -1
            for aKey in sorted(self.dataDescriptor.iterkeys()):
                if self.dataDescriptor[aKey][0].find("_chem_comp.id") !=-1:
                    nId = aKey
                if self.dataDescriptor[aKey][0].find("_chem_comp.name") !=-1:
                    nName = aKey
          
            # Header section 
            
            if nId !=-1:
                
                tOutFile.write(self.dataDescriptor[nId][1]+ "\n")
            else:
                tOutFile.write("UNL\n")
            if nName !=-1:
                tOutFile.write(self.dataDescriptor[nName][1]+ "\n")
            else:
                tOutFile.write("UNL\n")
            tOutFile.write("\n")
   
            # The Counts Line
            nA = str(len(self.atoms))
            nB = str(len(self.bonds))
            nC = ""
            if len(self.chirals):
                nC = "1"
            else:
                nC ="0"
            tOutFile.write("%s%s%s%s%s%s%s%s%s%s%s%s\n"%(nA.rjust(3), nB.rjust(3), "0".rjust(3), \
                                                         "0".rjust(3), nC.rjust(3), "0".rjust(3), \
                                                         " ".rjust(3), " ".rjust(3), " ".rjust(3), \
                                                         " ".rjust(3), "999".rjust(3), "V2000".rjust(6)))
            # Atom block

            # Re-arrange atoms: non-H atoms first, followed by H atoms
            print "Befor re-arrange, number of atoms is : ", len(self.atoms) 
            tNonHAtoms = []
            tHAtoms    = []
            for aAtom in self.atoms:
                id = ""
                if aAtom.has_key("_chem_comp_atom.type_symbol"):
                    tId = aAtom["_chem_comp_atom.type_symbol"].strip()
                    if len(tId) ==1:
                        id = tId.strip().upper()
                    elif len(tId) > 1:
                        id = tId[0].upper() + tId[1:].strip().lower()
                else:
                    print "Input file bug: no type_symbol for atoms!"
                    sys.exit()
                aAtom["type_symbol_in_mol"] = id
                if id != "H":
                    tNonHAtoms.append(aAtom)
                else :
                   tHAtoms.append(aAtom)
                
            self.atoms = []
            nAtm =0
            for aAtom in tNonHAtoms:
                self.atoms.append(aAtom)
                self.nameMapingCifMol["nonH"][nAtm] = aAtom["_chem_comp_atom.atom_id"]
                #print "NameMap ", aAtom["_chem_comp_atom.atom_id"]
                nAtm +=1
            for aAtom in tHAtoms:
                self.atoms.append(aAtom)

            # Set up atom seq match for bond section
            mapIdNum = {}
            nAtm =1
            for aAtom in self.atoms:
                if aAtom.has_key("_chem_comp_atom.atom_id"):
                    mapIdNum[aAtom["_chem_comp_atom.atom_id"]] = nAtm
                    nAtm +=1
                    #print "atom %s serial number in mol is %d "%(aAtom["_chem_comp_atom.atom_id"], mapIdNum[aAtom["_chem_comp_atom.atom_id"]])
                else:
                    print "Input file bug: no atom_id for atoms!"
                    sys.exit()
    
            # Now write out atom section in the Mol file
            hhh ="0"
            bbb ="0"
            vvv ="0"
            HHH ="0"
            rrr ="0"
            iii ="0"
            mmm ="0"
            nnn ="0"
            eee ="0"
            # Ignore all of the coordinates
            chargeAtomList = []
            idxAtom = 1
            for aAtom in self.atoms:
                #print "Atom ", aAtom["_chem_comp_atom.atom_id"], " with serial number ", idxAtom
                if aAtom.has_key("_chem_comp_atom.pdbx_model_Cartn_x_ideal") and\
                   aAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"].find("?")==-1:
                    x = aAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"]
                elif aAtom.has_key("_chem_comp_atom.model_Cartn_x"):
                    x = aAtom["_chem_comp_atom.model_Cartn_x"]
                elif aAtom.has_key("_chem_comp_atom.x"):
                    x = aAtom["_chem_comp_atom.x"]
                else:
                    x = "0.0000"

                if aAtom.has_key("_chem_comp_atom.pdbx_model_Cartn_y_ideal") and\
                   aAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"] ==-1:
                    y = aAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"]
                elif aAtom.has_key("_chem_comp_atom.model_Cartn_y"):
                    y = aAtom["_chem_comp_atom.model_Cartn_y"]
                elif aAtom.has_key("_chem_comp_atom.y"):
                    y = aAtom["_chem_comp_atom.y"]
                else:
                    y = "0.0000"

                if aAtom.has_key("_chem_comp_atom.pdbx_model_Cartn_z_ideal") and\
                   aAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"]==-1:
                    z = aAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"]
                elif aAtom.has_key("_chem_comp_atom.model_Cartn_z"):
                    z = aAtom["_chem_comp_atom.model_Cartn_z"]
                elif aAtom.has_key("_chem_comp_atom.z"):
                    z = aAtom["_chem_comp_atom.z"]
                else:
                    z = "0.0000"

                id = aAtom["type_symbol_in_mol"]
                    
                md =" 0 "
               
                # formal charges
 
                chargeMap = {}
                chargeMap[0]   = 0
                chargeMap[1]   = 3
                chargeMap[2]   = 2
                chargeMap[3]   = 1
                chargeMap[-1]  = 5
                chargeMap[-2]  = 6
                chargeMap[-3]  = 7

                chargeMap["doublet radical"]  = 4
                
                charge = " 0 "
                if aAtom.has_key("_chem_comp_atom.charge"):
                    #print "The charge is ", aAtom["_chem_comp_atom.charge"]
                    #print "Is that number digit ", aAtom["_chem_comp_atom.charge"].isdigit() 
                    if aAtom["_chem_comp_atom.charge"].find(".") !=-1:
                        aAtom["_chem_comp_atom.charge"] = aAtom["_chem_comp_atom.charge"].strip().split(".")[0]
                    nCharge =0
                    if aAtom["_chem_comp_atom.charge"].find("?") ==-1:
                        nCharge = int(aAtom["_chem_comp_atom.charge"])
                    #print "Atom ", aAtom["_chem_comp_atom.atom_id"], " has charge  ", aAtom["_chem_comp_atom.charge"]
                    #print " converted to charge symbol ", chargeMap[nCharge] 
                    if nCharge in chargeMap.keys():
                        charge  = " %d "%chargeMap[nCharge]
                    if nCharge !=0:
                        chargeAtomList.append([idxAtom, nCharge])       

                cha ="0"
                if aAtom.has_key("_chem_comp_atom.pdbx_stereo_config"):
                    if aAtom["_chem_comp_atom.pdbx_stereo_config"].find("S") !=-1:
                        cha = "1"
                    elif aAtom["_chem_comp_atom.pdbx_stereo_config"].find("R") !=-1:
                        cha = "2"

                tOutFile.write("%s%s%s %s%s%s%s%s%s%s%s%s%s%s%s%s\n"%(x.rjust(10), y.rjust(10), z.rjust(10), \
                                                                     id.ljust(3), md.rjust(2), charge.rjust(3), \
                                                                     cha.rjust(3), hhh.rjust(3), bbb.rjust(3), \
                                                                     vvv.rjust(3), HHH.rjust(3), rrr.rjust(3), \
                                                                     iii.rjust(3), mmm.rjust(3), nnn.rjust(3), \
                                                                     eee.rjust(3)))
                idxAtom +=1

            # Bond block
           
            self.DelocBondConvertor()
 
            sss ="0"
            xxx ="0"
            rrr ="0"
            ccc ="0"
            for aBond in self.bonds:
                id1 = aBond["_chem_comp_bond.atom_id_1"]
                id2 = aBond["_chem_comp_bond.atom_id_2"]
                #print "id1 %s "%id1 
                #print "id2 %s "%id2
                n1 = mapIdNum[id1]
                n2 = mapIdNum[id2]
                #print "id1 %s n1 %d "%(id1, n1)
                #print "id2 %s n2 %d "%(id2, n2)
                a1 = str(n1)
                a2 = str(n2)
                b4 = ""
                bt = ""
                if aBond.has_key("_chem_comp_bond.value_order") or  aBond.has_key("_chem_comp_bond.type"):         # mmcif in PDB or ccp4 monolib 
                    if aBond.has_key("_chem_comp_bond.value_order"):
                        aB4 = aBond["_chem_comp_bond.value_order"]
                    elif aBond.has_key("_chem_comp_bond.type"):             
                        aB4 = aBond["_chem_comp_bond.type"]
                    if len(aB4) >=4: 
                        b4 = aB4[:4].upper()
                    else:
                        b4 = aB4.upper()
                    # print b4
                    bt = self.bondTypeMmcifToMol[b4]
                else:
                    print "Input file bug: no bond type(order) for bonds!"        
                    sys.exit()
                aBL = "%s%s%s%s%s%s%s\n"%(a1.rjust(3), a2.rjust(3), bt.rjust(3), \
                       sss.rjust(3), xxx.rjust(3), rrr.rjust(3), ccc.rjust(3))    
                #print "The bond between %s of serial number %s and %s of serial number %s is : "\
                #      %(aBond["_chem_comp_bond.atom_id_1"], a1, aBond["_chem_comp_bond.atom_id_2"], a2)
                #print aBL
                tOutFile.write("%s%s%s%s%s%s%s\n"%(a1.rjust(3), a2.rjust(3), bt.rjust(3), \
                               sss.rjust(3), xxx.rjust(3), rrr.rjust(3), ccc.rjust(3)))

                elem1 = self.atoms[n1-1]["type_symbol_in_mol"]
                id1   = self.atoms[n1-1]["_chem_comp_atom.atom_id"]
                elem2 = self.atoms[n2-1]["type_symbol_in_mol"]
                id2   = self.atoms[n2-1]["_chem_comp_atom.atom_id"]
                #print "id1 %s elem1 %s "%(id1, elem1)
                #print "id2 %s elem2 %s "%(id2, elem2)
                if elem1 == "H":
                    if not self.nameMapingCifMol["H"].has_key(id2):
                        self.nameMapingCifMol["H"][id2] = [] 
                    self.nameMapingCifMol["H"][id2].append(id1)
                    #print "H atom %s is bonding to %s "%(id1, id2)
                if elem2 == "H":
                    if not self.nameMapingCifMol["H"].has_key(id1):
                        self.nameMapingCifMol["H"][id1] = [] 
                    self.nameMapingCifMol["H"][id1].append(id2)
                    #print "H atom %s is bonding to %s "%(id2, id1)
                   
            if len(chargeAtomList) != 0:
                aL = "M CHG  %d "%len(chargeAtomList)
                for aPair in chargeAtomList:
                    aL += " %d  %d "%(aPair[0], aPair[1])
                tOutFile.write(aL + "\n")          
   
            tOutFile.write("M  END\n")
          
    # Mol files related 
    def CheckElemSymbolsInMolFile(self, tInFileName, tOutFileName):

        if os.path.isfile(tInFileName):
            aMolF = open(tInFileName, "r")
            allMolLs = aMolF.readlines()
            aMolF.close()
            
            aMolSecs = {}
            iMol = 0
            aMolSecs[iMol] = {}
            aMolSecs[iMol]["Sec1"] = []
            aMolSecs[iMol]["Sec2"] = []
            aMolSecs[iMol]["Sec3"] = []
            nOneMolLines = 0
            nAtoms = 0
            for aL in allMolLs:
                if aL.find("$$$$") != -1:
                    aMolSecs[iMol]["Sec3"].append(aL)
                    iMol = iMol + 1
                    aMolSecs[iMol] = {}
                    aMolSecs[iMol]["Sec1"] = []
                    aMolSecs[iMol]["Sec2"] = []
                    aMolSecs[iMol]["Sec3"] = []
                    nOneMolLines = 0
                    nAtoms       = 0
                elif nOneMolLines < 4:
                    aMolSecs[iMol]["Sec1"].append(aL)
                    if nOneMolLines ==3:
                        tN = aL[:3].strip()
                        if tN.isdigit():
                            nAtoms = int(tN)
                        else:
                            print "Format error is input MOL/SDF file. The count line is : "
                            print aL 
                            sys.exit()
                elif nOneMolLines >= 4:
                     if nAtoms > 0:
                         if nOneMolLines >=4 and nOneMolLines < nAtoms:
                             aSym = aL[30:34].strip()
                             if len(aSym) !=2:  
                                 aMolSecs[iMol]["Sec2"].append(aL)
                             else:
                                 aSym2 = aSym[0] + aSym[1].lower()
                                 tL = aL[:30] + " " + aSym2 + " " + aL[34:]
                                 aMolSecs[iMol]["Sec2"].append(tL)
                         else: 
                             aMolSecs[iMol]["Sec3"].append(aL)
                     else:    
                            print "Format error is input MOL/SDF file. No count line found  "
                            sys.exit()
                nOneMolLines +=1

            if len(aMolSecs[0]["Sec2"]) > 0:
                # at least there is a molecule in the file.
                outF = open(tOutFileName, "w")
                for aMol in sorted(aMolSecs.iterkeys()):
                    for aL in aMolSecs[aMol]["Sec1"]:
                        outF.write(aL)
                    for aL in aMolSecs[aMol]["Sec2"]:
                        outF.write(aL)
                    for aL in aMolSecs[aMol]["Sec3"]:
                        outF.write(aL)
                outF.close()

    def MolToPDBFile(self, tOutFileName, idxMol, tMol, tDataDiscriptor=None, tMonoRoot="UNL", idxConf=0, tDelSign=""):

        try:
            tPDB = open(tOutFileName, "w")
        except IOError:
            print "%s (PDB format) can not be open for writing "%tOutFileName
            sys.exit()
        else:
            if not self.PdbForMols.has_key(idxMol):
                self.PdbForMols[idxMol] = []
            self.PdbForMols[idxMol].append(tOutFileName)
      
            self.PdbForMols    = {}
            # Head section 
            tClassification="UNL"
            tDate =time.strftime("%d/%m/%Y")
            tIdCode = str(tMonoRoot) 
           
            if tDataDiscriptor:
                for aIdx in tDataDiscriptor.keys():
                    if tDataDiscriptor[aIdx][0].find("_chem_comp.group") !=-1  \
                       or tDataDiscriptor[aIdx][0].find("_chem_comp.name") !=-1 \
                       or tDataDiscriptor[aIdx][0].find("_chem_comp.type") !=-1 :
                        tClassification= tDataDiscriptor[aIdx][1].strip()
                    elif tDataDiscriptor[aIdx][0].find("_chem_comp.id") !=-1:
                        tIdCode = tDataDiscriptor[aIdx][1].strip()
     
            tPDB.write("%s%s%s%s\n"%("HEADER".ljust(10), tClassification.ljust(40), tDate.ljust(9), tIdCode.rjust(7)))
            tPDB.write("CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1\n")
            # Atom section 
            # print idxConf
            aConf  =  tMol.GetConformer(idxConf)
            atoms      =  tMol.GetAtoms()
            for aAtom in atoms:
                idxA  = aAtom.GetIdx() 
                sIdxA  = str(aAtom.GetIdx() + 1)
                name  = aAtom.GetProp("Name").strip() 
                if len(name)<=2:
                    name = name + " "
                elif len(name)>3:
                    if name[0].find("\"")!=-1\
                       or name[0].find("\'")!=-1:
                        name = name[1:-1]
                rName = tIdCode[:3]
                seqNum = "1"
                empt1   = "   "
                pos = aConf.GetAtomPosition(idxA)
                posX ="%8.3f"%pos.x
                posY ="%8.3f"%pos.y
                posZ ="%8.3f"%pos.z
                ocp    = "%6.2f"%1.00
                b      = "%6.2f"%0.0
                empt2   = "         "
                elem   = aAtom.GetSymbol().strip()
                charge = aAtom.GetFormalCharge()
                sCharge = ""
                if charge > 0:
                   sCharge = str(charge) + "+"
                elif charge < 0:
                   sCharge = str(abs(charge)) + "-"
                if len(tDelSign)==0 and elem.find("H")==-1:
                    tPDB.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%("HETATM".ljust(6), sIdxA.rjust(5), name.rjust(5),  \
                                                            rName.rjust(4), seqNum.rjust(6), empt1.rjust(4), \
                                                            posX.rjust(8), posY.rjust(8),posZ.rjust(8),  \
                                                            ocp.rjust(6), b.rjust(6), empt2.rjust(10), elem.rjust(2), sCharge.rjust(2)))
                elif elem.find(tDelSign)==-1 and elem.find("H")==-1:
                    tPDB.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%("HETATM".ljust(6), sIdxA.rjust(5), name.rjust(5),  \
                                                            rName.rjust(4), seqNum.rjust(6), empt1.rjust(4), \
                                                            posX.rjust(8), posY.rjust(8),posZ.rjust(8),  \
                                                            ocp.rjust(6), b.rjust(6), empt2.rjust(10), elem.rjust(2), sCharge.rjust(2)))
       
                    
            # End section
            tPDB.write("END\n") 
            tPDB.close()

    def getAtomNamesInPDB(self, tInFileName):

        try:
            tPDB = open(tInFileName, "r")
        except IOError:
            print "%s (PDB format) can not be open for reading "%tInFileName
            sys.exit()
        else:
            allLines =tPDB.readlines()
            tPDB.close()
            
            for aL in allLines:
                head = aL[0:6].strip()
                if head.find("HETATM") !=-1 or head.find("ATOM") !=-1:
                    name = aL[12:16].strip() 
                    sNum = aL[6:11].strip()
                    if sNum.isdigit():
                        seriNum = int(sNum)-1
                        self.nameMapingPDBMol[seriNum] = name
                    else:
                        print "Errors in the input file %s format (PDB) "%tInFileName
                        print "Serial number of atom is ", sNum
                        print "in line %s "%aL
                        sys.exit()   
                if head.find("CONECT") !=-1 :
                    sNum = aL[6:11].strip()
                    if sNum.isdigit():
                        seriNum = int(sNum)-1
                        if not self.connInPDB.has_key(seriNum):
                            self.connInPDB[seriNum] = []
                    else:
                        print "Errors in the input file %s format (PDB) "%tInFileName
                        print "Serial number of atom is ", sNum
                        print "in line %s "%aL
                        sys.exit() 
                    if len(aL) > 16:
                        cNum = int(aL[11:16].strip())-1
                        self.connInPDB[seriNum].append(cNum)  
                    if len(aL) > 21:
                        cNum = int(aL[16:21].strip())-1
                        self.connInPDB[seriNum].append(cNum)  
                    if len(aL) > 26:
                        cNum = int(aL[21:26].strip())-1 
                        self.connInPDB[seriNum].append(cNum)  
                    if len(aL) >= 31:
                        cNum = int(aL[26:].strip())-1
                        self.connInPDB[seriNum].append(cNum)  
               
            # Check
            print "Atom names in PDB file %s are recorded here "%tInFileName 
            for aKey in sorted(self.nameMapingPDBMol.iterkeys()):
                print "Atom : ", aKey
                print "Name : ", self.nameMapingPDBMol[aKey]
                print "Its connection are : "
                for aConn in self.connInPDB[aKey]:
                    print "serial number %d and name %s "%(aConn, self.nameMapingPDBMol[aConn])
          
    def checkAndAddCryst1InPDB(self, tInFileName, tOutFileName):

        lCheck = False
        try:
            inPDB = open(tInFileName, "r")
        except IOError:
            print "%s (PDB format) can not be open for reading "%tInFileName
            sys.exit()
        else:
            allLs = inPDB.readlines()
            inPDB.close()
            for aL in allLs:
                if len(aL) > 6 :
                    if aL[0:6].find("CRYST1") != -1:
                        lCheck = True
                        break
       
            if lCheck:
                shutil.copyfile(tInFileName, tOutFileName)
            else:
                try: 
                    outPDB = open(tOutFileName, "w")
                except IOError:
                    print "%s (PDB format) can not be open for writing "%tOutFileName
                    sys.exit()
                else:
                    outPDB.write(self.cryst1 + "\n")
                    for aL in allLs:
                        outPDB.write(aL)
                    outPDB.close()

    def getResNameFromPDB(self, tInFileName):

        aResName = "UNL"
        try:
            inPDB = open(tInFileName, "r")
        except IOError:
            print "%s (PDB format) can not be open for reading "%tInFileName
            sys.exit()
        else:
            allLs = inPDB.readlines()
            inPDB.close()
            for aL in allLs:
                if len(aL) > 21 and (aL.find("ATOM") !=-1 or aL.find("HETATM") !=-1):
                    tStr = aL[17:21].strip()
                    if len (tStr) :
                        aResName = tStr
                        break

        return aResName
   
    def aLineToAlist(self, tL, tList):
   
        aTS = ""
        l0  = False
        l1  = False
        l2  = False   

        for aS in tL:
            aS = aS.strip()
            if len(aS) !=0:
                if aS.find("\"") !=-1:
                    if l2:
                        aTS = "\"" + aTS + "\""
                        tList.append(aTS.strip())
                        aTS = ""
                        l2  = False
                    elif l1:
                        aTS = aTS + aS
                    else:
                        l2 = True
                        l1 = False
                elif aS.find("\'") !=-1:
                    if l1:
                        aTS = "\'" + aTS + "\'"
                        tList.append(aTS)
                        aTS = ""
                        l1  = False
                    elif l2:
                        aTS = aTS + aS
                    else:
                        l1 = True
                        l2  = False
                else:
                    if not l0 and not l1 and not l2:
                        l0 = True
                    aTS = aTS + aS
            else:
                if l2 or l1:
                    aTS = aTS + aS
                elif l0:
                    if len(aTS):
                        tList.append(aTS.strip())
                        aTS = ""
                        l0  = False
                else:
                    l0  = True
        if aTS != "":
            tList.append(aTS.strip())        

class ChemCheck():

    def __init__( self):

        self.organicSec = ["AT", "At", "at", "B", "b", "BR", "Br", "br", "C", "c", "CL", "Cl", "cl", 
                   "F", "f", "H", "h", "I", "i", "N","n",  "O", "o", "P", "p", "S", "s", "SE", "Se", "se"]
        self.atomFileType = {}
        self.atomFileType["mmCif"]  = [11, 111, 16, 161, 51]
        self.atomFileType["simils"] = [12, 121, 52]
        self.atomFileType["mdl"]    = [13, 131, 14, 141, 33, 53, 54]
        self.atomFileType["mol2"]   = [15, 151, 55]

        self.pTable = rdchem.GetPeriodicTable()

        self.aminoAcids = ["ALA", "ARG", "ASN", "ASP", "CYS", \
                           "GLN", "GLU", "GLY", "HIS", "ILE", \
                           "LEU", "LYS", "MET", "PHE", "PRO", \
                           "SER", "THR", "TRP", "TYR", "VAL"]

        self.torsions   = []
        self.outChiSign  ={}
        self.tmpChiralSign = {}


    def isOrganicMol(self, tMol):

        pass 

    def isOrganic1(self, tElemList):

        organicOnly = True

        for aElem in tElemList:
            if not aElem.upper() in self.organicSec:
                organicOnly = False
                break;

        return organicOnly 

    def isOrganic(self, tInFileName, tTypeIdx):

        organicOnly = True

        allAtomElems = []
        self.getAtomElems(tInFileName, tTypeIdx, allAtomElems)
        
        for aAtomElem in allAtomElems:
            if not aAtomElem in self.organicSec:
                organicOnly = False
                break;
        if not organicOnly :
            print "The input ligands/molecules contains metal or other heavier atoms "
            print "Acedrg currently deals with ligands/molecules with following elements only "
            print "C, N, O, S, P, B, F, Cl, Br, I, H"
        if len(allAtomElems) ==0:
            print "The input ligands/molecules contains metal or other heavier atoms "
            print "Acedrg currently deals with ligands/molecules with following elements only "
            print "C, N, O, S, P, B, F, Cl, Br, I, H"
            #print "Can not get the element symbols of atoms. Check input file format"
            organicOnly = False

        return organicOnly 

    def isOrganicMol(self, tMol):

        pass 


    def addHs(self, tMol):

        # print "Number of atom before is ", tMol.GetNumAtoms()
        aMolT=Chem.MolFromSmiles(rdmolfiles.MolToSmiles(tMol))
        print rdmolfiles.MolToSmiles(tMol)
        aMol = Chem.AddHs(aMolT)
        rdmolfiles.MolToMolFile(tMol, "1.mol")
        print "Number of atom now is ", aMol.GetNumAtoms()
        tDoneAtms = []
        for i in range(len(tMol.GetAtoms())):
            elem = tMol.GetAtomWithIdx(i).GetSymbol()
            if elem != "H":   
                aBO = tMol.GetAtomWithIdx(i).GetTotalValence()
                aDV = self.pTable.GetDefaultValence(elem)
                nAddHs = aDV - aBO
                tMol.GetAtomWithIdx(i).SetFormalCharge(nAddHs) 

        aMolT = Chem.AddHs(tMol)
        print "Number of atom after is ", aMolT.GetNumAtoms()
        
    def getAtomElems(self, tInFileName, tTypeIdx, tAtomElems):
        if tTypeIdx in self.atomFileType["mmCif"]: 
            self.getAtomElemsFromMmcif(tInFileName, tAtomElems)
        elif tTypeIdx in self.atomFileType["simils"]:
                self.getAtomElemsFromSmiles(tInFileName, tAtomElems)
        elif tTypeIdx in self.atomFileType["mdl"]: 
            self.getAtomElemsFromMdl(tInFileName, tAtomElems)
        elif tTypeIdx in self.atomFileType["mol2"]: 
            self.getAtomElemsFromMol2(tInFileName, tAtomElems)
        else: 
            print "Could not recognize format for the input file. It should be one of cif, smiles, sdf/mol, mol2"
            sys.exit(1) 
        aLine = ""
        for aEl in tAtomElems:
             aLine += (aEl.strip()+"\t")
        if len(aLine.strip()) > 0:
            print "The system contains atoms of the following elements" 
            print aLine

    def getAtomElemsFromMmcif(self, tInFileName, tAtomElems):

        try :
            inFile = open(tInFileName, "r")
        except IOError :
             print "%s can not be opened for reading "%tInFileName
             sys.exit(1)
        else:
             allLines = inFile.readlines()
             inFile.close()
             lAtom = False
             iCol = 0
             colDict = {}
             for aL in allLines:
                if aL.find("loop") != -1 and lAtom:
                    lAtom = False
                    iCol =0
                    break
                elif lAtom and aL.find("_chem_comp_atom.") ==-1 :
                    strGrp = aL.strip().split()
                    if colDict.has_key("type_symbol") and len(strGrp) == len(colDict) \
                       and colDict["type_symbol"] < len(strGrp):
                        aAtomElem = strGrp[colDict["type_symbol"]]
                        if not aAtomElem in tAtomElems:
                            tAtomElems.append(aAtomElem)
                elif aL.find("_chem_comp_atom.") !=-1 :
                    strGrp = aL.strip().split(".")
                    if not lAtom:
                        lAtom = True
                    if len(strGrp) ==2:
                        colDict[strGrp[1]] = iCol
                        iCol +=1
                    else:
                        pass 
                        #print "Definition error in the input cif file %s"%tInFileName
                        #print "The entry is ", aL
                    
    def getAtomElemsFromMdl(self, tInFileName, tAtomElems):

        try :
            inFile = open(tInFileName, "r")
        except IOError :
            print "%s can not be opened for reading "%tInFileName
            sys.exit(1)
        else:
            allLines = inFile.readlines()
            inFile.close()
            if len(allLines) > 3 :
                nOneMolLines = 0
                nAtoms       = 0
                for aL in allLines:
                    if aL.find("$$$$") != -1:
                        nOneMolLines = 0 
                        nAtoms       = 0
                    elif nOneMolLines == 3:
                        tN = aL[:3].strip()
                        if tN.isdigit():
                            nAtoms = int(tN)
                        else:
                            print "Format error is input MOL/SDF file. The count line is : "
                            print aL 
                            sys.exit()
                    elif nOneMolLines > 3 and nOneMolLines < nAtoms:
                        aElem = aL[30:34].strip()
                        if not aElem in tAtomElems:
                            tAtomElems.append(aL[30:34].strip())  
                    nOneMolLines += 1

    def getAtomElemsFromSmiles(self, tFileName, tAtomElems):

        if os.path.isfile(tFileName):
            # SMILES string in a file
            try:
                fSmi = open(tFileName, "r")
            except IOError:
                print  "% can not be open for reading "%tFileName
                sys.exit()
            else:
                aSmiStr = fSmi.read()
                tSmiStr = aSmiStr.strip().split()
                if len(tSmiStr) >0 :
                    aSmiStr = tSmiStr[0].strip()
                else:
                    print "String format error"
                    sys.exit()    
                fSmi.close()
        else:
            # SMILES string in from a commandline  
            aSmiStr = tFileName.strip()
        if len(aSmiStr):
            checkMol=Chem.MolFromSmiles(aSmiStr)
  
            if checkMol:
                for aAtom in checkMol.GetAtoms():
                    aSym = aAtom.GetSymbol().strip()
                    if len(aSym) and not aSym in tAtomElems:
                        tAtomElems.append(aSym)
            else:
                print "Can not generate a molecule from the input SMILES string!"
                print "Check your SMILES or report the bug "

    """
    def getAtomElemsFromMol2(self, tFileName, tAtomElems):

        if os.path.isfile(tFileName):
            checkMol=Chem.MolFromMol2File(tFileName)
            if checkMol:
                for aAtom in checkMol.GetAtoms():
                    aSym = aAtom.GetSymbol().strip()
                    if len(aSym) and not aSym in tAtomElems:
                        tAtomElems.append(aSym)
            else:
                print "Can not generate a molecule from the input mol2 file!"
                print "Check your mol2 file or report the bug "
    """

    def getAtomElemsFromMol2(self, tFileName, tAtomElems):

        if os.path.isfile(tFileName):
            fM = open(tFileName, "r")
            allLs = fM.readlines()
            fM.close()
         
            lAtom = False
            for aL in allLs:
                if len(aL) > 0 and aL[0] != "#":
                    if aL.find("@<") !=-1 and aL.find("@<TRIPOS>ATOM") ==-1:
                        lAtom = False
                    elif lAtom :
                        strGrp = aL.strip().split()
                        if len(strGrp) >= 6:
                            aElem = ""
                            iPos =0
                            for aC in strGrp[1].strip():
                                if aC.isalpha() and not aC.isdigit(): 
                                    if iPos ==0:
                                        aElem = aElem + aC 
                                    else:
                                        if aC.lower() == aC:
                                            aElem = aElem + aC 
                                        else:
                                            break
                                else:
                                    break
                                iPos +=1
                            if len(aElem) and not aElem in tAtomElems:
                                tAtomElems.append(aElem)
                    elif aL.find("@<TRIPOS>ATOM") !=-1:
                        lAtom = True
        else: 
            print "%s can not be found for reading "%tFileName
                    
        
    def checkChiralCenters(self, tMol, tIdx):
   
        # RDKit misses some chiral centers such as N with 3 bonds.
        # Need to fix them here. May add other cases later.

        tCenAtm = tMol.GetAtomWithIdx(tIdx)

        aSymb = tCenAtm.GetSymbol().strip()
        if len(aSymb)==1 and aSymb.find("N") !=-1:
            aNBSet = tCenAtm.GetNeighbors()
            if len(aNBSet) == 3:
                lSP2 = False
                for aNB in aNBSet:
                    Hyb = aNB.GetHybridization()
                    if Hyb == Chem.rdchem.HybridizationType.SP2:
                        lSP2 = True
                        break

                if not lSP2:
                    # TEMP, re-calculated when the coordinates are available
                    #tCenAtm.SetChiralTag(rdchem.ChiralType.CHI_TETRAHEDRAL_CW)            
                    tCenAtm.SetProp("TmpChiral", "both")       
         
    def getChiralVolumeSign(self, cListCen, cList1, cList2, cList3):

        vSign = ""
        if len(cListCen)==3 and len(cList1)==3 and len(cList2)==3 and  len(cList3)==3:
            
            vect1 = []
            vect1.append(cList1[0]-cListCen[0])
            vect1.append(cList1[1]-cListCen[1])
            vect1.append(cList1[2]-cListCen[2])
       
            vect2 = []
            vect2.append(cList2[0]-cListCen[0])
            vect2.append(cList2[1]-cListCen[1])
            vect2.append(cList2[2]-cListCen[2])
       
            vect3 = [] 
            vect3.append(cList3[0]-cListCen[0])
            vect3.append(cList3[1]-cListCen[1])
            vect3.append(cList3[2]-cListCen[2])
       
            
            V =    vect1[0]*(vect2[1]*vect3[2]-vect2[2]*vect3[1]) \
                -  vect1[1]*(vect2[0]*vect3[2]-vect2[2]*vect3[0]) \
                +  vect1[2]*(vect2[0]*vect3[1]-vect2[1]*vect3[0]) \
         
            if V > 0.00001:
                vSign = "positive"
            elif V < -0.00001:
                vSign = "negative"
            else:
                vSign = "both"
    
            #print "Vol by acedrg :  ", V
            #print "Sign ", vSign 

        else :
            print "Can not calculate chiral volume. atom coordinates are not 3 dim "

        return vSign  
 


            
# Other supplement functions

def main():
    acedrgObj = Acedrg(sys.argv[1:])

if __name__ == "__main__":
    main()

