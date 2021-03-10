#!/usr/bin/env ccp4-python
# Python script
#
#
#     Copyright (C) 2014 --- 2020 Fei Long,  G. Murshudov
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Library.
#
#====================================================================
## The date of last modification: 21/02/2020
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
#from past.utils import old_div
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

from . exebase     import CExeCode

from . acedrgRDKit import AcedrgRDKit

from . chem        import ChemCheck

from . covLink     import CovLink
from . covLink     import CovLinkGenerator

from . filetools   import Ccp4MmCifObj
from . filetools   import FileTransformer

from . utility    import listComp
from . utility    import listComp2
from . utility    import listCompDes
from . utility    import listCompAcd
from . utility    import setBoolDict
from . utility    import splitLineSpa
from . utility    import splitLineSpa2
from . utility    import aLineToAlist
from . utility    import checkRepAtomTypes

if os.name != 'nt':
    import fcntl
import signal

class Acedrg(CExeCode ):

    def __init__( self, t_argvs):

        if len(sys.argv)==1:
            print("Look for help: %s -h "%(os.path.basename(sys.argv[0])))
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
        self.isAA             = False
        self.isNA             = False

        self.molGen           = False
        self.repCrds          = False
        self.neuDif           = False

        self.HMO              = False
        self.raS              = 0.5
        
        self.testMode         = False

        self.hasPDBxComp_descriptor = False
        self.outCifGlobSectHead    = []
        self.outCifGlobSect        = []

        self.allBondsAndAngles = {}

        self.upperSigForBonds   = 0.02
        self.lowSigForBonds     = 0.01
        self.upperSigForAngles  = 3.0
        self.lowSigForAngles    = 1.5
        

        self.allBondsAndAngles["atomClasses"] = {}
        self.allBondsAndAngles["bonds"]       = {}
        self.allBondsAndAngles["angles"]      = {}

        self.naTorsList                       = {}

        self.acedrg    = os.path.abspath(sys.argv[0])
        self.acedrgDir = os.path.dirname(os.path.dirname(self.acedrg))
        #self.acedrg    = ""
        #print "files ", glob.glob(sys.exec_prefix + "/*")
        #self.acedrgDir = sys.exec_prefix
        #print "files ", glob.glob(self.acedrgDir + "/*")
        self.qmInstructions   = ""
        self.qmSysDict        = {}
        inputOptionsP         = self.InputParser(t_argvs) 

        #if inputOptionsP.geneInFileName:
        #    self.checInputFormat()          

        self.runExitCode      = 0 

        self.setWorkMode(inputOptionsP)
        self.setInputProcPara(inputOptionsP)

        self.checkDependency()
        self.checkVersionInfo()
        self.showAcedrgPapers()

        if os.path.isfile(self.funcGroupTable):
            self.rdKit = AcedrgRDKit(self.funcGroupTable)
        else:
            self.rdKit = AcedrgRDKit()
        self.rdKit.setProcPara(self.inputPara)
        if self.useExistCoords or not self.workMode in [12, 121, 52]:
            self.rdKit.useExistCoords  = True        

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
        
    def showAcedrgPapers(self):

        print("=====================================================================")
        print("|                     Main reference                                |")
        print("| \"AceDRG: a stereochemical description generator for ligands\"      |")
        print("| Fei Long, Robert A. Nicholls, Paul Emsley, Saulius Grazulis,      |")
        print("| Andrius Merkys, Antanas Vaitkusb and Garib N. Murshudov,(2017)    |")
        print("| Acta Crystallogr. D73, 112-122                                    |")
        print("=====================================================================")

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
        self.inputParser.add_option("-a",  "--rechi", dest="ignoreInputChiral",  
                                    action="store_true",  default=False,
                                    help="igore chiral signs in input file, re-generate chiral signs by Acedrg")

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

        self.inputParser.add_option("-p",  "--coords", dest="useExistCoords",  
                                    action="store_true",  default=False,
                                    help="Using existing coordinates in the input file")

        self.inputParser.add_option("-r",  "--res", dest="monomRoot",  
                                    action="store", type="string", 
                                    help="The name of the chemical components users want to put into output files(e.g. PDB or MMCIF)")

        self.inputParser.add_option("-s",  "--sdf", dest="inSdfName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File of SDF format containing coordinates and bonds")

        self.inputParser.add_option("--bsu", dest="upperSigForBonds",  
                                    action="store", type="float", 
                                    help="Set upper bound for the sigma of bonds")

        self.inputParser.add_option("--bsl", dest="lowSigForBonds",  
                                    action="store", type="float", 
                                    help="Set low bound for the sigma of bonds")

        self.inputParser.add_option("--asu", dest="upperSigForAngles",  
                                    action="store", type="float", 
                                    help="Set upper bound for the sigma of angles")

        self.inputParser.add_option("--asl", dest="lowSigForAngles",  
                                    action="store", type="float", 
                                    help="Set low bound for the sigma of angles")

        self.inputParser.add_option("--fpar", dest="inParamFile", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File containing paramets that define upper and lower bounds for the sigma of bonds and angles")

        self.inputParser.add_option("--neu", dest="neuDif",  
                                    action="store_true",  default=False,  
                                    help="The option to look into structures (represented by cif files) determined by neutron diffraction")

        self.inputParser.add_option("--nucl", dest="protCol",
                                    action="store_true",  default=False,
                                    help="The option to add nucleus columns in mmcif (dictionary) files")

        self.inputParser.add_option("-t",  "--tab", dest="acedrgTables", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input path that stores all bond and angle tables (if no input, default CCP4 location will be used)")

        self.inputParser.add_option("-u",  "--hmo", dest="HMO",  
                                    action="store_true",  default=False,
                                    help="Calculate bond-orders and charges using HMO ")

        self.inputParser.add_option("-v",  "--version", dest="versionInfo",  
                                    action="store_true",  default=False,
                                    help="The option for checking version information of acedrg")

        self.inputParser.add_option("-x",  "--pdb", dest="inLigandPdbName", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File of PDB format containing coordinates of the ligand")

        self.inputParser.add_option("-y", "--repcrd",
                  action="store_true", dest="repCrd", default=False,
                  help="Use this keyword if you want to replace the atomic coordinates in the input mmCif with those in the input PDB")

        self.inputParser.add_option("-z",  "--noGeoOpt", dest="noGeoOpt",  
                                    action="store_true",  default=False,
                                    help="The option for not doing geometry optimization on coordinates")

        self.inputParser.add_option("-K",  "--noProt", dest="noProtonation",  
                                    action="store_true",  default=False,
                                    help="No further protonation/deprotonation will be done by Acedrg")

        self.inputParser.add_option("-L",  "--linkInstruction", dest="linkInstructions", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File that gives the instructions to build a link")

        self.inputParser.add_option("-Q",  "--qmInstruction", dest="qmInstructions", metavar="FILE", 
                                    action="store", type="string", 
                                    help="Input File that gives the instructions to do Quamtum Chemical calculations")

        self.inputParser.add_option("-T",  "--tempo", dest="testMode",  
                                    action="store_true",  default=False,  
                                    help="The mode temporarily exists in Acedrg, for some purposes")

        (inputOptionsP, inputOptionsU) = self.inputParser.parse_args(t_argvs)

        if inputOptionsU:    
            print("Line arguments for those keywords are missing ")
            for a_opt in inputOptionsU:
                print("Keyword %s ?"%a_opt)
            sys.exit()

        return inputOptionsP 

    def checkDependency(self):
    
        if self.acedrgDir !="" and os.path.isdir(self.acedrgDir):
            tLibmol = os.path.join(self.acedrgDir, "libexec","libmol")
            if platform.system()=="Windows": tLibmol += ".exe"
            if glob.glob(tLibmol):
                self.libmol = tLibmol

        if not self.libmol  and "CCP4" in os.environ:
            tLibmol = os.path.join(os.environ['CCP4'], "libexec","libmol")
            if platform.system()=="Windows": tLibmol += ".exe"
            if glob.glob(tLibmol):
                self.libmol = tLibmol

        if not self.libmol :
            if "LIBMOL_ROOT" in os.environ:
                tLibmol = os.path.join(os.environ['LIBMOL_ROOT'], "libexec","libmol")
                if platform.system()=="Windows": tLibmol += ".exe"
                if glob.glob(tLibmol):
                    self.libmol = tLibmol
                else:
                    print("libmol could not be found at %s"%tLibmol)
                    sys.exit()

        if not self.libmol: 
            print("can not find libmol at libexec/")
            sys.exit()
                                
        if "CCP4" in os.environ:
            tRefmac = os.path.join(os.environ['CBIN'], "refmac5")
            if platform.system()=="Windows": tRefmac += ".exe"
            if not glob.glob(tRefmac):
                print("refmac5 could not be found")
                sys.exit()
            else:
                self.refmac = tRefmac

            tLibcheck = os.path.join(os.environ['CBIN'], "libcheck")
            if platform.system()=="Windows": tLibcheck += ".exe"
            if not glob.glob(tLibcheck):
                print("libcheck could not be found")
                sys.exit()
            else:
                self.libcheck = tLibcheck
            
            if not self.libmol:
                tLibmol = os.path.join(os.environ['CCP4'], "libexec", "libmol")
                if platform.system()=="Windows": tLibmol += ".exe"
                if glob.glob(tLibmol):
                    self.libmol = tLibmol
        else :
            print("You need to install CCP4 suite")
            print("or activate ccp4.setup")
            sys.exit()
        if not self.acedrgTables:
            if "LIBMOL_ROOT" in os.environ:
                tAcedrgTables = os.path.join(os.environ['LIBMOL_ROOT'], "share","acedrg","tables")
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
                else:
                    print("Tables required could not be found at %s"%tAcedrgTables)
                    sys.exit()
        if not self.acedrgTables: 
            tAcedrgTables = os.path.join(self.acedrgDir, "share","acedrg","tables")
            # print tAcedrgTables
            if os.path.isdir(tAcedrgTables):
                self.acedrgTables = tAcedrgTables
        if not self.acedrgTables:
            tAcedrgTables = os.path.join(os.environ['CCP4'], "share","acedrg","tables")
            if glob.glob(tAcedrgTables):
                self.acedrgTables = tAcedrgTables
        if self.acedrgTables:
            tFuncGroupTable = os.path.join(self.acedrgTables, "funSmi.table")
            if os.path.isfile(tFuncGroupTable):
                self.funcGroupTable = tFuncGroupTable
            
        #print "The path to Acedrg tables is at ", self.acedrgTables
        #print "Libmol used is at ", self.libmol
        
    def checkVersionInfo(self):
  
        # Acedrg version info 
        self.versionInfo["man"] = os.path.join(self.acedrgTables, "manifest.txt")
        #print self.acedrgTables
        #print self.versionInfo["man"] 
        if not os.path.isfile(self.versionInfo["man"]):
            print("Version infomation is not available.")
        else:
            # print self.versionInfo["man"]
            try:
                vInfo = open(self.versionInfo["man"], "r")
            except IOError:
                print(self.versionInfo["man"], " Could not be opened for reading")
                print("Version infomation is not available.")
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
            print("Refmac version info is not available")

        if "REFMAC_VERSION" not in self.versionInfo:
            print("Refmac version info is not available")
           
             
    def setOutCifGlobSec(self):
        
        if not self.fileConv.hasStrDescriptors:
           for aL in self.fileConv.strDescriptors["defProps"]:
               self.outCifGlobSect.append(aL + "\n") 
        #self.outCifGlobSect.append("#loop_\n")
        #self.outCifGlobSect.append("#_software\n")
        #self.outCifGlobSect.append("#_version\n")
        #self.outCifGlobSect.append("#_purpose\n")
        if "ACEDRG_VERSION" in self.versionInfo:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), "acedrg".ljust(21), self.versionInfo["ACEDRG_VERSION"].strip().ljust(12), '\"dictionary generator\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), "acedrg".ljust(21), "?".ljust(12), '\"dictionary generator\"'.ljust(40)))
        
        if "DATABASE_VERSION" in self.versionInfo:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), "acedrg_database".ljust(21), self.versionInfo["DATABASE_VERSION"].strip().ljust(12), '\"data source\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), "acedrg_database".ljust(21), "?".ljust(12), '\"data source\"'.ljust(40)))

        self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17) ,"rdkit".ljust(21), rdBase.rdkitVersion.ljust(12), '\"Chemoinformatics tool\"' )) 
  
        if "REFMAC_NAME" in self.versionInfo:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), self.versionInfo["REFMAC_NAME"].ljust(21), self.versionInfo["REFMAC_VERSION"].ljust(12),\
                                       '\"optimization tool\"'.ljust(40)))
        else:
            self.outCifGlobSect.append("%s%s%s%s%s\n"%(self.monomRoot.ljust(4), "?".ljust(17), "refmac".ljust(21), "?".ljust(12), '\"optimization tool\"'.ljust(40)))
        
         
    def checInputFormat(self):
        
        try:
            tF =open(self.geneInFileName, "r")
        except IOError:
            print("%s can not be open for reading "%tFName)
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
            for aIdx in list(tDataDesc.keys()):
                if tDataDesc[aIdx][0].find("_chem_comp.id") !=-1:
                    self.monomRoot = tDataDesc[aIdx][1].strip()
                    return
        
        if self.inMdlName != "":
            try:
                tF =open(self.inMdlName, "r")
            except IOError:
                print("%s can not be open for reading "%self.inMdlName)
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
                print("%s can not be open for reading "%self.inMmCifName)
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

    def checkNAFromMmcif(self, tDataDesc=None):

        aRet = False
        if tDataDesc:
            for aK in tDataDesc.keys():
                if len(tDataDesc[aK])==2: 
                    if tDataDesc[aK][0].find("_chem_comp.type") !=-1 \
                       and (tDataDesc[aK][1].upper().find("DNA") !=-1 \
                            or tDataDesc[aK][1].upper().find("RNA") !=-1):
                        aRet = True 
                        break
        return aRet    

    def getNATors(self):

        torsFN = os.path.join(self.acedrgTables, "nucl_tors.table")
        if os.path.isfile(torsFN):
            torsF = open(torsFN, "r")
            allLs = torsF.readlines()
            torsF.close()
            
            if len(allLs) > 0:
                for aL in allLs:
                    aL = aL.strip()
                    print(aL)
                    if len(aL) > 0:
                        if aL[0].find("#") ==-1:
                            strgrp = aL.split()
                            if len(strgrp)==9:
                                print("here ")
                                torKey = strgrp[0].strip()
                                if not torKey in self.naTorsList.keys():
                                    self.naTorsList[torKey]    = {}
                                     
                                self.naTorsList[torKey][strgrp[1].strip()] = \
                                      [strgrp[2].strip(), strgrp[3].strip(), strgrp[4].strip(), strgrp[5].strip(),\
                                       strgrp[6].strip(), strgrp[7].strip(), strgrp[8].strip()]                           

        if len(self.naTorsList.keys()) > 0:
            print("special torsion angles for NA:")
            for aK in sorted(self.naTorsList.keys()):
                print("==========================")
                print("torsion key : %s "%aK)
                for aId in self.naTorsList[aK].keys():
                    print("ID          : %s "%aId)
                    print("atom Ids    :  ",self.naTorsList[aK][aId])
                    print("value       : %s "%self.naTorsList[aK][aId][4])
                    print("std_dev     : %s "%self.naTorsList[aK][aId][5])
                    print("period      : %s "%self.naTorsList[aK][aId][6])
        
    def setWorkMode(self, t_inputOptionsP = None):

        #print "acedrg is in ", self.acedrgDir
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
            if not self.acedrgTables or not os.path.isdir(self.acedrgTables):
                tAcedrgTables = os.path.join(self.acedrgDir, "tables")
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
            if  not self.acedrgTables or not os.path.isdir(self.acedrgTables):
                if "CCP4" in os.environ:
                    tAcedrgTables = os.path.join(os.environ['CCP4'], "share","acedrg","tables")
                    if os.path.isdir(tAcedrgTables):
                        self.acedrgTables = tAcedrgTables
                    else:
                        print("%s does not exist, check your installation of CCP4 "%tAcedrgTables)
                        sys.exit()
                else:
                    print("You need to install CCP4 suite")
                    sys.exit()
        else:
            if os.path.isdir(t_inputOptionsP.acedrgTables):
                self.acedrgTables = t_inputOptionsP.acedrgTables
 
        if not self.acedrgTables or not glob.glob(self.acedrgTables):
            if "LIBMOL_ROOT" in os.environ:
                tAcedrgTables = os.path.join(os.environ['LIBMOL_ROOT'], "share","acedrg","tables")
                print(tAcedrgTables)
                if os.path.isdir(tAcedrgTables):
                    self.acedrgTables = tAcedrgTables
        #print "Table is at ", self.acedrgTables

        if t_inputOptionsP.testMode :
            self.testMode = True

        if not t_inputOptionsP.molGen and not t_inputOptionsP.repCrd and not t_inputOptionsP.typeOut\
           and not t_inputOptionsP.HMO and not t_inputOptionsP.linkInstructions and\
           not t_inputOptionsP.qmInstructions and not t_inputOptionsP.testMode and not t_inputOptionsP.protCol: 
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
                if t_inputOptionsP.neuDif:
                    self.neuDif   = True
                    self.workMode = 211
                else:
                    self.workMode     = 21
            elif t_inputOptionsP.inStdCifDir: 
                self.inStdCifDir = t_inputOptionsP.inStdCifDir
                if t_inputOptionsP.neuDif:
                    self.neuDif   = True
                    self.workMode = 221
                else:
                    self.workMode    = 22
        elif t_inputOptionsP.typeOut:
            if t_inputOptionsP.inMmCifName:
                self.inMmCifName = t_inputOptionsP.inMmCifName
                self.workMode    = 311
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
        elif t_inputOptionsP.qmInstructions :
            self.qmInstructions = t_inputOptionsP.qmInstructions
            self.workMode = 70
        elif t_inputOptionsP.protCol :
            self.protCol = t_inputOptionsP.protCol
            self.workMode = 80
            if t_inputOptionsP.inMmCifName:
                self.inMmCifName = t_inputOptionsP.inMmCifName
        elif t_inputOptionsP.testMode :
            if t_inputOptionsP.inMmCifName:
                self.inMmCifName = t_inputOptionsP.inMmCifName
                self.workMode    = 1001
            
 
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
                    print("Could not put the input smiles in a file ")
                    sys.exit(1)
                else:
                    tF.write(self.inSmiName+ "\n")
                    tF.close()
                    self.inSmiName = tName
                    print(tName)


    def setInputProcPara(self, t_inputOptionsP = None):
       
        if t_inputOptionsP.upperSigForBonds:
            self.upperSigForBonds = t_inputOptionsP.upperSigForBonds
            print("upper bound for the bond sigma is set to ", self.upperSigForBonds)
        if t_inputOptionsP.lowSigForBonds:
            self.lowSigForBonds = t_inputOptionsP.lowSigForBonds
            print("low bound for the bond sigma is set to ", self.lowSigForBonds)
        if t_inputOptionsP.upperSigForAngles:
            self.upperSigForAngles = t_inputOptionsP.upperSigForAngles
            print("upper bound of the angle sigma is set to ", self.upperSigForAngles)
        if t_inputOptionsP.lowSigForAngles:
            self.lowSigForAngles = t_inputOptionsP.lowSigForAngles
            print("low bound of the angle sigma is set to ", self.lowSigForAngles)

        if t_inputOptionsP.inParamFile:
            self.inParamFile = t_inputOptionsP.inParamFile 
            print("Input paramet file is ", self.inParamFile)
            self.setSigmaBounds()


        if t_inputOptionsP.noProtonation:
            self.inputPara["noProtonation"]     = t_inputOptionsP.noProtonation

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

    def setSigmaBounds(self):
        
        if os.path.isfile(self.inParamFile):
            aPF = open(self.inParamFile, "r")
            allParaLs = aPF.readlines()
            aPF.close()
            for aLT in allParaLs:
                aL = aLT.upper().strip()
                if aL.find(":") != -1:
                    strs = aL.split(":") 
                    if len (strs) >= 2:
                        if strs[0].find("BSU") !=-1:
                            self.upperSigForBonds = float(strs[1])
                        if strs[0].find("BSL") !=-1:
                            self.lowSigForBonds = float(strs[1])
                        if strs[0].find("ASU") !=-1:
                            self.upperSigForAngles = float(strs[1])
                        if strs[0].find("ASL") !=-1:
                            self.lowSigForAngles = float(strs[1])
                    else:
                        print("Wrong format in the paramet file ",self.inParamFile)
                        print("Problem line at : ")
                        print(aLT)
                        sys.exit()
        else:
            print("Input paramet file %s does not exist"%self.inParamFile)
            sys.exit()        
        print("Upper bound for sigma of bonds ", self.upperSigForBonds)
        print("Lower bound for sigma of bonds ", self.lowSigForBonds)
        print("Upper bound for sigma of angles ", self.upperSigForAngles)
        print("Lower bound for sigma of angles ", self.lowSigForAngles)

    def printJobs(self):

        print("=====================================================================") 
        if "ACEDRG_VERSION" in self.versionInfo:
            print("| ACEDRG version:  %s|"%self.versionInfo["ACEDRG_VERSION"].ljust(49))
        else: 
            print("=====================================================================") 
            print("| ACEDRG version is not available                                   |")    
        if "DATABASE_VERSION" in self.versionInfo:
            print("| ACEDRG database: %s|"%self.versionInfo["DATABASE_VERSION"].ljust(49))
        else:
            print("| ACEDRG Database version is not available                          |")
        print("| RDKit version:  %s|"%rdBase.rdkitVersion.ljust(50))
        if "REFMAC_NAME" in self.versionInfo and "REFMAC_VERSION" in self.versionInfo:
            print("| %s  %s|"%((self.versionInfo["REFMAC_NAME"] + ":").ljust(15), self.versionInfo["REFMAC_VERSION"].ljust(49)))
        print("=====================================================================") 
        if self.workMode in [11, 12, 13, 14, 15, 16, 111, 121, 131, 141, 151, 161] :
        #if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 or self.workMode==15 \
        #   or self.workMode == 111 or self.workMode==121 or self.workMode ==131 or self.workMode==141 or self.workMode==151 :
            print("=====================================================================") 
            print("| Your job is  generating the dictionary (cif) and coord(pdb) files |")
            print("| for your ligand and/or monomer                                    |")
            print("=====================================================================") 
            if self.workMode==11 or self.workMode==111:
                print("Input file: %s"%os.path.basename(self.inMmCifName))
            if self.workMode==12 or self.workMode==121:
                print("Input file: %s"%os.path.basename(self.inSmiName))
            if self.workMode==13 or self.workMode==131:
                print("Input file: %s"%os.path.basename(self.inMdlName))
            if self.workMode==14 or self.workMode==141:
                print("Input file: %s"%os.path.basename(self.inSdfName))
            if self.workMode==15 or self.workMode==151:
                print("Input file: %s"%os.path.basename(self.inMol2Name))
            if self.workMode==16 or self.workMode==161:
                print("Input file: %s"%os.path.basename(self.inLigandPdbName))
            print("Output dictionary file: %s"%self.outRoot + ".cif")
            if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 or self.workMode==15:
                print("Output coordinate file: %s"%self.outRoot + ".pdb")

        if self.workMode in [21, 22]:
            print("=====================================================================") 
            print("| Your job is to generate molecules (sets of connected atoms), to   |")
            print("| get unique bonds and angles within the molecules and cluster them |") 
            print("| by specificted desigend atom types.                               |") 
            print("=====================================================================") 
            if self.workMode == 21:
                print("Input file: %s"%os.path.basename(self.inStdCifName))
                print("Output molecules: %s"%self.outRoot + "_all_mols.txt")
                print("Output bonds and angles : %s"%self.outRoot + "_unique_bond_and_angles.txt")
            elif self.workMode==22:
                print("Input directory where cif files are: %s"%self.inStdCifDir)
                print("Output bond and angle file : %s "%(self.outRoot + "_all_bonds_and_angles.table"))
        elif self.workMode in [211, 221]:
            print("=====================================================================") 
            print("| You would like to look into structures determined by neutron      |")
            print("| diffractions.                                                     |")
            print("| Your job is to generate molecules (sets of connected atoms), to   |")
            print("| get unique bonds and angles within the molecules and cluster them |") 
            print("| by specificted desigend atom types.                               |") 
            print("=====================================================================") 
            if self.workMode == 211:
                print("Input file: %s"%os.path.basename(self.inStdCifName))
                print("Output molecules: %s"%self.outRoot + "_all_mols.txt")
                print("Output bonds and angles : %s"%self.outRoot + "_unique_bond_and_angles.txt")
            elif self.workMode==221:
                print("Input directory where cif files are: %s"%self.inStdCifDir)
                print("Output bond and angle file : %s "%(self.outRoot + "_all_bonds_and_angles.table"))

        if self.workMode in [51, 52, 53, 54, 55] :
            print("=====================================================================") 
            print("| Your job is to calculate bond-orders and charges using HMO.      |")
            print("=====================================================================") 

        if self.workMode == 61:
            print("=====================================================================") 
            print("| Your job is to generate full descriptions for two monomers, their |")
            print("| modifications, and a link between them.                           |") 
            print("=====================================================================") 

    def runLibmol(self, tIn=None, tIdxMol=-1):
        self._cmdline = self.libmol
        if tIdxMol !=-1: 
            self._log_name       = os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(tIdxMol) + "_cod.log")
        else:
            self._log_name       = os.path.join(self.scrDir, self.baseRoot +  "_cod.log")

        if self.workMode == 11 or self.workMode==12 or self.workMode ==13 or self.workMode==14 \
           or self.workMode==15 or self.workMode ==16\
           or self.workMode == 111 or self.workMode == 121 or self.workMode == 131 \
           or self.workMode==141 or self.workMode==151 or self.workMode==161:
            self._cmdline += " -1 %f -2 %f -3 %f -4 %f "\
                            %(self.upperSigForBonds, self.lowSigForBonds,\
                              self.upperSigForAngles, self.lowSigForAngles)

            print("===================================================================") 
            print("| Generate the dictionary file using the internal database        |")
            print("===================================================================") 
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
           

        if self.workMode in [21, 211] :
            if tIn:
                self.inStdCifName = tIn
            self.outRstCifName  = self.outRoot + ".cif"
            self.outMolsName    = self.monomRoot + "_all_mols.txt"
            self.outBondsAndAnglesName  = self.monomRoot + "_unique_bond_and_angles.txt"

            self._cmdline += " -b %s  "%self.inStdCifName
            if self.workMode == 21:
                self._cmdline += " -m yes -r %s -o %s "%(self.monomRoot, self.outRoot)
            elif self.workMode == 211:
                self._cmdline += " -m yes -l yes -r %s -o %s "%(self.monomRoot, self.outRoot)
            #print self._cmdline
            #print "self.outRoot ", self.outRoot 
            self.runExitCode = self.subExecute()

        if self.workMode in [22, 221] :
            if os.path.isdir(self.inStdCifDir):
                bTable = self.outRoot + "_all_atoms_bonds_angles.table"
                self.workMode  == 21
                tAllBondInMols = []
                tExcludedCif   = []
                tempCifStr = os.path.join(self.inStdCifDir,"*.cif")
                for aCif in glob.glob(tempCifStr):
                    tMonomRoot = os.path.basename(aCif).split(".")[0]
                    print("Generate molecules, bonds and angles from ", aCif)
                    self._cmdline = self.libmol
                    self.outMolsName            = os.path.join(self.scrDir,tMonomRoot + "_all_mols.txt")
                    self.outBondsAndAnglesName  = os.path.join(self.scrDir, tMonomRoot + "_unique_bond_and_angles.txt")
                    self._log_name = os.path.join(self.scrDir,tMonomRoot + "_cod.log")
                    self._cmdline += " -b %s  "%aCif
                    tempMonRt = os.path.join(self.scrDir, tMonomRoot)
                    if self.workMode==22:
                        self._cmdline += " -m yes -r %s -o %s.cif "%(tMonomRoot, tempMonRt)
                    elif self.workMode==221:
                        self._cmdline += " -m yes -l yes -r %s -o %s.cif "%(tMonomRoot, tempMonRt)
                    #print self._cmdline
                    self.subExecute()
                    if os.path.isfile(self.outBondsAndAnglesName):
                        tAllBondInMols.append(self.outBondsAndAnglesName)
                    else:
                        print("No bonds and angles are generated from ", aCif)
                        tExcludedCif.append(aCif)
                if len(tAllBondInMols) !=0:
                    self.getBondsAndAngles(bTable, tAllBondInMols)
                if len(tExcludedCif) !=0:
                    print("The following cif files have been excluded from calculations because their large R factors")
                    for aExCif in tExcludedCif:
                        print(aExCif)

        if self.workMode == 31 or self.workMode ==311 or self.workMode == 32 or self.workMode == 33 or self.workMode == 34 or self.workMode == 35 : 
            if self.workMode == 31 or self.workMode == 311 or self.workMode == 32 or self.workMode==33:
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
                elif self.workMode == 311:
                    self.outAtmTypeName = self.outRoot + "_atomTypes.txt"
                    print("Output file name : %s "%self.outAtmTypeName)
                self._cmdline +=" -A yes -D %s -c %s  -o %s "%(self.acedrgTables, inFileName, self.outAtmTypeName)
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

        if self.workMode == 80 :
           
            self._cmdline += " -R yes -c %s -D %s  -o %s "%(tIn, self.acedrgTables, self.outProtDistTable)
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
                print("%s has not been found for reading"%aMonTab)
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
                         if tSets[0][1] not in self.allBondsAndAngles["bonds"]:         # class 1
                             self.allBondsAndAngles["bonds"][tSets[0][1]] = {}
                         if tSets[1][1] not in self.allBondsAndAngles["bonds"][tSets[0][1]]: # class 2 
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]] = {}
                         if tSets[0][0] not in self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]]: # id 1
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]] ={}
                         if tSets[1][0] not in self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]]: # id 2
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]] = {}
                         if strGrp[5] not in self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]]:
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]] = {}
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["observations"] = []
                             self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["stats"]        = {}

                         aTup = [tSets[0][2], tSets[1][2], float(strGrp[4]), aFileIdx]              # atom1_id, atom2_id, bond length and the file index
                         self.allBondsAndAngles["bonds"][tSets[0][1]][tSets[1][1]][tSets[0][0]][tSets[1][0]][strGrp[5]]["observations"].append(aTup) 
                

                         # Get sorted atom types
                         if tSets[0][1] not in self.allBondsAndAngles["atomClasses"]:            # atom 1 element type
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]] = {}
                         if tSets[0][0] not in self.allBondsAndAngles["atomClasses"][tSets[0][1]]: # atom 1 class type
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]][tSets[0][0]] = 1
                         else:
                             self.allBondsAndAngles["atomClasses"][tSets[0][1]][tSets[0][0]] +=1

                         if tSets[1][1] not in self.allBondsAndAngles["atomClasses"]:            # atom 2 element type
                              self.allBondsAndAngles["atomClasses"][tSets[1][1]] = {}
                         if tSets[1][0] not in self.allBondsAndAngles["atomClasses"][tSets[1][1]]: # atom 2 class type
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
                         if aSet[1] not in self.allBondsAndAngles["angles"]:         # center atom element  
                             self.allBondsAndAngles["angles"][aSet[1]] = {}
                         if tSets[0][1] not in self.allBondsAndAngles["angles"][aSet[1]]: # atom1 element 
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]] = {}
                         if tSets[1][1] not in self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]] : # atom2 elment 
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]] = {}
                         if aSet[0] not in self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]]: # center atom class
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]] = {}
                         if tSets[0][0] not in self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]]: #  atom1 class
                             self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]] = {}
                         if tSets[1][0] not in self.allBondsAndAngles["angles"][aSet[1]][tSets[0][1]][tSets[1][1]][aSet[0]][tSets[0][0]]: #  atom1 class
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
            print("%s has not been found for writing"%tFName)
        else:
            if len(self.allBondsAndAngles["atomClasses"]):
                f2.write("All unique atom types: \n")
                for aElem in sorted(self.allBondsAndAngles["atomClasses"].keys()):
                    for aCl in sorted(self.allBondsAndAngles["atomClasses"][aElem].keys()):
                        f2.write("%s        %s     \n"%(aElem,aCl))   
            if len(self.allBondsAndAngles["bonds"]):
                f2.write("All bond lengths: \n")  
                for aElem in sorted(self.allBondsAndAngles["bonds"].keys()):
                    for bElem in sorted(self.allBondsAndAngles["bonds"][aElem].keys()):
                        for aCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem].keys()):
                            for bCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl].keys()):
                                for rP in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl].keys()):
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
                     
                
                for aElem in sorted(self.allBondsAndAngles["bonds"].keys()):
                    for bElem in sorted(self.allBondsAndAngles["bonds"][aElem].keys()):
                        for aCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem].keys()):
                            for bCl in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl].keys()):
                                for rP in sorted(self.allBondsAndAngles["bonds"][aElem][bElem][aCl][bCl].keys()):
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
                for cenElem in sorted(self.allBondsAndAngles["angles"].keys()):
                    for aElem in sorted(self.allBondsAndAngles["angles"][cenElem].keys()):
                        for bElem in sorted(self.allBondsAndAngles["angles"][cenElem][aElem].keys()):
                            for cenC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem].keys()):
                                for aC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC].keys()):
                                    for bC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC].keys()):
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
                for cenElem in sorted(self.allBondsAndAngles["angles"].keys()):
                    for aElem in sorted(self.allBondsAndAngles["angles"][cenElem].keys()):
                        for bElem in sorted(self.allBondsAndAngles["angles"][cenElem][aElem].keys()):
                            for cenC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem].keys()):
                                for aC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC].keys()):
                                    for bC in sorted(self.allBondsAndAngles["angles"][cenElem][aElem][bElem][cenC][aC].keys()):
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
                    print("Bug. error in ", aTup)
                    sys.exit(1)

            tBondStatsDict["mean"] = old_div(sum,numObjs)

            sum_diff_sq = 0.0

            for aTup in tBondObjList:
                sum_diff_sq +=((aTup[2]- tBondStatsDict["mean"])*(aTup[2]- tBondStatsDict["mean"]))
                              
            if numObjs > 1:
                tBondStatsDict["sig"] = math.sqrt(old_div(sum_diff_sq,(numObjs-1)))
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
                    print("Bug. error in ", aTup)
                    sys.exit(1)

            tAngleStatsDict["mean"] = old_div(sum,numObjs)

            sum_diff_sq = 0.0

            for aTup in tAngleObjList:
                sum_diff_sq +=((aTup[3]- tAngleStatsDict["mean"])*(aTup[3]- tAngleStatsDict["mean"]))

            if numObjs > 1:
                tAngleStatsDict["sig"] = math.sqrt(old_div(sum_diff_sq,(numObjs-1)))
            else:
                tAngleStatsDict["sig"] = 0.0  
            
            
    def setLibcheckBat(self, tSmiName, tOutRoot):

        self._cmdline = self.libcheck 

        libcheckBatName = os.path.join(self.scrDir, self.baseRoot + "_libcheck.bat")
        try:
            libcheckBat = open(libcheckBatName,"w")
        except IOError:
            print(libcheckBatName, " could not be opened for write ")
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
                self._cmdline = "(ECHO ligand && ECHO ncyc 40 && ECHO make hout yes && ECHO make hydr full && ECHO end) | " + self._cmdline
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
      
        print("| Stage  %d        |"%tStage)
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
                                print("===================================================================") 
                            else:
                                print("Failed to produce %s after final geometrical optimization"%finPdb)
                        else:
                            print("Failed to produce the  coordinates at stage 3 optimization") 
                    else:
                        print("Failed to produce the coordinates at stage 2 optimization") 
                else:
                    print("Failed to produce the coordinates at stage 1 optimization") 
            else:
                print("No dictionary file produced ! ") 
        else:
            print("No dictionary file produced ! ") 
                        
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
                    print("Failed to produce the coordinates for input file %s in optimization"%tPdbIn) 
                else:
                    tFValue = -1.0 
                    tFValue = self.getRefmacFValue(self._log_name)
                    if tFValue >= -0.000001:
                        if len(tRoot) !=0:
                            if tRoot not in self.refmacXYSList:
                                 self.refmacXYSList[tRoot] = {}
                            self.refmacXYSList[tRoot]["log"]    = self._log_name 
                            self.refmacXYSList[tRoot]["xyz"]    = self.refmacXYZOUTName 
                            self.refmacXYSList[tRoot]["fvalue"] = tFValue 
                            # print "| %s|"%("Done: " + tRoot).ljust(64) 
                        aRefPair = [tFValue, self.refmacXYZOUTName]
                        self.refmacMinFValueList.append(aRefPair)
            else:
                print("No dictionary file for optimization ! ") 
        else:
            print("No input pdb and/or dictionary files for optimization ! ") 

    def runGeoOptOneMolFull(self, tIdxMol):

        print("Number of final output conformers for molecule %d is %d "%(tIdxMol+1, self.numConformers))
        tmpStr = ""
        if self.numConformers == 1:
            tmpStr = "_tmp"
        #nConf = self.rdKit.molecules[tIdxMol].GetNumConformers()
        #print "Number of intial conformers for the molecule  ", nConf
        print("Number of intial conformers for refmac geo-opt  ", len(self.rdKit.selecConformerIds))

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
            #print ("|%s%s|"%("Input XYZ : ".ljust(12), aPdbIn.ljust(53)))
            #print ("|%s%s|"%("Input LIB : ".ljust(12), aLibCifIn.ljust(53)))
            self.runGeoOpt(aFRoot, aPdbIn, aLibCifIn) 
            if  self.runExitCode :
                print("Geometrical optimization fails to produce the final coordinates for %s after geometrical optimization"%aPdbIn)
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
            print("Error %s can not be opened for reading"%tLogName)
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
                finTor = self.outRoot +  "_"+ tRoot +"_torsion_only.txt"
                tStr = tRoot.strip().split("_")[-1]
            else:
                finPdb = self.outRoot +  ".pdb"
                finRst = self.outRoot +  ".cif"
                finTor = self.outRoot +  "_torsion_only.txt"
            #print "tInPdb ", tInPdb 
            #print "finPdb ", finPdb

            shutil.copy(tInPdb, finPdb)

            if os.path.isfile(finPdb):
                self.inPdbName        = finPdb 
                self.inMmCifName      = tInCif
                self.outRstCifName    = finRst
                self.transCoordsPdbToCif(self.inPdbName, self.inMmCifName, self.outRstCifName, tMol, tDataDescriptor, tStrDescriptors, tDelocList)
                #self.outTorsionRestraints(self.outRstCifName, finTor)
        else:
            print("Failed to produce %s after final geometrical optimization"%tInPdb)

    def transCoordsPdbToCif(self, tPdbInName, tCifInName, tCifOutName, tMol=-1, tDataDescriptor=None, tStrDescriptors=None,tDelocList=None):

        cifCont = {}
        cifCont['head'] = []
        monoId = ""

        if self.monomRoot.find("UNL") ==-1:
            if len(self.monomRoot) >=3:
                monoId = self.monomRoot[:3]
            else:
                monoId = self.monomRoot
 
        elif tDataDescriptor:
            monoId = tDataDescriptor[-1].strip().split()[0]
        else:
            monoId = "UNL"

        #print monoId 

        if tDataDescriptor:
            cifCont['head']   = ["#\n", "data_comp_list\n", "loop_\n"]
            for aL in tDataDescriptor:
                if aL[0].find("_")==-1:
                    if len(monoId) < 3:
                        aNewL = "%s%s"%(monoId.ljust(6),monoId.ljust(6))+ aL[10:]
                    else:
                        aNewL = "%s%s"%(monoId.ljust(8),monoId.ljust(8))+ aL[15:]
                    cifCont['head'].append(aNewL+"\n")   
                else:
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
            print("%s can not be opened for reading"%tPdbInName)
            sys.exit()
        else:
            try: 
                tCifIn = open(tCifInName, "r")
            except IOError:
                print("%s can not be opened for reading"%tCifInName)
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
                            
                                if tName in pdbAtoms:
                                    bLine = "%s%s%s%s%s%s%s%s\n"%(monoId.ljust(8), tID.ljust(8), strGrp[2].ljust(8), \
                                                                strGrp[3].ljust(8), strGrp[4].ljust(8), \
                                                                pdbAtoms[tName][0].ljust(12), pdbAtoms[tName][1].ljust(12), \
                                                                pdbAtoms[tName][2].ljust(12)) 
                                    cifCont['atoms'].append(bLine)
                                else: 
                                    print("Bug. can not find atom %s in Pdb file %s "%(tName, tPdbInName)) 
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
                            
                            if tName in pdbAtoms:
                                bLine = "%s%s%s%s%s%s%s%s\n"%(strGrp[0].ljust(8), tID.ljust(8), strGrp[2].ljust(8), \
                                                              strGrp[3].ljust(8), strGrp[4].ljust(8), \
                                                              pdbAtoms[tName][0].ljust(12), pdbAtoms[tName][1].ljust(12), \
                                                              pdbAtoms[tName][2].ljust(12)) 
                                cifCont['atoms'].append(bLine)
                            else: 
                                print("Bug. can not find atom %s in Pdb file %s "%(tName, tPdbInName)) 
                                sys.exit()
                        elif lOther:
                            cifCont['others'].append(aLine)
 
                if "bonds" in cifCont and len(cifCont['bonds']) !=0:
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
                    print("Out cif name ", tCifOutName)
                except IOError:
                    print("%s can not be opened for reading"%tCifOutName)
                    sys.exit()
                else:
                            
                    for aL in cifCont['head']:
                        tOutCif.write(aL)
                    for aL in cifCont['atoms']:
                        tOutCif.write(aL)
                    if tDataDescriptor:
                        if not self.isNA:
                            for aL in cifCont['others1']:
                                strGrp = aL.strip().split()
                                if len(strGrp)==1:
                                    tOutCif.write(aL)
                                    print(aL)
                                else:
                                    aL1 = monoId + aL[3:]
                                    tOutCif.write(aL1)
                                    print(aL1)
                            if 'bonds' in cifCont:
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
                            self.naCorr(cifCont['others1'], tOutCif) 
                    else:
                        print("Here 2")
                        for aL in cifCont["others"]:
                            strGrp = aL.strip().split()
                            if len(strGrp)==1:
                                tOutCif.write(aL)
                            else:
                                aL1 = monoId + aL[3:]
                                print(aL1)
                                tOutCif.write(aL1)
                    
                    aPos = len(monoId)   
                    if tStrDescriptors and "props" in tStrDescriptors and len(tStrDescriptors["props"]) !=0:
                        aSetChars = ["#", ";", "\""]
                        tOutCif.write("loop_"+"\n")
                        for aL in tStrDescriptors["props"]:
                            tOutCif.write(aL+"\n")
                        for aL in tStrDescriptors["entries"]:
                            if len(aL.strip()) >0 and not aL[0] in aSetChars: 
                                aL1 = monoId + aL[aPos:]
                            else:
                                aL1 = aL
                            tOutCif.write(aL1+"\n")
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
                    if len(self.outCifGlobSect):
                        for aL in self.outCifGlobSect:
                            tOutCif.write(aL)

                    tOutCif.close()
    
    def naCorr(self, tCifCont, tOutF):

        lTor = False
        speTors   = {}
        otherTors = []
        for aL in tCifCont:
            if aL.find("_chem_comp_tor.period") != -1:
                tOutF.write(aL)
                lTor = True
            elif lTor and aL.find("loop_")==-1:
                if len(aL) >0:
                    if aL[0].find("#") == -1:
                        strs = aL.strip().split()
                        if len(strs)==9:
                            id1 = strs[2] + "_" + strs[3]  + "_" + strs[4] + "_" + strs[5]                      
                            id2 = strs[5] + "_" + strs[4]  + "_" + strs[3] + "_" + strs[2] 
                            id  = ""
                            if id1 in self.naTorsList.keys():
                                id = id1
                            elif id2 in self.naTorsList.keys():
                                id = id2
                            else:
                                otherTors.append(aL)
                            print("id1 is ", id1)
                            print("id2 is ", id2)
                            print("id is ",  id)
                            if len(id) > 0:
                                for aId in sorted(self.naTorsList[id].keys()):
                                    aTorL = "%s%s%s%s%s%s%s%s%s\n"%(strs[0].ljust(8), aId.ljust(16),\
                                        self.naTorsList[id][aId][0].ljust(10), self.naTorsList[id][aId][1].ljust(10),\
                                        self.naTorsList[id][aId][2].ljust(10), self.naTorsList[id][aId][3].ljust(10),\
                                        self.naTorsList[id][aId][4].ljust(14), self.naTorsList[id][aId][5].ljust(10),\
                                        self.naTorsList[id][aId][6].ljust(6))
                                    speTors[aId]= aTorL
                        else:
                            otherTors.append(aL)
            elif lTor and aL.find("loop_") !=-1:
                for aId in sorted(speTors.keys()):
                    tOutF.write(speTors[aId])
                for aTorL in otherTors:
                    tOutF.write(aTorL)
                print("Here 4")
                print(aL)
                lTor = False
                tOutF.write(aL)
            else:
                print("Here 5")
                print(aL)
                tOutF.write(aL)
                
    def outTorsionRestraints(self, tCifInName, tTorOutName):

        if os.path.isfile(tCifInName):
            aCifObj = Ccp4MmCifObj(tCifInName)
            if self.monomRoot in aCifObj["ccp4CifObj"]["comps"]:
                if "tors" in aCifObj["ccp4CifObj"]["comps"][self.monomRoot]\
                   and "atoms" in aCifObj["ccp4CifObj"]["comps"][self.monomRoot]\
                   and "atoms" in aCifObj["ccp4CifObj"]["comps"][self.monomRoot]:
                    try: 
                        aTorOut = open(tTorOutName, "w")
                    except IOError:
                        print("%s can not be opened for reading"%tTorOutName)
                    else:
                        aSetTor = []
                        for aTor in aCifObj["ccp4CifObj"]["comps"][self.monomRoot]["tors"]: 
                            id1 = aTor["atom_id_1"]
                            atom1 = aCifObj.getAtomById(id1, self.monomRoot)
                            id2 = aTor["atom_id_2"]
                            atom2 = aCifObj.getAtomById(id2, self.monomRoot)
                            id3 = aTor["atom_id_3"]
                            atom3 = aCifObj.getAtomById(id3, self.monomRoot)
                            id4 = aTor["atom_id_4"]
                            atom4 = aCifObj.getAtomById(id4, self.monomRoot)
                            if atom1 and atom2 and atom3 and atom4:
                                aBond = aCifObj.getBondByIds(id2, id3,  self.monomRoot)
                                if aBond:
                                    if aBond["type"].lower().find("sing") !=-1:
                                        pass 
          
 
    def outEnergyGeoMap(self, tIdxMol):

        aListFName = os.path.join(self.scrDir, self.baseRoot + "_energy_vs_Conformers.list")
        try: 
            aListF = open(aListFName, "w")
        except IOError:
            print("%s can not be opened for reading"%aListFName)
        else:
            tmpStr = ""
            if self.numConformers == 1:
                tmpStr = "_tmp"
            aListF.write("RDKit-Energy\tOriginal-Conformer-Id\tPDB-File-Name\n")
            nId = 1
            for aEng in sorted(self.rdKit.conformerEngMap.keys()):
                for aCId in self.rdKit.conformerEngMap[aEng]:
                    aPdbName = "mol_" + str(tIdxMol+1) + "_conf_" + str(nId) + tmpStr + "_init.pdb"
                    nId +=1   
                    aListF.write("%8.4f\t%d\t%s\n"%(aEng, aCId, aPdbName))

            aListF.write("REFMAC-Energy\tPDB-File-Name-With-Original-Conformer-Id\n")
            for aPair in self.refmacMinFValueList:
                aPair[1]= aPair[1].strip().split("/")[-1].strip()
                aListF.write("%8.4f\t%s\n"%(aPair[0], aPair[1]))
            aListF.close()

      
    def WriteSimplifiedMmcif(self, tMmcifName, tMonoName, tMonomer):
        
        # A simplified mmcif file contains:
        # (1) Header section      
        # (2) Description of atoms in the molecule
        # (3) Description of bonds in the molecule
    
        # This file is mainly used as an input file for Acedrg
        print("Ligand ID ", tMonoName)
       
        allAtoms = tMonomer["atoms"] 
        nAt  = len(allAtoms)
        print("number of atoms", nAt)
        allBonds = tMonomer["bonds"]
        nBo  = len(allBonds)
        print("number of bonds", nBo)

        tGroupName="non-polymer"    # place-holder, will not  be used 
        nHAt = nAt                  # place-holder, will not  be used
        try:
            aMmCif = open(tMmcifName, "w")
        except IOError:
            print(tMmcifName, " Could not be opened for reading")
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
        
            # Atom section
            lCP = False
            for aAtom in allAtoms: 
                if "partial_charge" in list(aAtom.keys()):
                    lCP = True
                break

            lCO = True
            for aAtom in allAtoms: 
                if not "x" in list(aAtom.keys()):
                    lCO = False
                break


            aMmCif.write("loop_\n")
            aMmCif.write("_chem_comp_atom.comp_id\n")
            aMmCif.write("_chem_comp_atom.atom_id\n")
            aMmCif.write("_chem_comp_atom.type_symbol\n")
            aMmCif.write("_chem_comp_atom.type_energy\n")
            if not lCP:
                aMmCif.write("_chem_comp_atom.charge\n")
            else:
                aMmCif.write("_chem_comp_atom.partial_charge\n")
            if lCO:
                aMmCif.write("_chem_comp_atom.x\n")
                aMmCif.write("_chem_comp_atom.y\n")
                aMmCif.write("_chem_comp_atom.z\n")
            #nTetraChi = 0 
            for aAtom in allAtoms:
                if not lCP :
                    if lCO:
                        aMmCif.write("%s         %s      %s    %s     %3.2f   %5.4f    %5.4f     %5.4f\n"
                                      %(tMonoName, aAtom["atom_id"], aAtom["type_symbol"],
                                      aAtom["type_energy"], float(aAtom["charge"]), 
                                      float(aAtom["x"]), float(aAtom["y"]), float(aAtom["z"]) ))
                    else:
                        aMmCif.write("%s         %s      %s    %s     %3.2f   \n"
                                      %(tMonoName, aAtom["atom_id"], aAtom["type_symbol"],
                                      aAtom["type_energy"], float(aAtom["charge"]) ))

                else :
                    if lCO:
                        aMmCif.write("%s         %s      %s    %s     %3.2f   %5.4f    %5.4f     %5.4f\n"
                                      %(tMonoName, aAtom["atom_id"], aAtom["type_symbol"],
                                        aAtom["type_energy"], float(aAtom["partial_charge"]), 
                                        float(aAtom["x"]), float(aAtom["y"]), float(aAtom["z"]) ))
                    else:
                        aMmCif.write("%s         %s      %s    %s     %3.2f   \n"
                                      %(tMonoName, aAtom["atom_id"], aAtom["type_symbol"],
                                      aAtom["type_energy"], float(aAtom["partial_charge"]) ))
                   
            # Bond section
            lAR = True
            for aBond in allBonds:
                if not "aromatic" in list(aBond.keys()):
                    lAR = False
                break

            aMmCif.write("#\n")
            aMmCif.write("_chem_comp_bond.comp_id\n")
            aMmCif.write("_chem_comp_bond.atom_id_1\n")
            aMmCif.write("_chem_comp_bond.atom_id_2\n")
            aMmCif.write("_chem_comp_bond.type\n")
            if lAR:
                aMmCif.write("_chem_comp_bond.aromatic\n")
            aMmCif.write("_chem_comp_bond.value_dist\n")
            aMmCif.write("_chem_comp_bond.value_dist_esd\n")
            for aBond in allBonds:
                name1 = aBond["atom_id_1"]
                name2 = aBond["atom_id_2"]

                bType = aBond["type"]
                isAro =""
                if lAR:
                    isAro = aBond["aromatic"]
                
                bLen  = float(aBond["value_dist"])
                dBlen = float(aBond["value_dist_esd"])
                if lAR:
                    aMmCif.write("%s       %s       %s       %s      %s     %5.4f     %5.4f\n" \
                                  %(tMonoName, name1, name2,  bType, \
                                    isAro, bLen, dBlen))
                else:
                    aMmCif.write("%s       %s       %s       %s       %5.4f     %5.4f\n" \
                                  %(tMonoName, name1, name2,  bType, bLen, dBlen))

            aMmCif.close()
            
    def addProtCols(self):
 
        aMmcifObj = Ccp4MmCifObj(self.inMmCifName)

        a3Blocks  = []
        self.get3Blocks(a3Blocks)    # tempo

        aNewCif   = self.baseRoot + ".cif"     
        if len(list(aMmcifObj["ccp4CifObj"]["comps"].keys()))==1:
            aMonoName = list(aMmcifObj["ccp4CifObj"]["comps"].keys())[0]
            aMonoObj  = aMmcifObj["ccp4CifObj"]["comps"][aMonoName]
            if len(aMonoObj["atoms"]) > 1:
                aSimpCifName =  os.path.join(self.scrDir, self.baseRoot + "_tmp.cif")
                print(aSimpCifName)
                self.WriteSimplifiedMmcif(aSimpCifName, aMonoName, aMonoObj)   
                if os.path.isfile(aSimpCifName):
                    self.outProtDistTable = os.path.join(self.scrDir, self.baseRoot + "_protDist.table")
                    self.runLibmol(aSimpCifName)
                
                if os.path.isfile(self.outProtDistTable):
                    aBondSet = {}
                    self.readProTab(aBondSet)
                    self.getNewCif(aNewCif, aBondSet, a3Blocks)
            else:
                shutil.copyfile(self.inMmCifName, aNewCif)
        else:
            shutil.copyfile(self.inMmCifName, aNewCif)

    def get3Blocks(self, t3B):

        fM = open(self.inMmCifName, "r")
        allLs = fM.readlines()
        fM.close()

        lB1 = True
        lB2 = False
        aBl = []
        for aL in allLs:
            if aL.find("_chem_comp_bond.value_dist") !=-1 and aL.find("_chem_comp_bond.value_dist_")==-1:
                t3B.append(aBl)
                lB1 = False
                lB2 = True
                aBl = []
                aBl.append(aL)
            elif  lB2 :
                if aL.find("loop_") != -1  :
                    t3B.append(aBl)
                    lB2= False
                    aBl = []
                    aBl.append(aL)
                else:
                    aBl.append(aL)
            else:    
                aBl.append(aL)
        t3B.append(aBl) 

    def readProTab(self, tBondSet):
 
        fP = open(self.outProtDistTable, "r")
        allLs = fP.readlines()
        fP.close()

        for aL in allLs:
            strGrp = aL.strip().split()
            if len(strGrp)== 6:
                combo_ID = strGrp[0] + "_" + strGrp[1]
                print(combo_ID)
                tBondSet[combo_ID] = {}
                tBondSet[combo_ID]["atom_id_1"] = strGrp[0]
                tBondSet[combo_ID]["atom_id_2"] = strGrp[1]
                tBondSet[combo_ID]["prot_h"]    = strGrp[2]
                tBondSet[combo_ID]["prot_h_s"]  = strGrp[3]
                tBondSet[combo_ID]["e_h"]       = strGrp[4]
                tBondSet[combo_ID]["e_h_s"]     = strGrp[5]
                print("EL = %s"%tBondSet[combo_ID]["e_h"])   

    def getNewCif(self, tCif, tBondSet, t3Bs):

        fO = open(tCif, "w")

        lFN = False
        for aL in t3Bs[0]:
            fO.write(aL)
            if aL.find("_chem_comp_bond.value_dist_nucleus") != -1:
                lFN = True  

        if not lFN:
            fO.write("_chem_comp_bond.value_dist_nucleus\n")
            fO.write("_chem_comp_bond.value_dist_nucleus_esd\n")

        aIDList = list(tBondSet.keys())
        print(aIDList)
        for aL in t3Bs[1]:
            outL = ""
            strs = aL.strip().split()
            nStrs = len(strs)
            if len(strs)==1:
                fO.write(aL)
            elif len(strs) > 0:
                if not lFN  and strs[0].find("#")==-1:
                    id1 = strs[1]
                    if id1.find("\'") !=-1 and id1.find("\"")==-1:
                        id1 = "\"" + id1 + "\""
                    id2 = strs[2]
                    if id2.find("\'") !=-1 and id2.find("\"")==-1:
                        id2 = "\"" + id2 + "\""
                    compId1 = id1 + "_" + id2
                    compId2 = id2 + "_" + id1
                    print("Id1 ", id1)
                    print("Id2 ", id2)

                    if compId1 in aIDList or compId2 in aIDList:
                        if compId1 in aIDList:
                            v  = "%4.3f"%float(tBondSet[compId1]["prot_h"])
                            vs = "%4.3f"%float(tBondSet[compId1]["prot_h_s"])
                            e  = "%4.3f"%float(tBondSet[compId1]["e_h"])
                            es = "%4.3f"%float(tBondSet[compId1]["e_h_s"])
                        elif compId2 in aIDList:
                            v  = "%4.3f"%float(tBondSet[compId2]["prot_h"])
                            vs = "%4.3f"%float(tBondSet[compId2]["prot_h_s"])
                            e  = "%4.3f"%float(tBondSet[compId2]["e_h"])
                            es = "%4.3f"%float(tBondSet[compId2]["e_h_s"])
                        for i in range(nStrs-2):
                            outL += "%s"%(strs[i].ljust(10))
                        outL+="%s%s%s%s\n"%(v.ljust(10), vs.ljust(10), e.ljust(10), es.ljust(10))
                        fO.write(outL)
                    else:
                        for i in range(nStrs):
                            outL += "%s"%(strs[i].ljust(10))
                        outL+="%s%s\n"%(strs[-2].ljust(10), strs[-1].ljust(10))
                        fO.write(outL)
                else :
                    for i in range(nStrs):
                        outL += "%s"%(strs[i].ljust(10))
                    fO.write(outL+"\n")
        
        if len(t3Bs) > 2:            
            for aL in t3Bs[2]:
                fO.write(aL)

        fO.close()

    def tmpCifInitMol(self):

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

    def execute(self):
        
        self.printJobs()
        if self.workMode == 11 or self.workMode == 111:

            # Stage 1: dictionary generation using a mmcif file 
            if os.path.isfile(self.inMmCifName):
                if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                    print("The input system contains metal or other heavier element")
                    print("The current version deals only with the atoms in the set of 'organic' elements") 
                    sys.exit()
                self.runLibmol(self.inMmCifName)
            else:
                print("The input %s does not exist"%self.inMmCifName)
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
                print("%s does not exist"%self.inSmiName)
                sys.exit()
            # Stage 2: dictionary generation using a mmcif file
            if os.path.isfile(self.libcheckLibName):
                if not self.chemCheck.isOrganic(self.libcheckLibName, 11):
                    print("The input system contains metal or other heavier element")
                    print("The current version deals only with the atoms in the set of 'organic' elements") 
                    sys.exit()
                self.runLibmol(self.libcheckLibName)
            else:
                print("The input smiles contain metal elements ")
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
                print("%s does not exist"%self.inSdfName)
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
                print("%s does not exist"%self.inMol2Name)
                sys.exit()
           
            if self.workMode == 15 and not self.runExitCode: 
                # Stage 2: optimization
                self.runGeoOpt()            
            else:
                self.printExitInfo() 

        if self.workMode in [21, 211]:
            
            # Stage 1: generate molecules and the associated bond and bond-angle values 
            # using a small molecule cif file
            if os.path.isfile(self.inStdCifName):
                self.runLibmol(self.inStdCifName)
            else:
                print("Can not find the input file ", self.inStdCifName) 
            
        if self.workMode in [22, 221]:
            
            # 1. Generate molecules using the small molecule cif files at the input directory. 
            # 2. Generate atom classes for atoms in the molecules.
            # 3. Obtain unique bond lengths and angles and cluster them according to their 
            #    component atoms in tables.
            if os.path.isdir(self.inStdCifDir):
                self.runLibmol()
            else:
                print("Can not find the input directory ", self.inStdCifDir)
                
        if self.workMode == 31 or self.workMode == 311:
            print ("===================================================================") 
            print ("| Generate atom types of Acedrg style                             |")
            print ("===================================================================") 
            if os.path.isfile(self.inMmCifName):
                self.runLibmol()    
        
        if self.workMode == 32:
            print("work mode ", self.workMode)
            if os.path.isfile(self.inSmiName):
                self.runLibcheck(self.inSmiName)
            else:
                print("%s does not exist"%self.inSmiName)
                sys.exit()

            if os.path.isfile(self.libcheckLibName):
                self.runLibmol(self.libcheckLibName)
            else:
                print("%s does not exist"%self.inMmCifName)
                sys.exit()
            
        if self.workMode == 33:
            print("work mode ", self.workMode)
            if os.path.isfile(self.inMdlName):
                self.runLibmol()    
        

        if self.workMode == 34:
            print("work mode ", self.workMode)
            if os.path.isfile(self.inSdfName):
                self.runLibmol()    

        if self.workMode == 35:
            print("work mode ", self.workMode)
            if os.path.isfile(self.inStdCifName):
                self.runLibmol()    

        if self.workMode ==111 or self.workMode ==121 or self.workMode ==131 or self.workMode ==141:
            if os.path.isfile(self.outRstCifName):
                tCif = self.outRoot + ".cif"
                #os.system("cp %s %s"%(self.outRstCifName, tCif)) 
                shutil.copy(self.outRstCifName, tCif) 
            else:
                print("acedrg failed to generate a dictionary file")       
    
    def getAAOut(self):
    
        ccpMonLib   = os.getenv("CLIBD_MON") 
        aSub        = self.monomRoot[0].lower() 
        aaDirC = os.path.join(ccpMonLib, aSub)       
        aaDir = os.path.join(self.acedrgTables, "AminoAcids")       
        iniPdb = os.path.join(aaDir, self.monomRoot + ".pdb")
        #print "iniPdb ", iniPdb
        finPdb = self.outRoot + ".pdb"
        if os.path.isfile(iniPdb):
            shutil.copy(iniPdb, finPdb)
        else:
            print("Error in dealing with the pdb of amino acid %s "%self.monomRoot) 
            sys.exit(1)
        iniCif = ""
        if os.path.isdir(aaDirC):
            iniCif = os.path.join(aaDirC, self.monomRoot + ".cif")
        else:
            iniCif = os.path.join(aaDir, self.monomRoot + ".cif")
        finCif = self.outRoot + ".cif"
        if os.path.isfile(iniCif):
            shutil.copy(iniCif, finCif)
        else:
            print("Can not find the cif file for  amino acid %s "%self.monomRoot) 
            sys.exit(1)
        if os.path.isfile(finPdb) and os.path.isfile(finCif):
           print("=====================================================================")
           print("|               Done                                                |")
           print("=====================================================================")


    def executeWithRDKit(self):
 
        self.printJobs()
       
        if self.useExistCoords or self.workMode==16 or self.workMode==161:
            print("One of output conformers will using input coordinates as initial ones")
        #elif self.workMode !=0 and self.workMode != 61 :
        #    print "Input coordinates will be ignored"
        print("workMode : ", self.workMode)
        # Stage 1: initiate a mol file for RDKit obj
        if self.workMode == 11 or self.workMode == 111:
            if self.monomRoot in self.chemCheck.aminoAcids:
                self.isAA = True
                self.getAAOut()
            elif os.path.isfile(self.inMmCifName) and self.chemCheck.isOrganic(self.inMmCifName, self.workMode)\
               and not self.isAA:
                tmpWorkMode = self.workMode
                self.workMode = 311
                self.runLibmol()
                self.workMode = tmpWorkMode
                # The input file is an mmcif file 
                self.fileConv.mmCifReader(self.inMmCifName)
                if len(self.fileConv.dataDescriptor):
                    self.setMonoRoot(self.fileConv.dataDescriptor)
                    self.isNA=self.checkNAFromMmcif(self.fileConv.dataDescriptor)
                    if self.isNA:
                        self.getNATors()
                    if self.monomRoot in self.chemCheck.aminoAcids:
                        self.isAA = True
                        self.getAAOut()
                if len(self.fileConv.atoms) !=0 and len(self.fileConv.bonds) !=0 and not self.isAA:
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
                    elif "props" in self.fileConv.strDescriptors \
                       and "entries" in self.fileConv.strDescriptors:

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
                            print("Smiles str  ", aSmi)
                            aIniSmiName = os.path.join(self.scrDir, self.baseRoot + "_init.smi")
                            print("Smiles file ", aIniSmiName)
                            print("Smiles str  ", aSmi)
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
                            if os.path.isfile(self.libmolAT1):
                                if not checkRepAtomTypes(self.libmolAT1, self.raS):
                                    self.tmpCifInitMol()
                                else:
                                    if os.path.isfile(self.libmolAT2):
                                        self.workMode = 900
                                        self.runLibmol()
                                        self.rdKit.reSetSmi = True
                                        self.rdKit.molecules = []
                                        self.rdKit.initMols("smi", aIniSmiName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)  
                                        self.fileConv.mergeAtomNames(self.libmolMatched, self.rdKit.molecules[0])
                                        self.fileConv.addAtomOrigChiralSign(self.rdKit.molecules[0])
                                        self.useCifCoords = True
                                    else:
                                        self.tmpCifInitMol()
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
                    """
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
                    """
                    self.fileConv.mol2Reader(self.inMol2Name)
                    self.rdKit.initMols("mol2", self.inMol2Name, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers, 0, self.fileConv.nameMapingMol2)
            else: 
                print(self.inMol2Name, " can not be found for reading ")
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
                    print("Can not add line with 'CRYST1' to the temp PDB file ", tPDBName)
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
                    print("Failed to generate initial dictionary file ", self.refmacXYZOUTName)
                    sys.exit()

        if not self.workMode in [80, 1001, 1002]:
            self.setOutCifGlobSec()
        if self.workMode in [11, 12, 13, 14, 15]:
            self.workMode = 11
        elif self.workMode in [111, 121, 131, 141, 151]:
            self.workMode = 111
        if self.workMode in [51, 52, 53, 54, 55]:
            self.workMode = 51
        if self.workMode in [11,  111, 51] and not self.isAA :
            #print len(self.rdKit.molecules)
            if len(self.rdKit.molecules):
                print("Ligand ID ", self.monomRoot)
                self.fileConv.getCCP4DataDescritor(self.rdKit.molecules[0],  self.chemCheck, self.monomRoot)
                self.rdKit.hasCCP4Type = self.fileConv.hasCCP4Type
            for iMol in range(len(self.rdKit.molecules)):
                self.inMmCifName =  os.path.join(self.scrDir, self.baseRoot + "_mol_" + str(iMol) + ".cif")
                self.initMmcifMolMap[iMol] = self.inMmCifName
                #if self.monomRoot.upper() in self.chemCheck.aminoAcids:
                #    self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[iMol], self.inMmCifName, self.chemCheck, self.monomRoot, "L-peptide")
                #else:
                if self.workMode in [11,  111]:
                    print("Using coords ", self.rdKit.useExistCoords)
                self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[iMol], self.inMmCifName, self.chemCheck, self.monomRoot, self.fileConv.chiralPre)
                if os.path.isfile(self.inMmCifName):
                    if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                        print("The input system contains metal or other heavier element")
                        print("The current version deals only with the atoms in the set of 'organic' elements") 
                        sys.exit()
                    self.runLibmol(self.inMmCifName, iMol)
                else:
                    print("The input %s does not exist"%self.inMmCifName)
                    sys.exit()

                if self.workMode == 11:
                    if  not self.runExitCode :
                        # Stage 2: optimization
                        print("===================================================================") 
                        print("| Geometrical Optimization                                        |")
                        print("===================================================================") 
                        
                        if len(self.rdKit.molecules) != 0 and os.path.isfile(self.outRstCifName):
                            inPdbNamesRoot = {} 
                            for idxMol in range(len(self.rdKit.molecules)): 
                                if idxMol not in inPdbNamesRoot:
                                    inPdbNamesRoot[idxMol] = []
                                print("Number of atoms in molecule %d is %d "%(idxMol+1, self.rdKit.molecules[idxMol].GetNumAtoms()))
                                #nConf = self.rdKit.molecules[idxMol].GetNumConformers()
                                self.runGeoOptOneMolFull(idxMol)
                                if not self.useExistCoords:
                                    self.outEnergyGeoMap(idxMol)
                    else:
                        print("Error: No dictionary cif file is generated by Acedrg ")

        elif self.workMode == 16 or self.workMode == 161:
            if os.path.isfile(self.inMmCifName):
                if not self.chemCheck.isOrganic(self.inMmCifName, self.workMode):
                    print("The input system contains metal or other heavier element")
                    print("The current version deals only with the atoms in the set of 'organic' elements") 
                    sys.exit()
                self.runLibmol()
            else:
                print("The input %s does not exist"%self.inMmCifName)
                sys.exit()
            if self.workMode == 16: 
                if not self.runExitCode :
                    # Stage 2: optimization
                    print("===================================================================") 
                    print("| Geometrical Optimization                                        |")
                    print("===================================================================") 
                        
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
                                print("Failed to produce %s after final geometrical optimization"%finPdb)

        if self.workMode in [21, 211]:
            
            # Stage 1: generate molecules and the associated bond and bond-angle values 
            # using a small molecule cif file
            if os.path.isfile(self.inStdCifName):
                self.runLibmol(self.inStdCifName)
            else:
                print("Can not find the input file ", self.inStdCifName) 
            
        if self.workMode in [22, 221]:
            
            # 1. Generate molecules using the small molecule cif files at the input directory. 
            # 2. Generate atom classes for atoms in the molecules.
            # 3. Obtain unique bond lengths and angles and cluster them according to their 
            #    component atoms in tables.
            if os.path.isdir(self.inStdCifDir):
                self.runLibmol()
            else:
                print("Can not find the input directory ", self.inStdCifDir)
                
        if self.workMode == 31 or self.workMode==311:
            print ("===================================================================") 
            print ("| Generate atom types of Acedrg style                             |")
            print ("===================================================================") 
            if os.path.isfile(self.inMmCifName):
                self.runLibmol()    
        
        if self.workMode == 32:
            print("work mode ", self.workMode)
            if os.path.isfile(self.inSmiName):
                self.runLibcheck(self.inSmiName)
            else:
                print("%s does not exist"%self.inSmiName)
                sys.exit()

            if os.path.isfile(self.libcheckLibName):
                self.runLibmol(self.libcheckLibName)
            else:
                print("%s does not exist"%self.inMmCifName)
                sys.exit()
            
        if self.workMode == 33:
            if os.path.isfile(self.inMdlName) and self.chemCheck.isOrganic(self.inMdlName, self.workMode):
                self.rdKit.initMols("mol", self.inMdlName, self.monomRoot, self.chemCheck, self.inputPara["PH"], self.numConformers)
                tMolMmcifName = os.path.join(self.scrDir, "tmpMol.cif")
                if len(self.rdKit.molecules)==1:
                    self.rdKit.MolToSimplifiedMmcif(self.rdKit.molecules[0], tMolMmcifName, self.chemCheck, self.monomRoot)
                    print("Number of molecules ", len(self.rdKit.molecules))
                    self.runLibmol(tMolMmcifName)
                else:
                    print("Do not process atom-types of multiple molecules at the same time, do nothing")

        if self.workMode == 34:
            if os.path.isfile(self.inSdfName):
                self.runLibmol()    

        if self.workMode == 35:
            if os.path.isfile(self.inStdCifName):
                self.runLibmol()    
        
        if self.workMode == 61:
            
            aAAdir      = os.path.join(self.acedrgTables, "AminoAcids")
            if not self.testMode: 
                aCLinkGenerator = CovLinkGenerator(aAAdir, self.linkInstructions, self.scrDir, self.outRoot, self.versionInfo)
            else:
                aCLinkGenerator = CovLinkGenerator(aAAdir, self.linkInstructions, self.scrDir, self.outRoot,\
                                                   self.versionInfo, self.testMode)
        if self.workMode ==111 or self.workMode ==121 or self.workMode ==131 or self.workMode ==141:
            if os.path.isfile(self.outRstCifName):
                tCif = self.outRoot + ".cif"
                #os.system("cp %s %s"%(self.outRstCifName, tCif)) 
                shutil.copy(self.outRstCifName, tCif) 
            else:
                print("acedrg failed to generate a dictionary file")     

        if self.workMode == 80:
            if os.path.isfile(self.inMmCifName):
                self.addProtCols()
            else:
                print("The mmcif file to be modified does not exist")

    def printExitInfo(self):

        print("Error: check log file at %s"%self._log_name)
    

            
# Other supplement functions

def main():
    acedrgObj = Acedrg(sys.argv[1:])

if __name__ == "__main__":
    main()

