#!/usr/bin/env ccp4-python
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

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from builtins import object
import os,os.path,sys
import platform
import glob,shutil
import re,string
import time
import math

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

if os.name != 'nt':
    import fcntl
import signal

from . exebase       import CExeCode

from . acedrgRDKit   import AcedrgRDKit
from . filetools     import FileTransformer
from . filetools     import Ccp4MmCifObj
from . chem          import ChemCheck

from . utility       import isInt
from . utility       import listComp
from . utility       import listComp2
from . utility       import listCompDes
from . utility       import listCompAcd
from . utility       import setBoolDict
from . utility       import splitLineSpa
from . utility       import BondOrderS2N
from . utility       import setNameByNumPrime

#################################################   

class CovLink(object):

    def __init__(self):
        
        #keyWordList = ["RES-NAME-1", "FILE-1"(optional), "ATOM-NAME-1",\
        #               "RES-NAME-2", "FILE-1"(optional), "ATOM-NAME-2", "BOND-TYPE",\
        #                "DELETE", ]

        self.describLevel                   = 1

        self.stdLigand1                     = {}
        self.stdLigand1["fromScr"]          = False
        self.stdLigand1["userIn"]           = False
        self.stdLigand1["compOut"]          = False
        self.stdLigand1["inCif"]            = False
        self.stdLigand1["dataBlock"]        = None
        self.stdLigand1["comp"]             = None
        self.stdLigand1["hAtom"]            = []
        self.stdLigand1["outComp"]          = False
        self.stdLigand1["remainAtoms"]      = []
        self.stdLigand1["remainIdxHs"]      = []
        self.stdLigand1["remainBonds"]      = []
        self.stdLigand1["remainAngs"]       = []
        self.stdLigand1["remainTors"]       = []
        self.stdLigand1["remainChirs"]      = []
        self.stdLigand1["linkChir"]         = []       # It will be [aChi, delAtmIdx]
        self.stdLigand1["remainPls"]        = []
        self.stdLigand2                     = {}
        self.stdLigand2["fromScr"]          = False
        self.stdLigand2["userIn"]           = False
        self.stdLigand2["compOut"]          = False
        self.stdLigand2["inCif"]            = False
        self.stdLigand2["dataBlock"]        = None
        self.stdLigand2["comp"]             = None
        self.stdLigand2["hAtom"]            = []
        self.stdLigand2["outComp"]          = False
        self.stdLigand2["remainAtoms"]      = []
        self.stdLigand2["remainIdxHs"]      = []
        self.stdLigand2["remainBonds"]      = []
        self.stdLigand2["remainAngs"]       = []
        self.stdLigand2["remainTors"]       = []
        self.stdLigand2["remainChirs"]      = []
        self.stdLigand2["linkChir"]         = []       # It will be [aChi, delAtmIdx]
        self.stdLigand2["remainPls"]        = []
        self.suggestBonds                   = []
        self.combLigand                     = {}
        self.combLigand["atoms"]            = []
        self.combLigand["hAtoms"]           = []
        self.combLigand["bonds"]            = []
        self.combLigand["chirs"]            = []
        self.outCombLigand                  = {}
        self.modLigand1                     = {}
        self.modLigand1["outComp"]            = True
        self.modLigand1["deleted"]            = {}
        self.modLigand1["deleted"]["atoms"]   = []
        self.modLigand1["deleted"]["bonds"]   = []
        self.modLigand1["deleted"]["angles"]  = []
        self.modLigand1["deleted"]["tors"]    = []
        self.modLigand1["deleted"]["chirs"]   = []
        self.modLigand1["deleted"]["planes"]  = []
        self.modLigand1["changed"]            = {}
        self.modLigand1["changed"]["atoms"]   = []
        self.modLigand1["changed"]["bonds"]   = []
        self.modLigand1["changed"]["angles"]  = []
        self.modLigand1["changed"]["tors"]    = []
        self.modLigand1["changed"]["chirs"]   = []
        self.modLigand1["changed"]["planes"]  = []
        self.modLigand1["changed"]["formal_charges"] = []
        self.modLigand1["added"]              = {}
        self.modLigand1["added"]["atoms"]     = []
        self.modLigand1["added"]["bonds"]     = []
        self.modLigand1["added"]["angles"]    = []
        self.modLigand1["added"]["tors"]      = []
        self.modLigand1["added"]["chirs"]     = []
        self.modLigand1["added"]["planes"]    = []
        self.modLigand2                       = {}
        self.modLigand2["deleted"]            = {}
        self.modLigand2["deleted"]["atoms"]   = []
        self.modLigand2["deleted"]["bonds"]   = []
        self.modLigand2["deleted"]["tors"]    = []
        self.modLigand2["deleted"]["angles"]  = []
        self.modLigand2["deleted"]["chirs"]   = []
        self.modLigand2["deleted"]["planes"]  = []
        self.modLigand2["changed"]            = {}
        self.modLigand2["changed"]["atoms"]   = []
        self.modLigand2["changed"]["bonds"]   = []
        self.modLigand2["changed"]["angles"]  = []
        self.modLigand2["changed"]["tors"]    = []
        self.modLigand2["changed"]["chirs"]   = []
        self.modLigand2["changed"]["planes"]  = []
        self.modLigand2["added"]              = {}
        self.modLigand2["added"]["atoms"]     = []
        self.modLigand2["added"]["bonds"]     = []
        self.modLigand2["added"]["angles"]    = []
        self.modLigand2["added"]["tors"]      = []
        self.modLigand2["added"]["chirs"]     = []
        self.modLigand2["added"]["planes"]    = []
        self.modLigand2["changed"]["formal_charges"] = []
        self.modLigand2["outComp"]            = True
        self.cLink                            = {}
        self.atomMap                          = {}
        self.bondMap                          = {}
        
        self.delSections                      = []

        self.errLevel                         = 0
        self.errMessage                       = []

    def checkInPara(self):
   
        aReturn = True

        if "name" not in self.stdLigand1 or "atomName" not in self.stdLigand1\
           or "name" not in self.stdLigand2 or "atomName" not in self.stdLigand2\
           or not os.path.isfile(self.stdLigand1["inCif"]) \
           or not os.path.isfile(self.stdLigand2["inCif"]):
            aReturn = False
            self.errLevel = 12
            if "name" not in self.stdLigand1:
                print("Residue 1 is not named ")
                self.errMessage.append("Residue 1 is not named\n")   
            if "atomName" not in self.stdLigand1:
                print("Linked atom in Residue 1 is not named ")
                self.errMessage.append("Linked atom in Residue 1 is not named\n")
            if "name" not in self.stdLigand2:
                print("Residue 2 is not named ")
                self.errMessage.append("Residue 2 is not named\n")
            if "atomName" not in self.stdLigand2:
                print("Linked atom in Residue 2 is not named ")
                self.errMessage.append("Linked atom in Residue 2 is not named ")
            if not os.path.isfile(self.stdLigand1["inCif"]):
                print("InCif 1", self.stdLigand1["inCif"])
                aBFName = os.path.basename(self.stdLigand1["inCif"])
                aMess = "The required file %s does not exist.\n"%aBFName
                aMess += "You need to provide this file.\n"
                print(aMess)
                self.errMessage.append(aMess)
            if not os.path.isfile(self.stdLigand2["inCif"]):
                print("InCif 2", self.stdLigand2["inCif"])
                aBFName = os.path.basename(self.stdLigand2["inCif"])
                aMess = "The required file %s does not exist.\n"%aBFName
                aMess += "You need to provide this file.\n"
                print(aMess)
                self.errMessage.append(aMess)
                #print "%s does not exist"%self.stdLigand2["inCif"]
                #self.errMessage.append("%s does not exist"%self.stdLigand2["inCif"])
        return aReturn 

    def setModiName(self):

        if self.stdLigand2["name"] != self.stdLigand1["name"]:
            self.modLigand2["name"]   = self.stdLigand2["name"] + "m1"
        else:
            self.modLigand2["name"]   = self.stdLigand2["name"] + "m2"

    def filterAtomsAndBonds(self, tFTool, tDS, tMonomer):

        for aAtom in tFTool.atoms:
            if "_chem_comp_atom.atom_id" in aAtom:
                if aAtom["_chem_comp_atom.atom_id"].upper() != tDS["ATOM-NAME"].upper():
                    tMonomer["remainAtoms"].append(aAtom)
                else:
                    # delete all bonds to the deleted atom
                    for aBond in tFTool.bonds:
                        if "_chem_comp_bond.atom_id_1" in aBond\
                           and "_chem_comp_bond.atom_id_2" in aBond:
                            if aBond["_chem_comp_bond.atom_id_1"].upper() != tDS["ATOM-NAME"].upper()\
                               and aBond["_chem_comp_bond.atom_id_2"].upper() != tDS["ATOM-NAME"].upper():
                                tMonomer["remainBonds"].append(aBond)


class CovLinkGenerator(CExeCode):

    def __init__(self, tAADir, tLinkInstructions, tScrDir, tOutRoot, tVerInfo=None, tTestMode=None):

        
        if tTestMode:
            self.testMode = True
        else:
            self.testMode = False
        self.verInfo          =  tVerInfo
        #print self.verInfo

        self.workMode         =  0   # Mode means
                                     # 0     generate a link without optimizing coordinates of 2 input monomers
                                     # 1     generate a linked system which optimizes everything. The input monomers 
                                     #       could be in form of mmcif and SMILES strings
        
        self.errMessage       = {}
        self.errLevel         = 0
        # 0 = OK  
        # 1 = No instructions file for building links
        # 2 = Instructions are not validated, no links will be built
        # 21 = "No component residue dictionary file (cif)"
        # 3 = Error in building individaul residues
        # 4 = Error in building the combined ligand
        # 5 = Error in extracting the final link info

        self.errFileName                      = ""

        # TEMPO
        if "CCP4" not in os.environ:
            print("You need to install or setup CCP4 first ")
            sys.exit() 
       
        self.allChemCombDir    = os.getenv("CLIBD_MON")
        #self.allChemCombDir   = "/Users/flong/DB/PDB_related/PDB_Ligands/Cif"

        self.aaDir            = tAADir
        self.scrDir           = tScrDir
        self.subRoot          = ""
        self.outRoot          = tOutRoot
        self.linkInstructions = tLinkInstructions
        self.ligSrcDir        = ""

        self.cLinks           = []

        # engs 

        self.chemCheck        = ChemCheck()
        self.fileTool         = FileTransformer()

        if os.path.isfile(self.linkInstructions):
            # input from a file
            self.getInstructionsForLink()
        else:
            print("%s can not be found for reading"%self.linkInstructions)
            self.errLevel = 1
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("Instruction file %s can not be found for reading"%self.linkInstructions)
        
        print("Number of links to be processed ", len(self.cLinks))

        if not self.errLevel:
            if len(self.cLinks):
                for aLink in self.cLinks:
                    self.processOneLink(aLink)
            else:
                self.errLevel  = 5
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Errors are found in the instruction file")
                self.errMessage[self.errLevel].append("Check your instruction file")

        if self.errLevel:
            print("=====================================================================")
            print("| Your job stoped because of the following errors:                  |")
            print("=====================================================================")

            self.errFileName = self.outRoot + "_errorInfo.txt"
            errF = open(self.errFileName, "w")
            for aKey in sorted(self.errMessage.keys()):
                for aL in self.errMessage[aKey]:
                    print(aL)
                    errF.write(aL)
            errF.close()
            time.sleep(0.5)
            sys.exit(self.errLevel)
            
    def getInstructionsForLink(self):

        #aFullName=os.path.basename(self.linkInstructions)
        #aExt     = ""
        #if aFullName.find(".") !=-1:
        #    aExt = aFullName.strip().split(".")[-1].strip()
        #    if aExt.find("cif") !=-1:
        #        self.getInstructionsForLinkFromCif()
        #    else:
        #        # Tempo comment off free format at the moment
        self.getInstructionsForLinkFreeFormat2()

    def checkInCompCif(self, tMonomer, tFName, tResName):
       
        if os.path.isfile(tFName):
            tMonomer["inCif"]    = tFName
            tMonomer["userIn"]   = True 
        else:
            if len(tResName) > 0:                    
                aNS = tResName[0].lower()
                aNL = tResName.upper()
                aFName = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                if os.path.isfile(aFName):
                    tMonomer["inCif"]  =  aFName
                else:
                    baseName =os.path.basename(tFName).strip()
                    self.errLevel    = 1
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Input cif %s does not exist"%baseName)
                    self.errMessage[self.errLevel].append("You need to provide an input dictionary file for %s "%aNL)
            else:
                self.errLevel    = 1
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Name for an input residue has not been given in the instruction file\n")
                                    
    def getInstructionsForLinkFromCif(self):
        
        aInsObj = Ccp4MmCifObj(self.linkInstructions)
        aLink = CovLink()
        print("Instruction file of cif format ", self.linkInstructions)
        print("Instruction keys ", list(aInsObj["ccp4CifObj"].keys()))
        if "instructs" in aInsObj["ccp4CifObj"] and "link" in aInsObj["ccp4CifObj"]["instructs"]:
            for aKey in list(aInsObj["ccp4CifObj"]["instructs"]["link"].keys()):
                if aKey =="atoms":
                    for aAtom in aInsObj["ccp4CifObj"]["instructs"]["link"]["atoms"]:
                        print("Atom in residue ", aAtom["comp_serial_num"])
                        resNum = int(aAtom["comp_serial_num"])
                        if resNum ==1:
                            aLink.stdLigand1["name"]      = aAtom["comp_id"]
                            aLink.stdLigand1["resNum"]    = resNum
                            self.checkInCompCif(aLink.stdLigand1, aAtom["comp_file"], aLink.stdLigand1["name"])    
                            aLink.stdLigand1["atomName"]  = aAtom["atom_id"]
                            aLink.modLigand1["name"]      = aLink.stdLigand1["name"] + "m1"
                        elif resNum ==2:
                            aLink.stdLigand2["name"]      = aAtom["comp_id"]
                            aLink.stdLigand2["resNum"]    = resNum
                            self.checkInCompCif(aLink.stdLigand2, aAtom["comp_file"], aLink.stdLigand2["name"])    
                            aLink.stdLigand2["inCif"]     = aAtom["comp_file"]
                            aLink.stdLigand2["atomName"]  = aAtom["atom_id"]
                            aLink.modLigand2["name"]      = aLink.stdLigand2["name"] + "m1"
                if aKey =="bonds":
                    for aBond in aInsObj["ccp4CifObj"]["instructs"]["link"]["bonds"]:
                         if aBond["type"].find(".") !=-1:
                             aBond["type"] = "single"
                         aBond["value_dist"]     = 0.0
                         aBond["value_dist_esd"] = 0.02
                         aLink.suggestBonds.append(aBond)
                if aKey =="datoms":
                    for aAtom in aInsObj["ccp4CifObj"]["instructs"]["link"]["datoms"]:
                        aDS = {}
                        aDS["atomName"]  = aAtom["atom_id"]
                        aDS["inRes"] = int(aAtom["comp_serial_num"])
                        aLink.delSections.append(aDS)

        print("Keys in the link object ", list(aInsObj["ccp4CifObj"]["instructs"]["link"].keys()))
        if aLink.checkInPara():
            self.cLinks.append(aLink)
            
    def getInstructionsForLinkFreeFormat(self):

        allLs = ""
        aList = []

        try :
            aInsF = open(self.linkInstructions, "r")
        except IOError:
            print("% can not be opened for reading ! "%self.linkInstructions)
            self.errLevel = 11
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("% can not be opened for reading ! "%self.linkInstructions)
        else:
            allLs = aInsF.readlines()
            aInsF.close()
            for aL in allLs:
                strs = aL.strip().split()
                if len(strs) !=0:
                    for aStr in strs:
                        if aStr.find("LINK:") == -1:
                            aList.append(aStr)

        if len(aList):
            lSt  = True
            lDel = False
            lCh  = False
            lAd  = False
            aDS  = {}
            i = 0
            idxAddAtm = 1
            aLink = CovLink()
           
            # Some default values. They will be changed at different stages
 
            aLink.stdLigand1["userIn"]    = False 
            aLink.stdLigand1["compOut"]   = False 
            aLink.stdLigand1["outComp"]   = False 
            aLink.stdLigand1["outMod"]    = True 
            aLink.stdLigand2["userIn"]    = False 
            aLink.stdLigand2["compOut"]   = False 
            aLink.stdLigand2["outComp"]   = False 
            aLink.stdLigand2["outMod"]    = True 

            aLink.stdLigand1["group"]     = "."
            aLink.stdLigand2["group"]     = "."

            nL = len(aList)
            while i < nL :
                #print aList[i]
                if lSt  :
                    if (i+1) < nL:
                        if aList[i].upper().find("RES-NAME-1") != -1:
                            aLink.stdLigand1["name"]   = aList[i+1]
                            aLink.stdLigand1["resNum"] = 1
                            aLink.modLigand1["name"]   = aLink.stdLigand1["name"] + "m1"
                            i +=2
                        elif aList[i].upper().find("FILE-1") != -1:
                            aLink.stdLigand1["inCif"]    = aList[i+1]
                            aLink.stdLigand1["compOut"]  =  True
                            i +=2
                        elif aList[i].upper().find("ATOM-NAME-1") != -1:
                            atm1Name = setNameByNumPrime(aList[i+1])
                            aLink.stdLigand1["atomName"] = atm1Name
                            i +=2
                        elif aList[i].upper().find("RES-NAME-2") != -1:
                            aLink.stdLigand2["name"] = aList[i+1]
                            aLink.stdLigand2["resNum"] = 2
                            if aLink.stdLigand2["name"] != aLink.stdLigand1["name"]:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m1"
                            else:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m2"
                            i +=2
                        elif aList[i].upper().find("FILE-2") != -1:
                            aLink.stdLigand2["inCif"] = aList[i+1]
                            aLink.stdLigand2["compOut"]  =  True
                            i +=2
                        elif aList[i].upper().find("ATOM-NAME-2") != -1:
                            atm2Name = setNameByNumPrime(aList[i+1])
                            aLink.stdLigand2["atomName"] =  atm2Name
                            i +=2
                        elif aList[i].upper().find("BOND-TYPE") != -1 :
                            aBond = {}
                            aBond["type"]           = "single"
                            aBond["type"]           = aList[i+1]
                            aBond["atom_id_1"]      = aLink.stdLigand1["atomName"]
                            aBond["comp_serial_num_1"] =  1 
                            aBond["atom_id_2"]      = aLink.stdLigand2["atomName"]
                            aBond["comp_serial_num_2"] =  2 
                            aBond["value_dist"]     = 0.0     
                            aBond["value_dist_esd"] = 0.02
                            aLink.suggestBonds.append(aBond) 
                            i+=2
                        elif aList[i].upper().find("DELETE") != -1:
                            lDel = True
                            lSt  = False
                            lCh  = False
                            lAd  = False
                            i +=1
                        elif aList[i].upper().find("CHANGE") != -1:
                            lCh  = True
                            lSt  = False
                            lDel = False
                            lAd  = False
                            i +=1
                        elif aList[i].upper().find("ADD") != -1:
                            lAd  = True
                            lSt  = False
                            lDel = False
                            lCh  = False
                            i +=1
                    else:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        aMess = "Format error in the instruction file. "
                        aMess+= "%s needs to be followed by a value.\n"% aList[i] 
                        self.errMessage[self.errLevel].append(aMess) 
                        break
                elif lDel :
                    if len(list(aDS.keys())) != 0:
                        aLink.delSections.append(aDS)
                        aDS = {}
                    if aList[i].upper().find("ATOM") != -1 and aList[i].upper().find("ATOM-NAME") == -1:
                        if i+2 < nL:
                            atmDName = setNameByNumPrime(aList[i+1])
                            aDS["atomName"]  = atmDName
                            aNum = int(aList[i+2])
                            aDS["inRes"]     = aNum
                            aAtom = {}
                            aAtom["atom_id"] = aList[i+1]
                            if aNum == 1:
                                aLink.modLigand1["deleted"]["atoms"].append(aAtom) 
                            elif aNum == 2:
                                aLink.modLigand2["deleted"]["atoms"].append(aAtom) 
                            else:
                                self.errLevel = 12
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE section of the instruction file. "
                                aMess+= "Residue number should be 1 or 2"
                                self.errMessage[self.errLevel].append(aMess)
                            i+=3
                        else:
                            self.errLevel = 12
                            if self.errMessage not in self.errMessage:
                                self.errMessage[self.errLevel] = []
                            aMess = "Format error in DELECTE ATOM section of the instruction file. "
                            aMess+= "%s needs to be followed by two values.\n"% aList[i]
                            self.errMessage[self.errLevel].append(aMess)
                            break
                    elif aList[i].upper().find("BOND") != -1 :
                        if i+3 < nL:
                            aBond = {}
                            aBond["atom_id_1"] = aList[i+1]
                            aBond["atom_id_2"] = aList[i+2]
                            aNum  = int(aList[i+3])
                            aBond["comp_serial_num_1"] = aNum 
                            aBond["comp_serial_num_2"] = aNum 
                            if aNum == 1:
                                aLink.modLigand1["deleted"]["bonds"].append(aBond)
                            elif aNum == 2:
                                aLink.modLigand2["deleted"]["bonds"].append(aBond)
                            else:
                                self.errLevel = 12
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE BOND section of the instruction file. "
                                aMess+= "Residue number should be 1 or 2"
                                self.errMessage[12].append(aMess)
                            i+=4
                        else:
                            self.errLevel = 12
                            if self.errLevel not in self.errMessage:
                                self.errMessage[self.errLevel] = []
                            aMess = "Format error in DELECTE BOND section of the instruction file. "
                            aMess+= "%s needs to be followed by two values.\n"% aList[i]
                            self.errMessage[self.errLevel].append(aMess)
                            break
                    elif aList[i].upper().find("CHANGE") != -1:
                        lCh  = True
                        lSt  = False
                        lDel = False
                        lAd  = False
                        i +=1
                    elif aList[i].upper().find("ADD") != -1:
                        lAd  = True
                        lSt  = False
                        lDel = False
                        lCh  = False
                        i +=1
                    else:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        aMess = "Format error in DELECTE section of the instruction file. "
                        aMess+= "Unknown key or value entry %s.\n"% aList[i]
                        self.errMessage[self.errLevel].append(aMess)
                        break
                elif lCh :
                    if aList[i].upper().find("BOND") != -1 :
                        if i+4 < nL:
                            aBond = {}
                            aBond["atom_id_1"]      = aList[i+1]
                            aBond["atom_id_2"]      = aList[i+2]
                            aBond["type"]           = aList[i+3]
                            aBond["comp_serial_num"]= int(aList[i+4])
                            aBond["value_dist"]     = 0.0     
                            aBond["value_dist_esd"] = 0.02
                            i+=5
                            if aBond["comp_serial_num"] ==1 :
                               aLink.modLigand1["changed"]["bonds"].append(aBond)
                            elif aBond["comp_serial_num"] ==2 :
                               aLink.modLigand2["changed"]["bonds"].append(aBond)
                            else:
                                print("Error in BOND CHANGE section: invalide monomer serial number (should be 1 or 2)")
                                print("Check your link instruction file")
                                self.errLevel = 1
                                break
                        else:
                            print("Error in CHANGE BOND section: check your link instruction file")
                            print(allLs)
                            self.errLevel = 2
                            break
                    elif aList[i].upper().find("DELETE") != -1:
                        lDel  = True
                        lSt   = False
                        lCh   = False
                        lAd   = False
                    elif aList[i].upper().find("ADD") != -1:
                        lAd  = True
                        lSt  = False
                        lDel = False
                        lCh  = False
                    else:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        aMess = "Format error in CHANGE section of the instruction file. "
                        aMess+= "Unknown key or value entry %s.\n"% aList[i]
                        self.errMessage[self.errLevel].append(aMess)
                        break
                elif lAd :
                    # Leave for discussion
                    if aList[i].upper().find("ATOM") != -1 :
                        if i+5 < nL:
                            aAtom = {}
                            aAtom["atom_id"]     = aList[i+1]
                            aAtom["type_symbol"] = aList[i+2]
                            aAtom["type_energy"] = aAtom["type_symbol"]
                            aNum = int(aList[i+3]) 
                            aAtom["res_idx"] =  aNum
                            aBond = {}
                            aBond["atom_id_1"] = aAtom["atom_id"]
                            aBond["atom_id_2"] = aList[i+4]
                            aBond["type"]      = aList[i+5]
                            i+=6
                            if aNum == 1:
                                aAtom["comp_id"] = aLink.stdLigand1["name"]
                                aLink.modLigand1["added"]["atoms"].append(aAtom)
                                aLink.modLigand1["added"]["bonds"].append(aBond)
                            elif aNum==2:
                                aAtom["comp_id"] = aLink.stdLigand2["name"]
                                aLink.modLigand2["added"]["atoms"].append(aAtom)
                                aLink.modLigand2["added"]["bonds"].append(aBond)
                            else:
                                print("Error in ADD ATOM section: Residue number should be 1 or 2 ! ")
                                print("Check your link instruction file")
                                print(allLs)
                                self.errLevel = 2
                        else:
                            print("Error in ADD ATOM section: Not enough infomation ")
                            print("Check your link instruction file")
                            print(allLs)
                            self.errLevel = 2
                    elif aList[i].upper().find("DELETE") != -1:
                        lDel  = True
                        lSt   = False
                        lCh   = False
                        lAd   = False
                    elif aList[i].upper().find("CHANGE") != -1:
                        lAd  = True
                        lSt  = False
                        lDel = False
                        lCh  = False
 
            if len(list(aDS.keys())) != 0:
                aLink.delSections.append(aDS)
            
            if len(aLink.suggestBonds)==0 and "atomName" in aLink.stdLigand1 and "atomName" in aLink.stdLigand2:
                aBond = {}
                aBond["type"]              = "single"
                aBond["atom_id_1"]         = aLink.stdLigand1["atomName"]
                aBond["comp_serial_num_1"] =  1 
                aBond["atom_id_2"]         = aLink.stdLigand2["atomName"]
                aBond["comp_serial_num_2"] =  2 
                aBond["value_dist"]        = 0.0     
                aBond["value_dist_esd"]    = 0.02
                aLink.suggestBonds.append(aBond) 

            aNS = aLink.stdLigand1["name"][0].lower()
            aNL = aLink.stdLigand1["name"].upper()
            if not aLink.stdLigand1["userIn"]:
                aLink.stdLigand1["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                aLink.stdLigand1["userIn"] = True
                #if aNL in self.chemCheck.aminoAcids :
                #    aLink.stdLigand1["inCif"] = os.path.join(self.aaDir, aNL + ".cif")

            aNS = aLink.stdLigand2["name"][0].lower()
            aNL = aLink.stdLigand2["name"].upper()
            if not aLink.stdLigand2["userIn"]:
                aLink.stdLigand2["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                aLink.stdLigand2["userIn"] = True
                #if aNL in self.chemCheck.aminoAcids :
                #    aLink.stdLigand2["inCif"] = os.path.join(self.aaDir, aNL + ".cif")

            if aLink.checkInPara() and not self.errLevel:
                self.cLinks.append(aLink)
                print("Instructions for build a link are  ")
                print("Link will happen between ")
                print("Atom : %s in Residue %s "%(aLink.stdLigand1["atomName"], aLink.stdLigand1["name"]))
                print(" and ") 
                print("Atom : %s in Residue %s "%(aLink.stdLigand2["atomName"], aLink.stdLigand2["name"]))
                if len(aLink.suggestBonds) >0:
                    print("The suggested bonds linking two residues are: ")
                    for aBond in aLink.suggestBonds:
                        print("Bond between atom %s in residue %s and atom %s in residue %s "\
                              %(aBond["atom_id_1"], aBond["comp_serial_num_1"],\
                                aBond["atom_id_2"], aBond["comp_serial_num_2"]))
                        print("Bond order is %s "%aBond["type"])
                       
                print("Two input cif files are : ")
                print("%s for comp 1 and "%aLink.stdLigand1["inCif"])
                print("%s for comp 2"%aLink.stdLigand2["inCif"])

                #if len(aLink.delSections):
                nda1 = len(aLink.modLigand1["deleted"]["atoms"])
                nda2 = len(aLink.modLigand2["deleted"]["atoms"])
                if nda1 >0 or nda2 >0:
                    print("The following atoms are deleted.")
                    if nda1 >0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aA in aLink.modLigand1["deleted"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])
                    if nda2 >0:
                        print("In residue %s: "%aLink.modLigand2["name"])
                        for aA in aLink.modLigand2["deleted"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])

                    #for aDS in aLink.delSections:
                    #    if aDS.has_key("inRes") and aDS.has_key("atomName"):
                    #        print "Atom %s in Residue %d "%(aDS["atomName"], aDS["inRes"]) 

                ndb1 = len(aLink.modLigand1["deleted"]["bonds"])
                ndb2 = len(aLink.modLigand2["deleted"]["bonds"])
                if ndb1 >0 or ndb2 >0:
                    print("The following bonds are deleted.")
                    if ndb1 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aB in aLink.modLigand1["deleted"]["bonds"]:
                            print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                    if ndb2 > 0:
                        print("In residue %s: "%aLink.modLigand2["name"])
                        for aB in aLink.modLigand2["deleted"]["bonds"]:
                            print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))

                naa1 = len(aLink.modLigand1["added"]["atoms"])
                naa2 = len(aLink.modLigand2["added"]["atoms"])
                if naa1 > 0 or naa2 > 0:
                    print("The following atoms are added.")
                    if naa1 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aA in aLink.modLigand1["added"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])
                    if naa2 > 0:
                        print("In residue %s: "%aLink.modLigand2["name"])
                        for aA in aLink.modLigand2["added"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])

                nab1 = len(aLink.modLigand1["added"]["bonds"])
                nab2 = len(aLink.modLigand2["added"]["bonds"])
                if nab1 >0 or nab2 >0:
                    print("The following bonds are added.")
                    if nab1 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aB in aLink.modLigand1["added"]["bonds"]:
                            print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                    if nab2 > 0:
                        print("In residue %s: "%aLink.modLigand2["name"])
                        for aB in aLink.modLigand2["added"]["bonds"]:
                            print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                        
                nca1 = len(aLink.modLigand1["changed"]["atoms"])
                nca2 = len(aLink.modLigand2["changed"]["atoms"])
                if nca1 > 0 or nca2 > 0:
                    print("The following atoms are changed.")
                    if nca1 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aA in aLink.modLigand1["changed"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])
                    if nca2 > 0:
                        print("In residue %s: "%aLink.modLigand2["name"])
                        for aA in aLink.modLigand2["changed"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])

                ncb1 = len(aLink.modLigand1["changed"]["bonds"])
                ncb2 = len(aLink.modLigand2["changed"]["bonds"])
                if ncb1 >0 or ncb2 >0:
                    print("The following bonds are changed.")
                    if ncb1 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aB in aLink.modLigand1["changed"]["bonds"]:
                            print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                            print("The bond-order is now ", aB["type"])
                    if ncb2 > 0:
                        print("In residue %s: "%aLink.modLigand1["name"])
                        for aB in aLink.modLigand2["changed"]["bonds"]:
                            print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                            print("The bond-order is now ", aB["type"])
            else:
                self.errLevel = aLink.errLevel
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                for aL in aLink.errMessage:
                    self.errMessage[self.errLevel].append(aL)
                self.errMessage[self.errLevel].append("Information in the instruction file is not correct/enough to build a link. \n")
        else:
            self.errLevel = 12
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("The input instruction file is empty\n")


    def getInstructionsForLinkFreeFormat2(self):

        allLs = ""
        aList = []

        try :
            aInsF = open(self.linkInstructions, "r")
        except IOError:
            print("% can not be opened for reading ! "%self.linkInstructions)
            self.errLevel = 11
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("% can not be opened for reading ! "%self.linkInstructions)
        else:
            allLs = aInsF.readlines()
            aInsF.close()
            for aL in allLs:
                strs = aL.strip().split()
                if len(strs) !=0:
                    for aStr in strs:
                        if aStr.find("LINK:") == -1:
                            aList.append(aStr)

        if len(aList):
            lSt  = True
            lDel = False
            lCh  = False
            lAd  = False
            aDS  = {}
            i = 0
            idxAddAtm = 1
            aLink = CovLink()
           
            # Some default values. They will be changed at different stages
 
            aLink.stdLigand1["userIn"]  = False 
            aLink.stdLigand1["outComp"] = False 
            aLink.stdLigand1["outMod"]  = True 
            aLink.stdLigand2["userIn"]  = False 
            aLink.stdLigand2["outComp"] = False 
            aLink.stdLigand2["outMod"]  = True 

            aLink.stdLigand1["group"]   = "."
            aLink.stdLigand2["group"]   = "."

            nL = len(aList)

            while i < nL :
                #print aList[i]
                if (i+1) < nL:
                    if aList[i].upper().find("RES-NAME-1") != -1:
                        aLink.stdLigand1["name"]   = aList[i+1]
                        aLink.stdLigand1["resNum"] = 1
                        aLink.modLigand1["name"]   = aLink.stdLigand1["name"] + "m1"
                        if "name" in aLink.stdLigand2:
                            if aLink.stdLigand2["name"] != aLink.stdLigand1["name"]:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m1"
                            else:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m2"
                        i +=2
                    elif aList[i].upper().find("FILE-1") != -1:
                        aLink.stdLigand1["inCif"] = aList[i+1]
                        aLink.stdLigand1["userIn"]   =  True
                        aLink.stdLigand1["compOut"]  =  True
                        i +=2
                    elif aList[i].upper().find("ATOM-NAME-1") != -1:
                        atm1Name = setNameByNumPrime(aList[i+1])
                        aLink.stdLigand1["atomName"] =  atm1Name
                        #aLink.stdLigand1["atomName"] = aList[i+1] 
                        i +=2
                    elif aList[i].upper().find("RES-NAME-2") != -1:
                        aLink.stdLigand2["name"] = aList[i+1]
                        aLink.stdLigand2["resNum"] = 2
                        if "name" in aLink.stdLigand1:
                            if aLink.stdLigand2["name"] != aLink.stdLigand1["name"]:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m1"
                            else:
                                aLink.modLigand2["name"]   = aLink.stdLigand2["name"] + "m2"
                        i +=2
                    elif aList[i].upper().find("FILE-2") != -1:
                        aLink.stdLigand2["inCif"] = aList[i+1]
                        aLink.stdLigand2["userIn"]   =  True
                        aLink.stdLigand2["compOut"]  =  True
                        i +=2
                    elif aList[i].upper().find("ATOM-NAME-2") != -1:
                        atm2Name = setNameByNumPrime(aList[i+1])
                        aLink.stdLigand2["atomName"] =  atm2Name
                        #aLink.stdLigand2["atomName"] = aList[i+1] 
                        i +=2
                    elif aList[i].upper().find("BOND-TYPE") != -1 :
                        aBond = {}
                        aBond["type"]           = "single"
                        aBond["type"]           = aList[i+1]
                        aBond["atom_id_1"]      = aLink.stdLigand1["atomName"]
                        aBond["comp_serial_num_1"] =  1 
                        aBond["atom_id_2"]      = aLink.stdLigand2["atomName"]
                        aBond["comp_serial_num_2"] =  2 
                        aBond["value_dist"]     = 0.0     
                        aBond["value_dist_esd"] = 0.02
                        aLink.suggestBonds.append(aBond) 
                        i+=2
                    elif aList[i].upper().find("DELETE") != -1:
                        lDel = True
                        lSt  = False
                        lCh  = False
                        lAd  = False
                        i +=1
                    elif aList[i].upper().find("CHANGE") != -1:
                        lCh  = True
                        lSt  = False
                        lDel = False
                        lAd  = False
                        i +=1
                    elif aList[i].upper().find("ADD") != -1:
                        lAd  = True
                        lSt  = False
                        lDel = False
                        lCh  = False
                        i +=1
                    elif lDel :
                        lDel = False
                        if len(list(aDS.keys())) != 0:
                            aLink.delSections.append(aDS)
                            aDS = {}
                        if aList[i].upper().find("ATOM") != -1 and aList[i].upper().find("ATOM-NAME") == -1:
                            if i+2 < nL:
                                if not aList[i+2].isdigit():
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Format error in DELECTE ATOM section of the instruction file. "
                                    aMess += "%s should be followed by residue number 1 or 2"%aList[i+1]
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                                atmDName = setNameByNumPrime(aList[i+1])
                                aDS["atomName"]  = atmDName
                                aNum = int(aList[i+2])
                                aDS["inRes"]     = aNum
                                aAtom = {}
                                aAtom["atom_id"] = atmDName
                                if aNum == 1:
                                    aLink.modLigand1["deleted"]["atoms"].append(aAtom) 
                                elif aNum == 2:
                                    aLink.modLigand2["deleted"]["atoms"].append(aAtom) 
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Format error in DELECTE section of the instruction file. "
                                    aMess+= "Residue number should be 1 or 2"
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                                i+=3
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE ATOM section of the instruction file. "
                                aMess+= "%s needs to be followed by two values.\n"% aList[i]
                                self.errMessage[self.errLevel].append(aMess)
                                break
                        elif aList[i].upper().find("BOND") != -1 :
                            if i+3 < nL:
                                aBond = {}
                                aBond["atom_id_1"] = aList[i+1]
                                aBond["atom_id_2"] = aList[i+2]
                                aNum  = int(aList[i+3])
                                aBond["comp_serial_num_1"] = aNum 
                                aBond["comp_serial_num_2"] = aNum 
                                if aNum == 1:
                                    aLink.modLigand1["deleted"]["bonds"].append(aBond)
                                elif aNum == 2:
                                    aLink.modLigand2["deleted"]["bonds"].append(aBond)
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Format error in DELECTE BOND section of the instruction file. "
                                    aMess+= "Residue number should be 1 or 2"
                                    self.errMessage[12].append(aMess)
                                    break
                                i+=4
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE BOND section of the instruction file.\n"
                                aMess+= "%s needs to be followed by two values.\n"% aList[i]
                                self.errMessage[self.errLevel].append(aMess)
                                break
                    elif lCh :
                        if aList[i].upper().find("BOND") != -1 :
                            if i+4 < nL:
                                aBond = {}
                                aBond["atom_id_1"]      = aList[i+1]
                                aBond["atom_id_2"]      = aList[i+2]
                                aBond["type"]           = aList[i+3]
                                aBond["comp_serial_num"]= int(aList[i+4])
                                aBond["value_dist"]     = 0.0     
                                aBond["value_dist_esd"] = 0.02
                                i+=5
                                if aBond["comp_serial_num"] ==1 :
                                    aLink.modLigand1["changed"]["bonds"].append(aBond)
                                elif aBond["comp_serial_num"] ==2 :
                                    aLink.modLigand2["changed"]["bonds"].append(aBond)
                                else:
                                    print("Error in BOND CHANGE section: invalide monomer serial number (should be 1 or 2)")
                                    print("Check your link instruction file")
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in BOND CHANGE section in the instruction file.\n"
                                    aMess+= "invalide monomer serial number (should be 1 or 2)\n"
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in BOND CHANGE section in the instruction file:\n "
                                aMess+= "Not enough information to define the change of the bond\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break

                        elif aList[i].upper().find("CHARGE") != -1 :
                            if i+3 < nL:
                                aAddCharge = {}
                                aList[i+1] = aList[i+1].strip()
                                if aList[i+1] in ["1", "2"]:
                                    aAddCharge["comp_serial_num"] = int(aList[i+1])
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in add CHARGE section in the instruction file:\n "
                                    aMess+= "Residue number must be 1 or 2.\n"
                                    aMess+= "Your input is %s.\n"%aList[i+1]
                                    self.errMessage[self.errLevel].append(aMess)
                                    break

                                aAddCharge["atom_id"]       = aList[i+2].strip()

                                aList[i+3] = aList[i+3].strip()
                                if isInt(aList[i+3]):
                                    aAddCharge["formal_charge"] = int(aList[i+3])
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in add CHARGE section in the instruction file:\n "
                                    aMess+= "Formal charge added must be an integer.\n"
                                    aMess+= "Your input is %s.\n"%aList[i+3]
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                                i+=4
                                if aAddCharge["comp_serial_num"] ==1 :
                                    aLink.modLigand1["changed"]["formal_charges"].append(aAddCharge)
                                elif aAddCharge["comp_serial_num"] ==2 :
                                    aLink.modLigand2["changed"]["formal_charges"].append(aAddCharge)
                            else :
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in the formal charge CHANGE section in the instruction file:\n "
                                aMess+= "Not enough information to define the change of the formal charge\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break


                    elif lAd :
                        # Leave for discussion
                        if aList[i].upper().find("ATOM") != -1 :
                            if i+4 < nL:
                                aAtom = {}
                                aAtom["atom_id"]     = aList[i+1]
                                aAtom["type_symbol"] = aList[i+2]
                                aAtom["charge"]      = aList[i+3]
                                aAtom["type_energy"] = aAtom["type_symbol"]
                                aAtom["is_added"]    = True
                                aNum = int(aList[i+4]) 
                                if aNum == 1:
                                    aAtom["comp_id"] = aLink.stdLigand1["name"]
                                    aLink.modLigand1["added"]["atoms"].append(aAtom)
                                elif aNum==2:
                                    aAtom["comp_id"] = aLink.stdLigand2["name"]
                                    aLink.modLigand2["added"]["atoms"].append(aAtom)
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in ADD ATOM section:\n "
                                    aMess += "Residue number for an atom should be 1 or 2 !\n "
                                    aMess += "Check your link instruction file.\n"
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                                i+=5
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in add ATOM section in the instruction file:\n "
                                aMess+= "Not enough information to define the added atom\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break
                        if aList[i].upper().find("BOND") != -1 :
                            if i+4 < nL:
                                aBond = {}
                                aBond["atom_id_1"] = aList[i+1]
                                aBond["atom_id_2"] = aList[i+2]
                                aBond["type"]      = aList[i+3]
                                aNum = int(aList[i+4]) 
                                aBond["comp_serial_num"]= aNum
                                aBond["value_dist"]     = 0.0     
                                aBond["value_dist_esd"] = 0.02
                                if aNum == 1:
                                    aLink.modLigand1["added"]["bonds"].append(aBond)
                                elif aNum==2:
                                    aLink.modLigand2["added"]["bonds"].append(aBond)
                                i+=5
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in add BOND section in the instruction file:\n "
                                aMess+= "Not enough information to define the added bond\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break
                    else :
                        self.errLevel = 2
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        aMess = "Unknown keyword  %s. Check your instruction file\n"%aList[i]
                        self.errMessage[self.errLevel].append(aMess)  
                        break
                else :
                    print("Errors begin at %s. Check your instruction file "%aList[i])
                    self.errLevel = 2
                    if 2 not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Errors begin at %s. Check your instruction file\n"%aList[i])             
                    break
            if self.errLevel==0:
                if len(list(aDS.keys())) != 0:
                    aLink.delSections.append(aDS)
            
                if len(aLink.suggestBonds)==0 and "atomName" in aLink.stdLigand1 and "atomName" in aLink.stdLigand2:
                    aBond = {}
                    aBond["type"]              = "single"
                    aBond["atom_id_1"]         = aLink.stdLigand1["atomName"]
                    aBond["comp_serial_num_1"] =  1 
                    aBond["atom_id_2"]         = aLink.stdLigand2["atomName"]
                    aBond["comp_serial_num_2"] =  2 
                    aBond["value_dist"]        = 0.0     
                    aBond["value_dist_esd"]    = 0.02
                    aLink.suggestBonds.append(aBond) 

                aNS = aLink.stdLigand1["name"][0].lower()
                aNL = aLink.stdLigand1["name"].upper()
                if not aLink.stdLigand1["userIn"]:
                    aLink.stdLigand1["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                    aLink.stdLigand1["userIn"] = True

                aNS = aLink.stdLigand2["name"][0].lower()
                aNL = aLink.stdLigand2["name"].upper()
                if not aLink.stdLigand2["userIn"]:
                    aLink.stdLigand2["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                    aLink.stdLigand2["userIn"] = True

                if aLink.checkInPara() and not self.errLevel:
                    self.cLinks.append(aLink)
                    print("Instructions for build a link are  ")
                    print("Link will happen between ")
                    print("Atom : %s in Residue %s "%(aLink.stdLigand1["atomName"], aLink.stdLigand1["name"]))
                    print(" and ") 
                    print("Atom : %s in Residue %s "%(aLink.stdLigand2["atomName"], aLink.stdLigand2["name"]))
                    if len(aLink.suggestBonds) >0:
                        print("The suggested bonds linking two residues are: ")
                        for aBond in aLink.suggestBonds:
                            print("Bond between atom %s in residue %s and atom %s in residue %s "\
                                  %(aBond["atom_id_1"], aBond["comp_serial_num_1"],\
                                    aBond["atom_id_2"], aBond["comp_serial_num_2"]))
                            print("Bond order is %s "%aBond["type"])

                    print("Two input cif files are : ")
                    print("%s for comp 1 and "%aLink.stdLigand1["inCif"])
                    print("%s for comp 2"%aLink.stdLigand2["inCif"])

                    #if len(aLink.delSections):
                    nda1 = len(aLink.modLigand1["deleted"]["atoms"])
                    nda2 = len(aLink.modLigand2["deleted"]["atoms"])
                    if nda1 >0 or nda2 >0:
                        print("The following atoms are deleted.")
                        if nda1 >0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aA in aLink.modLigand1["deleted"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])
                        if nda2 >0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aA in aLink.modLigand2["deleted"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])

                        #for aDS in aLink.delSections:
                        #    if aDS.has_key("inRes") and aDS.has_key("atomName"):
                        #        print "Atom %s in Residue %d "%(aDS["atomName"], aDS["inRes"]) 

                    ndb1 = len(aLink.modLigand1["deleted"]["bonds"])
                    ndb2 = len(aLink.modLigand2["deleted"]["bonds"])
                    if ndb1 >0 or ndb2 >0:
                        print("The following bonds are deleted.")
                        if ndb1 > 0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aB in aLink.modLigand1["deleted"]["bonds"]:
                                print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                        if ndb2 > 0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aB in aLink.modLigand2["deleted"]["bonds"]:
                                print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))


                    naa1 = len(aLink.modLigand1["added"]["atoms"])
                    naa2 = len(aLink.modLigand2["added"]["atoms"])
                    if naa1 > 0 or naa2 > 0:
                        print("The following atoms are added.")
                        if naa1 > 0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aA in aLink.modLigand1["added"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])
                        if naa2 > 0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aA in aLink.modLigand2["added"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])

                    nab1 = len(aLink.modLigand1["added"]["bonds"])
                    nab2 = len(aLink.modLigand2["added"]["bonds"])
                    if nab1 >0 or nab2 >0:
                        print("The following bonds are added.")
                        if nab1 > 0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aB in aLink.modLigand1["added"]["bonds"]:
                                print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                        if nab2 > 0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aB in aLink.modLigand2["added"]["bonds"]:
                                print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                        
                    nca1 = len(aLink.modLigand1["changed"]["atoms"])
                    ncc1 = len(aLink.modLigand1["changed"]["formal_charges"])
                    nca2 = len(aLink.modLigand2["changed"]["atoms"])
                    ncc2 = len(aLink.modLigand2["changed"]["formal_charges"])
                    if nca1 > 0 or ncc1 > 0 or nca2 > 0 or ncc2 > 0:
                        print("The following atoms are changed.")
                        if nca1 > 0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aA in aLink.modLigand1["changed"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])
                        if ncc1 > 0:
                            print("In residue %s: "%aLink.modLigand1["name"])
                            for aA in aLink.modLigand1["changed"]["formal_charges"]:
                                print("Atom %s "%aA["atom_id"])
                                print("Formal charge %d  "%aA["formal_charge"])
                        if nca2 > 0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aA in aLink.modLigand2["changed"]["atoms"]:
                                print("Atom %s "%aA["atom_id"])
                        if ncc2 > 0:
                            print("In residue %s: "%aLink.modLigand2["name"])
                            for aA in aLink.modLigand2["changed"]["formal_charges"]:
                                print("Atom %s "%aA["atom_id"])
                                print("Formal charge %d  "%aA["formal_charge"])


                    ncb1 = len(aLink.modLigand1["changed"]["bonds"])
                    ncb2 = len(aLink.modLigand2["changed"]["bonds"])
                    if ncb1 >0 or ncb2 >0:
                        print("The following bonds are changed.")
                        if ncb1 > 0:
                            for aB in aLink.modLigand1["changed"]["bonds"]:
                                print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                                print("The bond-order is now ", aB["type"])
                            print("In residue %s: "%aLink.modLigand1["name"])
                        if ncb2 > 0:
                            for aB in aLink.modLigand2["changed"]["bonds"]:
                                print("Bond between atoms %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                                print("The bond-order is now ", aB["type"])
                            print("In residue %s: "%aLink.modLigand1["name"])
                else:
                    self.errLevel = aLink.errLevel
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    for aL in aLink.errMessage:
                        self.errMessage[self.errLevel].append(aL)
                    self.errMessage[self.errLevel].append("Information in the instruction file is not correct/enough to build a link. \n")
            else:
                for aKey in list(self.errMessage.keys()):
                    for aLine in self.errMessage[aKey]:
                        print(aLine)
	#print "User input cif for L 1 ", aLink.stdLigand1["compOut"]
        #print "User input cif for L 2 ", aLink.stdLigand2["compOut"]

    def processOneLink(self, tLinkIns):
        if not self.errLevel:
            if not tLinkIns.stdLigand1["fromScr"]:
                self.setOneCompFromCif(tLinkIns.stdLigand1["inCif"], tLinkIns.stdLigand1)
                if not self.errLevel: 
                    print("Comp 1 ", tLinkIns.stdLigand1["name"], " contains ") 
                    if "atoms" in tLinkIns.stdLigand1["comp"]:
                        print(len(tLinkIns.stdLigand1["comp"]["atoms"]), " atoms.")
                    if "bonds" in tLinkIns.stdLigand1["comp"]:
                        print(len(tLinkIns.stdLigand1["comp"]["bonds"]), " bonds.")
                    if "chirs" in tLinkIns.stdLigand1["comp"]:
                        print(len(tLinkIns.stdLigand1["comp"]["chirs"]), " chiral centres.")
                    #tLinkIns.stdLigand1["ccp4MmCifObj"].printOneComp(tLinkIns.stdLigand1["name"])
                    if not self.errLevel:
                        self.selectHAtoms(tLinkIns.stdLigand1["comp"])
                        print("Number of H atoms in residue %s is %d "%(tLinkIns.stdLigand1["name"], len(tLinkIns.stdLigand1["comp"]["hAtoms"])))
            else:
                self.setOneMonomer(tLinkIns.stdLigand1)
                #print "output comp 1 ", tLinkIns.stdLigand1["outComp"] 
           
            if not self.errLevel:
                if not tLinkIns.stdLigand2["fromScr"]:
                    self.setOneCompFromCif(tLinkIns.stdLigand2["inCif"], tLinkIns.stdLigand2) 
                    if not self.errLevel:   
                        print("Comp 2 ", tLinkIns.stdLigand2["name"], " contains ") 
                        if "atoms" in tLinkIns.stdLigand2["comp"]:
                            print(len(tLinkIns.stdLigand2["comp"]["atoms"]), " atoms.")
                        if "bonds" in tLinkIns.stdLigand2["comp"]:
                            print(len(tLinkIns.stdLigand2["comp"]["bonds"]), " bonds.")
                        if "chirs" in tLinkIns.stdLigand2["comp"]:
                            print(len(tLinkIns.stdLigand2["comp"]["chirs"]), " chiral centres.")
                        #tLinkIns.stdLigand2["ccp4MmCifObj"].printOneComp(tLinkIns.stdLigand2["name"])
                        self.selectHAtoms(tLinkIns.stdLigand2["comp"])
                        print("Number of H atoms in residue %s is %d "%(tLinkIns.stdLigand2["name"], len(tLinkIns.stdLigand2["comp"]["hAtoms"])))
                else:
                    self.setOneMonomer(tLinkIns.stdLigand2)
                    #print "output comp 2 ", tLinkIns.stdLigand2["outComp"]

            if not self.errLevel:
                print("##################################################################")
                print("#                                                                #")
                print("#  Build a combo-ligand combining both monomers via the link     #")    
                print("#                                                                #")
                print("##################################################################")
                tLinkIns.cLink["name"] = tLinkIns.stdLigand1["name"] + "-" + tLinkIns.stdLigand2["name"]
                self.buildComboLigand(tLinkIns)
                print("##################################################################")
            if not self.errLevel:
                print("##################################################################")
                print("#                                                                #")
                print("#     Get information for all modifications and the link         #")    
                print("#                                                                #")
                print("##################################################################")
                self.extractOneLinkInfo(tLinkIns)

        # Finally print out 
        if not self.errLevel:
            print("##################################################################")
            print("#                                                                #")
            print("#          Output all information to the file                    #")    
            print("#                                                                #")
            print("##################################################################")
            self.outOneLinkInfo(tLinkIns)
            print("#                          Job done                              #")
            print("##################################################################")
        
    def setOneCompFromCif(self, tFileName, tMonomer):
     
        # Using the input file or the file in ccp4 monomer lib as it is
        aMmcifObj = Ccp4MmCifObj(tFileName)
        aMmcifObj.checkBlockCompsExist()
        if not aMmcifObj["errLevel"]:
            print(list(aMmcifObj["ccp4CifObj"].keys()))
            print(list(aMmcifObj["ccp4CifObj"]["comps"].keys()))
            if tMonomer["name"] in aMmcifObj["ccp4CifObj"]["comps"]:
                tMonomer["outComp"] = True
                tMonomer["comp"] = aMmcifObj["ccp4CifObj"]["comps"][tMonomer["name"]]
                if (not tMonomer["userIn"]) and (tMonomer["name"].upper() in self.chemCheck.aminoAcids):
                    self.chemCheck.tmpModiN_in_AA(tMonomer["name"].upper(), tMonomer["comp"])
                self.selectHAtoms(tMonomer["comp"])   
                print(aMmcifObj["ccp4CifObj"]["lists"]["comp"].keys())
                tMonomer["list"] = aMmcifObj["ccp4CifObj"]["lists"]["comp"][tMonomer["name"]]
                #print tMonomer["list"]
                dataHead = "data_comp_%s"%tMonomer["name"]
                #print "dataHead ", dataHead
                if dataHead in aMmcifObj["ccp4CifBlocks"]:
                    tMonomer["dataBlock"] = aMmcifObj["ccp4CifBlocks"][dataHead]
                    #for i in range(10): 
                    #    print tMonomer["dataBlock"][i].strip() 
            else:
                print("No %s in %s "%(tMonomer["name"], tFileName))
                print("Try CCP4 monomer lib ")
                #The input file does not contain the comp, try $CLIBD_MON
                if tMonomer["name"] !="":
                    aSub = tMonomer["name"][0].lower()
                    aNewCif = os.path.join(self.allChemCombDir, aSub, tMonomer["name"].upper() + ".cif")
                    print("aNewCif ", aNewCif) 
                    if os.path.isfile(aNewCif):
                        aMmcifObj = Ccp4MmCifObj(aNewCif)
                        if tMonomer["name"] in aMmcifObj["ccp4CifObj"]["comps"]:
                            tMonomer["comp"]  = aMmcifObj["ccp4CifObj"]["comps"][tMonomer["name"]]
                            self.selectHAtoms(tMonomer["comp"])   
                         
                            tMonomer["list"] = aMmcifObj["ccp4CifObj"]["lists"]["comp"][tMonomer["name"]]
                            tMonomer["inCif"] = aNewCif 
                            dataHead = "data_comp_%s"%tMonomer["name"]
                            if dataHead in aMmcifObj["ccp4CifBlocks"]:
                                tMonomer["dataBlock"] = aMmcifObj["ccp4CifBlocks"][dataHead]
                        else:
                            self.errLevel = 21
                            if self.errLevel not in self.errMessage:
                                   self.errMessage[self.errLevel] = []
                            self.errMessage[self.errLevel].append("Comp %s can not be found in both %s and %s\n"%(tMonomer["name"],\
                                                                   tMonomer["inCif"], aNewCif))
                    else:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        self.errMessage[self.errLevel].append("No %s in CCP4 monomer lib\n"%(aNewCif))
                else:        
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("One of residues is  without name\n")
              
            if not self.errLevel and "bonds" in tMonomer["comp"]:
                self.checkDelocAndAromaBonds(tMonomer["comp"], tMonomer["name"])

        else:
            self.errLevel = aMmcifObj["errLevel"]
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append(aMmcifObj["errMessage"])
          

    def checkDelocAndAromaBonds(self, tCompMonomer, tName):

        # check if a bond has the type of "deloc"

        allDelocBs = []
        if "bonds" in tCompMonomer and "atoms" in tCompMonomer:
            for iB in range(len(tCompMonomer["bonds"])):
                if "type" in tCompMonomer["bonds"][iB]: 
                    if tCompMonomer["bonds"][iB]["type"].upper().find("DELOC") != -1\
                       or tCompMonomer["bonds"][iB]["type"].upper().find("AROMATIC") != -1:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        aLine  = "The component %s.cif contains a deloc or aromatic bond. It needs to be kekulized.\n"%tName
                        aLine += "You could use acedrg to do kekulization of the bonds.\n"
                        aLine += "Try getting %s.cif from PDB/CCD and then use it as an input file\n"%tName
                        aLine += "to run acedrg to generate a cif file for the ligand description.\n"
                        aLine += "e.g. You get the file %s.cif from PDB/CCD and then,\n"%tName
                        aLine += "acedrg -c %s.cif -o  %s_fromAcedrg \n"%(tName, tName)
                        aLine += "After that, you can put the newly generated %s_fromAcedrg.cif\n"%tName
                        aLine += "into the instruction file %s, which contains(in one line), e.g. \n"%"instruc.txt"
                        aLine += "LINK: RES-NAME-1 MUR FILE-1 MUR_fromAcedrg.cif ATOM-NAME-1 C8 RES-NAME-2 HIS ATOM-NAME-2 N DELETE ATOM O9 1\n"
                        aLine += "Then run acedrg again to generate the link description, e.g. \n"
                        aLine += "acedrg -L %s   -o xxxxx(anyname) \n"%"instruct.txt" 
                        aLine += "The file generated for the link description is xxxxx_link.cif \n" 
                        self.errMessage[self.errLevel].append(aLine) 
                        break
                        """
                        if tCompMonomer["bonds"][iB]["type"].upper().find("DELOC") != -1:
                            self.modOneDelocBond(iB, tCompMonomer)
                        print "a deloc bond: "
                        print "atom 1 ",  tCompMonomer["bonds"][iB]["atom_id_1"], " element ",tCompMonomer["atoms"][idxAtm1]["type_symbol"]
                        print "atom 2 ",  tCompMonomer["bonds"][iB]["atom_id_2"], " element ",tCompMonomer["atoms"][idxAtm2]["type_symbol"]
                        print "bond type ", tCompMonomer["bonds"][iB]["type"]
                        """
    def modOneDelocBond(self, tIdxBond, tCompMonomer):

        lDone = False

        idxAtm1 = self.getAtomById(tCompMonomer["bonds"][iB]["atom_id_1"], tCompMonomer["atoms"])
        idxAtm2 = self.getAtomById(tCompMonomer["bonds"][iB]["atom_id_2"], tCompMonomer["atoms"])
        if idxAtm1 !=-1 and idxAtm2 !=-1:
            if tCompMonomer["atoms"][idxAtm1]["type_symbol"] == "O":
                 allDelocBs.append([iB,idxAtm2, idxAtm1])
            elif tCompMonomer["atoms"][idxAtm2]["type_symbol"] == "O": 
                 allDelocBs.append([iB,idxAtm1, idxAtm2])
            else:
                self.errLevel = 12
                self.errMessage[self.errLevel].append("%s or %s does not exist \n")\
                      %(tCompMonomer["bonds"][iB]["atom_id_1"], tCompMonomer["bonds"][iB]["atom_id_2"])

        

    def getAtomById(self, tId, atoms):

        aRet = -1
       
        for iA in range(len(atoms)):
            if tId == atoms[iA]["atom_id"] :
                aRet = iA
                break

        return aRet
        
           
    
    def selectHAtoms(self, tMonomer):

        tMonomer["hAtoms"] = []
        
        for aAtom in tMonomer["atoms"]:
            if aAtom["type_symbol"] == "H":
                tMonomer["hAtoms"].append(aAtom["atom_id"])
        
    def setOneMonomer(self, tMonomer):

        if os.path.isfile(tMonomer["inCif"]):
            #print(tMonomer["inCif"])
            #aNL = tMonomer["name"].upper()
            #if len(aNL) > 3:
            aNL = "UNL"
            self._log_name  = os.path.join(self.scrDir, aNL + "_for_link.log")
            self.subRoot    = os.path.join(self.scrDir, aNL + "_for_link")
            self._cmdline   = "acedrg -c %s  -r %s -o %s "%(tMonomer["inCif"], aNL, self.subRoot)   
            #print self._cmdline
            self.runExitCode = self.subExecute()
            if not self.runExitCode :
                aOutLigCif = self.subRoot + ".cif"
                if os.path.isfile(aOutLigCif): 
                    tMonomer["outCif"] = aOutLigCif
                else: 
                    self.errLevel = 31
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Run time error : no result dictionary for %s \n"%aNL) 
                
            else: 
                self.errLevel = 31
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Run time error : generating a dictionary for %s failed\n"%aNL) 

    #def adjustAtomsAndOthersForComboLigand(self, tLinkedObj):

        #self.setDeletedAtomsForModification(tLinkedObj) 

        #if len(tLinkedObj.modLigand1["deleted"]["atoms"]) !=0:
        #self.setDeletedInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
        #if len(tLinkedObj.modLigand2["deleted"]["atoms"]) !=0:
        #self.setDeletedInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)
        #self.setAddedInOneResForModification(tLinkedObj)   

    def adjustAtomsAndOthersForComboLigand(self, tLinkedObj):
     
        if len(tLinkedObj.modLigand1["changed"]["formal_charges"]) > 0:
            self.addjustFormalChargeInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1)
        elif len(tLinkedObj.modLigand2["changed"]["formal_charges"]) > 0:
            self.addjustFormalChargeInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2)
        
        self.setAddedInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
        if not self.errLevel:
            self.setDeletedInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
            if len(tLinkedObj.stdLigand1["linkChir"]) ==2:
                self.setLinkChir(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.stdLigand2)
            if not self.errLevel:
                self.setDeletedInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)
                if len(tLinkedObj.stdLigand2["linkChir"]) ==2:
                    self.setLinkChir(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.stdLigand1)
                #if not self.errLevel:
                #    self.setChargeInLinkAtom(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
                    #if not self.errLevel:
                    #    self.setChargeInLinkAtom(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)

        if self.errLevel:
            print("Error level is ", self.errLevel)
            print("Error message is : ")
            for aL in self.errMessage[self.errLevel]:
                print(aL)

    def  setLinkChir(self, tRes, tMod, tOthRes):
        
         #print("Here 1") 
         idKey = "atom_id_" + tRes["linkChir"][1]
         print(idKey)
         print("center ", tRes["linkChir"][0]["atom_id_centre"])
         print("before ", tRes["linkChir"][0][idKey])
         tRes["linkChir"][0][idKey]=tOthRes["atomName"]
         print("after ", tRes["linkChir"][0][idKey])
         print("deleted chiral center")
         print ("center atom ", tMod["deleted"]["chirs"][-1]["atom_id_centre"], " in chir ", tMod["deleted"]["chirs"][-1]['id'])
         print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_1"], " in chir ", tMod["deleted"]["chirs"][-1]['id'])
         print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_2"], " in chir ", tMod["deleted"]["chirs"][-1]['id'])
         print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_3"], " in chir ", tMod["deleted"]["chirs"][-1]['id'])
         
    def addjustFormalChargeInOneResForModification(self, tRes, tMod):

        changeAtms = []
        if len(tMod["changed"]["atoms"]) > 0:
            for aA in tMod["changed"]["atoms"]:
                changeAtms.append(aA["atom_id"].strip().upper())
        for aFC in tMod["changed"]["formal_charges"]:
            aId = aFC["atom_id"].strip().upper()
            for aAt in tRes["comp"]["atoms"]:
                if aAt["atom_id"].strip().upper() ==aId:
                    aAt["charge"] = aFC["formal_charge"]
                    if not aId in changeAtms:
                        tMod["changed"]["atoms"].append(aAt)
                    break
        
        if len(tMod["changed"]["atoms"]) > 0:
            print("The following atoms are changed: ")
            for aAt in tMod["changed"]["atoms"]:
                print("Atom ",aAt["atom_id"]) 
                print("Formal charge ", aAt["charge"])
        
    def setChargeInLinkAtom(self, tRes, tMod, tLinkBonds):

        for idxA in range(len(tRes["remainAtoms"])):
            if tRes["remainAtoms"][idxA]["atom_id"] == tRes["atomName"]:
               if tRes["remainAtoms"][idxA]["type_symbol"] != "N":
                   if "charge" in tRes["remainAtoms"][idxA]:
                       nC = int(tRes["remainAtoms"][idxA]["charge"])
                       if nC !=0:
                           if len(tLinkBonds) > 0:
                               aOr = BondOrderS2N(tLinkBonds[0]["type"])                   
                               tRes["remainAtoms"][idxA]["charge"] = str(nC+aOr)
  
    def setDeletedAtomsForModification(self, tLinkedObj):
        # Not used anymore 
        # Delete atoms as instructed 
        lDelRes1 = False
        for aDS in tLinkedObj.delSections:
            print("A delete section has the following keys : ", list(aDS.keys()))
            if "inRes" in aDS:
                print(aDS["inRes"])
                if aDS["inRes"] == 1:
                    lDelRes1  = True
                    self.deleteOneAtomAndConnectedHAtoms(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, aDS)
                    print("Number of atoms remained in residue %s is %d "%(tLinkedObj.stdLigand1["name"], len(tLinkedObj.stdLigand1["remainAtoms"])))  
                    print("Number of atoms deleted in this residue is", len(tLinkedObj.modLigand1["deleted"]["atoms"]))
                    print("They are : ")
                    for aAtom in tLinkedObj.modLigand1["deleted"]["atoms"]:
                        print(aAtom["atom_id"])  
                elif aDS["inRes"] == 2:
                    lDelRes2  = True
                    self.deleteOneAtomAndConnectedHAtoms(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, aDS)
                    print("Number of atoms remained in residue %s is %d "%(tLinkedObj.stdLigand2["name"], len(tLinkedObj.stdLigand2["remainAtoms"])))  
                    print("Number of atoms deleted in this residue is", len(tLinkedObj.modLigand2["deleted"]["atoms"]))
                    print("They are : ")
                    for aAtom in tLinkedObj.modLigand2["deleted"]["atoms"]:
                        print(aAtom["atom_id"])  
                else:
                    self.errLevel    = 21 
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Error in residue number for for DELETED atoms in the instruction file\n")


        if not lDelRes1:
            for aAtom in tLinkedObj.stdLigand1["comp"]["atoms"]:
                tLinkedObj.stdLigand1["remainAtoms"].append(aAtom)
        if not lDelRes2:
            for aAtom in tLinkedObj.stdLigand2["comp"]["atoms"]:
                tLinkedObj.stdLigand2["remainAtoms"].append(aAtom)
       
    def deleteOneAtomAndConnectedHAtoms(self, tStdMonomer, tModMonomer, tId):

        aName = tId.upper()
        connectedH = []
        self.getHAtomConnected(aName, tStdMonomer, connectedH)
        for aAtom in tStdMonomer["comp"]["atoms"]:
            if aAtom["atom_id"].upper() == aName or aAtom["atom_id"] in connectedH:
                tModMonomer["deleted"]["atoms"].append(aAtom)
 
    def getHAtomConnected(self, tNonHAtomName, tMonomer, tHNames):
        
        for aBond in tMonomer["comp"]["bonds"]:
            if aBond["atom_id_1"].upper()==tNonHAtomName and aBond["atom_id_2"].upper() in tMonomer["comp"]["hAtoms"]:
                tHNames.append(aBond["atom_id_2"])
            if aBond["atom_id_2"].upper()==tNonHAtomName and aBond["atom_id_1"].upper() in tMonomer["comp"]["hAtoms"]:
                tHNames.append(aBond["atom_id_1"])
                
    def setAddedInOneResForModification(self, tRes, tMod, tLinkBonds):

        print("For residue ", tRes["name"])
        
        if len(tMod["added"]["atoms"]) > 0:
            for aAtom in tMod["added"]["atoms"]:
                print("Add atom ", aAtom["atom_id"])
                aBool = self.checkDubAtomNameInOneRes(tRes, aAtom)
                if not aBool : 
                    self.addOneAtomAndConnectedBonds(tRes, tMod, tLinkBonds, aAtom)
                else:
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    aLine = "Atom Name Dublication : added atom %s already exists in residue %s\n"\
                            %(aAtom["atom_id"], tRes["name"])
                    print(aLine) 
                    self.errMessage[self.errLevel].append(aLine)
                    break

    def checkDubAtomNameInOneRes(self, tRes, tAtom):

        aRet = False

        existAtmIds = []
        for aAtom in tRes["comp"]["atoms"]:
            existAtmIds.append(aAtom["atom_id"])
       
        if tAtom["atom_id"] in existAtmIds:
            aRet = True

        return aRet          

    def addOneAtomAndConnectedBonds(self, tRes, tMod, tLinkBonds, tAtom):

        tRes["remainAtoms"].append(tAtom)
        tmpBandA = []
        tmpBonds = []

        nAtoms = len(tRes["comp"]["atoms"])
        for aBond in tMod["added"]["bonds"]:
            tRes["remainBonds"].append(aBond)
            if aBond["atom_id_1"]==tAtom["atom_id"]:
                atmIdx = self.getOneAtomSerialById(aBond["atom_id_2"], tRes["comp"]["atoms"])
                if atmIdx > 0 and atmIdx < nAtoms:
                    tmpBandA.append([aBond, tRes["comp"]["atoms"][atmIdx]])
                    tmpBonds.append(aBond)
                else:
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    aLine = "Can not find idx for atom ", aBond["atom_id_2"] 
                    self.errMessage[self.errLevel].append(aLine)
                    break
            elif aBond["atom_id_2"]==tAtom["atom_id"]:
                atmIdx = self.getOneAtomSerialById(aBond["atom_id_1"], tRes["comp"]["atoms"])
                if atmIdx > 0 and atmIdx < nAtoms:
                    tmpBandA.append([aBond, tRes["comp"]["atoms"][atmIdx]])
                    tmpBonds.append(aBond)
                else:
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    aLine = "Can not find idx for atom ", aBond["atom_id_1"] 
                    self.errMessage[self.errLevel].append(aLine)
                    break


        for aBond in tLinkBonds :
            if aBond["atom_id_1"]==tAtom["atom_id"]:
                tmpBonds.append(aBond)
            elif aBond["atom_id_2"]==tAtom["atom_id"]:
                tmpBonds.append(aBond)

        if len(tmpBonds) > 0:
            # Check bonds around the added atoms

            print("Check added atom %s and around bonds "%tAtom["atom_id"])
            aPair = self.chemCheck.valideBondOrderForOneOrgAtom(tAtom, tmpBonds)
            if not aPair[0]:
                self.errLevel = 12
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append(aPair[1])
                print(aPair[1])
            else:
                print("Bonds connected to atom %s are OK "%tAtom["atom_id"])
        else :
            aLine = "No bonds are defined to connected to atom %s, check your input file! \n"%tAtom["atom_id"]
            self.errLevel = 12
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append(aLine)

        if not self.errLevel:
            if len(tmpBandA) >0:
                for aPair in tmpBandA:
                    print("Added bond between atom %s and %s "%(aPair[0]["atom_id_1"], aPair[0]["atom_id_2"]))
                    print("Need to adjust bonds connected to atom %s "%aPair[1]["atom_id"])
            
            tmpDelIds =[]
            for [aBond, aAtom] in tmpBandA:
                aBandASet = self.getBondSetForOneLinkedAtom(aAtom["atom_id"], tRes["comp"]["atoms"],tRes["comp"]["bonds"], tmpDelIds)
                aBandASet[0].append(tAtom)
                aBandASet[1].append(aBond)
                if len(aBandASet[0]) > 0:
                    print("Check the bonds connected atom %s "%aAtom["atom_id"])
                    print("These bonds are : ")
                    for aB in aBandASet[1]:
                        print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                    print("Atoms involved are: ")
                    for aA in aBandASet[0]:
                        print("Atom %s "%aA["atom_id"])
                    [noErr, errLine] =self.chemCheck.adjustNBForOneAddedAtom(aAtom, aBandASet[0], aBandASet[1], tRes, tMod, tmpDelIds)
                    if not noErr:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        self.errMessage[self.errLevel].append(errLine)  
                        break

    def setDeletedInOneResForModification(self, tRes, tMod, tLinkBonds):

        print("For residue ", tRes["name"])
        delAtomIdSet =[]
        for aAtom in tMod["deleted"]["atoms"]:
            delAtomIdSet.append(aAtom["atom_id"])

        if len(delAtomIdSet) > 0:
            tMod["deleted"]["atoms"] = []
            for aId in delAtomIdSet:
                self.deleteOneAtomAndConnectedHAtoms(tRes, tMod, aId)            
        print("Number of deleted atoms according to the instruction file is ", len(delAtomIdSet))
        self.resetChargeForLinkedNAtom(tRes, tMod, tLinkBonds, delAtomIdSet) 
        self.adjustAtomsAroundOneAtom(tRes["atomName"], tRes, tMod, tLinkBonds, delAtomIdSet, 1)
        extraAtomSet = []
        if len(tMod["changed"]["bonds"]) > 0:
            for aB in tMod["changed"]["bonds"]:
                if aB["atom_id_1"] == tRes["atomName"]:
                    extraAtomSet.append(aB["atom_id_2"])
                elif aB["atom_id_2"] == tRes["atomName"]:
                    extraAtomSet.append(aB["atom_id_1"])
        for atmId in extraAtomSet:
            self.adjustAtomsAroundOneAtom(atmId, tRes, tMod, tLinkBonds, delAtomIdSet, 2)
  
        delAtomIdSet =[]
        for aAtom in tMod["deleted"]["atoms"]:
            print(list(aAtom.keys()))
            delAtomIdSet.append(aAtom["atom_id"])
 
        if len(delAtomIdSet) > 0:
           print("After atom deleting procedure, the following atoms are deleted :")
           for aId in delAtomIdSet:
               print("Atom ", aId)
        aTmpRemain = []
        for aA in tRes["remainAtoms"]:
            aTmpRemain.append(aA)
        tRes["remainAtoms"] = []
        for aAtom in tRes["comp"]["atoms"]:
            if not aAtom["atom_id"] in delAtomIdSet:
                tRes["remainAtoms"].append(aAtom)
        for aA in aTmpRemain:
            tRes["remainAtoms"].append(aA)
        print("Number of total atoms is ", len(tRes["comp"]["atoms"]))
        print("Number of remained atoms is ", len(tRes["remainAtoms"])) 
        print("Those atoms are :")
        for aAtom in  tRes["remainAtoms"]:
            print("Atom ", aAtom["atom_id"])
   
        # Delete all bonds that contains the deleted atom, change the bond-order for those changed bonds
        tmpRemBs = []
        for aB in tRes["remainBonds"]:
            tmpRemBs.append(aB)
        tRes["remainBonds"] = []
        i = 0
        chBondIdMap = {}
        print(" changed bonds  ", len(tMod["changed"]["bonds"]))
        for chBond in tMod["changed"]["bonds"]:
            aList = [chBond["atom_id_1"], chBond["atom_id_2"]]
            aList.sort()
            aStr = aList[0] + "_" + aList[1]
            chBondIdMap[aStr] = i  
            i = i+1
        for aBond in tRes["comp"]["bonds"]:
            if (aBond["atom_id_1"].upper() in delAtomIdSet) or (aBond["atom_id_2"].upper() in delAtomIdSet):
                tMod["deleted"]["bonds"].append(aBond)
            else:
                bList = [aBond["atom_id_1"], aBond["atom_id_2"]]
                bList.sort()
                bStr = bList[0] + "_" + bList[1]
                if bStr in list(chBondIdMap.keys()):
                    tRes["remainBonds"].append(tMod["changed"]["bonds"][chBondIdMap[bStr]])
                else:
                    tRes["remainBonds"].append(aBond)

        for aB in tmpRemBs:
            tRes["remainBonds"].append(aB)

        print("Number of remained bonds : ", len(tRes["remainBonds"]))
        print("They are : ")
        for aBond in tRes["remainBonds"]:
            print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))
            print("Bond-order is %s "%aBond["type"])     
        print("Number of deleted bonds : ", len(tMod["deleted"]["bonds"]))
        if len(tMod["deleted"]["bonds"]):
            print("They are : ")
            for aBond in tMod["deleted"]["bonds"]: 
                print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))      
                print("Bond-order is %s "%aBond["type"])  
 
        # Delete all angles that contain the deleted atom
        for aAng in tRes["comp"]["angles"]:
            if (aAng["atom_id_1"].upper() in delAtomIdSet)\
               or (aAng["atom_id_2"].upper() in delAtomIdSet)\
               or (aAng["atom_id_3"].upper() in delAtomIdSet):
                tMod["deleted"]["angles"].append(aAng)
            else:
                tRes["remainAngs"].append(aAng)
        print("Number of remained Angles : ", len(tRes["remainAngs"]))
        print("They are : ")
        for aAng in tRes["remainAngs"]:
            print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
        print("Number of deleted angles ", len(tMod["deleted"]["angles"]))
        print("They are : ")
        for aAng in tMod["deleted"]["angles"]:
            print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
                  
        # Delete all tors that contain the deleted atom
        if "tors" in tRes["comp"]:
            for aTor in tRes["comp"]["tors"]:
                if (aTor["atom_id_1"].upper() in delAtomIdSet)\
                   or (aTor["atom_id_2"].upper() in delAtomIdSet)\
                   or (aTor["atom_id_3"].upper() in delAtomIdSet)\
                   or (aTor["atom_id_4"].upper() in delAtomIdSet):
                    tMod["deleted"]["tors"].append(aTor)
                else:
                    tRes["remainTors"].append(aTor)
                       
        # Delete all chiral centers that contain the deleted atom
        # as the center atom.
        # Change signs of chiral volumes of chiral centers
        # where one of NB atoms replaced.
        if "chirs" in tRes["comp"]:
            for aChi in tRes["comp"]["chirs"]:
                if aChi["atom_id_centre"] != tRes["atomName"]:
                    if aChi["atom_id_centre"].upper() in delAtomIdSet\
                       or aChi["atom_id_1"].upper() in delAtomIdSet\
                       or aChi["atom_id_2"].upper() in delAtomIdSet\
                       or aChi["atom_id_3"].upper() in delAtomIdSet:
                        tMod["deleted"]["chirs"].append(aChi)
                        if aChi["atom_id_centre"].upper() in delAtomIdSet:
                            print("atom ", aChi["atom_id_centre"], " in chir ", aChi['id'], " is deleted ")
                        elif aChi["atom_id_1"].upper() in delAtomIdSet:
                            print("atom ", aChi["atom_id_1"], " in chir ", aChi['id'], " is deleted ")
                        elif aChi["atom_id_2"].upper() in delAtomIdSet:
                            print("atom ", aChi["atom_id_2"], " in chir ", aChi['id'], " is deleted ")
                        elif aChi["atom_id_3"].upper() in delAtomIdSet:
                            print("atom ", aChi["atom_id_3"], " in chir ", aChi['id'], " is deleted ")
                        print("Chiral center ", aChi['id'], " is deleted ")
                    else:
                        tRes["remainChirs"].append(aChi)
                        #print("remain chiral center")
                        #print ("center atom ", aChi["atom_id_centre"], " in chir ", aChi['id'])
                        #print ("atom ", aChi["atom_id_1"], " in chir ", aChi['id'])
                        #print ("atom ", aChi["atom_id_2"], " in chir ", aChi['id'])
                        #print ("atom ", aChi["atom_id_3"], " in chir ", aChi['id'])
                else:
                    delIdSet = []
                    if aChi["atom_id_1"].upper() in delAtomIdSet:
                         delIdSet.append("1")
                    if aChi["atom_id_2"].upper() in delAtomIdSet:
                         delIdSet.append("2")
                    if aChi["atom_id_3"].upper() in delAtomIdSet:
                         delIdSet.append("3")
                    if len(delIdSet) > 1 :
                        # Original chiral center destroyed. Let combo-conformer to determine if a chir exists.
                        tMod["deleted"]["chirs"].append(aChi)
                    elif len(delIdSet)==1 :
                        tMod["deleted"]["chirs"].append(aChi)
                        tRes["linkChir"].append(self.copyChi(aChi))
                        tRes["linkChir"].append(delIdSet[0])
                        print("deleted chiral center")
                        print ("center atom ", tMod["deleted"]["chirs"][-1]["atom_id_centre"], " in chir ", aChi['id'])
                        print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_1"], " in chir ", aChi['id'])
                        print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_2"], " in chir ", aChi['id'])
                        print ("atom ", tMod["deleted"]["chirs"][-1]["atom_id_3"], " in chir ", aChi['id'])
                    else:
                        tRes["remainChirs"].append(aChi)
                        tRes["linkChir"].append(self.copyChi(aChi))

        #print("Here 0")
        # Delete the deleted atom from a plane and even delete a plane 
        # if number of atoms in a plane smaller then 3 (after deleting the assigned atom)
        if "planes" in tRes["comp"]:
            for aPl in list(tRes["comp"]["planes"].keys()):
                lDel = False
                aPlGrp = []
                for aPAtm in tRes["comp"]["planes"][aPl]:
                    if aPAtm["atom_id"] in delAtomIdSet:
                        lDel = True
                    else:
                        aPlGrp.append(aPAtm)
                if lDel:
                    tMod["deleted"]["planes"].append(tRes["comp"]["planes"][aPl])
                else:
                    tRes["remainPls"].append(aPlGrp)

    def copyChi(self, tChi):

        aChi = {}
        for aKey in list(tChi.keys()):
            aChi[aKey] = tChi[aKey]
        return aChi

    def resetChargeForLinkedNAtom(self,  tRes, tMod, tLinkBonds, tDelAtomIds):
        # For atoms of type_symbol "N"  only at the moment.
        aLAtmId = tRes["atomName"]
        aLAtmSerial = self.getOneAtomSerialById(aLAtmId, tRes["comp"]["atoms"])
        if aLAtmSerial != -1:
            aLAtmElem = tRes["comp"]["atoms"][aLAtmSerial]["type_symbol"]
            if aLAtmElem == "N":
                aLAAtoms    = []
                aTmpLABonds = []
                [aLAAtoms,aTmpLABonds] =  self.getBondSetForOneLinkedAtom(aLAtmId, tRes["comp"]["atoms"], tRes["comp"]["bonds"], tDelAtomIds)
                print("Number of connected atoms for N atom", len(aTmpLABonds))
                aName = aLAtmId.upper()
                connectedH = []
                self.getHAtomConnected(aName, tRes, connectedH)
                if len(aTmpLABonds)==4 and len(connectedH)==3:
                    if "charge" in tRes["comp"]["atoms"][aLAtmSerial]:
                        tRes["comp"]["atoms"][aLAtmSerial]["charge"] = "0"
                    elif "form_charge" in tRes["comp"]["atoms"][aLAtmSerial]:
                        tRes["comp"]["atoms"][aLAtmSerial]["form_charge"] = "0"
                    print("Charges in atom %s is set to 0 "%aLAtmId)  
           
    def adjustAtomsAroundOneAtom(self, tCenAtomId, tRes, tMod, tLinkBonds, tDelAtomIds, tMode):
        #aLAtmId = tRes["atomName"]
        aLAtmId  = tCenAtomId
        aLAtmSerial = self.getOneAtomSerialById(aLAtmId, tRes["comp"]["atoms"])
        tExtraChAtms = []
        #print "center atom: ", aLAtmSerial
        print("center atom: ", aLAtmId) 
        addAtmIds = []
        if len(tMod["added"]["atoms"]) !=0:
            for aA in tMod["added"]["atoms"]:
                addAtmIds.append(aA["atom_id"])
        if aLAtmSerial != -1 and not aLAtmId in addAtmIds:
            aLAtmElem = tRes["comp"]["atoms"][aLAtmSerial]["type_symbol"]
            # Atoms and bonds around the linked atom
            [aLAAtoms,aTmpLABonds] =  self.getBondSetForOneLinkedAtom(aLAtmId, tRes["comp"]["atoms"], tRes["comp"]["bonds"], tDelAtomIds)
            print("connected atoms", len(aTmpLABonds))
            aLABonds = []    
            if len(tMod["changed"]["bonds"]) > 0:
                # Consider the effect of bond-order changes for some bonds
                aChBonds = {}
                for aB in tMod["changed"]["bonds"]:
                    tBIdList1 = [aB["atom_id_1"], aB["atom_id_2"]]
                    tBIdList1.sort()
                    aStr = tBIdList1[0] + "_" + tBIdList1[1]
                    aChBonds[aStr] = []
                    aChBonds[aStr].append(aB)
                for aB in aTmpLABonds:
                    tBIdList2 = [aB["atom_id_1"], aB["atom_id_2"]]
                    tBIdList2.sort()
                    bStr = tBIdList2[0] + "_" + tBIdList2[1]
                    print(bStr)
                    if bStr in list(aChBonds.keys()):
                        aLABonds.append(aChBonds[bStr][0])
                    else:
                        aLABonds.append(aB) 
            else:
                for aB in aTmpLABonds:
                    aLABonds.append(aB)
        
            if tMode == 1:
                # Add the linked bond
                aLABonds.append(tLinkBonds[0])
            print("Number of bonds ", len(aLABonds))
            print("tMode ", tMode) 
            print("The linked atom %s in residue %s "%(aLAtmId, tRes["name"]))
            if len(aLABonds):
                print("It appears in the following bonds now ")
                for aB in aLABonds:
                    print("Bond between %s and %s of order %s "%(aB["atom_id_1"], aB["atom_id_2"], aB["type"]))

                nTotalVa=self.getTotalBondOrderInOneMmcifAtom(aLAtmId, aLABonds)

                print("total Valence is ", nTotalVa)
                print("atom ", aLAtmElem.upper())
                print("Default Valence is ", self.chemCheck.defaultBo[aLAtmElem.upper()])
                if aLAtmElem in self.chemCheck.orgVal:
                    if "charge" in tRes["comp"]["atoms"][ aLAtmSerial]:
                        allowedBO = self.chemCheck.orgVal[aLAtmElem][0] + int(tRes["comp"]["atoms"][ aLAtmSerial]["charge"])
                    else:
                        allowedBO = self.chemCheck.orgVal[aLAtmElem][0]
                    print("Allowed order is ", allowedBO)
                    if nTotalVa != allowedBO:
                        lUseOther = False
                        if len(self.chemCheck.orgVal[aLAtmElem]) > 1:
                            for i in range(1, len(self.chemCheck.orgVal[aLAtmElem])):
                                if nTotalVa == self.chemCheck.orgVal[aLAtmElem][i]:
                                    lUseOther = True
                                    break

                        if not lUseOther:
                            allHIds = []
                            aLAHIds = []
                            for aA in tRes["comp"]["atoms"]:
                                if aA["type_symbol"] == "H" :
                                    allHIds.append(aA["atom_id"])
                            for aA in aLAAtoms:
                                if aA["type_symbol"] == "H" and not aA["atom_id"] in tDelAtomIds:
                                    aLAHIds.append(aA["atom_id"])
                                    print("H atom ", aA["atom_id"])
                            if nTotalVa > allowedBO:
                                aN = nTotalVa-allowedBO
                                if aN <=len(aLAHIds) :
                                    print("%d H atom will be deleted "%aN)
                                    print("%s connects %d H atom "%(tRes["comp"]["atoms"][aLAtmSerial]["atom_id"], len(aLAHIds)))
                                    tmpDelAtomSet = []
                                    aLAHIds.sort()
                                    idxD = -1
                                    for i in range(aN):
                                        aHName = aLAHIds[idxD]
                                        tmpDelAtomSet.append(aHName)
                                        idxD = idxD -1
                                    
                                    print("The following H atom is deleted ")
                                    for aId in tmpDelAtomSet:
                                        print("Atom ", aId) 
 
                                    for aAtom in tRes["comp"]["atoms"]:
                                        if aAtom["atom_id"] in tmpDelAtomSet:
                                            tMod["deleted"]["atoms"].append(aAtom)

                                elif len(self.chemCheck.orgVal[aLAtmElem]) > 1:
                                    nV = len(self.chemCheck.orgVal[aLAtmElem])
                                    lAdded = False
                                    for i in range(1,nV):
                                        if "charge" in tRes["comp"]["atoms"][ aLAtmSerial]:
                                            allowedBO1 = self.chemCheck.orgVal[aLAtmElem][i] + int(tRes["comp"]["atoms"][ aLAtmSerial]["charge"])
                                        else:
                                            allowedBO1 = self.chemCheck.orgVal[aLAtmElem][i]
                                        aN = allowedBO1 - nTotalVa
                                        if aN==1:
                                            print("H connected atom id ", tRes["comp"]["atoms"][aLAtmSerial]["atom_id"])
                                            self.addOneHInRes(tRes["comp"]["atoms"][aLAtmSerial], aLAAtoms, allHIds, tRes, tMod)
                                            lAdded = True
                                            break
                                    if not lAdded:
                                        self.errLevel    = 22 
                                        if self.errLevel not in self.errMessage:
                                            self.errMessage[self.errLevel] = []
                                            aL = "Currently it is not allowed to do such a change for atom %s\n"\
                                                 %(tRes["comp"]["atoms"][aLAtmSerial]["atom_id"])
                                            self.errMessage[self.errLevel].append(aL)

                            elif nTotalVa < self.chemCheck.orgVal[aLAtmElem][0]:
                                aN = allowedBO - nTotalVa
                                if aN == 1:
                                    print("H connected atom id ", tRes["comp"]["atoms"][aLAtmSerial]["atom_id"])
                                    self.addOneHInRes(tRes["comp"]["atoms"][aLAtmSerial], aLAAtoms, allHIds, tRes, tMod)
        else :
            if not aLAtmId in addAtmIds:
                self.errLevel    = 22 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Can not find the linked atom in residue %s "%tRes["name"])
                self.errMessage[self.errLevel].append("Check your input instruction file ")
        
    def getBondSetForOneLinkedAtom(self, tAtomId, tAtoms, tBonds, tDelAtomIds):

        aAtomSet = []
        aBondSet = []
        #print "atom ", tAtomId
        #print "Number of bonds ", len(tBonds)
        for aBond in tBonds:
            aId1 = aBond["atom_id_1"].strip()
            aId2 = aBond["atom_id_2"].strip()
            #print "bond between atom ", aId1, " and atom  ", aId2
            if aId1 == tAtomId and not aId2 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print "a bond found "
            if aId2 == tAtomId and not aId1 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print "a bond found "
        #print "Number of bonds found ", len(aBondSet)
        return [aAtomSet, aBondSet]      

    def getBondSetForOneAtomByAlias(self, tAtomId, tAtoms, tBonds, tDelAtomIds):

        aAtomSet = []
        aBondSet = []
        print("atom ", tAtomId)
        print("Number of bonds ", len(tBonds))
        for aBond in tBonds:
            aId1 = aBond["atom_id_1_alias"].strip()
            aId2 = aBond["atom_id_2_alias"].strip()
            if aId1 == tAtomId and not aId2 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                print("bond between atom ", aId1, " and atom  ", aId2)
                print("a bond found ")
            if aId2 == tAtomId and not aId1 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                print("bond between atom ", aId1, " and atom  ", aId2)
                print("a bond found ")
        print("Number of bonds found ", len(aBondSet))
        return [aAtomSet, aBondSet]      
        
    def getTotalBondOrderInOneMmcifAtom(self, tAtomId, tBonds):
        # tBonds is prefilted in getBondSetForOneMmcifAtom() and other
        totalOr = 0
        for aBond in tBonds:
            aOr = BondOrderS2N(aBond["type"])
            if aOr != -1:
                totalOr += aOr
            else:
                totalOr = -1
                break

        return totalOr
   

    def addOneHInRes(self, tHConnAtom, tOtheAtmSet, tAllAtmIds, tRes, tMod):

        aAtom = {}
        aAtom["comp_id"]     = tRes["name"]
        aAtom["atom_id"]     = self.chemCheck.setHName(tHConnAtom, tOtheAtmSet, tAllAtmIds)         
        aAtom["type_symbol"] = "H"
        aAtom["type_energy"] = "H"
        aAtom["charge"]      = "0"
        tRes["remainAtoms"].append(aAtom)
        tMod["added"]["atoms"].append(aAtom)
        print("an H atom %s is added into %s"%(aAtom["atom_id"], aAtom["comp_id"]))                 
        aBond = {}
        aBond["comp_id"]         = tRes["name"]
        aBond["atom_id_1"]       = tHConnAtom["atom_id"]
        aBond["atom_id_2"]       = aAtom["atom_id"]
        aBond["type"]            = "single"
        aBond["value_dist"]      = "0.860"
        aBond["value_dist_esd"]  = "0.02"
        tRes["remainBonds"].append(aBond)
        tMod["added"]["bonds"].append(aBond)
        print("a bond between  %s and %s is added into bonds in %s"%(aBond["atom_id_1"], aBond["atom_id_2"], aBond["comp_id"]))                 
    
    def modOneChir(self, tLinkBonds, tCenAtmID, tKw, tChir, tResChirs):

        for aBond in tLinkBonds:
            if aBond["atom_id_1"] ==tCenAtmID:
                tChir[tKw]=aBond["atom_id_2"]
                aChir["volume_sign"] = "both"
                tResChirs.append(aChir)
                break
            elif aBond["atom_id_2"] ==tCenAtmID:
                tChir[tKw]=aBond["atom_id_1"]
                tChir["volume_sign"] = "both"
                tResChirs.append(tChir)
                break

    def setInitComboLigand(self, tLinkedObj):

        # All atoms and bonds in the combo-ligand in the linkedObj
 
        tLinkedObj.combLigand["atoms"] = []
        idxAtmComb = 0

        idxA1      = 0
        for aAtom in tLinkedObj.stdLigand1["remainAtoms"]:
            aAtom["idx"]             = idxA1
            aAtom["idxCombo"]        = idxAtmComb
            if aAtom["type_symbol"]=="H":
                tLinkedObj.stdLigand1["remainIdxHs"].append(idxA1)
                tLinkedObj.combLigand["hAtoms"].append(idxAtmComb)
            idxA1 +=1
            idxAtmComb +=1
            aAtmName = aAtom["type_symbol"] + str(idxAtmComb)   
            aAtom["atom_id_alias"] = aAtmName  
            tLinkedObj.atomMap[aAtmName] = [1, aAtom["atom_id"]] 
            tLinkedObj.combLigand["atoms"].append(aAtom)
            if aAtom["atom_id"] == tLinkedObj.stdLigand1["atomName"]:
                tLinkedObj.stdLigand1["atomName_alias"] =aAtmName
 
        idxA2      = 0
        for aAtom in tLinkedObj.stdLigand2["remainAtoms"]:
            aAtom["idx"]             = idxA2
            aAtom["idxCombo"]        = idxAtmComb
            if aAtom["type_symbol"]=="H":
                tLinkedObj.stdLigand2["remainIdxHs"].append(idxA2)
                tLinkedObj.combLigand["hAtoms"].append(idxAtmComb)
            idxA2 +=1
            idxAtmComb +=1   
            aAtmName = aAtom["type_symbol"] + str(idxAtmComb)   
            aAtom["atom_id_alias"] = aAtmName 
            tLinkedObj.atomMap[aAtmName] = [2, aAtom["atom_id"]] 
            tLinkedObj.combLigand["atoms"].append(aAtom)
            if aAtom["atom_id"] == tLinkedObj.stdLigand2["atomName"]:
                tLinkedObj.stdLigand2["atomName_alias"] =aAtmName 
       
         
        tLinkedObj.combLigand["bonds"] = []
        for aBond in tLinkedObj.stdLigand1["remainBonds"]:
            aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
            aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand1["remainAtoms"])
            if aBond["atom_id_1_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine = "An error happens in residue 1:\n"
                aLine += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"]
                self.errMessage[self.errLevel].append(aLine) 
            if aBond["atom_id_2_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine = "An error happens in residue 1:\n"
                aLine+= "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
                self.errMessage[self.errLevel].append(aLine)
            if not self.errLevel:
                tLinkedObj.combLigand["bonds"].append(aBond)
        for aBond in tLinkedObj.stdLigand2["remainBonds"]:
            aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand2["remainAtoms"])
            aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
            if aBond["atom_id_1_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine = "An error happens in residue 2:\n"
                aLine += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"] 
                self.errMessage[self.errLevel].append(aLine)
            if aBond["atom_id_2_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine  = "An error happens in residue 2:\n"
                aLine += "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
                self.errMessage[self.errLevel].append(aLine)
            if not self.errLevel:
                tLinkedObj.combLigand["bonds"].append(aBond)

        # The added bond as input instruction file
        for aBond in tLinkedObj.suggestBonds:
            aBond["comp_id"]         = tLinkedObj.combLigand["name"]
            aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
            print("bond ", aBond["comp_id"])
            print("atom 1 alias ", aBond["atom_id_1_alias"])
            if aBond["atom_id_1_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine = "An error happens in residue 1:\n"
                aLine += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"] 
                self.errMessage[self.errLevel].append(aLine)
            aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
            print("atom 2 alias ", aBond["atom_id_2_alias"])
            if aBond["atom_id_2_alias"] == "":
                self.errLevel    = 32 
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                aLine = "An error happens in residue 2:\n"
                aLine += "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
                self.errMessage[self.errLevel].append(aLine)
            aBond["value_dist"]      =  0.0             
            aBond["value_dist_esd"]  =  0.20 
            if not self.errLevel:            
                tLinkedObj.combLigand["bonds"].append(aBond)

        # Chiral centers 
        #tLinkedObj.combLigand["chirs"]
        for aChir in  tLinkedObj.stdLigand1["remainChirs"]:
            aChir["atom_id_centre_alias"] = self.getAtomAlias(aChir["atom_id_centre"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_1_alias"]      = self.getAtomAlias(aChir["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_2_alias"]      = self.getAtomAlias(aChir["atom_id_2"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_3_alias"]      = self.getAtomAlias(aChir["atom_id_3"], tLinkedObj.stdLigand1["remainAtoms"])
            if aChir["atom_id_centre_alias"] and aChir["atom_id_1_alias"] and aChir["atom_id_2_alias"] and aChir["atom_id_3_alias"]:
                tLinkedObj.combLigand["chirs"].append(aChir)
        if len(tLinkedObj.stdLigand1["linkChir"]) ==2:
            tLinkedObj.stdLigand1["linkChir"][0]["atom_id_centre_alias"] =\
                             self.getAtomAlias(tLinkedObj.stdLigand1["atomName"], tLinkedObj.stdLigand1["remainAtoms"])
            print ("link atom, centre alias : ", tLinkedObj.stdLigand1["linkChir"][0]["atom_id_centre_alias"])
            otherId = "atom_id_" + tLinkedObj.stdLigand1["linkChir"][1]
            idKey = "atom_id_" + tLinkedObj.stdLigand1["linkChir"][1] + "_alias"
            nN =0
            for aC in ["1","2", "3"]:
                curKey1 = "atom_id_" + aC
                curKey2 = curKey1 + "_alias"
                if curKey2 != idKey:
                    tLinkedObj.stdLigand1["linkChir"][0][curKey2]  = \
                                  self.getAtomAlias(tLinkedObj.stdLigand1["linkChir"][0][curKey1], tLinkedObj.stdLigand1["remainAtoms"])   
                else:
                    tLinkedObj.stdLigand1["linkChir"][0][curKey2]  = \
                                  self.getAtomAlias(tLinkedObj.stdLigand1["linkChir"][0][curKey1], tLinkedObj.stdLigand2["remainAtoms"])   
                print(curKey2,  tLinkedObj.stdLigand1["linkChir"][0][curKey2])
            tLinkedObj.combLigand["chirs"].append(tLinkedObj.stdLigand1["linkChir"][0])
        for aChir in  tLinkedObj.stdLigand2["remainChirs"]:
            aChir["atom_id_centre_alias"] = self.getAtomAlias(aChir["atom_id_centre"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_1_alias"] = self.getAtomAlias(aChir["atom_id_1"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_2_alias"] = self.getAtomAlias(aChir["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_3_alias"] = self.getAtomAlias(aChir["atom_id_3"], tLinkedObj.stdLigand2["remainAtoms"])
            if aChir["atom_id_centre_alias"] and aChir["atom_id_1_alias"] and aChir["atom_id_2_alias"] and aChir["atom_id_3_alias"]:
                tLinkedObj.combLigand["chirs"].append(aChir)

        if len(tLinkedObj.stdLigand2["linkChir"]) ==2:
            tLinkedObj.stdLigand2["linkChir"][0]["atom_id_centre_alias"] =\
                             self.getAtomAlias(tLinkedObj.stdLigand2["atomName"], tLinkedObj.stdLigand2["remainAtoms"])
            print ("link atom, centre alias : ", tLinkedObj.stdLigand2["linkChir"][0]["atom_id_centre_alias"])
            otherId = "atom_id_" + tLinkedObj.stdLigand2["linkChir"][1]
            idKey = "atom_id_" + tLinkedObj.stdLigand2["linkChir"][1] + "_alias"
            for aC in ["1","2", "3"]:
                curKey1 = "atom_id_" + aC
                curKey2 = curKey1 + "_alias"
                
                if curKey2 != idKey:
                    tLinkedObj.stdLigand2["linkChir"][0][curKey2]  =  \
                                        self.getAtomAlias(tLinkedObj.stdLigand2["linkChir"][0][curKey1], tLinkedObj.stdLigand2["remainAtoms"])  
                else:
                    tLinkedObj.stdLigand2["linkChir"][0][curKey2]  = \
                                        self.getAtomAlias(tLinkedObj.stdLigand2["linkChir"][0][curKey1], tLinkedObj.stdLigand1["remainAtoms"])
                print(curKey2,  tLinkedObj.stdLigand2["linkChir"][0][curKey2])
            tLinkedObj.combLigand["chirs"].append(tLinkedObj.stdLigand2["linkChir"][0])
  
        self.outAtomNameMapToJSon(tLinkedObj)

    def outAtomNameMapToJSon(self, tLinkedObj):

        #aOutFN = os.path.join(self.scrDir, tLinkedObj.combLigand["name"] + "_AtomNameMapping.json")
        aOutFN = os.path.join(self.scrDir, "UNL_AtomNameMapping.json")
        try: 
            aOutF = open(aOutFN, "w")
        except IOError:
            print("%s can not be open for writing "%aOutFN)
            self.errLevel    = 32
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for write "%tOutFName)
        else:       
            resNames =["UNL", tLinkedObj.stdLigand1["name"], tLinkedObj.stdLigand2["name"]]
            aOutF.write("{\n\"ResidueNameMapping\":\n")
            aOutF.write("    [{ \"RESIDUE_1\":\"%s\"},\n"%tLinkedObj.stdLigand1["name"])
            aOutF.write("     { \"RESIDUE_2\":\"%s\"}],\n"%tLinkedObj.stdLigand2["name"])
            aOutF.write("\"AtomNameMapping\":\n")
            nLs = len(list(tLinkedObj.atomMap.keys()))
            idxL = 0
            aOutF.write("    [\n")
            for aCombAtmName in list(tLinkedObj.atomMap.keys()):
                aCombAtmNameT  = "\"" + aCombAtmName + "\","
                aAtmName  = "\"" + tLinkedObj.atomMap[aCombAtmName][1] + "\","
                aResIndex = tLinkedObj.atomMap[aCombAtmName][0]
                aLine = ""
                if idxL != nLs -1:
                    aLine = "    { \"ALIAS\":%s \"ATOM_NAME\":%s \"RESIDUE_INDEX\":%d },\n"\
                            %((aCombAtmNameT).ljust(8), (aAtmName).ljust(8), aResIndex)
                else:
                    aLine = "    { \"ALIAS\":%s \"ATOM_NAME\":%s \"RESIDUE_INDEX\":%d }\n"\
                            %((aCombAtmNameT).ljust(8), (aAtmName).ljust(8), aResIndex)
                    
                aOutF.write(aLine)
                idxL = idxL +1
            aOutF.write("    ]\n")
            aOutF.write("}\n")

            aOutF.close()
             
    def getOneAtomSerialById(self, tId, tAtoms):

        nReturn = -1
        #print tId
        for i in range(len(tAtoms)):
            #print tAtoms[i]["atom_id"]
            if tAtoms[i]["atom_id"]== tId:
                nReturn = i
                break

        return nReturn            
            
    def getAtomAlias(self, tAtomName, tAtoms):

        aReturnAlias = ""

        for aAtom in tAtoms:
            if aAtom["atom_id"]==tAtomName:
                if "atom_id_alias" in aAtom:
                    aReturnAlias = aAtom["atom_id_alias"]
                else:
                    self.errLevel    = 33 
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Atom %s does not have alias name\n"%aAtom["atom_id"]) 
                break
        return aReturnAlias

    def setTorIdsInOneLink(self, tTors):

        if len(tTors) !=0:

           aTorIdMap = {}
           newStrs = []
           aSeriNum = 0 
           for aTor in tTors:
               idStrs = aTor["id"].strip().split("_")
               if len(idStrs)==3:
                   newStrs = idStrs[:2]
               elif len(idStrs)==4:
                   newStrs = idStrs[1:3]
    
               if len(newStrs)==2:
                   newId = newStrs[0] + "_" + newStrs[1]
                   if not newId in list(aTorIdMap.keys()):
                       aTorIdMap[newId] = []
                   aTorIdMap[newId].append(aSeriNum)
               aSeriNum = aSeriNum + 1       

           if len(aTorIdMap) >0:
               for aKey in list(aTorIdMap.keys()):
                   aNum = 1
                   for aIdx in aTorIdMap[aKey]:
                       tTors[aIdx]["id"] = aKey + "_" + str(aNum)
                       aNum = aNum + 1

    def outTmpComboLigandMap(self, tLinkObj):

        # Output the mapping between 2 monomers and the combo-ligand for checking
 
        tmpFName = os.path.join(self.scrDir, "tmpComboLigandMap.txt")

        try:
            tmpF = open(tmpFName, "w")
        except IOError:
            self.errLevel    = 33 
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for reading "%tmpFName)
        else:
            tmpF.write("The mapping between original atom names and their alias\n")
            tmpF.write("Monomer 1 \n")
            for aAtom in tLinkObj.stdLigand1["remainAtoms"]:
                tmpF.write("%10d1%10d%s%s\n"%(aAtom["idxCombo"], aAtom["idx"], aAtom["atom_id"].ljust(10), aAtom["atom_id_alias"].ljust(10)))
            tmpF.write("\nMonomer 2 \n")
            for aAtom in tLinkObj.stdLigand2["remainAtoms"]:
                tmpF.write("%10d1%10d%s%s\n"%(aAtom["idxCombo"], aAtom["idx"], aAtom["atom_id"].ljust(10), aAtom["atom_id_alias"].ljust(10)))

            tmpF.write("The mapping of bonds: original atom names vs their alias\n")

            for aBond in tLinkedObj.combLigand["bonds"]:
                tmpF.write("-------------------------\n")
                tmpF.write("%s%s\n"%(aBond["atom_id_1"].ljust(), aBond["atom_id_2"].ljust()))
                tmpF.write("%s%s\n"%(aBond["atom_id_1_alias"].ljust(), aBond["atom_id_2_alias"].ljust()))
            tmpF.write("-------------------------\n")
            tmpF.close()

    def checkChemInMonomer(self, tMonomer, tMode):

        tDelAtomIds = []
        for aAtom in tMonomer["atoms"]:
            atmId = aAtom["atom_id"]
            atmElm = aAtom["type_symbol"]
            if tMode == 2:
                atmId = aAtom["atom_id_alias"]
            [bondAtomSet, aLABonds] =  self.getBondSetForOneAtomByAlias(atmId, tMonomer["atoms"], tMonomer["bonds"], tDelAtomIds)
            nTotalVa=self.getTotalBondOrderInOneMmcifAtom(atmId, aLABonds) 
            aCharge = 0
            if "charge" in aAtom:
                aCharge = int(aAtom["charge"])
            if aCharge != 0: 
                nTotalVa = nTotalVa - aCharge
            print("atom ", atmId, "charge ", aCharge, " equiv bond-order ", nTotalVa)
            if atmElm in self.chemCheck.orgVal:
                if not nTotalVa in self.chemCheck.orgVal[atmElm]:
                    self.errLevel    = 45
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("atom %s in monomer %s has a total valence of %d, which is not allowed!\n"\
                                                          %(aAtom["atom_id"], aAtom["comp_id"], nTotalVa))
                    break 

    def comboLigToSimplifiedMmcif(self, tMonomer, tOutFName):
      
        try: 
            aOutF = open(tOutFName, "w")
        except IOError:
            print("%s can not be open for reading "%tOutFName)
            self.errLevel    = 34 
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for reading "%tOutFName)
        else:

            # A list block
            aOutF.write("data_comp_list\n")
            aOutF.write("loop_\n")
            aOutF.write("_chem_comp.id\n")
            aOutF.write("_chem_comp.three_letter_code\n")
            aOutF.write("_chem_comp.name\n")
            aOutF.write("_chem_comp.group\n")
            aOutF.write("_chem_comp.number_atoms_all\n")
            aOutF.write("_chem_comp.number_atoms_nh\n")
         
            tCombId = "UNL"   
            numAts = len(tMonomer["atoms"]) 
            numH   = len(tMonomer["hAtoms"]) 
            numNonH = numAts - numH
            aOutF.write("%s%s%s%s%6d%6d\n"%(tCombId.ljust(8), tCombId.ljust(8), "\'.           \'".ljust(20),\
                                           "non-polymer".ljust(20), numAts, numNonH))
            aOutF.write("data_comp_UNL\n")
            aOutF.write("#\n")
            if numAts != 0: 
                aOutF.write("loop_\n")
                aOutF.write("_chem_comp_atom.comp_id\n")
                aOutF.write("_chem_comp_atom.atom_id\n")
                aOutF.write("_chem_comp_atom.type_symbol\n")
                aOutF.write("_chem_comp_atom.type_energy\n")
                aOutF.write("_chem_comp_atom.charge\n")
                aOutF.write("_chem_comp_atom.x\n")
                aOutF.write("_chem_comp_atom.y\n")
                aOutF.write("_chem_comp_atom.z\n")

                    
                for aAtom in tMonomer["atoms"]:
                    aCharge = "0"
                    if "charge" in aAtom:
                        aAtom["charge"] = str(aAtom["charge"])
                        if aAtom["charge"].find(".")==-1:
                            aCharge = aAtom["charge"]
                        else:
                            aCharge = aAtom["charge"].strip().split(".")[0]
                    if "x" in aAtom and "y" in aAtom and "z" in aAtom: 
                        x = aAtom["x"]
                        y = aAtom["y"]
                        z = aAtom["z"]
                    else:
                        x = "0.000"
                        y = "0.000"
                        z = "0.000"

                    aOutF.write("%s%s%s%s%s%s%s%s\n"%(tCombId.ljust(10), aAtom["atom_id_alias"].ljust(8),\
                                aAtom["type_symbol"].ljust(6), aAtom["type_energy"].ljust(8),\
                                aCharge.ljust(10), x.ljust(10), y.ljust(10), z.ljust(10)))
                
            else: 
                self.errLevel    = 35
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("The combo-ligand has no atoms\n")
            if not self.errLevel:
                if len(tMonomer["bonds"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_bond.comp_id\n")
                    aOutF.write("_chem_comp_bond.atom_id_1\n")
                    aOutF.write("_chem_comp_bond.atom_id_2\n")
                    aOutF.write("_chem_comp_bond.type\n")
                    aOutF.write("_chem_comp_bond.value_dist\n")
                    aOutF.write("_chem_comp_bond.value_dist_esd\n")
                    for aBond in tMonomer["bonds"]:
                        aOutF.write("%s%s%s%s%s%s\n"%(tCombId.ljust(10), aBond["atom_id_1_alias"].ljust(8),\
                                     aBond["atom_id_2_alias"].ljust(8), aBond["type"].ljust(20), \
                                     str(aBond["value_dist"]).ljust(10), str(aBond["value_dist_esd"]).ljust(10)))
                 
                else: 
                    self.errLevel    = 35
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("The combo-ligand has no bonds\n")

                if not self.errLevel and len(tMonomer["chirs"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_chir.comp_id\n")
                    aOutF.write("_chem_comp_chir.id\n")
                    aOutF.write("_chem_comp_chir.atom_id_centre\n")
                    aOutF.write("_chem_comp_chir.atom_id_1\n")
                    aOutF.write("_chem_comp_chir.atom_id_2\n")
                    aOutF.write("_chem_comp_chir.atom_id_3\n")
                    aOutF.write("_chem_comp_chir.volume_sign\n")
                    i = 1
                    for aChi in tMonomer["chirs"]:
                        aChiId = "chir_" + str(i)
                        aOutF.write("%s%s%s%s%s%s%s\n"%(tCombId.ljust(10), aChiId.ljust(10), aChi["atom_id_centre_alias"].ljust(6),\
                                                        aChi["atom_id_1_alias"].ljust(6), aChi["atom_id_2_alias"].ljust(6),\
                                                        aChi["atom_id_3_alias"].ljust(6), aChi["volume_sign"]))
                        i +=1
                aOutF.close()   
                        
    def buildComboLigand(self, tLinkedObj):

        if not self.errLevel:

            self.adjustAtomsAndOthersForComboLigand(tLinkedObj)
            tLinkedObj.combLigand["name"] = tLinkedObj.stdLigand1["name"].strip() + "-" + tLinkedObj.stdLigand2["name"].strip()
            print("The name of combo-ligand : %s "%tLinkedObj.combLigand["name"])
            self.setInitComboLigand(tLinkedObj)
           
            if not self.errLevel: 
                print("Number of atoms in the combo-ligand is ", len(tLinkedObj.combLigand["atoms"]))
                print("They are : ")
                for aAtom in tLinkedObj.combLigand["atoms"]:
                    print("%s%s%s"%(aAtom["atom_id"].ljust(10), aAtom["atom_id_alias"].ljust(10), aAtom["type_symbol"]))
                for aBond in tLinkedObj.combLigand["bonds"]:
                    print("%s%s%s%s%s\n"%(aBond["atom_id_1_alias"].ljust(10), aBond["atom_id_2_alias"].ljust(10),
                                   ("("+aBond["atom_id_1"]).ljust(10), (aBond["atom_id_2"]+ ")").ljust(10),
                                    aBond["type"].ljust(1) )) 
                #self.outTmpComboLigandMap(tLinkedObj)    # Check 
                self.checkChemInMonomer(tLinkedObj.combLigand, 2)
                if not self.errLevel:
                    tLinkedObj.combLigand["inCif"] = os.path.join(self.scrDir, tLinkedObj.combLigand["name"] + "_comboIn.cif")
                    print("The cif file of the combo-ligand for input ", tLinkedObj.combLigand["inCif"])
                    self.comboLigToSimplifiedMmcif(tLinkedObj.combLigand, tLinkedObj.combLigand["inCif"])
                    self.setOneMonomer(tLinkedObj.combLigand)
                    if not self.errLevel:
                        tLinkedObj.outCombLigand["name"] = tLinkedObj.combLigand["name"]
                        tLinkedObj.outCombLigand["cifObj"] = Ccp4MmCifObj(tLinkedObj.combLigand["outCif"])["ccp4CifObj"]            
                        print("output comboLigand name ", tLinkedObj.outCombLigand["name"])
                        #print("Number of atoms in the comboLigand : ", len(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]))
                        comboScrDN = self.subRoot + "_TMP"
                        if os.path.isdir(comboScrDN) and not self.testMode:
                            #print "Delete the tempo dir ", comboScrDN
                            shutil.rmtree(comboScrDN)
                        comboLigMolFileName = self.subRoot + ".mol"
                        self.fileTool.MmCifToMolFile(tLinkedObj.combLigand["outCif"], comboLigMolFileName)
            else:
                print(self.errLevel)            
       
    def getChangesInModificationFromCombLigand(self, tLinkedObj):

  
        #for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]:
        #    print(aAtom.keys())
        #    print aAtom["atom_id"]
        #    print "Atom %s is in residue %s "%(aAtom["atom_id"], aAtom["res_idx"])
        #    print "Is it added ? ", aAtom["is_added"]

        addedSet1 = []
        addedSet2 = []
        # Atoms
 
        existChangeAtmIdsRes1 = []
        if len(tLinkedObj.modLigand1["changed"]["atoms"]) > 0:
            for aA in tLinkedObj.modLigand1["changed"]["atoms"]:
                existChangeAtmIdsRes1.append(aA["atom_id"])
        existChangeAtmIdsRes2 = []
        if len(tLinkedObj.modLigand2["changed"]["atoms"]) > 0:
            for aA in tLinkedObj.modLigand2["changed"]["atoms"]:
                existChangeAtmIdsRes2.append(aA["atom_id"])
         
        for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]:
            #print("Atom id ",aAtom["atom_id"])
            #print("Atom id alias ", aAtom["atom_id_alias"]) 
            if aAtom["res_idx"] == 1:
                if not aAtom["is_added"]:
                    #print("It is in changed section ")
                    self.checkAtomMod(tLinkedObj.stdLigand1["remainAtoms"], aAtom, tLinkedObj.modLigand1["changed"]["atoms"], existChangeAtmIdsRes1)
                else:
                    #print("It is in added section ")
                    addedSet1.append(aAtom["atom_id"])
                    tLinkedObj.modLigand1["added"]["atoms"].append(aAtom)
            elif aAtom["res_idx"] == 2:
                if not aAtom["is_added"]:
                    self.checkAtomMod(tLinkedObj.stdLigand2["remainAtoms"], aAtom, tLinkedObj.modLigand2["changed"]["atoms"], existChangeAtmIdsRes2)
                else:
                    addedSet2.append(aAtom["atom_id"])
                    tLinkedObj.modLigand2["added"]["atoms"].append(aAtom)
            else:
                self.errLevel    = 36
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Bug: Atom %s does not attach to any residue\n"%aAtom["atom_id"])

        print("Num of mod atoms in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["atoms"]))
        print("Num of add atoms in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["atoms"]))
        print("Num of deleted atoms in residue 1 is %d "%(len(tLinkedObj.modLigand1["deleted"]["atoms"])))
        print("Num of mod atoms in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["atoms"]))
        print("Num of add atoms in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["atoms"]))
        print("Num of deleted atoms in residue 2 is %d "%(len(tLinkedObj.modLigand2["deleted"]["atoms"])))
        # Bonds
        aTmpChBonds_1 =  {}
        i1 = 0
        for aB in tLinkedObj.modLigand1["changed"]["bonds"]:
            aList1 = [aB["atom_id_1"], aB["atom_id_2"]]
            aList1.sort()
            aStr1 = aList1[0] + "_" + aList1[1]
            aTmpChBonds_1[aStr1] = [aB]
            i1= i1 + 1
        tLinkedObj.modLigand1["changed"]["bonds"] = []
        aTmpChBonds_2 =  {}
        i2 = 0
        for aB in tLinkedObj.modLigand2["changed"]["bonds"]:
            aList2 = [aB["atom_id_1"], aB["atom_id_2"]]
            aList2.sort()
            aStr2 = aList2[0] + "_" + aList2[1]
            aTmpChBonds_2[aStr2] = [aB]
            i2= i2 + 1
        tLinkedObj.modLigand2["changed"]["bonds"] = []
        
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            #print("A bond ")
            #print("atom %s in residue %d "%(aBond["atom_id_1"], aBond["atom_id_1_resNum"]))
            #print("atom %s in residue %d "%(aBond["atom_id_2"], aBond["atom_id_2_resNum"]))
            if aBond["atom_id_1_resNum"]==1 and aBond["atom_id_2_resNum"]==1:
                # In residue 1
                if not aBond["atom_id_1"] in addedSet1 and not aBond["atom_id_2"] in addedSet1:
                    aList = [aBond["atom_id_1"], aBond["atom_id_2"]]
                    aList.sort()
                    aStr = aList[0] + "_" + aList[1]
                    if not aStr in list(aTmpChBonds_1.keys()):
                        self.checkBondMod(tLinkedObj.stdLigand1["remainBonds"], aBond, tLinkedObj.modLigand1["changed"]["bonds"])
                    else:
                        tLinkedObj.modLigand1["changed"]["bonds"].append(aBond) 
                else:
                    tLinkedObj.modLigand1["added"]["bonds"].append(aBond)

            elif aBond["atom_id_1_resNum"]==2 and aBond["atom_id_2_resNum"]==2:
                # In residue 2
                if not aBond["atom_id_1"] in addedSet2 and not aBond["atom_id_2"] in addedSet2:
                    aList = [aBond["atom_id_1"], aBond["atom_id_2"]]
                    aList.sort()
                    aStr = aList[0] + "_" + aList[1]
                    if not aStr in list(aTmpChBonds_2.keys()):
                        self.checkBondMod(tLinkedObj.stdLigand2["remainBonds"], aBond, tLinkedObj.modLigand2["changed"]["bonds"])
                    else:
                        tLinkedObj.modLigand2["changed"]["bonds"].append(aBond) 
                else:
                    tLinkedObj.modLigand2["added"]["bonds"].append(aBond)

        print("Number of changed bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["bonds"]))  
        print("Number of added bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["bonds"]))  
        print("Number of deleted bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["bonds"]))  
        print("Number of changed bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["bonds"]))  
        print("Number of added bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["bonds"]))  
        print("Number of deleted bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["bonds"]))  

        if len(tLinkedObj.modLigand1["added"]["atoms"]) > 0:
            for aAtom in tLinkedObj.modLigand1["added"]["atoms"]:
                addedSet1.append(aAtom["atom_id"])
        if len(tLinkedObj.modLigand2["added"]["atoms"]) > 0:
            for aAtom in tLinkedObj.modLigand2["added"]["atoms"]:
                addedSet2.append(aAtom["atom_id"])
        # Angles
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            #print("A angle ")
            #print("atom %s in residue %d "%(aAng["atom_id_1"], aAng["atom_id_1_resNum"]))
            #print("atom %s in residue %d "%(aAng["atom_id_2"], aAng["atom_id_2_resNum"]))
            #print("atom %s in residue %d "%(aAng["atom_id_3"], aAng["atom_id_3_resNum"]))
            if aAng["atom_id_1_resNum"]==1 and aAng["atom_id_2_resNum"]==1 and aAng["atom_id_3_resNum"]==1:
                if not aAng["atom_id_1"] in addedSet1 and not aAng["atom_id_2"] in addedSet1 and not aAng["atom_id_3"] in addedSet1:
                    self.checkAngMod(tLinkedObj.stdLigand1["remainAngs"], aAng, tLinkedObj.modLigand1["changed"]["angles"]) 
                else:
                    tLinkedObj.modLigand1["added"]["angles"].append(aAng) 
            if aAng["atom_id_1_resNum"]==2 and aAng["atom_id_2_resNum"]==2 and aAng["atom_id_3_resNum"]==2:
                if not aAng["atom_id_1"] in addedSet2 and not aAng["atom_id_2"] in addedSet2 and not aAng["atom_id_3"] in addedSet2:
                    self.checkAngMod(tLinkedObj.stdLigand2["remainAngs"], aAng, tLinkedObj.modLigand2["changed"]["angles"]) 
                else:
                    tLinkedObj.modLigand2["added"]["angles"].append(aAng) 
         
        print("Number of changed angles in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["angles"]))  
        print("Number of added angles in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["angles"]))  
        print("Number of deleted angles in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["angles"]))  
        print("Number of changed angles in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["angles"]))  
        print("Number of added angles in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["angles"]))  
        print("Number of deleted angles in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["angles"])) 
        # Torsions
        if "tors" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
                #print("A torsion ")
                #print("atom %s in residue %d "%(aTor["atom_id_1"], aTor["atom_id_1_resNum"]))   
                #print("atom %s in residue %d "%(aTor["atom_id_2"], aTor["atom_id_2_resNum"]))   
                #print("atom %s in residue %d "%(aTor["atom_id_3"], aTor["atom_id_3_resNum"]))  
                #print("atom %s in residue %d "%(aTor["atom_id_4"], aTor["atom_id_4_resNum"]))  
                if aTor["atom_id_1_resNum"]==1 and aTor["atom_id_2_resNum"]==1\
                   and aTor["atom_id_3_resNum"]==1 and aTor["atom_id_4_resNum"]==1:
                    if not aTor["atom_id_1"] in addedSet1 and not aTor["atom_id_2"] in addedSet1\
                       and not aTor["atom_id_3"] in addedSet1 and not aTor["atom_id_4"] in addedSet1: 
                        self.checkTorMod(tLinkedObj.stdLigand1["remainTors"], aTor, tLinkedObj.modLigand1["changed"]["tors"])
                    else: 
                        tLinkedObj.modLigand1["added"]["tors"].append(aTor)
                elif aTor["atom_id_1_resNum"]==2 and aTor["atom_id_2_resNum"]==2\
                   and aTor["atom_id_3_resNum"]==2 and aTor["atom_id_4_resNum"]==2:
                    if not aTor["atom_id_1"] in addedSet2 and not aTor["atom_id_2"] in addedSet2\
                       and not aTor["atom_id_3"] in addedSet2 and not aTor["atom_id_4"] in addedSet2: 
                        self.checkTorMod(tLinkedObj.stdLigand2["remainTors"], aTor, tLinkedObj.modLigand2["changed"]["tors"])
                    else: 
                        tLinkedObj.modLigand2["added"]["tors"].append(aTor)
    
            print("Number of changed tors in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["tors"]))  
            print("Number of added tors in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["tors"]))  
            print("Number of deleted tors in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["tors"]))  
            print("Number of changed tors in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["tors"]))  
            print("Number of added tors in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["tors"]))  
            print("Number of deleted tors in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["tors"])) 
         
        # Chirs
        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
                print("A chiral center ")
                print("atom %s in residue %d "%(aChi["atom_id_centre"], aChi["atom_id_centre_resNum"]))  
                print("atom %s in residue %d "%(aChi["atom_id_1"], aChi["atom_id_1_resNum"]))   
                print("atom %s in residue %d "%(aChi["atom_id_2"], aChi["atom_id_2_resNum"]))   
                print("atom %s in residue %d "%(aChi["atom_id_3"], aChi["atom_id_3_resNum"]))  
                if aChi["atom_id_centre_resNum"] == 1 and aChi["atom_id_1_resNum"] == 1 and\
                   aChi["atom_id_2_resNum"] == 1 and aChi["atom_id_3_resNum"] == 1: 
                    if not aChi["atom_id_centre"] in addedSet1 and not aChi["atom_id_1"] in addedSet1\
                       and not aChi["atom_id_2"] in addedSet1 and not aChi["atom_id_3"] in addedSet1:
                        self.checkChiMod(tLinkedObj.stdLigand1["remainChirs"], aChi, tLinkedObj.modLigand1["changed"]["chirs"])
                    else:
                        tLinkedObj.modLigand1["added"]["chirs"].append(aChi)
                elif aChi["atom_id_centre_resNum"] == 2 and  aChi["atom_id_1_resNum"] == 2 and\
                     aChi["atom_id_2_resNum"] == 2 and  aChi["atom_id_3_resNum"] == 2: 
                    if not aChi["atom_id_centre"] in addedSet2 and not aChi["atom_id_1"] in addedSet2\
                       and not aChi["atom_id_2"] in addedSet2 and not aChi["atom_id_3"] in addedSet2:
                        self.checkChiMod(tLinkedObj.stdLigand2["remainChirs"], aChi, tLinkedObj.modLigand2["changed"]["chirs"])
                    else:
                        tLinkedObj.modLigand2["added"]["chirs"].append(aChi)

            print("Number of changed chirs in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["chirs"]))  
            print("Number of added chirs in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["chirs"]))  
            print("Number of deleted chirs in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["chirs"]))  
            print("Number of changed chirs in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["chirs"]))  
            print("Number of added chirs in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["chirs"]))  
            print("Number of deleted chirs in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["chirs"])) 
         
        # Planes
        self.checkPlModFromCombo(tLinkedObj)
        self.checkPlModFromOrig(tLinkedObj)
     
        nDP1 = len(tLinkedObj.modLigand1["deleted"]["planes"])   
        nAP1 = len(tLinkedObj.modLigand1["added"]["planes"])
        nDP2 = len(tLinkedObj.modLigand2["deleted"]["planes"])   
        nAP2 = len(tLinkedObj.modLigand2["added"]["planes"])
        print("For residue 1 : ")
        print("Number of deleted planes ", nDP1)
        if nDP1:
            print("They are : ")
            for aPl in tLinkedObj.modLigand1["deleted"]["planes"]:
                print("----------- A plane -------------")
                for aPlAtm in aPl:
                    print("%s    %s  "%(aPlAtm["plane_id"], aPlAtm["atom_id"]))
            print("------------------ -------------")
        print("Number of added planes ", nAP1)
        if nAP1:
            print("They are : ")
            for aPl in tLinkedObj.modLigand1["added"]["planes"]:
                print("----------- A plane -------------")
                for aPlAtm in aPl:
                    print("%s    %s  "%(aPlAtm["plane_id"], aPlAtm["atom_id"]))
            print("------------------ -------------")
        print("\nFor residue 2 : ")
        print("Number of deleted planes ", nDP2)
        if nDP2:
            print("They are : ")
            for aPl in tLinkedObj.modLigand2["deleted"]["planes"]:
                print("----------- A plane -------------")
                for aPlAtm in aPl:
                    print("%s    %s  "%(aPlAtm["plane_id"], aPlAtm["atom_id"]))
            print("------------------ -------------")
        print("Number of added planes ", nAP2)
        if nAP2:
            print("They are : ")
            for aPl in tLinkedObj.modLigand2["added"]["planes"]:
                print("----------- A plane -------------")
                for aPlAtm in aPl:
                    print("%s    %s  "%(aPlAtm["plane_id"], aPlAtm["atom_id"]))
            print("------------------ -------------")

    def checkAtomMod(self, tOriAtoms, tAtom, tModAtoms, tExistChangedAtoms):

        aId = tAtom["atom_id"]
 
        for aAtm in tOriAtoms:
            if aAtm["atom_id"] == aId:
                if self.compare2Atoms(aAtm, tAtom) and not aId in tExistChangedAtoms:
                    tModAtoms.append(tAtom)
   
    def compare2Atoms(self, tOriAtom, tAtom):

        lChange = False
        if "type_energy" in tOriAtom and "type_energy" in tAtom:
            if tOriAtom["type_energy"] !=tAtom["type_energy"]: 
                tOriAtom["type_energy"] =tAtom["type_energy"]
                lChange = True
        if "charge" in tOriAtom and "charge" in tAtom:
            if tOriAtom["charge"] !=tAtom["charge"]:
                lCharge = True

        return lChange
        
    def checkBondMod(self, tOrigBonds, tBond, tModBonds):

        id1 = tBond["atom_id_1"]
        id2 = tBond["atom_id_2"]
        
        for aB in tOrigBonds:
            if (aB["atom_id_1"]==id1 and aB["atom_id_2"]==id2) or\
               (aB["atom_id_2"]==id1 and aB["atom_id_1"]==id2):
                if self.compare2Bonds(aB, tBond):
                    tModBonds.append(tBond)
                break

    def compare2Bonds(self, tOriBond, tBond):
   
        lChanged = False
        if tOriBond["type"].upper()[:3] !=tBond["type"].upper()[:3]:
            lChanged = True
        if float(tOriBond["value_dist"]) != float(tBond["value_dist"]):
            lChanged = True
        return lChanged
   
    def checkAngMod(self, tOrigAngs, tAng, tModAngs):

        id1 = tAng["atom_id_1"]
        id2 = tAng["atom_id_2"]
        id3 = tAng["atom_id_3"]
        
        for aAng in tOrigAngs:
            if aAng["atom_id_2"]==id2:
                if (aAng["atom_id_1"]==id1 and aAng["atom_id_3"]==id3) or\
                   (aAng["atom_id_1"]==id3 and aAng["atom_id_3"]==id1):
                    if self.compare2Angs(aAng, tAng):
                        tModAngs.append(tAng)
                    break

    def compare2Angs(self, tOrigAng, tAng):
   
        lChanged = False

        if float(tOrigAng["value_angle"]) != float(tAng["value_angle"]):
            lChanged = True
        return lChanged
   
    def checkTorMod(self, tOrigTors, tTor, tModTors):

        id1 = tTor["atom_id_1"]
        id2 = tTor["atom_id_2"]
        id3 = tTor["atom_id_3"]
        id4 = tTor["atom_id_4"]
        
        for aTor in tOrigTors:
            if (aTor["atom_id_2"]==id2 and aTor["atom_id_3"]==id3):
                if (aTor["atom_id_1"]==id1 and aTor["atom_id_4"]==id4):
                    if self.compare2Tors(aTor, tTor):
                        tModTors.append(tTor)
                    break
            elif (aTor["atom_id_2"]==id3 and aTor["atom_id_3"]==id2):
                if (aTor["atom_id_1"]==id4 and aTor["atom_id_4"]==id1):
                    if self.compare2Tors(aTor, tTor):
                        tModTors.append(tTor)
                    break

    def compare2Tors(self, tOrigTor, tTor):
   
        lChanged = False
        if float(tOrigTor["value_angle"]) !=float(tTor["value_angle"]):
            lChanged = True
        return lChanged
   
    def checkChiMod(self, tOrigChirs, tChi, tModChi):

        pass

    def compare2Chis(self, tOrigChi, tChi):

        pass

    def checkPlModFromCombo(self, tLinkedObj):

        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys()):
                print("For plane ", aPl)
                inRes1 = []
                inRes2 = []
                nAtmInPl = len(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl])
                for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                    print("atom %s in residue %d "%(aAtom["atom_id"], aAtom["atom_id_resNum"]))
                    if aAtom["atom_id_resNum"]== 1:
                        inRes1.append(aAtom["atom_id"])
                    elif aAtom["atom_id_resNum"]== 2:
                        inRes2.append(aAtom["atom_id"])
                if len(inRes1) == nAtmInPl:
                    self.checkPlMod(tLinkedObj.stdLigand1["remainPls"], tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl],\
                                    tLinkedObj.modLigand1["deleted"]["planes"], tLinkedObj.modLigand1["added"]["planes"]) 
                elif len(inRes2) == nAtmInPl:
                    self.checkPlMod(tLinkedObj.stdLigand2["remainPls"], tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl],\
                                tLinkedObj.modLigand2["deleted"]["planes"], tLinkedObj.modLigand2["added"]["planes"])
     
    def checkPlModFromOrig(self, tLinkedObj):

        for aPl in tLinkedObj.stdLigand1["remainPls"]:
            if not "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
                tLinkedObj.modLigand1["deleted"]["planes"].append(aPl)  
            else:
                self.checkPlMod2(aPl, tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"], tLinkedObj.modLigand1["deleted"]["planes"])

        for aPl in tLinkedObj.stdLigand2["remainPls"]:
            if not "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
                tLinkedObj.modLigand2["deleted"]["planes"].append(aPl)  
            else:
                self.checkPlMod2(aPl, tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"], tLinkedObj.modLigand2["deleted"]["planes"])

    def checkPlMod(self, tOrigPls,  tPl, tModDelPls, tModAddPls):

        tPlAtmIds = []
        print("Atoms in tPl : ")
        for aAtm in tPl:
            #if aAtm["atom_id"][0] != "H":
            tPlAtmIds.append(aAtm["atom_id"])
            print(aAtm["atom_id"])
        nAtoms = len(tPlAtmIds)
        for aPl in tOrigPls:
            print("Atoms in aPl : ")
            tOrigPlAtmIds = []
            for aPlAtm in aPl:
                #if aPlAtm["atom_id"][0] != "H":
                tOrigPlAtmIds.append(aPlAtm["atom_id"])
                print(aPlAtm["atom_id"])
            nOrigAtoms = len(tOrigPlAtmIds)
            if nAtoms > nOrigAtoms: 
                print("here1")
                overlapAtms = []
                for aOId in tOrigPlAtmIds:
                    if aOId in tPlAtmIds:
                        overlapAtms.append(aOId)
                print ("overlaped atoms ",  len(overlapAtms))
                if len(overlapAtms) == nOrigAtoms:
                    tModDelPls.append(aPl) 
                    tModAddPls.append(tPl)
                    print("orignal plane deleted")
                    break
            elif nAtoms < nOrigAtoms: 
                print("here2")
                overlapAtms = []
                for aOId in tPlAtmIds:
                    if aOId in tOrigPlAtmIds:
                        overlapAtms.append(aOId)
                if len(overlapAtms) == nAtoms:
                    print ("overlaped atoms ",  len(overlapAtms))
                    for aId in overlapAtms:
                        print ("overlaped atom : %s "%aId)
                    tModDelPls.append(aPl) 
                    tModAddPls.append(tPl)
                    print("orignal plane deleted")
                    break

        """
        for aPl in tModDelPls:
            tOrigPlAtmIds = []
            for aPlAtm in aPl:
                if aPlAtm["atom_id"][0] != "H":
                    tOrigPlAtmIds.append(aPlAtm["atom_id"])
            nOrigAtoms = len(tOrigPlAtmIds)
            if nAtoms > nOrigAtoms: 
                overlapAtms = []
                for aOId in tOrigPlAtmIds:
                    if aOId in tPlAtmIds:
                        overlapAtms.append(aOId)
                if len(overlapAtms) == nOrigAtoms:
                    tModAddPls.append(tPl)
                    break
        """
    

    def checkPlMod2(self, tPl, tComboPls, tModDelPls):

        tPlAtmIds = []
        for aAtm in tPl:
            if aAtm["atom_id"][0] != "H":
                tPlAtmIds.append(aAtm["atom_id"])

        lOverlapPls = False
        for aPl in tComboPls:
            overlapAtms = []
            for aPlAtm in tComboPls[aPl]:
                if aPlAtm["atom_id"][0] != "H":
                    if aPlAtm["atom_id"] in tPlAtmIds:
                        overlapAtms.append(aPlAtm["atom_id"])
            if len(overlapAtms) == len(tPlAtmIds):
                lOverlapPls = True
                break

        if not lOverlapPls:
            tModDelPls.append(tPl)

 
            
    def extractOneLinkInfo(self, tLinkedObj):

        self.reIndexCombLigand(tLinkedObj) 
        self.outComboPdb(tLinkedObj)
        self.getLinkInfo(tLinkedObj)
        self.getChangesInModificationFromCombLigand(tLinkedObj)
 
    def getOneOrigAtomIdFromAlias(self, tAtomMap, tAlias):

        aReturn = [-1, ""]

        if tAlias in tAtomMap:
            aReturn = [tAtomMap[tAlias][0], tAtomMap[tAlias][1]]

        return aReturn

    def getOneOrigAtomFromAlias(self, tAtomMap, tAtom):
   
        # tAtoms come from self.combLigand["atoms"]
        # outCombLigand["cifObj"] contains ids which are actually alias_id

        aReturn = False

        if tAtom["atom_id_alias"] in tAtomMap:
            tAtom["res_idx"] = tAtomMap[tAtom["atom_id_alias"]][0]
            tAtom["atom_id"] = tAtomMap[tAtom["atom_id_alias"]][1]
            aReturn = True
        return aReturn

    def getResidueIdxFromBonding(self, tBonds, tAtmMap, tAtom):

        lFind = False
        for aBond in tBonds:
            if aBond["atom_id_1"] == tAtom["atom_id_alias"]:
                if aBond["atom_id_2"] in tAtmMap:
                    tAtom["res_idx"] = tAtmMap[aBond["atom_id_2"]][0]
                    lFind = True
                    break
            elif aBond["atom_id_2"] == tAtom["atom_id_alias"]:
                if aBond["atom_id_1"] in tAtmMap:
                    tAtom["res_idx"] = tAtmMap[aBond["atom_id_1"]][0]
                    lFind = True
                    break

        if not lFind:
            self.errLevel    = 37
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("Atom %s does not attach any atoms in residue 1 or 2 ")
                    
    def reIndexCombLigand(self, tLinkedObj):
  
        addedAtoms = []

        # Atoms 
        addAtoms1 =[]
        addAtoms2 =[]
        for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]:
            aAtom["is_added"] = False 
            aAtom["atom_id_alias"] = aAtom["atom_id"]
            lSet = self.getOneOrigAtomFromAlias(tLinkedObj.atomMap, aAtom)
            if not lSet:
                aAtom["is_added"] = True 
                addedAtoms.append(aAtom)

        if len(addedAtoms) !=0:
            for aAtom in addedAtoms:
                self.getResidueIdxFromBonding(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"], \
                                              tLinkedObj.atomMap, aAtom)
            
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            aBond["atom_id_1_alias"] = aBond["atom_id_1"]
            [aBond["atom_id_1_resNum"], aBond["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aBond["atom_id_1_alias"]) 
            aBond["atom_id_2_alias"] = aBond["atom_id_2"]
            [aBond["atom_id_2_resNum"], aBond["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aBond["atom_id_2_alias"]) 

        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            aAng["atom_id_1_alias"] = aAng["atom_id_1"]
            [aAng["atom_id_1_resNum"], aAng["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_1_alias"]) 
            aAng["atom_id_2_alias"] = aAng["atom_id_2"]
            [aAng["atom_id_2_resNum"], aAng["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_2_alias"]) 
            aAng["atom_id_3_alias"] = aAng["atom_id_3"]
            [aAng["atom_id_3_resNum"], aAng["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_3_alias"]) 

        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
            aTor["atom_id_1_alias"] = aTor["atom_id_1"]
            [aTor["atom_id_1_resNum"], aTor["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_1_alias"]) 
            aTor["atom_id_2_alias"] = aTor["atom_id_2"]
            [aTor["atom_id_2_resNum"], aTor["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_2_alias"]) 
            aTor["atom_id_3_alias"] = aTor["atom_id_3"]
            [aTor["atom_id_3_resNum"], aTor["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_3_alias"]) 
            aTor["atom_id_4_alias"] = aTor["atom_id_4"]
            [aTor["atom_id_4_resNum"], aTor["atom_id_4"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_4_alias"]) 

        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
                aChi["atom_id_centre_alias"] = aChi["atom_id_centre"]
                [aChi["atom_id_centre_resNum"], aChi["atom_id_centre"]]=\
                                  self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_centre_alias"])  
                aChi["atom_id_1_alias"] = aChi["atom_id_1"]
                [aChi["atom_id_1_resNum"], aChi["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_1_alias"]) 
                aChi["atom_id_2_alias"] = aChi["atom_id_2"]
                [aChi["atom_id_2_resNum"], aChi["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_2_alias"]) 
                aChi["atom_id_3_alias"] = aChi["atom_id_3"]
                [aChi["atom_id_3_resNum"], aChi["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_3_alias"]) 
        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys()):
                for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                    aAtom["atom_id_alias"] = aAtom["atom_id"]
                    [aAtom["atom_id_resNum"], aAtom["atom_id"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAtom["atom_id_alias"]) 

    def outComboPdb(self, tLinkedObj):

        import gemmi
        aPdb = self.subRoot + ".pdb"
        if os.path.isfile(aPdb):
            aMol = gemmi.read_pdb(aPdb)
            if len(aMol) > 0:
                aMod =  aMol[0]
                if len(aMod) > 0:
                    aChain = aMod[0]
                    if len(aChain) > 0:
                        aRes = aChain[0]
                        if len(aRes) > 0:
                            allAtoms = aRes
                            for aAtom in allAtoms:
                                print("atom name before: %s"%aAtom.name)
                                if aAtom.name in tLinkedObj.atomMap:
                                    aAtom.name = tLinkedObj.atomMap[aAtom.name][1]
                                    if aAtom.name.count("\"")==2:
                                        aAtom.name=aAtom.name.strip("\"")
                                print("atom name after: %s"%aAtom.name)
                        else:
                            print("Residues in the combo-ligand have no atom ")
                            self.errLevel = 32
                            if self.errLevel not in self.errMessage:
                                self.errMessage[self.errLevel] = []
                            self.errMessage[self.errLevel].append("Residues in the combo-ligand have no atom ")
                    else:
                        print("Chains in the combo-ligand are empty ")
                        self.errLevel = 32
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        self.errMessage[self.errLevel].append("Chains in the combo-ligand are empty")

            else:
                print("Models in the combo-ligand are empty ")
                self.errLevel = 32
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Models in the combo-ligand are empty")

            aOutPdbName = self.outRoot + "_link.pdb"
            aMol.write_minimal_pdb(aOutPdbName)

        else:
            print("File %s doest not exist"%aPdb)
            self.errLevel = 32
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("Run time error :  %s does not exist\n"%aPdb)

    def getLinkInfo(self, tLinkedObj):
        # 
        atm1 = tLinkedObj.stdLigand1["atomName_alias"]
        atm2 = tLinkedObj.stdLigand2["atomName_alias"]
        tLinkedObj.cLink["bonds"] = []
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            if aBond["atom_id_1_alias"] == atm1 and aBond["atom_id_2_alias"] == atm2:
                aBond["atom_1_comp_id"] =  tLinkedObj.stdLigand1["name"]
                aBond["atom_2_comp_id"] =  tLinkedObj.stdLigand2["name"]
                tLinkedObj.cLink["bonds"].append(aBond)
            elif aBond["atom_id_1_alias"] == atm2 and aBond["atom_id_2_alias"] == atm1:
                aBond["atom_1_comp_id"] =  tLinkedObj.stdLigand2["name"]
                aBond["atom_2_comp_id"] =  tLinkedObj.stdLigand1["name"]
                tLinkedObj.cLink["bonds"].append(aBond)
        print("Number of Link bonds ", len(tLinkedObj.cLink["bonds"]))
        print("They are :")
        for aBond in tLinkedObj.cLink["bonds"]:
            print("Bond between atom %s in residue %d and atom %s in residue %d"\
                  %(aBond["atom_id_1"], aBond["atom_id_1_resNum"], aBond["atom_id_2"], aBond["atom_id_2_resNum"]))

        tLinkedObj.cLink["angles"] = []
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            if aAng["atom_id_2_alias"]==atm1:
                if aAng["atom_id_1_alias"]==atm2 or aAng["atom_id_3_alias"]==atm2:
                    tLinkedObj.cLink["angles"].append(aAng) 
            elif aAng["atom_id_2_alias"]==atm2:
                if aAng["atom_id_1_alias"]==atm1 or aAng["atom_id_3_alias"]==atm1:
                    tLinkedObj.cLink["angles"].append(aAng) 
        print("Number of Link angles ", len(tLinkedObj.cLink["angles"]))
        print("They are :")
        for aAng in tLinkedObj.cLink["angles"]:
            print("Angle between atom %s in residue %d, atom  %s in residue %d and atom %s in residue %d "\
                  %(aAng["atom_id_1"], aAng["atom_id_1_resNum"], aAng["atom_id_2"], aAng["atom_id_2_resNum"],\
                    aAng["atom_id_3"], aAng["atom_id_3_resNum"]))
   
        tLinkedObj.cLink["torsions"] = []
        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
            if (aTor["atom_id_2_alias"] == atm1 and aTor["atom_id_3_alias"] == atm2) or\
               (aTor["atom_id_2_alias"] == atm2 and aTor["atom_id_3_alias"] == atm1):
                tLinkedObj.cLink["torsions"].append(aTor)
        print("Number of Link Torsions ", len(tLinkedObj.cLink["torsions"]))
        print("They are :")
        for aTor in tLinkedObj.cLink["torsions"]:
            print("Torsion between atom %s in residue %d, atom %s in residue %d, atom %s in residue %d, and atom %s in residue %d "\
                 %(aTor["atom_id_1"],aTor["atom_id_1_resNum"], aTor["atom_id_2"], aTor["atom_id_2_resNum"],\
                   aTor["atom_id_3"],aTor["atom_id_3_resNum"], aTor["atom_id_4"], aTor["atom_id_4_resNum"]))

        
        tLinkedObj.cLink["chirals"] = [] 
        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
                if aChi["atom_id_centre_alias"] ==atm1:
                    if len(tLinkedObj.stdLigand1["linkChir"])==0:
                        aChi["volume_sign"] = "both"
                    tLinkedObj.cLink["chirals"].append(aChi)
                elif aChi["atom_id_centre_alias"] ==atm2:
                    if len(tLinkedObj.stdLigand2["linkChir"])==0:
                        aChi["volume_sign"] = "both"
                    tLinkedObj.cLink["chirals"].append(aChi)
            print("Number of Link chirals ", len(tLinkedObj.cLink["chirals"]))
            print("They are :")
            for aChi in tLinkedObj.cLink["chirals"]:
                print("Chiral centre on atom %s in residue %d "%(aChi["atom_id_centre"], aChi["atom_id_centre_resNum"]))
                print("other atoms are : ")
                print("atom %s in residue %d"%(aChi["atom_id_1"], aChi["atom_id_1_resNum"]))
                print("atom %s in residue %d"%(aChi["atom_id_2"], aChi["atom_id_2_resNum"]))
                print("atom %s in residue %d"%(aChi["atom_id_3"], aChi["atom_id_3_resNum"]))
                print("Volume sign ", aChi["volume_sign"])

        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]:
            if len(list(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys())) !=0:
                tLinkedObj.cLink["planes"] = {}
                for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys()): 
                    inAtmIds = []
                    for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                        if aPlAtm["atom_id_alias"] == atm1 or  aPlAtm["atom_id_alias"]==atm2:
                            inAtmIds.append(aPlAtm["atom_id_alias"])
                    if atm1 in inAtmIds and atm2 in inAtmIds:      
                        tLinkedObj.cLink["planes"][aPl] = [] 
                        for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                            tLinkedObj.cLink["planes"][aPl].append(aPlAtm)
        if "planes" in tLinkedObj.cLink:
            print("Number of Link planes ", len(list(tLinkedObj.cLink["planes"].keys())))
            for aPl in sorted(tLinkedObj.cLink["planes"].keys()):
                print("Plane %s contains %d atoms. They are: "%(aPl, len(tLinkedObj.cLink["planes"][aPl])))
                for aAtm in tLinkedObj.cLink["planes"][aPl]:
                    print("Plane %s, atom  %s in residue %d, dist_esd  %s "%(aAtm["plane_id"], aAtm["atom_id"], aAtm["atom_id_resNum"], aAtm["dist_esd"]))

       
    def outOneLinkInfo(self, tLinkedObj):

        #aOutCifName = tLinkedObj.combLigand["name"] + "_link.cif"
        aOutCifName = self.outRoot + "_link.cif"
        
        try: 
            aOutCif = open(aOutCifName, "w")
        except IOError:
            self.errLevel    = 41
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for writting "%aOutCifName)
        else:
            self.outVerInfo(aOutCif)
            self.outCompList(aOutCif, tLinkedObj)
            self.outModList(aOutCif, tLinkedObj)
            self.outLinkList(aOutCif, tLinkedObj)
            self.outAllComps(aOutCif,  tLinkedObj) 
            self.outAllMods(aOutCif,  tLinkedObj) 
            self.outAllLinks(aOutCif, tLinkedObj.cLink) 
      
            aOutCif.close()

        aOutPDBName = self.outRoot + "_fullMol.pdb"

        #for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]:
            
  
    def outVerInfo(self, tOutFile):

        tOutFile.write("#%s%s\n"%("ACEDRG_VERSION".ljust(30),    self.verInfo["ACEDRG_VERSION"].ljust(20)))
        tOutFile.write("#%s%s\n"%("ACEDRG_DB_VERSION".ljust(30), self.verInfo["DATABASE_VERSION"].ljust(20)))
        tOutFile.write("#%s%s\n"%("RDKit_VERSION".ljust(30),    self.verInfo["RDKit_VERSION"].ljust(20)))
        tOutFile.write("#%s%s\n"%("REFMAC_VERSION".ljust(30),   self.verInfo["REFMAC_VERSION"].ljust(20)))
        tOutFile.write("#\n\n")
        

    def outCompList(self, tOutFile, tLinkedObj):

        if tLinkedObj.stdLigand1["compOut"] or tLinkedObj.stdLigand2["compOut"]: 
            tOutFile.write("data_comp_list\n\n")
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_comp.id\n")
            tOutFile.write("_chem_comp.three_letter_code\n")
            tOutFile.write("_chem_comp.name\n")
            tOutFile.write("_chem_comp.group\n")
            tOutFile.write("_chem_comp.number_atoms_all\n")
            tOutFile.write("_chem_comp.number_atoms_nh\n")
        
        if tLinkedObj.stdLigand1["compOut"] :    
            aNL = tLinkedObj.stdLigand1["name"]
            aN3 =""
            if len(aNL) >= 3:
                aN3 = tLinkedObj.stdLigand1["name"][:3]
            else:
                aN3 = aNL
            aName   = tLinkedObj.stdLigand1["list"]["name"]
            aNameL  = len(aName) + 6
            aGrp    = tLinkedObj.stdLigand1["list"]["group"]
            NA      = tLinkedObj.stdLigand1["list"]["number_atoms_all"]
            NAH     = tLinkedObj.stdLigand1["list"]["number_atoms_nh"]
            aL="%s%s%s%s%s%s\n"%(aNL.ljust(10), aN3.ljust(10), aName.ljust(aNameL), aGrp.ljust(20), NA.ljust(10), NAH)
            tOutFile.write(aL) 
        
        if tLinkedObj.stdLigand2["compOut"]:          
            aNL = tLinkedObj.stdLigand2["name"]
            aN3 =""
            if len(aNL) >= 3:
                aN3 = tLinkedObj.stdLigand2["name"][:3]
            else:
                aN3 = aNL
        
            aName   = tLinkedObj.stdLigand2["list"]["name"]
            aNameL  = len(aName) + 6
            aGrp    = tLinkedObj.stdLigand2["list"]["group"]
            NA      = tLinkedObj.stdLigand2["list"]["number_atoms_all"]
            NAH     = tLinkedObj.stdLigand2["list"]["number_atoms_nh"]
            aL="%s%s%s%s%s%s\n"%(aNL.ljust(10), aN3.ljust(10), aName.ljust(aNameL), aGrp.ljust(20), NA.ljust(10), NAH)
            tOutFile.write(aL) 
        if tLinkedObj.stdLigand1["compOut"] or tLinkedObj.stdLigand2["compOut"]: 
            tOutFile.write("\n")
                 
    def outModList(self, tOutFile, tLinkedObj):

        tOutFile.write("data_mod_list\n\n")
        tOutFile.write("loop_\n")
        tOutFile.write("_chem_mod.id\n")
        tOutFile.write("_chem_mod.name\n")
        tOutFile.write("_chem_mod.comp_id\n")
        tOutFile.write("_chem_mod.group_id\n")
            
        aMN     = tLinkedObj.modLigand1["name"]
        aName   = tLinkedObj.stdLigand1["list"]["name"]
        aNameL  = len(aName) + 6
        aLN = tLinkedObj.stdLigand1["name"]
        aGrp    = tLinkedObj.stdLigand1["list"]["group"]
        aL="%s%s%s%s\n"%(aMN.ljust(10), aName.ljust(aNameL), aLN.ljust(10), ".".ljust(20))
        tOutFile.write(aL) 
        
        aMN     = tLinkedObj.modLigand2["name"]
        aName   = tLinkedObj.stdLigand2["list"]["name"]
        aNameL  = len(aName) + 6
        aLN     = tLinkedObj.stdLigand2["name"]
        aGrp    = tLinkedObj.stdLigand2["list"]["group"]
        aL="%s%s%s%s\n"%(aMN.ljust(10), aName.ljust(aNameL), aLN.ljust(10), ".".ljust(20))
        tOutFile.write(aL) 
        
        tOutFile.write("\n")

    def outLinkList(self, tOutFile, tLinkedObj):

        tOutFile.write("data_link_list\n\n")
        tOutFile.write("loop_\n")
        tOutFile.write("_chem_link.id\n")
        tOutFile.write("_chem_link.comp_id_1\n")
        tOutFile.write("_chem_link.mod_id_1\n")
        tOutFile.write("_chem_link.group_comp_1\n")
        tOutFile.write("_chem_link.comp_id_2\n")
        tOutFile.write("_chem_link.mod_id_2\n")
        tOutFile.write("_chem_link.group_comp_2\n")
        tOutFile.write("_chem_link.name\n")

        aLID  = tLinkedObj.cLink["name"]
        aLN1  = tLinkedObj.stdLigand1["name"]
        aLMN1 = tLinkedObj.modLigand1["name"]
        aG1   = self.setGroupId(tLinkedObj.stdLigand1["list"]["group"])
        print("Ligand 1 group id %s "%aG1)
        #aG1   = "."
        aLN2  = tLinkedObj.stdLigand2["name"]
        aLMN2 = tLinkedObj.modLigand2["name"]
        aG2   = self.setGroupId(tLinkedObj.stdLigand2["list"]["group"])
        print("Ligand 2 group id %s "%aG2)
        #aG2   = "."
        aL="%s%s%s%s%s%s%s%s\n"%(aLID.ljust(15), aLN1.ljust(10), aLMN1.ljust(12), aG1.ljust(20),\
                                 aLN2.ljust(10), aLMN2.ljust(12), aG2.ljust(20), aLID.ljust(15))
        tOutFile.write(aL) 
         
        tOutFile.write("\n")

    def setGroupId(self, tId):

        retId = tId.lower()

        if retId.find("peptide") != -1 and retId.find("-") !=-1:
            strs = retId.strip().split("-")
            if len(strs)==2:
                retId = strs[0].upper() + "-" + strs[1]

        return retId

        

    def outAllComps(self, tOutFile, tLinkedObj):
        if tLinkedObj.stdLigand1["dataBlock"] and tLinkedObj.stdLigand1["outComp"] and tLinkedObj.stdLigand1["compOut"]:
            self.outOneComp(tOutFile, tLinkedObj.stdLigand1)
        if tLinkedObj.stdLigand2["dataBlock"] and tLinkedObj.stdLigand2["outComp"] and tLinkedObj.stdLigand2["compOut"]:
            self.outOneComp(tOutFile, tLinkedObj.stdLigand2)

    def outOneComp(self, tOutFile, tMonomer):
        
        tOutFile.write("data_comp_%s\n\n"%tMonomer["name"])
        for aL in tMonomer["dataBlock"]:
            tOutFile.write(aL)
        tOutFile.write("\n")

    def outAllMods(self, tOutFile, tLinkedObj):

        aPoolAtoms1 = []
        if "atomName" in tLinkedObj.stdLigand1:
            aPoolAtoms1.append(tLinkedObj.stdLigand1["atomName"])
        self.outOneMod(tOutFile, tLinkedObj.modLigand1, tLinkedObj.describLevel, aPoolAtoms1)
        
        aPoolAtoms2 = []
        if "atomName" in tLinkedObj.stdLigand2:
            aPoolAtoms2.append(tLinkedObj.stdLigand2["atomName"])
        self.outOneMod(tOutFile, tLinkedObj.modLigand2, tLinkedObj.describLevel, aPoolAtoms2)

        
    def outOneMod(self, tOutFile, tModLigand, tLevel, tPoolAtoms):

        tOutFile.write("data_mod_%s\n\n"%tModLigand["name"])

        # Atoms, including all changed atoms whatever the description level is
        nDA = len(tModLigand["deleted"]["atoms"]) 
        nCA = len(tModLigand["changed"]["atoms"])
        nAA = len(tModLigand["added"]["atoms"])
 
        if nDA != 0 or nCA != 0 or nAA !=0 :
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_atom.mod_id\n")
            tOutFile.write("_chem_mod_atom.function\n")
            tOutFile.write("_chem_mod_atom.atom_id\n")
            tOutFile.write("_chem_mod_atom.new_atom_id\n")
            tOutFile.write("_chem_mod_atom.new_type_symbol\n")
            tOutFile.write("_chem_mod_atom.new_type_energy\n")
            tOutFile.write("_chem_mod_atom.new_charge\n")

            if nDA !=0:
                for aAtom in tModLigand["deleted"]["atoms"]:
                    tPoolAtoms.append(aAtom["atom_id"])
                    aC = "."
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    #for aKey in aAtom.keys():
                    #    print aKey, " : ", aAtom[aKey] 
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), aAtom["atom_id"].ljust(10),\
                                             ".".ljust(10), aAtom["type_symbol"].ljust(10), aAtom["type_energy"].ljust(10),\
                                             aC.ljust(10))
                    tOutFile.write(aL)

            if nCA !=0:
                for aAtom in tModLigand["changed"]["atoms"]:
                    aC = "."
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aAtom["atom_id"].ljust(10),\
                                             ".".ljust(10), aAtom["type_symbol"].ljust(10), aAtom["type_energy"].ljust(10),\
                                             aC.ljust(10))
                    tOutFile.write(aL)

            if nAA !=0:
                for aAtom in tModLigand["added"]["atoms"]:
                    aC = "."
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), aAtom["atom_id"].ljust(10),\
                                             ".".ljust(10), aAtom["type_symbol"].ljust(10), aAtom["type_energy"].ljust(10),\
                                             aC.ljust(10))
                    tOutFile.write(aL)
           
            tOutFile.write("\n")
                
        # Bonds, depending on  the description level 

        nDB  = len(tModLigand["deleted"]["bonds"])
        nCB  = 0
        CB_Bonds = []
        nAB = len(tModLigand["added"]["bonds"])

        if len(tModLigand["changed"]["bonds"]) !=0:
            for aBond in tModLigand["changed"]["bonds"]:
                if tLevel ==1:
                    if aBond["atom_id_1"] in tPoolAtoms and not aBond["atom_id_2"] in tPoolAtoms:
                        CB_Bonds.append(aBond)
                    elif aBond["atom_id_2"] in tPoolAtoms and not aBond["atom_id_1"] in tPoolAtoms:
                        CB_Bonds.append(aBond)
            nCB = len(CB_Bonds)
         
        if nDB !=0 or nCB !=0 or nAB !=0:
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_bond.mod_id\n")
            tOutFile.write("_chem_mod_bond.function\n")
            tOutFile.write("_chem_mod_bond.atom_id_1\n")
            tOutFile.write("_chem_mod_bond.atom_id_2\n")
            tOutFile.write("_chem_mod_bond.new_type\n")
            tOutFile.write("_chem_mod_bond.new_value_dist\n")
            tOutFile.write("_chem_mod_bond.new_value_dist_esd\n")               
         
            if nDB !=0:
                for aBond in tModLigand["deleted"]["bonds"]:
                    aBT = aBond["type"].lower()
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), ".".ljust(15),\
                                             ".".ljust(10))
                    tOutFile.write(aL)
                    
            
            if nCB !=0:
                #print tModLigand["name"]
                for aBond in CB_Bonds:
                    aBT = aBond["type"].lower()
                    #print aBond["atom_id_1"]
                    #print aBond["atom_id_2"]
                    #print aBond["value_dist"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), aBond["value_dist"].ljust(15),\
                                             aBond["value_dist_esd"].ljust(10))
                    tOutFile.write(aL)
                    
            if nAB !=0:
                for aBond in tModLigand["added"]["bonds"]:
                    aBT = aBond["type"].lower()
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), str(aBond["value_dist"]).ljust(15),\
                                             str(aBond["value_dist_esd"]).ljust(10))
                    tOutFile.write(aL)
                    
            
            tOutFile.write("\n")

        # Angles, depending on  the description level 

        nDAngs  = len(tModLigand["deleted"]["angles"])
        nCAngs  = 0
        CA_Angs = []
        nAAngs  = len(tModLigand["added"]["angles"])
        
        if len(tModLigand["changed"]["angles"]) != 0:
            for aAng in tModLigand["changed"]["angles"]:
                if tLevel ==1:
                    if aAng["atom_id_2"] in tPoolAtoms and not aAng["atom_id_1"] in tPoolAtoms\
                        and not aAng["atom_id_3"] in tPoolAtoms:
                         CA_Angs.append(aAng)
            nCAngs = len(CA_Angs)
        
        if nDAngs !=0 or nCAngs !=0 or nAAngs !=0:
            
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_angle.mod_id\n")
            tOutFile.write("_chem_mod_angle.function\n")
            tOutFile.write("_chem_mod_angle.atom_id_1\n")
            tOutFile.write("_chem_mod_angle.atom_id_2\n")
            tOutFile.write("_chem_mod_angle.atom_id_3\n")
            tOutFile.write("_chem_mod_angle.new_value_angle\n")
            tOutFile.write("_chem_mod_angle.new_value_angle_esd\n")
   
            if nDAngs !=0 :
                for aAng in tModLigand["deleted"]["angles"]:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), aAng["atom_id_1"].ljust(10),\
                                            aAng["atom_id_2"].ljust(10), aAng["atom_id_3"].ljust(10),\
                                            ".".ljust(15), ".".ljust(15))
                    tOutFile.write(aL)
            
            if nCAngs !=0 :
                for aAng in CA_Angs:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aAng["atom_id_1"].ljust(10),\
                                            aAng["atom_id_2"].ljust(10), aAng["atom_id_3"].ljust(10),\
                                            aAng["value_angle"].ljust(15), aAng["value_angle_esd"].ljust(15))
                    tOutFile.write(aL)
            
            if nAAngs !=0 :
                for aAng in tModLigand["added"]["angles"]:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), aAng["atom_id_1"].ljust(10),\
                                            aAng["atom_id_2"].ljust(10), aAng["atom_id_3"].ljust(10),\
                                            aAng["value_angle"].ljust(15), aAng["value_angle_esd"].ljust(15))
                    tOutFile.write(aL)
              
            tOutFile.write("\n")
            
        # Torsions, depending on  the description level 
        
        nDTors  = len(tModLigand["deleted"]["tors"])
        nCTors  = 0
        CT_Tors = []
        nATors  = len(tModLigand["added"]["tors"])

        if len(tModLigand["changed"]["tors"]) !=0:
            for aTor in tModLigand["changed"]["tors"]:
                if tLevel == 1:
                    if aTor["atom_id_2"] in tPoolAtoms and not aTor["atom_id_1"] in tPoolAtoms\
                       and not aTor["atom_id_3"] in tPoolAtoms and not aTor["atom_id_4"] in tPoolAtoms:
                        CT_Tors.append(aTor)
                    elif aTor["atom_id_3"] in tPoolAtoms and not aTor["atom_id_1"] in tPoolAtoms\
                       and not aTor["atom_id_2"] in tPoolAtoms and not aTor["atom_id_4"] in tPoolAtoms:
                        CT_Tors.append(aTor)
            nCTors = len(CT_Tors)

             
        if nDTors !=0 or nCTors != 0 or nATors !=0:

            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_tor.mod_id\n")
            tOutFile.write("_chem_mod_tor.function\n")
            tOutFile.write("_chem_mod_tor.atom_id_1\n")
            tOutFile.write("_chem_mod_tor.atom_id_2\n")
            tOutFile.write("_chem_mod_tor.atom_id_3\n")
            tOutFile.write("_chem_mod_tor.atom_id_4\n")
            tOutFile.write("_chem_mod_tor.new_value_angle\n")
            tOutFile.write("_chem_mod_tor.new_value_angle_esd\n")
   
            if nDTors !=0: 
                for aTor in tModLigand["deleted"]["tors"]:
                    aL ="%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              ".".ljust(15), ".".ljust(15))                  
                    tOutFile.write(aL)
           
            if nCTors !=0: 
                for aTor in CT_Tors:
                    aL ="%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15))                  
                    tOutFile.write(aL)
   
            if nATors !=0: 
                for aTor in tModLigand["added"]["tors"]:
                    aL ="%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15))                  
                    tOutFile.write(aL)
   
            tOutFile.write("\n")
    
        # Chiral centers, depending on  the description level 

        nDChirs  = len(tModLigand["deleted"]["chirs"])
        nCChirs  = 0
        CT_Chirs = []
        nAChirs  = len(tModLigand["added"]["chirs"])

        if len(tModLigand["changed"]["chirs"]) !=0:
            for aChi in tModLigand["changed"]["chirs"]:
                if tLevel == 1:
                    if not aChi["atom_id_centre"] in tPoolAtoms and (aChi["atom_id_1"] in tPoolAtoms\
                       or  aChi["atom_id_2"] in tPoolAtoms or aChi["atom_id_3"] in tPoolAtoms):
                        CT_Chirs.append(aChi)
            nCChirs = len(CT_Chirs)

        if nDChirs !=0 or nCChirs != 0 or nAChirs !=0:

            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_chir.mod_id\n")
            tOutFile.write("_chem_mod_chir.function\n")
            tOutFile.write("_chem_mod_chir.atom_id_centre\n")
            tOutFile.write("_chem_mod_chir.atom_id_1\n")
            tOutFile.write("_chem_mod_chir.atom_id_2\n")
            tOutFile.write("_chem_mod_chir.atom_id_3\n")
            tOutFile.write("_chem_mod_chir.new_volume_sign\n")
   
            if nDChirs !=0: 
                for aChi in tModLigand["deleted"]["chirs"]:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), \
                                            aChi["atom_id_centre"].ljust(10), aChi["atom_id_1"].ljust(10),\
                                            aChi["atom_id_2"].ljust(10), aChi["atom_id_3"].ljust(10),\
                                              ".".ljust(15))                  
                    tOutFile.write(aL)
           
            if nCChirs !=0: 
                for aChi in CT_Chirs:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), \
                                            aChi["atom_id_centre"].ljust(10), aChi["atom_id_1"].ljust(10),\
                                            aChi["atom_id_2"].ljust(10), aChi["atom_id_3"].ljust(10),\
                                            aChi["volume_sign"].ljust(15))                  
                    tOutFile.write(aL)
           
            if nAChirs !=0: 
                for aChi in tModLigand["added"]["chirs"]:
                    aL ="%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), \
                                            aChi["atom_id_centre"].ljust(10), aChi["atom_id_1"].ljust(10),\
                                            aChi["atom_id_2"].ljust(10), aChi["atom_id_3"].ljust(10),\
                                            aChi["volume_sign"].ljust(15))                  
                    tOutFile.write(aL)
           
            tOutFile.write("\n")

        # Planes

        nDPls  = len(tModLigand["deleted"]["planes"])
        nCPls  = len(tModLigand["changed"]["planes"])
        nAPls  = len(tModLigand["added"]["planes"])
         
        if nDPls or nCPls or nAPls:
     
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod_plane.mod_id\n")
            tOutFile.write("_chem_mod_plane.function\n")
            tOutFile.write("_chem_mod_plane.plane_id\n")
            tOutFile.write("_chem_mod_plane.atom_id\n")
            tOutFile.write("_chem_mod_plane.dist_esd\n")

            if nDPls: 
                for aPl in tModLigand["deleted"]["planes"]:
                    for aPlAtm in aPl:
                        aL ="%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15),\
                            aPlAtm["plane_id"].ljust(15), aPlAtm["atom_id"].ljust(15), aPlAtm["dist_esd"]) 
                             
                        tOutFile.write(aL)

        tOutFile.write("\n")

    def outAllLinks(self, tOutFile, tLink):

        tOutFile.write("data_link_%s\n\n"%tLink["name"])

        if len(tLink["bonds"]) !=0:
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_link_bond.link_id\n")
            tOutFile.write("_chem_link_bond.atom_1_comp_id\n")
            tOutFile.write("_chem_link_bond.atom_id_1\n")
            tOutFile.write("_chem_link_bond.atom_2_comp_id\n")
            tOutFile.write("_chem_link_bond.atom_id_2\n")
            tOutFile.write("_chem_link_bond.type\n")
            tOutFile.write("_chem_link_bond.value_dist\n")
            tOutFile.write("_chem_link_bond.value_dist_esd\n")
            for aBond in tLink["bonds"]:
                aL = "%s%s%s%s%s%s%s%s\n"%(tLink["name"].ljust(10), str(aBond["atom_id_1_resNum"]).ljust(10), aBond["atom_id_1"].ljust(10),\
                                           str(aBond["atom_id_2_resNum"]).ljust(10), aBond["atom_id_2"].ljust(10), aBond["type"].ljust(10),\
                                           aBond["value_dist"].ljust(12), aBond["value_dist_esd"].ljust(12))
                tOutFile.write(aL)
            tOutFile.write("\n")

        if "angles" in tLink and len(tLink["angles"]) !=0:
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_link_angle.link_id\n")
            tOutFile.write("_chem_link_angle.atom_1_comp_id\n")
            tOutFile.write("_chem_link_angle.atom_id_1\n")
            tOutFile.write("_chem_link_angle.atom_2_comp_id\n")
            tOutFile.write("_chem_link_angle.atom_id_2\n")
            tOutFile.write("_chem_link_angle.atom_3_comp_id\n")
            tOutFile.write("_chem_link_angle.atom_id_3\n")
            tOutFile.write("_chem_link_angle.value_angle\n")
            tOutFile.write("_chem_link_angle.value_angle_esd\n")

            for aAng in tLink["angles"]: 
                aL="%s%s%s%s%s%s%s%s%s\n"%(tLink["name"].ljust(10),str(aAng["atom_id_1_resNum"]).ljust(10),  aAng["atom_id_1"].ljust(10),\
                                           str(aAng["atom_id_2_resNum"]).ljust(10),  aAng["atom_id_2"].ljust(10),\
                                           str(aAng["atom_id_3_resNum"]).ljust(10),  aAng["atom_id_3"].ljust(10),\
                                           aAng["value_angle"].ljust(15), aAng["value_angle_esd"].ljust(10))
                tOutFile.write(aL)
            tOutFile.write("\n")

        if "torsions" in tLink and len(tLink["torsions"]) !=0:
  
            self.setTorIdsInOneLink(tLink["torsions"])
 
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_link_tor.link_id\n")
            tOutFile.write("_chem_link_tor.id\n")
            tOutFile.write("_chem_link_tor.atom_1_comp_id\n")
            tOutFile.write("_chem_link_tor.atom_id_1\n")
            tOutFile.write("_chem_link_tor.atom_2_comp_id\n")
            tOutFile.write("_chem_link_tor.atom_id_2\n")
            tOutFile.write("_chem_link_tor.atom_3_comp_id\n")
            tOutFile.write("_chem_link_tor.atom_id_3\n")
            tOutFile.write("_chem_link_tor.atom_4_comp_id\n")
            tOutFile.write("_chem_link_tor.atom_id_4\n")
            tOutFile.write("_chem_link_tor.value_angle\n")
            tOutFile.write("_chem_link_tor.value_angle_esd\n")
            tOutFile.write("_chem_link_tor.period\n")
          
            for aTor in tLink["torsions"]: 
                aL="%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%(tLink["name"].ljust(10), aTor["id"].ljust(16),\
                                           str(aTor["atom_id_1_resNum"]).ljust(10),  aTor["atom_id_1"].ljust(10),\
                                           str(aTor["atom_id_2_resNum"]).ljust(10),  aTor["atom_id_2"].ljust(10),\
                                           str(aTor["atom_id_3_resNum"]).ljust(10),  aTor["atom_id_3"].ljust(10),\
                                           str(aTor["atom_id_4_resNum"]).ljust(10),  aTor["atom_id_4"].ljust(10),\
                                           aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(10), aTor["period"])
                tOutFile.write(aL)
            tOutFile.write("\n")

        if "chirals" in tLink and  len(tLink["chirals"]) !=0:
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_link_chir.link_id\n")
            tOutFile.write("_chem_link_chir.atom_centre_comp_id\n")
            tOutFile.write("_chem_link_chir.atom_id_centre\n")
            tOutFile.write("_chem_link_chir.atom_1_comp_id\n")
            tOutFile.write("_chem_link_chir.atom_id_1\n")
            tOutFile.write("_chem_link_chir.atom_2_comp_id\n")
            tOutFile.write("_chem_link_chir.atom_id_2\n")
            tOutFile.write("_chem_link_chir.atom_3_comp_id\n")
            tOutFile.write("_chem_link_chir.atom_id_3\n")
            tOutFile.write("_chem_link_chir.volume_sign\n")
          
            for aChi in tLink["chirals"]:
                aL="%s%s%s%s%s%s%s%s%s%s\n"%(tLink["name"].ljust(10),str(aChi["atom_id_centre_resNum"]).ljust(10),\
                                           aChi["atom_id_centre"].ljust(10), str(aChi["atom_id_1_resNum"]).ljust(10),\
                                           aChi["atom_id_1"].ljust(10), str(aChi["atom_id_2_resNum"]).ljust(10),\
                                           aChi["atom_id_2"].ljust(10), str(aChi["atom_id_3_resNum"]).ljust(10),\
                                           aChi["atom_id_3"].ljust(10), aChi["volume_sign"].ljust(10))
                tOutFile.write(aL)
            tOutFile.write("\n")

        if "planes" in tLink and  len(list(tLink["planes"].keys())) !=0:
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_link_plane.link_id\n")
            tOutFile.write("_chem_link_plane.plane_id\n")
            tOutFile.write("_chem_link_plane.atom_comp_id\n")
            tOutFile.write("_chem_link_plane.atom_id\n")
            tOutFile.write("_chem_link_plane.dist_esd\n")
            for aPl in list(tLink["planes"].keys()):
                for aAtm in tLink["planes"][aPl]:
                    aL="%s%s%s%s%s\n"%(tLink["name"].ljust(10), aAtm["plane_id"].ljust(10),\
                                       str(aAtm["atom_id_resNum"]).ljust(10),\
                                           aAtm["atom_id"].ljust(10), aAtm["dist_esd"].ljust(10))
                    tOutFile.write(aL)
            tOutFile.write("\n")
                     
            
            

