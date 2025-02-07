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
import time
import math
import shutil

#if os.name != 'nt':
#    import fcntl
#import signal

from . exebase       import CExeCode

from . filetools     import FileTransformer
from . filetools     import Ccp4MmCifObj
from . chem          import ChemCheck
from . metalMode     import metalMode

from . utility       import isInt
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
        self.stdLigand1["origCharge"]       = {}
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
        self.stdLigand2["origCharge"]       = {}
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
        self.delSectionsBonds                 = []
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
        print("error info", self.errMessage)
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

class MondRes(object):

    def __init__(self):
        
        #keyWordList = ["RES-NAME", "FILE"(optional), "ADD",
        #                "DELETE", "CHANGE"]

        self.describLevel                   = 1

        self.stdLigand                     = {}
        self.stdLigand["fromScr"]          = False
        self.stdLigand["userIn"]           = False
        self.stdLigand["compOut"]          = False
        self.stdLigand["inCif"]            = False
        self.stdLigand["dataBlock"]        = None
        self.stdLigand["comp"]             = None
        self.stdLigand["hAtom"]            = []
        self.stdLigand["outComp"]          = False
        self.stdLigand["remainAtoms"]      = []
        self.stdLigand["remainIdxHs"]      = []
        self.stdLigand["remainBonds"]      = []
        self.stdLigand["remainAngs"]       = []
        self.stdLigand["remainTors"]       = []
        self.stdLigand["remainChirs"]      = []
        self.stdLigand["linkChir"]         = []       # It will be [aChi, delAtmIdx]
        self.stdLigand["remainPls"]        = []

        self.newLigand                       = {}
        
        
        self.modLigand                       = {}
        self.modLigand["outComp"]            = True
        self.modLigand["deleted"]            = {}
        self.modLigand["deleted"]["atoms"]   = []
        self.modLigand["deleted"]["bonds"]   = []
        self.modLigand["deleted"]["angles"]  = []
        self.modLigand["deleted"]["tors"]    = []
        self.modLigand["deleted"]["chirs"]   = []
        self.modLigand["deleted"]["planes"]  = []
        self.modLigand["changed"]            = {}
        self.modLigand["changed"]["atoms"]   = []
        self.modLigand["changed"]["bonds"]   = []
        self.modLigand["changed"]["angles"]  = []
        self.modLigand["changed"]["tors"]    = []
        self.modLigand["changed"]["chirs"]   = []
        self.modLigand["changed"]["planes"]  = []
        #self.modLigand["changed"]["formal_charges"] = []
        self.modLigand["changed"]["charges"] = []
        self.modLigand["added"]              = {}
        self.modLigand["added"]["atoms"]     = []
        self.modLigand["added"]["bonds"]     = []
        self.modLigand["added"]["angles"]    = []
        self.modLigand["added"]["tors"]      = []
        self.modLigand["added"]["chirs"]     = []
        self.modLigand["added"]["planes"]    = []
        
        self.delSectionsAtoms                = []
        self.delSectionsBonds                = []
        self.errLevel                        = 0
        self.errMessage                      = []
        
        self.dataHead                        = {}


class CovLinkGenerator(CExeCode):

    def __init__(self, tSetParas, tInstructions, tScrDir, tOutRoot, tTabLoc, tVerInfo=None, tTestMode=None):

        
        if tTestMode:
            self.testMode = True
        else:
            self.testMode = False
        self.verInfo          = tVerInfo
        #print self.verInfo

        self.workMode         = 0   # Mode means
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
        
        self.lMode = 0                        # 0 No job, 
                                              # 1 modificate a monomer
                                              # 2 generate a link
        self.lNoMeta = True

        # TEMPO
        if "CCP4" not in os.environ:
            print("You need to install or setup CCP4 first ")
            sys.exit() 
       
        self.allChemCombDir    = os.getenv("CLIBD_MON")
        #self.allChemCombDir   = "/Users/flong/DB/PDB_related/PDB_Ligands/Cif"

        self.setParas         = tSetParas 
        self.scrDir           = tScrDir
        self.subRoot          = ""
        self.outRoot          = tOutRoot
        
        self.exeAcedrg        = "acedrg " 
        #self.exeAcedrg       = "/lmb/home/flong/workplace/LIBMOL/bin/acedrg"
        
        self.instructions     = tInstructions
        
        # Links related 
        self.linkInstructionsContent   = []
        self.ligSrcDir        = ""

        self.cLinks           = []

        #Modification of monomers related 
        
        self.modInstructionsContent = []
        self.modSrcDir        = ""
        self.modRess          = []

        self.chemCheck        = ChemCheck()
        self.fileTool         = FileTransformer()
        self.metalMode        = metalMode(tTabLoc)


        if os.path.isfile(self.instructions):
            # input from a file
            self.getInstructions()
        else:
            print("%s can not be found for reading"%self.instructions)
            self.errLevel = 1
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("Instruction file %s can not be found for reading"%self.instructions)
        
        if self.lMode ==1:
            if not self.errLevel:
                if len(self.modRess):
                    print("Number of modifcations to be processed ", len(self.modRess))
                    for aMod in self.modRess:
                        self.processOneModRes(aMod)
            else:
                self.errLevel  = 5
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append("Errors are found in the instruction file")
                self.errMessage[self.errLevel].append("Check your instruction file")
                
        elif self.lMode ==2:
            
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
            
    def getInstructions(self):

        #aFullName=os.path.basename(self.linkInstructions)
        #aExt     = ""
        #if aFullName.find(".") !=-1:
        #    aExt = aFullName.strip().split(".")[-1].strip()
        #    if aExt.find("cif") !=-1:
        #        self.getInstructionsForLinkFromCif()
        #    else:
        #        # Tempo comment off free format at the moment
        try :
            aInsF = open(self.instructions, "r")
        except IOError:
            print("% can not be opened for reading ! "%self.instructions)
            self.errLevel = 11
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("% can not be opened for reading ! "%self.instructions)
        else:
            allLs = aInsF.readlines()
            aInsF.close()
            insLines = []
            
            for aL in allLs:
                aL = aL.strip()
                if len(aL) !=0:
                    strGrp = aL.split()
                    if len(strGrp) > 0 and strGrp[0] !="#":
                        if strGrp[0].upper().find("MOD:")!=-1:
                            self.lMode = 1
                            insLines.append(aL) 
                        elif strGrp[0].upper().find("LINK:") !=-1:
                            self.lMode = 2
                            insLines.append(aL)
            
            if self.lMode == 1:
                self.getInstructionsForModification(insLines)
            elif self.lMode ==2:
                
                self.getInstructionsForLinkFreeFormat2(insLines)
            
             
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
        
        aInsObj = Ccp4MmCifObj(self.instructions)
        aLink = CovLink()
        print("Instruction file of cif format ", self.instructions)
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
            aInsF = open(self.instructions, "r")
        except IOError:
            print("% can not be opened for reading ! "%self.instructions)
            self.errLevel = 11
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("% can not be opened for reading ! "%self.instructions)
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
                print (aList[i])
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


    def getInstructionsForLinkFreeFormat2(self, tInLines):

        aList = []
           
        for aL in tInLines:
            strs = aL.strip().split()
            for aStr in strs:
                if aStr.find("LINK:") == -1:
                    aList.append(aStr)

        self.linkInstructionsContent = tInLines

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
                #print(aList[i])
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
                                aBond["comp_serial_num"] = aNum 
                                aBond["type"] = "SING"
                                aLink.delSectionsBonds.append(aBond)
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

                # Set Ligand 1 file  
                aNS = aLink.stdLigand1["name"][0].lower()
                aNL = aLink.stdLigand1["name"].upper()
                if not aLink.stdLigand1["userIn"]:
                    aLink.stdLigand1["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                    aLink.stdLigand1["userIn"] = True
                    
                # Look into the file to see if it contains metal element and then dealt with them
                aMmcifObj1 = Ccp4MmCifObj(aLink.stdLigand1["inCif"])
                aMmcifObj1.checkBlockCompsExist()
                #if not aMmcifObj1["errLevel"]:
                #    if aLink.stdLigand1["name"] in aMmcifObj1["ccp4CifObj"]["comps"]:
                #        self.lNoMetal=self.chemCheck.isOrganicInCif(aMmcifObj1["ccp4CifObj"]["comps"][aLink.stdLigand1["name"]]["atoms"]) 
                #        if not self.lNoMetal:    # contain metal 
                #            print(aLink.stdLigand1["name"], " contains metal atoms ")
                #            self.subRoot = aLink.stdLigand1["name"] + "_intmedia"
                #            self.setOneMonomer(aLink.stdLigand1, 4)
                #            aIntmedFN = self.subRoot + "_tmp.cif"
                #            if os.path.isfile(aIntmedFN):
                #                aLink.stdLigand1["inCif"]=aIntmedFN  
                #                print("one component cif is now",  aLink.stdLigand1["inCif"])
                #            #aIntMediatCif = os.path.join(self.subRoot + "_TMP", )
                #            #print("the intermediate file is ", )
                
                # Ligand 2 file and elements
                aNS = aLink.stdLigand2["name"][0].lower()
                aNL = aLink.stdLigand2["name"].upper()
                if not aLink.stdLigand2["userIn"]:
                    aLink.stdLigand2["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                    aLink.stdLigand2["userIn"] = True
                  
                aMmcifObj2 = Ccp4MmCifObj(aLink.stdLigand2["inCif"])
                aMmcifObj2.checkBlockCompsExist()
                """
                if not aMmcifObj2["errLevel"]:
                    if aLink.stdLigand2["name"] in aMmcifObj2["ccp4CifObj"]["comps"]:
                        self.lNoMetal=self.chemCheck.isOrganicInCif(aMmcifObj2["ccp4CifObj"]["comps"][aLink.stdLigand2["name"]]["atoms"]) 
                        if not self.lNoMetal:    # contain metal 
                            print(aLink.stdLigand2["name"], " contains metal atoms ")
                            self.subRoot = aLink.stdLigand2["name"] + "_intmedia"
                            self.setOneMonomer(aLink.stdLigand2, 4)
                            aIntmedFN = self.subRoot + "_tmp.cif"
                            if os.path.isfile(aIntmedFN):
                                aLink.stdLigand2["inCif"]=aIntmedFN  
                                print("one component cif is now",  aLink.stdLigand2["inCif"])
                
                """
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
        
        
    def getInstructionsForModification(self, tInLines):
        
        insList = []
        for aL in tInLines:
            strs = aL.strip().split()
            if len(strs) !=0 :
                for aStr in strs:
                    if aStr.upper().find("MOD:") == -1:
                        insList.append(aStr)

        self.modInstructionsContent = tInLines
        
        if len(insList):
            lDel = False
            lCh  = False
            lAd  = False
            i    = 0
            
            aModRes = MondRes()
            
            # Some default values. They will be changed at different stages
 
            aModRes.stdLigand["userIn"]  = False 
            aModRes.stdLigand["outComp"] = False 
            aModRes.stdLigand["outMod"]  = True 
            aModRes.stdLigand["group"]   = "."
            
            
            nL = len(insList)
            
            aKList = ["ATOM", "BOND", "CHARGE", "ADD", "DELETE", "CHANGE"]
            
            while i < nL :
                #print(insList[i])                    
                #print(lCh)
                #print(lDel)
                #print(lAd)
                if (i+1) < nL:
                    if insList[i].upper().find("RES-NAME") != -1:
                        aModRes.stdLigand["name"]   = insList[i+1]
                        aModRes.modLigand["name"]   = aModRes.stdLigand["name"] + "m"               
                        i +=2
                    elif insList[i].upper().find("FILE") != -1:
                        aModRes.stdLigand["inCif"]   = insList[i+1]
                        aModRes.stdLigand["userIn"]  =  True
                        aModRes.stdLigand["compOut"] =  True
                        i +=2
                    elif insList[i].upper().find("DELETE") != -1:
                        lDel = True
                        lCh  = False
                        lAd  = False
                        i +=1
                    elif insList[i].upper().find("CHANGE") != -1:
                        lCh  = True
                        lDel = False
                        lAd  = False
                        i +=1
                    elif insList[i].upper().find("ADD") != -1:
                        lAd  = True
                        lDel = False
                        lCh  = False
                        i +=1
                    elif lDel :
                        if insList[i].upper().find("ATOM") != -1 :
                            if i+1 < nL:
                                atmDName = setNameByNumPrime(insList[i+1])
                                aAtom = {}
                                aAtom["atom_id"] = atmDName
                                aModRes.modLigand["deleted"]["atoms"].append(aAtom) 
                                i+=2 
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE ATOM section of the instruction file. "
                                aMess+= "%s needs to be followed by  the atom id.\n"%insList[i]
                                self.errMessage[self.errLevel].append(aMess)
                                break
                        elif insList[i].upper().find("BOND") != -1 :
                            if i+2 < nL:
                                aBond = {}
                                aBond["atom_id_1"] = insList[i+1]
                                aBond["atom_id_2"] = insList[i+2]
                                aModRes.modLigand["deleted"]["bonds"].append(aBond)
                                i+=3
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Format error in DELECTE BOND section of the instruction file.\n"
                                aMess+= "%s needs to be followed by id of two atoms.\n"% insList[i]
                                self.errMessage[self.errLevel].append(aMess)
                                break
                    elif lCh :
                        if insList[i].upper().find("BOND") != -1 :
                            if i+3 < nL:
                                if not  insList[i+3] in aKList:
                                    aBond = {}
                                    aBond["atom_id_1"]      = insList[i+1]
                                    aBond["atom_id_2"]      = insList[i+2]
                                    aBond["type"]           = insList[i+3]
                                    aBond["value_order"]    = insList[i+3]
                                    aBond["value_dist"]     = 0.0     
                                    aBond["value_dist_esd"] = 0.02
                                    i+=4
                                    aModRes.modLigand["changed"]["bonds"].append(aBond)
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in BOND CHANGE section in the instruction file:\n "
                                    aMess+= "Not enough information to define the change of the bond\n"
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
                        elif insList[i].upper().find("CHARGE") != -1 :
                            if i+2 < nL:
                                aAddCharge = {}
                                aAddCharge["atom_id"]       = insList[i+1].strip()
                                insList[i+2] = insList[i+2].strip()
                                if isInt(insList[i+2]):
                                    aAddCharge["charge"] = int(insList[i+2])
                                    print("Charge is %d"%aAddCharge["charge"])
                                    aModRes.modLigand["changed"]["charges"].append(aAddCharge)
                                  
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in add CHARGE section in the instruction file:\n "
                                    aMess+= "Formal charge added must be an integer.\n"
                                    aMess+= "Your input is %s.\n"%insList[i+2]
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                                i+=3
                            else :
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in the formal charge CHANGE section in the instruction file:\n "
                                aMess+= "Not enough information to define the change of the formal charge\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break
                        
                    elif lAd :
                        if insList[i].upper().find("ATOM") != -1 :
                            if i+3 < nL:
                                if not insList[i+3] in aKList:
                                    
                                    aAtom = {}
                                    aAtom["atom_id"]     = insList[i+1]
                                    aAtom["type_symbol"] = insList[i+2]
                                    aAtom["charge"]      = int(insList[i+3])
                                    aAtom["type_energy"] = aAtom["type_symbol"]
                                    aAtom["comp_id"]     = aModRes.modLigand["name"]
                                    aAtom["is_added"]    = True
                                    aModRes.modLigand["added"]["atoms"].append(aAtom)
                                    i+=4
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in ADD ATOM section in the instruction file:\n "
                                    aMess+= "Not enough information to define added atom\n"
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
                            else:
                                self.errLevel = 2
                                if self.errLevel not in self.errMessage:
                                    self.errMessage[self.errLevel] = []
                                aMess = "Error in add ATOM section in the instruction file:\n "
                                aMess+= "Not enough information to define the added atom\n"
                                self.errMessage[self.errLevel].append(aMess)
                                break
                        elif insList[i].upper().find("BOND") != -1 :
                            if i+3 < nL:
                                if not insList[i+3] in aKList:
                                    aBond = {}
                                    aBond["atom_id_1"] = insList[i+1]
                                    aBond["atom_id_2"] = insList[i+2]
                                    aBond["type"]      = insList[i+3]
                                    aBond["value_order"] = insList[i+3]
                                    aBond["value_dist"]     = 0.0     
                                    aBond["value_dist_esd"] = 0.02
                                    aModRes.modLigand["added"]["bonds"].append(aBond)
                                    i+=4
                                else:
                                    self.errLevel = 2
                                    if self.errLevel not in self.errMessage:
                                        self.errMessage[self.errLevel] = []
                                    aMess = "Error in ADD BOND section in the instruction file:\n "
                                    aMess+= "Not enough information to define added bond\n"
                                    self.errMessage[self.errLevel].append(aMess)
                                    break
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
                        aMess = "Unknown keyword  %s. Check your instruction file\n"%insList[i]
                        self.errMessage[self.errLevel].append(aMess)
                        break
                
                else :
                    
                    print("Errors begin at %s. Check your instruction file "%insList[i])
                    self.errLevel = 2
                    if 2 not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("Errors begin at %s. Check your instruction file\n"%insList[i])
                    break
            
            if not self.errLevel:
                
                aNS = aModRes.stdLigand["name"][0].lower()
                aNL = aModRes.stdLigand["name"].upper()
                if not aModRes.stdLigand["userIn"]:
                    aModRes.stdLigand["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")
                    aModRes.stdLigand["userIn"] = True
                
                #self.checkMetalAndAetOneMonomer(aModRes.stdLigand)
                
                
                self.modRess.append(aModRes)
                if aModRes.stdLigand["userIn"]:
                    print("The input file is %s"%aModRes.stdLigand["inCif"])
                print("The modifcations to the monomer %s are : "%aModRes.stdLigand["name"])
                nda = len(aModRes.modLigand["deleted"]["atoms"])
                if nda > 0:
                    print("The following atoms are deleted.")
                    for aA in aModRes.modLigand["deleted"]["atoms"]:
                        print("Atom %s "%aA["atom_id"])
                ndb = len(aModRes.modLigand["deleted"]["bonds"])
                if ndb >0:
                    print("The following bonds are deleted.")
                    if ndb > 0:
                        for aB in aModRes.modLigand["deleted"]["bonds"]:
                            print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                naa = len(aModRes.modLigand["added"]["atoms"])
                if naa > 0 :
                    print("The following atoms are added.")
                    if naa > 0:
                        for aA in aModRes.modLigand["added"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])
                nab = len(aModRes.modLigand["added"]["bonds"])
                if nab > 0:
                    print("The following bonds are added.")
                    for aB in aModRes.modLigand["added"]["bonds"]:
                        print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                        print("Its bond order is ", aB["value_order"])
                nca = len(aModRes.modLigand["changed"]["atoms"])
                ncc = len(aModRes.modLigand["changed"]["charges"])
                
                if nca > 0 or ncc > 0 :
                    print("The following atoms are changed.")
                    if nca > 0:
                        for aA in aModRes.modLigand["changed"]["atoms"]:
                            print("Atom %s "%aA["atom_id"])
                    if ncc > 0:
                        for aA in aModRes.modLigand["changed"]["charges"]:
                            print("Atom %s with charge %d "%(aA["atom_id"], aA["charge"]))
                
        
    def processOneModRes(self, tModRes):
        
        if not self.errLevel:
            print("Atoms are read from ", tModRes.stdLigand["inCif"])
            self.setOneCompFromCifModRes(tModRes.stdLigand["inCif"], tModRes.stdLigand)
            
            print("Number of the original atoms ", len(tModRes.stdLigand["comp"]["atoms"]))
            for aAt in tModRes.stdLigand["comp"]["atoms"]:
                print("Atom ", aAt["atom_id"], " with charge ", aAt["charge"])
            if not self.errLevel:
                print("##################################################################")
                print("#                                                                #")
                print("#  Build a modification of the monomers                          #")    
                print("#                                                                #")
                print("##################################################################")
                self.buildOneModMonomer(tModRes)
                
                if not self.errLevel:
                    print("##################################################################")
                    print("#                                                                #")
                    print("#     Get information for the modifications                      #")    
                    print("#                                                                #")
                    print("##################################################################")
                    self.extractOneModMonomerInfo(tModRes)
                    # Finally print out 
                    if not self.errLevel:
                        print("##################################################################")
                        print("#                                                                #")
                        print("#          Output all information to the file                    #")    
                        print("#                                                                #")
                        print("##################################################################")
                        self.outOneModResInfo(tModRes)
                        print("#                          Job done                              #")
                        print("##################################################################")
    
    def buildOneModMonomer(self, tModRes):
        
        if not self.errLevel:
            self.adjustAtomsAndOthersForOneMod(tModRes)
            #for aAtom in tModRes.stdLigand["remainAtoms"]:
            #    print("XXXX atom ", aAtom["atom_id"])
            
            
            if not self.errLevel:
                tModRes.newLigand["name"] = tModRes.stdLigand["name"]
                tModRes.newLigand["inCif"] = os.path.join(self.scrDir, tModRes.newLigand["name"] + "_In.cif")
                self.newLigToSimplifiedMmcif(tModRes,  tModRes.newLigand["inCif"])
            
                lKeku = self.checkKekuInModRes(tModRes)
                if not self.errLevel:
                    if lKeku : 
                        print("re-keku needed")
                        self.setOneMonomer(tModRes.newLigand, 21)
                    else:
                        print("No need for re-keku")
                        self.setOneMonomer(tModRes.newLigand, 2)
                    if not self.errLevel:
                        #print("new ligand outcif ", tModRes.newLigand["outCif"])
                        tModRes.newLigand["cifObj"] = Ccp4MmCifObj(tModRes.newLigand["outCif"])["ccp4CifObj"]
                        
                        if os.path.isfile(tModRes.newLigand["outCif"]):
                            newCifName = self.outRoot + ".cif"
                            shutil.copy(tModRes.newLigand["outCif"], newCifName)
        
        
    def newLigToSimplifiedMmcif(self, tModRes, tOutFName):
        
        try: 
            aOutF = open(tOutFName, "w")
        except IOError:
            print("%s can not be open for writting "%tOutFName)
            self.errLevel    = 34 
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for writing "%tOutFName)
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
            
            tModRes.stdLigand["remainHAtoms"] = []
            for aAt in tModRes.stdLigand["remainAtoms"]:
                if aAt["type_symbol"] == "H":
                    tModRes.stdLigand["remainHAtoms"].append(aAt)
            
            aNewLigId =  tModRes.stdLigand["name"] 
            numAts = len(tModRes.stdLigand["remainAtoms"]) 
            numH   = len(tModRes.stdLigand["remainHAtoms"]) 
            numNonH = numAts - numH
            
            print("n atoms ", numAts)
            print("numH ", numH)
            print("numNonH ", numNonH)
            
            #aOutF.write("%s%s%s%s%6d%6d\n"%(tModRes.stdLigand["name"].ljust(8), 
            #                                tModRes.stdLigand["name"].ljust(8), 
            #                                ".".ljust(20),\
            #                               "non-polymer".ljust(20), numAts, numNonH))
            aName = tModRes.stdLigand["name"] + "m"
            aDisc = tModRes.stdLigand["list"]["name"]
            aGroup = tModRes.stdLigand["list"]["group"]
            aOutF.write("%s%s%s%s%6d%6d\n"%(aName.ljust(8), 
                                            aName.ljust(8), 
                                            aDisc.ljust(20),\
                                            aGroup.ljust(20), numAts, numNonH))
                    
            aOutF.write("data_comp_%s\n"%tModRes.stdLigand["name"])
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
                for aAtom in tModRes.stdLigand["remainAtoms"]:
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

                    aOutF.write("%s%s%s%s%s%s%s%s\n"%(aNewLigId.ljust(10), aAtom["atom_id"].ljust(8),\
                                aAtom["type_symbol"].ljust(6), aAtom["type_energy"].ljust(8),\
                                aCharge.ljust(10), x.ljust(10), y.ljust(10), z.ljust(10)))
                        
                        
            if not self.errLevel:
                if len(tModRes.stdLigand["remainBonds"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_bond.comp_id\n")
                    aOutF.write("_chem_comp_bond.atom_id_1\n")
                    aOutF.write("_chem_comp_bond.atom_id_2\n")
                    aOutF.write("_chem_comp_bond.type\n")
                    aOutF.write("_chem_comp_bond.value_dist\n")
                    aOutF.write("_chem_comp_bond.value_dist_esd\n")
                    for aBond in tModRes.stdLigand["remainBonds"]:
                        
                        aOutF.write("%s%s%s%s%s%s\n"%(aNewLigId.ljust(10), aBond["atom_id_1"].ljust(8),\
                                     aBond["atom_id_2"].ljust(8), aBond["type"].ljust(20), \
                                     str(aBond["value_dist"]).ljust(10), str(aBond["value_dist_esd"]).ljust(10)))
                else: 
                    self.errLevel    = 35
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    self.errMessage[self.errLevel].append("The combo-ligand has no bonds\n")
                
                if not self.errLevel and len(tModRes.stdLigand["remainChirs"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_chir.comp_id\n")
                    aOutF.write("_chem_comp_chir.id\n")
                    aOutF.write("_chem_comp_chir.atom_id_centre\n")
                    aOutF.write("_chem_comp_chir.atom_id_1\n")
                    aOutF.write("_chem_comp_chir.atom_id_2\n")
                    aOutF.write("_chem_comp_chir.atom_id_3\n")
                    aOutF.write("_chem_comp_chir.volume_sign\n")
                    i = 1
                    for aChi in tModRes.stdLigand["remainChirs"]:
                        aChiId = "chir_" + str(i)
                        aOutF.write("%s%s%s%s%s%s%s\n"%(aNewLigId.ljust(10), aChiId.ljust(10), aChi["atom_id_centre"].ljust(6),\
                                                        aChi["atom_id_1"].ljust(6), aChi["atom_id_2"].ljust(6),\
                                                        aChi["atom_id_3"].ljust(6), aChi["volume_sign"]))
                        i +=1
                aOutF.close()   
            
            
    def checkKekuInModRes(self, tModRes):
        
        lK = False
        
        for aAtom in tModRes.stdLigand["remainAtoms"]:
            aElem = aAtom["type_symbol"].upper()
            aId   = aAtom["atom_id"]
            connBs = self.getBondSetOnlyForOneAtomById(aId, tModRes.stdLigand["remainBonds"], tModRes.stdLigand["remainAtoms"])        
            if aElem in  self.chemCheck.orgVal and len(connBs) > 0:
                curVal = self.getTotalBondOrderInOneMmcifAtom(aId, tModRes.stdLigand["remainAtoms"], connBs)
                #print("atom ", aId, " has charge ", aAtom["charge"])
                aC = float(aAtom["charge"]) 
                if  aC !=0:
                    curVal -= aC 
                if not curVal in self.chemCheck.orgVal[aElem]:
                    lK = True
                    print("BondOrder inconsistent atom %s : actual val %s "
                          %(aId, curVal))
        return lK 
    
    def adjustAtomsAndOthersForOneMod(self, tModRes):
        
        #if len(tModRes.modLigand["changed"]["charges"]) > 0:
        self.addjustFormalChargeInOneModRes(tModRes.stdLigand, tModRes.modLigand)
        
        
        self.setAddedInOneResForModRes(tModRes.stdLigand, tModRes.modLigand)
        
        
        if not self.errLevel:
            self.setDeletedInOneResForModRes(tModRes.stdLigand, tModRes.modLigand)
         
        tempDelList =[]
        for aEnt in tModRes.modLigand["deleted"]["atoms"]:
            if "type_symbol" in aEnt.keys():
                tempDelList.append(aEnt)
        tModRes.modLigand["deleted"]["atoms"].clear()
        existIds = []
        for aAtm in tempDelList:
            if not aAtm["atom_id"] in existIds:
                tModRes.modLigand["deleted"]["atoms"].append(aAtm)
                existIds.append(aAtm["atom_id"])
        for aAtm in tModRes.modLigand["deleted"]["atoms"]:
            print("deleted ", aAtm["atom_id"])
        
    def extractOneModMonomerInfo(self, tModRes):
        
        
        #for aA in tModRes.modLigand["changed"]["atoms"]:
        #    print("here ", aA["atom_id"])
            
        addedSet = []
        existChangeAtmIdsRes = []
        if len(tModRes.modLigand["changed"]["atoms"]) > 0:
            for aA in tModRes.modLigand["changed"]["atoms"]:
                existChangeAtmIdsRes.append(aA["atom_id"])
        
        existAddAtmIdsRes = []
        tmpAddedAtoms = []
        if len(tModRes.modLigand["added"]["atoms"]) > 0:
            for aA in tModRes.modLigand["added"]["atoms"]:
                existAddAtmIdsRes.append(aA["atom_id"].strip())
                print("existing added ", aA["atom_id"] )
        
        compId = tModRes.newLigand["name"]
        #print("Check")
        for aAtom in tModRes.newLigand["cifObj"]["comps"][compId]["atoms"]:
            #print("Atom id ",aAtom["atom_id"])
            #print("ccp4-type ", aAtom["type_energy"])
            #print("charge ", aAtom["charge"])
            if aAtom["atom_id"].strip() in existAddAtmIdsRes: 
                #if not aAtom["is_added"]:
                #    self.checkAtomMod(tModRes.stdLigand["remainAtoms"], aAtom, tModRes.modLigand["changed"]["atoms"], existChangeAtmIdsRes)
                #    print("existChangeAtmIdsRes2 ", existChangeAtmIdsRes)
                #else:
                #print("It is in added section ")
                addedSet.append(aAtom["atom_id"])
                tmpAddedAtoms.append(aAtom)
                #tModRes.modLigand["added"]["atoms"].append(aAtom)
            else:
                #print("Its in changed section")
                #self.checkAtomMod(tModRes.stdLigand["remainAtoms"], aAtom, tModRes.modLigand["changed"]["atoms"], existChangeAtmIdsRes)
                self.checkAtomMod(tModRes.stdLigand["remainAtoms"], aAtom, tModRes.modLigand["changed"]["atoms"], existChangeAtmIdsRes)
                
        
        tModRes.modLigand["added"]["atoms"].clear()
        for aA in tmpAddedAtoms:
            tModRes.modLigand["added"]["atoms"].append(aA)
            
        print("Num of mod atoms in the monomer is %d "%len(tModRes.modLigand["changed"]["atoms"]))
        #if len(tModRes.modLigand["changed"]["atoms"]) >0:
        #    print("They are : ")
        #    for aA in tModRes.modLigand["changed"]["atoms"]:
        #        print(aA["atom_id"])
        print("Num of add atoms in the monomer is %d "%len(tModRes.modLigand["added"]["atoms"]))
        #if len(tModRes.modLigand["added"]["atoms"]) > 0:
        #    print("They are : ")
        #    for aA  in tModRes.modLigand["added"]["atoms"]:
        #        print(aA["atom_id"])
        print("Num of deleted atoms in the monomer is %d "%(len(tModRes.modLigand["deleted"]["atoms"])))
        #if len(tModRes.modLigand["deleted"]["atoms"]) > 0:
        #    print("They are : ")
        #    for aA  in tModRes.modLigand["deleted"]["atoms"]:
        #        print(aA["atom_id"])
    
        # Bonds
        aTmpChBonds =  {}
        i1 = 0
        for aB in tModRes.modLigand["changed"]["bonds"]:
            aList1 = [aB["atom_id_1"], aB["atom_id_2"]]
            aList1.sort()
            aStr1 = aList1[0] + "_" + aList1[1]
            aTmpChBonds[aStr1] = [aB]
            i1= i1 + 1
        
        if len(tModRes.modLigand["added"]["atoms"]) > 0:
            for aAtom in tModRes.modLigand["added"]["atoms"]:
                addedSet.append(aAtom["atom_id"])
                
        tModRes.modLigand["changed"]["bonds"] = []
        tModRes.modLigand["added"]["bonds"]   = []
        for aBond in tModRes.newLigand["cifObj"]["comps"][compId]["bonds"]:
            #print("A bond ")
            #print("atom %s in residue %d "%(aBond["atom_id_1"], aBond["atom_id_1_resNum"]))
            #print("atom %s in residue %d "%(aBond["atom_id_2"], aBond["atom_id_2_resNum"]))
            if not aBond["atom_id_1"] in addedSet and not aBond["atom_id_2"] in addedSet:
                    aList = [aBond["atom_id_1"], aBond["atom_id_2"]]
                    aList.sort()
                    aStr = aList[0] + "_" + aList[1]
                    if not aStr in list(aTmpChBonds.keys()):
                        self.checkBondMod(tModRes.stdLigand["remainBonds"], aBond, tModRes.modLigand["changed"]["bonds"])
                    else:
                        tModRes.modLigand["changed"]["bonds"].append(aBond) 
            else:
                tModRes.modLigand["added"]["bonds"].append(aBond)
        
        print("Number of changed bonds in modified residue  is %d "%len(tModRes.modLigand["changed"]["bonds"]))  
        print("Number of added bonds in modified residue is %d "%len(tModRes.modLigand["added"]["bonds"]))  
        print("Number of deleted bonds in modified residue  is %d "%len(tModRes.modLigand["deleted"]["bonds"]))
        
        # Angles
        for aAng in tModRes.newLigand["cifObj"]["comps"][compId]["angles"]:
            #print("A angle between atom %s, atom %s and atom %s "%(aAng["atom_id_1"], aAng["atom_id_2"], aAng["atom_id_3"]))
            if not aAng["atom_id_1"] in addedSet and not aAng["atom_id_2"] in addedSet and not aAng["atom_id_3"] in addedSet:
                self.checkAngMod(tModRes.stdLigand["remainAngs"], aAng, tModRes.modLigand["changed"]["angles"]) 
            else:
                tModRes.modLigand["added"]["angles"].append(aAng) 
        print("Number of changed angles is %d "%len(tModRes.modLigand["changed"]["angles"]))  
        print("Number of added angles is %d "%len(tModRes.modLigand["added"]["angles"]))  
        print("Number of deleted angles is %d "%len(tModRes.modLigand["deleted"]["angles"]))
        
        # Torsions
        torIdNumLig = {}
        torIdNumLig["SP2_SP2"]=0
        torIdNumLig["SP2_SP3"]=0
        torIdNumLig["SP3_SP3"]=0
        self.getTorIdSets(torIdNumLig, tModRes.stdLigand["remainTors"])
        print("torIdNumLig: ", torIdNumLig)
        
        if "tors" in tModRes.newLigand["cifObj"]["comps"][compId]:
            for aTor in tModRes.newLigand["cifObj"]["comps"][compId]["tors"]:
                #print("A torsion between : ")
                #print("atom %s,  %s, %s and %s "%(aTor["atom_id_1"], aTor["atom_id_2"], aTor["atom_id_3"],aTor["atom_id_4"]))   
                if not aTor["atom_id_1"] in addedSet and not aTor["atom_id_2"] in addedSet\
                   and not aTor["atom_id_3"] in addedSet and not aTor["atom_id_4"] in addedSet: 
                    self.checkTorMod(tModRes.stdLigand["remainTors"], aTor, tModRes.modLigand["changed"]["tors"], torIdNumLig)
                else: 
                    tModRes.modLigand["added"]["tors"].append(aTor)
        
        print("Number of changed tors is %d "%len(tModRes.modLigand["changed"]["tors"]))  
        print("Number of added tors is %d "%len(tModRes.modLigand["added"]["tors"]))  
        print("Number of deleted tors is %d "%len(tModRes.modLigand["deleted"]["tors"])) 
        
        # Chirs
        if "chirs" in tModRes.newLigand["cifObj"]["comps"][compId]:
            for aChi in tModRes.newLigand["cifObj"]["comps"][compId]["chirs"]:
                #print("A chiral center by the following atoms")
                #print("atom %s "%aChi["atom_id_centre"])  
                #print("atom %s "%aChi["atom_id_1"])   
                #print("atom %s "%aChi["atom_id_2"])   
                #print("atom %s "%aChi["atom_id_3"])  
                
                if not aChi["atom_id_centre"] in addedSet and not aChi["atom_id_1"] in addedSet\
                   and not aChi["atom_id_2"] in addedSet and not aChi["atom_id_3"] in addedSet:
                    self.checkChiMod(tModRes.stdLigand["remainChirs"], aChi, tModRes.modLigand["changed"]["chirs"])
                else:
                    tModRes.modLigand["added"]["chirs"].append(aChi)
        print("Number of changed chirs is %d "%len(tModRes.modLigand["changed"]["chirs"]))  
        print("Number of added chirs is %d "%len(tModRes.modLigand["added"]["chirs"]))  
        print("Number of deleted chirs is %d "%len(tModRes.modLigand["deleted"]["chirs"]))
        
        # Planes
        
        # checkPlModFromNewLigand
        print("Number of the orig planes  : ", len(tModRes.stdLigand["remainPls"]))
        aLDelPlIds = []
        if "planes" in tModRes.newLigand["cifObj"]["comps"][compId]:
            for aPl in list(tModRes.newLigand["cifObj"]["comps"][compId]["planes"].keys()):
                #print("For plane ", aPl)
                inModRes = []
                nAtmInPl = len(tModRes.newLigand["cifObj"]["comps"][compId]["planes"][aPl])
                for aAtom in tModRes.newLigand["cifObj"]["comps"][compId]["planes"][aPl]:
                    #print("atom %s in residue %d "%(aAtom["atom_id"], aAtom["atom_id_resNum"]))
                    inModRes.append(aAtom["atom_id"])
                    
                if len(inModRes) == nAtmInPl:
                    self.checkPlMod(tModRes.stdLigand["remainPls"], tModRes.newLigand["cifObj"]["comps"][compId]["planes"][aPl],\
                                    tModRes.modLigand["deleted"]["planes"], tModRes.modLigand["added"]["planes"], aLDelPlIds)
        # checkPlModFromOrig
        for aPl in tModRes.stdLigand["remainPls"]:
            if not "planes" in tModRes.newLigand["cifObj"]["comps"][compId]:
                tModRes.modLigand["deleted"]["planes"].append(aPl)  
            else:
                self.checkPlMod2(aPl, tModRes.newLigand["cifObj"]["comps"][compId]["planes"], 
                                 tModRes.modLigand["deleted"]["planes"])
        nDP = len(tModRes.modLigand["deleted"]["planes"])   
        nAP = len(tModRes.modLigand["added"]["planes"])
        print("Number of deleted planes ", nDP)
        print("Number of added planes ", nAP)
        
    def outOneModResInfo(self, tModRes):
        
        aOutCifName = self.outRoot + "_modres.cif"
        #print(aOutCifName) 
        try: 
            aOutCif = open(aOutCifName, "w")
        except IOError:
            self.errLevel    = 41
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for writting "%aOutCifName)
        else:
            self.outVerInfo(aOutCif)
            self.outCompListModRes(aOutCif, tModRes)
            self.outModListModRes(aOutCif, tModRes)   
            if tModRes.stdLigand["dataBlock"] and tModRes.stdLigand["outComp"]\
                and tModRes.stdLigand["compOut"]:
                self.outOneComp(aOutCif, tModRes.stdLigand)
            
            self.outOneModModRes(aOutCif, tModRes.modLigand)
            aOutCif.close()
    
    def outCompListModRes(self, tOutFile, tModRes):
        
        if tModRes.stdLigand["compOut"]: 
            tOutFile.write("data_comp_list\n\n")
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_comp.id\n")
            tOutFile.write("_chem_comp.three_letter_code\n")
            tOutFile.write("_chem_comp.name\n")
            tOutFile.write("_chem_comp.group\n")
            tOutFile.write("_chem_comp.number_atoms_all\n")
            tOutFile.write("_chem_comp.number_atoms_nh\n")
            
            aNL = tModRes.stdLigand["name"]
            aN3 =""
            if len(aNL) >= 3:
                aN3 = tModRes.stdLigand["name"][:3]
            else:
                aN3 = aNL
            aName   = tModRes.stdLigand["list"]["name"]
            aNameL  = len(aName) + 6
            aGrp    = tModRes.stdLigand["list"]["group"]
            NA      = tModRes.stdLigand["list"]["number_atoms_all"]
            NAH     = tModRes.stdLigand["list"]["number_atoms_nh"]
            aL="%s%s%s%s%s%s\n"%(aNL.ljust(10), aN3.ljust(10), aName.ljust(aNameL), aGrp.ljust(20), NA.ljust(10), NAH)
            tOutFile.write(aL) 
    
    def outModListModRes(self, tOutFile, tModRes):
        
        tOutFile.write("data_mod_list\n\n")
        tOutFile.write("loop_\n")
        tOutFile.write("_chem_mod.id\n")
        tOutFile.write("_chem_mod.name\n")
        tOutFile.write("_chem_mod.comp_id\n")
        tOutFile.write("_chem_mod.group_id\n")
      
        aMN     = tModRes.modLigand["name"]
        aName   = tModRes.stdLigand["list"]["name"]
        aNameL  = len(aName) + 6
        aLN     = tModRes.stdLigand["name"]
        aGrp    = tModRes.stdLigand["list"]["group"]
        aL="%s%s%s%s\n"%(aMN.ljust(10), aName.ljust(aNameL), aLN.ljust(10), aGrp.ljust(20))
        tOutFile.write(aL) 
        
    def outOneModModRes(self, tOutFile, tModLigand):
        
        tPoolAtoms = []
        
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
                    aC = "0.0"
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
                    aC = "0.0"
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aAtom["atom_id"].ljust(10),\
                                             ".".ljust(10), aAtom["type_symbol"].ljust(10), aAtom["type_energy"].ljust(10),\
                                             aC.ljust(10))
                    tOutFile.write(aL)

            if nAA !=0:
                for aAtom in tModLigand["added"]["atoms"]:
                    aC = "0.0"
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), ".".ljust(10),\
                                             aAtom["atom_id"].ljust(10), aAtom["type_symbol"].ljust(10),\
                                             aAtom["type_energy"].ljust(10), aC.ljust(10))
                    tOutFile.write(aL)
           
            tOutFile.write("\n")
        
        
            # Bonds

            nDB  = len(tModLigand["deleted"]["bonds"])
            nCB  = len(tModLigand["changed"]["bonds"])
            nAB  = len(tModLigand["added"]["bonds"])
            
            print("Number of deleted bonds ", nDB)
            print("number of changed bonds ", nCB)
            print("Number of added bonds ", nAB)
            
            if nDB !=0 or nCB !=0 or nAB !=0:
                tOutFile.write("loop_\n")
                tOutFile.write("_chem_mod_bond.mod_id\n")
                tOutFile.write("_chem_mod_bond.function\n")
                tOutFile.write("_chem_mod_bond.atom_id_1\n")
                tOutFile.write("_chem_mod_bond.atom_id_2\n")
                tOutFile.write("_chem_mod_bond.new_type\n")
                tOutFile.write("_chem_mod_bond.new_value_dist\n")
                tOutFile.write("_chem_mod_bond.new_value_dist_esd\n")               
                tOutFile.write("_chem_mod_bond.new_value_dist_nucleus\n")
                tOutFile.write("_chem_mod_bond.new_value_dist_nucleus_esd\n")               
         
                if nDB !=0:
                    for aBond in tModLigand["deleted"]["bonds"]:
                        aBT = aBond["type"].lower()
                        aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), ".".ljust(15),\
                                             ".".ljust(10), ".".ljust(10), ".".ljust(10))
                        tOutFile.write(aL)
                
                if nCB !=0:
                    #print tModLigand["name"]
                    for aBond in tModLigand["changed"]["bonds"]:
                        aBT = aBond["type"].lower()
                        #print aBond["atom_id_1"]
                        #print aBond["atom_id_2"]
                        #print aBond["value_dist"]
                        aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aBond["atom_id_1"].ljust(10),\
                                                     aBond["atom_id_2"].ljust(10), aBT.ljust(15), aBond["value_dist"].ljust(15),\
                                                     aBond["value_dist_esd"].ljust(10), aBond["value_dist_nucleus"].ljust(15),\
                                                     aBond["value_dist_nucleus_esd"].ljust(10))
                        tOutFile.write(aL)
                
                if nAB !=0:
                    for aBond in tModLigand["added"]["bonds"]:
                        #print(aBond.keys())
                        aBT = aBond["type"].lower()
                        if not "value_dist_nucleus" in aBond.keys():
                            print("No n-dist for bond between %s and %s "%(aBond["atom_id_1"].ljust(10), aBond["atom_id_2"].ljust(10)))
                            print("in ligand %s\n"%tModLigand["name"])
                            sys.exit()
                        else:
                            aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), aBond["atom_id_1"].ljust(10),\
                                                         aBond["atom_id_2"].ljust(10), aBT.ljust(15), str(aBond["value_dist"]).ljust(15),\
                                                         aBond["value_dist_esd"].ljust(10), aBond["value_dist_nucleus"].ljust(15),\
                                                         aBond["value_dist_nucleus_esd"].ljust(10))
                     
                            tOutFile.write(aL)
                
        
                tOutFile.write("\n")
                
                
                # Angles, depending on  the description level 

                nDAngs  = len(tModLigand["deleted"]["angles"])
                nCAngs  = len(tModLigand["changed"]["angles"])
                nAAngs  = len(tModLigand["added"]["angles"])
                
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
                        for aAng in tModLigand["changed"]["angles"]:
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
                nCTors  = len(tModLigand["changed"]["tors"])
                nATors  = len(tModLigand["added"]["tors"])
                
                if nDTors !=0 or nCTors != 0 or nATors !=0:
                    tOutFile.write("loop_\n")
                    tOutFile.write("_chem_mod_tor.mod_id\n")
                    tOutFile.write("_chem_mod_tor.function\n")
                    tOutFile.write("_chem_mod_tor.atom_id_1\n")
                    tOutFile.write("_chem_mod_tor.atom_id_2\n")
                    tOutFile.write("_chem_mod_tor.atom_id_3\n")
                    tOutFile.write("_chem_mod_tor.atom_id_4\n")
                    tOutFile.write("_chem_mod_tor.id\n")
                    tOutFile.write("_chem_mod_tor.new_value_angle\n")
                    tOutFile.write("_chem_mod_tor.new_value_angle_esd\n")
                    tOutFile.write("_chem_mod_tor.new_period\n")
                
                
                    if nDTors !=0: 
                        for aTor in tModLigand["deleted"]["tors"]:
                            aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), \
                                                          aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                                          aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                                          ".".ljust(12), ".".ljust(12), ".".ljust(12), str(aTor["period"]).ljust(15))                  
                            tOutFile.write(aL)
           
                    if nCTors !=0: 
                        for aTor in tModLigand["changed"]["tors"]:
                            aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), \
                                                          aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                                               aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                                          aTor["id"].ljust(18),\
                                                          aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15),\
                                                          aTor["period"].ljust(6))                  
                            tOutFile.write(aL)
   
                    if nATors !=0: 
                        for aTor in tModLigand["added"]["tors"]:
                            aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), \
                                                      aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                                      aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                                      aTor["id"].ljust(18),\
                                                      aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15),
                                                      aTor["period"].ljust(6))                  
                            tOutFile.write(aL)
   
                    tOutFile.write("\n")
                    
                    
                # Chiral centers 

                nDChirs  = len(tModLigand["deleted"]["chirs"])
                nCChirs  = len(tModLigand["changed"]["chirs"])
                nAChirs  = len(tModLigand["added"]["chirs"])
                    
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
                        for aChi in tModLigand["changed"]["chirs"]:
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
                print("number of deleted plane ", nDPls)
                print("number of changed plane ", nCPls)
                print("number of added plane ", nAPls)
                if nDPls or nCPls or nAPls:
                    tOutFile.write("loop_\n")
                    tOutFile.write("_chem_mod_plane_atom.mod_id\n")
                    tOutFile.write("_chem_mod_plane_atom.function\n")
                    tOutFile.write("_chem_mod_plane_atom.plane_id\n")
                    tOutFile.write("_chem_mod_plane_atom.atom_id\n")
                    tOutFile.write("_chem_mod_plane_atom.new_dist_esd\n")
                if nDPls: 
                    for aPl in tModLigand["deleted"]["planes"]:
                        for aPlAtm in aPl:
                            aL ="%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15),\
                                                aPlAtm["plane_id"].ljust(15), aPlAtm["atom_id"].ljust(15), aPlAtm["dist_esd"]) 
                            tOutFile.write(aL)
                if nAPls:
                    for aPl in tModLigand["added"]["planes"]:
                        for aPlAtm in aPl:
                            aL ="%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15),\
                                                aPlAtm["plane_id"].ljust(15), aPlAtm["atom_id"].ljust(15), aPlAtm["dist_esd"]) 
                            tOutFile.write(aL)

                tOutFile.write("\n")

            
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
                self.setOneMonomer(tLinkIns.stdLigand1,3)
                #print "output comp 1 ", tLinkIns.stdLigand1["outComp"] 
                
            if not self.errLevel:
                if not tLinkIns.stdLigand2["fromScr"]:
                    print("Here inCif ", tLinkIns.stdLigand2["inCif"])
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
                        
                        #print("Number of H atoms in residue %s is %d "%(tLinkIns.stdLigand2["name"], len(tLinkIns.stdLigand2["comp"]["hAtoms"])))
                else:
                    self.setOneMonomer(tLinkIns.stdLigand2,3)
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
        
        if not aMmcifObj["errLevel"] :
            #print(list(aMmcifObj["ccp4CifObj"].keys()))
            #print(list(aMmcifObj["ccp4CifObj"]["comps"].keys()))
            if tMonomer["name"] in aMmcifObj["ccp4CifObj"]["comps"]:
                tMonomer["outComp"] = True
                tMonomer["comp"] = aMmcifObj["ccp4CifObj"]["comps"][tMonomer["name"]]
                if (not tMonomer["userIn"]) and (tMonomer["name"].upper() in self.chemCheck.aminoAcids):
                    self.chemCheck.tmpModiN_in_AA(tMonomer["name"].upper(), tMonomer["comp"])
                self.selectHAtoms(tMonomer["comp"])   
                #print(aMmcifObj["ccp4CifObj"]["lists"]["comp"].keys())
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
                for iB in range(len(tMonomer["comp"]["bonds"])):
                    if "type" in tMonomer["comp"]["bonds"][iB] and not "value_order" in tMonomer["comp"]["bonds"][iB]:
                        tMonomer["comp"]["bonds"][iB]["value_order"] = tMonomer["comp"]["bonds"][iB]["type"]
                    elif not "type" in tMonomer["comp"]["bonds"][iB] and "value_order" in tMonomer["comp"]["bonds"][iB]:
                        tMonomer["comp"]["bonds"][iB]["type"] = tMonomer["comp"]["bonds"][iB]["value_order"] 
                        
                self.checkDelocAndAromaBonds(tMonomer["comp"], tMonomer["name"])

        else:
            self.errLevel = aMmcifObj["errLevel"]
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append(aMmcifObj["errMessage"])

    def setOneCompFromCifModRes(self, tFileName, tMonomer):
     
        # Using the input file or the file in ccp4 monomer lib as it is
        aMmcifObj = Ccp4MmCifObj(tFileName)
        aMmcifObj.checkBlockCompsExist()
        
        
        if not aMmcifObj["errLevel"]:
            #print(list(aMmcifObj["ccp4CifObj"].keys()))
            #print(list(aMmcifObj["ccp4CifObj"]["comps"].keys()))
            if tMonomer["name"] in aMmcifObj["ccp4CifObj"]["comps"]:
                tMonomer["outComp"] = True
                tMonomer["comp"] = aMmcifObj["ccp4CifObj"]["comps"][tMonomer["name"]]
                if (not tMonomer["userIn"]) and (tMonomer["name"].upper() in self.chemCheck.aminoAcids):
                    self.chemCheck.tmpModiN_in_AA(tMonomer["name"].upper(), tMonomer["comp"])
                self.selectHAtoms(tMonomer["comp"])   
                #print(aMmcifObj["ccp4CifObj"]["lists"]["comp"].keys())
                tMonomer["list"] = aMmcifObj["ccp4CifObj"]["lists"]["comp"][tMonomer["name"]]
                #print(tMonomer["list"])
            
                
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
                for iB in range(len(tMonomer["comp"]["bonds"])):
                    if "type" in tMonomer["comp"]["bonds"][iB] and not "value_order" in tMonomer["comp"]["bonds"][iB]:
                        tMonomer["comp"]["bonds"][iB]["value_order"] = tMonomer["comp"]["bonds"][iB]["type"]
                    elif not "type" in tMonomer["comp"]["bonds"][iB] and "value_order" in tMonomer["comp"]["bonds"][iB]:
                        tMonomer["comp"]["bonds"][iB]["type"] = tMonomer["comp"]["bonds"][iB]["value_order"] 
                        
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
        
    def setOneMonomer(self, tMonomer, tMode=1):             # tMode ==1 link generation
                                                            # tMode ==2 generate a modification of a  monomer
        #print("tMode==", tMode)
        if os.path.isfile(tMonomer["inCif"]):
            #print(tMonomer["inCif"])
            #aNL = tMonomer["name"].upper()
            #if len(aNL) > 3:
            if tMode ==1:
                aNL = "LIG"
                self._log_name  = os.path.join(self.scrDir, aNL + "_for_link.log")
                self.subRoot    = os.path.join(self.scrDir, aNL + "_for_link")
                self._cmdline   = "%s -c %s  -r %s -o %s "%(self.exeAcedrg, tMonomer["inCif"], aNL, self.subRoot)   
                if "modifiedPlanes" in self.setParas.keys():
                    self._cmdline   += " -M "
            elif tMode ==2:
                aNL = tMonomer["name"]
                self._log_name  = os.path.join(self.scrDir, aNL + "_for_Mod.log")
                #self.subRoot    = os.path.join(self.scrDir, aNL + "_for_Mod")
                self.subRoot     = aNL + "_for_Mod"
                self._cmdline   = "%s -c %s  -r %s -o %s -K "%(self.exeAcedrg, tMonomer["inCif"], aNL, self.subRoot)
            elif tMode ==21:
                aNL = tMonomer["name"]
                self._log_name  = os.path.join(self.scrDir, aNL + "_for_Mod.log")
                #self.subRoot    = os.path.join(self.scrDir, aNL + "_for_Mod")
                self.subRoot  = aNL + "_for_Mod"
                self._cmdline   = "%s -c %s  -r %s -o %s --keku -K "%(self.exeAcedrg, tMonomer["inCif"], aNL, self.subRoot)
            elif tMode ==3:
                aNL = tMonomer["name"]
                self._log_name  = os.path.join(self.scrDir, aNL + ".log")
                self.subRoot    = os.path.join(self.scrDir, aNL)
                self._cmdline   = "%s -c %s  -r %s -o %s "%(self.exeAcedrg, tMonomer["inCif"], aNL, self.subRoot)
            elif tMode ==4:
                aNL = tMonomer["name"]
                self._log_name  = os.path.join(self.scrDir, aNL + ".log")
                self._cmdline   = "%s -c %s  -r %s -o %s "%(self.exeAcedrg, tMonomer["inCif"], aNL, self.subRoot)
            if tMode ==1 :
                print ("Link generation")
            elif tMode==2 or tMode==21:
                print("modification generation")
            print(self._cmdline)
            self.runExitCode = self.subExecute()
            print("tMode=", tMode)
            print("self.subRoot ", self.subRoot)
            #self.runExitCode = os.system(self._cmdline)
            if not self.runExitCode :
                if tMode !=4:
                     aOutLigCif1 = self.subRoot + ".cif"
                     aOutLigCif2 = self.subRoot + "_final.cif"
                     if os.path.isfile(aOutLigCif1):
                         aOutLigCif = aOutLigCif1
                     elif os.path.isfile(aOutLigCif2):
                         aOutLigCif = aOutLigCif2
                     print("1. intermediate cif is ",aOutLigCif)
                else:
                    aOutLigCif = self.subRoot + "_tmp.cif"
                    print("2. intermediate cif is ",aOutLigCif)
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
    
    def checkMetalAndAetOneMonomer(self, tMon):                                 # tMon = .stdLigand1 or .stdLigand2 or .stdLigand
    
        # Look into the file to see if it contains metal element and then dealt with them
        aMmcifObj = Ccp4MmCifObj(tMon["inCif"])
        aMmcifObj.checkBlockCompsExist()
        if not aMmcifObj["errLevel"]:
            if tMon["name"] in aMmcifObj["ccp4CifObj"]["comps"]:
                self.lNoMetal=self.chemCheck.isOrganicInCif(aMmcifObj["ccp4CifObj"]["comps"][tMon["name"]]["atoms"]) 
                if not self.lNoMetal:    # contain metal 
                    print(tMon["name"], " contains metal atoms ")
                    self.subRoot = tMon["name"] + "_intmedia"
                    self.setOneMonomer(tMon, 4)
                    #aIntmedFN = self.subRoot + "_tmp.cif"
                    aIntmedFN = self.subRoot + "_tmp.cif"
                    if os.path.isfile(aIntmedFN):
                        tMon["inCif"]=aIntmedFN  
                        print("one component cif is now",  tMon["inCif"])
                        
        
        
    def adjustAtomsAndOthersForComboLigand(self, tLinkedObj):
         
        self.setRemainBonds(tLinkedObj.stdLigand1, tLinkedObj.modLigand1)
        self.setRemainBonds(tLinkedObj.stdLigand2, tLinkedObj.modLigand2)
        
        if len(tLinkedObj.modLigand1["changed"]["formal_charges"]) > 0:
            self.addjustFormalChargeInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1)
        elif len(tLinkedObj.modLigand2["changed"]["formal_charges"]) > 0:
            self.addjustFormalChargeInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2)
            
        self.setAddedInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
        
        self.setAddedInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)
        
        
        
        if not self.errLevel:
            self.setDeletedInOneResForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
            if len(tLinkedObj.stdLigand1["linkChir"]) ==2:
                self.setLinkChir(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.stdLigand2)
            if not self.errLevel:
                self.setDeletedInOneResForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)
                
                if len(tLinkedObj.stdLigand2["linkChir"]) ==2:
                    self.setLinkChir(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.stdLigand1)
                #if not self.errLevel:
                    #self.setChargeInLinkAtom(tLinkedObj.stdLigand1, tLinkedObj.modLigand1, tLinkedObj.suggestBonds)
                    #if not self.errLevel:
                        #self.setChargeInLinkAtom(tLinkedObj.stdLigand2, tLinkedObj.modLigand2, tLinkedObj.suggestBonds)
        
        if self.errLevel:
            print("Error level is ", self.errLevel)
            print("Error message is : ")
            for aL in self.errMessage[self.errLevel]:
                print(aL)
                
    def setRemainBonds(self, tRes, tMod):
        # changed band 
        tmpAllBs = []
        for aB in tRes["comp"]["bonds"]:
            tmpAllBs.append(aB)
            
        i = 0
        chBondIdMap = {}
        print(" changed bonds  ", len(tMod["changed"]["bonds"]))
        for chBond in tMod["changed"]["bonds"]:
            aList = [chBond["atom_id_1"], chBond["atom_id_2"]]
            aList.sort()
            aStr = aList[0] + "_" + aList[1]
            chBondIdMap[aStr] = i  
            i = i+1
        
        for aBond in tmpAllBs:
            bList = [aBond["atom_id_1"], aBond["atom_id_2"]]
            bList.sort()
            bStr = bList[0] + "_" + bList[1]
            if bStr in list(chBondIdMap.keys()):
                tRes["remainBonds"].append(tMod["changed"]["bonds"][chBondIdMap[bStr]])
            else:
                tRes["remainBonds"].append(aBond)
            

    def  setLinkChir(self, tRes, tMod, tOthRes):
        
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
    
    def addjustFormalChargeInOneModRes(self, tRes, tMod):
   
        
        
        changeAtms = []
        if len(tMod["changed"]["atoms"]) > 0:
            for aA in tMod["changed"]["atoms"]:
                changeAtms.append(aA["atom_id"].strip().upper())
        
        addedAtms = []
        if len(tMod["added"]["atoms"]) > 0:
            for aA in tMod["added"]["atoms"]:
                addedAtms.append(aA["atom_id"].strip().upper())
                
        for aBond in tMod["added"]["bonds"]:
            atm1 = aBond["atom_id_1"]
            atm2 = aBond["atom_id_2"]
            print("Bond between ", atm1,  " and ", atm2, " is Aadded ")
            if atm1 in addedAtms and not atm2 in addedAtms:
                idxAtm2 = self.getAtomById(atm2, tRes["comp"]["atoms"])
                if not idxAtm2==-1 and int(tRes["comp"]["atoms"][idxAtm2]["charge"]) !=0:
                    tRes["comp"]["atoms"][idxAtm2]["charge"]=0
                    tMod["changed"]["atoms"].append(tRes["comp"]["atoms"][idxAtm2])
                    print("Atom ",  tRes["comp"]["atoms"][idxAtm2]["atom_id"], " charge : ",  tRes["comp"]["atoms"][idxAtm2]["charge"])
            elif atm2 in addedAtms and not atm1 in addedAtms:
                idxAtm1 = self.getAtomById(atm1, tRes["comp"]["atoms"])
                if not idxAtm1==-1 and int(tRes["comp"]["atoms"][idxAtm1]["charge"]) !=0:
                    tRes["comp"]["atoms"][idxAtm1]["charge"]=0
                    tMod["changed"]["atoms"].append(tRes["comp"]["atoms"][idxAtm1])         
                    print("Atom ",  tRes["comp"]["atoms"][idxAtm1]["atom_id"], " charge : ",  tRes["comp"]["atoms"][idxAtm1]["charge"])
        
        changeId = []
        for aTmpAtm in tMod["changed"]["atoms"]:
            changeId.append(aTmpAtm["atom_id"])
         
        for aFC in tMod["changed"]["charges"]:
            aId = aFC["atom_id"].strip().upper()
            for aAt in tRes["comp"]["atoms"]:
                if aAt["atom_id"].strip().upper() ==aId:
                    aAt["charge"] = aFC["charge"]
                    if not aAt["atom_id"] in  changeId:
                        tMod["changed"]["atoms"].append(aAt)
                        #if not aId in changeAtms:
                        #    addedHs=self.checkAssocHAtoms(tRes, tMod, aAt)
                        break
    
        
        #for aAt in tRes["comp"]["atoms"]:
        #    print("aId ", aAt["atom_id"])
        #    print("Its charge ", aAt["charge"])
        if len(tMod["changed"]["atoms"]) > 0:
            print("2. The following atoms are changed: ")
            for aAt in tMod["changed"]["atoms"]:
                print("Atom ",aAt["atom_id"]) 
                print("Formal charge ", aAt["charge"])
        
         
    def addjustFormalChargeInOneResForModification(self, tRes, tMod):
        
        
   
    
        changeAtms = []
        if len(tMod["changed"]["atoms"]) > 0:
            for aA in tMod["changed"]["atoms"]:
                changeAtms.append(aA["atom_id"].strip().upper())
        for aAt in tRes["remainAtoms"]:
            print("Remain id ", aAt["atom_id"], " charge ", aAt["charge"]), " "
            if aAt["atom_id"] == tRes["atomName"]:
                tRes["origCharge"][aAt["atom_id"]] = aAt["charge"]
                aAt["charge"] = "0"
        
        for aFC in tMod["changed"]["formal_charges"]:
            aId = aFC["atom_id"].strip().upper()
            for aAt in tRes["comp"]["atoms"]:
                if aAt["atom_id"].strip().upper() ==aId:
                    
                    tRes["origCharge"][aAt["atom_id"]] = aAt["charge"]
                    aAt["charge"] = aFC["formal_charge"]
                    
                    if not aId in changeAtms and aAt["atom_id"] != tRes["atomName"]:
                        addedHs=self.checkAssocHAtoms(tRes, tMod, aAt)
                    break
        
       
        if len(tMod["changed"]["atoms"]) > 0:
            print("The following atoms are changed: ")
            for aAt in tMod["changed"]["atoms"]:
                print("Atom ",aAt["atom_id"]) 
                print("Charge ", aAt["charge"])
    
    def checkAssocHAtoms(self, tRes, tMod, tAt):
        
        tDelAtomIds = []
        allHIds     = []
        iH          = 0
        for idxA in range(len(tRes["remainAtoms"])): 
            if tRes["remainAtoms"][idxA]["type_symbol"] == "H":        
                allHIds.append(tRes["remainAtoms"][idxA]["atom_id"])
        
        #[bondAtomSet, aLABonds] =  self.getBondSetForOneAtomById(tAt["atom_id"], tRes["comp"]["atoms"], tRes["comp"]["bonds"], tDelAtomIds)
        [bondAtomSet, aLABonds] =  self.getBondSetForOneAtomById(tAt["atom_id"], tRes["comp"]["atoms"], tRes["remainBonds"], tDelAtomIds)
        nTotalVa=self.getTotalBondOrderInOneMmcifAtom(tAt["atom_id"], tRes["remainAtoms"], aLABonds)  
        aCharge = 0
    
        if "charge" in tAt:
            aCharge = int(tAt["charge"])
           
        if aCharge != 0: 
            nTotalVa = nTotalVa - aCharge
        if tAt["type_symbol"] in self.chemCheck.orgVal:
            if not nTotalVa in self.chemCheck.orgVal[tAt["type_symbol"]]:
                nVDiff = nTotalVa - self.chemCheck.orgVal[tAt["type_symbol"]][0]
                nAbs   = abs(nVDiff)
                if nVDiff < 0: 
                    while iH < nAbs:
                        self.addOneHInRes(tAt, bondAtomSet, allHIds, tRes, tMod) 
                        iH+=1
        
        return iH            
                    
    def checkAssocHAtoms2(self, tRes, tMod, tAt, tDelBAtmIds):
        
        allHIds     = []
        iH          = 0
        for idxA in range(len(tRes["comp"]["atoms"])): 
            if tRes["comp"]["atoms"][idxA]["type_symbol"] == "H":        
                allHIds.append(tRes["comp"]["atoms"][idxA]["atom_id"])
                
        [bondAtomSet, aLABonds] =  self.getBondSetForOneAtomById2(tAt["atom_id"], 
                                   tRes["comp"]["atoms"], tRes["comp"]["bonds"], tDelBAtmIds)
        nTotalVa=self.getTotalBondOrderInOneMmcifAtom(tAt["atom_id"], tRes["comp"]["atoms"], aLABonds)  
        aCharge = 0
        if "charge" in tAt:
            aCharge = int(tAt["charge"])
        if aCharge != 0:
            nTotalVa = nTotalVa - aCharge
        print("atom ",tAt["atom_id"] , "charge ", aCharge, " equiv bond-order ", nTotalVa)
        
        if tAt["type_symbol"] in self.chemCheck.orgVal:
            if not nTotalVa in self.chemCheck.orgVal[tAt["type_symbol"]]:
                nVDiff = nTotalVa - self.chemCheck.orgVal[tAt["type_symbol"]][0]
                nAbs   = abs(nVDiff)
                if nVDiff < 0:
                    
                    while iH < nAbs:
                        self.addOneHInRes2(tAt, bondAtomSet, allHIds, tRes, tMod) 
                        iH+=1
        return iH        
    
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
               #else:
               #    if "charge" in tRes["remainAtoms"][idxA]:
               #        tRes["remainAtoms"][idxA]["charge"] = "0"
  
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

        #print("For residue ", tRes["name"])
        
        if len(tMod["added"]["atoms"]) > 0:
            for aAtom in tMod["added"]["atoms"]:
                #print("Add atom ", aAtom["atom_id"])
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
    
    def setAddedInOneResForModRes(self, tRes, tMod):

        print("For residue ", tRes["name"])
    
        aSetDelIds = []
        if len(tMod["added"]["atoms"]) > 0:
            initAddedAtoms = []
            for aAtom in tMod["added"]["atoms"]:
                initAddedAtoms.append(aAtom)
                print(aAtom["atom_id"], " is added ")
            
            initAddedBonds = []
            for aBond in tMod["added"]["bonds"]:
                initAddedBonds.append(aBond)
            
            initDelAtoms = []
            
            for aAtom in tMod["deleted"]["atoms"]:
                aSetDelIds.append(aAtom["atom_id"])
            #print("Initially deleted atoms are: ", aSetDelIds)
            
            
            
            for aAtom in initAddedAtoms:
                print("Add atom ", aAtom["atom_id"])
                aBool = self.checkDubAtomNameInOneRes(tRes, aAtom)
                if not aBool : 
                    self.addOneAtomAndConnectedBondsToModRes(tRes, tMod, aAtom,aSetDelIds, initAddedAtoms, initAddedBonds, aSetDelIds)
            
                else:
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    aLine = "Atom Name Dublication : added atom %s already exists in residue %s\n"\
                            %(aAtom["atom_id"], tRes["name"])
                    #print(aLine) 
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
    
    def addOneAtomAndConnectedBondsToModRes(self, tRes, tMod, tAtom, tDelIds, tInitAddedAtoms, tInitAddedBonds, tInitDelAtomIds):

        tRes["remainAtoms"].append(tAtom)
        
        tmpBandA = []
        tmpAtoms = []
        tmpBonds = []
        
        addAtmMap = {}
        i=0
        for aAtm in tInitAddedAtoms:
            addAtmMap[aAtm["atom_id"]] = i
            i+=1

        nAtoms = len(tRes["comp"]["atoms"])
        #for aBond in tMod["added"]["bonds"]:
        for aBond in tInitAddedBonds:
            tRes["remainBonds"].append(aBond)
            print("check bond between ",  aBond["atom_id_1"], " and ",  aBond["atom_id_2"])
            if aBond["atom_id_1"]==tAtom["atom_id"]:
                atmIdx = self.getOneAtomSerialById(aBond["atom_id_2"], tRes["comp"]["atoms"])
                if atmIdx > 0 and atmIdx < nAtoms:
                    if not tRes["comp"]["atoms"][atmIdx]["atom_id"] in tInitDelAtomIds:
                        tmpBandA.append([aBond, tRes["comp"]["atoms"][atmIdx]])
                        tmpAtoms.append(tRes["comp"]["atoms"][atmIdx])
                        tmpBonds.append(aBond)
                elif aBond["atom_id_2"] in addAtmMap:
                    idxA=addAtmMap[aBond["atom_id_2"]]
                    tmpBandA.append([aBond, tInitAddedAtoms[idxA]])
                    tmpAtoms.append(tInitAddedAtoms[idxA])
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
                    if not tRes["comp"]["atoms"][atmIdx]["atom_id"] in tInitDelAtomIds:
                        tmpBandA.append([aBond, tRes["comp"]["atoms"][atmIdx]])
                        tmpAtoms.append(tRes["comp"]["atoms"][atmIdx])
                        tmpBonds.append(aBond)
                elif aBond["atom_id_1"] in addAtmMap:
                    idxA=addAtmMap[aBond["atom_id_1"]]
                    tmpBandA.append([aBond, tInitAddedAtoms[idxA]])
                    tmpAtoms.append(tInitAddedAtoms[idxA])
                    tmpBonds.append(aBond)
                else:
                    self.errLevel = 12
                    if self.errLevel not in self.errMessage:
                        self.errMessage[self.errLevel] = []
                    aLine = "Can not find idx for atom ", aBond["atom_id_1"] 
                    self.errMessage[self.errLevel].append(aLine)
                    break
       
        
        if len(tmpBonds) > 0:
            # Check bonds around the added atoms
            print("Check added atom %s and around bonds "%tAtom["atom_id"])
            #aPair = self.chemCheck.valideBondOrderForOneOrgAtom2(tAtom, tmpBonds)
            self.chemCheck.adjustNBForOneAddedAtom(tAtom, tmpAtoms, tmpBonds, tRes, tMod, tDelIds)
        
            """           
                self.errLevel = 12
                if self.errLevel not in self.errMessage:
                    self.errMessage[self.errLevel] = []
                self.errMessage[self.errLevel].append(aPair[1])
                print(aPair[1])
            else:
                print("Bonds connected to atom %s are OK "%tAtom["atom_id"])
            """
        else :
            aLine = "No bonds are defined to connected to atom %s, check your input file! \n"%tAtom["atom_id"]
            self.errLevel = 12
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append(aLine)
    
        #print("tmpDelIds ", tDelIds)
        if not self.errLevel:
            #if len(tmpBandA) >0:
                #for aPair in tmpBandA:
                #    print("Added bond between atom %s and %s "%(aPair[0]["atom_id_1"], aPair[0]["atom_id_2"]))
                #    print("Need to adjust bonds connected to atom %s "%aPair[1]["atom_id"])
            
            
            for [aBond, aAtom] in tmpBandA:
                
                aBandASet = self.getBondSetForOneLinkedAtomModRes(aAtom["atom_id"], tRes["comp"]["atoms"],
                                                            tRes["comp"]["bonds"], tDelIds, tInitAddedAtoms, tInitAddedBonds)
                #aBandASet[0].append(tAtom)
                #aBandASet[1].append(aBond)
                if len(aBandASet[0]) > 0:
                    #print("Check the bonds connected atom %s "%aAtom["atom_id"])
                    #print("These bonds are : ")
                    #for aB in aBandASet[1]:
                    #    print("Bond between atom %s and %s "%(aB["atom_id_1"], aB["atom_id_2"]))
                    #print("Atoms involved are: ")
                    #for aA in aBandASet[0]:
                    #    print("Atom %s "%aA["atom_id"])
                    [noErr, errLine] =self.chemCheck.adjustNBForOneAddedAtom(aAtom, aBandASet[0], aBandASet[1], tRes, tMod, tDelIds)
                    if not noErr:
                        self.errLevel = 12
                        if self.errLevel not in self.errMessage:
                            self.errMessage[self.errLevel] = []
                        self.errMessage[self.errLevel].append(errLine)  
                        break
        
        
        
        
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
        
        # Delecte bonds as instructed directly by the input file
    
    
        self.setDeletedBondInOneResForModification(tRes, tMod)
    
    
        # Delecte atoms 
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
        #print("Those atoms are :")
        #for aAtom in  tRes["remainAtoms"]:
        #    if "charge" in aAtom.keys():
        #        print("Atom %s with charge %s "%(aAtom["atom_id"], aAtom["charge"]))
   
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
        #for aBond in tRes["comp"]["bonds"]:
        for aBond in tmpRemBs:
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
                    
       
        print("Number of remained bonds : ", len(tRes["remainBonds"]))
        #print("They are : ")
        #for aBond in tRes["remainBonds"]:
        #    print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))
        #    print("Bond-order is %s "%aBond["type"])     
        print("Number of deleted bonds : ", len(tMod["deleted"]["bonds"]))
        #if len(tMod["deleted"]["bonds"]):
        #    print("They are : ")
        #    for aBond in tMod["deleted"]["bonds"]: 
        #        print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))      
        #        print("Bond-order is %s "%aBond["type"])  


        
 
        # Delete all angles that contain the deleted atom
        for aAng in tRes["comp"]["angles"]:
            if (aAng["atom_id_1"].upper() in delAtomIdSet)\
               or (aAng["atom_id_2"].upper() in delAtomIdSet)\
               or (aAng["atom_id_3"].upper() in delAtomIdSet):
                tMod["deleted"]["angles"].append(aAng)
            else:
                tRes["remainAngs"].append(aAng)
        print("Number of remained Angles : ", len(tRes["remainAngs"]))
        #print("They are : ")
        #for aAng in tRes["remainAngs"]:
        #    print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
        print("Number of deleted angles ", len(tMod["deleted"]["angles"]))
        #print("They are : ")
        #for aAng in tMod["deleted"]["angles"]:
        #    print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
                  
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
         
        self.setDeletedAngInOneResForModification(tRes, tMod)   
        self.setDeleteTorsInOneResForModification(tRes, tMod) 
        self.setDeleteChirsInOneResForModification(tRes, tMod)

    def setDeletedInOneResForModRes(self, tRes, tMod):
        
        # Delecte bonds as instructed directly by the input file
        
        self.setDeletedBondInOneResForModification(tRes, tMod)
        
        # Delecte atoms 
        print("For residue ", tRes["name"])
        delAtomIdSet =[]
        for aAtom in tMod["deleted"]["atoms"]:
            delAtomIdSet.append(aAtom["atom_id"])

        if len(delAtomIdSet) > 0:
            tMod["deleted"]["atoms"] = []
            for aId in delAtomIdSet:
                self.deleteOneAtomAndConnectedHAtoms(tRes, tMod, aId)            
        print("Number of deleted atoms according to the instruction file is ", len(delAtomIdSet))
        #if len(delAtomIdSet) > 0:
        #    print("They are : ")
        #    for aId in delAtomIdSet: 
        #        print("atom ", aId)
    
        
        extraAtomSet = []
        if len(tMod["changed"]["bonds"]) > 0:
            for aB in tMod["changed"]["bonds"]:
                if aB["atom_id_1"] == tRes["atomName"]:
                    extraAtomSet.append(aB["atom_id_2"])
                elif aB["atom_id_2"] == tRes["atomName"]:
                    extraAtomSet.append(aB["atom_id_1"])
        tLinkBonds = None  
        for atmId in extraAtomSet:
            self.adjustAtomsAroundOneAtom(atmId, tRes, tMod, tLinkBonds, delAtomIdSet, 2)
      
        delAtomIdSet =[]
        for aAtom in tMod["deleted"]["atoms"]:
            #print(list(aAtom.keys()))
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
        #for aBond in tRes["comp"]["bonds"]:
        for aBond in tmpRemBs:
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
                    
       
        print("Number of remained bonds : ", len(tRes["remainBonds"]))
        #print("They are : ")
        #for aBond in tRes["remainBonds"]:
        #    print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))
        #    print("Bond-order is %s "%aBond["type"])     
        print("Number of deleted bonds : ", len(tMod["deleted"]["bonds"]))
        #if len(tMod["deleted"]["bonds"]):
        #    print("They are : ")
        #     for aBond in tMod["deleted"]["bonds"]: 
        #        print("Bond between atom %s and %s "%(aBond["atom_id_1"], aBond["atom_id_2"]))      
        #        print("Bond-order is %s "%aBond["type"])  
        
 
        # Delete all angles that contain the deleted atom
        for aAng in tRes["comp"]["angles"]:
            if (aAng["atom_id_1"].upper() in delAtomIdSet)\
               or (aAng["atom_id_2"].upper() in delAtomIdSet)\
               or (aAng["atom_id_3"].upper() in delAtomIdSet):
                tMod["deleted"]["angles"].append(aAng)
            else:
                tRes["remainAngs"].append(aAng)
        print("Number of remained Angles : ", len(tRes["remainAngs"]))
        #print("They are : ")
        #for aAng in tRes["remainAngs"]:
        #    print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
        print("Number of deleted angles ", len(tMod["deleted"]["angles"]))
        #print("They are : ")
        #for aAng in tMod["deleted"]["angles"]:
        #    print("Angle formed by atoms %s, %s and %s "%(aAng["atom_id_1"], aAng["atom_id_2"],aAng["atom_id_3"]))
                  
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
                #if aChi["atom_id_centre"] != tRes["atomName"]:
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
         
        self.setDeletedAngInOneResForModification(tRes, tMod)   
        self.setDeleteTorsInOneResForModification(tRes, tMod) 
        self.setDeleteChirsInOneResForModification(tRes, tMod)     
        
                    
    def setDeletedBondInOneResForModification(self, tRes, tMod):
        
        tmpRemBs = []
        for aB in tRes["comp"]["bonds"]:
            tmpRemBs.append(aB)
        
        if len(tMod["added"]["bonds"]) > 0:
            for aB in tMod["added"]["bonds"]:
                tmpRemBs.append(aB)
            
        #print("Number of tmpB is ", len(tmpRemBs))
        
        tRes["remainBonds"] = []
        
        dBIDs = []
        for aDB in tMod["deleted"]["bonds"]:
            
            atm1DB = aDB["atom_id_1"]
            atm2DB = aDB["atom_id_2"]
            print("deleted bond between atom  %s and %s "%(atm1DB, atm2DB))
            aList = [atm1DB, atm2DB]
            aList.sort()
            aStr = aList[0] + "_" + aList[1]
            dBIDs.append(aStr)
        
        if len(dBIDs) > 0:
            print ("deleted bond IDs :")
            print (dBIDs)
        
        for aB in tmpRemBs:
            atm1 = aB["atom_id_1"]
            atm2 = aB["atom_id_2"]
            bList = [atm1, atm2]
            bList.sort()
            bStr = bList[0] + "_" + bList[1]
            
            if not bStr in dBIDs:
                #print("Bond id ", bStr)
                #print("Bond included")
                tRes["remainBonds"].append(aB)
        
        #print("The following bonds are kept:")  
        
        #for aB in tRes["remainBonds"]:
        #    print("Bond between %s and %s of order %s "%(aB["atom_id_1"], aB["atom_id_2"], aB["type"]))
        
        
        # check and add H atom to those atoms involved in deleted bonds
        for aDB in tMod["deleted"]["bonds"]:    
            atm1DB = aDB["atom_id_1"]
            atm2DB = aDB["atom_id_2"] 
            self.addjustAtomInDeletedBondInOneResForModification(tRes, tMod, atm1DB, dBIDs)
            self.addjustAtomInDeletedBondInOneResForModification(tRes, tMod, atm2DB, dBIDs)
    
    
    def setDeletedAngInOneResForModification(self, tRes, tMod):
        
        # check and add H atom to those atoms involved in deleted bonds
        tmpRemAngs = []
        for aAng in tRes["remainAngs"]:
            tmpRemAngs.append(aAng)
        
        tRes["remainAngs"] = []
        
        if len(tMod["deleted"]["bonds"]) > 0:
            for aDB in tMod["deleted"]["bonds"]:    
                atm1DB = aDB["atom_id_1"]
                atm2DB = aDB["atom_id_2"] 
                aDList = [aDB["atom_id_1"], aDB["atom_id_2"]]
                for aAng in tmpRemAngs:
                    #print("Ang among %s %s %s"%(aAng["atom_id_1"], aAng["atom_id_2"], aAng["atom_id_3"]))
                    if aAng["atom_id_2"] in aDList:
                        aIdList = [aAng["atom_id_1"], aAng["atom_id_2"], aAng["atom_id_3"]]
                        if atm1DB in aIdList and atm2DB in aIdList:
                            tMod["deleted"]["angles"].append(aAng)
                        else:
                            tRes["remainAngs"].append(aAng)
                    else:
                        tRes["remainAngs"].append(aAng)
        else:
            for aAng in tmpRemAngs:
                tRes["remainAngs"].append(aAng)
                    
    def setDeleteTorsInOneResForModification(self, tRes, tMod):
        
        # check and add H atom to those atoms involved in deleted bonds
        tmpRemTors = []
        for aTor in tRes["remainTors"]:
            tmpRemTors.append(aTor)
        
        tRes["remainTors"] = []
        if len(tMod["deleted"]["bonds"]) > 0:
            for aDB in tMod["deleted"]["bonds"]:    
                atm1DB = aDB["atom_id_1"]
                atm2DB = aDB["atom_id_2"] 
                for aTor in tmpRemTors:
                    aIdList = [aTor["atom_id_1"], aTor["atom_id_2"], aTor["atom_id_3"],  aTor["atom_id_4"]]
                    if atm1DB in aIdList and atm2DB in aIdList:
                        tMod["deleted"]["tors"].append(aTor)
                    else:
                        tRes["remainTors"].append(aTor) 
        else:
            for aTor in tmpRemTors:
                tRes["remainTors"].append(aTor)
                
    
    def setDeleteChirsInOneResForModification(self, tRes, tMod):
            
        # check and add H atom to those atoms involved in deleted bonds
        tmpRemChirs = []
        for aCh in tRes["remainChirs"] :
            tmpRemChirs.append(aCh)
            
        tRes["remainChirs"] = []
        if len(tMod["deleted"]["bonds"]) > 0:    
            for aDB in tMod["deleted"]["bonds"]:    
                atm1DB = aDB["atom_id_1"]
                atm2DB = aDB["atom_id_2"]
                aDList = [aDB["atom_id_1"], aDB["atom_id_2"]]
                print("deleted bond between %s and %s "%(atm1DB, atm2DB))
                for aCh in tmpRemChirs:
                    aIdList = [aCh["atom_id_1"], aCh["atom_id_2"], aCh["atom_id_3"], aCh["atom_id_centre"]]
                    print("For chir %s  %s %s %s "%(aCh["atom_id_centre"], aCh["atom_id_1"], aCh["atom_id_2"], aCh["atom_id_3"]) )
                    if aCh["atom_id_centre"] in aDList:
                        if atm1DB in aIdList and atm2DB in aIdList:
                            tMod["deleted"]["chirs"].append(aCh)
                            print("It is deleted")
                        else:
                            tRes["remainChirs"].append(aCh)
                            print("It is included")
                    else:
                        tRes["remainChirs"].append(aCh)
                        print("It is included")
        else:
            for aCh in tmpRemChirs:
                tRes["remainChirs"].append(aCh)
                
    def addjustAtomInDeletedBondInOneResForModification(self, tRes, tMod, tAtmId, tDelBAtmIds):
        
        
        changeAtms = []
        if len(tMod["changed"]["atoms"]) > 0:
            for aA in tMod["changed"]["atoms"]:
                changeAtms.append(aA["atom_id"].strip().upper())
        print("Check atom ", tAtmId)
        for aAt in tRes["comp"]["atoms"]:
            if aAt["atom_id"].strip().upper() ==tAtmId:
                
                if aAt["atom_id"] != tRes["atomName"]:
                    
                    addedHs=self.checkAssocHAtoms2(tRes, tMod, aAt, tDelBAtmIds)
                    
                    if addedHs >0 and not tAtmId in changeAtms :
                        print("number of added H is ", addedHs)
                        tMod["changed"]["atoms"].append(aAt)
                        print("the State of Atom %s is changed because of bond deletion"%tAtmId)
    
    
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
                
                if len(aTmpLABonds)==4 and len(connectedH) >0:
                    
                    if "charge" in tRes["comp"]["atoms"][aLAtmSerial]:
                        #tRes["origCharge"][aLAtmId] = tRes["comp"]["atoms"][aLAtmSerial]["charge"]
                        tRes["comp"]["atoms"][aLAtmSerial]["charge"] = "0"
                    elif "form_charge" in tRes["comp"]["atoms"][aLAtmSerial]:
                        #tRes["origCharge"][aLAtmId] =tRes["comp"]["atoms"][aLAtmSerial]["formal_charge"]
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
            
            if len(tMod["changed"]["bonds"]) > 0 or len(tMod["deleted"]["bonds"]) > 0:
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
                if len(tMod["deleted"]["bonds"]) > 0:
                    # Consider the effect of bond-order changes for some bonds
                    aDelBondIds = []
                    for aB in tMod["deleted"]["bonds"]:
                        tBIdList1 = [aB["atom_id_1"], aB["atom_id_2"]]
                        tBIdList1.sort()
                        aStr = tBIdList1[0] + "_" + tBIdList1[1]
                        aDelBondIds.append(aStr)
                    for aB in aTmpLABonds:
                        tBIdList2 = [aB["atom_id_1"], aB["atom_id_2"]]
                        tBIdList2.sort()
                        bStr = tBIdList2[0] + "_" + tBIdList2[1]
                        print(bStr)
                        if not bStr in aDelBondIds:
                            aLABonds.append(aB)
                    
            else:
                for aB in aTmpLABonds:
                    aLABonds.append(aB)
        
            if tMode == 1:
                # Add the linked bond
                aLABonds.append(tLinkBonds[0])
            #print("Number of bonds ", len(aLABonds))
            print("tMode ", tMode) 
            print("===============")
            print("The linked atom %s in residue %s "%(aLAtmId, tRes["name"]))
            if len(aLABonds):
                print("It appears in the following bonds now ")
                for aB in aLABonds:
                    print("Bond between %s and %s of order %s "%(aB["atom_id_1"], aB["atom_id_2"], aB["type"]))

                nTotalVa=self.getTotalBondOrderInOneMmcifAtom(aLAtmId, tRes["comp"]["atoms"], aLABonds)

                print("total Valence is ", nTotalVa)
                print("atom ", aLAtmElem.upper())
                print("Default Valence is ", self.chemCheck.defaultBo[aLAtmElem.upper()])
                aC = 0
                if "charge" in tRes["comp"]["atoms"][ aLAtmSerial].keys():
                    aC = tRes["comp"]["atoms"][ aLAtmSerial]["charge"]
                print("The initial charge is  ", aC)
                if aLAtmElem in self.chemCheck.orgVal:
                    if "charge" in tRes["comp"]["atoms"][ aLAtmSerial]:
                        print(tRes["comp"]["atoms"][ aLAtmSerial]["charge"])
                        allowedBO = self.chemCheck.orgVal[aLAtmElem][0] + int(float(tRes["comp"]["atoms"][ aLAtmSerial]["charge"]))
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
    
    def getBondSetForOneLinkedAtomModRes(self, tAtomId, tAtoms, tBonds, tDelAtomIds, tAddedAtoms, tAddedBonds):

        aAtomSet = []
        aBondSet = []
        #print("atom ", tAtomId)
        #print "Number of bonds ", len(tBonds)
        #print("addedBonds ", tAddedBonds)
        allBonds = []
        for aB in tBonds:
            allBonds.append(aB)
        for aB in tAddedBonds:
            allBonds.append(aB)
        for aBond in allBonds:
            aId1 = aBond["atom_id_1"].strip()
            aId2 = aBond["atom_id_2"].strip()
            #print "bond between atom ", aId1, " and atom  ", aId2
            if aId1 == tAtomId and not aId2 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                if aIdx !=-1:
                    aAtomSet.append(tAtoms[aIdx])
                else:
                    aIdx = self.getAtomById(aId2, tAddedAtoms)
                    aAtomSet.append(tAddedAtoms[aIdx])
                aBondSet.append(aBond)
                print ("a bond found between %s and %s "%(aId1, aId2))
            if aId2 == tAtomId and not aId1 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                if aIdx !=-1:
                    aAtomSet.append(tAtoms[aIdx])
                else:
                    aIdx = self.getAtomById(aId1, tAddedAtoms)
                    aAtomSet.append(tAddedAtoms[aIdx])
                aBondSet.append(aBond)
                print ("a bond found between %s and %s "%(aId1, aId2))
                #print "a bond found "
        #print "Number of bonds found ", len(aBondSet)
        return [aAtomSet, aBondSet] 

    def getBondSetForOneAtomByAlias(self, tAtomId, tAtoms, tBonds, tDelAtomIds):

        aAtomSet = []
        aBondSet = []
        for aBond in tBonds:
            aId1 = aBond["atom_id_1_alias"].strip()
            aId2 = aBond["atom_id_2_alias"].strip()
            if aId1 == tAtomId and not aId2 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
            if aId2 == tAtomId and not aId1 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
        #print("Number of bonds found ", len(aBondSet))
        return [aAtomSet, aBondSet]   
    
    
    def getBondSetOnlyForOneAtomById(self, tAtomId, tBonds, tAtoms):
        
        aBondSet = []
        #print("atom ", tAtomId)
        #print("Number of bonds ", len(tBonds))
        for aBond in tBonds:
            aId1 = aBond["atom_id_1"].strip()
            aIdx1 = self.getOneAtomSerialById(aId1, tAtoms)
            aId2 = aBond["atom_id_2"].strip()
            aIdx2 = self.getOneAtomSerialById(aId2, tAtoms)
            if tAtoms[aIdx1]["type_symbol"].upper() in  self.chemCheck.orgVal\
                and  tAtoms[aIdx2]["type_symbol"].upper() in  self.chemCheck.orgVal:
                if aId1 == tAtomId :
                    aBondSet.append(aBond)
                    #print("bond between atom ", aId1, " and atom  ", aId2)
                    #print("a bond found ")
                if aId2 == tAtomId :
                    aBondSet.append(aBond)
                    #print("bond between atom ", aId1, " and atom  ", aId2)
                    #print("a bond found ")
        return aBondSet  
        
    def getBondSetForOneAtomById(self, tAtomId, tAtoms, tBonds, tDelAtomIds):

        aAtomSet = []
        aBondSet = []
        print("atom ", tAtomId)
        print("Number of bonds ", len(tBonds))
        for aBond in tBonds:
            aId1 = aBond["atom_id_1"].strip()
            aId2 = aBond["atom_id_2"].strip()
            if aId1 == tAtomId and not aId2 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
            if aId2 == tAtomId and not aId1 in tDelAtomIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
        #print("Number of bonds found ", len(aBondSet))
        return [aAtomSet, aBondSet]      
        
    def getBondSetForOneAtomById2(self, tAtomId, tAtoms, tBonds, tDelBAtmIds):

        aAtomSet = []
        aBondSet = []
        #print("atom ", tAtomId)
        #print("Number of bonds ", len(tBonds))
        for aBond in tBonds:
            aId1 = aBond["atom_id_1"].strip()
            aId2 = aBond["atom_id_2"].strip()
            aList = [aId1, aId2]
            aList.sort()
            aStr = aList[0] + "_" + aList[1] 
            if aId1 == tAtomId and not aStr in tDelBAtmIds:
                aIdx = self.getOneAtomSerialById(aId2, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
            if aId2 == tAtomId and not aStr in tDelBAtmIds:
                aIdx = self.getOneAtomSerialById(aId1, tAtoms)
                aAtomSet.append(tAtoms[aIdx])
                aBondSet.append(aBond)
                #print("bond between atom ", aId1, " and atom  ", aId2)
                #print("a bond found ")
        #print("Number of bonds found ", len(aBondSet))
        return [aAtomSet, aBondSet]  
    
    def getTotalBondOrderInOneMmcifAtom(self, tAtomId, tAtoms, tBonds):
        # tBonds is prefilted in getBondSetForOneMmcifAtom() and other
        totalOr = 0
        for aBond in tBonds:
            aId1  = aBond["atom_id_1"].strip()
            aIdx1 = self.getOneAtomSerialById(aId1, tAtoms)
            aId2  = aBond["atom_id_2"].strip()
            aIdx2 = self.getOneAtomSerialById(aId2, tAtoms)
            elem1 = tAtoms[aIdx1]["type_symbol"] 
            elem2 = tAtoms[aIdx2]["type_symbol"]
            if elem1 in self.chemCheck.organicSec and\
                elem2 in self.chemCheck.organicSec: 
                aOr = BondOrderS2N(aBond["type"])
                if aOr != -1:
                    totalOr += aOr
                else:
                    totalOr = 0
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


    def addOneHInRes2(self, tHConnAtom, tOtheAtmSet, tAllAtmIds, tRes, tMod):

        aAtom = {}
        aAtom["comp_id"]     = tRes["name"]
        aAtom["atom_id"]     = self.chemCheck.setHName2(tHConnAtom, tOtheAtmSet, tAllAtmIds)         
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
                tChir["volume_sign"] = "both"
                tResChirs.append(tChir)
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
            #print("Atom id ", aAtom["atom_id"], "   Charge ", aAtom["charge"])
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
            #print("bond ", aBond["comp_id"])
            #print("atom 1 alias ", aBond["atom_id_1_alias"])
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
            #print ("link atom, centre alias : ", tLinkedObj.stdLigand1["linkChir"][0]["atom_id_centre_alias"])
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
                #print(curKey2,  tLinkedObj.stdLigand1["linkChir"][0][curKey2])
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
            #print ("link atom, centre alias : ", tLinkedObj.stdLigand2["linkChir"][0]["atom_id_centre_alias"])
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
                #print(curKey2,  tLinkedObj.stdLigand2["linkChir"][0][curKey2])
            tLinkedObj.combLigand["chirs"].append(tLinkedObj.stdLigand2["linkChir"][0])
        
        self.outAtomNameMapToJSon(tLinkedObj)
    
        
    def outAtomNameMapToJSon(self, tLinkedObj):

        #aOutFN = os.path.join(self.scrDir, tLinkedObj.combLigand["name"] + "_AtomNameMapping.json")
        aOutFN = os.path.join(self.scrDir, "LIG_AtomNameMapping.json")
        try: 
            aOutF = open(aOutFN, "w")
        except IOError:
            print("%s can not be open for writing "%aOutFN)
            self.errLevel    = 32
            if self.errLevel not in self.errMessage:
                self.errMessage[self.errLevel] = []
            self.errMessage[self.errLevel].append("%s can not be open for write "%aOutFN)
        else:       
            resNames =["LIG", tLinkedObj.stdLigand1["name"], tLinkedObj.stdLigand2["name"]]
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

    def getOneAtomSymbById(self, tId, tAtoms):
        
        aRet = ""
        print("seach atom", tId)
        for aAtm in tAtoms:
            #print("current Id", aAtm["atom_id"])
            #print("its symbol ", aAtm["type_symbol"])
            if aAtm["atom_id"]== tId:
                aRet = aAtm["type_symbol"]
                
                break
        print("found symbol is ", aRet)
        return aRet
    
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
               if len(idStrs)==2:
                   newStrs = idStrs
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

            for aBond in tLinkObj.combLigand["bonds"]:
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
            print("-==========")
            print("Atom id is ", atmId)
            print("Atom elem is ", atmElm)
            if tMode == 2:
                atmId = aAtom["atom_id_alias"]
            [bondAtomSet, aLABonds] =  self.getBondSetForOneAtomByAlias(atmId, tMonomer["atoms"], tMonomer["bonds"], tDelAtomIds)
            nTotalVa=self.getTotalBondOrderInOneMmcifAtom(atmId, tMonomer["atoms"], aLABonds) 
            aCharge = 0
            if "charge" in aAtom:
                aCharge = int(float(aAtom["charge"]))
            if aCharge != 0: 
                nTotalVa = nTotalVa - aCharge
            print("its charge is ", aCharge)
            print("nTotalVa is ", nTotalVa)
            if atmElm in self.chemCheck.orgVal and atmElm !="H":
                print("default=", self.chemCheck.orgVal[atmElm][0])
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
         
            tCombId = "LIG"   
            numAts = len(tMonomer["atoms"]) 
            numH   = len(tMonomer["hAtoms"]) 
            numNonH = numAts - numH
            aOutF.write("%s%s%s%s%6d%6d\n"%(tCombId.ljust(8), tCombId.ljust(8), "\'.           \'".ljust(20),\
                                           "non-polymer".ljust(20), numAts, numNonH))
            aOutF.write("data_comp_LIG\n")
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
                    if "charge" in aAtom.keys():
                        aC = str(aAtom["charge"])
                    else:
                        aC = "0"
                    print("%s%s%s%s"%(aAtom["atom_id"].ljust(10), aAtom["atom_id_alias"].ljust(10), 
                                      aAtom["type_symbol"].ljust(6), aC.ljust(6)))
                for aBond in tLinkedObj.combLigand["bonds"]:
                    print("%s%s%s%s%s\n"%(aBond["atom_id_1_alias"].ljust(10), aBond["atom_id_2_alias"].ljust(10),
                                   ("("+aBond["atom_id_1"]).ljust(10), (aBond["atom_id_2"]+ ")").ljust(10),
                                    aBond["type"].ljust(1))) 
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
                        #print("Number of atoms in the comboLigand : ", len(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["atoms"]))
                        comboScrDN = self.subRoot + "_TMP"
                        if os.path.isdir(comboScrDN) and not self.testMode:
                            #print "Delete the tempo dir ", comboScrDN
                            #shutil.rmtree(comboScrDN)
                            pass
                        comboLigMolFileName = self.subRoot + ".mol"
                        self.fileTool.MmCifToMolFile(tLinkedObj.combLigand["outCif"], comboLigMolFileName)
            else:
                print(self.errLevel) 
            
       
    def getChangesInModificationFromCombLigand(self, tLinkedObj):

  
        #for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["atoms"]:
        #    print(aAtom.keys())
        #    print aAtom["atom_id"]
        #    print "Atom %s is in residue %s "%(aAtom["atom_id"], aAtom["res_idx"])
        #    print "Is it added ? ", aAtom["is_added"]
        for aAt in tLinkedObj.stdLigand1["remainAtoms"]:
            if aAt["atom_id"] in tLinkedObj.stdLigand1["origCharge"]:
                aAt["charge"] = tLinkedObj.stdLigand1["origCharge"][aAt["atom_id"]] 
        for aAt in tLinkedObj.stdLigand2["remainAtoms"]:
            if aAt["atom_id"] in tLinkedObj.stdLigand2["origCharge"]:
                aAt["charge"] = tLinkedObj.stdLigand2["origCharge"][aAt["atom_id"]]
             
            
        
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
         
        for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["atoms"]:
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

        if len(tLinkedObj.modLigand1["added"]["atoms"]) > 0:
            for aAtom in tLinkedObj.modLigand1["added"]["atoms"]:
                addedSet1.append(aAtom["atom_id"])
        if len(tLinkedObj.modLigand2["added"]["atoms"]) > 0:
            for aAtom in tLinkedObj.modLigand2["added"]["atoms"]:
                addedSet2.append(aAtom["atom_id"])

        tLinkedObj.modLigand2["changed"]["bonds"] = []
        tLinkedObj.modLigand1["added"]["bonds"]   = []
        tLinkedObj.modLigand2["added"]["bonds"]   = []

        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["bonds"]:
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
                        print("added bond between atom1 %s atom2 %s "%(aBond["atom_id_1"], aBond["atom_id_2"])) 
                        tLinkedObj.modLigand2["added"]["bonds"].append(aBond)
                
        print("Number of changed bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["bonds"]))  
        print("Number of added bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["bonds"]))  
        print("Number of deleted bonds in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["bonds"]))  
        print("Number of changed bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["bonds"]))  
        print("Number of added bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["bonds"]))  
        print("Number of deleted bonds in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["bonds"]))  

        # Angles
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["angles"]:
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
        #for aAng in tLinkedObj.modLigand2["changed"]["angles"]:
        #    id1 = aAng["atom_id_1"]
        #    id2 = aAng["atom_id_2"]
        #    id3 = aAng["atom_id_3"]
        #    print("combo angle between %s - %s -%s "%(id1, id2, id3))
        print("Number of added angles in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["angles"]))  
        print("Number of deleted angles in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["angles"])) 
        
        # Torsions
        torIdNumLig1 = {}
        torIdNumLig1["SP2_SP2"]=0
        torIdNumLig1["SP2_SP3"]=0
        torIdNumLig1["SP3_SP3"]=0
        self.getTorIdSets(torIdNumLig1, tLinkedObj.stdLigand1["remainTors"])
        print("torIdNumLig1: ", torIdNumLig1)
        
        torIdNumLig2 = {}
        torIdNumLig2["SP2_SP2"]=0
        torIdNumLig2["SP2_SP3"]=0
        torIdNumLig2["SP3_SP3"]=0
        self.getTorIdSets(torIdNumLig2, tLinkedObj.stdLigand2["remainTors"])
        print("torIdNumLig2: ", torIdNumLig2)
        
        if "tors" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["tors"]:
                print("A torsion ")
                print("atom %s in residue %d "%(aTor["atom_id_1"], aTor["atom_id_1_resNum"]))   
                print("atom %s in residue %d "%(aTor["atom_id_2"], aTor["atom_id_2_resNum"]))   
                print("atom %s in residue %d "%(aTor["atom_id_3"], aTor["atom_id_3_resNum"]))  
                print("atom %s in residue %d "%(aTor["atom_id_4"], aTor["atom_id_4_resNum"]))  
                if aTor["atom_id_1_resNum"]==1 and aTor["atom_id_2_resNum"]==1\
                   and aTor["atom_id_3_resNum"]==1 and aTor["atom_id_4_resNum"]==1:
                    if not aTor["atom_id_1"] in addedSet1 and not aTor["atom_id_2"] in addedSet1\
                       and not aTor["atom_id_3"] in addedSet1 and not aTor["atom_id_4"] in addedSet1: 
                        self.checkTorMod(tLinkedObj.stdLigand1["remainTors"], aTor, tLinkedObj.modLigand1["changed"]["tors"], torIdNumLig1)
                    else: 
                        tLinkedObj.modLigand1["added"]["tors"].append(aTor)
                elif aTor["atom_id_1_resNum"]==2 and aTor["atom_id_2_resNum"]==2\
                   and aTor["atom_id_3_resNum"]==2 and aTor["atom_id_4_resNum"]==2:
                    if not aTor["atom_id_1"] in addedSet2 and not aTor["atom_id_2"] in addedSet2\
                       and not aTor["atom_id_3"] in addedSet2 and not aTor["atom_id_4"] in addedSet2: 
                        self.checkTorMod(tLinkedObj.stdLigand2["remainTors"], aTor, tLinkedObj.modLigand2["changed"]["tors"], torIdNumLig2)
                    else: 
                        tLinkedObj.modLigand2["added"]["tors"].append(aTor)
                elif aTor["atom_id_1_resNum"]==1 and aTor["atom_id_2_resNum"]==1\
                   and aTor["atom_id_3_resNum"]==1 and aTor["atom_id_4_resNum"]==2:
                       self.checkTorDelId(tLinkedObj.stdLigand1["remainTors"], aTor, tLinkedObj.modLigand1["deleted"]["tors"])
                elif aTor["atom_id_1_resNum"]==2 and aTor["atom_id_2_resNum"]==1\
                   and aTor["atom_id_3_resNum"]==1 and aTor["atom_id_4_resNum"]==1:
                       self.checkTorDelId(tLinkedObj.stdLigand1["remainTors"], aTor, tLinkedObj.modLigand1["deleted"]["tors"])
                elif aTor["atom_id_1_resNum"]==2 and aTor["atom_id_2_resNum"]==2\
                   and aTor["atom_id_3_resNum"]==2 and aTor["atom_id_4_resNum"]==1:
                       self.checkTorDelId(tLinkedObj.stdLigand2["remainTors"], aTor, tLinkedObj.modLigand2["deleted"]["tors"])
                elif aTor["atom_id_1_resNum"]==1 and aTor["atom_id_2_resNum"]==2\
                   and aTor["atom_id_3_resNum"]==2 and aTor["atom_id_4_resNum"]==2:
                       print("Here")
                       self.checkTorDelId(tLinkedObj.stdLigand2["remainTors"], aTor, tLinkedObj.modLigand2["deleted"]["tors"])
    
            print("Number of changed tors in residue 1 is %d "%len(tLinkedObj.modLigand1["changed"]["tors"]))  
            print("Number of added tors in residue 1 is %d "%len(tLinkedObj.modLigand1["added"]["tors"]))  
            print("Number of deleted tors in residue 1 is %d "%len(tLinkedObj.modLigand1["deleted"]["tors"]))  
            print("Number of changed tors in residue 2 is %d "%len(tLinkedObj.modLigand2["changed"]["tors"]))  
            print("Number of added tors in residue 2 is %d "%len(tLinkedObj.modLigand2["added"]["tors"]))  
            print("Number of deleted tors in residue 2 is %d "%len(tLinkedObj.modLigand2["deleted"]["tors"])) 
         
        # Chirs
        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["chirs"]:
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
        delPlanIds = {}
        delPlanIds["1"] = []
        delPlanIds["2"] = []
        self.checkPlModFromCombo(tLinkedObj, delPlanIds)
        
        self.checkPlModFromOrig(tLinkedObj, delPlanIds)
     
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
        print("comp atom id ", aId)
        for aAtm in tOriAtoms:
            if aAtm["atom_id"] == aId:
                print("comp id ", aAtm["atom_id"])
                if self.compare2Atoms(aAtm, tAtom) and not aId in tExistChangedAtoms:
                    print(aId, " changed")
                    tModAtoms.append(tAtom)
                    if not aId in tExistChangedAtoms:
                        tExistChangedAtoms.append(aId)
        
                
   
    def compare2Atoms(self, tOriAtom, tAtom):

        lChange = False
        if "type_energy" in tOriAtom and "type_energy" in tAtom:
            if tOriAtom["type_energy"] !=tAtom["type_energy"]: 
                tOriAtom["type_energy"] =tAtom["type_energy"]
                lChange = True
        print("O: Is charge a key ", "charge" in tOriAtom)
        print("C: Is charge a key ", "charge" in tAtom)
        if "charge" in tOriAtom and "charge" in tAtom:
            print("charge1 ", tOriAtom["charge"])
            print("charge2 ", tAtom["charge"])
            if tOriAtom["charge"] !=tAtom["charge"]:
                lChange = True
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
        #print(tOriBond.keys())
        #print(tBond.keys())
        if tOriBond["type"].upper()[:3] !=tBond["type"].upper()[:3]:
            lChanged = True
        #print("Bond changes between %s and %s"%(tOriBond["atom_id_1"], tOriBond["atom_id_2"]))
        v1 = float(tOriBond["value_dist"]) 
        v2 = float(tBond["value_dist"])
        dV = math.fabs(v1-v2)
        s1 = float(tOriBond["value_dist_esd"])
        s2 = float(tBond["value_dist_esd"])
        dS = math.fabs(s1-s2)
        #print("XXX delta is ", dS)
        #print("0.3*s1 is ", 0.3*s1)
        if  dV > s1 :
            lChanged = True
        elif dS > 0.3*s1:
            lChanged = True
                
        return lChanged
   
    def checkAngMod(self, tOrigAngs, tAng, tModAngs):

        id1 = tAng["atom_id_1"]
        id2 = tAng["atom_id_2"]
        id3 = tAng["atom_id_3"]
        #print("combo angle between %s - %s -%s "%(id1, id2, id3))
        for aAng in tOrigAngs:
            if aAng["atom_id_2"]==id2:
                if (aAng["atom_id_1"]==id1 and aAng["atom_id_3"]==id3) or\
                   (aAng["atom_id_1"]==id3 and aAng["atom_id_3"]==id1):
                    if self.compare2Angs(aAng, tAng):
                        #print("the angle is included")
                        tModAngs.append(tAng)
                    break

    def compare2Angs(self, tOrigAng, tAng):
   
        lChanged = False
        v1 = float(tOrigAng["value_angle"])
        v2 = float(tAng["value_angle"])
        dV = math.fabs(v1-v2)
        s1 = float(tOrigAng["value_angle_esd"])
        s2 = float(tAng["value_angle_esd"])
        dS = math.fabs(s1-s2)
        if  dV > s1:
            lChanged = True
        #elif dS > 0.3*s1:
        #    lChanged = True
            
        return lChanged
    
    def getTorIdSets(self, tIdMaps, tOrigTors):
        
        for aTor in tOrigTors:
            if aTor["id"].upper().find("CONST")==-1:
                if aTor["id"].upper().find("SP2_SP2") !=-1:
                    tIdMaps["SP2_SP2"] +=1
                elif aTor["id"].upper().find("SP2_SP3") !=-1:
                    tIdMaps["SP2_SP3"] +=1
                elif aTor["id"].upper().find("SP3_SP3") !=-1:
                    tIdMaps["SP3_SP3"] +=1  
   
    def checkTorMod(self, tOrigTors, tTor, tModTors, tTorIdNum):

        id1 = tTor["atom_id_1"]
        id2 = tTor["atom_id_2"]
        id3 = tTor["atom_id_3"]
        id4 = tTor["atom_id_4"]
        
        
        
        speTorIds = ["chi1", "chi2", "chi3", "chi4", "chi5", "hh", "hh1", "hh2"]
        
        
        for aTor in tOrigTors:
            
            if not aTor["id"] in speTorIds:
                if (aTor["atom_id_2"]==id2 and aTor["atom_id_3"]==id3):
                    if (aTor["atom_id_1"]==id1 and aTor["atom_id_4"]==id4):
                        if self.compare2Tors(aTor, tTor, tTorIdNum):
                            tModTors.append(tTor)
                        break
                elif (aTor["atom_id_2"]==id3 and aTor["atom_id_3"]==id2):
                    if (aTor["atom_id_1"]==id4 and aTor["atom_id_4"]==id1):
                        if self.compare2Tors(aTor, tTor, tTorIdNum):
                            tModTors.append(tTor)
                        break

    def compare2Tors(self, tOrigTor, tTor, tTorIdNum):
   
        lChanged = False
        if float(tOrigTor["value_angle"]) !=float(tTor["value_angle"]):
            lChanged = True
        #elif float(tOrigTor["value_angle_esd"]) !=float(tTor["value_angle_esd"]):
        #    lChanged = True
        elif tOrigTor["period"] != tTor["period"]:
            lChanged = True
        
        if tOrigTor["id"] != tTor["id"]:
            lIdChanged = self.checkTorId(tOrigTor["id"], tTor["id"])
            if lIdChanged :
                lChanged = lIdChanged
                #print("id before: ", tTor["id"])
                self.setTorId(tTor, tTorIdNum)
                #print("id after : ", tTor["id"])
            
        return lChanged
   
    def checkTorId(self, tId1, tId2):
        
        
        aRet = False
        if (tId1.upper().find("CONST") !=-1 and \
            tId2.upper().find("CONST") ==-1) or \
           (tId2.upper().find("CONST") !=-1 and \
            tId1.upper().find("CONST") ==-1) : 
            aRet = True       
        elif tId1.upper().find("CONST") ==-1 and \
            tId2.upper().find("CONST") ==-1:
            strs1 = tId1.strip().split("_")
            strs2 = tId2.strip().split("_")
            if len(strs1) > 2 and len(strs2) > 2:
                if strs1[0] != strs2[0] or strs1[1] !=strs2[1]:
                    aRet = True
                    
        return aRet
            
    def setTorId(self, tTor,tTorIdNum):
        
        if tTor["id"].upper().find("CONST")==-1:
            if tTor["id"].upper().find("SP2_SP2") !=-1:
                aStr = "%d"%(tTorIdNum["SP2_SP2"]+1)
                tTor["id"] = "sp2_sp2_" + aStr
            elif  tTor["id"].upper().find("SP2_SP3") !=-1:
                aStr = "%d"%(tTorIdNum["SP2_SP3"] +1)
                tTor["id"] = "sp2_sp3_" + aStr
            elif  tTor["id"].upper().find("SP3_SP3") !=-1:
                aStr = "%d"%(tTorIdNum["SP3_SP3"]+1)
                tTor["id"] = "sp3_sp3_" + aStr
                
    def checkTorDelId(self, tOrigTors, tTor, tModTors):
        
        id2 = tTor["atom_id_2"]
        id3 = tTor["atom_id_3"]
        for aTor in tOrigTors:
            
            if ((aTor["atom_id_2"]==id2 and aTor["atom_id_3"]==id3)
                or (aTor["atom_id_2"]==id3 and aTor["atom_id_3"]==id2)):
                # delet the torsion where its bond appearing in the link section
                tModTors.append(aTor)
        
    def checkChiMod(self, tOrigChirs, tChi, tModChi):

        pass

    def compare2Chis(self, tOrigChi, tChi):

        pass

    def checkPlModFromCombo(self, tLinkedObj, tDelPlanIds):
        
        print("Number of the orig planes for lig 1 : ", len(tLinkedObj.stdLigand1["remainPls"]))
        print("Number of the orig planes for lig 2 : ", len(tLinkedObj.stdLigand2["remainPls"]))
        
        

        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"].keys()):
                #print("For plane ", aPl)
                inRes1 = []
                inRes2 = []
                nAtmInPl = len(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl])
                for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl]:
                    #print("atom %s in residue %d "%(aAtom["atom_id"], aAtom["atom_id_resNum"]))
                    if aAtom["atom_id_resNum"]== 1:
                        inRes1.append(aAtom["atom_id"])
                    elif aAtom["atom_id_resNum"]== 2:
                        inRes2.append(aAtom["atom_id"])
                #print("Check plane combo plane ", aPl)
                #print("nAtmInPl=", nAtmInPl)
                #print("inRes1=", len(inRes1))
                #print("inRes2=", len(inRes2))
                if len(inRes1) == nAtmInPl:
                    self.checkPlMod(tLinkedObj.stdLigand1["remainPls"], tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl],\
                                    tLinkedObj.modLigand1["deleted"]["planes"], tLinkedObj.modLigand1["added"]["planes"], tDelPlanIds["1"]) 
                elif len(inRes2) == nAtmInPl:
                    self.checkPlMod(tLinkedObj.stdLigand2["remainPls"], tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl],\
                                    tLinkedObj.modLigand2["deleted"]["planes"], tLinkedObj.modLigand2["added"]["planes"], tDelPlanIds["2"])
     
    def checkPlModFromOrig(self, tLinkedObj, tDelPlaIds):

        #print("check from original ")
        
        for aPl in tLinkedObj.stdLigand1["remainPls"]:
            if not "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
                tLinkedObj.modLigand1["deleted"]["planes"].append(aPl)  
            else:
                if not aPl[0]["plane_id"] in tDelPlaIds["1"]:
                    self.checkPlMod2(aPl, tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"], tLinkedObj.modLigand1["deleted"]["planes"])

        for aPl in tLinkedObj.stdLigand2["remainPls"]:
            if not "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
                tLinkedObj.modLigand2["deleted"]["planes"].append(aPl)  
            else:
                if not aPl[0]["plane_id"] in tDelPlaIds["2"]:
                    self.checkPlMod2(aPl, tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"], tLinkedObj.modLigand2["deleted"]["planes"])

    def checkPlMod(self, tOrigPls,  tPl, tModDelPls, tModAddPls, tDelPlIds):

        
        tPlAtmIds = []
        #print("Atoms in tPl : ")
        for aAtm in tPl:
            #if aAtm["atom_id"][0] != "H":
            tPlAtmIds.append(aAtm["atom_id"])
            #print(aAtm["atom_id"])
        #print("==========================")
        #print("Check the combo pl ", tPl[0]["plane_id"])
        #print("==========================")
        nAtoms = len(tPlAtmIds)
        for aPl in tOrigPls:
            #print("Atoms in aPl : ")
            #print("compare Orig pl ", aPl[0]["plane_id"])
            tOrigPlAtmIds = []
            for aPlAtm in aPl:
                #if aPlAtm["atom_id"][0] != "H":
                tOrigPlAtmIds.append(aPlAtm["atom_id"])
                #print(aPlAtm["atom_id"])
            nOrigAtoms = len(tOrigPlAtmIds)
            if nAtoms > nOrigAtoms: 
                overlapAtms = []
                for aOId in tOrigPlAtmIds:
                    if aOId in tPlAtmIds:
                        overlapAtms.append(aOId)
                #print ("overlaped atoms ",  len(overlapAtms))
                if len(overlapAtms) == nOrigAtoms:
                    aPlId = aPl[0]["plane_id"]
                    if not aPlId in tDelPlIds:
                        tModDelPls.append(aPl) 
                        tModAddPls.append(tPl)
                        tDelPlIds.append(aPlId)
                        #print("Num atom in combPl ", nAtoms )
                        #print("Num atom in origPl ", nOrigAtoms)
                        #print("!!!!orignal %s is deleted"%aPlId)
                        break
            elif nAtoms < nOrigAtoms: 
                overlapAtms = []
                for aOId in tPlAtmIds:
                    if aOId in tOrigPlAtmIds:
                        overlapAtms.append(aOId)
                #print("nOrigAtm ", nOrigAtoms)
                #print("nOverlap ", len(overlapAtms))
                if len(overlapAtms) == nAtoms:
                    #for aId in overlapAtms:
                    #    print ("overlaped atom : %s "%aId)
                    aPlId = aPl[0]["plane_id"]
                    if not aPlId in tDelPlIds:
                        tModDelPls.append(aPl) 
                        tModAddPls.append(tPl)
                        tDelPlIds.append(aPlId)
                        #print("Num atom in combPl ", nAtoms )
                        #print("Num atom in origPl ", nOrigAtoms)
                        #print("!!!!orignal %s is deleted"%aPlId)
                        break
                    else:
                        tModAddPls.append(tPl)
                        print("!!!!combo plan %s is added"%aPlId)
                        break
        """
        print("2. tModDelPls")
        print(tModDelPls)
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
                    #tModAddPls.append(tPl)
                    print("HereP")
                    break
        """
    def checkPlMod2(self, tPl, tComboPls, tModDelPls):

        tPlAtmIds = []
        #print("==========================")
        #print("Check the orig pl ", tPl[0]["plane_id"])
        #print("==========================")
        for aAtm in tPl:
            if aAtm["atom_id"][0] != "H":
                tPlAtmIds.append(aAtm["atom_id"])

        lOverlapPls = False
        for aPl in tComboPls:
            #print("compare combo pl ", aPl[0]["plane_id"])
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
            #print("The orig pl ", tPl[0]["plane_id"])
            
 
            
    def extractOneLinkInfo(self, tLinkedObj):

        self.reIndexCombLigand(tLinkedObj) 
        #self.outComboPdb(tLinkedObj)
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
        for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["atoms"]:
            aAtom["is_added"] = False 
            aAtom["atom_id_alias"] = aAtom["atom_id"]
            lSet = self.getOneOrigAtomFromAlias(tLinkedObj.atomMap, aAtom)
            if not lSet:
                aAtom["is_added"] = True 
                addedAtoms.append(aAtom)

        if len(addedAtoms) !=0:
            for aAtom in addedAtoms:
                self.getResidueIdxFromBonding(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["bonds"], \
                                              tLinkedObj.atomMap, aAtom)
            
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["bonds"]:
            aBond["atom_id_1_alias"] = aBond["atom_id_1"]
            [aBond["atom_id_1_resNum"], aBond["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aBond["atom_id_1_alias"]) 
            aBond["atom_id_2_alias"] = aBond["atom_id_2"]
            [aBond["atom_id_2_resNum"], aBond["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aBond["atom_id_2_alias"]) 

        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["angles"]:
            aAng["atom_id_1_alias"] = aAng["atom_id_1"]
            [aAng["atom_id_1_resNum"], aAng["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_1_alias"]) 
            aAng["atom_id_2_alias"] = aAng["atom_id_2"]
            [aAng["atom_id_2_resNum"], aAng["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_2_alias"]) 
            aAng["atom_id_3_alias"] = aAng["atom_id_3"]
            [aAng["atom_id_3_resNum"], aAng["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aAng["atom_id_3_alias"]) 

        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["tors"]:
            aTor["atom_id_1_alias"] = aTor["atom_id_1"]
            [aTor["atom_id_1_resNum"], aTor["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_1_alias"]) 
            aTor["atom_id_2_alias"] = aTor["atom_id_2"]
            [aTor["atom_id_2_resNum"], aTor["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_2_alias"]) 
            aTor["atom_id_3_alias"] = aTor["atom_id_3"]
            [aTor["atom_id_3_resNum"], aTor["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_3_alias"]) 
            aTor["atom_id_4_alias"] = aTor["atom_id_4"]
            [aTor["atom_id_4_resNum"], aTor["atom_id_4"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aTor["atom_id_4_alias"]) 

        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["chirs"]:
                aChi["atom_id_centre_alias"] = aChi["atom_id_centre"]
                [aChi["atom_id_centre_resNum"], aChi["atom_id_centre"]]=\
                                  self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_centre_alias"])  
                aChi["atom_id_1_alias"] = aChi["atom_id_1"]
                [aChi["atom_id_1_resNum"], aChi["atom_id_1"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_1_alias"]) 
                aChi["atom_id_2_alias"] = aChi["atom_id_2"]
                [aChi["atom_id_2_resNum"], aChi["atom_id_2"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_2_alias"]) 
                aChi["atom_id_3_alias"] = aChi["atom_id_3"]
                [aChi["atom_id_3_resNum"], aChi["atom_id_3"]]=self.getOneOrigAtomIdFromAlias(tLinkedObj.atomMap, aChi["atom_id_3_alias"]) 
        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"].keys()):
                for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl]:
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
                                #print("atom name before: %s"%aAtom.name)
                                if aAtom.name in tLinkedObj.atomMap:
                                    aAtom.name = tLinkedObj.atomMap[aAtom.name][1]
                                    if aAtom.name.count("\"")==2:
                                        aAtom.name=aAtom.name.strip("\"")
                                #print("atom name after: %s"%aAtom.name)
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
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["bonds"]:
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
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["angles"]:
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
        self.selectTors(atm1, atm2, tLinkedObj.cLink["torsions"], tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["tors"])
        #for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["tors"]:
        #    if (aTor["atom_id_2_alias"] == atm1 and aTor["atom_id_3_alias"] == atm2) or\
        #       (aTor["atom_id_2_alias"] == atm2 and aTor["atom_id_3_alias"] == atm1):
        #        tLinkedObj.cLink["torsions"].append(aTor)
        print("Number of Link Torsions ", len(tLinkedObj.cLink["torsions"]))
        print("They are :")
        for aTor in tLinkedObj.cLink["torsions"]:
            print("Torsion between atom %s in residue %d, atom %s in residue %d, atom %s in residue %d, and atom %s in residue %d "\
                 %(aTor["atom_id_1"],aTor["atom_id_1_resNum"], aTor["atom_id_2"], aTor["atom_id_2_resNum"],\
                   aTor["atom_id_3"],aTor["atom_id_3_resNum"], aTor["atom_id_4"], aTor["atom_id_4_resNum"]))
            print("dist_esd ", aTor["value_angle_esd"])
            

        
        tLinkedObj.cLink["chirals"] = [] 
        if "chirs" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["chirs"]:
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

        if "planes" in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]:
            if len(list(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"].keys())) !=0:
                tLinkedObj.cLink["planes"] = {}
                for aPl in list(tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"].keys()): 
                    inAtmIds = []
                    for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl]:
                        if aPlAtm["atom_id_alias"] == atm1 or  aPlAtm["atom_id_alias"]==atm2:
                            inAtmIds.append(aPlAtm["atom_id_alias"])
                    if atm1 in inAtmIds and atm2 in inAtmIds:      
                        tLinkedObj.cLink["planes"][aPl] = [] 
                        for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["planes"][aPl]:
                            tLinkedObj.cLink["planes"][aPl].append(aPlAtm)
        if "planes" in tLinkedObj.cLink:
            print("Number of Link planes ", len(list(tLinkedObj.cLink["planes"].keys())))
            for aPl in sorted(tLinkedObj.cLink["planes"].keys()):
                print("Plane %s contains %d atoms. They are: "%(aPl, len(tLinkedObj.cLink["planes"][aPl])))
                for aAtm in tLinkedObj.cLink["planes"][aPl]:
                    print("Plane %s, atom  %s in residue %d, dist_esd  %s "%(aAtm["plane_id"], aAtm["atom_id"], aAtm["atom_id_resNum"], aAtm["dist_esd"]))

    def selectTors(self, tAtm1Id, tAtm2Id, tLinkTors, tOutCombLigandTors):
       
        doneTB = [] 

        lD1 = False
        lD2 = False
        lD3 = False

        # Two extract loops to keep tors enter in order
                
        for aTor in tOutCombLigandTors:
            if (aTor["atom_id_2_alias"] == tAtm1Id and aTor["atom_id_3_alias"] == tAtm2Id) or\
               (aTor["atom_id_2_alias"] == tAtm2Id and aTor["atom_id_3_alias"] == tAtm1Id):
                if aTor["atom_id_1_alias"][0] !="H" or aTor["atom_id_4_alias"][0] !="H":
                    if not lD1:                  
                        tLinkTors.append(aTor)
                        lD1 = True
        for aTor in tOutCombLigandTors:
            if (aTor["atom_id_1_alias"] == tAtm1Id and aTor["atom_id_2_alias"] == tAtm2Id):
                if aTor["atom_id_3_alias"][0] !="H" or aTor["atom_id_4_alias"][0] !="H":
                    if not lD2:                  
                        tLinkTors.append(aTor)
                        lD2 = True
            elif (aTor["atom_id_4_alias"] == tAtm1Id and aTor["atom_id_3_alias"] == tAtm2Id):
                if aTor["atom_id_2_alias"][0] !="H" and aTor["atom_id_1_alias"][0] !="H":
                    if not lD2:                  
                        tLinkTors.append(aTor)
                        lD2 = True
               
        for aTor in tOutCombLigandTors:
            #print("Torsion angle between %s, %s, %s, %s "%(aTor["atom_id_1"], aTor["atom_id_2"], aTor["atom_id_3"], aTor["atom_id_4"]))
            #print("Torsion id is ", aTor["id"])
            if (aTor["atom_id_1_alias"] == tAtm2Id and aTor["atom_id_2_alias"] == tAtm1Id):
                if aTor["atom_id_3_alias"][0] !="H" or aTor["atom_id_4_alias"][0] !="H":
                    if not lD3:                  
                        tLinkTors.append(aTor)
                        lD3 = True
                    tLinkTors.append(aTor)
            elif (aTor["atom_id_4_alias"] == tAtm2Id and aTor["atom_id_3_alias"] == tAtm1Id):
                if aTor["atom_id_2_alias"][0] !="H" or aTor["atom_id_1_alias"][0] !="H":
                    #print("Torsion angle between %s and %s "%(aTor["atom_id_1"], aTor["atom_id_2"]))
                    if not lD3:                  
                        tLinkTors.append(aTor)
                        lD3 = True
               
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

        #for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["LIG"]["atoms"]:
            
  
    def outVerInfo(self, tOutFile):
       
        tOutFile.write("data_program_info\n\n") 
        tOutFile.write("%s%s\n"%("_acedrg_version".ljust(30),    self.verInfo["ACEDRG_VERSION"].ljust(20)))
        tOutFile.write("%s%s\n"%("_acedrg_db_version".ljust(30), self.verInfo["DATABASE_VERSION"].ljust(20)))
        tOutFile.write("%s%s\n"%("_rdkit_version".ljust(30),    self.verInfo["RDKit_VERSION"].ljust(20)))
        tOutFile.write("%s%s\n"%("_servalcat_version".ljust(30),   self.verInfo["SERVALCAT_VERSION"].ljust(20)))
        if len(self.linkInstructionsContent) > 0:
            tOutFile.write("_CCP4_AceDRG_link_generation.instruction\n")
            tOutFile.write(";\n")
            for aL in self.linkInstructionsContent:
                aL = aL.strip()
                if len(aL):
                    tOutFile.write(aL + "\n")
            tOutFile.write(";\n\n")
        elif len(self.modInstructionsContent) > 0:
            tOutFile.write("_CCP4_AceDRG_Mod_generation.instruction\n")
            tOutFile.write(";\n")
            for aL in self.modInstructionsContent:
                aL = aL.strip()
                if len(aL):
                    tOutFile.write(aL + "\n")
            tOutFile.write(";\n\n")
        
        

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
            aGrpL   = len(aGrp)  + 4
            NA      = tLinkedObj.stdLigand1["list"]["number_atoms_all"]
            NAH     = tLinkedObj.stdLigand1["list"]["number_atoms_nh"]
            aL="%s%s%s%s%s%s\n"%(aNL.ljust(10), aN3.ljust(10), aName.ljust(aNameL), aGrp.ljust(aGrpL), NA.ljust(10), NAH)
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
            aGrpL   = len(aGrp)  + 4
            NA      = tLinkedObj.stdLigand2["list"]["number_atoms_all"]
            NAH     = tLinkedObj.stdLigand2["list"]["number_atoms_nh"]
            aL="%s%s%s%s%s%s\n"%(aNL.ljust(10), aN3.ljust(10), aName.ljust(aNameL), aGrp.ljust(aGrpL), NA.ljust(10), NAH)
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
        aL="%s%s%s%s\n"%(aMN.ljust(10), aName.ljust(aNameL), aLN.ljust(10), aGrp.ljust(20))
        tOutFile.write(aL) 
        
        aMN     = tLinkedObj.modLigand2["name"]
        aName   = tLinkedObj.stdLigand2["list"]["name"]
        aNameL  = len(aName) + 6
        aLN     = tLinkedObj.stdLigand2["name"]
        aGrp    = tLinkedObj.stdLigand2["list"]["group"]
        aL="%s%s%s%s\n"%(aMN.ljust(10), aName.ljust(aNameL), aLN.ljust(10), aGrp.ljust(20))
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
        aG1   = tLinkedObj.stdLigand1["list"]["group"]
        print("Ligand 1 group id %s "%aG1)
        #aG1   = "."
        aLN2  = tLinkedObj.stdLigand2["name"]
        aLMN2 = tLinkedObj.modLigand2["name"]
        aG2   = tLinkedObj.stdLigand2["list"]["group"]
        print("Ligand 2 group id %s "%aG2)
        #aG2   = "."
        aL="%s%s%s%s%s%s%s%s\n"%(aLID.ljust(15), aLN1.ljust(10), aLMN1.ljust(12), aG1.ljust(20),\
                                 aLN2.ljust(10), aLMN2.ljust(12), aG2.ljust(20), aLID.ljust(15))
        tOutFile.write(aL) 
         
        tOutFile.write("\n")


    def outAllComps(self, tOutFile, tLinkedObj):
        if tLinkedObj.stdLigand1["dataBlock"] and tLinkedObj.stdLigand1["outComp"] and tLinkedObj.stdLigand1["compOut"]:
            self.outOneComp(tOutFile, tLinkedObj.stdLigand1)
        if tLinkedObj.stdLigand1["compOut"] and tLinkedObj.stdLigand2["compOut"]:
            if tLinkedObj.stdLigand1["inCif"] == tLinkedObj.stdLigand2["inCif"]:
                tLinkedObj.stdLigand2["compOut"] = False
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
            for aBond in tLinkedObj.stdLigand1["remainBonds"]:
                if aBond["atom_id_1"]==tLinkedObj.stdLigand1["atomName"]:
                    if not aBond["atom_id_2"] in aPoolAtoms1:
                        aPoolAtoms1.append(aBond["atom_id_2"])
                elif  aBond["atom_id_2"]==tLinkedObj.stdLigand1["atomName"]:
                    if not aBond["atom_id_1"] in aPoolAtoms1:
                        aPoolAtoms1.append(aBond["atom_id_1"])
        self.outOneMod(tOutFile, tLinkedObj.modLigand1, tLinkedObj.describLevel, aPoolAtoms1)
        
        aPoolAtoms2 = []
        if "atomName" in tLinkedObj.stdLigand2:
            aPoolAtoms2.append(tLinkedObj.stdLigand2["atomName"])
            for aBond in tLinkedObj.stdLigand2["remainBonds"]:
                if aBond["atom_id_1"]==tLinkedObj.stdLigand2["atomName"]:
                    if not aBond["atom_id_2"] in aPoolAtoms2:
                        aPoolAtoms1.append(aBond["atom_id_2"])
                elif  aBond["atom_id_2"]==tLinkedObj.stdLigand2["atomName"]:
                    if not aBond["atom_id_1"] in aPoolAtoms2:
                        aPoolAtoms2.append(aBond["atom_id_1"])
        
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
                    aC = "0.0"
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
                    aC = "0.0"
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aAtom["atom_id"].ljust(10),\
                                             ".".ljust(10), aAtom["type_symbol"].ljust(10), aAtom["type_energy"].ljust(10),\
                                             aC.ljust(10))
                    tOutFile.write(aL)

            if nAA !=0:
                for aAtom in tModLigand["added"]["atoms"]:
                    aC = "0.0"
                    if "charge" in aAtom:
                        aC = aAtom["charge"]
                    aL = "%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), ".".ljust(10),\
                                             aAtom["atom_id"].ljust(10), aAtom["type_symbol"].ljust(10),\
                                             aAtom["type_energy"].ljust(10), aC.ljust(10))
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
                    if aBond["atom_id_1"] in tPoolAtoms or aBond["atom_id_2"] in tPoolAtoms:
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
            tOutFile.write("_chem_mod_bond.new_value_dist_nucleus\n")
            tOutFile.write("_chem_mod_bond.new_value_dist_nucleus_esd\n")               
         
            if nDB !=0:
                for aBond in tModLigand["deleted"]["bonds"]:
                    aBT = aBond["type"].lower()
                    aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), ".".ljust(15),\
                                             ".".ljust(10), ".".ljust(10), ".".ljust(10))
                    tOutFile.write(aL)
                    
            
            if nCB !=0:
                #print tModLigand["name"]
                for aBond in CB_Bonds:
                    aBT = aBond["type"].lower()
                    #print aBond["atom_id_1"]
                    #print aBond["atom_id_2"]
                    #print aBond["value_dist"]
                    aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), aBond["value_dist"].ljust(15),\
                                             aBond["value_dist_esd"].ljust(10), aBond["value_dist_nucleus"].ljust(15),\
                                             aBond["value_dist_nucleus_esd"].ljust(10))
                    tOutFile.write(aL)
                    
            if nAB !=0:
                for aBond in tModLigand["added"]["bonds"]:
                    #print(aBond.keys())
                    aBT = aBond["type"].lower()
                    if not "value_dist_nucleus" in aBond.keys():
                        print("No n-dist for bond between %s and %s "%(aBond["atom_id_1"].ljust(10), aBond["atom_id_2"].ljust(10)))
                        print("in ligand %s\n"%tModLigand["name"])
                        sys.exit()
                    else:
                        aL = "%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), aBond["atom_id_1"].ljust(10),\
                                             aBond["atom_id_2"].ljust(10), aBT.ljust(15), str(aBond["value_dist"]).ljust(15),\
                                             aBond["value_dist_esd"].ljust(10), aBond["value_dist_nucleus"].ljust(15),\
                                             aBond["value_dist_nucleus_esd"].ljust(10))
                      
                    tOutFile.write(aL)
                    
            
            tOutFile.write("\n")

        # Angles, depending on  the description level 

        nDAngs  = len(tModLigand["deleted"]["angles"])
        nCAngs  = len(tModLigand["changed"]["angles"])
        CA_Angs = []
        nAAngs  = len(tModLigand["added"]["angles"])
        
        #if len(tModLigand["changed"]["angles"]) != 0:
        #    for aAng in tModLigand["changed"]["angles"]:
        #        if tLevel ==1:
        #            if aAng["atom_id_1"] in tPoolAtoms or aAng["atom_id_2"] in tPoolAtoms or aAng["atom_id_3"] in tPoolAtoms:
        #                 CA_Angs.append(aAng)
        #    nCAngs = len(CA_Angs)
        
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
                for aAng in tModLigand["changed"]["angles"]:
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
                    if aTor["atom_id_2"] in tPoolAtoms or aTor["atom_id_3"] in tPoolAtoms:
                       #and not aTor["atom_id_3"] in tPoolAtoms and not aTor["atom_id_4"] in tPoolAtoms:
                       # CT_Tors.append(aTor)
                       #elif aTor["atom_id_3"] in tPoolAtoms and not aTor["atom_id_1"] in tPoolAtoms\
                       #and not aTor["atom_id_2"] in tPoolAtoms and not aTor["atom_id_4"] in tPoolAtoms:
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
            tOutFile.write("_chem_mod_tor.id\n")
            tOutFile.write("_chem_mod_tor.new_value_angle\n")
            tOutFile.write("_chem_mod_tor.new_value_angle_esd\n")
            tOutFile.write("_chem_mod_tor.new_period\n")
   
            if nDTors !=0: 
                for aTor in tModLigand["deleted"]["tors"]:
                    aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              ".".ljust(12), ".".ljust(12), ".".ljust(12), str(aTor["period"]).ljust(15))                  
                    tOutFile.write(aL)
           
            if nCTors !=0: 
                for aTor in CT_Tors:
                    aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "change".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              aTor["id"].ljust(18),\
                                              aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15),\
                                              aTor["period"].ljust(6))                  
                    tOutFile.write(aL)
   
            if nATors !=0: 
                for aTor in tModLigand["added"]["tors"]:
                    aL ="%s%s%s%s%s%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15), \
                                              aTor["atom_id_1"].ljust(10), aTor["atom_id_2"].ljust(10),\
                                              aTor["atom_id_3"].ljust(10), aTor["atom_id_4"].ljust(10),\
                                              aTor["id"].ljust(18),\
                                              aTor["value_angle"].ljust(15), aTor["value_angle_esd"].ljust(15),
                                              aTor["period"].ljust(6))                  
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
            tOutFile.write("_chem_mod_plane_atom.mod_id\n")
            tOutFile.write("_chem_mod_plane_atom.function\n")
            tOutFile.write("_chem_mod_plane_atom.plane_id\n")
            tOutFile.write("_chem_mod_plane_atom.atom_id\n")
            tOutFile.write("_chem_mod_plane_atom.new_dist_esd\n")

            if nDPls: 
                for aPl in tModLigand["deleted"]["planes"]:
                    for aPlAtm in aPl:
                        aL ="%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "delete".ljust(15),\
                            aPlAtm["plane_id"].ljust(15), aPlAtm["atom_id"].ljust(15), aPlAtm["dist_esd"]) 
                        tOutFile.write(aL)
            if nAPls:
                for aPl in tModLigand["added"]["planes"]:
                    for aPlAtm in aPl:
                        aL ="%s%s%s%s%s\n"%(tModLigand["name"].ljust(15), "add".ljust(15),\
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
                     
            
            

