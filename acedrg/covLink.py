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

from exebase       import CExeCode

from acedrgRDKit   import AcedrgRDKit
from filetools     import FileTransformer
from filetools     import Ccp4MmCifObj
from chem          import ChemCheck

from utility       import listComp
from utility       import listComp2
from utility       import listCompDes
from utility       import listCompAcd
from utility       import setBoolDict
from utility       import splitLineSpa

#################################################   

class CovLink:

    def __init__(self):
        
        #keyWordList = ["RES-NAME-1", "FILE-1"(optional), "ATOM-NAME-1",\
        #               "RES-NAME-2", "FILE-1"(optional), "ATOM-NAME-2", "BOND-TYPE",\
        #                "DELETE", ]

        self.stdLigand1                     = {}
        self.stdLigand1["fromScr"]          = False
        self.stdLigand1["useIn"]            = False
        self.stdLigand1["comp"]             = None
        self.stdLigand1["hAtom"]            = []
        self.stdLigand1["outComp"]          = False
        self.stdLigand1["remainAtoms"]      = []
        self.stdLigand1["remainNumHs"]      = []
        self.stdLigand1["remainBonds"]      = []
        self.stdLigand1["remainAngs"]       = []
        self.stdLigand1["remainTors"]       = []
        self.stdLigand1["remainChirals"]    = []
        self.stdLigand2["remainPls"]        = []
        self.stdLigand2                     = {}
        self.stdLigand2["fromScr"]          = False
        self.stdLigand2["useIn"]            = False
        self.stdLigand2["comp"]             = None
        self.stdLigand2["hAtom"]            = []
        self.stdLigand2["outComp"]          = False
        self.stdLigand2["remainAtoms"]      = []
        self.stdLigand2["remainNumHs"]      = []
        self.stdLigand2["remainBonds"]      = []
        self.stdLigand2["remainAngs"]       = []
        self.stdLigand2["remainTors"]       = []
        self.stdLigand2["remainChirals"]    = []
        self.stdLigand2["remainPls"]        = []
        self.suggestBondOrder               = "SINGLE"
        self.combLigand                     = {}
        self.combLigand["atoms"]            = {}
        self.combLigand["bonds"]            = {}
        self.combLigand["chirs"]            = {}
        self.outCombLigand                  = {}
        self.modLigand1                     = {}
        self.modLigand1["outComp"]          = True
        self.modLigand1["deleted"]          = {}
        self.modLigand1["deleted"]["atoms"]   = []
        self.modLigand1["deleted"]["bonds"]   = []
        self.modLigand1["deleted"]["angles"]  = []
        self.modLigand1["deleted"]["tors"]   = []
        self.modLigand1["deleted"]["chirs"]   = []
        self.modLigand1["deleted"]["planes"]  = []
        self.modLigand1["changed"]["atoms"]   = []
        self.modLigand1["changed"]["bonds"]   = []
        self.modLigand1["changed"]["angles"]  = []
        self.modLigand1["changed"]["tors"]   = []
        self.modLigand1["changed"]["chirs"]   = []
        self.modLigand1["changed"]["planes"]  = []
        self.modLigand1["added"]["atoms"]   = []
        self.modLigand1["added"]["bonds"]   = []
        self.modLigand1["added"]["angles"]  = []
        self.modLigand1["added"]["tors"]   = []
        self.modLigand1["added"]["chirs"]   = []
        self.modLigand1["added"]["planes"]  = []
        self.modLigand2                       = {}
        self.modLigand2["deleted"]            = {}
        self.modLigand2["deleted"]["atoms"]   = []
        self.modLigand2["deleted"]["bonds"]   = []
        self.modLigand2["deleted"]["tors"]   = []
        self.modLigand2["deleted"]["angles"]  = []
        self.modLigand2["deleted"]["chirs"]   = []
        self.modLigand2["deleted"]["planes"]  = []
        self.modLigand2["changed"]["atoms"]   = []
        self.modLigand2["changed"]["bonds"]   = []
        self.modLigand2["changed"]["angles"]  = []
        self.modLigand2["changed"]["tors"]    = []
        self.modLigand2["changed"]["chirs"]   = []
        self.modLigand2["changed"]["planes"]  = []
        self.modLigand2["added"]["atoms"]   = []
        self.modLigand2["added"]["bonds"]   = []
        self.modLigand2["added"]["angles"]  = []
        self.modLigand2["added"]["tors"]   = []
        self.modLigand2["added"]["chirs"]   = []
        self.modLigand2["added"]["planes"]  = []
        self.modLigand2["outComp"]            = True
        self.cLink                            = {}
        self.atomMap                          = {}
        self.bondMap                          = {}
        
        self.delSections                      = []

    def checkInPara(self):
   
        aReturn = True

        if not self.stdLigand1.has_key("name") or not self.stdLigand1.has_key("atomName")\
           or not self.stdLigand2.has_key("name") or not self.stdLigand2.has_key("atomName"):
            aReturn = False
            
            if not self.stdLigand1.has_key("name"):
                print "Residue 1 is not named "
            if not self.stdLigand1.has_key("atomName"):
                print "Linked atom in Residue 1 is not named "
            if not self.stdLigand2.has_key("name"):
                print "Residue 2 is not named "
            if not self.stdLigand2.has_key("atomName"):
                print "Linked atom in Residue 2 is not named "
   
        return aReturn 

    def filterAtomsAndBonds(self, tFTool, tDS, tMonomer):

        for aAtom in tFTool.atoms:
            if aAtom.has_key("_chem_comp_atom.atom_id"):
                if aAtom["_chem_comp_atom.atom_id"].upper() != tDS["ATOM-NAME"].upper():
                    tMonomer["remainAtoms"].append(aAtom)
                else:
                    # delete all bonds to the deleted atom
                    for aBond in tFTool.bonds:
                        if aBond.has_key("_chem_comp_bond.atom_id_1")\
                           and aBond.has_key("_chem_comp_bond.atom_id_2"):
                            if aBond["_chem_comp_bond.atom_id_1"].upper() != tDS["ATOM-NAME"].upper()\
                               and aBond["_chem_comp_bond.atom_id_2"].upper() != tDS["ATOM-NAME"].upper():
                                tMonomer["remainBonds"].append(aBond)


class CovLinkGenerator(CExeCode):

    def __init__(self, tLinkInstructions, tScrDir, tOutRoot, tVerInfo=None):

        self.verInfo          =  tVerInfo
        #print self.verInfo

        self.workMode         =  0   # Mode means
                                     # 0     generate a link without optimizing coordinates of 2 input monomers
                                     # 1     generate a linked system which optimizes everything. The input monomers 
                                     #       could be in form of mmcif and SMILES strings
        
        self.errMessage       = []
        self.errLevel         = 0
        # 0 = OK  
        # 1 = No instructions file for building links
        # 2 = Instructions are not validated, no links will be built
        # 3 = Error in building individaul residues
        # 4 = Error in building the combined ligand
        # 5 = Error in extracting the final link info

        # TEMPO
        if not os.environ.has_key("CCP4"):
            print "You need to install or setup CCP4 first "
            sys.exit() 
       
        self.allChemCombDir    = os.getenv("CLIBD_MON")
        #self.allChemCombDir   = "/Users/flong/DB/PDB_related/PDB_Ligands/Cif"

        self.scrDir           = tScrDir
        self.subRoot          = ""
        self.outRoot          = tOutRoot
        self.linkInstructions = tLinkInstructions
        self.ligSrcDir        = ""

        self.cLinks           = []

        # engs 
        self.chemCheck        = ChemCheck()

        if os.path.isfile(self.linkInstructions):
            # input from a file
            self.getInstructionsForLink()
        else:
            print "%s can not be found for reading"%self.linkInstructions
            self.errLevel = 1

        if not self.errLevel:
            if len(self.cLinks):
                for aLink in self.cLinks:
                    self.processOneLink(aLink)
            else:
                self.errMessage.append("No links to be found in the instruction file")
                print self.errMessage
                self.errLevel  = 2
        else:
            print self.errMessage

    def getInstructionsForLink(self):

        aFullName=os.path.basename(self.linkInstructions)
        aExt     = ""
        if aFullName.find(".") !=-1:
            aExt = aFullName.strip().split(".")[-1].strip()
            if aExt.find("cif") !=-1:
                self.getInstructionsForLinkFromCif()
            else:
                # Tempo comment off free format at the moment
                #self.getInstructionsForLinkFreeFormat()
                print "Please use cif format for link instruction at the moment"
                self.errMassage += "Please use cif format for link instruction at the moment"
                self.errLevel = 1

            
    def getInstructionsForLinkFromCif(self):
        
        sys.exit()
        aInsObj = Ccp4MmCifObj(self.linkInstructions)

               

    def getInstructionsForLinkFreeFormat(self):

        allLs =""
        aList = []

        try :
            aInsF = open(self.linkInstructions, "r")
        except IOError:
            print "% can not be open for reading ! "%self.linkInstructions
        else:
            allLs = aInsF.readlines()
            aInsF.close()
            for aL in allLs:
                strs = aL.strip().split()
                if len(strs) !=0:
                    for aStr in strs:
                        if aStr.find("LINK:") == -1:
                            aList.append(aStr)
        print len(aList)
        print aList

        if len(aList):
            lDel = False
            aDS  = {}
            i = 0
            aLink = CovLink()
           
            # Some default values. They will be changed at different stages
 
            aLink.stdLigand1["userIn"]  = False 
            aLink.stdLigand1["OutComp"] = False 
            aLink.stdLigand1["OutModi"] = True 
            aLink.stdLigand2["userIn"]  = False 
            aLink.stdLigand2["OutComp"] = False 
            aLink.stdLigand2["OutModi"] = True 

	    aLink.stdLigand1["Group"]   = "."
            aLink.stdLigand2["Group"]   = "."

            nL = len(aList)
            while i < nL :
                if not lDel:
                    if (i+1) < nL:
                        if aList[i].find("RES-NAME-1") != -1:
                            aLink.stdLigand1["name"] = aList[i+1]
                            i +=2
                        elif aList[i].find("FILE-1") != -1:
                            aLink.stdLigand1["inCif"] = aList[i+1]
                            aLink.stdLigand1["userIn"]     = True 
                            i +=2
                        elif aList[i].find("ATOM-NAME-1") != -1:
                            aLink.stdLigand1["atomName"] = aList[i+1] 
                            i +=2
                        elif aList[i].find("RES-NAME-2") != -1:
                            aLink.stdLigand2["name"] = aList[i+1]
                            idxR = 2 
                            i +=2
                        elif aList[i].find("FILE-2") != -1:
                            aLink.stdLigand2["inCif"] = aList[i+1]
                            aLink.stdLigand2["useIn"]     = True 
                        elif aList[i].find("ATOM-NAME-2") != -1:
                            aLink.stdLigand2["atomName"] = aList[i+1] 
                            i +=2
                        elif aList[i].find("BOND-TYPE") != -1 :
                            aLink.suggestBondOrder = aList[i+1]
                            i+=2
                        elif aList[i].find("DELETE") != -1:
                            lDel = True
                            if len(aDS.keys()) != 0:
                                aLink.delSections.apppend(aDS)
                            i +=1
                    else:
                        print "Error : check your link instruction file"
                        print allLs
                        self.errLevel = 2
                        break
                else:
                    if aList[i].find("ATOM") != -1 and aList[i].find("ATOM-NAME") == -1:
                        if i+2 < nL:
                            if len(aDS.keys()) != 0:
                                aLink.delSections.apppend(aDS)
                                aDS = {}
                            aDS["atomName"] = aList[i+1]
                            aDS["inRes"]     = int(aList[i+2])
                            i+=3
                        else:
                            print "Error Delete Atom section: check your link instruction file"
                            print allLs
                            self.errLevel = 2
                            break

                    elif aList[i].find("BOND-TYPE") != -1 :
                        if i+1 < nL:
                            aDS["bondOrder"] = aList[i+1]
                            i+=2
                    else:
                        print "Error Delete bond section: check your link instruction file"
                        print allLs
                        self.errLevel = 1
                        break

            if len(aDS.keys()) != 0:
                aLink.delSections.append(aDS)
             
            aNS = aLink.stdLigand1["name"][0].lower()
            aNL = aLink.stdLigand1["name"].upper()
            aLink.modLigand1["name"] = aNL + "mod1"
            if not aLink.stdLigand1["userIn"]:
                aLink.stdLigand1["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")

            aNS = aLink.stdLigand2["name"][0].lower()
            aNL = aLink.stdLigand2["name"].upper()
            aLink.modLigand2["name"] = aNL + "mod2"
            if not aLink.stdLigand2["userIn"]:
                aLink.stdLigand2["inCif"] = os.path.join(self.allChemCombDir, aNS, aNL + ".cif")

            if aLink.checkInPara():
                self.cLinks.append(aLink)
                print "Instructions for build a link are  "
                print "Link will happen between "
                print "Atom : %s in Residue %s "%(aLink.stdLigand1["atomName"], aLink.stdLigand1["name"])
                print " and " 
                print "Atom : %s in Residue %s "%(aLink.stdLigand2["atomName"], aLink.stdLigand2["name"])
                if aLink.suggestBondOrder != "":
                    print "The suggested bond-order between those two atom is %s "%aLink.suggestBondOrder 
                print "Two input cif files are : "
                print "%s and "%aLink.stdLigand1["inCif"]
                print "%s"%aLink.stdLigand2["inCif"]
                if len(aLink.delSections):
                    print "The following items are deleted."
                    for aDS in aLink.delSections:
                        if aDS.has_key("inRes") and aDS.has_key("atomName"):
                            print "Atom %s in Residue %d "%(aDS["atomName"], aDS["inRes"]) 
                        if aDS.has_key("bondOrder"):
                            print "%s bond for the above atom "%aDS["BOND-TYPE"] 
            else:
                self.errMessage += "Information in the instruction file is not enough to build a link. \n"
                self.errMessage += "Please check your instruction file to see names of monomers, atoms etc are there\n"
                self.errLevel    = 2
      
        else:
            print "A empty set of input instructions for link building"  
   
    def processOneLink(self, tLinkIns):

        if not self.errLevel:
            if not tLinkIns.stdLigand1["fromScr"]:
                self.setOneCompFromCif(tLinkIns.stdLigand1["inCif"], tLinkIns.stdLigand1)    
                print "Comp 1 ", tLinkIns.stdLigand1["name"], " contains " 
                tLinkIns.stdLigand1["ccp4MmCifObj"].printOneComp(tLinkIns.stdLigand1["name"])
            else:
                self.setOneMonomer(tLinkIns.stdLigand1)

            if not tLinkIns.stdLigand2["fromScr"]:
                self.setOneCompFromCif(tLinkIns.stdLigand2["inCif"], tLinkIns.stdLigand2)    
                print "Comp 2 ", tLinkIns.stdLigand2["name"], " contains " 
                tLinkIns.stdLigand2["ccp4MmCifObj"].printOneComp(tLinkIns.stdLigand2["name"])
            else:
                self.setOneMonomer(tLinkIns.stdLigand2)

            if not self.errLevel: 
                self.buildComboLigand(tLinkIns)

            if not self.errLevel:
                self.extractOneLinkInfo(tLinkIns)

        # Finally print out 
        if not self.errLevel:
            self.outLinkInfo(tLinkIns)

    def setOneCompFromCif(self, tFileName, tMonomer):
     
        # Using the input file or the file in ccp4 monomer lib as it is
        aMmcifObj = Ccp4MmCifObj(tFileName) 
        if aMmcifObj["comps"].has_key(tMonomer["name"]):
            tMonomer["comp"] = aMmcifObj["ccp4CifObj"]["comps"][tMonomer["name"]]
            self.selectHAtoms(tMonomer["comp"])    
        else:
            #The input file does not contain the comp, try $CLIBD_MON
            if tMonomer["useIn"]:
                aSub = tMonomer["name"][0].lower()
                aNewCif = os.path.join(self.allChemCombDir, aSub, tMonomer["name"].upper() + ".cif")
                if os.path.isfile(aNewCif):
                    tMonomer["comp"]  = Ccp4MmCifObj(aNewCif)["ccp4CifObj"]["comps"][tMonomer["name"]]
                    self.selectHAtoms(tMonomer["comp"])    
                    tMonomer["inCif"] = aNewCif 
                else:
                    self.errMassage+="Comp %s can not be found in both %s and %s\n"%(tLinkIns.stdLigand1["name"],\
                                      tLinkIns.stdLigand1["inCif"], aNewCif)
                    self.errLevel = 3
            else:        
                self.errMassage+="No file in ccp4 monomer lib contains  %s \n"%tLinkIns.stdLigand1["name"]
                self.errLevel = 3
           
    def selectHAtoms(self, tMonomer):

        tMonomer["hAtoms"] = []
        
        for aAtom in tMonomer["atoms"]:
            if aAtom["type_symbol"] == "H":
                tMonomer["hAtoms"].append(aAtom["atom_id"])
        
    def setOneMonomer(self, tMonomer):

        if os.path.isfile(tMonomer["inCif"]):
            aNL = tMonomer["name"].upper()
            self._log_name  = os.path.join(self.scrDir, aNL + "_for_link.log")
            self.subRoot    = os.path.join(self.scrDir, aNL + "_for_link")
            self._cmdline   = "acedrg -c %s  -r %s -o %s "%(tMonomer["inCif"], aNL, self.subRoot)   
            self.runExitCode = self.subExecute()
            if not self.runExitCode :
                aOutLigCif = self.subRoot + ".cif"
                if os.path.isfile(aOutLigCif): 
                    tMonomer["outCif"] = aOutLigCif
                else: 
                    self.errMessage.append("Error : The output dictionary for %s does not exist\n"%aNL) 
                    self.errLevel = 3
            else: 
                self.errMessage.append("Error (run time): generate a dictionary for %s does failed\n"%aNL) 
                self.errLevel = 3

    def adjustAtomsAndOthersForComboLigand(self, tLinkedObj):

        self.selectHAtoms(tLinkedObj.stdLigand1["comp"])
        self.selectHAtoms(tLinkedObj.stdLigand2["comp"])

        self.setDeletedAtomsForModification(tLinkedObj) 

        if len(tLinkedObj.modLigand1["deleted"]["atoms"]) !=0:
             self.setDeletedOthersOneResidueForModification(tLinkedObj.stdLigand1, tLinkedObj.modLigand1)    
        if len(tLinkedObj.modLigand2["deleted"]["atoms"]) !=0:
             self.setDeletedOthersOneResidueForModification(tLinkedObj.stdLigand2, tLinkedObj.modLigand2)    
        
    def setDeletedAtomsForModification(self, tLinkedObj):

        # Delete atoms as instructed 
        for aDS in tLinkedObj.delSections:
            if aDS.has_key("InRes"):
                aRes = {}
                if aDB["InRes"] == 1:
                    aRes = tLinkedObj.stdLigand1
                    aMod = tLinkedObj.modLigand1
                elif aDB["InRes"] == 2:
                    aRes = tLinkedObj.stdLigand2
                    aMod = tLinkedObj.modLigand2
                else:
                    errMessage += "Error in setup residue number in the section for DELETE in the input file\n"
                    errLevel    = 2 
                
                if not errLevel:
                    if aDS.has_key("ATOM-NAME"):
                        aName = aDS["ATOM-NAME"].upper()
                        connectedH = []
                        self.getHAtomConnected(aName, aRes, connectedH)
                        for aAtom in aRes["comp"]["atoms"]:
                            if aAtom.has_key("atom_id"):
                                if aAtom["atom_id"].upper() != aName and not aAtom["atom_id"] in connectedH:
                                    aRes["remainAtoms"].append(aAtom)
                                else:
                                    aMod["deleted"]["atoms"].append(aAtom)
                        
    def getHAtomConnected(self, tNonHAtomName, tMonomer, tHNames):

          
        for aBond in tMonomer["comp"]["bonds"]:
            if aBond["atom_id_1"].upper()==tNonHAtomName and aBond["atom_id_2"].upper() in tMonomer["hAtoms"]:
                tHNames.append(aBond["atom_id_2"])
            if aBond["atom_id_2"].upper()==tNonHAtomName and aBond["atom_id_1"].upper() in tMonomer["hAtoms"]:
                tHNames.append(aBond["atom_id_1"])
                
    def setDeletedOthersOneResForModification(self, tLinkedObj, tRes, tMod):

        delAtomIdSet =[]
        for aAtom in tMod["deleted"]["atoms"]:
            delAtomIdSet.append(aAtom["atom_id"])
        print "Deleted atoms in one residue ", delAtomIdSet

        # Delete all bonds that contains the deleted atom
        for aBond in tRes["comp"]["bonds"]:
            if (aBond["atom_id_1"].upper() in delAtomIdSet) or (aBond["atom_id_2"].upper() in delAtomIdSet):
                tMod["deletedBonds"].append(aBond)
            else:
                tRes["remainBonds"].append(aBond)
                      
        # Delete all angles that contain the deleted atom
        for aAng in tRes["comp"]["angles"]:
            if (aAng["atom_id_1"].upper() in delAtomIdSet)\
               or (aAng["atom_id_2"].upper() in delAtomIdSet)\
               or (aAng["atom_id_3"].upper() in delAtomIdSet):
                tMod["deletedAngs"].append(aAng)
            else:
                tRes["remainAngs"].append(aAng)
        
        # Delete all tors that contain the deleted atom
        for aTor in tRes["comp"]["tors"]:
            if (aTor["atom_id_1"].upper() in delAtomIdSet)\
               or (aTor["atom_id_2"].upper() in delAtomIdSet)\
               or (aTor["atom_id_3"].upper() in delAtomIdSet)\
               or (aTor["atom_id_4"].upper() in delAtomIdSet):
                tMod["deletedTors"].append(aTor)
            else:
                tRes["remainTors"].append(aTor)
                       
        # Delete all chiral centers that contain the deleted atom
        for aChi in tRes["comp"]["chirs"]:
            if (aChi["atom_id_centre"].upper() in delAtomIdSet)\
               or (aChi["atom_id_1"].upper() in delAtomIdSet)\
               or (aChi["atom_id_2"].upper() in delAtomIdSet)\
               or (aChi["atom_id_3"].upper() in delAtomIdSet):
                tMod["deletedChirs"].append(aChi)
            else:
                tRes["remainChirs"].append(aChi)
                       
        # Delete the deleted atom from a plane and even delete a plane 
        # if number of atoms in a plane smaller then 3 (after deleting the assigned atom)
        for aPl in tRes["comp"]["planes"].keys():
            aPlGrp = []
            for aPAtm in tRes["comp"]["planes"][aPl]:
                if not aPAtm["atom_id"] in delAtomIdSet:
                    aPlGrp.append(aPAtm)
            if len(aPlGrp) > 2:
                tRes["remainPls"].append(aPlGrp)
            else:
                tMod["deletedPls"].append(aPlGrp)    
        
    def adjustChiralCenterForComboLigand(self, tLinkedObj):
        
        for aDS in tLinkedObj.delSections:
            if aDS.has_key("InRes"):
                aRes = {}
                if aDB["InRes"] == 1:
                    aRes = tLinkedObj.stdLigand1
                    otherRes = tLinkedObj.stdLigand2
                else:
                    aRes = tLinkedObj.stdLigand2
                    otherRes = tLinkedObj.stdLigand1

                if aDS.has_key("ATOM-NAME"):
                    AName = aDS["ATOM-NAME"].upper()    
                    # Filter the chiral center of the deleted atom, or
                    # replace the delected atom with the link one in the other monomer
                    for aChi in  aRes["comp"]["chirs"]:
                        
                        if aChi["atom_id_centre"] != AName:
                            if aChi["atom_id_centre"] == aRes["AtomName"]:
                                lChange = False
                                if aChi["atom_id_1"] == AName:    
                                   aChi["atom_id_1"] = otherRes["AtomName"]
                                   lChange = True
                                elif aChi["atom_id_2"] == AName:    
                                   aChi["atom_id_2"] = otherRes["AtomName"]
                                   lChange = True
                                elif aChi["atom_id_3"] == AName:    
                                   aChi["atom_id_3"] = otherRes["AtomName"]
                                   lChange = True
                                if lChange:
                                    aChi["volume_sign"] = "both"
                                aRes["remainChirs"].append(aChir)
                            else:
                                aRes["remainChirs"].append(aChir)

        
    def setInitComboLigand(self, tLinkedObj):

        # All atoms and bonds in the combo-ligand in the linkedObj
 
        tLinkedObj.combLigand["atoms"] = []
        idxAtmComb = 0
        idxA1      = 0
        idxA2      = 0

        for aAtom in tLinkedObj.stdLigand1["remainAtoms"]:
            self.atomMap[idxAtmComb] = [1, idxA1]
            aAtom["idx"]             = idxA1
            aAtom["idxCombo"]        = idxAtmComb
            if aAtom["type_symbol"]=="H":
                tLinkedObj.stdLigand1["remainIdxHs"].append(idxA1)
            idxA1 +=1
            idxAtmComb +=1
            aAtom["atom_id_alias"] = aAtom["type_symbol"] + str(idxAtmComb)   
            tLinkedObj.combLigand["atoms"].append(aAtom)

        for aAtom in tLinkedObj.stdLigand2["remainAtoms"]:
            self.atomMap[idxAtmComb] = [2, idxA2]
            aAtom["idx"]             = idxA2
            aAtom["idxCombo"]        = idxAtmComb
            if aAtom["type_symbol"]=="H":
                tLinkedObj.stdLigand2["remainIdxHs"].append(idxA2)
            idxA2 +=1
            idxAtmComb +=1   
            aAtom["atom_id_alias"] = aAtom["type_symbol"] + str(idxAtmComb)   
            tLinkedObj.combLigand["atoms"].append(aAtom)
        
        tLinkedObj.combLigand["bonds"] = []
        for aBond in tLinkedObj.stdLigand1["remainBonds"]:
            aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
            aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand1["remainAtoms"])
            if aBond["atom_id_1_alias"] == "":
                self.errMassage += "An error happens in residue 1:\n"
                self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"] 
                self.errLevel    = 3
            if aBond["atom_id_2_alias"] == "":
                self.errMassage += "An error happens in residue 1:\n"
                self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
                self.errLevel    = 3
            tLinkedObj.combLigand["bonds"].append(aBond)
        for aBond in tLinkedObj.stdLigand2["remainBonds"]:
            aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand2["remainAtoms"])
            aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
            if aBond["atom_id_1_alias"] == "":
                self.errMassage += "An error happens in residue 2:\n"
                self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"] 
                self.errLevel    = 3
            if aBond["atom_id_2_alias"] == "":
                self.errMassage += "An error happens in residue 2:\n"
                self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
                self.errLevel    = 3
            tLinkedObj.combLigand["bonds"].append(aBond)

        # The added bond as input instruction file
        aBond = {}
        aBond["comp_id"]         = tLinkedObj.combLigand["name"]
        aBond["atom_id_1"]       = tLinkedObj.stdLigand1["atomName"] 
        aBond["atom_id_1_alias"] = self.getAtomAlias(aBond["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
        aBond["atom_id_2"]       = tLinkedObj.stdLigand2["atomName"] 
        aBond["atom_id_2_alias"] = self.getAtomAlias(aBond["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
        aBond["type"]            = tLinkedObj.suggestBondOrder
        aBond["value_dist"]      =  0.0             
        aBond["value_dist_esd"]  =  0.20             
        tLinkedObj.combLigand["bonds"].append(aBond)
        if aBond["atom_id_1_alias"] == "":
            self.errMassage += "An error happens in residue 1:\n"
            self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_1"] 
            self.errLevel    = 3
        if aBond["atom_id_2_alias"] == "":
            self.errMassage += "An error happens in residue 2:\n"
            self.errMassage += "The alias name for atom %s can not be found\n"%aBond["atom_id_2"] 
            self.errLevel    = 3

        # Chiral centers 
        for aChir in  tLinkedObj.stdLigand1["remainChirs"]:
            aChir["atom_id_centre_alias"] = self.getAtomAlias(aChir["atom_id_centre"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_1_alias"] = self.getAtomAlias(aChir["atom_id_1"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_2_alias"] = self.getAtomAlias(aChir["atom_id_2"], tLinkedObj.stdLigand1["remainAtoms"])
            aChir["atom_id_3_alias"] = self.getAtomAlias(aChir["atom_id_3"], tLinkedObj.stdLigand1["remainAtoms"])
            if aChir["atom_id_centre_alias"] and aChir["atom_id_1_alias"] and aChir["atom_id_2_alias"] and aChir["atom_id_3_alias"]:
                tLinkedObj.combLigand["chirs"].append(aChir)
            
        for aChir in  tLinkedObj.stdLigand2["remainChirs"]:
            aChir["atom_id_centre_alias"] = self.getAtomAlias(aChir["atom_id_centre"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_1_alias"] = self.getAtomAlias(aChir["atom_id_1"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_2_alias"] = self.getAtomAlias(aChir["atom_id_2"], tLinkedObj.stdLigand2["remainAtoms"])
            aChir["atom_id_3_alias"] = self.getAtomAlias(aChir["atom_id_3"], tLinkedObj.stdLigand2["remainAtoms"])
            if aChir["atom_id_centre_alias"] and aChir["atom_id_1_alias"] and aChir["atom_id_2_alias"] and aChir["atom_id_3_alias"]:
                tLinkedObj.combLigand["chirs"].append(aChir)
            
    def getAtomAlias(self, tAtomName, tAtoms):

        aReturnAlias = ""

        for aAtom in tAtoms:
            if aAtom["atom_id"]==tAtomName:
                if aAtom.has_key("atom_id_alias"):
                    aReturnAlias = aAtom["atom_id_alias"]
                else:
                   self.errMassage += "Atom %s does not have alias name\n"%aAtom["atom_id"] 
                   self.errLevel    = 3
                break
        return aReturnAlias

    def outTmpComboLigandMap(self, tLinkObj):

        # Output the mapping between 2 monomers and the combo-ligand for checking
 
        tmpFName = os.path.join(self.scrDir, "tmpComboLigandMap.txt")

        try:
            tmpF = open(tmpFName, "w")
        except IOError:
            print "%s can not be open for reading "%tmpFName
            self.errMessage +="%s can not be open for reading "%tmpFName
            self.errLevel = 3
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

    def comboLigToSimplifiedMmcif(self, tMonomer, tOutFName):

        try: 
            aOutF = open(tOutFName, "w")
        except IOError:
            print "%s can not be open for reading "%tOutFName
            self.errMessage +="%s can not be open for reading "%tOutFName
            self.errLevel = 3
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
         
            tComboId = "UNL"   
            numAts = len(tLinkedObj.combLigand["atoms"]) 
            numH   = len(tLinkedObj.stdLigand1["remainIdxHs"]) + len(tLinkedObj.stdLigand2["remainIdxHs"]) 
            aOutF.write("%s%s%s%s%6d%6d\n"%(tComboId.ljust(8), tComboId.ljust(8), "\'.           \'".just(20),\
                                           "non-polymer".ljust(20), numAts, numH))

            aOutF.write("data_comp_UNL\n")
            aOutF.write("#\n")
            if numAts != 0: 
                aOutF.write("loop_\n")
                aOutF.write("_chem_comp_atom.comp_id\n")
                aOutF.write("_chem_comp_atom.atom_id\n")
                aOutF.write("_chem_comp_atom.type_symbol\n")
                aOutF.write("_chem_comp_atom.type_energy\n")
                aOutF.write("_chem_comp_atom.partial_charge\n")
                aOutF.write("_chem_comp_atom.x\n")
                aOutF.write("_chem_comp_atom.y\n")
                aOutF.write("_chem_comp_atom.z\n")
                for aAtom in tLinkedObj.combLigand["atoms"]:
                    aOutF.write("%s%s%s%s%s%s%s%s\n"%(tCombId.ljust(10), aAtom["atom_id_alias"].ljust(8),\
                                aAtom["type_symbol"].ljust(6), aAtom["type_energy"].ljust(8), aAtom["partial_charge"].ljust(10),\
                                aAtom["x"].ljust(10), aAtom["y"].ljust(10), aAtom["z"].ljust(10)))
                
            else: 
                print "The combo-ligand has no atoms"
                self.errMessage += "The combo-ligand has no atoms\n"
                self.errLevel = 3
           
            if not self.errLevel:
                if len(tLinkedObj.combLigand["bonds"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_bond.comp_id\n")
                    aOutF.write("_chem_comp_bond.atom_id_1\n")
                    aOutF.write("_chem_comp_bond.atom_id_2\n")
                    aOutF.write("_chem_comp_bond.type\n")
                    aOutF.write("_chem_comp_bond.value_dist\n")
                    aOutF.write("_chem_comp_bond.value_dist_esd\n")
                    for aBond in tLinkedObj.combLigand["bonds"]:
                        aOutF.write("%s%s%s%s%s%s\n"%(tCombId.ljust(10), aBond["atom_id_1_alias"].ljust(8),\
                                     aBond["atom_id_2_alias"].ljust(8), aBond["type"].ljust(20), \
                                     aBond["value_dist"].ljust(10), aBond["value_dist_esd"].ljust(10)))
                 
                else: 
                    print "The combo-ligand has no bonds"
                    self.errMessage += "The combo-ligand has no bonds\n"
                    self.errLevel = 3

                if not self.errLevel and len(tLinkedObj.combLigand["chirs"]) !=0:
                    aOutF.write("loop_\n")
                    aOutF.write("_chem_comp_chir.comp_id\n")
                    aOutF.write("_chem_comp_chir.id\n")
                    aOutF.write("_chem_comp_chir.atom_id_centre\n")
                    aOutF.write("_chem_comp_chir.atom_id_1\n")
                    aOutF.write("_chem_comp_chir.atom_id_2\n")
                    aOutF.write("_chem_comp_chir.atom_id_3\n")
                    aOutF.write("_chem_comp_chir.volume_sign\n")
                    for aChi in tLinkedObj.combLigand["chirs"]:
                        aOutF.write("%s%s%s%s%s%s%s\n"%(tCombId.ljust(10), aChi["atom_id_centre_alias"], aChi["atom_id_1_alias"],\
                                                      aChi["atom_id_2_alias"], aChi["atom_id_3_alias"], aChi["volume_sign"]))
                aOutF.close()   
                        
    def buildComboLigand(self, tLinkedObj):

        # Change to mode 1 whether it is 0 or 1 initially 
 
        if not self.errLevel:

            tLinkedObj.combLigand["name"] = tLinkedObj.stdLigand1["name"].strip() + "-" + tLinkedObj.stdLigand2["name"].strip()
            print "The name of combo-ligand : %s "%tLinkedObj.combLigand["name"]

            self.adjustAtomsAndOthersForComboLigand(tLinkedObj)
       
            self.setInitComboLigand(tLinkedObj)

            self.outTmpComboLigandMap(tLinkedObj)    # Check 
    
            tLinkedObj.combLigand["inCif"] = os.path.join(self.scrDir, tLinkedObj.combLigand["name"] + "_comboIn.cif")

            self.comboLigToSimplifiedMmcif(tLinkedObj.combLigand, tLinkedObj.combLigand["inCif"])
  
            if self.errLevel: 
                self.setOneMonomer(tLinkedObj.combLigand)
                if self.errLevel:
                    tLinkedObj.outCombLigand["name"] = tLinkedObj.combLigand["name"]
                    tLinkedObj.outCombLigand["cifObj"] = Ccp4MmCifObj(tLinkedObj.combLigand["outCif"])            

    def getOneOrigAtomIdFromAlias(self, tAtoms, tAlias):
   
        # tAtoms come from self.combLigand["atoms"]
        # outCombLigand["cifObj"] contains ids which are actually alias_id
        aReturn = ""
        for aAtom in tAtoms:
            if aAtom["atom_id_alias"] == tAlias:
                aReturn = aAtom["atom_id"]
                break

        return aReturn
     
    def getChangesInModificationFromComboLigand(self, tLinkedObj):

        idSetL1 = []
        for aAtom in tLinkedObj.stdLigand1["remainAtoms"]:
            idSet1.append(aAtom["atom_id_alias"])

        idSetL2 = []
        for aAtom in tLinkedObj.stdLigand2["remainAtoms"]:
            idSet2.append(aAtom["atom_id_alias"])

        # Atoms 
        for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["atoms"]: 
            aAtom["atom_id_alias"] = aAtom["atom_id"]
            aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aAtom["atom_id_alias"])
            if aID != "":
                aAtom["atom_id"] = aID
            else:
                if not tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"].has_key("add_atoms"):
                    tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["add_atoms"] = []
                tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["add_atoms"].append(aAtom) 
                         
        # Bonds
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            aBond["atom_id_1_alias"] = aBond["atom_id_1"]
            aID1   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aBond["atom_id_1_alias"])
            if aID1 != "":
                aBond["atom_id_1"] = aID1
            aBond["atom_id_2_alias"] = aBond["atom_id_2"]
            aID2   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aBond["atom_id_2_alias"])
            if aID2 != "":
                aBond["atom_id_2"] = aID2
            if aBond["atom_id_1_alias"] in idSetL1 and aBond["atom_id_2_alias"] in idSetL1:
                self.checkBondMod(tLinkedObj.stdLigand1["remainBonds"], aBond, tLinkedObj.modLigand1["changed"]["bonds"]) 
            elif aBond["atom_id_1_alias"] in idSetL2 and aBond["atom_id_2_alias"] in idSetL2:
                self.checkBondMod(tLinkedObj.stdLigand2["remainBonds"], aBond, tLinkedObj.modLigand2["changed"]["bonds"]) 
            elif aBond["atom_id_1_alias"] in idSetL1 and aID2=="":
                tLinkedObj.modLigand1["added"]["bonds"].append(aBond)
            elif aBond["atom_id_2_alias"] in idSetL1 and aID1=="":
                tLinkedObj.modLigand1["added"]["bonds"].append(aBond)
            elif aBond["atom_id_1_alias"] in idSetL2 and aID2=="":
                tLinkedObj.modLigand2["added"]["bonds"].append(aBond)
            elif aBond["atom_id_2_alias"] in idSetL2 and aID1=="":
                tLinkedObj.modLigand2["added"]["bonds"].append(aBond)

            
        # Angles
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            aAng["atom_id_1_alias"] = aAng["atom_id_1"]
            aID1   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aAng["atom_id_1_alias"])
            if aID1 != "":
                aAng["atom_id_1"] = aID1
            aAng["atom_id_2_alias"] = aAng["atom_id_2"]
            aID2   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aAng["atom_id_2_alias"])
            if aID2 != "":
                aAng["atom_id_2"] = aID2
            aAng["atom_id_3_alias"] = aAng["atom_id_3"]
            aID3   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aAng["atom_id_3_alias"])
            if aID3 != "":
                aAng["atom_id_3"] = aID3
            if aAng["atom_id_1_alias"] in idSetL1 and aAng["atom_id_2_alias"] in idSetL1\
               and aAng["atom_id_3_alias"] in idSetL1:
                self.checkAngMod(tLinkedObj.stdLigand1["remainAngs"], aAng, tLinkedObj.modLigand1["changed"]["angles"])  
            elif aAng["atom_id_1_alias"] in idSetL2 and aAng["atom_id_2_alias"] in idSetL2\
               and aAng["atom_id_3_alias"] in idSetL2:
                self.checkAngMod(tLinkedObj.stdLigand2["remainAngs"], aAng, tLinkedObj.modLigand2["changed"]["angles"])  
            elif  aAng["atom_id_1_alias"] in idSetL1 and aAng["atom_id_2_alias"] in idSetL1\
               and AD3 =="":
                LinkedObj.modLigand1["added"]["angles"].append(aAng)
            elif  aAng["atom_id_3_alias"] in idSetL1 and aAng["atom_id_2_alias"] in idSetL1\
               and AD1 =="":
                LinkedObj.modLigand1["added"]["angles"].append(aAng)
            elif  aAng["atom_id_1_alias"] in idSetL2 and aAng["atom_id_2_alias"] in idSetL2\
               and AD3 =="":
                LinkedObj.modLigand2["added"]["angles"].append(aAng)
            elif  aAng["atom_id_3_alias"] in idSetL2 and aAng["atom_id_2_alias"] in idSetL2\
               and AD1 =="":
                LinkedObj.modLigand2["added"]["angles"].append(aAng)
                 
        # Torsions
        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
            aTor["atom_id_1_alias"] = aTor["atom_id_1"]
            aID1   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aTor["atom_id_1_alias"])
            if aID1 != "":
                aTor["atom_id_1"] = aID1
            aTor["atom_id_2_alias"] = aTor["atom_id_2"]
            aID2   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aTor["atom_id_2_alias"])
            if aID2 != "":
                aTor["atom_id_2"] = aID2
            aTor["atom_id_3_alias"] = aTor["atom_id_3"]
            aID3   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aTor["atom_id_3_alias"])
            if aID3 != "":
                aTor["atom_id_3"]   = aID3
            aTor["atom_id_4_alias"] = aTor["atom_id_4"]
            aID4   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aTor["atom_id_4_alias"])
            if aID4 != "":
                aTor["atom_id_4"]   = aID4
            if aTor["atom_id_2"] in idSetL1 and aTor["atom_id_3"] in idSetL1:
                if aTor["atom_id_1"] in idSetL1 and aTor["atom_id_4"] in idSetL1:
                    self.checkTorMod(tLinkedObj.stdLigand1["remainTors"], aTor, tLinkedObj.modLigand1["changed"]["tors"])
                elif (aTor["atom_id_1"] in idSetL1 and aID4=="") or (aTor["atom_id_4"] in idSetL1 and aID1==""):
                    LinkedObj.modLigand1["added"]["tors"].append(aTor)
            elif aTor["atom_id_2"] in idSetL2 and aTor["atom_id_3"] in idSetL2:
                if aTor["atom_id_1"] in idSetL2 and aTor["atom_id_4"] in idSetL2:
                    self.checkTorMod(tLinkedObj.stdLigand2["remainTors"], aTor, tLinkedObj.modLigand2["changed"]["tors"])
                elif (aTor["atom_id_1"] in idSetL2 and aID4=="") or (aTor["atom_id_4"] in idSetL2 and aID1==""):
                    LinkedObj.modLigand2["added"]["tors"].append(aTor)
            
        # Chirs
        for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
            aChi["atom_id_centre_alias"] = aChi["atom_id_centre"]
            aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aChi["atom_id_4_alias"])
            if aID != "":
                aChi["atom_id_centre"]   = aID
            aChi["atom_id_1_alias"] = aChi["atom_id_1"]
            aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aChi["atom_id_1_alias"])
            if aID != "":
                aChi["atom_id_1"] = aID
            aChi["atom_id_2_alias"] = aChi["atom_id_2"]
            aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aChi["atom_id_2_alias"])
            if aID != "":
                aChi["atom_id_2"] = aID
            aChi["atom_id_3_alias"] = aChi["atom_id_3"]
            aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aChi["atom_id_3_alias"])
            if aID != "":
                aChi["atom_id_3"]   = aID
            if aChi["atom_id_centre"] in idSetL1:
                if aChi["atom_id_1"] in idSetL1 and aChi["atom_id_2"] in idSetL1 and aChi["atom_id_3"] in idSetL1:
                    self.checkChiMod(tLinkedObj.stdLigand1["remainChirs"], aChi, tLinkedObj.modLigand1["changed"]["chirs"])
            if aChi["atom_id_centre"] in idSetL2:
                if aChi["atom_id_1"] in idSetL2 and aChi["atom_id_2"] in idSetL2 and aChi["atom_id_3"] in idSetL2:
                    self.checkChiMod(tLinkedObj.stdLigand2["remainChirs"], aChi, tLinkedObj.modLigand2["changed"]["chirs"])

        # Planes
        for aPl in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys():
            for aAtom in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                aAtom["atom_id_alias"] = aAtom["atom_id"]
                aID   = getOneOrigAtomIdFromAlias(tLinkedObj.combLigand["atoms"], aAtom["atom_id_alias"])
                if aID != "":
                    aAtom["atom_id"] = aID
   
    def checkBondMod(self, tOriBonds, tBond, tModBonds):

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
        if float(tOriBond["value_dist"]) != float(tBond["type"]):
            lChanged = True
        return lChanged
   
    def checkAngMod(self, tOriAngs, tAng, tModAngs):

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

    def compare2Angs(self, tOriAng, tAng):
   
        lChanged = False
        if tOriAng["type"].upper()[:3] !=tAng["type"].upper()[:3]:
            lChanged = True
        if float(tOriAng["value_dist"]) != float(tAng["type"]):
            lChanged = True
        return lChanged
   
    def checkTorMod(self, tOriTors, tTor, tModTors):

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

    def compare2Tors(self, tOriTor, tTor):
   
        lChanged = False
        if float(tOriTor["value_angle"]) !=float(tTor["value_angle"]):
            lChanged = True
        return lChanged
   

    def extractOneLinkInfo(self, tLinkedObj):
 
        self.getLinkInfo(tLinkedObj)
        self.getChangesInModification(tLinkedObj) 

    def getLinkInfo(self, tLinkedObj):

        atm1 = tLinkedObj.stdLigand1["atomName_alias"]
        atm2 = tLinkedObj.stdLigand2["atomName_alias"]
        tLinkedObj.clink["bonds"] = []
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            if (aBond["atom_id_1"] == atm1 and aBond["atom_id_2"] == atm2)\
                or (aBond["atom_id_1"] == atm2 and aBond["atom_id_2"] == atm1):
                tLinkedObj.clink["bonds"].append(aBond)

        tLinkedObj.clink["angles"] = []
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            if aAng["atom_id_2"]==atm1:
                if aAng["atom_id_1"]==atm2 or aAng["atom_id_3"]==atm2:
                    tLinkedObj.clink["angles"].append(aAng) 
            elif aAng["atom_id_2"]==atm2:
                if aAng["atom_id_1"]==atm1 or aAng["atom_id_3"]==atm1:
                    tLinkedObj.clink["angles"].append(aAng) 
   
        tLinkedObj.clink["torsions"] = []
        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
            if (aTor["atom_id_2"] == atm1 and aTor["atom_id_3"] == atm2) or\
               (aTor["atom_id_2"] == atm2 and aTor["atom_id_3"] == atm1):
                tLinkedObj.clink["torsions"].append(aTor)

        tLinkedObj.clink["chirals"] = [] 
        for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
            if aChi["atom_id_centre"] ==atm1:
                tLinkedObj.clink["chirals"].append(aChi)
            elif aChi["atom_id_centre"] ==atm2:
                tLinkedObj.clink["chirals"].append(aChi)

        if tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"].has_key("planes"):
            if len(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys()) !=0:
                tLinkedObj.clink["planes"] = []
                for aPl in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys():   
                    for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                        if aPlAtm[atom_id] == atm1 or aPlAtm[atom_id]==atm2:
                            tLinkedObj.clink["planes"].append(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl])
                            break
    
       

    def outOneLinkInfo(self, tLinkedObj):

        aOutCifName = tLinkedObj.combLigand["name"] + "_link.cif"
        
        try: 
            aOutCif = open(aOutCifName, "w")
        except IOError:
            self.errMessage.append("%s can not be open for writting "%aOutCifName)
            self.errLevel  = 6
        
        else:
            self.outVerInfo(aOutCif)
            self.outCompList(aOutCif, tLinkedObj)
            self.outModiList(aOutCif, tLinkedObj)
            self.outLinkList(aOutCif, tLinkedObj)
            self.outAllComps(aOutCif,  tLinkedObj) 
            self.outAllModis(aOutCif,  tLinkedObj) 
            self.outAllLinks(aOutCif, tLinkedObj) 
      
            aOutCif.close()

    def outVerInfo(self, tOutFile):

        pass

    def outCompList(self, tOutFile, tLinkedObj):

        if tLinkedObj.stdLigand1["OutComp"] or tLinkedObj.stdLigand2["OutComp"] :
            tOutFile.write("data_comp_list\n\n")
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_comp.id\n")
            tOutFile.write("_chem_comp.three_letter_code\n")
            tOutFile.write("_chem_comp.name\n")
            tOutFile.write("_chem_comp.group\n")
            tOutFile.write("_chem_comp.number_atoms_all\n")
            tOutFile.write("_chem_comp.number_atoms_nh\n")
            

    def extractOneLinkInfo(self, tLinkedObj):
 
        self.getLinkInfo(tLinkedObj)
        self.getChangesInModificationFromCombLigand(tLinkedObj) 

    def getLinkInfo(self, tLinkedObj):

        atm1 = tLinkedObj.stdLigand1["atomName_alias"]
        atm2 = tLinkedObj.stdLigand2["atomName_alias"]
        tLinkedObj.clink["bonds"] = []
        for aBond in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["bonds"]:
            if (aBond["atom_id_1"] == atm1 and aBond["atom_id_2"] == atm2)\
                or (aBond["atom_id_1"] == atm2 and aBond["atom_id_2"] == atm1):
                tLinkedObj.clink["bonds"].append(aBond)

        tLinkedObj.clink["angles"] = []
        for aAng in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["angles"]:
            if aAng["atom_id_2"]==atm1:
                if aAng["atom_id_1"]==atm2 or aAng["atom_id_3"]==atm2:
                    tLinkedObj.clink["angles"].append(aAng) 
            elif aAng["atom_id_2"]==atm2:
                if aAng["atom_id_1"]==atm1 or aAng["atom_id_3"]==atm1:
                    tLinkedObj.clink["angles"].append(aAng) 
   
        tLinkedObj.clink["torsions"] = []
        for aTor in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["tors"]:
            if (aTor["atom_id_2"] == atm1 and aTor["atom_id_3"] == atm2) or\
               (aTor["atom_id_2"] == atm2 and aTor["atom_id_3"] == atm1):
                tLinkedObj.clink["torsions"].append(aTor)

        tLinkedObj.clink["chirals"] = [] 
        for aChi in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["chirs"]:
            if aChi["atom_id_centre"] ==atm1:
                tLinkedObj.clink["chirals"].append(aChi)
            elif aChi["atom_id_centre"] ==atm2:
                tLinkedObj.clink["chirals"].append(aChi)

        if tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"].has_key("planes"):
            if len(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys()) !=0:
                tLinkedObj.clink["planes"] = []
                for aPl in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"].keys():   
                    for aPlAtm in tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl]:
                        if aPlAtm[atom_id] == atm1 or aPlAtm[atom_id]==atm2:
                            tLinkedObj.clink["planes"].append(tLinkedObj.outCombLigand["cifObj"]["comps"]["UNL"]["planes"][aPl])
                            break
    
       
    def outOneLinkInfo(self, tLinkedObj):

        aOutCifName = tLinkedObj.combLigand["name"] + "_link.cif"
        
        try: 
            aOutCif = open(aOutCifName, "w")
        except IOError:
            self.errMessage.append("%s can not be open for writting "%aOutCifName)
            self.errLevel  = 6
        
        else:
            self.outVerInfo(aOutCif)
            self.outCompList(aOutCif, tLinkedObj)
            self.outModiList(aOutCif, tLinkedObj)
            self.outLinkList(aOutCif, tLinkedObj)
            self.outAllComps(aOutCif,  tLinkedObj) 
            self.outAllModis(aOutCif,  tLinkedObj) 
            self.outAllLinks(aOutCif, tLinkedObj) 
      
            aOutCif.close()

    def outVerInfo(self, tOutFile):

        pass

    def outCompList(self, tOutFile, tLinkedObj):

        if tLinkedObj.stdLigand1["OutComp"] or tLinkedObj.stdLigand2["OutComp"] :
            tOutFile.write("data_comp_list\n\n")
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_comp.id\n")
            tOutFile.write("_chem_comp.three_letter_code\n")
            tOutFile.write("_chem_comp.name\n")
            tOutFile.write("_chem_comp.group\n")
            tOutFile.write("_chem_comp.number_atoms_all\n")
            tOutFile.write("_chem_comp.number_atoms_nh\n")
            
            if tLinkedObj.stdLigand1["OutComp"]:
                aNL = tLinkedObj.stdLigand1["Name"]
                if len(aNL) >= 3:
                    aN3 = tLinkedObj.stdLigand1["Name"][:3]
                else:
                    aN3 = aNL
                aLD = ""
                if tLinkedObj.stdLigand1["LongDis"].find("\"") != -1:
                    aLD = "\'" + tLinkedObj.stdLigand1["LongDis"] + "\'"
                elif tLinkedObj.stdLigand1["LongDis"].find("\'") != -1:
                    aLD = "\"" + tLinkedObj.stdLigand1["LongDis"] + "\""
                else:
                    aLD = "\"" + tLinkedObj.stdLigand1["LongDis"] + "\""
                aG   = tLinkedObj.stdLigand1["Group"]
                aAN  = len(tLinkedObj.stdLigand1["Atoms"])
                aAHN = len(tLinkedObj.stdLigand1["HAtoms"])
                aL="%s%s%s%s%20d%20d\n"%(aNL.ljust(10), aN3.ljust(10), aLD.ljust(len(aLD)+3), aG.ljust(10), aAN, aAHN)
                tOutFile.write(aL) 
                 
            if tLinkedObj.stdLigand2["OutComp"]:
                aNL = tLinkedObj.stdLigand2["Name"]
                if len(aNL) >= 3:
                    aN3 = tLinkedObj.stdLigand2["Name"][:3]
                else:
                    aN3 = aNL
                aLD = ""
                if tLinkedObj.stdLigand2["LongDis"].find("\"") != -1:
                    aLD = "\'" + tLinkedObj.stdLigand2["LongDis"] + "\'"
                elif tLinkedObj.stdLigand2["LongDis"].find("\'") != -1:
                    aLD = "\"" + tLinkedObj.stdLigand2["LongDis"] + "\""
                else:
                    aLD = "\"" + tLinkedObj.stdLigand2["LongDis"] + "\""
                aG   = tLinkedObj.stdLigand2["Group"]

                aAN  = len(tLinkedObj.stdLigand2["Atoms"])
                aAHN = len(tLinkedObj.stdLigand2["HAtoms"])
                aL="%s%s%s%s%20d%20d\n"%(aNL.ljust(10), aN3.ljust(10), aLD.ljust(len(aLD)+3), aG.ljust(10), aAN, aAHN)
                tOutFile.write(aL) 

            tOutFile.write("\n")
                 
    def outModiList(self, tOutFile, tLinkedObj):

        if tLinkedObj.stdLigand1["OutModi"] or tLinkedObj.stdLigand2["OutModi"] :
            tOutFile.write("data_mod_list\n\n")
            tOutFile.write("loop_\n")
            tOutFile.write("_chem_mod.id\n")
            tOutFile.write("_chem_mod.name\n")
            tOutFile.write("_chem_mod.comp_id\n")
            tOutFile.write("_chem_mod.group_id\n")
            
            if tLinkedObj.stdLigand1["OutModi"]:
                aMN = tLinkedObj.stdLigand1["ModName"]
                aLD = tLinkedObj.stdLigand1["LongDis"]
                aNL = tLinkedObj.stdLigand1["Name"]
                aG   = tLinkedObj.stdLigand1["Group"]
                aL="%s%s%s%s\n"%(aMN.ljust(10), aLD.ljust(len(aLD)+3), aNL.ljust(10), aG.ljust(10))
                tOutFile.write(aL) 
                
            if tLinkedObj.stdLigand2["OutModi"]:
                aMN = tLinkedObj.stdLigand2["ModName"]
                aLD = tLinkedObj.stdLigand2["LongDis"]
                aNL = tLinkedObj.stdLigand2["Name"]
                aG   = tLinkedObj.stdLigand2["Group"]
                aL="%s%s%s%s\n"%(aMN.ljust(10), aLD.ljust(len(aLD)+3), aNL.ljust(10), aG.ljust(10))
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
        aLN1  = tLinkedObj.stdLigand1["Name"]
        aLMN1 = tLinkedObj.stdLigand1["ModName"]
        aG1   = tLinkedObj.stdLigand1["Group"]
        aLN2  = tLinkedObj.stdLigand2["Name"]
        aLMN2 = tLinkedObj.stdLigand2["ModName"]
        aG2   = tLinkedObj.stdLigand2["Group"]
        aL="%s%s%s%s%s%s%s%s\n"%(aLID.ljust(15), aLN1.ljust(10), aLMN1.ljust(12), aG1.ljust(8),\
                                 aLN2.ljust(10), aLMN2.ljust(12), aG2.ljust(8), aLID.ljust(15))
        tOutFile.write(aL) 
         
    def outAllComps(self, tOutFile, tLinkedObj):

        if tLinkedObj.stdLigand1["OutComp"]:
            self.outOneComp(tOutFile, tLinkedObj.stdLigand1)
        if tLinkedObj.stdLigand2["OutComp"]:
            self.outOneComp(tOutFile, tLinkedObj.stdLigand2)

    def outOneComp(self, tOutFile, tMonomer):

        tOutFile.write("data_comp_%s\n\n"%tMonomer["Name"])
        tOutFile.write("loop_\n")
    
        # Atom section 
        
        
