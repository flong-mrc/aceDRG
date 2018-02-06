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

from utility  import listComp
from utility  import listComp2
from utility  import listCompDes
from utility  import listCompAcd
from utility  import setBoolDict
from utility  import splitLineSpa
from utility  import splitLineSpa2

class FileTransformer :

    def __init__(self):

        self.dataDescriptor  = {} 
        self.strDescriptors  = {} 
        self.strDescriptors["defProps"]    = ["loop_", "_pdbx_chem_comp_descriptor.comp_id", "_pdbx_chem_comp_descriptor.type", \
                               "_pdbx_chem_comp_descriptor.program", "_pdbx_chem_comp_descriptor.program_version",\
                               "_pdbx_chem_comp_descriptor.descriptor"]
        self.hasStrDescriptors = False
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

        self.mmCifHasCoords   = False

        self.ccp4MmCifDataMap = {}

        self.nameMapingCifMol = {}             # e.g. self.nameMapingCifMol[1]   = name 
                                               # where 1 : atom serial number   name : atom name
        self.nameMapingCifMol = {}            
        self.nameMapingCifMol["nonH"] = {}            
        self.nameMapingCifMol["H"]    = {}  

        self.nameMapingMol2 = {}
        self.nameMapingMol2["nonH"] = {}
        self.nameMapingMol2["H"]    = {}          
                                             

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
       
            if len(self.allBlockLs):
                for aBlk in self.allBlockLs:
                    self.parseOneMmCifBlk(aBlk)

            if len(all2ColLines) !=0:
                self.parserAll2Cols(all2ColLines)

            self.TmpChemCheck()
            
            self.selectAtomCoordinates()
 
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
                    self.hasStrDescriptors = True
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
                        # print "aAtom['_chem_comp_atom.charge'] ", aAtom["_chem_comp_atom.charge"] 
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

        #if len(colIdx):
        #    if tProp =="atom":
        #        for aAtom in self.atoms:
        #            for i in range(len(colIdx)):
        #                if colIdx[i].find("chem_comp_atom.atom_id") != -1:
        #                    print "Prop %s is %s "%(colIdx[i], aAtom[colIdx[i]].ljust(8))

        # Check
        if len(self.inputCharge) !=0:
            print "The following atoms have charges "
            for aName in self.inputCharge.keys():
                print "Name : ", aName, " charge : ", self.inputCharge[aName]

        if len(colIdx):
            #for i in range(len(colIdx)):
            #    print colIdx[i]
            if tProp =="atom":
                for aAtom in self.atoms:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aAtom[colIdx[i]].ljust(8)
                    #print aStr
            elif tProp =="bond":
                for aBond in self.bonds:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aBond[colIdx[i]].ljust(8)
                    #print aStr
            elif tProp =="chiral":
                for aChiral in self.chirals:
                    aStr = ""
                    for i in range(len(colIdx)):
                        aStr +="%s"%aChiral[colIdx[i]].ljust(8)
                    #print aStr

        #if tProp =="strDescriptor":
        #    if len(self.strDescriptors["props"]):
        #        for aProp in self.strDescriptors["props"]:
        #            print aProp
        #        for aEn in self.strDescriptors["entries"]:
        #            print aEn.strip()

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
                         
    def selectAtomCoordinates(self):

        # Select coordinates as .x, .y and .z when there are several sets of
        # coordinates for atoms 
        lq1 = False
        lq2 = False
        lq3 = False

        self.mmCifHasCoords = False  
        if len(self.atoms) > 0:
            if self.atoms[0].has_key("_chem_comp_atom.model_Cartn_x") and\
               self.atoms[0].has_key("_chem_comp_atom.model_Cartn_y") and\
               self.atoms[0].has_key("_chem_comp_atom.model_Cartn_z"):
                for aAtom in self.atoms:
                    if aAtom["_chem_comp_atom.model_Cartn_x"].find("?") != -1\
                       or aAtom["_chem_comp_atom.model_Cartn_y"].find("?") != -1\
                       or aAtom["_chem_comp_atom.model_Cartn_z"].find("?") != -1:
                        lq1 = True
                        break
                if not lq1 : 
                    print "Using _chem_comp_atom.pdbx_model_Cartn_x"
                    self.mmCifHasCoords   = True
                    for aAtom in self.atoms:
                        aAtom["_chem_comp_atom.x"] = aAtom["_chem_comp_atom.model_Cartn_x"]
                        aAtom["_chem_comp_atom.y"] = aAtom["_chem_comp_atom.model_Cartn_y"]
                        aAtom["_chem_comp_atom.z"] = aAtom["_chem_comp_atom.model_Cartn_z"]
                       
            if not self.mmCifHasCoords :
                print "Not _chem_comp_atom.pdbx_model_Cartn_x"
                if self.atoms[0].has_key("_chem_comp_atom.pdbx_model_Cartn_x_ideal")\
                   and self.atoms[0].has_key("_chem_comp_atom.pdbx_model_Cartn_y_ideal")\
                   and self.atoms[0].has_key("_chem_comp_atom.pdbx_model_Cartn_z_ideal"):
                    for aAtom in self.atoms:           
                        if aAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"].find("?") !=-1\
                           or aAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"].find("?") !=-1\
                           or aAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"].find("?") !=-1:
                            lq2 = True
                            break
                    if not lq2:
                        self.mmCifHasCoords   = True
                        print "Using _chem_comp_atom.pdbx_model_Cartn_x_ideal"
                        for aAtom in self.atoms:           
                            aAtom["_chem_comp_atom.x"] = aAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"]
                            aAtom["_chem_comp_atom.y"] = aAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"]
                            aAtom["_chem_comp_atom.z"] = aAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"]

            if not self.mmCifHasCoords :
                print "not _chem_comp_atom.model_Cartn_x_ideal"
                if self.atoms[0].has_key("_chem_comp_atom.x") and\
                   self.atoms[0].has_key("_chem_comp_atom.y") and\
                   self.atoms[0].has_key("_chem_comp_atom.z"):
                    for aAtom in self.atoms:
                        if aAtom["_chem_comp_atom.x"].find("?") !=-1\
                           or aAtom["_chem_comp_atom.y"].find("?") !=-1\
                           or aAtom["_chem_comp_atom.z"].find("?") !=-1:
                            lq3 = True
                            break
                    if not lq3:
                        print "Using _chem_comp_atom.x,y,z"
                        self.mmCifHasCoords   = True
                    else:
                        print "not _chem_comp_atom.x"
                        print "No original coordinates will be used for conformer generations"
                        self.mmCifHasCoords = False  

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
            # print "Befor re-arrange, number of atoms is : ", len(self.atoms) 
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
                print "NameMap ", nAtm, " : ", self.nameMapingCifMol["nonH"][nAtm]
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
                if self.mmCifHasCoords :
                    x = aAtom["_chem_comp_atom.x"]
                    y = aAtom["_chem_comp_atom.y"]
                    z = aAtom["_chem_comp_atom.z"]
                else:
                    x = "0.0000"
                    y = "0.0000"
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

    def mol2Reader(self, tFileName):
        """Read a file of mol2 for name matching purpose"""

        self.atoms = []
        self.bonds = []
        try:
            aFile = open(tFileName, "r")
        except IOError:
            print "%s can not be open for reading "%tFileName
            sys.exit()
        else:
            allLs = aFile.readlines()
            aFile.close()

        lAtom = False
        lBond = False
  
        idxNonHA = 0
        idxHA    = 0
        idxBo    = 0
       
        tNonHAtmOldNewMap = {} 
        tHAtoms = []
        for aL in allLs:
            aL = aL.strip()
            if len(aL) > 0 and aL[0].find("#") == -1:
                if aL.find("@<TRIPOS>ATOM") !=-1:
                    lAtom = True
                    lBond = False
                elif aL.find("@<TRIPOS>BOND") !=-1:
                    lAtom = False
                    lBond = True
                elif lAtom or lBond :
                    strs = aL.split()
                    if lAtom:
                        if len(strs) >= 6:
                            aAtom = {}
                            aAtom["mol2Seri"] = strs[0].strip()
                            aAtom["name"]     = strs[1].strip()
                            elem = ""
                            if strs[5].find(".")==-1:
                                elem = strs[5].strip()
                            else:
                                elem = strs[5].strip().split(".")[0].strip()
                            aAtom["element"] = elem
                            if elem.upper() != "H":
                                aAtom["seriNum"]  = idxNonHA
                                self.nameMapingMol2["nonH"][aAtom["seriNum"]] = aAtom["name"]
                                tNonHAtmOldNewMap[aAtom["mol2Seri"]]  = aAtom["seriNum"] 
                                idxNonHA +=1
                                self.atoms.append(aAtom)
                            else:
                                tHAtoms.append(aAtom)
            
                    elif lBond:
                        if len(strs) >= 4:
                            aBond = {}
                            aBond["mol2Seri"] = strs[0].strip()     
                            aBond["1stAtm"]   = strs[1].strip()
                            aBond["2ndAtm"]   = strs[2].strip()
                            aBond["mol2Ord"]  = strs[3].strip()
                            self.bonds.append(aBond)
      
        if len(tHAtoms) > 0:
           for aAtm in tHAtoms:
               curN = len(self.atoms)
               tNonHAtmOldNewMap[aAtm["mol2Seri"]] = curN
               self.atoms.append(aAtm)

        self.setNameMap(tNonHAtmOldNewMap)
        
    def setNameMap(self, tNonHAtmOldNewMap):

        for aBond in self.bonds:                    
             idxH     = -1
             idxOther = -1
             if tNonHAtmOldNewMap.has_key(aBond["1stAtm"])\
                 and tNonHAtmOldNewMap.has_key(aBond["2ndAtm"]):
                  idx1 = tNonHAtmOldNewMap[aBond["1stAtm"]]
                  elem1 = self.atoms[idx1]["element"]
                  idx2 = tNonHAtmOldNewMap[aBond["2ndAtm"]]
                  elem2 = self.atoms[idx2]["element"]
                  if elem1=="H" and elem2 != "H":
                      idxH      = idx1
                      idxOther  = idx2
                  elif elem2=="H" and elem1 != "H":
                      idxH      = idx2
                      idxOther  = idx1
                  if idxH !=-1 and idxOther !=-1:
                      otherId  = self.atoms[idxOther]["name"]
                      if not self.nameMapingMol2["H"].has_key(otherId):
                          self.nameMapingMol2["H"][otherId] = []
                      self.nameMapingMol2["H"][otherId].append(self.atoms[idxH]["name"])
        print "The follow are information on H connections "
        for aK in sorted(self.nameMapingMol2["H"].iterkeys()):
            for aH in self.nameMapingMol2["H"][aK]:
                print "%s  connects to %s "%(aK, aH)
        
class Ccp4MmCifObj (dict) :

    def __init__(self, tInFileName):

        self["errMessage"] = ""
        self["errLevel"]   = 0

        self["ccp4CifBlocks"] = {}
        self["ccp4CifBlocks"]["head_sec"] = []

        self["ccp4CifObj"]               = {}
        self["ccp4CifObj"]["lists"]      = {}
        self["ccp4CifObj"]["comps"]      = {}
        self["ccp4CifObj"]["mods"]       = {}
        self["ccp4CifObj"]["links"]      = {}
        self["ccp4CifObj"]["instructs"]  = {}

        self["inCif"] = tInFileName
        print self["inCif"]
        if os.path.isfile(self["inCif"]):
            self.getCCP4CifObj()
        else:
            self["errMessage"] = "%s does not exist "%self["inCif"]
            self["errLevel"]   = 1
            print self["errMessage"]
            
       
    def getCCP4CifObj(self):
       
        try : 
            aCCP4Cif = open(self["inCif"], "r")
        except IOError:
            self["errMessage"] =  "% can not be open for reading ! "%self["inCif"]
            self["errLevel"]   = 1
        else:
            self["ccp4CifBlocks"] = {}
            self["ccp4CifBlocks"]["head_sec"] = []
            allLs = aCCP4Cif.readlines()
            lCifBegin = False
            aDataHead = ""
            for aL in allLs:
                if len(aL.strip()) !=0:
                    if aL.find("data_") !=-1:
                        aDataHead = aL.strip()
                        # print aDataHead
                        self["ccp4CifBlocks"][aDataHead] = []
                        lCifBegin = True
                    elif not lCifBegin:
                        self["ccp4CifBlocks"]["head_sec"].append(aL) 
                    elif aDataHead !="":
                        self["ccp4CifBlocks"][aDataHead].append(aL)
                    else:
                        print "Bug in analyse the line : "
                        print aL
                        sys.exit()

            for aKey in self["ccp4CifBlocks"].keys():        
                if aKey != "head_sec":
                    self.setupOneBlock(aKey, self["ccp4CifBlocks"])
            """
            if len(self["ccp4CifObj"]["comps"].keys()):
                print "Number of comps ", len(self["ccp4CifObj"]["comps"].keys())
                for aKey in sorted(self["ccp4CifObj"]["comps"].iterkeys()):
                    self.printOneComp(aKey)
            """

 
    def setupOneBlock(self, tKey, tBlock):
        
        print "Block ", tKey
        allLoops = []
        aLoop = []
        for aL in tBlock[tKey]:
            if aL.find("loop_") !=-1:
                if len(aLoop) != 0:
                    allLoops.append(aLoop)
                aLoop = []
            else:
                aL = aL.strip()
                if len(aL):
                    aLoop.append(aL)
 
        if len(aLoop) != 0:
            allLoops.append(aLoop)
        print "It contains ", len(allLoops), " loops "

        lB = {}
        lB["LC"]  =False
        lB["LM"]  =False
        lB["LL"]  =False
        lB["INS"] =False
        lB["COMP"]=False
        lB["MOD"] =False
        lB["LINK"]=False

        strsKey = tKey.strip().split("_")
           
        if len(strsKey) == 3:
            strsKey[2] = strsKey[2].strip() 
            if strsKey[2].find("list") !=-1:
                if strsKey[1].find("comp")   !=-1: 
                    setBoolDict("LC", True, lB, False)                
                elif strsKey[1].find("mod")  !=-1:                 
                    setBoolDict("LM", True, lB, False)                
                elif strsKey[1].find("link") !=-1:                 
                    setBoolDict("LL", True, lB, False)                
            else:
                if strsKey[1].find("comp") !=-1:
                    setBoolDict("COMP", True, lB, False)                
                if strsKey[1].find("mod") !=-1:
                    setBoolDict("MOD", True, lB, False)                
                if strsKey[1].find("link") !=-1:
                    setBoolDict("LINK", True, lB, False)                
                if strsKey[1].find("instruct") !=-1:
                    setBoolDict("INS", True, lB, False)                
                    
            for aLoop in allLoops:
                if lB["LC"]:
                    #print "set comp list loop "
                    self["ccp4CifObj"]["lists"]["comp"] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["lists"]["comp"], aLoop, "list")
                elif lB["LM"]:
                    #print "set mod list loop "
                    self["ccp4CifObj"]["lists"]["mod"] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["lists"]["mod"],  aLoop, "list")
                elif lB["LL"]:
                    #print "set link list loop "
                    self["ccp4CifObj"]["lists"]["link"] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["lists"]["link"], aLoop, "list")
                elif lB["COMP"]:
                    #print "set  comp %s loop "%strsKey[2]
                    if not self["ccp4CifObj"]["comps"].has_key(strsKey[2]):
                        self["ccp4CifObj"]["comps"][strsKey[2]] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["comps"][strsKey[2]],aLoop)
                elif lB["MOD"]:
                    #print "set mod %s loop "%strsKey[2]
                    if not self["ccp4CifObj"]["mods"].has_key(strsKey[2]):
                        self["ccp4CifObj"]["mods"][strsKey[2]] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["mods"][strsKey[2]], aLoop)
                elif lB["LINK"]:
                    #print "set list %s loop "%strsKey[2]
                    if not self["ccp4CifObj"]["links"].has_key(strsKey[2]):
                        self["ccp4CifObj"]["links"][strsKey[2]] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["links"][strsKey[2]], aLoop)
                elif lB["INS"]:
                    if not self["ccp4CifObj"]["instructs"].has_key(strsKey[2]):
                        self["ccp4CifObj"]["instructs"][strsKey[2]] = {}
                    self.setupOneLoop(self["ccp4CifObj"]["instructs"][strsKey[2]], aLoop)
        
            # Re-group plane-atoms in comps according to planes:
            if self["ccp4CifObj"].has_key("comps"):
                for aComp in self["ccp4CifObj"]["comps"]:
                    self.setPlaneAtmGroupInOneComp(aComp)
 
    def setupOneLoop(self, tDict, tLoop, tKW=""):

        aKeySet   = []
        aValueSet = []
        for aL in tLoop:
            aL = aL.strip()
            if len(aL) !=0:
                aStrSet = [] 
                if tKW =="list":
                    aStrSet = splitLineSpa2(aL)
                    print aL
                    print aStrSet
                else:
                    aStrSet = splitLineSpa(aL)
                for aStr in aStrSet:
                    if aStr.strip()[0] =="_":
                        aKeySet.append(aStr)
                        if tKW =="list":
                            print aKeySet
                    else:
                        aValueSet.append(aStr)

        idKeys =[]
        for aKey in aKeySet:
            k = aKey.strip().split(".")[-1]
            idKeys.append(k)

        nK = len(aKeySet)
        nV = len(aValueSet) 
        print "Number of keys in the loop is ", nK
        print "Number of values in the loop is ", nV
        if nK > 0 and nV > 0:
            if tKW =="list":
                idxV = 0
                remV = nV
                while remV >= nK:
                    idx = idxV
                    compId = aValueSet[idx]
                    tDict[compId] = {}
                    for i in range(nK):
                        k = idKeys[i]
                        v = aValueSet[idxV].strip()
                        tDict[compId][k] = v
                        idxV +=1
                    remV -=nK
                    print "For list entry  %s, it has the following properties "%compId
                    for aPro in tDict[compId].keys():
                        print "%s : %s "%(aPro, tDict[compId][aPro])
            else:
                aType    = ""
                parts = aKeySet[0].strip().split(".")[0].strip().split("_")
                if len(parts) ==5:
                    aType = parts[3] + "_" + parts[4] + "s"
                else:
                    aType = parts[3] + "s"
                if not tDict.has_key(aType):
                    tDict[aType] = []

                idxV = 0
                remV = nV
                print aType
                print idKeys
                while remV >= nK:
                    aObj = {}
                    for i in range(nK):
                        k = idKeys[i]
                        v = aValueSet[idxV].strip()
                        aObj[k] = v
                        idxV +=1
                    remV -=nK
                    tDict[aType].append(aObj)
                    print "-----------------"
                    for aKey in aObj.keys():
                        print "%s   :   %s "%(aKey, aObj[aKey])
                print "Number of %s is %d "%(aType, len(tDict[aType])) 

    def setPlaneAtmGroupInOneComp(self, tName):

        if self["ccp4CifObj"]["comps"].has_key(tName):
            if self["ccp4CifObj"]["comps"][tName].has_key("plane_atoms"):
                self["ccp4CifObj"]["comps"][tName]["planes"] = {}
                for aPlAtm in self["ccp4CifObj"]["comps"][tName]["plane_atoms"]:
                    if not self["ccp4CifObj"]["comps"][tName]["planes"].has_key(aPlAtm["plane_id"]):
                        self["ccp4CifObj"]["comps"][tName]["planes"][aPlAtm["plane_id"]] = []
                    self["ccp4CifObj"]["comps"][tName]["planes"][aPlAtm["plane_id"]].append(aPlAtm)
   

    def printOneComp(self, tName):

        if self["ccp4CifObj"]["comps"].has_key(tName):

            print "Information on compound %s : \n"%tName

            if self["ccp4CifObj"]["comps"][tName].has_key("atoms"):
                print "The following is information for atoms :"
                print "number of atom is ", len(self["ccp4CifObj"]["comps"][tName]["atoms"])
                i = 0
                for aAtom in  self["ccp4CifObj"]["comps"][tName]["atoms"]:
                    print "For atom ", i
                    i+=1
                    for aKey in aAtom.keys():
                        print "%s   :   %s  "%(aKey, aAtom[aKey])
 
            if self["ccp4CifObj"]["comps"][tName].has_key("bonds"):
                print "The following is information for bonds :"
                i = 0
                for aBond in  self["ccp4CifObj"]["comps"][tName]["bonds"]:
                    print "For bond ", i
                    i+=1
                    for aKey in aBond.keys():
                        print "%s   :   %s  "%(aKey, aBond[aKey])
 

            if self["ccp4CifObj"]["comps"][tName].has_key("angles"):
                print "The following is information for angles :"
                i = 0
                for aAngle in  self["ccp4CifObj"]["comps"][tName]["angles"]:
                    print "For angle ", i
                    i+=1
                    for aKey in aAngle.keys():
                        print "%s   :   %s  "%(aKey, aAngle[aKey])
 

            if self["ccp4CifObj"]["comps"][tName].has_key("tors"):
                print "The following is information for torsions :"
                i = 0
                for aTor in  self["ccp4CifObj"]["comps"][tName]["tors"]:
                    print "For torsion angle ", i
                    i+=1
                    for aKey in aTor.keys():
                        print "%s   :   %s  "%(aKey, aTor[aKey])
 
            if self["ccp4CifObj"]["comps"][tName].has_key("chirs"):
                print "The following is information for chiral centers :"
                i = 0
                for aChi in  self["ccp4CifObj"]["comps"][tName]["chirs"]:
                    print "For chiral center ", i
                    i+=1
                    for aKey in achi.keys():
                        print "%s   :   %s  "%(aKey, aChi[aKey])
 
            if self["ccp4CifObj"]["comps"][tName].has_key("plane_atoms"):
                print "The following is information for plane atoms :"
                i = 0
                for aPA in  self["ccp4CifObj"]["comps"][tName]["plane_atoms"]:
                    print "For plane atom ", i
                    i+=1
                    for aKey in aPA.keys():
                        print "%s   :   %s  "%(aKey, aPA[aKey])

    def printOneModOrLink(self, tName):

        if self["ccp4CifObj"]["mods"].has_key(tName):

            print "Information on modification %s : \n"%tName

            for aKey in self["ccp4CifObj"]["mods"][tName].keys():
                print "\nProp %s"%aKey
                for aProp in self["ccp4CifObj"]["mods"][tName][aKey]:
                    for aEntry in aProp.keys():
                        print "%s   :   %s  "%(aKey, aPA[aKey])
                         

    def getAtomById(self, tId, tName):

        aReturn = None
        for aAtom in self["ccp4CifObj"]["comps"][tName]["atoms"]:
            if aAtom["atom_id"]==tId:
                aReturn = aAtom
                break
        return aReturn

    def getBondByIds(self, tId1, tId2, tName):

        aReturn = None
        for aBond in self["ccp4CifObj"]["comps"][tName]["bonds"]:
            if (aBond["atom_id_1"] == tId1 and aBond["atom_id_2"] == tId2) \
              or (aBond["atom_id_1"] == tId2 and aBond["atom_id_2"] == tId1) :
               aReturn = aBond
               break
        return aReturn
         
       
