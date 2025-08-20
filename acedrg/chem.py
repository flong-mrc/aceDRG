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
from builtins   import str
from builtins   import range
from builtins   import object
import os,os.path,sys

import numpy    as np


from rdkit      import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles


from . utility  import BondOrderS2N
from . utility  import aLineToAlist
from . utility  import getLinkedGroups2
from . utility  import aMolToAGraph
from . utility  import MatchTwoGraphs
from . utility  import getAUnitDirVec
from . utility  import setOneAtomCoordsOutCB
from . utility  import getAUnitDirVec2P
from . utility  import getLengBetween2Pos
from . utility  import getAAngFrom3Ps
from . utility  import getPlaneNormal
from . utility  import replaceAtomCoords
from . utility  import DataStrTransferAtomsAndBonds

from . filetools import outputOneMolInACif
from . filetools import fromACrysToMolCifsGemmi
from . filetools import outputOneMolFullDictInACif

class ChemCheck(object):

    def __init__( self):

        self.organicSec = ["AS", "As", "as", "AT", "At", "at", "B", "b", "BR", "Br", "br", "C", "c", "CL", "Cl", "cl",
                   "D", "d", "F", "f", "Ge", "Ge", "ge", "H", "h", "I", "i", "N","n",  "O", "o", "P", "p", "S", "s",
                   "SE", "Se", "se", "SI", "Si", "si"]
        self.atomFileType = {}
        self.atomFileType["mmCif"]  = [11, 111, 112, 16, 161, 51]
        self.atomFileType["simils"] = [12, 121, 52]
        self.atomFileType["mdl"]    = [13, 131, 14, 141, 33, 53, 54]
        self.atomFileType["mol2"]   = [15, 151, 55]

        self.pTable = rdchem.GetPeriodicTable()

        self.aminoAcids = ["ALA", "ARG", "ASN", "ASP", "CYS", \
                           "GLN", "GLU", "GLY", "HIS", "ILE", \
                           "LEU", "LYS", "MET", "PHE", "PRO", \
                           "SER", "THR", "TRP", "TYR", "VAL"]

        self.sugar = {}
        self.sugar["pyranose"]  = []
        self.sugar["pyranose"].append("OC1C(C(C(C(O1)O)O)O)O")        # Pyranose
        self.sugar["pyranose"].append("C1C(C(C(C(O1)O)O)O)O")         # 0MK, A6P, XYS
        self.sugar["pyranose"].append("C1C(C(C(C(O1)*)O)O)*")         # B1H
        self.sugar["pyranose"].append("C1C(C(C(C(O1)*)*)O)O")         # AGL
        self.sugar["furanoses"] = []
        self.sugar["furanoses"].append("C1CCCO1")

        # Just for checking of amino acids in the old ccp4 monomer lib
        self.defaultBo ={}
        self.defaultBo["H"]  = 1
        self.defaultBo["F"]  = 1
        self.defaultBo["CL"] = 1
        self.defaultBo["BR"] = 1
        self.defaultBo["I"]  = 1
        self.defaultBo["AT"] = 1
        self.defaultBo["O"]  = 2
        self.defaultBo["S"]  = 2
        self.defaultBo["SE"] = 2
        self.defaultBo["N"]  = 3
        self.defaultBo["AS"] = 3
        self.defaultBo["B"]  = 3
        self.defaultBo["C"]  = 4
        self.defaultBo["GE"] = 4
        self.defaultBo["SI"] = 4
        self.defaultBo["P"]  = 5

        # Copy some data from libmol to here
        self.orgVal          = {}
        self.orgVal["AS"]    = [3]
        self.orgVal["AT"]    = [1]
        self.orgVal["B"]     = [3]
        self.orgVal["BR"]    = [1]
        self.orgVal["C"]     = [4]
        self.orgVal["CL"]    = [1]
        self.orgVal["F"]     = [1]
        self.orgVal["GE"]    = [4]
        self.orgVal["H"]     = [1]
        self.orgVal["D"]     = [1]
        self.orgVal["I"]     = [1]
        self.orgVal["N"]     = [3]
        self.orgVal["O"]     = [2]
        self.orgVal["P"]     = [5]    # should be [5, 3]. Use one at the moment
        self.orgVal["S"]     = [2, 4, 6]
        self.orgVal["SE"]    = [2]
        self.orgVal["SI"]    = [4]
        
        
        self.maxVal          = {}
        self.maxVal["AS"]    = 3
        self.maxVal["AT"]    = 1
        self.maxVal["B"]     = 4
        self.maxVal["BR"]    = 1
        self.maxVal["C"]     = 4
        self.maxVal["CL"]    = 1
        self.maxVal["F"]     = 1
        self.maxVal["GE"]    = 4
        self.maxVal["H"]     = 1
        self.maxVal["D"]     = 1
        self.maxVal["I"]     = 1
        self.maxVal["N"]     = 5
        self.maxVal["O"]     = 3
        self.maxVal["P"]     = 5    # should be [5, 3]. Use one at the moment
        self.maxVal["S"]     = 6
        self.maxVal["SE"]    = 4
        self.maxVal["SI"]    = 4
        self.maxVal["OTHERS"]    = 1  # metal elemets 

        self.torsions        = []
        self.outChiSign      ={}
        self.tmpChiralSign   = {}
        
        self.tmpSeriNum      =-1

        # for Carbone-related applications
        self.CBGrapList      = {}
        self.tmpBrokenBondIds = {}
        self.tmpBrokenBonds   = []
        self.tmpHIdMap        = {}
        
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
        """
        if not organicOnly :
            print("The input ligands/molecules contains metal or other heavier atoms ")
            print("Acedrg currently deals with ligands/molecules with following elements only ")
            print("As, B, Br, C, Cl, F, Ge, H, I, N, O, P, S, Si")
            print("The job finishes without the output cif file")
        if len(allAtomElems) ==0:
            print ("Can not get the element symbols of atoms. Check input file format")
            print("Error : The job stops because of errors")
            organicOnly = False
        """

        return organicOnly


    def isOrganicInCif(self, tCifAtoms):

        aRet = True
        lElem = True
        nonOrgSet =[]
        for aAtm in tCifAtoms:
            if "_chem_comp_atom.type_symbol" in aAtm.keys():
                 if not aAtm["_chem_comp_atom.type_symbol"] in self.organicSec :
                     if not aAtm["_chem_comp_atom.type_symbol"] in nonOrgSet:
                         nonOrgSet.append(aAtm["_chem_comp_atom.type_symbol"])
            elif "type_symbol" in aAtm.keys():
                 if not aAtm["type_symbol"] in self.organicSec :
                     if not aAtm["type_symbol"] in nonOrgSet:
                         nonOrgSet.append(aAtm["type_symbol"])
            else:
                 lElem = False

        if len(nonOrgSet)> 0:
            aRet = False
        #if not aRet :
        #    print("The input ligands/molecules contains metal or other heavier atoms ")
        #    print("Acedrg currently deals with ligands/molecules with following elements only ")
        #    print("As, B, Br, C, Cl, F, Ge, H, I, N, O, P, S, Si")
        #    print ("Your molecule contains atoms of elements of the following type_symbol:")
        #    for aE in nonOrgSet:
        #        print(aE)
        #print("The job finishes without the output cif file.")

        if not lElem :
            print ("Can not get the element symbols of atoms. Check input file format")
            print("Error : The job stops because of errors")
            organicOnly = False

        return aRet

    def isOrganicMol(self, tMol):

        pass

    def getMetalAtoms(self, tCifAtoms, tMetalAtoms):

        for aAtm in tCifAtoms:
            if "_chem_comp_atom.type_symbol" in aAtm.keys():
                 if not aAtm["_chem_comp_atom.type_symbol"] in self.organicSec :
                     tMetalAtoms.append(aAtm)

    def addHs(self, tMol):

        # print "Number of atom before is ", tMol.GetNumAtoms()
        aMolT=Chem.MolFromSmiles(rdmolfiles.MolToSmiles(tMol))
        print(rdmolfiles.MolToSmiles(tMol))
        aMol = Chem.AddHs(aMolT)
        rdmolfiles.MolToMolFile(tMol, "1.mol")
        print("Number of atom now is ", aMol.GetNumAtoms())
        tDoneAtms = []
        for i in range(len(tMol.GetAtoms())):
            elem = tMol.GetAtomWithIdx(i).GetSymbol()
            if elem != "H":
                aBO = tMol.GetAtomWithIdx(i).GetTotalValence()
                aDV = self.pTable.GetDefaultValence(elem)
                nAddHs = aDV - aBO
                tMol.GetAtomWithIdx(i).SetFormalCharge(nAddHs)

        aMolT = Chem.AddHs(tMol)
        print("Number of atom after is ", aMolT.GetNumAtoms())

    def getNumNBHAtoms(self, tAtom):

        nRet =0

        nbAtoms = tAtom.GetNeighbors()
        for aNBAtm in nbAtoms:
            if aNBAtm.GetSymbol()=="H":
                nRet+=1
        return nRet

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
            print("Could not recognize format for the input file. It should be one of cif, smiles, sdf/mol, mol2")
            sys.exit(1)
        aLine = ""
        for aEl in tAtomElems:
             aLine += (aEl.strip()+"\t")
        if len(aLine.strip()) > 0:
            print("The system contains atoms of the following elements")
            print(aLine)

    def getAtomElemsFromMmcif(self, tInFileName, tAtomElems):
        #print(tInFileName)
        try :
            inFile = open(tInFileName, "r")
        except IOError :
             print("%s can not be opened for reading "%tInFileName)
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
                    #strGrp = aL.strip().split()
                    strGrp = []
                    aLineToAlist(aL, strGrp)
                    if "type_symbol" in colDict and len(strGrp) == len(colDict) \
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
            print("%s can not be opened for reading "%tInFileName)
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
                            print("Format error is input MOL/SDF file. The count line is : ")
                            print(aL)
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
                print("% can not be open for reading "%tFileName)
                sys.exit()
            else:
                aSmiStr = fSmi.read()
                tSmiStr = aSmiStr.strip().split()
                if len(tSmiStr) >0 :
                    aSmiStr = tSmiStr[0].strip()
                else:
                    print("String format error")
                    sys.exit()
                fSmi.close()
        else:
            # SMILES string in from a commandline
            aSmiStr = tFileName.strip()
        if len(aSmiStr):
            checkMol=Chem.MolFromSmiles(aSmiStr, sanitize=False)

            if checkMol:
                for aAtom in checkMol.GetAtoms():
                    aSym = aAtom.GetSymbol().strip()
                    if len(aSym) and not aSym in tAtomElems:
                        tAtomElems.append(aSym)
            else:
                print("Can not generate a molecule from the input SMILES string!")
                print("Check your SMILES or report the bug ")
    
    def checkCarbonMolFromSmi(self, tFileName):
        
        aRet = False
        if os.path.isfile(tFileName):
            # SMILES string in a file
            try:
                fSmi = open(tFileName, "r")
            except IOError:
                print("% can not be open for reading "%tFileName)
                sys.exit()
            else:
                aSmiStr = fSmi.read()
                tSmiStr = aSmiStr.strip().split()
                if len(tSmiStr) >0 :
                    aSmiStr = tSmiStr[0].strip()
                else:
                    print("String format error")
                    sys.exit()
                fSmi.close()
        else:
            # SMILES string in from a commandline
            aSmiStr = tFileName.strip()
        if len(aSmiStr):
            checkMol=Chem.MolFromSmiles(aSmiStr, sanitize=False)
            
            allNonHAtoms = checkMol.GetAtoms()
            allNonHBonds = checkMol.GetBonds()
            aElmMap = {}
            connBMap = {}
            for aAt in allNonHAtoms:
                aElmMap[aAt.GetIdx()] = aAt.GetSymbol()
                connBMap[aAt.GetIdx()] = []
                
            for aBo in allNonHBonds:
                idx1 = aBo.GetBeginAtomIdx()
                idx2 = aBo.GetEndAtomIdx()
                if aElmMap[idx1] == "B":
                    connBMap[idx2].append(idx1)
                if aElmMap[idx2] == "B":
                    connBMap[idx1].append(idx2)
            
            for aIdx in connBMap:
                if len(connBMap[aIdx]) >=4:
                    aRet = True
                    break
        return aRet 
    
    def GetAtomsAndBondsFromSmi(self, tFileName, tAtoms, tBonds):
        
        # Transfer rdkit data structures to a dictionary-style (e.g. gemmi atoms)
        
        
        
        if os.path.isfile(tFileName):
            # SMILES string in a file
            try:
                fSmi = open(tFileName, "r")
            except IOError:
                print("% can not be open for reading "%tFileName)
                sys.exit()
            else:
                aSmiStr = fSmi.read()
                tSmiStr = aSmiStr.strip().split()
                if len(tSmiStr) >0 :
                    aSmiStr = tSmiStr[0].strip()
                else:
                    print("String format error")
                    sys.exit()
                fSmi.close()
                print(aSmiStr)
        else:
            # SMILES string in from a commandline
            aSmiStr = tFileName.strip()
            print(aSmiStr)
            
        initMol=Chem.MolFromSmiles(aSmiStr, sanitize=False)
        sys.exit()
        aMol= Chem.AddHs(initMol)
            
        print(aMol.GetNumAtoms())
            
        aIdSet = {}
        for aAt in aMol.GetAtoms():
            aDAt = {}
            aDAt["_chem_comp_atom.type_symbol"] = aAt.GetSymbol()
            aName = ""
            if not aDAt["_chem_comp_atom.type_symbol"] in aIdSet:
                aIdSet[aDAt["_chem_comp_atom.type_symbol"]] = 1
            else:
                aIdSet[aDAt["_chem_comp_atom.type_symbol"]] +=1
            if aIdSet[aDAt["_chem_comp_atom.type_symbol"]] == 1:
                aName = aDAt["_chem_comp_atom.type_symbol"] 
            else:
                aName = aDAt["_chem_comp_atom.type_symbol"] + str(aIdSet[aDAt["_chem_comp_atom.type_symbol"]])
            aDAt["_chem_comp_atom.atom_id"] = aName
                
            aDAt["_chem_comp_atom.type_energy"] = aDAt["_chem_comp_atom.type_symbol"]
            aDAt["_chem_comp_atom.atom_serial_number"] = aAt.GetIdx()
            
            tAtoms.append(aDAt)
                
        for aBo in aMol.GetBonds():
            aDB = {}
            aDB["_chem_comp_bond.bond_serial_number"] = aBo.GetIdx()
            idx1 = aBo.GetBeginAtomIdx()
            id1  = tAtoms[idx1]["_chem_comp_atom.atom_id"]
            idx2 = aBo.GetEndAtomIdx()
            id2  = tAtoms[idx2]["_chem_comp_atom.atom_id"]
            aDB["_chem_comp_bond.atom_serial_number_1"] = idx1
            aDB["_chem_comp_bond.atom_serial_number_2"] = idx2
            aDB["_chem_comp_bond.atom_id_1"] = id1
            aDB["_chem_comp_bond.atom_id_2"] = id2
            tAtoms[idx1]["_chem_comp_atom.atom_conn"].append(idx2)
            tAtoms[idx2]["_chem_comp_atom.atom_conn"].append(idx1)
            tBonds.append(aDB)
            
            print("Initial atoms from smi : ", len(tAtoms))
            print("They are : ")
                
                
    
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
            print("%s can not be found for reading "%tFileName)


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

    def checkLocalStructAndCharge(self, tAtomId, tAtoms,  tBondSet1, tBondSet2):

        # For a linked atom, the environment mostly changes.
        # Formal charge should be corrected respectively

        pass

    def getAllowBondOrderForOneOrgAtom(self, tAtom):

        allowedBo = []
        aT = tAtom["type_symbol"]
        if aT in self.orgVal:
            for aOrd in self.orgVal[aT]:
                if "charge" in tAtom:
                    allowedBo.append(aOrd + int(tAtom["charge"]))
                else:
                    allowedBo.append(aOrd)

        return allowedBo

    def valideBondOrderForOneOrgAtom(self, tAtom, tBonds):

        aBool = True
        aLine = ""

        nOrder = 0
        for aB in tBonds:
            nBT = BondOrderS2N(aB["type"])
            if nBT != -1 :
                if "charge" in tAtom:
                    nOrder += (nBT + int(tAtom["charge"]))
                else:
                    nOrder += nBT
            else:
                #aLine = "Unknown type %s for the added bond between %s and %s \n"\
                #        %(aPair[0]["type"], aPair[0]["atom_id_1"], aPair[0]["atom_id_2"])
                aBool = False
                break
        if aBool:
            print("Total Bond order (including charges) for atom %s is %d"%(tAtom["atom_id"], nOrder))
            curAllowOrds = self.getAllowBondOrderForOneOrgAtom(tAtom)
            if len(curAllowOrds) !=0:
                print("Total allowed bond orders are : ")
                for aOrd in curAllowOrds:
                    print(aOrd)
                if not nOrder in curAllowOrds:
                     aLine ="Bond order (valence) %d around added atom %s are wrong.\n"%(nOrder, tAtom["atom_id"])
                     aLine += "It should be %d, check your input file(bond order or charge)\n"%curAllowOrds[0]
                     aBool = False
            else:
                aLine = "Bug, can not find the allowed valence for added atom %s\n"%tAtom["atom_id"]
                aBool = False

        return [aBool, aLine]

    def valideBondOrderForOneOrgAtom2(self, tAtom, tBonds):

        aBool = True
        aLine = ""

        nOrder = 0
        for aB in tBonds:
            nBT = BondOrderS2N(aB["type"])
            if nBT != -1 :
                if "charge" in tAtom:
                    nOrder += (nBT + int(tAtom["charge"]))
                else:
                    nOrder += nBT
            else:
                #aLine = "Unknown type %s for the added bond between %s and %s \n"\
                #        %(aPair[0]["type"], aPair[0]["atom_id_1"], aPair[0]["atom_id_2"])
                aBool = False
                break
        if aBool:
            print("Total Bond order (including charges) for atom %s is %d"%(tAtom["atom_id"], nOrder))
            curAllowOrds = self.getAllowBondOrderForOneOrgAtom(tAtom)
            if len(curAllowOrds) !=0:
                print("Total allowed bond orders are : ")
                for aOrd in curAllowOrds:
                    print(aOrd)
                if not nOrder in curAllowOrds:
                     difBo =  nOrder - curAllowOrds[0]
                     aLine ="Bond order (valence) %d around added atom %s are wrong.\n"%(nOrder, tAtom["atom_id"])
                     aLine += "It should be %d, check your input file(bond order or charge)\n"%curAllowOrds[0]
                     aBool = False
            else:
                aLine = "Bug, can not find the allowed valence for added atom %s\n"%tAtom["atom_id"]
                aBool = False

        return [aBool, aLine, difBo]

    def adjustNBForOneAddedAtom(self, tAtom, tNBAtoms, tBonds, tRes, tMod, tDelAtmIds):

        aBool = True
        aLine = ""


        nOrder = 0
        for aB in tBonds:
            nBT = -1
            if "type" in aB:
                nBT = BondOrderS2N(aB["type"])
                #print("Bond-type between atom %s and %s is %d "%(aB["atom_id_1"], aB["atom_id_2"], nBT))
            if nBT != -1 :
                nOrder += nBT
            else:
                aLine = "Unknown type for the added bond between %s and %s \n"\
                        %(aB["atom_id_1"], aB["atom_id_2"])
                aBool = False
                break

        if aBool:
            print("Total Bond order (including charges) for atom %s is %d"%(tAtom["atom_id"], nOrder))
            curAllowOrds = self.getAllowBondOrderForOneOrgAtom(tAtom)
            if len(curAllowOrds) !=0:
                print("Total allowed bond orders are : ")
                for aOrd in curAllowOrds:
                    print(aOrd)
                #print("Current bond order is : ", nOrder)
                if not nOrder in curAllowOrds:
                    # If total valence does not equal to one of allowed valences. Use curAllowOrds[0]
                    # The procedures are:
                    # (1) if the atom has formal charge, dealt with it first
                    # (2) adjust H atom around the atom
                    # (3) error info : tell the user to re-define the added bonds
                    defOrd = curAllowOrds[0] - nOrder
                    if defOrd == 1:
                        if tAtom["charge"]=="1":
                            # (1a) nOrder = curAllowOrds[0] - 1, cancel the formal charge in atom
                            tAtom["charge"]="0"
                            tMod["changed"]["atoms"].append(tAtom)
                            print("The formal charge at atom %s is reduced to %s"%(tAtom["atom_id"], tAtom["charge"]))
                        else:
                            # (1b) to add H atoms and therefore the valence of the atom.
                            aAtom = {}

                            aId = ""
                            if len(tAtom["atom_id"]) > 2:
                                aId = tAtom["atom_id"][1:]
                            else :
                                aId = tAtom["atom_id"]
                            aAtom["atom_id"]     = "H" + aId
                            aAtom["type_symbol"] = "H"
                            aAtom["type_energy"] = aAtom["type_symbol"]
                            aAtom["comp_id"]     = tAtom["comp_id"]
                            aAtom["charge"]      = "0"
                            tMod["added"]["atoms"].append(aAtom)
                            tRes["remainAtoms"].append(aAtom)
                            print("One H atom %s is added into the residue %s "%(aAtom["atom_id"], aAtom["comp_id"]))
                            aBond = {}
                            aBond["atom_id_1"] = aAtom["atom_id"]
                            aBond["atom_id_2"] = tAtom["atom_id"]
                            aBond["type"]      = "SINGLE"
                            aBond["value_dist"] = "1.0"
                            aBond["value_dist_esd"] = "0.01"
                            tMod["added"]["bonds"].append(aBond)
                            tRes["remainBonds"].append(aBond)
                    elif defOrd > 1:
                        for i in range(defOrd):
                            aAtom = {}
                            aId = ""
                            if len(tAtom["atom_id"]) > 1:
                                aId = tAtom["atom_id"][1:]
                            else :
                                aId = tAtom["atom_id"]
                            aAtom["atom_id"]     = "H" + aId + str(i+1)         # tempo
                            aAtom["type_symbol"] = "H"
                            aAtom["type_energy"] = aAtom["type_symbol"]
                            aAtom["comp_id"]     = tRes["name"]
                            aAtom["charge"]      = "0"
                            tMod["added"]["atoms"].append(aAtom)
                            tRes["remainAtoms"].append(aAtom)
                            print("One H atom %s is added into the residue %s "%(aAtom["atom_id"], aAtom["comp_id"]))
                            aBond = {}
                            aBond["atom_id_1"] = aAtom["atom_id"]
                            aBond["atom_id_2"] = tAtom["atom_id"]
                            aBond["type"]      = "SINGLE"
                            aBond["value_dist"] = "1.0"
                            aBond["value_dist_esd"] = "0.01"
                            tMod["added"]["bonds"].append(aBond)
                            tRes["remainBonds"].append(aBond)
                    elif nOrder > curAllowOrds[0]:
                        if len(curAllowOrds)==1:
                            aN = nOrder-curAllowOrds[0]
                            print("%d H atom will be deleted "%aN)
                            hAtomIds =[]
                            for aAtom in tNBAtoms:
                                if aAtom["type_symbol"]=="H" and not aAtom["atom_id"] in tDelAtmIds:
                                    hAtomIds.append(aAtom["atom_id"])
                            print("%s connects %d H atom "%(tAtom["atom_id"], len(hAtomIds)))
                            if aN <=len(hAtomIds) :
                                hAtomIds.sort()
                                idxD = -1
                                for i in range(aN):
                                    aHName = hAtomIds[idxD]
                                    if not aHName in tDelAtmIds:
                                        tDelAtmIds.append(aHName)
                                        print("H atom %s is deleted "%aHName)
                                    idxD = idxD -1
                                for aAtom in tRes["comp"]["atoms"]:
                                    if aAtom["atom_id"] in tDelAtmIds:
                                        tMod["deleted"]["atoms"].append(aAtom)
                            else :
                                 aLine = "Currently it is not allowed to do such a change for atom %s %s\n"%tAtom["atom_id"]
                                 aBool = False
                        else:

                            nAdded = 1 # TEMPO
                            lAdded = False
                            for i in range(1,len(curAllowOrds)):
                                if nOrder == curAllowOrds[i] -1:
                                    self.addOneHAtomAndRelatedBond(tMod["added"]["atoms"], tMod["added"]["bonds"], tAtom, nAdded)
                                    nAdded +=1
                                    lAdded = True
                                    break
                            if not lAdded:
                                aLine = "Currently it is not allowed to do such a change for atom %s %s\n"%tAtom["atom_id"]
                                aBool = False
            else:
                aLine = "Bug, can not find the allowed valence for added atom %s\n"%tAtom["atom_id"]
                aBool = False
        return [aBool, aLine]

    def addOneHAtomAndRelatedBond(self, tAtoms, tBonds, tConnAtom, tIdxH):
        # Add a H atom that connect tConnAtom
        aAtom = {}
        aStr = str(tIdxH)
        aAtom["atom_id"]     = "HX" + aStr          # tempo
        aAtom["type_symbol"] = "H"
        aAtom["type_energy"] = aAtom["type_symbol"]
        aAtom["comp_id"]     = tConnAtom["comp_id"]
        tAtoms.append(aAtom)
        print("One H atom %s is added into the residue %s "%(aAtom["atom_id"], aAtom["comp_id"]))
        aBond = {}
        aBond["atom_id_1"] = aAtom["atom_id"]
        aBond["atom_id_2"] = tConnAtom["atom_id"]
        aBond["type"]      = "SINGLE"
        tBonds.append(aBond)

    def addOneHAtomAndRelatedBond2(self, tAtoms, tBonds, tConnAtom, tIdxH):
        # Add a H atom that connect tConnAtom
        aAtom = {}
        aStr = str(tIdxH)
        aAtom["_chem_comp_atom.atom_id"]     = "HX" + aStr          # tempo
        aAtom["_chem_comp_atom.type_symbol"] = "H"
        aAtom["_chem_comp_atom.comp_id"]     = tConnAtom["_chem_comp_atom.comp_id"]
        aAtom["_chem_comp_atom.x"] = "0.0"
        aAtom["_chem_comp_atom.y"] = "0.0"
        aAtom["_chem_comp_atom.z"] = "0.0"
        tAtoms.append(aAtom)

        print("One H atom %s is added into the residue %s "%(aAtom["_chem_comp_atom.atom_id"], aAtom["_chem_comp_atom.comp_id"]))
        aBond = {}
        aBond["_chem_comp_bond.atom_id_1"] = aAtom["_chem_comp_atom.atom_id"]
        aBond["_chem_comp_bond.atom_id_2"] = tConnAtom["_chem_comp_atom.atom_id"]
        aBond["_chem_comp_bond.value_order"]      = "SING"
        tBonds.append(aBond)

    def tmpModiN_in_AA(self, tCompId, tCifMonomer):

        # Tempo function to add a H atom connected to N in AA coming from the current $CLIBD_MON
        aAtom = {}
        aAtom["comp_id"]     = tCompId
        aAtom["atom_id"]     = "H2"
        aAtom["type_symbol"] = "H"
        aAtom["type_energy"] = "HNH2"
        #aAtom["charge"]      = "0.0"
        tCifMonomer["atoms"].append(aAtom)
        aBond = {}
        aBond["comp_id"]         = tCompId
        aBond["atom_id_1"]       = "N"
        aBond["atom_id_2"]       = "H2"
        aBond["type"]            = "single"
        aBond["value_dist"]      = "0.860"
        aBond["value_dist_esd"]  = "0.02"
        tCifMonomer["bonds"].append(aBond)

        for aAtm in tCifMonomer["atoms"]:
            if aAtm["atom_id"].strip()=="N":
                aAtm["type_energy"] = "NH2"
            elif aAtm["atom_id"].strip()=="H":
                aAtm["type_energy"] = "HNH2"

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
                nHs = 0
                for aNB in aNBSet:
                    symb = aNB.GetSymbol().strip()
                    if len(symb)==1 and symb.find("H") !=-1:
                        nHs +=1
                #print("atom symb ", aSymb)
                #print("number of H NB ", nHs)
                if not lSP2 and nHs < 2:
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
            print("Can not calculate chiral volume. atom coordinates are not 3 dim ")

        return vSign


    def setHName(self, tHConnAtom, tOtheAtmSet, tAllHIds):

        nAllH =len(tAllHIds)
        reName = "H" + str(nAllH +1)
        hIds = []
        for aA in tOtheAtmSet:
            if aA["type_symbol"].strip() == "H":
                 hIds.append(aA["atom_id"])

        hIds.sort()
        nL = len(tHConnAtom["atom_id"])
        nH = len(hIds)
        if nH ==0:
           aId = "H" + str(nAllH)
           if aId in tAllHIds:
               reName = "H" + str(nAllH +1)
           else:
               if nL==1:
                   aId1 = "H"
                   aId2 = "H" + tHConnAtom["type_symbol"]
                   if not aId1 in tAllHIds:
                       reName = aId1
                   elif not aId2 in tAllHIds:
                       reName = aId2
        elif nH==1:
             aId = hIds[0] + "2"
             if not aId in tAllHIds:
                 reName = aId
        else:

            if hIds[-1][-1].isdigit():
                aN = int(hIds[-1][-1]) + 1
                aId = hIds[-1][0:-1] + str(aN)
                if not aId in tAllHIds:
                    reName = aId

        return reName

    def setHName2(self, tHConnAtom, tOtheAtmSet, tAllHIds):

        nAllH =len(tAllHIds)
        reName = "H" + str(nAllH +1)
        print("HHHH ", reName)
        return reName

    def confirmAAandNames(self, tAtoms, tBonds):

        lAA = False
        atomLinks = {}
        baseAtomIds = []
        atomDicts   = {}
        hAtomNameMap = {}
        hBondNameMap = {}
        # Check AA atoms exist
        i=0
        for aAtom in tAtoms:
            atomDicts[aAtom["_chem_comp_atom.atom_id"]] =i
            i = i+1
            if "_chem_comp_atom.atom_id" in aAtom.keys() and\
               "_chem_comp_atom.type_symbol" in aAtom.keys():
                # Check CA
                if aAtom["_chem_comp_atom.atom_id"].strip().upper()== "CA"\
                   and aAtom["_chem_comp_atom.type_symbol"].strip().upper()=="C":
                    baseAtomIds.append("CA")
                elif aAtom["_chem_comp_atom.atom_id"].strip().upper()== "C"\
                   and aAtom["_chem_comp_atom.type_symbol"].strip().upper()=="C":
                    baseAtomIds.append("C")
                elif aAtom["_chem_comp_atom.atom_id"].strip().upper()== "O"\
                   and aAtom["_chem_comp_atom.type_symbol"].strip().upper()=="O":
                    baseAtomIds.append("O")
                elif aAtom["_chem_comp_atom.atom_id"].strip().upper()== "OXT"\
                   and aAtom["_chem_comp_atom.type_symbol"].strip().upper()=="O":
                    baseAtomIds.append("OXT")
                elif aAtom["_chem_comp_atom.atom_id"].strip().upper()== "N"\
                   and aAtom["_chem_comp_atom.type_symbol"].strip().upper()=="N":
                    baseAtomIds.append("N")
        # Check connections
        if len(baseAtomIds)==5 and "CA" in baseAtomIds and "C" in baseAtomIds\
           and "O" in baseAtomIds and "OXT" in baseAtomIds and "N" in baseAtomIds:
            for aAId in baseAtomIds:
                atomLinks[aAId] = []
                for aBond in tBonds:
                    if "_chem_comp_bond.atom_id_1" in aBond.keys() and\
                       "_chem_comp_bond.atom_id_2" in aBond.keys():
                        if aBond["_chem_comp_bond.atom_id_1"].strip().upper()==aAId:
                            atomLinks[aAId].append(aBond["_chem_comp_bond.atom_id_2"])
                        elif aBond["_chem_comp_bond.atom_id_2"].strip().upper()==aAId:
                            atomLinks[aAId].append(aBond["_chem_comp_bond.atom_id_1"])
            if "N" in atomLinks["CA"] and "C" in atomLinks["CA"] and len(atomLinks["CA"])==4\
               and "O" in atomLinks["C"] and "OXT" in atomLinks["C"] and len(atomLinks["C"])==3\
               and len(atomLinks["N"])<=4:
                lAA = True
                # Check and change if required
                hAtomNameMap["CA"]= []
                for aId in atomLinks["CA"]:
                    if tAtoms[atomDicts[aId]]["_chem_comp_atom.type_symbol"]=="H":
                        hAtomNameMap["CA"].append(tAtoms[atomDicts[aId]]["_chem_comp_atom.atom_id"])
                if len(hAtomNameMap["CA"]) !=1:
                    lAA = False
                elif hAtomNameMap["CA"][0] != "HA":
                    lAA = False
                hAtomNameMap["N"]= []
                for aId in atomLinks["N"]:
                    if tAtoms[atomDicts[aId]]["_chem_comp_atom.type_symbol"]=="H":
                        hAtomNameMap["N"].append(tAtoms[atomDicts[aId]]["_chem_comp_atom.atom_id"])
                if len(atomLinks["N"])==3:
                    if len(hAtomNameMap["N"]) ==0 or len(hAtomNameMap["N"])> 3 :
                        lAA = False
                    else:
                        if len(hAtomNameMap["N"]) ==1 :
                            if not "H2" in hAtomNameMap["N"] :
                                lAA = False
                        else:
                            if not "H" in hAtomNameMap["N"] or not "H2" in hAtomNameMap["N"]:
                                lAA = False
                elif len(atomLinks["N"])==4:
                    if len(hAtomNameMap["N"]) != 3  :
                        lAA = False
                    else:
                        if not "H" in hAtomNameMap["N"] or not "H2" in hAtomNameMap["N"] or not "H3" in hAtomNameMap["N"]:
                            lAA = False
            else:
                lAA = False
        else:
            lAA = False
        return lAA

    def containAROMA(self, tBonds):

        aRet = False

        for aB in tBonds:

            if "_chem_comp_bond.value_order" in aB.keys():
                if aB["_chem_comp_bond.value_order"].upper()[:4]=="AROM":
                    aRet = True
                    break


        return aRet

    def addjustAtomsAndBonds2(self, tAtoms, tBonds, tBOHFile):

        tBOHFile = tBOHFile + "_BOs_and_Hs.list"
        if os.path.isfile(tBOHFile):
            f = open(tBOHFile, "r")
            allLs = f.readlines()
            f.close()

            modBonds = {}
            addHAtms = {}
            for aL in allLs:
                strs = aL.strip().split()
                if len(strs) >0:
                    if strs[0].find("BOND-ORDER:") !=-1:
                        combId = strs[1] + "_" + strs[2]
                        if not combId in modBonds.keys():
                            modBonds[combId] = strs[3]
                    elif strs[0].find("HATOM:") !=-1:
                        nH = int(strs[2])
                        if nH > 0:
                            addHAtms[strs[1]] = nH
            print(addHAtms)

            for aB in tBonds:
                if "_chem_comp_bond.value_order" in aB.keys():
                    atm1Id = aB["_chem_comp_bond.atom_id_1"].strip()
                    atm2Id = aB["_chem_comp_bond.atom_id_2"].strip()
                    combId1 = atm1Id + "_" + atm2Id
                    combId2 = atm1Id + "_" + atm2Id
                    combId  = ""
                    if combId1 in modBonds.keys():
                        aB["_chem_comp_bond.value_order"]=modBonds[combId1]
                    elif combId2 in modBonds.keys():
                        aB["_chem_comp_bond.value_order"]=modBonds[combId2]
            # Add H atoms and related bpnds
            atmIdMap = {}
            for i in range(len(tAtoms)):
                aId = tAtoms[i]["_chem_comp_atom.atom_id"]
                atmIdMap[aId] = i
            compId = tAtoms[0]["_chem_comp_atom.comp_id"]
            for aAtmId in addHAtms.keys():
                for j in range(addHAtms[aAtmId]):
                    aHAtom = {}
                    aHAtom["_chem_comp_atom.comp_id"] = compId
                    idxA   = atmIdMap[aAtmId]
                    elemA  = tAtoms[idxA]["_chem_comp_atom.type_symbol"]
                    part2Name = '"'
                    if elemA=="C":
                        part2Name = tAtoms[idxA]["_chem_comp_atom.atom_id"][1:]
                    else:
                        part2Name = tAtoms[idxA]["_chem_comp_atom.atom_id"]
                    sj = str(j+1)
                    if addHAtms[aAtmId] ==1:
                        aHAtom["_chem_comp_atom.atom_id"] = "H"+part2Name
                    else:
                        if len(part2Name) < 3:
                            aHAtom["_chem_comp_atom.atom_id"] = "H"+part2Name + sj
                        else:
                            aHAtom["_chem_comp_atom.atom_id"] = "H"+part2Name[0:-2] + sj
                    aHAtom["_chem_comp_atom.type_symbol"] = "H"
                    aHAtom["_chem_comp_atom.charge"]      = "0"
                    aHAtom["_chem_comp_atom.x"]           = "0.0"
                    aHAtom["_chem_comp_atom.y"]           = "0.0"
                    aHAtom["_chem_comp_atom.z"]           = "0.0"
                    tAtoms.append(aHAtom)

                    aHBond = {}
                    aHBond["_chem_comp_bond.comp_id"]   = compId
                    aHBond["_chem_comp_bond.atom_id_1"] = tAtoms[idxA]["_chem_comp_atom.atom_id"]
                    aHBond["_chem_comp_bond.atom_id_2"] = aHAtom["_chem_comp_atom.atom_id"]
                    aHBond["_chem_comp_bond.type"]      = "SING"
                    aHBond["_chem_comp_bond.aromatic"]  = "n"
                    # the following values are just place holders
                    aHBond["_chem_comp_bond.value_dist_nucleus"]     = "1.0"
                    aHBond["_chem_comp_bond.value_dist_nucleus_esd"] = "0.1"
                    aHBond["_chem_comp_bond.value_dist"]             = "1.0"
                    aHBond["_chem_comp_bond.value_dist_esd"]         = "0.1"
                    tBonds.append(aHBond)

    def addjustAtomsAndBonds3(self, tAtoms, tBonds, tBCFName):


        tBCFName  = tBCFName + "_ac.txt"
        if os.path.isfile(tBCFName):
            f = open(tBCFName, "r")
            allLs = f.readlines()
            f.close()

        aAtmMap = {}
        idxA =0
        for aA in tAtoms:
            idStr = aA["_chem_comp_atom.atom_id"]
            aAtmMap[idStr] = idxA
            idxA+=1

        aBondMap = {}
        idxB =0;
        for aB in tBonds:
            idStr = aB["_chem_comp_bond.atom_id_1"].strip() + "_" +\
                     aB["_chem_comp_bond.atom_id_2"].strip()
            aBondMap[idStr] =  idxB
            idxB+=1

        for aL in allLs:
            if aL.upper().find("CHARGE:") !=-1:
                strGrp = aL.strip().split()
                if len(strGrp) ==3:
                    if strGrp[1] in aAtmMap.keys():
                        idxAt = aAtmMap[strGrp[1]]
                        #aCha = int(strGrp[2].strip())
                        tAtoms[idxAt]["_chem_comp_atom.charge"] = strGrp[2].strip()
                        #if aCha !=0:
                        #    tCharges[tAtoms[idxAt]["_chem_comp_atom.atom_id"]] = aCha
                    else:
                        print("Bug: atom %s can not be found"%strGrp[1])
            elif aL.upper().find("BOND:") !=-1:
                strGrp = aL.strip().split()
                if len(strGrp) ==4:
                    idStr =""
                    idStr1 =strGrp[1].strip() + "_" + strGrp[2].strip()
                    idStr2 =strGrp[2].strip() + "_" + strGrp[1].strip()
                    if idStr1 in aBondMap.keys():
                        idStr= idStr1
                    elif idStr2 in aBondMap.keys():
                        idStr= idStr2

                    if idStr:
                        idxBo = aBondMap[idStr]
                        tBonds[idxBo]["_chem_comp_bond.type"] = strGrp[3].strip()
                    else:
                        print("Bug: bond between atoms %s and %s can not be found"
                              %(strGrp[1], strGrp[2]))

        print("After keku: ")
        for aAt in tAtoms:
            print("Atom %s has charge of %s"%(aAt["_chem_comp_atom.atom_id"],  aAt["_chem_comp_atom.charge"]))
        for aB in tBonds:
            print("Bond-order between atoms %s and %s is %s"
                  %(aB["_chem_comp_bond.atom_id_1"], aB["_chem_comp_bond.atom_id_2"],
                    aB["_chem_comp_bond.type"]))


    def addjustAtomsAndBonds(self, tAtoms, tBonds):

        # For molecules contain bonds called explicitly aromatic bonds

        aCurVaMap   = {}
        adefVaMap   = {}
        aChargeMap  = {}
        aElemMap    = {}
        aIdAtmMap   = {}
        aHMap       = {}
        doneAtmIds   = []

        for aAtm in tAtoms:
            atmId = aAtm["_chem_comp_atom.atom_id"]
            aElemMap [atmId]  = aAtm["_chem_comp_atom.type_symbol"]
            aIdAtmMap[atmId]  = [aAtm]
            aCurVaMap[atmId]  = 0
            aChargeMap[atmId] = 0
            if aAtm["_chem_comp_atom.type_symbol"].upper() in self.defaultBo.keys():
                adefVaMap[atmId] = self.defaultBo[aAtm["_chem_comp_atom.type_symbol"]]
                aCurVaMap[atmId] = adefVaMap[atmId]
            if "_chem_comp_atom.charge" in aAtm.keys():
                aChargeMap[atmId] = (float(aAtm["_chem_comp_atom.charge"]))
                aCurVaMap[atmId] -=aChargeMap[atmId]

            else:
                print("The input molecule contains atom %s of element %s, which acedrg can not deal with at the moment.")
                sys.exit(1)


        aromAtmMap = {}
        allAtmBondingMap = {}
        idxB = 0
        for aB in tBonds:
            if "_chem_comp_bond.value_order" in aB.keys():
                atm1Id = aB["_chem_comp_bond.atom_id_1"].strip()
                atm2Id = aB["_chem_comp_bond.atom_id_2"].strip()
                aBO    = aB["_chem_comp_bond.value_order"].strip().upper()
                #print("Bond between %s and %s"%(atm1Id,atm2Id))
                if len(aBO) > 4:
                    aBO = aBO[:4]
                aBON = BondOrderS2N(aBO)

                if not atm1Id in allAtmBondingMap.keys():
                    allAtmBondingMap[atm1Id] = {}
                allAtmBondingMap[atm1Id][atm2Id] = [aBON, idxB, False]
                if not atm2Id in allAtmBondingMap.keys():
                    allAtmBondingMap[atm2Id] = {}
                allAtmBondingMap[atm2Id][atm1Id] = [aBON, idxB, False]

                if aBO.find("AROM") !=-1:
                    if not atm1Id in aromAtmMap.keys():
                        aromAtmMap[atm1Id] = []
                    if not atm2Id in aromAtmMap.keys():
                        aromAtmMap[atm2Id] = []
                    aromAtmMap[atm1Id].append(atm2Id)
                    aromAtmMap[atm2Id].append(atm1Id)
            idxB+=1


        #print(aromAtmMap.keys())
        #print(allAtmBondingMap.keys())

        for aAtm in tAtoms:
            self.setInitNumH(aAtm, allAtmBondingMap,
                             aromAtmMap, aHMap, adefVaMap,aCurVaMap)
            aAtmId = aAtm["_chem_comp_atom.atom_id"]
            if not aAtmId in aromAtmMap.keys():
                doneAtmIds.append(aAtmId)

        for aAtmId in aHMap.keys():
            print("Atom %s needs to add %d Hs"%(aAtmId, aHMap[aAtmId]))
        for aAtmId in doneAtmIds:
            print("Atom %s is set"%aAtmId)



        # New parts to use linear programming to solve the problem
        # The application may be more general (replacing the last one)
        if len(aromAtmMap.keys()) > 0:
            nextAtoms = []


            self.initRoundOfSetBO(tAtoms,  tBonds, allAtmBondingMap,
                                  aromAtmMap, doneAtmIds, nextAtoms,
                                  aHMap, aCurVaMap)


            self.setCurVal(tAtoms, tBonds, doneAtmIds, aCurVaMap, allAtmBondingMap)


            """
            for aAtm in tAtoms:
                aAtmId = aAtm["_chem_comp_atom.atom_id"]
                if not aAtmId in doneAtmIds:
                    print("Current val for %s is %d "%(aAtmId, aCurVaMap[aAtmId]))

            print("Atoms left : ")
            for aAtm in tAtoms:
                aAtmId = aAtm["_chem_comp_atom.atom_id"]
                if not aAtmId in doneAtmIds:
                    print(aAtmId)
            """
            self.setCPRound(tAtoms, tBonds, allAtmBondingMap, aElemMap,
                           doneAtmIds, aCurVaMap)

            print("Done atoms after 1st round :")
            print(doneAtmIds)

            self.setOtherRound(tAtoms, tBonds, allAtmBondingMap, aElemMap,
                           doneAtmIds, aCurVaMap, adefVaMap, aChargeMap)
            print("Done atoms after 2nd round :")
            print(doneAtmIds)
            #self.setNumHOnN(tAtoms, allAtmBondingMap, aElemMap,aHMap,
            #                adefVaMap)

            #print("Done atoms after 3rd round :")
            #print(doneAtmIds)
            if len(doneAtmIds) > 0:
                print("The following atoms are assigned the bonds ")
                for aAtm1 in doneAtmIds:
                    print("Atom %s has following bonds:"%aAtm1)
                    for aAtm2 in allAtmBondingMap[aAtm1].keys():
                       print("Bonding with atom %s, bond-order %d"%(aAtm2, allAtmBondingMap[aAtm1][aAtm2][0]))
                    if aAtm1 in aHMap.keys():
                        print("Plus it connect %d H atoms"%aHMap[aAtm1])

        print("Before add H atoms, number of atoms is ", len(tAtoms))
        idxH = 0
        for aAtmId in aHMap.keys():
            for i in range(int(aHMap[aAtmId])):
                self.addOneHAtomAndRelatedBond2(tAtoms, tBonds, aIdAtmMap[aAtmId][0],  idxH)
                idxH+=1

        print("After add H atoms, number of atoms is ", len(tAtoms))

        for aB in tBonds:
            print("Bond-order between atoms %s and %s is %s "
                  %(aB["_chem_comp_bond.atom_id_1"],aB["_chem_comp_bond.atom_id_2"],
                    aB["_chem_comp_bond.value_order"]))
        sys.exit()

    def setInitNumH(self, tAtm, tAllAtmBondingMap, tAromAtmMap,tHMap,
                    tDefVaMap, tCurVaMap):

        aSumBO = 0

        aAtm1  = tAtm["_chem_comp_atom.atom_id"]
        if not  aAtm1 in tAromAtmMap.keys():
            for aAtm2 in tAllAtmBondingMap[aAtm1].keys():
                aSumBO += tAllAtmBondingMap[aAtm1][aAtm2][0]
            aNumH = tDefVaMap[aAtm1] -float(tAtm["_chem_comp_atom.charge"])- aSumBO
            aNumH = aNumH
            if aNumH < 0:
                print("Bug, sum of bond-order around atom %s is wrong!"%tAtm)
            elif aNumH > 0:
                tHMap[aAtm1] = aNumH
        elif tAtm["_chem_comp_atom.type_symbol"] == "C":
            nB  = len(tAllAtmBondingMap[aAtm1].keys())
            nBA = len(tAromAtmMap[aAtm1])
            if nB==2 and nBA==2:
                tHMap[aAtm1] = 1
        if aAtm1 in tHMap.keys():
            tCurVaMap[aAtm1] -=tHMap[aAtm1]

    def setNumHOnN(self, tAtoms, tAllAtmBondingMap, tElemMap,tHMap,
                    tDefVaMap):


        for aAtm in tAtoms:
            aSumBO = 0

            aAtm1  = aAtm["_chem_comp_atom.atom_id"]
            if tElemMap[aAtm1] =="N" :
                #print("For ", aAtm1)
                for aAtm2 in tAllAtmBondingMap[aAtm1].keys():
                    #print("Bonding with ", aAtm2)
                    aSumBO += tAllAtmBondingMap[aAtm1][aAtm2][0]
                    #print(aSumBO)
                #print("aSumBo=", aSumBO)
                #print("tDefVaMap[aAtm1]=",tDefVaMap[aAtm1])
                aNumH = tDefVaMap[aAtm1] -float(aAtm["_chem_comp_atom.charge"])- aSumBO
                if aNumH < 0:
                    print("Bug, sum of bond-order around atom %s is wrong!"%aAtm)
                elif aNumH > 0:
                    tHMap[aAtm1] = aNumH


    def initRoundOfSetBO(self, tAtoms, tBonds, tAllAtmBondingMap, tAromAtmMap,
                         tDoneAtoms, tNextAtoms, tHMap, tCurVaMap):


        for aAtm in tAtoms:
            aAtm1Id  = aAtm["_chem_comp_atom.atom_id"]
            if not aAtm1Id in tDoneAtoms:
                if not aAtm1Id in tAromAtmMap.keys():
                    aSumBo = 0
                    print("Atm1", aAtm1Id )
                    for aAtm2Id in tAllAtmBondingMap[aAtm1Id].keys():
                        print("aAtm2 ", aAtm2Id)
                        print(tAllAtmBondingMap[aAtm1Id][aAtm2Id][0])
                        aSumBo+=tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]
                    print("Cur Val ", tCurVaMap[aAtm1Id])
                    aDeff = tCurVaMap[aAtm1Id] - aSumBo
                    print("aDeff ", aDeff)
                    if aAtm1Id in tHMap.keys():
                        aDeff -= tHMap[aAtm1Id]
                    print("again, aDeff ", aDeff)
                    if aDeff ==0:
                        if not aAtm1Id in tDoneAtoms:
                            for aAtm2Id in tAllAtmBondingMap[aAtm1Id].keys():
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                            tDoneAtoms.append(aAtm1Id)

        print(tDoneAtoms)
        if len(tDoneAtoms) > 0:
            print("The following atoms are assigned the bonds ")
            for aAtm1 in tDoneAtoms:
                print("Atom %s has following bonds:"%aAtm1)
                for aAtm2 in tAllAtmBondingMap[aAtm1].keys():
                   print("Bonding with atom %s, bond-order %d"%(aAtm2, tAllAtmBondingMap[aAtm1][aAtm2][0]))
                if aAtm1 in tHMap.keys():
                    print("Plus it connect %d H atoms"%tHMap[aAtm1])



        """
                elif aDeff==1:
                    aBOr = "SINGLE"
                    tBonds[tAllAtmBondingMap[aAtm1Id][aAtm2Id][0][1]]["_chem_comp_bond.value_order"]
                         = aBor
                    self.setABond(aSetReBA[0], tBonds, aBOr, tNAtms, tDBo)

         """

    def setCurVal(self, tAtoms, tBonds, tDoneAtmIds, tCurVaMap, tAllAtmBondingMap):

        print("Before :")
        for aAtm in tCurVaMap.keys():
            print ("Current val for atom ", aAtm, "  is ", tCurVaMap[aAtm])

        idxB = 0
        for aB in tBonds:
            if "_chem_comp_bond.value_order" in aB.keys():
                atm1Id = aB["_chem_comp_bond.atom_id_1"].strip()
                atm2Id = aB["_chem_comp_bond.atom_id_2"].strip()
                aBO    = aB["_chem_comp_bond.value_order"].strip().upper()
                if len(aBO) > 4:
                    aBO = aBO[:4]
                if aBO.upper().find("A") !=-1:
                    aBON = 1
                    aB["_chem_comp_bond.value_order"]="SING"
                else:
                    aBON = BondOrderS2N(aBO)

                tAllAtmBondingMap[atm1Id][atm2Id][0] = aBON
                tAllAtmBondingMap[atm2Id][atm1Id][0] = aBON
                #if not atm1Id in tDoneAtmIds:
                #    if atm2Id in tDoneAtmIds:

                tCurVaMap[atm1Id]-=aBON
                #if not atm2Id in tDoneAtmIds:
                #    if atm1Id in tDoneAtmIds:
                tCurVaMap[atm2Id]-=aBON
            idxB+=1
        print("After :")
        for aAtm in tCurVaMap.keys():
            print ("Current val for atom ", aAtm, "  is ", tCurVaMap[aAtm])


    def setCPRound(self, tAtoms, tBonds, tAllAtmBondingMap, tElemMap,
                   tDoneAtmIds, tCurVaMap):

        for aAtm in tAtoms:
            aAtm1Id = aAtm["_chem_comp_atom.atom_id"]
            if tElemMap[aAtm1Id]=="C" or tElemMap[aAtm1Id]=="P":
                if not aAtm1Id in tDoneAtmIds and tCurVaMap[aAtm1Id] !=0:
                    lDone = False
                    for aAtm2Id in tAllAtmBondingMap[aAtm1Id].keys() :
                        if tElemMap[aAtm2Id]=="C" and not aAtm2Id in tDoneAtmIds:
                            idxB = tAllAtmBondingMap[aAtm1Id][aAtm2Id][1]
                            if not lDone:

                                tBonds[idxB]["_chem_comp_bond.value_order"] = "DOUB"
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]=2
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][0]=2
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                                #print("Bond order between %s and %s is %d "%(aAtm1Id, aAtm2Id, tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]))
                                print("Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                        tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                        tBonds[idxB]["_chem_comp_bond.value_order"]))
                                lDone = True
                            else:
                                tBonds[idxB]["_chem_comp_bond.value_order"] = "SING"
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]=1
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                                #print("Bond order between %s and %s is %d "%(aAtm1Id, aAtm2Id, tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]))
                                print("Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                        tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                        tBonds[idxB]["_chem_comp_bond.value_order"]))

                            tCurVaMap[aAtm2Id]-=1
                            for aAtm3Id in tAllAtmBondingMap[aAtm2Id].keys():
                                if not aAtm3Id in  tDoneAtmIds and aAtm3Id !=aAtm1Id:
                                    idxB2 = tAllAtmBondingMap[aAtm2Id][aAtm3Id][1]
                                    tBonds[idxB2]["_chem_comp_bond.value_order"] = "SING"
                                    tAllAtmBondingMap[aAtm2Id][aAtm3Id][0]=1
                                    tAllAtmBondingMap[aAtm3Id][aAtm2Id][0]=1
                                    tAllAtmBondingMap[aAtm2Id][aAtm3Id][2]=True
                                    tAllAtmBondingMap[aAtm3Id][aAtm2Id][2]=True
                                    #print("Bond order between %s and %s is %d "%(aAtm1Id, aAtm2Id, tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]))
                                    print("Bond order between %s  and %s is %s"%(tBonds[idxB2]["_chem_comp_bond.atom_id_1"],
                                                                tBonds[idxB2]["_chem_comp_bond.atom_id_2"],
                                                                tBonds[idxB2]["_chem_comp_bond.value_order"]))

                            tDoneAtmIds.append(aAtm2Id)
                            break
                    tDoneAtmIds.append(aAtm1Id)
                    tCurVaMap[aAtm1Id]-=1

        # Further clean those current val ==0 atoms
        for aAtm in tAtoms:
            aAtm1Id = aAtm["_chem_comp_atom.atom_id"]

            if tElemMap[aAtm1Id]=="C" or tElemMap[aAtm1Id]=="P":
                if not aAtm1Id in tDoneAtmIds and tCurVaMap[aAtm1Id] ==0:
                    for aAtm2Id in tAllAtmBondingMap[aAtm1Id].keys():
                        idxB = tAllAtmBondingMap[aAtm1Id][aAtm2Id][1]
                        if tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]==1:
                            tBonds[idxB]["_chem_comp_bond.value_order"] = "SING"
                        elif tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]==2:
                            tBonds[idxB]["_chem_comp_bond.value_order"] = "DOUB"
                        elif tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]==3:
                            tBonds[idxB]["_chem_comp_bond.value_order"] = "TRIP"
                        print("Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                tBonds[idxB]["_chem_comp_bond.value_order"]))
                        tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                        tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                    tDoneAtmIds.append(aAtm1Id)

    def setOtherRound(self, tAtoms, tBonds, tAllAtmBondingMap, tElemMap,
                   tDoneAtmIds, tCurVaMap, tDefVaMap, tChargeMap):

        for aAtm in tAtoms:
            aAtm1Id = aAtm["_chem_comp_atom.atom_id"]
            if not aAtm1Id in tDoneAtmIds:
                print("For ", aAtm1Id)
                for aAtm2Id in tAllAtmBondingMap[aAtm1Id].keys():

                    if not tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]:
                        print("NB ", aAtm2Id)
                        idxB = tAllAtmBondingMap[aAtm1Id][aAtm2Id][1]
                        aSumBo = 0

                        for aAtm22Id in tAllAtmBondingMap[aAtm1Id]:
                            print("For NB2 ", aAtm22Id)
                            if aAtm22Id!=aAtm2Id:
                                aSumBo+= tAllAtmBondingMap[aAtm1Id][aAtm22Id][0]
                        print("sum =", aSumBo)
                        tDeff = tDefVaMap[aAtm1Id]-aSumBo-tChargeMap[aAtm1Id]
                        print("tDeff=", tDeff)
                        if tDeff > 1:
                            aSumBo2=0
                            for aAtm23Id in tAllAtmBondingMap[aAtm2Id]:
                                if aAtm23Id!=aAtm1Id:
                                    print("For NB23 ", aAtm23Id)
                                    aSumBo2+= tAllAtmBondingMap[aAtm2Id][aAtm23Id][0]
                            print("sum23 =", aSumBo2)
                            tDeff23 = tDefVaMap[aAtm2Id]-aSumBo2-tChargeMap[aAtm2Id]
                            if tDeff23 > 1:
                                print("doub")
                                tBonds[idxB]["_chem_comp_bond.value_order"] = "DOUB"
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]=2
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][0]=2
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                                print("2 Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                    tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                    tBonds[idxB]["_chem_comp_bond.value_order"]))
                            else:
                                print("sing")
                                tBonds[idxB]["_chem_comp_bond.value_order"] = "SING"
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]=1
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][0]=1
                                tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                                tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                                print("1 Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                        tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                        tBonds[idxB]["_chem_comp_bond.value_order"]))

                        else:
                            print("sing")
                            tBonds[idxB]["_chem_comp_bond.value_order"] = "SING"
                            tAllAtmBondingMap[aAtm1Id][aAtm2Id][0]=1
                            tAllAtmBondingMap[aAtm2Id][aAtm1Id][0]=1
                            tAllAtmBondingMap[aAtm1Id][aAtm2Id][2]=True
                            tAllAtmBondingMap[aAtm2Id][aAtm1Id][2]=True
                            print("1 Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                    tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                    tBonds[idxB]["_chem_comp_bond.value_order"]))
                        tCurVaMap[aAtm2Id]-=1
                        for aAtm3Id in tAllAtmBondingMap[aAtm2Id].keys():
                            if not aAtm3Id in  tDoneAtmIds and aAtm3Id !=aAtm1Id:
                                idxB2 = tAllAtmBondingMap[aAtm2Id][aAtm3Id][1]
                                tBonds[idxB2]["_chem_comp_bond.value_order"] = "SING"
                                tAllAtmBondingMap[aAtm2Id][aAtm3Id][0]=1
                                tAllAtmBondingMap[aAtm3Id][aAtm2Id][0]=1
                                tAllAtmBondingMap[aAtm2Id][aAtm3Id][2]=True
                                tAllAtmBondingMap[aAtm3Id][aAtm2Id][2]=True
                                print("3 Bond order between %s  and %s is %s"%(tBonds[idxB]["_chem_comp_bond.atom_id_1"],
                                                        tBonds[idxB]["_chem_comp_bond.atom_id_2"],
                                                        tBonds[idxB]["_chem_comp_bond.value_order"]))
                        if not aAtm2Id in tDoneAtmIds:
                            tDoneAtmIds.append(aAtm2Id)

                if not aAtm1Id in  tDoneAtmIds:
                    tDoneAtmIds.append(aAtm1Id)
                tCurVaMap[aAtm1Id]-=1


    def getAndSetStartAtm(self, tAromAtoms, tAtmBondingMap, tAtoms, tBonds,
                          tDAtms, tNAtms, tDBo, tdefVaMap, tCurVaMap, tHMap):

        set1Atms = []
        setUAtms = []

        for aAtm in tAromAtoms.keys():
            for aAtm2 in tAtmBondingMap[aAtm]:
                if tAtmBondingMap[aAtm][aAtm2][0] ==1:
                    set1Atms.append(aAtm)
                else:
                    setUAtms.append(aAtm)


        if len(set1Atms) > 0:
            for aAtm in set1Atms:
                #print("Bond-order 1 atom : ", aAtm)
                print("\n***************")
                print("deal with ", aAtm)
                print("***************")
                self.setAtmBondOrder(aAtm, tAtmBondingMap,tAtoms, tBonds,
                                          tDAtms, tNAtms, tDBo, tdefVaMap,
                                          tCurVaMap)
                if not aAtm in tDAtms:
                    tDAtms.append(aAtm)

    def setAtmBondOrder(self, tAtm, tAtmBondingMap,tAtoms, tBonds,
                        tDAtms, tNAtms, tDBo, tdefVaMap, tCurVaMap):


        aSetBA   =[]
        aSetReBA =[]


        for aAtm2 in tAtmBondingMap[tAtm].keys():
            print("aAtm2 ", aAtm2)
            print(tAtmBondingMap[tAtm][aAtm2][0])
            if tAtmBondingMap[tAtm][aAtm2][0] != 1.5\
                and tAtmBondingMap[aAtm2][tAtm][0] !=1.5:
                tCurVaMap[tAtm]-=tAtmBondingMap[tAtm][aAtm2][0]
                print("cur val ",tCurVaMap[tAtm])
                aSetBA.append([aAtm2,tAtmBondingMap[tAtm][aAtm2]])
            else:
                aSetReBA.append([aAtm2, tAtmBondingMap[tAtm][aAtm2]])


        reVal = tCurVaMap[tAtm]
        nReA  = len(aSetReBA)
        print(reVal)
        print(nReA)
        if reVal ==2.0:
            if nReA ==1:
                aBOr = "SINGLE"
                self.setABond(aSetReBA[0], tBonds, aBOr, tNAtms, tDBo)
                tCurVaMap[tAtm]          -=1.0
                tCurVaMap[aSetReBA[0][0]]-=1.0
                tAtmBondingMap[tAtm][aSetReBA[0][0]][0]=1.0

            if nReA ==2:
                aBOr = "SINGLE"
                for aPair in aSetReBA:
                    self.setABond(aPair, tBonds, aBOr, tNAtms, tDBo)
                    print("Here %s  and %f"%(aPair[0], tCurVaMap[aPair[0]]))
                    tCurVaMap[tAtm]     -=1.0
                    tCurVaMap[aPair[0]] -=1.0
                    tAtmBondingMap[tAtm][aPair[0]][0]=1.0

        elif reVal ==3.0:
            if nReA ==1:
                aBOr = "DOUBLE"
                self.setABond(aSetReBA[0], tBonds, aBOr, tNAtms, tDBo)
                tCurVaMap[tAtm]          -=2.0
                tCurVaMap[aSetReBA[0][0]]-=2.0
                tAtmBondingMap[tAtm][aSetReBA[0][0]][0]=2.0
            elif nReA ==2:
                aBOr = "SINGLE"
                self.setABond(aSetReBA[0], tBonds, aBOr, tNAtms, tDBo)
                tCurVaMap[tAtm]        -=1.0
                tCurVaMap[aSetReBA[0][0]] -=1.0
                tAtmBondingMap[tAtm][aSetReBA[0][0]][0]=1.0
                aBOr = "DOUBLE"
                self.setABond(aSetReBA[1], tBonds, aBOr, tNAtms, tDBo)
                tCurVaMap[tAtm]    -=2.0
                tCurVaMap[aSetReBA[1][0]]-=2.0
                tAtmBondingMap[tAtm][aSetReBA[0][0]][0]=2.0

        print("Bond orders around Atom %s has been set, they are:  "%tAtm)
        for aAtm2 in tAtmBondingMap[tAtm]:
            if not aAtm2 in tNAtms and not aAtm2 in tDAtms:
                tNAtms.append(aAtm2)
            idxB = tAtmBondingMap[tAtm][aAtm2][1]
            print("Bond between atom %s and %s is %s"%(tAtm, aAtm2,
                      tBonds[idxB]["_chem_comp_bond.value_order"]))
            #print("curV for ", aAtm2, " is ", tCurVaMap[aAtm2])

    def setABond(self, tPair, tBonds, tBOr, tNAtms, tDBo):

        #print("Bond order should be  ", tBOr)
        tBonds[tPair[1][1]]["_chem_comp_bond.value_order"] = tBOr
        tDBo.append(tPair[1][1])
        tNAtms.append(tPair[0])
        #print("Bond order is ", tBonds[tPair[1][1]]["_chem_comp_bond.value_order"])
    
    
    # A set of methods initially for Carbone molecules
    
    def addAtomSeri(self, tAtoms, tBonds):
    
        nIdx =0
        atIdMap = {}
        for aAt in tAtoms:
            aAt["_chem_comp_atom.atom_serial_number"] = nIdx
            aId = aAt["_chem_comp_atom.atom_id"] 
            atIdMap[aId] = aAt["_chem_comp_atom.atom_serial_number"]
            nIdx+=1
    
        for aB in tBonds:
            aId1 = aB["_chem_comp_bond.atom_id_1"]
            aId2 = aB["_chem_comp_bond.atom_id_2"]
            if aId1 in atIdMap:
                aB["_chem_comp_bond.atom_serial_number_1"] = atIdMap[aId1]
            else:
                print("Bug. can not find atom ", aId1)
        
            if aId2 in atIdMap:
                aB["_chem_comp_bond.atom_serial_number_2"] = atIdMap[aId2]
            else:
                print("Bug. can not find atom ", aId2)
    
    def setAtomConn(self, tAtoms, tBonds):
    
        for aAt in tAtoms:
            aAt["_chem_comp_atom.atom_conn"] = []
    
        if "_chem_comp_bond.atom_serial_number_1" in tBonds[0]:
            for aB in tBonds:
                               
                idxA1 = int(aB["_chem_comp_bond.atom_serial_number_1"])
                idxA2 = int(aB["_chem_comp_bond.atom_serial_number_2"])
                tAtoms[idxA1]["_chem_comp_atom.atom_conn"].append(idxA2) 
                tAtoms[idxA2]["_chem_comp_atom.atom_conn"].append(idxA1) 
    
    def checkCarbonMol(self, tAtoms, tBonds):
        
        aRet = False
        
        aElmMap = {}
        aIdxMap  = {}
        aIdx =0
        for aAt in tAtoms:
            aElmMap[aAt["_chem_comp_atom.atom_id"]] = aAt["_chem_comp_atom.type_symbol"]
            aIdxMap[aAt["_chem_comp_atom.atom_id"]] = aIdx 
            aIdx +=1
            aAt["connB"] = []
        
        for aBo in tBonds:
            id1 = aBo["_chem_comp_bond.atom_id_1"]
            idx1 = aIdxMap[id1]
            id2 = aBo["_chem_comp_bond.atom_id_2"]
            idx2 = aIdxMap[id2]
            if  aElmMap[id2] == "B":
                 tAtoms[idx1]["connB"].append(idx2)
            if  aElmMap[id1] == "B":
                 tAtoms[idx2]["connB"].append(idx1)
       
        for aAt in tAtoms:
            if len(aAt["connB"]) >=4:
                aRet = True
                break
        
        
        return aRet
    
    def getNewConns(self, tAtoms, tBonds, tNewBonds, tNewConns, tBrokenConns):
        
        notBroken = ["H", "CL", "F", "I", "D", "BR", "AT"]
        idxB =0
        for aB in tBonds:
            idxA1 = int(aB["_chem_comp_bond.atom_serial_number_1"])
            id1   = tAtoms[idxA1]["_chem_comp_atom.atom_id"]
            elem1 = tAtoms[idxA1]["_chem_comp_atom.type_symbol"].upper()
            conn1 = len(tAtoms[idxA1]["_chem_comp_atom.atom_conn"])
            
            val1  = 0
            if elem1.upper() in self.maxVal:
                val1 = self.maxVal[elem1]
            else:
                # metal elements
                val1 = self.maxVal["OTHERS"]
                
            #print("============")
            #print("id1=", id1)
            #print("elem1=", elem1)
            #print("conn1=", conn1)
            #print("val1=", val1)
            #print("len(tAtoms[connB])=", len(tAtoms[idxA1]["connB"]))    
            idxA2 = int(aB["_chem_comp_bond.atom_serial_number_2"])
            id2   = tAtoms[idxA2]["_chem_comp_atom.atom_id"]
            elem2 = tAtoms[idxA2]["_chem_comp_atom.type_symbol"].upper()
            conn2 = len(tAtoms[idxA2]["_chem_comp_atom.atom_conn"])
            val2  = 0
            if elem2 in self.maxVal:
                val2 = self.maxVal[elem2]
            else:
                # metal elements
                val2 = self.maxVal["OTHERS"]
            #print("id2=", id2)
            #print("elem2=", elem2)
            #print("conn2=", conn2)
            #print("val2=", val2)
            #print("len(tAtoms[connB])=", len(tAtoms[idxA2]["connB"]))
            # break bonds if need
            if len(tAtoms[idxA1]["connB"]) ==0 and  len(tAtoms[idxA2]["connB"])==0:
                tNewBonds.append(aB) 
            elif len(tAtoms[idxA1]["connB"])==0 and not elem1 in self.maxVal:
                if  conn2 <= val2 :
                     tNewBonds.append(aB)
            elif not tAtoms[idxA2]["connB"] and not elem2 in self.maxVal:
                if conn1 <= val1 :
                    tNewBonds.append(aB)
            elif conn1 > val1 :
                if conn2 > val2:
                    tNewBonds.append(aB)
                elif conn2==val2 and len(tAtoms[idxA2]["connB"]) >=2:
                    tNewBonds.append(aB)
                else:
                    if not elem2 in notBroken:
                        print("1. Bond between atom %s of elem %s and atom %s of elem %s is broken "%(id1, elem1, id2, elem2))
                        tBrokenConns.append([idxA1, idxA2, idxB])
                    else:
                        tNewBonds.append(aB)
            elif conn2 > val2 :    #and elem2 in tVals:
                if conn1 > val1:
                    tNewBonds.append(aB)
                elif conn1==val1 and len(tAtoms[idxA1]["connB"]) >=2:
                    
                    tNewBonds.append(aB)
                else:
                    if not elem1 in notBroken:
                        print("2. Bond between atom %s of elem %s and atom %s of elem %s is broken "%(id1, elem1, id2, elem2))
                        tBrokenConns.append([idxA1, idxA2, idxB])
                    else:
                        tNewBonds.append(aB)
            else:   
                 tNewBonds.append(aB)
            idxB+=1
        
        # Get new set of conn (id not idx) for new mols  
        
        for aB in tNewBonds:
            idxA1 = int(aB["_chem_comp_bond.atom_serial_number_1"])
            id1   =  tAtoms[idxA1]["_chem_comp_atom.atom_id"]
            #print("a1 ", idxA1)
            #print("id1 ", tAtoms[idxA1]["_chem_comp_atom.atom_id"])
            idxA2 = int(aB["_chem_comp_bond.atom_serial_number_2"])
            id2   = tAtoms[idxA2]["_chem_comp_atom.atom_id"]
            #print("a2 ", idxA2)
            #print("id2 ", tAtoms[idxA2]["_chem_comp_atom.atom_id"])
            #print(tAtoms[idxA1]["_chem_comp_atom.atom_conn"])
            if not id1 in tNewConns :
                tNewConns[id1] = []
            tNewConns[id1].append(id2) 
            if not id2 in tNewConns:
                tNewConns[id2] = []
            tNewConns[id2].append(id1)
        
    def addHAtomsToMols(self,  tAtoms, tMol, tBrokenBs):
        # H atoms and bonds here are replace the atoms at broken bonds
        # tAtoms and tBonds are original atoms > tMols["atoms"]
        
        tMol["tmpHAtoms"] = []
        
        aIdxPool = []
        for aAt in tMol["atoms"]:
            aIdxPool.append(aAt["_chem_comp_atom.atom_serial_number"])
        
        
        # This method add tempo Hs to Temp mols from CB mols
        for aPair in tBrokenBs:
            if not aPair[0] in aIdxPool and  aPair[1] in aIdxPool:
                #print("Bond between %s and %s is broken"%(tAtoms[aPair[0]]["_chem_comp_atom.atom_id"], tAtoms[aPair[1]]["_chem_comp_atom.atom_id"]))
                self.addOneHAtom(aPair[0], aPair[1], tAtoms, tMol, self.tmpSeriNum)
                self.tmpSeriNum +=1 
            elif not aPair[1] in aIdxPool and  aPair[0] in aIdxPool: 
                self.addOneHAtom(aPair[1], aPair[0], tAtoms, tMol, self.tmpSeriNum)
                self.tmpSeriNum +=1
            
    def addOneHAtom(self, tIdx1, tIdx2, tAtoms, tMol, tNumAts):
        
        aHAt ={}
        
        aHAt["_chem_comp_atom.atom_id"] = "H" + str(tNumAts+1)
        
        aHAt["_chem_comp_atom.type_symbol"]   = "H"
        aHAt["_chem_comp_atom.model_Cartn_x"] = tAtoms[tIdx1]["_chem_comp_atom.model_Cartn_x"] 
        aHAt["_chem_comp_atom.model_Cartn_y"] = tAtoms[tIdx1]["_chem_comp_atom.model_Cartn_y"]
        aHAt["_chem_comp_atom.model_Cartn_z"] = tAtoms[tIdx1]["_chem_comp_atom.model_Cartn_z"]   
        
        aHB = {}
        aHB["_chem_comp_bond.atom_id_1"] = aHAt["_chem_comp_atom.atom_id"] 
        aHB["_chem_comp_bond.atom_id_2"] = tAtoms[tIdx2]["_chem_comp_atom.atom_id"] 
        aHB["_chem_comp_bond.value_order"] = "SINGLE"
        
        tMol["tmpHAtoms"].append(aHAt["_chem_comp_atom.atom_id"])
        # Add H atom to atoms and bonds of the mol
        tMol["atoms"].append(aHAt)
        tMol["bonds"].append(aHB)
        
        self.tmpHIdMap[aHAt["_chem_comp_atom.atom_id"]] = tAtoms[tIdx1]["_chem_comp_atom.atom_id"]
        
    def getNewMolsFromConns(self, tAtoms, tBonds, tBrokenBs, tOutDir, tRoot):
        
        aSetMols = []
        self.tmpSeriNum = len(tAtoms)
    
        idxM = 1
        for aKey in self.newLinkedGroups:
            
            aMol = {}
            aMol["atoms"] = []
            aMol["bonds"] = []
            aMol["fileIdx"] = tRoot + "_" + str(idxM) 
            idxM +=1
            idxAtms = []
            for aB in tBonds:
                idx1 = aB["_chem_comp_bond.atom_serial_number_1"]
                idx2 = aB["_chem_comp_bond.atom_serial_number_2"]
                if idx1 in self.newLinkedGroups[aKey] and idx2 in self.newLinkedGroups[aKey]:
                    aMol["bonds"].append(aB)
                    if not idx1 in idxAtms:
                        idxAtms.append(idx1)
                        aMol["atoms"].append(tAtoms[idx1])
                        
                    if not idx2 in idxAtms:
                        idxAtms.append(idx2)
                        aMol["atoms"].append(tAtoms[idx2])
            #print(idxAtms)
            for aIdx in self.newLinkedGroups[aKey]:
                if not aIdx in idxAtms:
                    idxAtms.append(aIdx)
                    aMol["atoms"].append(tAtoms[aIdx])
                    
            
            if len(aMol["atoms"]) > 0 and len(aMol["bonds"]) > 0:
                if self.checkCarbonMol(aMol["atoms"], aMol["bonds"]):
                    aMol["isCBMol"] = True
                    self.reSetAtomConnsForCBMol(tAtoms, aMol)
                else:
                    aMol["isCBMol"] = False
                
                if not aMol["isCBMol"] : 
                    self.addHAtomsToMols(tAtoms, aMol, tBrokenBs)
                
                aSetMols.append(aMol) 
                
                
                aCifName = os.path.join(tOutDir, aMol["fileIdx"] +".cif")
                try:
                    aCif = open(aCifName, "w")
                except IOError :
                    print("%s can not be opened for writing "%aCifName)
                    sys.exit(1)
                else:
                    outputOneMolInACif(aCif, tRoot, aMol)
                    aCif.close()
                    
         
        return  aSetMols
    
    def reSetAtomConnsForCBMol(self, tAllAtoms, tMol):
        
        idMap = {}
        idxA =0
        for aAt in tMol["atoms"]:
            idMap[aAt["_chem_comp_atom.atom_id"]] = idxA
            idxA+=1 
        
        for aAt in tMol["atoms"]:
            aTmpConn = []
            for aIdxNAt in aAt["_chem_comp_atom.atom_conn"]:
                aId = tAllAtoms[aIdxNAt]["_chem_comp_atom.atom_id"]
                if aId in idMap:
                    aTmpConn.append(idMap[aId])
            
            aAt["_chem_comp_atom.atom_conn"] = []
            #print("Now atom ", aAt["_chem_comp_atom.atom_id"], " conn : ")
            for aIdxNew in aTmpConn:
                aAt["_chem_comp_atom.atom_conn"].append(aIdxNew)
                #print(tMol["atoms"][aIdxNew]["_chem_comp_atom.atom_id"])
        
    def reSetAtomConnNonCBMol(self, tAtoms, tBonds, tExcHs):
    
        idMap = {}
        idxA =0
        for aAt in tAtoms:
            aAt["_chem_comp_atom.atom_conn"] = []
            idMap[aAt["_chem_comp_atom.atom_id"]] = idxA
            idxA+=1 
            
        for aB in tBonds:
            id1 = aB["_chem_comp_bond.atom_id_1"]
            id2 = aB["_chem_comp_bond.atom_id_2"]
            if not id1 in tExcHs and not id2 in tExcHs:
                if id1 in idMap and id2 in idMap:
                    idx1 = idMap[id1]
                    idx2 = idMap[id2]
                    tAtoms[idx1]["_chem_comp_atom.atom_conn"].append(idx2) 
                    tAtoms[idx2]["_chem_comp_atom.atom_conn"].append(idx1)            
            
    
    def splitOneCBMol(self, tMol):
        
        tMol["storedHs"] = []
        tMol["storedHBs"] = []
        tMol["remainAtoms"] = []
        tMol["remainBonds"] = []
        tMol["remainAngs"] = []
        
        hIdList = []
    
        
        for aAt in tMol["atoms"]:
            if aAt["_chem_comp_atom.type_symbol"] == "H":
                tMol["storedHs"].append(aAt)
                hIdList.append(aAt["_chem_comp_atom.atom_id"])
            else:
                tMol["remainAtoms"].append(aAt)
        #print("storedHs=", tMol["storedHs"])        
        # Modified atom idx in new mole but keep the record of old atom idx
        idMap = {}
        
        for aIdx, aAt in enumerate(tMol["remainAtoms"]):
            if "_chem_comp_atom.atom_serial_number" in aAt:
                aAt["chem_comp_atom.atom_tmp_serial_number"]   = aAt["_chem_comp_atom.atom_serial_number"]
            aAt["chem_comp_atom.atom_serial_number"]           = aIdx
            aId = aAt["_chem_comp_atom.atom_id"] 
            idMap[aId] = aIdx 
            
        for aB in tMol["bonds"]:
            id1 = aB["_chem_comp_bond.atom_id_1"]
            id2 = aB["_chem_comp_bond.atom_id_2"]
            if "_chem_comp_bond.atom_serial_number_1" in aB:
                aB["_chem_comp_bond.atom_tmp_serial_number_1"] = aB["_chem_comp_bond.atom_serial_number_1"]      
            if "_chem_comp_bond.atom_serial_number_2" in aB:
                aB["_chem_comp_bond.atom_tmp_serial_number_2"] = aB["_chem_comp_bond.atom_serial_number_2"]
            if id1 in hIdList or id2 in hIdList:
                tMol["storedHBs"].append(aB)
            else:
                aB["_chem_comp_bond.atom_serial_number_1"] = idMap[id1]
                aB["_chem_comp_bond.atom_serial_number_2"] = idMap[id2]
                tMol["remainBonds"].append(aB)
             
    def getClassIds(self, tMol):
        
        elmMap = {}
        elmConnMap = {}
        
        
        for aNIdx in tMol["CB_Graph"].nodes:    
            aElm  = tMol["CB_Graph"].nodes[aNIdx]["_chem_comp_atom.type_symbol"]
            if not aElm in elmMap:
                elmMap[aElm] = []
            elmMap[aElm].append(tMol["CB_Graph"].nodes[aNIdx]["_chem_comp_atom.atom_id"])
            nConn = len(tMol["CB_Graph"].nodes[aNIdx]["_chem_comp_atom.atom_conn"]) 
            if not aElm in elmConnMap:
                elmConnMap[aElm] = {}
            if not nConn in elmConnMap[aElm]:
                elmConnMap[aElm][nConn] = []
            elmConnMap[aElm][nConn].append(tMol["CB_Graph"].nodes[aNIdx]["_chem_comp_atom.atom_id"])
        
        
        
        classId0 = ""
        classId1 = ""
        
        if "B" in elmMap:
            classId0 +="[B" + str(len(elmMap["B"])) + "]"
        if "C" in elmMap:
            classId0 +="[C" + str(len(elmMap["C"])) + "]"
        for aElm in sorted(elmMap.keys()):
            if aElm !="B" and aElm !="C":
                classId0 +="[" + aElm  + str(len(elmMap[aElm])) + "]"   
             
        if "B" in elmConnMap:
            for aC in sorted(elmConnMap["B"].keys()):
                classId1 += "[B" + str(aC) + "_" + str(len(elmConnMap["B"][aC])) + "]"
        if "C" in elmConnMap:
            for aC in sorted(elmConnMap["C"].keys()):
                classId1 += "[C" + str(aC) + "_" + str(len(elmConnMap["C"][aC])) + "]"
        for aElm in elmConnMap:
            if aElm !="B" and aElm != "C":
                for aC in sorted(elmConnMap[aElm].keys()):
                    classId1 += "[" + aElm + str(aC) + "_" + str(len(elmConnMap[aElm][aC])) + "]"
               
        return classId0, classId1
    
    def serchGraphDB(self, tMol):
        
        aCId0, aCId1 = self.getClassIds(tMol)
        print("The graph cid0 is %s and cid1 is %s "%(aCId0, aCId1))
        if aCId0 in self.CBGrapList:
            if aCId1 in self.CBGrapList[aCId0]:
                for aFN in self.CBGrapList[aCId0][aCId1]:
                    print("Match DB graph %s with target %s"%(aFN, tMol["CB_Graph"].graph["name"]))
                    getMatched, matchedIdxs, aDB_G=self.graphMatchTwoMols(aFN, tMol)
                    if getMatched :
                        print("%s and %s is isomorphic"%(tMol["CB_Graph"].graph["name"], aFN))
                        self.recoverAllAtomsAndBonds(tMol, aDB_G,  matchedIdxs)
                        break
                    #else:       
                    #    print("%s and %s is not isomorphic"%(tMol["CB_Graph"].graph["name"], aFN))
    
    def recoverAllAtomsAndBonds(self, tMol, tDB_G, tIdxMap):
        
        self.reSetCoordsByIsomMatch(tMol, tDB_G, tIdxMap)
       
        print("Before add, number of atoms", len(tMol["remainAtoms"]))
        print("Number of H detached ", len(tMol["storedHs"]))
        print("Number of bonds is ", len(tMol["remainBonds"]))
        print("Number of Broken H bonds", len(tMol["storedHBs"]))
        self.setHCoordsAndBondsOnCBBall(tMol)
        print("After adding: ")
        print("Number of atoms", len(tMol["remainAtoms"]))
        print("Number of bonds is ", len(tMol["remainBonds"]))
        self.setBondLengForCBBalls(tMol)
    
        
    def setHCoordsAndBondsOnCBBall(self, tMol):
        
        aHIdMap = {}
        idxH = 0
        for aH in tMol["storedHs"]:
            aId = aH["_chem_comp_atom.atom_id"]
            aHIdMap[aId] = idxH
            idxH+=1
        
        aNonHMap = {}
        idxNonH =0
        for aAt in tMol["remainAtoms"]:
            aId = aAt["_chem_comp_atom.atom_id"]
            aNonHMap[aId] =  idxNonH
            idxNonH +=1  
        
        for aB in tMol["storedHBs"]:
            idxH    =-1
            idxNonH =-1
            id1 = aB["_chem_comp_bond.atom_id_1"]
            
            id2  = aB["_chem_comp_bond.atom_id_2"]
            if id1 in aHIdMap and id2 in aNonHMap:
                idxH = aHIdMap[id1]
                idxNonH = aNonHMap[id2]
            elif id2 in aHIdMap and id1 in aNonHMap:
                idxH = aHIdMap[id2]
                idxNonH =aNonHMap[id1]
             
            if idxH !=-1 and idxNonH !=-1:
                #print("NonH atom ", tMol["remainAtoms"][idxNonH]["_chem_comp_atom.atom_id"])
                #print("H atom ", tMol["storedHs"][idxH]["_chem_comp_atom.atom_id"])
                aLengH = 1.1
                self.setOneHAtomCoordsAndOneBondAndAllAngs(idxH, idxNonH, tMol, aLengH, aB)
                #print("atom ", tMol["remainAtoms"][idxNonH]["_chem_comp_atom.atom_id"], "conns : ")
                #for aIdx in tMol["remainAtoms"][idxNonH]["_chem_comp_atom.atom_conn"]:
                #    print("atom ", tMol["remainAtoms"][aIdx]["_chem_comp_atom.atom_id"])
                
    def setBondLengForCBBalls(self, tMol):
        
        connMap  = {}
        idMap = {}
        idxAt =0
        for aAt in tMol["remainAtoms"]:
            
            idMap[aAt["_chem_comp_atom.atom_id"]] = idxAt
            idxAt+=1
        
        
        for aB in tMol["remainBonds"]:
            id1 = aB["_chem_comp_bond.atom_id_1"]
            id2 = aB["_chem_comp_bond.atom_id_2"]
            
            idx1 = idMap[id1]
            x1 = float(tMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_x"])
            y1 = float(tMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_y"])
            z1 = float(tMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_z"])
            pos1 = [x1, y1, z1]
            
            idx2 = idMap[id2]
            x2 = float(tMol["remainAtoms"][idx2]["_chem_comp_atom.model_Cartn_x"])
            y2 = float(tMol["remainAtoms"][idx2]["_chem_comp_atom.model_Cartn_y"])
            z2 = float(tMol["remainAtoms"][idx2]["_chem_comp_atom.model_Cartn_z"])
            pos2 = [x2, y2, z2]
            aL = getLengBetween2Pos(pos1, pos2)
            aB["_chem_comp_bond.value_dist"] = aL
            aB["_chem_comp_bond.value_dist_nucleus"] = aL
            if not "_chem_comp_bond.value_dist_esd" in aB:
                aB["_chem_comp_bond.value_dist_esd"] = 0.01
                aB["_chem_comp_bond.value_dist_nucleus_esd"] = 0.01
            
            #print("Bond leng between %s and %s is %6.4f"%(id1, id2, aB["_chem_comp_bond.value_dist"]))
            #print("Its esd is ", aB["_chem_comp_bond.value_dist_esd"])
        
        
         
    def setOneHAtomCoordsAndOneBondAndAllAngs(self, tIdxH, tIdxNonH, tMol, tLeng, tBond):
        
        if len(tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.atom_conn"]) > 4:
            
            x1 = float(tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.model_Cartn_x"])
            y1 = float(tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.model_Cartn_y"])
            z1 = float(tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.model_Cartn_z"])
            aPosNonH = [x1, y1, z1]
            
            aSetPos = []
            
            aSetCombP = []
            for aIdx in tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.atom_conn"]:
                aId = tMol["remainAtoms"][aIdx]["_chem_comp_atom.atom_id"]
                x = float(tMol["remainAtoms"][aIdx]["_chem_comp_atom.model_Cartn_x"])
                y = float(tMol["remainAtoms"][aIdx]["_chem_comp_atom.model_Cartn_y"])
                z = float(tMol["remainAtoms"][aIdx]["_chem_comp_atom.model_Cartn_z"])
                aPos = [x,y,z]
                aSetPos.append(aPos)
                aSetCombP.append([aPos, aId])
                
            aDirV= getAUnitDirVec(aPosNonH, aSetPos)

            tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_x"] = x1 + tLeng*aDirV[0]
            tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_y"] = y1+ tLeng*aDirV[1]
            tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_z"] = z1 + tLeng*aDirV[2]
            aPosH = [tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_x"],tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_y"],
                    tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_z"]]
            tMol["remainAtoms"].append(tMol["storedHs"][tIdxH])
            tBond["_chem_comp_bond.value_dis"] = tLeng
            tMol["remainBonds"].append(tBond)
            
    
            for aPair in aSetCombP:
                aAng = {}
                aAng["_chem_comp_angle.comp_id"]     = "LIG"
                aAng["_chem_comp_angle.atom_id_1"]   =  tMol["storedHs"][tIdxH]["_chem_comp_atom.atom_id"]
                aAng["_chem_comp_angle.atom_id_2"]   =  tMol["remainAtoms"][tIdxNonH]["_chem_comp_atom.atom_id"]
                aAng["_chem_comp_angle.atom_id_3"]   =  aPair[1]
                aAng["_chem_comp_angle.value_angle"] =  118.0        # getAAngFrom3Ps(aPosNonH, aPosH, aPair[0])
                #print("an angle between %s-%s-%s is %6.4f" %(aAng["_chem_comp_angle.atom_id_1"], aAng["_chem_comp_angle.atom_id_2"], 
                #                                            aAng["_chem_comp_angle.atom_id_3"], aAng["_chem_comp_angle.value_angle"]))
                aAng["_chem_comp_angle.value_angle_esd"] = 3.00
                tMol["remainAngs"].append(aAng)       
            
            
    def addOneHAtomCoordsAndOneBondAndAllAngs(self, tAt, tMol, tLeng, tCurNumAtoms, tExcH):
        
        
        if len(tAt["_chem_comp_atom.atom_conn"]) > 4:
            x1 = float(tAt["_chem_comp_atom.model_Cartn_x"])
            y1 = float(tAt["_chem_comp_atom.model_Cartn_y"])
            z1 = float(tAt["_chem_comp_atom.model_Cartn_z"])
            aPosNonH = [x1, y1, z1]
            
            aSetPos = []
            
            aSetCombP = []
            for aIdx in tAt["_chem_comp_atom.atom_conn"]:
                aId = tMol["atoms"][aIdx]["_chem_comp_atom.atom_id"]
                x = float(tMol["atoms"][aIdx]["_chem_comp_atom.model_Cartn_x"])
                y = float(tMol["atoms"][aIdx]["_chem_comp_atom.model_Cartn_y"])
                z = float(tMol["atoms"][aIdx]["_chem_comp_atom.model_Cartn_z"])
                aPos = [x,y,z]
                aSetPos.append(aPos)
                aSetCombP.append([aPos, aId])
                
            aDirV= getAUnitDirVec(aPosNonH, aSetPos)
             
            aHAtom = {}
            aHAtom["_chem_comp_atom.comp_id"] = tAt["_chem_comp_atom.comp_id"]
            aHAtom["_chem_comp_atom.atom_id"] = "H" + str(tCurNumAtoms)
            if aHAtom["_chem_comp_atom.atom_id"] in tExcH:
                aHAtom["_chem_comp_atom.atom_id"] =aHAtom["_chem_comp_atom.atom_id"]+"a"
            tCurNumAtoms+=1
            aHAtom["_chem_comp_atom.alt_atom_id"] = aHAtom["_chem_comp_atom.atom_id"]
            aHAtom["_chem_comp_atom.type_symbol"] = "H"
            aHAtom["_chem_comp_atom.charge"]      = "0"
            aHAtom["_chem_comp_atom.pdbx_align"]  = "1"
            aHAtom["_chem_comp_atom.pdbx_aromatic_flag"] = "N"
            aHAtom["_chem_comp_atom.pdbx_leaving_atom_flag"] = "N"
            aHAtom["_chem_comp_atom.pdbx_stereo_config"]  = "N"
            
            aHAtom["_chem_comp_atom.model_Cartn_x"] = x1 + tLeng*aDirV[0]
            aHAtom["_chem_comp_atom.model_Cartn_y"] = y1+ tLeng*aDirV[1]
            aHAtom["_chem_comp_atom.model_Cartn_z"] = z1 + tLeng*aDirV[2]
            aHAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"] = aHAtom["_chem_comp_atom.model_Cartn_x"]
            aHAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"] = aHAtom["_chem_comp_atom.model_Cartn_y"]
            aHAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"] = aHAtom["_chem_comp_atom.model_Cartn_z"]
            aHAtom["_chem_comp_atom.x"] = aHAtom["_chem_comp_atom.model_Cartn_x"]
            aHAtom["_chem_comp_atom.y"] = aHAtom["_chem_comp_atom.model_Cartn_y"]
            aHAtom["_chem_comp_atom.z"] = aHAtom["_chem_comp_atom.model_Cartn_z"]
            #aPosH = [tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_x"],tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_y"],
            #        tMol["storedHs"][tIdxH]["_chem_comp_atom.model_Cartn_z"]]
            tMol["remainAtoms"].append(aHAtom)
            print("H atom ", aHAtom["_chem_comp_atom.alt_atom_id"], " is added ")
            aBond = {}
            aBond["_chem_comp_bond.comp_id"]   = tAt["_chem_comp_atom.comp_id"]
            aBond["_chem_comp_bond.atom_id_1"] = tAt["_chem_comp_atom.atom_id"]
            aBond["_chem_comp_bond.atom_id_2"] = aHAtom["_chem_comp_atom.atom_id"]
            aBond["_chem_comp_bond.value_order"] = "SING" 
            aBond["_chem_comp_bond.value_dist"] = tLeng
            aBond["_chem_comp_bond.value_dist_esd"] = 0.01
            aBond["_chem_comp_bond.pdbx_aromatic_flag"] = "N"
            aBond["_chem_comp_bond.value_dist_nucleus"] = tLeng
            aBond["_chem_comp_bond.value_dist_nucleus_esd"] = 0.01
            
            tMol["remainBonds"].append(aBond)
            
    
            for aPair in aSetCombP:
                aAng = {}
                aAng["_chem_comp_angle.comp_id"]     =  tAt["_chem_comp_atom.comp_id"]
                aAng["_chem_comp_angle.atom_id_1"]   =  aHAtom["_chem_comp_atom.atom_id"]
                aAng["_chem_comp_angle.atom_id_2"]   =  tAt["_chem_comp_atom.atom_id"]
                aAng["_chem_comp_angle.atom_id_3"]   =  aPair[1]
                aAng["_chem_comp_angle.value_angle"] =  118.0        # getAAngFrom3Ps(aPosNonH, aPosH, aPair[0])
                #print("an angle between %s-%s-%s is %6.4f" %(aAng["_chem_comp_angle.atom_id_1"], aAng["_chem_comp_angle.atom_id_2"], 
                #                                            aAng["_chem_comp_angle.atom_id_3"], aAng["_chem_comp_angle.value_angle"]))
                aAng["_chem_comp_angle.value_angle_esd"] = 3.00
                tMol["remainAngs"].append(aAng)                 
        
        return tCurNumAtoms
        
    def reSetCoordsByIsomMatch(self, tTargetMol, tG_DB, tMapIdxs):
        
        """
        idMap = {}
        idxA =0
        for aAt in tTargetMol["atoms"]:
            idMap[aAt["_chem_comp_atom.atom_id"]] = idxA
            idxA+=1
        """    
        for idx1, idx2 in tMapIdxs:
            #id1 = tTargetMol["remainAtoms"][idx1]["_chem_comp_atom.atom_id"]
            #oriIdx = idMap[id1]
            tTargetMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_x"] = tG_DB.nodes[idx2]["_chem_comp_atom.model_Cartn_x"]
            tTargetMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_y"] = tG_DB.nodes[idx2]["_chem_comp_atom.model_Cartn_y"]
            tTargetMol["remainAtoms"][idx1]["_chem_comp_atom.model_Cartn_z"] = tG_DB.nodes[idx2]["_chem_comp_atom.model_Cartn_z"]
            
        
                    
    def graphMatchTwoMols(self, tCif, tMol):
        
        lM = False
        aIdxMap = {}
        aCifName = os.path.join(self.acedrgTables, "CarboroneSamples", tCif + ".cif")
        if os.path.isfile(aCifName):
            tmpMols = []
            fromACrysToMolCifsGemmi(aCifName, tCif, tmpMols)
            self.addAtomSeri(tmpMols[0]["atoms"], tmpMols[0]["bonds"])
            aGFromGB = aMolToAGraph(tmpMols[0]["atoms"], tmpMols[0]["bonds"], tCif)
            lM, aIdxMap =  MatchTwoGraphs(tMol["CB_Graph"], aGFromGB)
            
            
        return lM, aIdxMap , aGFromGB
        
                
    def processOneCBMol(self, tMol):
        
        self.splitOneCBMol(tMol)
        self.setAtomConn(tMol["remainAtoms"], tMol["remainBonds"])
        tMol["CB_Graph"] = aMolToAGraph(tMol["remainAtoms"], tMol["remainBonds"], tMol["fileIdx"])
        #print(tMol["CB_Graph"])
        self.serchGraphDB(tMol)
        
        
        
        
    def getCB_Graph_DB(self, tTabDir):
        
        self.acedrgTables = tTabDir
        
        sampDir = os.path.join(self.acedrgTables, "CarboroneSamples")
        if not os.path.isdir(sampDir):
            print("Install bug : %s does not exist"%sampDir)
            sys.exit()
        graphListN = os.path.join(self.acedrgTables, "CarboroneSamples.list")
        if not os.path.isfile(graphListN):
            print("Install bug : %s does not exist"%sampDir)
            sys.exit()
        graphListF = open(graphListN, "r")
        allLs      = graphListF.readlines()
        graphListF.close()
        for aL in allLs:
            strs = aL.strip().split()
            if len(strs)==3:
                strs[0] = strs[0].strip().upper()
                strs[1] = strs[1].strip().upper()
                strs[2] = strs[2].strip()
                if not strs[0] in self.CBGrapList:
                    self.CBGrapList[strs[0]] = {}
                if not  strs[1] in self.CBGrapList[strs[0]]:
                    self.CBGrapList[strs[0]][strs[1]] = []
                self.CBGrapList[strs[0]][strs[1]].append(strs[2]) 
         
                    
    def getNewMols(self, tAtoms, tBonds, tOutDir, tRoot):
        
        newBonds = []
        newConns = {}
        
        tmpBrokenBonds = []
        
        self.addAtomSeri(tAtoms, tBonds)
        self.setAtomConn(tAtoms, tBonds)
        self.getNewConns(tAtoms, tBonds, newBonds, newConns, tmpBrokenBonds)
        self.newLinkedGroups = getLinkedGroups2(tAtoms, newConns)
        
        print(self.newLinkedGroups)
        self.aSetNewMols = self.getNewMolsFromConns(tAtoms, newBonds,tmpBrokenBonds,  tOutDir, tRoot)
        
        #for aMol in self.aSetNewMols:
        #    if aMol["isCBMol"]:
        #        for aAt in aMol["atoms"]:
        #            print("2:Now atom ", aAt["_chem_comp_atom.atom_id"], " conn : ")
        #            for aIdxNew in aAt["_chem_comp_atom.atom_conn"]:
        #                print(aMol["atoms"][aIdxNew]["_chem_comp_atom.atom_id"])
                        
        if len(tmpBrokenBonds) > 0:
            self.tmpBrokenBonds = []
            for aIdxB in tmpBrokenBonds:
                id1 = tBonds[aIdxB[2]]["_chem_comp_bond.atom_id_1"]
                id2 = tBonds[aIdxB[2]]["_chem_comp_bond.atom_id_2"]
                self.tmpBrokenBondIds[id1] = id2
                self.tmpBrokenBondIds[id2] = id1
                self.tmpBrokenBonds.append(aIdxB)
        
        print("Numbber of Mols is ", len(self.aSetNewMols))
        for aMol in self.aSetNewMols:
            if aMol["isCBMol"]:
                print(aMol["fileIdx"], " is a CB molecule")
        
        
        """    
        print("Those bonds are temporarily broken:")
        for aIdxB in tmpBrokenBonds:
            id1 = tAtoms[aIdxB[0]]["_chem_comp_atom.atom_id"]
            id2 = tAtoms[aIdxB[1]]["_chem_comp_atom.atom_id"]
            
            print("Bond between %s and %s is broken"%(id1, id2))
        """
        
    def assembleNewMods(self, tAllAtoms, tAllBonds, tNonCBMols, tCBMols, tOutDir, tRoot, tComboMol):
        
        atomIdMap = {}
        print("Assember the mols ")
        """
        for aMol in self.aSetNewMols:
            idxA=0
            for aAt in aMol["atoms"] :
                atomIdMap[aAt["_chem_comp_atom.atom_id"]] = [idxM, idxA]
                idxA+=1
            idxM+=1 
        """
        
        tmpHs = []
        for aMol in self.aSetNewMols:
            if "tmpHAtoms" in aMol:
                for aHId in aMol["tmpHAtoms"]:
                    tmpHs.append(aHId)
        
        atomIdMapNonCB = {}
        idxM =0
        for aMol in tNonCBMols:
            idxA=0
            for aAt in aMol["atoms"] :
                atomIdMapNonCB[aAt["_chem_comp_atom.atom_id"]] = [idxM, idxA]
                idxA+=1
            idxM+=1
            self.reSetAtomConnNonCBMol(aMol["atoms"], aMol["bonds"], tmpHs)
        
        atomIdMapCB = {}
        idxM =0
        for aMol in tCBMols:
            idxA=0
            for aAt in aMol["atoms"] :
                #print("Now 6.1 atom ", aAt["_chem_comp_atom.atom_id"], " conn :")
                aAt["_chem_comp_atom.atom_conn"] = []
                for aN in aAt["_chem_comp_atom.tmp_atom_conn"]:
                    aAt["_chem_comp_atom.atom_conn"].append(aN)
                    #print(aMol["atoms"][aN]["_chem_comp_atom.atom_id"])
                atomIdMapCB[aAt["_chem_comp_atom.atom_id"]] = [idxM, idxA]
                idxA+=1
            idxM+=1
        
        extraBonds = []
        extraAngs  = []
        extraTors  = []
        
        #print(self.tmpBrokenBonds)
        
        for aIdxB in self.tmpBrokenBonds:
            id1 = tAllAtoms[aIdxB[0]]["_chem_comp_atom.atom_id"]
            id2 = tAllAtoms[aIdxB[1]]["_chem_comp_atom.atom_id"]
            #print("Broken bond id1 ", id1, " id2 ", id2)
            
            if id1 in atomIdMapNonCB and id2 in atomIdMapCB:
                self.TRNonCBMols(tNonCBMols, tCBMols, atomIdMapNonCB[id1], atomIdMapCB[id2], aIdxB[2], tAllBonds, extraBonds, extraAngs)
                self.RotateNonCBMols(tNonCBMols, tCBMols, atomIdMapNonCB[id1], atomIdMapCB[id2])
            elif id2 in atomIdMapNonCB and id1 in atomIdMapCB:
                self.TRNonCBMols(tNonCBMols, tCBMols, atomIdMapNonCB[id2], atomIdMapCB[id1], aIdxB[2], tAllBonds, extraBonds, extraAngs)
                self.RotateNonCBMols(tNonCBMols, tCBMols, atomIdMapNonCB[id2], atomIdMapCB[id1])
        
        tComboMol["fileIdx"]  = "comboMol_LIG"
        tComboMol["atoms"]    = []
        tComboMol["bonds"]    = []
        tComboMol["angles"]   = []
        tComboMol["torsions"] = []
        tComboMol["chirs"]    = []
        tComboMol["planes"]   = {}
        tComboMol["rings"]    = {}
        numAllAtoms = len(tAllAtoms)
        
        for aCBMol in tCBMols:
            self.reSetAtomConnInMol(aCBMol)
            self.checkAndRepairHAtoms(aCBMol, numAllAtoms, tmpHs)
        
        self.setAComboMol(tNonCBMols, tCBMols,  extraBonds, extraAngs, extraTors, tmpHs, tComboMol)
        
            
        #self.outIniCifForAssemMol(extraBonds, extraAngs)
    def reSetAtomConnInMol(self, tMol):
        
        aIdMap = {}
        idxA =0
        for aAt in tMol["atoms"]:
            aId = aAt["_chem_comp_atom.atom_id"]
            aAt["_chem_comp_atom.atom_serial_number"] = idxA
            aIdMap[aId] = idxA 
            idxA+=1
            aAt["_chem_comp_atom.atom_conn"] = []
        
        for aB in tMol["bonds"]:
            id1 = aB["_chem_comp_bond.atom_id_1"]
            idx1 = aIdMap[id1]
            id2 = aB["_chem_comp_bond.atom_id_2"]
            idx2 = aIdMap[id2]
            
            tMol["atoms"][idx1]["_chem_comp_atom.atom_conn"].append(idx2)
            tMol["atoms"][idx2]["_chem_comp_atom.atom_conn"].append(idx1)
            
        
    def checkAndRepairHAtoms(self, tCBMol, tCurNumAtoms, tExcHs):
        
        aLeng = 1.1
        for aAt in tCBMol["atoms"]:
            if aAt["_chem_comp_atom.type_symbol"] =="B" or aAt["_chem_comp_atom.type_symbol"] =="C":
                #print("atom ", aAt["_chem_comp_atom.atom_id"], " conn ", len(aAt["_chem_comp_atom.atom_conn"]))
                
                if len(aAt["_chem_comp_atom.atom_conn"])==5 and not aAt["_chem_comp_atom.atom_id"] in self.tmpBrokenBondIds:
                    #print("it needs a H")
                    #print("Its index is ", aAt["_chem_comp_atom.atom_serial_number"])
                    #print("Confirm: its id is ", tCBMol["atoms"][aAt["_chem_comp_atom.atom_serial_number"]]["_chem_comp_atom.atom_id"])
                    if not "angles" in tCBMol:
                        tCBMol["angles"] = []
                    tCurNumAtoms=self.addOneHAtomCoordsAndOneBondAndAllAngs(aAt, tCBMol, aLeng, tCurNumAtoms, tExcHs)
                    
        
        
    def setAComboMol(self, tNonCBMols, tCBMols, tExtraBonds, tExtraAngs, tExtraTors, tExcHs, tComboMol):
        
        aTorIdxs = {}
        aTorIdxs["const"]    = 0
        aTorIdxs["sp1_sp1"]  = 1
        aTorIdxs["sp1_sp2"]  = 1
        aTorIdxs["sp1_sp3"]  = 1
        aTorIdxs["sp2_sp2"]  = 1
        aTorIdxs["sp2_sp3"]  = 1
        aTorIdxs["sp3_sp3"]  = 1
        
        numChirs = 1
        numPlans = 1
        numRings = 1
        
        #print(self.tmpBrokenBondIds)
        
        allMols = []
        for aMol in tNonCBMols:
            allMols.append(aMol)
        for aMol in   tCBMols:
            allMols.append(aMol)
            
        for aMol in tNonCBMols:
            
            for aAt in aMol["atoms"]:
                
                if not aAt["_chem_comp_atom.atom_id"] in tExcHs:
                    tComboMol["atoms"].append(aAt)
                    
            for aBo in aMol["bonds"]:
                id1 = aBo["_chem_comp_bond.atom_id_1"]
                id2 = aBo["_chem_comp_bond.atom_id_2"]
                if not id1 in tExcHs and not id2 in tExcHs:
                    tComboMol["bonds"].append(aBo)
                   
            if "angles" in aMol:
                for aAng in aMol["angles"]:
                    id1 = aAng["_chem_comp_angle.atom_id_1"]
                    id2 = aAng["_chem_comp_angle.atom_id_2"]
                    id3 = aAng["_chem_comp_angle.atom_id_3"]
                    if id1 in self.tmpHIdMap:
                        aAng["_chem_comp_angle.atom_id_1"] = self.tmpHIdMap[id1]
                        tComboMol["angles"].append(aAng)
                    elif id3 in self.tmpHIdMap:
                        aAng["_chem_comp_angle.atom_id_3"] = self.tmpHIdMap[id3]
                        tComboMol["angles"].append(aAng)
                    else:
                        tComboMol["angles"].append(aAng)
                        
            if "torsions" in aMol:
                for aTor in aMol["torsions"]:
                    id1 = aTor["_chem_comp_tor.atom_id_1"]
                    id2 = aTor["_chem_comp_tor.atom_id_2"]
                    id3 = aTor["_chem_comp_tor.atom_id_3"]
                    id4 = aTor["_chem_comp_tor.atom_id_4"]
                    if id1 in self.tmpHIdMap :
                        aTor["_chem_comp_tor.atom_id_1"] = self.tmpHIdMap[id1]
                        aTor["_chem_comp_tor.id"]=self.reAdjustTorId(aTor["_chem_comp_tor.id"], aTorIdxs)
                        tComboMol["torsions"].append(aTor)
                    elif id4 in self.tmpHIdMap :
                        aTor["_chem_comp_tor.atom_id_4"] = self.tmpHIdMap[id4]
                        aTor["_chem_comp_tor.id"]=self.reAdjustTorId(aTor["_chem_comp_tor.id"], aTorIdxs)
                        tComboMol["torsions"].append(aTor)
                    else:
                        aTor["_chem_comp_tor.id"]=self.reAdjustTorId(aTor["_chem_comp_tor.id"], aTorIdxs)
                        tComboMol["torsions"].append(aTor)
            if "chirs" in aMol:
                for aCh in aMol["chirs"]:
                    id1 = aCh["_chem_comp_chir.atom_id_1"]
                    id2 = aCh["_chem_comp_chir.atom_id_2"]
                    id3 = aCh["_chem_comp_chir.atom_id_3"]
                    if not id1 in tExcHs and not id2 in tExcHs and not id3 in tExcHs:
                        aCh["_chem_comp_chir.id"] = str(numChirs)
                        numChirs +=1
                        tComboMol["chirs"].append(aCh)
            if "planes" in aMol:
                for aPId in aMol["planes"]:
                    aNewId = "plan-" + str(numPlans)
                    numPlans +=1
                    tComboMol["planes"][aNewId] = []
                    for aAt in aMol["planes"][aPId]:
                        aAt["_chem_comp_plane_atom.plane_id"] = aNewId 
                        if aAt["_chem_comp_plane_atom.atom_id"] in self.tmpHIdMap:
                            aAt["_chem_comp_plane_atom.atom_id"] = self.tmpHIdMap[aAt["_chem_comp_plane_atom.atom_id"]] 
                        tComboMol["planes"][aNewId].append(aAt)
            if "rings" in aMol:
                for aRIdx in aMol["rings"]:
                    aNewRId = "ring-" + str(numRings)
                    numRings +=1 
                    tComboMol["rings"][aNewRId] = []
                    for aAt in aMol["rings"][aRIdx]:
                        aAt["_chem_comp_ring_atom.ring_serial_number"] = aNewRId 
                        tComboMol["rings"][aNewRId].append(aAt)
                        
        for aMol in tCBMols:
            
            for aAt in aMol["remainAtoms"]:
                
                if not aAt["_chem_comp_atom.atom_id"] in tExcHs:
                # Cancel the charges at CB ball
                    aAt["_chem_comp_atom.charge"] = "0.0"
                    tComboMol["atoms"].append(aAt)
                    
            for aBo in aMol["remainBonds"]:
                id1 = aBo["_chem_comp_bond.atom_id_1"]
                id2 = aBo["_chem_comp_bond.atom_id_2"]
                if not id1 in tExcHs and not id2 in tExcHs:
                    tComboMol["bonds"].append(aBo)
            
            if "remainAngs" in aMol:
                for aAng in aMol["remainAngs"]:
                    tComboMol["angles"].append(aAng)
            
                    
        for aBo in tExtraBonds:
            tComboMol["bonds"].append(aBo)
        for aAng in tExtraAngs:
            tComboMol["angles"].append(aAng)
        for aTor in tExtraTors:
            self.reAdjustTorId(aTor["_chem_comp_tor.id"], aTorIdxs)
            #print("aTor id is ", aTor["_chem_comp_tor.id"])
            tComboMol["torsions"].append(aTor)
        
        
        
    def reAdjustTorId(self, tId, tIdIdx):
        
        
        if tId.find("const") !=-1:
            tId = "const_" + str(tIdIdx["const"])
            tIdIdx["const"] +=1
        elif tId.find("sp1_sp1") !=-1:
            tId = "sp1_sp1_" + str(tIdIdx["sp1_sp1"])
            tIdIdx["sp1_sp1"] +=1
        elif tId.find("sp1_sp2") !=-1:
            tId = "sp1_sp2_" + str(tIdIdx["sp1_sp2"]) 
            tIdIdx["sp1_sp2"] +=1
        elif tId.find("sp1_sp3") !=-1:
            tId = "sp1_sp3_" + str(tIdIdx["sp1_sp3"]) 
            tIdIdx["sp1_sp3"] +=1
        elif tId.find("sp2_sp2") !=-1:
            tId = "sp2_sp2_" + str(tIdIdx["sp2_sp2"]) 
            tIdIdx["sp2_sp2"] +=1
        elif tId.find("sp2_sp3") !=-1:
            tId = "sp2_sp3_" + str(tIdIdx["sp2_sp3"]) 
            tIdIdx["sp2_sp3"] +=1
        elif tId.find("sp3_sp3") !=-1:
            tId = "sp3_sp3_" + str(tIdIdx["sp3_sp3"]) 
            tIdIdx["sp3_sp3"] +=1
        
        #print(tIdIdx)
        return tId
            
    
    def TRNonCBMols(self, tNonCBMols, tCBMols,tNonCBIdxMA, tCBIdxMA, tIdxB, tAllBonds, tExtraBonds, tExtraAngs):
        
        # 
        
        idxCBM    = -1
        idxNonCBM = -1
        
        
        idxNonCBM = tNonCBIdxMA[0]
        idxNonCBA = tNonCBIdxMA[1] 
        idxCBM    = tCBIdxMA[0]
        idxCBA    = tCBIdxMA[1]
        
        
        
        if idxCBM  !=-1 and idxNonCBM !=-1:
            
            # Tempo for B element
            aLeng = 1.57  
            if tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.type_symbol"] == "C":
                aLeng = 1.51     # Temp for CB 
            
            transV= setOneAtomCoordsOutCB(idxCBA, idxNonCBA, tCBMols[idxCBM]["atoms"], tNonCBMols[idxNonCBM]["atoms"], aLeng)
            #print(self.aSetNewMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.atom_id"])
            #print(self.aSetNewMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_id"])        
            tAllBonds[tIdxB]["_chem_comp_bond.value_dist"] = aLeng
            #print("again")
            #print(tAllBonds[tIdxB]["_chem_comp_bond.atom_id_1"])
            #print(tAllBonds[tIdxB]["_chem_comp_bond.atom_id_2"])
            tExtraBonds.append(tAllBonds[tIdxB])
            
            for aAt in tNonCBMols[idxNonCBM]["atoms"]:
                if aAt["_chem_comp_atom.atom_id"] != tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_id"]:
                    
                    if "_chem_comp_atom.model_Cartn_x" in aAt:
                        oldX = float(aAt["_chem_comp_atom.model_Cartn_x"]) 
                        oldY = float(aAt["_chem_comp_atom.model_Cartn_y"])
                        oldZ = float(aAt["_chem_comp_atom.model_Cartn_z"])
                    elif "_chem_comp_atom.pdbx_model_Cartn_x_ideal" in aAt:
                        oldX= float(aAt["_chem_comp_atom.pdbx_model_Cartn_x_ideal"])
                        oldY= float(aAt["_chem_comp_atom.pdbx_model_Cartn_y_ideal"])
                        oldZ= float(aAt["_chem_comp_atom.pdbx_model_Cartn_z_ideal"])
                    elif "_chem_comp_atom.x" in aAt:
                        oldX = float(aAt["_chem_comp_atom.x"])
                        oldY = float(aAt["_chem_comp_atom.y"])
                        oldZ = float(aAt["_chem_comp_atom.z"])
                    else:
                        print("No coordinates ?. Program at this section needs coordinates")
                        sys.exit(1)
                        
                    newX = oldX + transV[0]
                    newY = oldY + transV[1]
                    newZ = oldZ + transV[2]
                     
                    replaceAtomCoords(aAt, newX,  newY, newZ)
            
        
            
            for a3rdIdx in tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.atom_conn"]:
                if tCBMols[idxCBM]["atoms"][a3rdIdx]["_chem_comp_atom.atom_id"] !=\
                    tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_id"]:
                    aAng = {}
                    aAng["_chem_comp_angle.comp_id"]     = "LIG"
                    aAng["_chem_comp_angle.atom_id_1"]   =  tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_id"]
                    aAng["_chem_comp_angle.atom_id_2"]   =  tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.atom_id"]
                    aAng["_chem_comp_angle.atom_id_3"]   =  tCBMols[idxCBM]["atoms"][a3rdIdx]["_chem_comp_atom.atom_id"]
                    aAng["_chem_comp_angle.value_angle"] =  118.0
                    aAng["_chem_comp_angle.value_angle_esd"] = 3.00
                    tExtraAngs.append(aAng)
                    #print("2: angle between %s--%s--%s is added "%(aAng["_chem_comp_angle.atom_id_1"], aAng["_chem_comp_angle.atom_id_2"],
                    #                                             aAng["_chem_comp_angle.atom_id_3"]))

    def RotateNonCBMols(self, tNonCBMols, tCBMols,tNonCBIdxMA, tCBIdxMA):
        
        idxCBM    = -1
        idxNonCBM = -1


        idxNonCBM = tNonCBIdxMA[0]
        idxNonCBA = tNonCBIdxMA[1] 
        idxCBM    = tCBIdxMA[0]
        idxCBA    = tCBIdxMA[1]

        if idxCBM  !=-1 and idxNonCBM !=-1:
            # x1 = y1  corresponds to idxNonCBA
            # y2 = idxCBA, x2->idxNonCBA2 needs to be found out
            for aAt in tNonCBMols[idxNonCBM]["atoms"]:
                aAt["_chem_comp_atom.model_Cartn_x"] = float(aAt["_chem_comp_atom.model_Cartn_x"])
                aAt["_chem_comp_atom.model_Cartn_y"] = float(aAt["_chem_comp_atom.model_Cartn_y"])
                aAt["_chem_comp_atom.model_Cartn_z"] = float(aAt["_chem_comp_atom.model_Cartn_z"])
            for aAt in tCBMols[idxCBM]["atoms"]:
                aAt["_chem_comp_atom.model_Cartn_x"] = float(aAt["_chem_comp_atom.model_Cartn_x"])
                aAt["_chem_comp_atom.model_Cartn_y"] = float(aAt["_chem_comp_atom.model_Cartn_y"])
                aAt["_chem_comp_atom.model_Cartn_z"] = float(aAt["_chem_comp_atom.model_Cartn_z"])
                
            idxNonCBA2 = -1 
            idY1  = tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_id"]
            print("x1=y1=", idY1)
            idY2 =  tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.atom_id"]
            print("y2=", idY2)
            
            for aIdx in tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_conn"]:
                idX3 = tNonCBMols[idxNonCBM]["atoms"][aIdx]["_chem_comp_atom.atom_id"]
                elm3 = tNonCBMols[idxNonCBM]["atoms"][aIdx]["_chem_comp_atom.type_symbol"]
                if elm3 != "H" and idX3 != idY2 :
                    idxNonCBA2 = aIdx 
                    break
            
            if idxNonCBA2 != -1 :
                
                aPosC = [tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_x"],
                         tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_y"], 
                         tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_z"]]
                aPos1 = [tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.model_Cartn_x"],
                         tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.model_Cartn_y"],
                         tCBMols[idxCBM]["atoms"][idxCBA]["_chem_comp_atom.model_Cartn_z"]]
                aPos2 = [tNonCBMols[idxNonCBM]["atoms"][idxNonCBA2]["_chem_comp_atom.model_Cartn_x"],
                         tNonCBMols[idxNonCBM]["atoms"][idxNonCBA2]["_chem_comp_atom.model_Cartn_y"],
                         tNonCBMols[idxNonCBM]["atoms"][idxNonCBA2]["_chem_comp_atom.model_Cartn_z"]]
                
                aTheta = getAAngFrom3Ps(aPosC, aPos1, aPos2) 
                nConn = len(tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_conn"])
                beta = 120.0
                if  (nConn +1) == 4:
                    beta = 109.4
                elif  (nConn +1) == 2:
                    beta =  180.0
                
                aPhi = -beta + aTheta 
                
                sPhi = np.sin(aPhi * np.pi / 180.)
                cPhi = np.cos(aPhi * np.pi / 180.)
                    
                idX2 = tNonCBMols[idxNonCBM]["atoms"][idxNonCBA2]["_chem_comp_atom.atom_id"]
                print("The current angle betweeen %s--%s--%s is %8.4f"%(idY2, idY1, idX2, aTheta))
                print("Should rotate %8.4f angle "%aPhi)
                print("atom %s conn %d atoms. "%(idY1, len(tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.atom_conn"])))
                vecN = -getPlaneNormal(aPosC, aPos1, aPos2)
                #a1   = getAUnitDirVec2P(aPosC, aPos1)
                #a2   = getAUnitDirVec2P(aPosC, aPos2)
                #print("n_a1=", np.dot(vecN, a1)) 
                #print("n_a2=", np.dot(vecN, a2))
                #print("a1_a2=", np.dot(a2, a1))
                #print("vecN=",vecN)
                #print("Its length ", np.linalg.norm(vecN))
                
                K   = np.array([[0.0, -vecN[2], vecN[1]],
                                [vecN[2],  0.0, -vecN[0]],
                                [-vecN[1], vecN[0], 0]])
                
                #print("K=", K)
                
                KSq = K@K
                
                """
                n1n2   = vecN[0]*vecN[1]
                n1n3   = vecN[0]*vecN[2]
                n2n3   = vecN[1]*vecN[2]
                
                n1n2sq = -vecN[0]*vecN[0] - vecN[1]*vecN[1]
                n1n3sq = -vecN[0]*vecN[0] - vecN[2]*vecN[2]
                n2n3sq = -vecN[1]*vecN[1] - vecN[2]*vecN[2]
            
                KSq = np.array([[n2n3sq, n1n2,  n1n3],
                                [n1n2,  n1n3sq, n2n3],
                                [n1n3,  n2n3, n1n2sq]])
                """
                
                I   = np.identity(3)
                
               
                
                
                R = I + K*sPhi + KSq*(1-cPhi) 
                
                #print("R=", R)
                #print("R_det=", np.linalg.det(R))
                x = tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_x"]
                y = tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_y"]
                z = tNonCBMols[idxNonCBM]["atoms"][idxNonCBA]["_chem_comp_atom.model_Cartn_z"]
                for aAt in tNonCBMols[idxNonCBM]["atoms"]:
                    if aAt["_chem_comp_atom.atom_id"] != idY1:
                        aDef = np.array([[aAt["_chem_comp_atom.model_Cartn_x"]-x],
                                       [aAt["_chem_comp_atom.model_Cartn_y"]-y],
                                       [aAt["_chem_comp_atom.model_Cartn_z"]-z]])
                        aNewCoords =np.matmul(R,aDef) + np.array([[x], [y], [z],])
                        #print("atom : ", aAt["_chem_comp_atom.atom_id"])
                        #print("aNewCoords=", aNewCoords)
                        aAt["_chem_comp_atom.model_Cartn_x"] = aNewCoords[0][0]
                        aAt["_chem_comp_atom.model_Cartn_y"] = aNewCoords[1][0]
                        aAt["_chem_comp_atom.model_Cartn_z"] = aNewCoords[2][0]