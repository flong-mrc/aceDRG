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
from optparse import OptionParser 
import time
import math
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

from . utility  import listComp
from . utility  import listComp2
from . utility  import listCompDes
from . utility  import listCompAcd
from . utility  import BondOrderS2N
from . utility  import splitLineSpa
from . utility  import splitLineSpa2
from . utility  import aLineToAlist

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
        self.orgVal["H"]     = [1]
        self.orgVal["D"]     = [1]
        self.orgVal["F"]     = [1]
        self.orgVal["CL"]    = [1]
        self.orgVal["BR"]    = [1]
        self.orgVal["I"]     = [1]
        self.orgVal["AT"]    = [1]
        self.orgVal["O"]     = [2]
        self.orgVal["S"]     = [2, 4, 6]
        self.orgVal["SE"]    = [2]
        self.orgVal["N"]     = [3]
        self.orgVal["B"]     = [3]
        self.orgVal["C"]     = [4]
        self.orgVal["P"]     = [5]    # should be [5, 3]. Use one at the moment


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
            print("The input ligands/molecules contains metal or other heavier atoms ")
            print("Acedrg currently deals with ligands/molecules with following elements only ")
            print("As, B, Br, C, Cl, F, Ge, H, I, N, O, P, S, Si")
            print("The job finishes without the output cif file")
        if len(allAtomElems) ==0:
            print ("Can not get the element symbols of atoms. Check input file format")
            print("Error : The job stops because of errors")
            organicOnly = False

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
                     
            else:
                 lElem = False

        if len(nonOrgSet)> 0: 
            aRet = False
        if not aRet :
            print("The input ligands/molecules contains metal or other heavier atoms ")
            print("Acedrg currently deals with ligands/molecules with following elements only ")
            print("As, B, Br, C, Cl, F, Ge, H, I, N, O, P, S, Si")
            print ("Your molecule contains atoms of elements of the following type_symbol:")
            for aE in nonOrgSet:
                print(aE)
            print("The job finishes without the output cif file.")
       
        if not lElem :
            print ("Can not get the element symbols of atoms. Check input file format")
            print("Error : The job stops because of errors")
            organicOnly = False

        return aRet

    def isOrganicMol(self, tMol):

        pass 


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
            checkMol=Chem.MolFromSmiles(aSmiStr)
  
            if checkMol:
                for aAtom in checkMol.GetAtoms():
                    aSym = aAtom.GetSymbol().strip()
                    if len(aSym) and not aSym in tAtomElems:
                        tAtomElems.append(aSym)
            else:
                print("Can not generate a molecule from the input SMILES string!")
                print("Check your SMILES or report the bug ")

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
                aLine = "Unknown type %s for the added bond between %s and %s \n"\
                        %(aPair[0]["type"], aPair[0]["atom_id_1"], aPair[0]["atom_id_2"])
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

    def adjustNBForOneAddedAtom(self, tAtom, tNBAtoms, tBonds, tRes, tMod, tDelAtmIds):

        aBool = True
        aLine = ""

        nOrder = 0
        for aB in tBonds:
            nBT = -1
            if "type" in aB:
                nBT = BondOrderS2N(aB["type"])
                print("Bond-type between atom %s and %s is %d "%(aB["atom_id_1"], aB["atom_id_2"], nBT))
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
                print("Current bond order is : ", nOrder)
                if not nOrder in curAllowOrds:
                    # If total valence does not equal to one of allowed valences. Use curAllowOrds[0]
                    # The procedures are:
                    # (1) if the atom has formal charge, dealt with it first
                    # (2) adjust H atom around the atom
                    # (3) error info : tell the user to re-define the added bonds

                    if nOrder == (curAllowOrds[0] - 1):
                        if tAtom["charge"]=="1":
                            # (1a) nOrder = curAllowOrds[0] - 1, cancel the formal charge in atom
                            tAtom["charge"]="0"
                            tMod["changed"]["atoms"].append(tAtom)
                            print("The formal charge at atom %s is reduced to %s"%(tAtom["atom_id"], tAtom["charge"]))
                        else:
                            # (1b) to add H atoms and therefore the valence of the atom.
                            aAtom = {}
                            aAtom["atom_id"]     = "HXX"          # tempo
                            aAtom["type_symbol"] = "H"
                            aAtom["type_energy"] = aAtom["type_symbol"]
                            aAtom["comp_id"]     = tAtom["comp_id"] 
                            tMod["added"]["atoms"].append(aAtom)
                            print("One H atom %s is added into the residue %s "%(aAtom["atom_id"], aAtom["comp_id"]))
                            aBond = {}
                            aBond["atom_id_1"] = aAtom["atom_id"]
                            aBond["atom_id_2"] = tAtom["atom_id"]
                            aBond["type"]      = "SINGLE"
                            tMod["added"]["bonds"].append(aBond)
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
                            print("Here ")
                            sys.exit()
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
               and (len(atomLinks["N"])==2 or len(atomLinks["N"])==3):
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
                    if len(hAtomNameMap["N"]) ==0 or len(hAtomNameMap["N"])> 2 :
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
    
    def addjustAtomsAndBonds2(self, tAtoms, tBonds):
        
        pass 
        
    def addjustAtomsAndBonds(self, tAtoms, tBonds):
        
        # For molecules contain bonds called explicitly aromatic bonds
        
        aCurVaMap   = {}
        adefVaMap   = {}
        aChargeMap  = {}
        aElemMap    = {}
        aHMap       = {}
        
        for aAtm in tAtoms:
            atmId = aAtm["_chem_comp_atom.atom_id"]
            aElemMap [atmId]  = aAtm["_chem_comp_atom.type_symbol"]
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
            print()
            if "_chem_comp_bond.value_order" in aB.keys():
                atm1Id = aB["_chem_comp_bond.atom_id_1"].strip()
                atm2Id = aB["_chem_comp_bond.atom_id_2"].strip()
                aBO    = aB["_chem_comp_bond.value_order"].strip().upper()
                print("Bond between %s and %s"%(atm1Id,atm2Id))
                if len(aBO) > 4:
                    aBO = aBO[:4]
                aBON = BondOrderS2N(aBO)
                
                if not atm1Id in allAtmBondingMap.keys():
                    allAtmBondingMap[atm1Id] = {}
                allAtmBondingMap[atm1Id][atm2Id] = [aBON, idxB] 
                if not atm2Id in allAtmBondingMap.keys():
                    allAtmBondingMap[atm2Id] = {}
                allAtmBondingMap[atm2Id][atm1Id] = [aBON, idxB]
                
                if aBO.find("AROM") !=-1:
                    if not atm1Id in aromAtmMap.keys():
                        aromAtmMap[atm1Id] = [] 
                    if not atm2Id in aromAtmMap.keys():
                        aromAtmMap[atm2Id] = []
                    aromAtmMap[atm1Id].append(atm2Id)
                    aromAtmMap[atm2Id].append(atm1Id)
            idxB+=1
            
        
        for aAtm in tAtoms:
            
            self.setInitNumH(aAtm, allAtmBondingMap, 
                             aromAtmMap, aHMap, adefVaMap)
                                                                    
                
        for aAtmId in aHMap.keys():
            print("Atom %s needs to add %d Hs"%(aAtmId, aHMap[aAtmId]))
        
        
        # New parts to use linear programming to solve the problem
        # The application may be more general (replacing the last one)
        if len(aromAtmMap.keys()) > 0:
            doneAtoms = []
            nextAtoms = []
            doneBonds = []
            # 
                
            self.initRoundOfSetBO(tAtoms,  tBonds, allAtmBondingMap,
                                  aromAtmMap, doneAtoms, nextAtoms,
                                  aHMap, aCurVaMap)
            
            
        
       
            self.getAndSetStartAtm(aromAtmMap, allAtmBondingMap, tAtoms, 
                                   tBonds,doneAtoms, nextAtoms, doneBonds,
                                   adefVaMap, aCurVaMap, aHMap)
            
            
            print("atoms to be done", nextAtoms)
    
            
            
            i=1
            while len(nextAtoms) > 0:
                
                for aAtm in nextAtoms:
                    print("i=",i)
                    if not aAtm in doneAtoms:
                        if aAtm in aromAtmMap.keys() and\
                           not aromAtmMap[aAtm] in doneAtoms:
                               print("\n***************")
                               print("deal with ", aAtm)
                               print("***************")
                               self.setAtmBondOrder(aAtm, allAtmBondingMap,
                                               tAtoms, tBonds, 
                                               doneAtoms, nextAtoms, doneBonds, 
                                               adefVaMap,aCurVaMap)
                               nextAtoms.remove(aAtm)
                               doneAtoms.append(aAtm)
                               sys.exit(1)
                    i+=1    
                               
        
        
        # Now check if H atoms should be added
        #addedHAtms = []
        #for aAtm in tAtoms: 
        
        
    def setInitNumH(self, tAtm, tAllAtmBondingMap, tAromAtmMap,tHMap, tDefVaMap):
        
        aSumBO = 0

        aAtm1  = tAtm["_chem_comp_atom.atom_id"]
        if not  aAtm1 in tAromAtmMap.keys():
            for aAtm2 in tAllAtmBondingMap[aAtm1].keys():
                aSumBO += tAllAtmBondingMap[aAtm1][aAtm2][0]
            aNumH = tDefVaMap[aAtm1] -float(tAtm["_chem_comp_atom.charge"])- aSumBO 
            if aNumH < 0:
                print("Bug, sum of bond-order around atom %s is wrong!"%tAtm)
            elif aNumH > 0:
                tHMap[aAtm1] = aNumH        
        elif tAtm["_chem_comp_atom.type_symbol"] == "C":
            nB  = len(tAllAtmBondingMap[aAtm1].keys())
            nBA = len(tAromAtmMap[aAtm1])
            if nB==2 and nBA==2:
                tHMap[aAtm1] = 1
                

    def initRoundOfSetBO(self, tAtoms, tBonds, tAllAtmBondingMap, tAromAtmMap,
                         tDoneAtoms, tNextAtoms, tHMap, tCurVaMap):
        
        
        for aAtm in tAtoms:
            aAtm1Id  = aAtm["_chem_comp_atom.atom_id"]
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
    
