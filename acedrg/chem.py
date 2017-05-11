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
 

