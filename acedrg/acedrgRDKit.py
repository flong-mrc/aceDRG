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
from builtins import range
from builtins import object
import os,os.path,sys
import platform
import glob,shutil
import re,string
from optparse import OptionParser 
import time
import math
import select
import random

from functools  import cmp_to_key

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

from . chem          import ChemCheck
from . periodicTable import PeriodicTab

from . utility  import listComp
from . utility  import listComp2
from . utility  import listCompDes
from . utility  import listCompAcd
from . utility  import BondOrderS2N

   
class AcedrgRDKit(object):

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
        self.noProtonation           = False
        self.conformerEngMap         = {}
        self.numSelectForRefConfs    = 25
        self.selecConformerIds       = []

        self.hasCCP4Type     = False
   
        self.chemCheck   = ChemCheck()
        self.periodicTab = PeriodicTab()
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
        
        if "useExistCoords" in tProcessParaSet:
            self.useExistCoords = tProcessParaSet["useExistCoords"]
        else:
            self.useExistCoords = False
  
        if "noProtonation" in tProcessParaSet:
            self.noProtonation = tProcessParaSet["noProtonation"]

        if "numRDKitOptmSteps" in tProcessParaSet:
            self.numRDKitOptmSteps = tProcessParaSet["numRDKitOptmSteps"]
        else:
            self.numRDKitOptmSteps = 1000

        if "numInitConformers" in tProcessParaSet:
            self.numInitConformers = tProcessParaSet["numInitConformers"]
        else:
            if self.useExistCoords:
                self.numInitConformers = 1
            else:
                self.numInitConformers = 20

        if "numConformers" in tProcessParaSet:
            self.numConformers = tProcessParaSet["numConformers"]
        else:
            self.numConformers = 1

        if self.numConformers > self.numInitConformers:
            self.numInitConformers = self.numConformers
  
        if self.numInitConformers > self.maxNumInitConformers:
            self.numInitConformers = self.maxNumInitConformers
 
        #self.numSelectForRefConfs = self.numInitConformers/20
        #if self.numSelectForRefConfs < 25:
        #    self.numSelectForRefConfs = 25
        self.numSelectForRefConfs = 20
        

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

        #print("Initial subUnit ", tUniTmp)
        tNewStr =""
        tStrs =[]
        if tUniTmp.find("N@@") !=-1:
           tStrs =tUniTmp.strip().split('N@@')
           tRepStr  = "N@+1" 
        elif tUniTmp.find("N@") !=-1 and tUniTmp.find("N@@")==-1:
           tStrs =tUniTmp.strip().split('N@')
           tRepStr  = "N@@+1" 
        if len(tStrs) >=2:
            tNewStrs = tStrs[0] + tRepStr + tStrs[1] + "(%s)"%self.repSign
            self.smiMod +=tNewStrs
        #print "New repUnit ", tRepStr
        #print "New Unit ", tNewStrs
           
    def initMols(self, tFileType, tFileName, tMonoRoot, tChemCheck, tPH, tNumConf, tMode=0, tNameMap=None, tChargeList=None):
        aTmpMode = 0
        aMolT = None
        aMolName = ""
        print("File type is ", tFileName)
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
                print("File %s does not exist "%tFileName)
                sys.exit()

        elif tFileType == "mol2" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                     aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                     aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                print("Molecule name:  ", aMolName)
                aMolT     = Chem.MolFromMol2File(tFileName)
                aTmpMode  = 1 
                #self.chemCheck.addHs(aMolT)
            else:
                print("File %s does not exist "%tFileName)
                sys.exit()

        elif tFileType =="smi" :
            if os.path.isfile(tFileName):
                # SMILES string in a file
                try:
                    fSmi = open(tFileName, "r")
                except IOError:
                    print("% can not be open for reading "%tFileName)
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
                        print("SMILES string formation error ")
                        sys.exit()
                              
            else:
                # SMILES string in from a commandline  
                aSmiStr = tFileName.strip()
            if len(aSmiStr):
                self.smiOrig = aSmiStr
                print(self.smiOrig)
                #if self.reSetSmi:
                #    self.modifySmiTmp()
                #    self.smiMod = self.smiMod.strip()
                #    aMolT     = Chem.MolFromSmiles(self.smiMod.strip())
                #else:
                aMolT     = Chem.MolFromSmiles(self.smiOrig)
                if aMolT:
                    aMolT.SetProp("SmilesIn", self.smiOrig)
                else:
                    print("No molecule is generated using SMILES str ", self.smiOrig) 
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
                print("File %s does not exist "%tFileName)

        elif tFileType == "pdb" :
            if os.path.isfile(tFileName):
                if platform.system()=="Windows":
                    aMolName = tFileName.strip().split("\\")[-1].strip().split(".")[0]
                else:
                    aMolName = tFileName.strip().split("/")[-1].strip().split(".")[0]
                print("Molecule name:  ", aMolName)
                # Test. H atoms removed as usual. Maybe change that later
                aMolT     = Chem.MolFromPDBFile(tFileName)
                if aMolT:
                    self.removeWater(aMolT)
        if not aMolT:
            print("Molecules can not generated  from file %s ! "%tFileName)
            print("Check your file format ")
            sys.exit()
        else:
            aMolT.SetProp("ResidueName", tMonoRoot)
            #print("file %s aTmpMode %d"%(tFileName,  aTmpMode))
            self.setOtherMolInfo(aMolT, tNumConf, tChemCheck, tPH, tNameMap, tMode, tChargeList, aTmpMode)
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
                if aElem not in dictAtomTypes:
                    dictAtomTypes[aElem] = []
                dictAtomTypes[aElem].append(aA.GetIdx())
            for aElem in list(dictAtomTypes.keys()):
                i = 1
                for aIdx in dictAtomTypes[aElem]:
                    aName = aElem + str(i)
                    dictAtomNames[aIdx] = aName
                    i = i+1

            for aAtom in tMol.GetAtoms():
                print(dictAtomNames[aAtom.GetIdx()])
                aS = str(dictAtomNames[aAtom.GetIdx()])
                aAtom.SetProp("Name", aS) 
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

    def setNamesForAtomsInMol2(self, tMol, tChemCheck, tNameMap, tStage=0):
   
        dictAtomTypes = {}
        dictAtomNames = {}

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
                print("H Atom  of serial number ", aAtom.GetIdx())
                idxB  = -1
                idxE  = -1
                idxC  = -1
                aSetBonds = aAtom.GetBonds()
                if len(aSetBonds)==1: 
                    idxB = aSetBonds[0].GetBeginAtomIdx()
                    idxE = aSetBonds[0].GetEndAtomIdx()
                    print("Bond atom idx 1 ", idxB)
                    print("Bond atom idx 2 ", idxE)
                    if idxH == idxB:
                        idxC = idxE
                    elif idxH==idxE: 
                        idxC = idxB
                    #print "Bonding  to ", " atom ", tMol.GetAtomWithIdx(idxC).GetProp("Name")
                    if idxC !=-1 :
                        nonH_Id = tNameMap["nonH"][idxC] 
                        if nonH_Id not in HConns:
                            HConns[nonH_Id] = []
                        HConns[nonH_Id].append(idxH)  
                    else:
                        print("Can not find the non-H atom that connects to H atom of serial number ", idxH)   
                else:
                    print("One H atom connect more than one atoms, check")
                    sys.exit()
        print("Total number of H atoms is ", nh)
        #for aKey in HConns.keys():
        #    print "Atom ", aKey, " bonds to ", len(HConns[aKey]), " H atoms "


        if len(list(HConns.keys())) !=0:  
            # Check if total numbers of H atoms are different between the original mol and current mol  
            # Names for H atoms in the original file
            origHNames =[]
            for aNonH in list(HConns.keys()):
                if  aNonH in tNameMap["H"]:
                    c3 = len(tNameMap["H"][aNonH])
                    for i in range(c3):
                        origHNames.append(tNameMap["H"][aNonH][i])
            numOrigH = len(origHNames)
            nExtra = numOrigH + 1
            for aK in list(HConns.keys()):
                c1 =  len(HConns[aK])
                if aK in tNameMap["H"]:
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
                        print("H rootName is ", hRootId)         
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
                                print("An added H is named as ", aNewHName) 
                                origHNames.append(aNewHName)
                            else:
                                nExtra = 2
                                tHName = hRootId + str(idxMax+nExtra)
                                while True:
                                    if not tHName in origHNames:
                                        tMol.GetAtomWithIdx(HConns[aK][c2+j]).SetProp("Name",tHName)
                                        origHNames.append(tHName)
                                        print("An added H is named as ", tHName) 
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
                                print("An added H is named as ", tHName) 
                                nExtra +=1
                                break
                            else:
                                nExtra +=1
                                tHName = hRootId + str(nExtra)
        
        #for aAtom in allAtoms:
        #    if aAtom.GetSymbol() =="H":
        #        aSetBonds = aAtom.GetBonds()
        #        if aAtom.HasProp("Name"):
        #            #print "\nH atom idx ", aAtom.GetIdx(), "  its matched Name ", aAtom.GetProp("Name")
        #            #print "It is in the following bonds: "
        #            #for aB in aSetBonds:
        #            #    print "Bond ", aB.GetIdx()
        #            #    print "Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name"))
        #            #    print "Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name"))
        #        else:
        #            print "\nH atom without name, its idx is ", aAtom.GetIdx()

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
                    print("Bug: a H atom of index %d does not bond to any atom "%aIdxH)
                    sys.exit()
                elif len(aB) >1 : 
                    print("Bug: a H atom of index %d bond to more than one atom "%aIdxH)
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
                        print("Bug: a H atom of index %d is not in the bond of index %d "%(aB.GetIdx()))
                        sys.exit()

                    # print "It bonds atom ", tIdxBH
   
                    if tIdxBH not in tIdxHs:
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
         
            for aIdxBH in list(tIdxHs.keys()):
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
                       print("Atom name %s can not be handled "%name1) 
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
                    print("Geometry of molecule %d can not be optimized within %d circles "%(len(self.allMols), self.nMaxMult*self.nMaxIters))
            elif tNConfId > 1:
                #  multiple conformers
                confIds =AllChem.EmbedMultipleConfs(aMol, tNConfId)
                for aId in confIds:
                    aFailure = self.optOneConformer(aMol, aId)
                    if not aFailure:
                        rdmolops.AssignAtomChiralTagsFromStructure(aMol, aId)
                    else:
                        print("Conformer: ", aId, " is not optimized")
            tMols.append(aMol)

        if len(tMols) != len(self.molecules):
            print("Bug in setInitConformersAllMols ")
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
                print(tReq)
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
            confIds = self.generateMultiComformersByRDKit(tMol)
        if confIds:
            print("Number of initial conformers requested", self.numInitConformers)
            print("Number of number of opt step requested for each conformer ", self.numRDKitOptmSteps)
            print("Number of new conformers ", len(confIds))
            nConf = tMol.GetNumConformers()
            print("Number of initial conformers obtained", nConf)
            allConfs = tMol.GetConformers()
            allAtoms = tMol.GetAtoms()
            formalE = -len(allConfs)    # Formal energy
            #print Chem.MolToMolBlock(tMol)
            iFormalE = 0
            for aConf in allConfs:
                aCIdx =  aConf.GetId()
                aWCoordList = []
                lNorm = self.checkH_Abnormal(aConf, allAtoms, aWCoordList)
                if not lNorm:
                    print("Conformer ", aCIdx, " has abormal coordinates.")
                    self.forceH_coords(aConf, allAtoms, aWCoordList)
                try:
                    aFailure = self.optOneConformer(tMol, self.numRDKitOptmSteps, self.nMaxIters, aCIdx)
                except :
                    # If RDKit can not optimize the conformers, then use the existing conformers and use the formal energies
                    print("Conformer ", aCIdx, " not optimized ")
                    #aForceField = AllChem.UFFGetMoleculeForceField(tMol, confId=aCIdx)
                    #aEng        = aForceField.CalcEnergy()
                    if iFormalE==0:
                        if formalE not in self.conformerEngMap:
                            self.conformerEngMap[formalE] = []
                        self.conformerEngMap[formalE].append(aCIdx)
                        formalE +=1
                        iFormalE +=1
                        #print "Conf : ", aCIdx, " Engergy : ", formalE
                        rdmolops.AssignAtomChiralTagsFromStructure(tMol, aCIdx)
                else:
                    #print "opted id ", aId
                    aForceField = AllChem.UFFGetMoleculeForceField(tMol, confId=aCIdx)
                    aEng        = aForceField.CalcEnergy()
                    if aEng not in self.conformerEngMap:
                        self.conformerEngMap[aEng] = []
                    self.conformerEngMap[aEng].append(aCIdx)
                    #print "Conf : ", aCIdx, " Engergy : ", aEng
                    rdmolops.AssignAtomChiralTagsFromStructure(tMol, aCIdx)
            if len(self.conformerEngMap):
                #print "Current conformers have %d energy levels from UFF force field "%len(self.conformerEngMap)
                #print "They are : "
                #for aEng in sorted(self.conformerEngMap.iterkeys()):
                #    print "Energy ", aEng
                #    for aCid in self.conformerEngMap[aEng]:
                #        print aCid

                #nSelect = 2
                #if not self.useExistCoords: 
                #    nSelect = self.numSelectForRefConfs
                #else:
                #    if tMol.GetNumConformers() <2:
                #        nSelect = tMol.GetNumConformers()

                print("The following conformers are selected for refinement: ")
                nID= 0
                for aEng in sorted(self.conformerEngMap.keys()):
                    for aCId in self.conformerEngMap[aEng]:
                        if nID < self.numSelectForRefConfs:
                            self.selecConformerIds.append(aCId)
                            print("Conformer ID: ", aCId, " UFF energy : ", aEng)
                            nID +=1
                        else:
                            break
            print("Number of conformers selected for refinement is ",  len(self.selecConformerIds))

    def generateMultiComformersByRDKit(self, tMol):

        retConfIds = AllChem.EmbedMultipleConfs(tMol, self.numInitConformers, maxAttempts=0, randomSeed=-1, clearConfs=True)
        if len(retConfIds)==0:
            retConfIds = AllChem.EmbedMultipleConfs(tMol, self.numInitConformers, maxAttempts=0,\
                         randomSeed=-1, clearConfs=True, useRandomCoords=True)

        if len(retConfIds)==0:
            retConfIds = AllChem.EmbedMultipleConfs(tMol, self.numInitConformers, maxAttempts=0,\
                         randomSeed=-1, clearConfs=True, useRandomCoords=True, boxSizeMult=2.0, \
                         randNegEig=True, numZeroFail=1, pruneRmsThresh=-1.0,  coordMap={}, forceTol=0.01)

        return retConfIds

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
                print("for atom ", aIdx)
                print("Its position is :")
                print("x = ", aPos.x)
                print("Is it nan ? ", math.isnan(aPos.x))
                print("Is it inf ", math.isinf(aPos.x))
                print("y = ", aPos.y)
                print("Is it nan ? ", math.isnan(aPos.y))
                print("Is it inf ", math.isinf(aPos.y))
                print("z = ", aPos.z)
                print("Is it nan ? ", math.isnan(aPos.z))
                print("Is it inf ", math.isinf(aPos.z))
               
        return lN
              
    def forceH_coords(self, tConf, tAtoms, tWrongCoordList, tReadInCoordMap=None):
        
        if tReadInCoordMap: 
            for aIdx in tWrongCoordList:
                aPos   = rdGeometry.Point3D()
                if aIdx in tReadInCoordMap:
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
            print("For atom ", name)
            print("Its coordinates are : ")
            print("x:  ", aPos.x)           
            print("y:  ", aPos.y)           
            print("z:  ", aPos.z)           
            
    def setOtherMolInfo(self, tMol, tNumConf, tChemCheck, tPH, tNameMap, tMode=0, tChargeList=None, tMapMode=0):
        
        print("A molecule with residue name %s is generated"%tMol.GetProp("ResidueName"))
        nAtoms      =  tMol.GetNumAtoms()
        print("Number of atoms in the molecule is ", nAtoms)
        initAtoms   =  tMol.GetAtoms()
        # check Atom elements
        elemList = []
        for aAtom in initAtoms:
            elemList.append(aAtom.GetSymbol())
        print("number of props in the mol ", len(tMol.GetPropNames()))
        if not len(elemList) :
            print("No atoms in from your file, check your file format")
        elif not tChemCheck.isOrganic1(elemList):
            print("Your molecule contains METAL or other NON-ORGANIC elements ")
            print("The molecule contains atoms of the following elements ")
            aLine = ""
            for aElem in elemList:
                aLine.append(aElem+ "   ")
            print(aLine)
            sys.exit()
        #print "Number of atoms in this molecule is initially ", nAtoms
        # self.showInfoAboutAtomsAndBonds(aMol, 0)

        if tMapMode == 1:
            self.setNamesForAtomsInMol2(tMol, ChemCheck,  tNameMap, 0)
        else:
            self.setNamesForAtomsInMol(tMol, ChemCheck,  tNameMap, 0)
        # self.showInfoAboutAtomsAndBonds(aMol, 1)
       
        Chem.SanitizeMol(tMol)
        Chem.Kekulize(tMol)

        # Make sure an atom in the molecule has the same charge in the input file. 
        # RDKit sometimes change it when initiating the molecule
        if tChargeList:
            for aAtom in tMol.GetAtoms():
                if aAtom.GetSymbol() != "H":  
                    name  = aAtom.GetProp("Name")
                    if name in tChargeList: 
                        aAtom.SetFormalCharge(int(tChargeList[name]))
                    else:
                        aAtom.SetFormalCharge(0)

        # Extra check added because of bugs in RDKit
        aErrDict = {}
        aErrDict["notExist"]   = []
        aErrDict["wrongOrder"] = []
        aErrDict["needMod"]    = []
        aErrDict["unMod"]      = []
        self.furtherCheckValForAllAtoms(tMol, aErrDict)
        if len(aErrDict["wrongOrder"]) > 0:
            for aErr in aErrDict["wrongOrder"]:
                print(aErr)
            sys.exit()
        if len(aErrDict["needMod"]) > 0:
            print("Number of atoms need to modified ", len(aErrDict["needMod"]))
            for aGrp in aErrDict["needMod"]:
                aModAtom = tMol.GetAtomWithIdx(aGrp[1])
                print("Need to extra val %d to atom %s of %s"%(aGrp[0], aGrp[1], aModAtom.GetSymbol()))    
                newCharge = aModAtom.GetFormalCharge() + aGrp[0]
                oldExHs     = aAtom.GetNumExplicitHs()
                print("before : num of explicit Hs ", oldExHs)
                tMol.GetAtomWithIdx(aGrp[1]).SetNumExplicitHs(aGrp[0]+oldExHs) 
                tMol.GetAtomWithIdx(aGrp[1]).UpdatePropertyCache()
        tMol.UpdatePropertyCache()
        # self.showInfoAboutAtomsAndBonds(aMol, 1)

        if not self.noProtonation:
            if tPH[0] :
                self.setAllFormalChargeFuncGroupAtoms(tMol, tPH[1])
            else:
                self.setAllFormalChargeFuncGroupAtoms(tMol)
            tMol.UpdatePropertyCache()

        if self.useExistCoords:
            aMol = Chem.AddHs(tMol, explicitOnly=False, addCoords=True)
        else:
            aMol = Chem.AddHs(tMol)
        tMol.UpdatePropertyCache()
        #self.showInfoAboutAtomsAndBonds(aMol, 1)
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
 
        # Further: give names to those newly added H atoms
        if tMapMode == 1:
            self.setNamesForAtomsInMol2(aMol,tChemCheck,  tNameMap, 1)
        else:
            self.setNamesForAtomsInMol(aMol, tChemCheck, tNameMap, 1)
        # self.showInfoAboutAtomsAndBonds(aMol, 2)

        #if self.reSetChirals:
        #    self.reAssignChirals(aMol) 
    
        allAtoms = aMol.GetAtoms()
        for aAtom in allAtoms:
            aIdx = aAtom.GetIdx()
            tChemCheck.checkChiralCenters(aMol,aIdx)
            #print "Atom ", aAtom.GetProp("Name")
            #print "Is it a temporal chiral center ", aAtom.HasProp("TmpChiral") 
        # self.showInfoAboutAtomsAndBonds(aMol, 2)
        self.setInitConformersOneMol(aMol)
      
        if len(aMol.GetConformers()) !=0: 
            #rdmolfiles.MolToPDBFile(aMol, "Test_2.pdb")
            aSetTorsions = []
            self.assignTorsions(aMol, aSetTorsions)

            print("Number of torsions in this molecule is ", len(aSetTorsions))
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
                rdmolops.AssignStereochemistry(aMol, cleanIt=False, force=False, flagPossibleStereoCenters=True)
                self.showInfoAboutAtomsAndBonds(aMol, 3)
        else:
            print("RDKit failed to produce any conformers. Acedrg needs at least one conformer ")
            print("Acedrg stops ")

    def showInfoAboutAtomsAndBonds(self, tMol, tLev=0):

        allAtoms = tMol.GetAtoms()
        print("Number of atoms in the molecule is ", len(allAtoms))
        aConf    = 0 
        lConf = False
       
        if tLev ==3: 
            conformers =  tMol.GetConformers()
            if len(conformers) > 0:
                lConf = True
                aConf = conformers[0]

        for aAtom in allAtoms:
            idxA  = aAtom.GetIdx()
            elemA = aAtom.GetSymbol()
            if tLev==1:
                if elemA !="H":
                    name  = aAtom.GetProp("Name") 
            elif tLev > 1:
                name  = aAtom.GetProp("Name") 
            exHs    = aAtom.GetNumExplicitHs()
            imHs    = aAtom.GetNumImplicitHs()
            charge = aAtom.GetFormalCharge()
            val    = aAtom.GetTotalValence()
            print("\nFor atom of index ", idxA)
            print("Its element symbol  ", elemA)
            if tLev==1:
                if elemA !="H":
                    print("Its name ", name)
            elif tLev > 1:
                print("Its name ", name)
            print("Its explicit Hs ", exHs)
            print("Its mplicit Hs ", imHs)
            print("Its formal charge ", charge)
            print("Its total valence ", val)
            if tLev > 2 and aConf:
                pos = aConf.GetAtomPosition(idxA)
                if elemA.find("H") ==-1:
                    print("the coordinates are:")
                    print("X : %7.4f."%pos.x)
                    print("Y : %7.4f"%pos.y)
                    print("Z : %7.4f"%pos.z)

            if tLev > 0:
                print("It is in the following bonds: ")
                aSetBonds = aAtom.GetBonds()
                for aB in aSetBonds:
                    print("Bond ", aB.GetIdx())
                    if tLev > 1:
                        print("Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetProp("Name")))
                        print("Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetProp("Name")))
                    else:
                        print("Its begin atom  %d of %s"%(aB.GetBeginAtomIdx(), allAtoms[aB.GetBeginAtomIdx()].GetSymbol()))
                        print("Its end atom %d of %s "%(aB.GetEndAtomIdx(), allAtoms[aB.GetEndAtomIdx()].GetSymbol()))
                        
                    print("Bond type is ", aB.GetBondType())


    def furtherCheckValForAllAtoms(self, tMol, tErrDict):
        # Temp for bugs appear in RDKit
        # Should be used after the molecule has been kekulized
        checkElems1 = ["I", "P"]
        checkElems2 = ["S", "Se"]
        allAtoms = tMol.GetAtoms() 
        for aAtom in allAtoms:
            # Temporally, just for element S and Se
            aElem = aAtom.GetSymbol().strip()
            if aElem in checkElems1 or aElem in checkElems2 :
                elm    = self.periodicTab.standlizeSymbol(aAtom.GetSymbol().strip())
                if not elm or elm not in self.periodicTab:
                    tErrDict["notExist"].append(aAtom)
                else:
                    aTotalOrder = aAtom.GetTotalValence()
                    self.checkBondOrderAndAllowedVal(aAtom, elm, aTotalOrder, tErrDict)                 
     
    def checkBondOrderAndAllowedVal(self, tAtom, tElm, tOrd, tErrDict):
        checkElems1 = ["I", "P"]
        checkElems2 = ["S", "Se"]
        val    = self.periodicTab[tElm]["val"]
        charge = tAtom.GetFormalCharge()
        idx    = tAtom.GetIdx()
        print("Atom ", tElm)
        actV = val 
        if val > 4:
            actV = 8-val
        allowOrd = actV+charge 
        print("Total bond order calculated from th bonds ", tOrd)
        print("Default allowed bond order ", allowOrd)  
        if tOrd != allowOrd :
            if "extraVal" in self.periodicTab[tElm]:
                aExtraAllow = []
                lMod        = False
                for aV in self.periodicTab[tElm]["extraVal"]:
                    #if aV > 4:
                    aV = 8-aV
                    aExtraAllow.append(aV + charge) 
                    print("ExtraAllowed total bond order ",  aV+charge)
                aExtraAllow.sort()
                if not tOrd in aExtraAllow:
                    if tElem in checkElems2:
                        for aMV in aExtraAllow:
                            if aMV > tOrd:
                                tErrDict["needMod"].append([(aMV-tOrd),idx])
                                lMod = True
                                break
                    elif tElem in checkElems1:
                        tErrDict["wrongOrder"].append("Check, wrong valence on atom "+ tAtom.GetProp("Name") + " !")       
                if not lMod:
                    tErrDict["unMod"].append(tAtom)
        
            else:
                tErrDict["unMod"].append(tAtom)
        
 
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
                print("Chiral center ", name)
                print("CIP rank %s : Stero code %s"%(aATom.GetProp("_CIPRank"), aAtom.GetProp("_CIPCode"))) 

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
        allBonds  = tMol.GetBonds()       
        print("Element to be deleted in the mol is ", self.repSign)
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
            print("Those atoms are linked the del atoms : ")
            for aIdx in tAtomsBondedToDel:
                print(tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
            
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
        print("Ligand ID ", tMonoName)
        print("Group Name ", tGroupName) 
        
        allAtoms         = []
        allAtoms1        = tMol.GetAtoms()
        delAtomIdxs      = []
        atomsBondedToDel = []

        allBonds   = []
        allChirals = []

        nAt  = 0
        nHAt = 0

        #if self.reSetSmi:
        #    self.modifyMol(tMol, allAtoms, allBonds, allChirals, delAtomIdxs, atomsBondedToDel)
        #    nAt  = len(allAtoms)
        #    nHAt = tMol.GetNumHeavyAtoms() - len(delAtomIdxs)
        #else:
        allAtoms = tMol.GetAtoms()
        allBonds = tMol.GetBonds()
        nAt  = len(allAtoms)
        nHAt = tMol.GetNumHeavyAtoms()

        #print "number of atoms with pseudo-atoms is ", tMol.GetNumAtoms()
        print("number of atoms in the initial smiles is  ", len(allAtoms))
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
            #nTetraChi = 0 
            for aAtom in allAtoms:
                pos = aConformer.GetAtomPosition(aAtom.GetIdx())
                #aChi = aAtom.GetChiralTag()
                #if aChi != rdchem.ChiralType.CHI_UNSPECIFIED:
                    #print "Atom ", aAtom.GetProp("Name")
                    #print "Chiral center ? ", aChi
                    #nTetraChi +=1
                #elif aAtom.HasProp("TmpChiral") !=0:
                #    nTetraChi +=1
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
            print("N Bonds ", len(allBonds))
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
                #if self.reSetSmi:
                #    if symb1.find(self.repSign)==-1 and symb2.find(self.repSign)==-1:
                #        bLen  = rdMolTransforms.GetBondLength(aConformer, idx1, idx2)
                #        dBlen = 0.20
                #        aMmCif.write("%s       %s       %s       %s      %s     %5.4f     %5.4f\n" \
                #                     %(tMonoName, name1, name2,  bType, \
                #                       isAro, bLen, dBlen))
                #else:
                bLen  = rdMolTransforms.GetBondLength(aConformer, idx1, idx2)
                dBlen = 0.20
                aMmCif.write("%s       %s       %s       %s      %s     %5.4f     %5.4f\n" \
                              %(tMonoName, name1, name2,  bType, \
                                isAro, bLen, dBlen))

            # chiral center section
            atomNBCIPMap = self.setCIPCodeSerialForNBAtoms(tMol, delAtomIdxs)
            aChiralSignMap = self.setChiralsByMultiConformers(tChemCheck, tMol, atomNBCIPMap)
            self.doubleCheckRDKitChiralCenters(aChiralSignMap, atomNBCIPMap)
            for aId in sorted(aChiralSignMap.keys()):
                if aChiralSignMap[aId]["isChiraled"] and "finalChiVolSign" in aChiralSignMap[aId]:
                    print("==============================================")
                    print("| Centered atom :     %s"%aId)
                    print("----------------------------------------------")
                    for aPair in atomNBCIPMap[aId]:
                        print("NB Atom %s : CIPRank %s "%(aPair[0].GetProp("Name"), aPair[1]))
                    print("----------------------------------------------")
                    print("Its output chiral volume sign   %s"%aChiralSignMap[aId]["finalChiVolSign"])
                    for aCid in list(aChiralSignMap[aId].keys()):
                        if aCid != "isChiraled" and aCid != "finalChiVolSign" and "confSign" in aChiralSignMap[aId][aCid]:
                            print("In conformer %d, its chiral volume sign is %s "%(aCid, aChiralSignMap[aId][aCid]["confSign"]))
                    print("==============================================\n")

            chiCenAtmIds1     = []
            nChiPre = 0
            if tChiDes:
                nChiPre = len(tChiDes)
            print("number of chiral center predefined ", nChiPre)
            if nChiPre !=0:
                for aChiral in tChiDes:
                    chiStrs = aChiral.strip().split()
                    if len(chiStrs) ==7:
                        chiCenAtmIds1.append(chiStrs[2])
       
            #chiCenAtmIds2 = []
            #for aAtom in allAtoms:
            #    aCT = aAtom.GetChiralTag()
            #    if aCT != rdchem.ChiralType.CHI_UNSPECIFIED:
            #        aAtmName=aAtom.GetProp("Name")
                    
            #        if not aAtmName in chiCenAtmIds1:
            #            if aChiralSignMap[aId]["isChiraled"] and "finalChiVolSign" in aChiralSignMap[aId]:
            #                print("Chiral center %s is not in predefined chiral centers"%aAtmName)
            #                chiCenAtmIds2.append(aAtmName)
            #nTetraChi = len(chiCenAtmIds2)
            #print(" Number of chiral centers get from the conformer ", nTetraChi)
            chiCenAtms3 = []
            for aAtom in allAtoms:
                aElem = aAtom.GetSymbol()
                aName=aAtom.GetProp("Name")
                aHyb  = aAtom.GetHybridization()
                nNB   = len(aAtom.GetNeighbors())
                #print("for atom ", aName)
                #print("its hyb is ", aHyb)
                #print("it connects to ",nNB, " NBs")
                if aHyb == rdchem.HybridizationType.SP3 and aElem != "H" and aElem != "O":
                    nConnHs =  self.chemCheck.getNumNBHAtoms(aAtom)
                    #print("number of H atoms connected  ", nConnHs)
                    if not aName in chiCenAtmIds1 and nConnHs < 2 and nNB > 2:
                        chiCenAtms3.append(aAtom)                        
            nChiBoth = len(chiCenAtms3)
            print("nChi ",nChiBoth)
            if nChiPre !=0 or nChiBoth !=0 :
                # The molecule contain chiral centers
                aMmCif.write("#\n")
                aMmCif.write("_chem_comp_chir.comp_id\n")
                aMmCif.write("_chem_comp_chir.id\n")
                aMmCif.write("_chem_comp_chir.atom_id_centre\n")
                aMmCif.write("_chem_comp_chir.atom_id_1\n")
                aMmCif.write("_chem_comp_chir.atom_id_2\n")
                aMmCif.write("_chem_comp_chir.atom_id_3\n")
                aMmCif.write("_chem_comp_chir.volume_sign\n")

            chiralIdx = 1
            if nChiPre:
                for aChiral in tChiDes:
                    aMmCif.write(aChiral+"\n")
                    print("aChiral ")
                    print(aChiral)
                print("Predefined chiral centres ", chiCenAtmIds1)
            #if nTetraChi : 
            #    # output all chiral centers in form of mmCif
            #    chiralIdx = nChiPre + 1
            #    for aId in chiCenAtmIds2 :
            #        if aChiralSignMap[aId]["isChiraled"] and "finalChiVolSign" in aChiralSignMap[aId]:
            #            #print("Chiral center %s is not in predefined chiral centers"%aId)
            #            aCTName = "chir_" + str(chiralIdx)
            #            aLine = "%s%s%s%s%s%s%s\n"\
            #                    %(tMonoName.ljust(12), aCTName.ljust(12),\
            #                      aId.ljust(12),\
            #                      atomNBCIPMap[aId][0][0].GetProp("Name").ljust(12),\
            #                      atomNBCIPMap[aId][1][0].GetProp("Name").ljust(12),\
            #                      atomNBCIPMap[aId][2][0].GetProp("Name").ljust(12),\
            #                      aChiralSignMap[aId]["finalChiVolSign"])
            #            aMmCif.write(aLine)
            #            print(aLine)
            #            chiralIdx +=1

            if nChiBoth :
                for aAtom in chiCenAtms3 :
                    aId = aAtom.GetProp("Name")  
                    aCTName = "chir_" + str(chiralIdx)
                    aLine = "%s%s%s%s%s%s%s\n"\
                            %(tMonoName.ljust(12), aCTName.ljust(12),\
                              aId.ljust(12),\
                              atomNBCIPMap[aId][0][0].GetProp("Name").ljust(12),\
                              atomNBCIPMap[aId][1][0].GetProp("Name").ljust(12),\
                              atomNBCIPMap[aId][2][0].GetProp("Name").ljust(12),\
                              aChiralSignMap[aId]["finalChiVolSign"])
                    aMmCif.write(aLine)
                    print(aLine)
                    chiralIdx +=1
            aMmCif.close()
        

    def checkUncertainChirals(self, tNBAtoms):
        
        tPass = True
        tOneBondOs = 0
        tNBRanks   = {}
        for aPair in tNBAtoms:
            if aPair[0].GetSymbol().find("O") !=-1 and len(aPair[0].GetBonds())==1:
                tOneBondOs +=1 
            if not aPair[1] in tNBRanks:
                tNBRanks[aPair[1]] = []
            tNBRanks[aPair[1]].append(1)

        if tOneBondOs > 1 :
            tPass = False
        else :
            for aKey in list(tNBRanks.keys()):
                if len(tNBRanks[aKey]) > 1:
                    tPass = False

        return tPass

   
    def setCIPCodeSerialForNBAtoms(self, tMol, delAtomIdxs):

        reNameSet = {}
        
        rdmolops.AssignStereochemistry(tMol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

        allAtoms = tMol.GetAtoms()

        for aAtom in allAtoms:
            aIdx  = aAtom.GetIdx()
            aId   = aAtom.GetProp("Name")
            reNameSet[aId] = []
            aSetBonds = aAtom.GetBonds()
            nNB = len(aAtom.GetNeighbors())
            #print("nNB %d"%nNB)
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
                                #print name2
                                reNameSet[aId].append([allAtoms[atmIdx2], allAtoms[atmIdx2].GetProp("_CIPRank")])
                            elif atmIdx2 == aIdx:
                                #print name1
                                reNameSet[aId].append([allAtoms[atmIdx1], allAtoms[atmIdx1].GetProp("_CIPRank")])
                            else:
                                print("Bug! atom %s is not in bonds obtained by aAtom.GetBonds()"%(aAtom.GetProp("Name")))
                                break
                else:
                    if atmIdx1 == aIdx:                      
                        reNameSet[aId].append([allAtoms[atmIdx2], allAtoms[atmIdx2].GetProp("_CIPRank")])
                    elif atmIdx2 == aIdx:
                        reNameSet[aId].append([allAtoms[atmIdx1], allAtoms[atmIdx1].GetProp("_CIPRank")])
                    else:
                        print("Bug! atom is not in bonds obtained by aAtom.GetBonds()"%(aAtom.GetProp("Name")))
                        break
            #reNameSet[aId].sort(key=listCompDes) 
            reNameSet[aId]=sorted(reNameSet[aId], key=cmp_to_key(listCompDes))
            #if len(reNameSet[aId]) > 2:
            #    print("center atom ",aId)
            #    print(len(reNameSet[aId]))
            #    for aPair in reNameSet[aId]:
            #        print(aPair[0].GetProp("Name"), " : ", aPair[1])
        for aId in reNameSet.keys():
            if aId.find("H")==-1:
                print("Atom ",aId)
                print("Its NB CIP is arranged as follow ")
                for aPair in reNameSet[aId]:
                    print ("Atom %s : CIP RANK %s "%(aPair[0].GetProp("Name"), aPair[1]))
        return reNameSet 

    def setChiralsByMultiConformers(self, tChemCheck, tMol, tNBCIPMap):

        """
           Using several lowest energy conformers to check if the chiral atoms
           should have "both" as the volume signs. 
        """
        
        aChiralSetMap = {}       #  aChiralSetMap[cenAtomIdx][confId][chiralVolSign]
 
        allAtoms  = tMol.GetAtoms()  
 
        if len(self.selecConformerIds):
            for aCid in self.selecConformerIds:
                rdmolops.AssignAtomChiralTagsFromStructure(tMol, confId=aCid)
                aConformer =  tMol.GetConformer(aCid)
                #print "In conformer ", aCid
                for aAtom in allAtoms:
                    aCT = aAtom.GetChiralTag()
                    #print "Atom ", aAtom.GetProp("Name")
                    #print "RDKit sign ", aCT
                    aTmpCT = aAtom.HasProp("TmpChiral")
                    aId  = aAtom.GetProp("Name")
                    aIdx = aAtom.GetIdx()
                    if aId not in aChiralSetMap:
                        aChiralSetMap[aId] = {}
                        if aCT != rdchem.ChiralType.CHI_UNSPECIFIED or aTmpCT !=0 :
                            aChiralSetMap[aId]["isChiraled"] = True
                        else:
                            aChiralSetMap[aId]["isChiraled"] = False
                            
                    if aChiralSetMap[aId]["isChiraled"] :
                        aChiralSetMap[aId][aCid] = {}
                        aPass = self.checkUncertainChirals(tNBCIPMap[aId])
                        if not aPass:
                            aChiralSetMap[aId]["finalChiVolSign"] = "both"
                        else:
                            if "finalChiVolSign" in aChiralSetMap[aId] and aChiralSetMap[aId]["finalChiVolSign"]=="both":
                                # Already in "both" sign, no need to further check
                                pass
                            else:
                                if aId in tNBCIPMap and len(tNBCIPMap[aId]) > 2: # the second one should be for insurance
                                    posCen = aConformer.GetAtomPosition(aIdx)
                                    coordsCen = [float(posCen.x), float(posCen.y), float(posCen.z)] 
                                    pos1 = aConformer.GetAtomPosition(tNBCIPMap[aId][0][0].GetIdx())
                                    coords1 = [float(pos1.x), float(pos1.y), float(pos1.z)] 
                                    pos2 = aConformer.GetAtomPosition(tNBCIPMap[aId][1][0].GetIdx())
                                    coords2 = [float(pos2.x), float(pos2.y), float(pos2.z)] 
                                    pos3 = aConformer.GetAtomPosition(tNBCIPMap[aId][2][0].GetIdx())
                                    coords3 = [float(pos3.x), float(pos3.y), float(pos3.z)]
                                    aChiralSetMap[aId][aCid]["confSign"] = tChemCheck.getChiralVolumeSign(posCen, pos1, pos2, pos3)
                                    if "finalChiVolSign" not in aChiralSetMap[aId]:
                                        aChiralSetMap[aId]["finalChiVolSign"] = aChiralSetMap[aId][aCid]["confSign"]
                                    #if aChiralSetMap[aId][aCid]["confSign"] != aChiralSetMap[aId]["finalChiVolSign"]:
                                    #    aChiralSetMap[aId]["finalChiVolSign"] = "both"
                                    
        return aChiralSetMap 

    def doubleCheckRDKitChiralCenters(self, tChiralSetMap, tNBCIPMap):

        # some RDKit assignes a chiral center to an atom which connects 3 or 4 NB atoms, with at least 2 of them 
        # have the same CIP-RANK numbers. e.g. the atom of id "P" in mononer "AFS"
        # Need to exclude those atoms in the list of chiral centers    
        for aIdx in list(tChiralSetMap.keys()):
            if "isChiraled" in tChiralSetMap[aIdx] and tChiralSetMap[aIdx]["isChiraled"]:
                aCIPMap = {}         # aCIPMap[aCIP-Rank] = [NB_atomId1,NB_atomId2, .....]
                for aPair in tNBCIPMap[aIdx]:
                    if aPair[1] not in aCIPMap:
                        aCIPMap[aPair[1]] = []
                    aCIPMap[aPair[1]].append(aPair[0].GetProp("Name"))
         
                for aCIP in list(aCIPMap.keys()):
                    if len(aCIPMap[aCIP]) > 1: 
                        tChiralSetMap[aIdx]["isChiraled"] = False
  
                    
    def setupFuncGroupTab(self, tFuncFileName):

        try:
            tFuncF = open(tFuncFileName, "r")
        except IOError:
            print("%s can not be open for reading "%tFuncFileName)
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
        for aKey in sorted(self.funcGroupTab.keys()):
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
                if aTri[0] not in self.funcGroups:
                    self.funcGroups[aTri[0]] = []
                for oneAtmIdxGrp in atmGrpIdxs:
                    #print aTri[0], "add func atoms"
                    self.funcGroups[aTri[0]].append(oneAtmIdxGrp)
                #print ""

    # The following methods are migrated from Acedrg c++ section

    def setAllFormalChargeFuncGroupAtoms(self, tMol, tPH=7.0):
      
        self.getAllFuncGroupAtoms(tMol)
        if len (self.funcGroups):
            for aFuncG in list(self.funcGroups.keys()):
                # print "check ", aFuncG 
                if aFuncG.find("CARBOXY-AMINO-TERS") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeC_A_T(tMol, aFuncG,  self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-TER")!=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeC_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-ASP")!=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeC_ASP(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("AMINO-TER") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeA_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("CARBOXY-AMINO-GLY") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeGLY(tMol, aFuncG,  self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("VAR-A-TER") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeV_A_T(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-LYS") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeNH_LYS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-HIS") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeNH_HIS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-PRO") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeNH_PRO(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NH-ARG") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeNH_ARG(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SH-CYS") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeSH_CYS(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-TYR") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeBZ_TYR(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-NH") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeBZ_NH(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("BZ-NHA") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeBZ_NHA(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SO3") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeSO3(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("SO4") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargeSO4(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO4") !=-1:
                    print("Doing ", aFuncG)
                    self.setFormalChargePO4(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO3R") !=-1:
                    print("Doing ", aFuncG)
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargePO3R(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("PO3N") !=-1:
                    print("Doing ", aFuncG)
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargePO3N(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NITROMETHANE0") !=-1:
                    print("Doing ", aFuncG)
                    #print self.funcGroups[aFuncG] 
                    self.setFormalChargeNITROMETHANE(tMol, aFuncG, self.funcGroups[aFuncG], tPH)
                elif aFuncG.find("NITROMETHANE1") !=-1:
                    print("Doing ", aFuncG)
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
                        #print "Total H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()
                        #print "Exp   H ", tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                            and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==0:     
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
     

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
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
         
    def setFormalChargeC_ASP(self, tMol, tFunG, tAtomIdxs, tPH):

        tPka = self.funcGroupTab[tFunG][1]
        print "Pka ", tPka
        print "PH  ", tPH
        if tPH > tPka :
            for aSetIdxs in tAtomIdxs:
                #print "a set of atoms: "
                #print aSetIdxs
                for aIdx in aSetIdxs:
                    print "Atom : ", tMol.GetAtomWithIdx(aIdx).GetProp("Name")
                    print "exH  : ", tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()
                    print "imH  : ", tMol.GetAtomWithIdx(aIdx).GetNumImplicitHs()
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="O":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 0\
                            and (tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs() ==1 
                                 or tMol.GetAtomWithIdx(aIdx).GetNumImplicitHs()==1):  
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNumExplicitHs(0)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
     
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
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
                    if tMol.GetAtomWithIdx(aIdx).GetSymbol()=="N":
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() >= 2:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                   tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
         
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
                        print("atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
                        print("total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence())
                        print("number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds()))
                        print("num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()) 
                        print("num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors()))
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
                        print("atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
                        print("total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence())
                        print("number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds()))
                        print("num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()) 
                        print("num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors()))
                        if tMol.GetAtomWithIdx(aIdx).GetTotalValence() == 3 \
                            and tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() > 1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache()
                            print("atom %s has a charge %d "%(tMol.GetAtomWithIdx(aIdx).GetProp("Name"), \
                                  tMol.GetAtomWithIdx(aIdx).GetFormalCharge()))
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
                        print("atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
                        print("total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence())
                        print("number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds()))
                        print("num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs()) 
                        print("num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors()))
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() == 1     \
                           and tMol.GetAtomWithIdx(aIdx).GetTotalValence()==3 \
                           and len(tMol.GetAtomWithIdx(aIdx).GetBonds())==1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                            tMol.GetAtomWithIdx(aIdx).UpdatePropertyCache() 
                            print("Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
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
                        print("\nHere are details of O: ")           
                        print("atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
                        print("total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence())
                        print("number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds()))
                        print("num of H ", tMol.GetAtomWithIdx(aIdx).GetTotalNumHs())
                        print("num of neighbor ", len(tMol.GetAtomWithIdx(aIdx).GetNeighbors())) 
                        if tMol.GetAtomWithIdx(aIdx).GetTotalNumHs() != 0\
                           and tMol.GetAtomWithIdx(aIdx).GetNumExplicitHs()==0 :
                            # need to deprotonation for the O atom 
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(-1)
                            tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                            print("Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge())
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
                        print("\nHere are details of O: ")
                        print("atom %s :"%tMol.GetAtomWithIdx(aIdx).GetProp("Name"))
                        print("total Val ", tMol.GetAtomWithIdx(aIdx).GetTotalValence())
                        print("number of bonds ", len(tMol.GetAtomWithIdx(aIdx).GetBonds()))
                        if tMol.GetAtomWithIdx(aIdx).GetFormalCharge() !=1:
                            tMol.GetAtomWithIdx(aIdx).SetFormalCharge(1)
                        tMol.GetAtomWithIdx(aIdx).SetNoImplicit(True)
                        print("Its formal charge ", tMol.GetAtomWithIdx(aIdx).GetFormalCharge())  
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
        print("\nHere are details of O: ")           
        print("atom %s :"%tMol.GetAtomWithIdx(tIdx).GetProp("Name"))
        print("total Val ", tMol.GetAtomWithIdx(tIdx).GetTotalValence())
        print("number of bonds ", len(tMol.GetAtomWithIdx(tIdx).GetBonds()))
        print("num of H ", tMol.GetAtomWithIdx(tIdx).GetTotalNumHs())
        print("num of neighbor ", len(tMol.GetAtomWithIdx(tIdx).GetNeighbors())) 
        if tMol.GetAtomWithIdx(tIdx).GetTotalNumHs() != 0\
            and tMol.GetAtomWithIdx(tIdx).GetNumExplicitHs()==0 :
            tMol.GetAtomWithIdx(tIdx).SetFormalCharge(-1)
            tMol.GetAtomWithIdx(tIdx).SetNoImplicit(True)
            print("Its formal charge ", tMol.GetAtomWithIdx(tIdx).GetFormalCharge())
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
       
        print("Number of atoms ", tMol.GetNumAtoms()) 
        # Remove water molecules read from the input files (usually pdb format)
        aEdMol = rdchem.EditableMol(tMol)
        for aAt in tMol.GetAtoms():
            if aAt.GetSymbol() =="O":
                numBs = aAt.GetTotalDegree()
                numHs = aAt.GetTotalNumHs()
                if numBs == numHs:
                    print("Atom's element type ", aAt.GetSymbol())
                    print("Number of bonds it has ", numBs)
                    print("Number of H  it connects ", numHs)
                    aEdMol.RemoveAtom(aAt.GetIdx())
                    print("atom ", aAt.GetIdx(), " is removed ")
        tMol = aEdMol.GetMol()
        print("Number of atoms after removing water", tMol.GetNumAtoms()) 
