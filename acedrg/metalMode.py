#!/usr/bin/env python
# Python script

import os,os.path,sys


from functools  import cmp_to_key


from .  exebase    import CExeCode

from . filetools   import FileTransformer

from .  utility    import BondOrderS2N
from .  utility    import getLinkedGroups


class metalMode(CExeCode):

    def __init__( self, tTabLoc):
        
        self.acedrgTables         = tTabLoc
        
        self.remainAtoms          = []
        self.metalAtoms           = []
        self.metalConnAtoms       = []
        self.metalConnAtomsMap    = {}
        self.metalConnFullAtoms   = []
        self.nonMAtmConnsMap      = {}
        self.connMAMap            = {}
        self.atmNonHMap           = {}
        
        self.addedHs          = []
        
        self.remainBonds      = []
        self.metalBonds       = []
        self.mcAtomBonds      = {}
        self.addedHBs         = []
        
        self.speAngs          = []
        self.metalPA          = []
        self.atmHybr          = {}
        self.simpP            = {}
        self.simpR            = {}   # simple rings (full fledge rings are in c++ libmol )
        self.metalInPs        = {}
        
        self.dataDiescriptor  = []
        
        self.outRstCifName    = ""
        
        self.mols             = []
        
        self.monomRoot        = ""
        self.outRoot          = ""

        self.getRadii() 
        
        """
        check
        """
        
    def execute(self, tAtoms, tBonds, tMonomRoot, tOutRoot, tFileConv, tChem, tVersionInfo):
        
        aRet = ""
        self.monomRoot = tMonomRoot
        self.outRoot   = tOutRoot 
        
        print("Stage 1 ")
        
        
        self.getNewMolWithoutMetal(tAtoms, tBonds, tChem)
        
        tmpMmcifName = self.monomRoot + "_tmpIn.cif"
    
        self.setCCP4MonDataDescritor(tFileConv.dataDescriptor)
        
        self.outInitModifiedCif(tmpMmcifName)
        
        nBonds = len(self.remainBonds) 
        
        if os.path.isfile(tmpMmcifName) and nBonds:
            print("Stage 2 ")
            print("acedrg is running")
            outTmpRootName = self.outRoot + "_tmp"
            self.runAcedrg(tmpMmcifName, outTmpRootName)
            outTmpCifName  =  outTmpRootName + ".cif"
            
            
            #for aA in self.metalAtoms:
            #    print("Metal atom ", aA['_chem_comp_atom.atom_id'], 
            #          " carries ", aA['_chem_comp_atom.charge'], 
            #          " formal charge ")
    
            if os.path.isfile(outTmpCifName):
                outTmp2CifName = self.outRoot   + "_tmp2.cif"
                self.modiCif(outTmpCifName, outTmp2CifName)
                if os.path.isfile(outTmp2CifName):
                    aRet = outTmp2CifName
        elif nBonds==0:
            outInterMCifName = self.outRoot   + "_final.cif"
            self.outInterMCif(outInterMCifName, tAtoms, tBonds, tVersionInfo)
            print("==============================================")
            print("No non-metal related bonds exist in the molecule")
            #print("The final output file is %s"%outInterMCifName)
            print("==============================================")
        elif not os.path.isfile(tmpMmcifName):
            print("==============================================")
            print("Bug: no input file %s "%outTmpCifName)
            print("==============================================")
                
        return aRet         
    
    def executeStage1Only(self, tAtoms, tBonds, tMonomRoot, tOutRoot, tFileConv, tChem):
        
        aRet = ""
        self.monomRoot = tMonomRoot
        self.outRoot   = tOutRoot 
        
        print("Stage 1 : prepare input file for ", )
        
        
        self.getNewMolWithoutMetal(tAtoms, tBonds, tChem)
        
        tmpMmcifName = self.monomRoot + "_tmpIn.cif"
    
        self.setCCP4MonDataDescritor(tFileConv.dataDescriptor)
        
        self.outInitModifiedCif(tmpMmcifName)
        
        nBonds = len(self.remainBonds) 
        
        if os.path.isfile(tmpMmcifName) and nBonds:
            print("Stage 2 ")
            print("acedrg is running")
            outTmpRootName = self.outRoot + "_tmp"
            self.runAcedrg(tmpMmcifName, outTmpRootName)
            outTmpCifName  =  outTmpRootName + ".cif"
            if os.path.isfile(outTmpCifName):
                aRet = outTmpCifName
                
        return aRet

    def getRadii(self):

        aFN = os.path.join(self.acedrgTables, "radii.table") 
        if os.path.isfile(aFN):
            aF = open(aFN, "r")
            allLs = aF.readlines()
            aF.close()
            
            self.radii = {}
            for aL in allLs:
                strs = aL.strip().split()
                if len (aL) > 0 and len(strs)==3 and aL[0].find("#")==-1:
                    strs[0] = strs[0].upper()
                    if not strs[0] in self.radii.keys():
                        self.radii[strs[0]] = {}
                    self.radii[strs[0]][strs[1]] = float(strs[2].strip())
    
    def setBondValueFromRadii(self, tAtom1, tAtom2, tBond):
        
        aElem1 = tAtom1['_chem_comp_atom.type_symbol'].upper()
        aElem2 = tAtom2['_chem_comp_atom.type_symbol'].upper() 
        #print ("Atom1 ", tAtom1['_chem_comp_atom.atom_id'], " and Atom2 ", tAtom2['_chem_comp_atom.atom_id'])
        #print("Elem1 ", aElem1, " and Elem2 ", aElem2)
        if aElem1 in self.radii and aElem2 in self.radii:
            if "cova" in self.radii[aElem1] and "cova" in self.radii[aElem2]:
                aVal = self.radii[aElem1]["cova"] + self.radii[aElem2]["cova"]      
                tBond["_chem_comp_bond.value_dist"] = "%4.3f"%aVal
            else:
                tBond["_chem_comp_bond.value_dist"] = "2.0"
        else:
            tBond["_chem_comp_bond.value_dist"]  = "2.0"
        #print("Bond length %s "%tBond["_chem_comp_bond.value_dist"])
            
    
    def setMMBond(self, tAtoms, tBonds):
        
        for aBond in tBonds:
            if not "_chem_comp_bond.value_dist" in aBond.keys():
                atm1 = self.getAtomById(tAtoms, aBond['_chem_comp_bond.atom_id_1'])
                atm2 = self.getAtomById(tAtoms, aBond['_chem_comp_bond.atom_id_2'])
                self.setBondValueFromRadii(atm1, atm2, aBond)
                aBond['_chem_comp_bond.value_dist_esd'] = "0.04"
                aBond['_chem_comp_bond.value_dist_nucleus']     = aBond['_chem_comp_bond.value_dist']
                aBond['_chem_comp_bond.value_dist_nucleus_esd'] = "0.01"
                
                
        
    def getNewMolWithoutMetal(self, tAtoms, tBonds, tChem):
        
        
        
        
        metalIds = []
        for aAtom in tAtoms:
            if "_chem_comp_atom.type_symbol" in aAtom.keys():
                print("Atom ", aAtom['_chem_comp_atom.atom_id'], 
                      " is ", aAtom["_chem_comp_atom.type_symbol"])
                if not aAtom['_chem_comp_atom.type_symbol'] in tChem.organicSec:
                    #print("atom %s is a metal atom of elem %s "
                    #      %(aAtom['_chem_comp_atom.atom_id'], aAtom['_chem_comp_atom.type_symbol'])) 
                    self.metalAtoms.append(aAtom)
                    metalIds.append(aAtom['_chem_comp_atom.atom_id'])
                else:
                    self.remainAtoms.append(aAtom)
                    
                    
        print("number of all atoms is ", len(tAtoms))
        print("number of metal atoms is ", len(self.metalAtoms))
        print("number of org atoms is ", len(self.remainAtoms))
        #for aA in tAtoms:
        #    print(aA.keys())
        
        
        atmConnsMap = {} 
    
        for aBond in tBonds:
            atm1 = self.getAtomById(tAtoms, aBond['_chem_comp_bond.atom_id_1'])
            atm2 = self.getAtomById(tAtoms, aBond['_chem_comp_bond.atom_id_2'])
            
            if aBond['_chem_comp_bond.atom_id_1'] in metalIds:
                if not aBond['_chem_comp_bond.atom_id_1'] in self.metalConnAtomsMap:
                    self.metalConnAtomsMap[aBond['_chem_comp_bond.atom_id_1']] = []
                self.metalConnAtomsMap[aBond['_chem_comp_bond.atom_id_1']].append(aBond['_chem_comp_bond.atom_id_2'])
                
                self.setBondValueFromRadii(atm1, atm2, aBond)
                self.metalBonds.append(aBond)
                if not aBond['_chem_comp_bond.atom_id_2'] in metalIds:
                    if not aBond['_chem_comp_bond.atom_id_2'] in self.connMAMap:
                        self.connMAMap[aBond['_chem_comp_bond.atom_id_2']] = []
                    self.connMAMap[aBond['_chem_comp_bond.atom_id_2']].append(aBond['_chem_comp_bond.atom_id_1'])
            if aBond['_chem_comp_bond.atom_id_2'] in metalIds:
                if not aBond['_chem_comp_bond.atom_id_2'] in self.metalConnAtomsMap:
                    self.metalConnAtomsMap[aBond['_chem_comp_bond.atom_id_2']] = []
                self.metalConnAtomsMap[aBond['_chem_comp_bond.atom_id_2']].append(aBond['_chem_comp_bond.atom_id_1'])
                self.setBondValueFromRadii(atm1, atm2, aBond)
                self.metalBonds.append(aBond)
                if not aBond['_chem_comp_bond.atom_id_1'] in metalIds:
                    if not aBond['_chem_comp_bond.atom_id_1'] in self.connMAMap:
                        self.connMAMap[aBond['_chem_comp_bond.atom_id_1']] = []
                    self.connMAMap[aBond['_chem_comp_bond.atom_id_1']].append(aBond['_chem_comp_bond.atom_id_2'])
            
                
            if not aBond['_chem_comp_bond.atom_id_1'] in metalIds and not aBond['_chem_comp_bond.atom_id_2'] in metalIds:
                self.remainBonds.append(aBond)
                
                if not aBond['_chem_comp_bond.atom_id_1'] in atmConnsMap.keys():
                    atmConnsMap[aBond['_chem_comp_bond.atom_id_1']] = []
                if not aBond['_chem_comp_bond.atom_id_2'] in atmConnsMap.keys():
                    atmConnsMap[aBond['_chem_comp_bond.atom_id_2']] = []    
                if not aBond['_chem_comp_bond.atom_id_2'] in atmConnsMap[aBond['_chem_comp_bond.atom_id_1']]:
                    atmConnsMap[aBond['_chem_comp_bond.atom_id_1']].append(aBond['_chem_comp_bond.atom_id_2'])
                if not aBond['_chem_comp_bond.atom_id_1'] in atmConnsMap[aBond['_chem_comp_bond.atom_id_2']]:
                    atmConnsMap[aBond['_chem_comp_bond.atom_id_2']].append(aBond['_chem_comp_bond.atom_id_1'])
                
                if atm1 and atm2:
                    aElem1 = atm1['_chem_comp_atom.type_symbol'] 
                    aElem2 = atm2['_chem_comp_atom.type_symbol']          
                    aId1   = atm1['_chem_comp_atom.atom_id']
                    aId2   = atm2['_chem_comp_atom.atom_id'] 
                    
                    
                    if not aId1 in self.atmNonHMap and aElem1 !="H":
                        self.atmNonHMap[aId1] = []
                    if not aId2 in self.atmNonHMap and aElem2 !="H":
                        self.atmNonHMap[aId2] = []
                    if aElem1 !="H":
                        if aElem2 != "H":
                            self.atmNonHMap[aId1].append(aId2)
                            
                    if aElem2 !="H":
                        if aElem1 != "H":
                            self.atmNonHMap[aId2].append(aId1)
                            
                    
        # Check 
        #for aMAtm in self.metalConnAtomsMap:
        #    print("%s is a metal atom "%aMAtm)
        #    if len(self.metalConnAtomsMap[aMAtm]) > 0:
        #        print("It connects to the following atoms")
        #        for aCAtm in self.metalConnAtomsMap[aMAtm]:
        #            print("atom : ", aCAtm)
        
        
        
                
        for aAtm in self.remainAtoms:
            aId = aAtm['_chem_comp_atom.atom_id']
            if not aId in atmConnsMap:
                # Isolated atoms
                atmConnsMap[aId] = {}
            #aN  = len(atmConnsMap[aId])
            #print("Atom %s has %d non metal atom connections."%(aId, aN))
            #if aN > 0:
                #print("They are: ")
                #for aCId in atmConnsMap[aId]:
                #    print("atom ", aCId)
        
        
        # Get all fragments 
        aSetFrags = {}
        getLinkedGroups(self.remainAtoms, atmConnsMap, aSetFrags)
        
        
        idxR = {}
        for aMA in self.metalConnAtomsMap.keys():
            i =0
            for aBond in self.remainBonds:
                id1 = aBond['_chem_comp_bond.atom_id_1']
                id2 = aBond['_chem_comp_bond.atom_id_2']
                if id1 in self.metalConnAtomsMap[aMA]:
                    if not id1 in  idxR:       #  self.mcAtomBonds.keys():
                        idxR[id1]=[]
                    if not i in idxR[id1]:
                        idxR[id1].append(i)
                    if not id1 in self.mcAtomBonds.keys():
                        self.mcAtomBonds[id1] = []
                    #self.mcAtomBonds[id1].append(aBond)
                
                if id2 in self.metalConnAtomsMap[aMA]:
                    if not id2 in  idxR:       #  self.mcAtomBonds.keys():
                        idxR[id2]=[]
                    if not i in idxR[id2]:
                        idxR[id2].append(i)
                    if not id2 in self.mcAtomBonds.keys():
                       self.mcAtomBonds[id2] = []
                    #self.mcAtomBonds[id2].append(aBond)
                i+=1
        
        for aK in idxR:
            for aI in idxR[aK]:
                self.mcAtomBonds[aK].append(self.remainBonds[aI])
                
        #for aMC in self.mcAtomBonds:
        #    print("a metal bonding atom ", aMC)
        #    print("It has ", len(self.mcAtomBonds[aMC]), " non metal bonds")
        #    for aB in self.mcAtomBonds[aMC]:
        #        print(aB['_chem_comp_bond.atom_id_1'])
        #        print(aB['_chem_comp_bond.atom_id_2'])
        #        if '_chem_comp_bond.value_order' in aB:
        #            print(aB['_chem_comp_bond.value_order'])
        #        elif '_chem_comp_bond.type' in aB:
        #            print(aB['_chem_comp_bond.type'])
                    
        if len(self.metalConnAtomsMap):
            #print("These atoms are : ")
            #for aA in self.metalConnAtomsMap.keys():
            #    print(aA)
            
            cKeys = ['_chem_comp_atom.atom_id',  '_chem_comp_atom.type_symbol',
                     '_chem_comp_atom.charge', '_chem_comp_atom.x', 
                     '_chem_comp_atom.y', '_chem_comp_atom.z']
            otherKeys = []
            for aK in tAtoms[0].keys():
                if not aK in cKeys:
                    otherKeys.append(aK)
                    
            bKeys = ['_chem_comp_bond.atom_id_1',  '_chem_comp_bond.atom_id_2',
                     '_chem_comp_bond.value_order', '_chem_comp_bond.type',
                     '_chem_comp_bond.value_dist_nucleus', '_chem_comp_bond.value_dist_nucleus_esd',
                     '_chem_comp_bond.value_dist', '_chem_comp_bond.value_dist_esd']
            bOtherKeys = []
            for aBK in tBonds[0].keys():
                if not aBK in bKeys:
                    bOtherKeys.append(aBK)
                    
            idxH =0
            doneOrgAtms = []
            
            for aMA in self.metalConnAtomsMap:
                print("For a metal atom : ", aMA)
                for aMCId in self.metalConnAtomsMap[aMA]:
                    print("For metal bonding atom ", aMCId)
                    aMCAtom = self.getAtomById(tAtoms, aMCId)
                    if aMCAtom['_chem_comp_atom.type_symbol'] in tChem.organicSec\
                       and not aMCId in doneOrgAtms:
                        doneOrgAtms.append(aMCId)
                        self.metalConnFullAtoms.append(aMCAtom)
                        aMCAtom['_chem_comp_atom.charge'] = 0
                        if aMCAtom !=None:
                            #numC   = -self.checkVal(aMCAtom)
                            #print("It should carry ", numC, " ext charge")
                            #aMCAtom['_chem_comp_atom.charge'] = numC
                            aMCAtom['_chem_comp_atom.charge'] = -self.checkVal(aMCAtom, tChem.defaultBo) 
                            print("It carries ", aMCAtom['_chem_comp_atom.charge'], " charges")
                            
        
        for aNonM in atmConnsMap.keys():
            if not aNonM in self.nonMAtmConnsMap.keys():
                self.nonMAtmConnsMap[aNonM] = []
            for bNonM in atmConnsMap[aNonM]:
                self.nonMAtmConnsMap[aNonM].append(bNonM)
                
        # Check
        self.setAtomHybr(tAtoms)
        
        #for aA in tAtoms:
        #    print(aA['_chem_comp_atom.atom_id'], " carries ", aA['_chem_comp_atom.charge'])
         
    def setAtomHybr(self, tAtoms):
        
        #print(self.nonMAtmConnsMap)
        # First round
        for aAtm in tAtoms:
            aId = aAtm['_chem_comp_atom.atom_id']
            #print(aId)
            self.atmHybr[aId] =0
            if aId in self.nonMAtmConnsMap.keys():
                aM=0
                if aId in self.connMAMap:
                    aM = len(self.connMAMap[aId])
                aL = len(self.nonMAtmConnsMap[aId])
                #print("conn ", aL)
                if aL > 4 :
                    self.atmHybr[aId] = aL
                #elif aL < 2:
                #    self.atmHybr[aId] = 1
                elif aAtm['_chem_comp_atom.type_symbol']=="C" or \
                     aAtm['_chem_comp_atom.type_symbol']=="SI" or\
                     aAtm['_chem_comp_atom.type_symbol']=="Si" or\
                     aAtm['_chem_comp_atom.type_symbol']=="GE" or\
                     aAtm['_chem_comp_atom.type_symbol']=="Ge":
                     if aL==4:
                         self.atmHybr[aId] = 3
                     elif aL==3:
                         if aM==1:
                             self.atmHybr[aId] = 3
                         else:
                             self.atmHybr[aId] = 2
                     elif aL==2 :
                         if aM==1:
                             self.atmHybr[aId] = 2
                         else:
                             self.atmHybr[aId] = 1
                     elif aL==1:
                         self.atmHybr[aId] = 1
                elif aAtm['_chem_comp_atom.type_symbol']=="N" or \
                     aAtm['_chem_comp_atom.type_symbol']=="AS" or\
                     aAtm['_chem_comp_atom.type_symbol']=="As":
                     if aL==4 or aL==3:
                         self.atmHybr[aId] = 3
                     elif aL==2:
                         if aId in self.connMAMap:
                             self.atmHybr[aId] = 3
                         else:
                             self.atmHybr[aId] = 1
                     elif aL==1:
                         self.atmHybr[aId] = 1
                elif aAtm['_chem_comp_atom.type_symbol']=="B":
                     if aL==4:
                         self.atmHybr[aId] = 3
                     elif aL==3:
                         self.atmHybr[aId] = 2
                     elif aL==2:
                         if int(aAtm['_chem_comp_atom.charge']) ==1:
                             self.atmHybr[aId] = 1  
                         else:
                             self.atmHybr[aId] = 2
                     elif aL==1:
                         self.atmHybr[aId] = 1
                elif aAtm['_chem_comp_atom.type_symbol']=="O":
                     if aL==2 or aL ==1:
                         self.atmHybr[aId] = 3
                elif aAtm['_chem_comp_atom.type_symbol']=="P":
                     if aL > 1 and aL <6:
                         self.atmHybr[aId] = 3     
                elif aAtm['_chem_comp_atom.type_symbol']=="S":
                     if aL==4 or aL==3 or aL==2 or aL==1:
                         self.atmHybr[aId] = 3
                     elif aL==6:
                         self.atmHybr[aId] = 5  
                elif aAtm['_chem_comp_atom.type_symbol']=="SE" or\
                    aAtm['_chem_comp_atom.type_symbol']=="Se":
                     if aL==3:
                         self.atmHybr[aId] = 2
                     else:
                         self.atmHybr[aId] = 3
                
                         
        # second round
        for aAtm in tAtoms:
            if aAtm['_chem_comp_atom.type_symbol']=="N":
                aId = aAtm['_chem_comp_atom.atom_id']
                lSp = False
                for aNM in self.nonMAtmConnsMap[aId]:
                    if self.atmHybr[aNM]== 2:
                        lSp=True
                        break
                if lSp:
                    self.atmHybr[aId]=2
                
                        
                
                         
    def setASpeAng(self, tMA, tMN, tNN, tSP, tAngSum):
        
        aAng = {}
        aAng["_chem_comp_angle.atom_id_1"] = tMA 
        aAng["_chem_comp_angle.atom_id_2"] = tMN
        aAng["_chem_comp_angle.atom_id_3"] = tNN
        if tSP==1:
            aAng["_chem_comp_angle.value_angle"] = "180.00"
        elif tSP==2:
            aVal = (360.0-tAngSum[tMN])/2.0
            aSV  = "%8.4f"%aVal
            aAng["_chem_comp_angle.value_angle"] = aSV
        elif tSP==3:
            aAng["_chem_comp_angle.value_angle"] = "109.47"
        else:
            aAng["_chem_comp_angle.value_angle"] = "0.0"
        
        aAng["_chem_comp_angle.value_angle_esd"] = "5.0"
        #print("add metal related angle:", aAng)
    
        self.speAngs.append(aAng)
        
    def setMetalPA(self):
        
        rAndPMap1 = {}
        rAndPMap2 = {} 
        
        """
        for aR in self.simpR:
            rAndPMap1[aR] = {}
            for aOA in self.simpR[aR]: 
                for aP in self.simpP:
                    if aOA in self.simpP[aP]:
                        if not aP in rAndPMap1[aR]:
                            rAndPMap1[aR][aP]=1
                        else:
                            rAndPMap1[aR][aP]=rAndPMap1[aR][aP]+1
        for aR in rAndPMap1:
            for aP in rAndPMap1[aR]:
                if rAndPMap1[aR][aP]== len(self.simpR[aR]):
                    rAndPMap2[aR]=aP
                    
        print(rAndPMap2)
        """
            
        for aMA in self.metalConnAtomsMap:
            #print("a MA ", aMA)
            for aOA in self.metalConnAtomsMap[aMA]:
                #print(" a OA ", aOA)
                for aP in self.simpP:
                    #print("aP ", aP)
                    #print(self.simpP[aP])
                    if aOA in self.simpP[aP]:
                        #print(aOA, ' is in ', aP)
                        if not aMA in self.metalInPs:
                            self.metalInPs[aMA] = []
                        if not aP in self.metalInPs[aMA]:
                            self.metalInPs[aMA].append(aP)
        
        #print(self.metalInPs)
        
    def getAtomById(self, tAtoms, tId):
        
        aRetA = None 
        
        for aA in tAtoms:
            if aA['_chem_comp_atom.atom_id'] == tId:
                aRetA = aA
                break
        
        return aRetA
            
        
    def checkVal(self, tAtom, tDefBo):
        
        retH = 0
        
        numBO = 0
        aId = tAtom['_chem_comp_atom.atom_id']
        #print(" Here For atom ", aId)
        if aId in self.mcAtomBonds:
            for aB in self.mcAtomBonds[aId]:
                if "_chem_comp_bond.type" in aB.keys():
                    numBO+=BondOrderS2N(aB['_chem_comp_bond.type']) 
                elif '_chem_comp_bond.value_order' in aB.keys():
                    #print("Bond order is ", aB['_chem_comp_bond.value_order'])
                    numBO+=BondOrderS2N(aB['_chem_comp_bond.value_order']) 
                    #print("numBO is ", numBO)
        #print("def val is ", self.defaultBo[tAtom['_chem_comp_atom.type_symbol']])     
        #print("charge is ", int(tAtom['_chem_comp_atom.charge']))
            retH = tDefBo[tAtom['_chem_comp_atom.type_symbol']]\
                   -numBO + int(tAtom['_chem_comp_atom.charge'])
            #print("number H need ", retH)
        else :
            retH = tDefBo[tAtom['_chem_comp_atom.type_symbol']]
        return retH
    
    def setCCP4MonDataDescritor(self, tDataDescriptor):
        
        # Get a CCP4 monomer lib data descriptor
    
        s1 = self.monomRoot
        s2 = self.monomRoot
        s3 =".         "
        s4 ="non-polymer"
        s4 ="non-polymer"
        s5 = str(len(self.remainAtoms) + len(self.addedHs))
        numHEAVA =0
        for aA in self.remainAtoms:
            if aA['_chem_comp_atom.type_symbol']=="H":
                 numHEAVA +=1
                 
        s6 = str(numHEAVA)
        s7 ="."
        
        for aKey in sorted(tDataDescriptor.keys()):
            if tDataDescriptor[aKey][0].find("_chem_comp.id") !=-1:
                s1 = tDataDescriptor[aKey][1]  
                s2 = tDataDescriptor[aKey][1]  
            elif tDataDescriptor[aKey][0].find("_chem_comp.name") !=-1:
                temStrs = tDataDescriptor[aKey][1].strip().split()
                tmpDD   = tDataDescriptor[aKey][1].strip()
                if len(temStrs)>1 and tmpDD[0] !="\"" and tmpDD[-1] !="\"":
                    s3 = "\"" + tDataDescriptor[aKey][1] + "\""
                else:
                    s3 = tDataDescriptor[aKey][1]
            elif tDataDescriptor[aKey][0].find("_chem_comp.group") !=-1:
                s4 = tDataDescriptor[aKey][1]  
            elif tDataDescriptor[aKey][0].find("_chem_comp.type") !=-1:
                
                if tDataDescriptor[aKey][1].upper().find("L-PEPTIDE") !=-1 :
                    s4 = "L-PEPTIDE"
                elif tDataDescriptor[aKey][1].upper().find("D-PEPTIDE") !=-1 :
                    s4 = "D-PEPTIDE"
                elif tDataDescriptor[aKey][1].upper().find("M-PEPTIDE") !=-1 :
                    s4 = "M-PEPTIDE"
                elif tDataDescriptor[aKey][1].upper().find("PEPTIDE-") !=-1 :
                    s4 = "PEPTIDE"
                elif tDataDescriptor[aKey][1].upper().find("NON-POLYMER") !=-1 :
                    s4 = "NON-POLYMER"
                elif tDataDescriptor[aKey][1].upper().find("DNA ") !=-1 :
                    s4 = "DNA"
                elif tDataDescriptor[aKey][1].upper().find("RNA ") !=-1 :
                    s4 = "RNA"
                else:
                    s4 = tDataDescriptor[aKey][1]  
    
        aLine = "%s%s%s%s%s%s%s"%(s1.ljust(len(s1)+5), s2.ljust(len(s2)+5), \
                            s3.ljust(len(s3)+5), s4.ljust(len(s4)+5), \
                            s5.ljust(len(s5)+5), s6.ljust(len(s6)+5), \
                            s7.ljust(len(s7)+5))
            
        self.dataDiescriptor.append(aLine)
        
    def outInitModifiedCif(self, tMmcifName):
        
        try:
            aMmCif = open(tMmcifName, "w")
        except IOError:
            print(tMmcifName, " Could not be opened for writing")
        else:
            # Header section

            aMmCif.write("global_\n")
            aMmCif.write("_lib_name         ?\n")
            aMmCif.write("_lib_version      ?\n")
            aMmCif.write("_lib_update       ?\n")
            aMmCif.write(
                "# ------------------------------------------------\n")
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
            aMmCif.write(self.dataDiescriptor[0]+"\n")
            
            aLigId = ""
            if len(self.monomRoot) < 3:
                aLigId = self.monomRoot
            else:
                aLigId = self.monomRoot[:3]                
            aMmCif.write("data_comp_%s\n"%aLigId)
            # Atom section
            aMmCif.write("loop_\n")
            aMmCif.write("_chem_comp_atom.comp_id\n")
            aMmCif.write("_chem_comp_atom.atom_id\n")
            aMmCif.write("_chem_comp_atom.type_symbol\n")
            aMmCif.write("_chem_comp_atom.charge\n")
            aMmCif.write("_chem_comp_atom.x\n")
            aMmCif.write("_chem_comp_atom.y\n")
            aMmCif.write("_chem_comp_atom.z\n")
            
            for aAtom in self.remainAtoms:
                posX =0.0
                posY =0.0
                posZ =0.0
                
                if "_chem_comp_atom.model_Cartn_x" in aAtom.keys() and \
                    aAtom['_chem_comp_atom.model_Cartn_x'].find("?")==-1:
                    posX = float(aAtom['_chem_comp_atom.model_Cartn_x'])
                    posY = float(aAtom['_chem_comp_atom.model_Cartn_y'])
                    posZ = float(aAtom['_chem_comp_atom.model_Cartn_z'])
                elif '_chem_comp_atom.x' in aAtom.keys():
                    posX = float(aAtom['_chem_comp_atom.x'])
                    posY = float(aAtom['_chem_comp_atom.y'])
                    posZ = float(aAtom['_chem_comp_atom.z'])
                
                aMmCif.write("%s%s%s%10.2f%10.4f%10.4f%10.4f\n"
                             % (self.monomRoot.ljust(8), aAtom['_chem_comp_atom.atom_id'].ljust(10),
                                aAtom['_chem_comp_atom.type_symbol'].ljust(6), 
                                float(aAtom['_chem_comp_atom.charge']), posX, posY, posZ))
                
            for aAtom in self.addedHs:
                posX =0.0
                posY =0.0
                posZ =0.0
                
                if "_chem_comp_atom.model_Cartn_x" in aAtom.keys() and \
                    aAtom['_chem_comp_atom.model_Cartn_x'].find("?")==-1:
                    posX = float(aAtom['_chem_comp_atom.model_Cartn_x'])
                    posY = float(aAtom['_chem_comp_atom.model_Cartn_y'])
                    posZ = float(aAtom['_chem_comp_atom.model_Cartn_z'])
                elif '_chem_comp_atom.x' in aAtom.keys():
                    posX = float(aAtom['_chem_comp_atom.x'])
                    posY = float(aAtom['_chem_comp_atom.y'])
                    posZ = float(aAtom['_chem_comp_atom.z'])
                
                aMmCif.write("%s%s%s%10.2f%10.4f%10.4f%10.4f\n"
                             % (self.monomRoot.ljust(8), aAtom['_chem_comp_atom.atom_id'].ljust(10),
                                aAtom['_chem_comp_atom.type_symbol'].ljust(6), 
                                float(aAtom['_chem_comp_atom.charge']), posX, posY, posZ)) 
                
                
            # Bond section
            if len(self.remainBonds) or len(self.addedHBs):
                aMmCif.write("#\n")
                aMmCif.write("loop_\n")
                aMmCif.write("_chem_comp_bond.comp_id\n")
                aMmCif.write("_chem_comp_bond.atom_id_1\n")
                aMmCif.write("_chem_comp_bond.atom_id_2\n")
                aMmCif.write("_chem_comp_bond.type\n")
                aMmCif.write("_chem_comp_bond.value_dist\n")
                aMmCif.write("_chem_comp_bond.value_dist_esd\n")
                
                for aBond in self.remainBonds:
                    aType = ""
                    if '_chem_comp_bond.value_order' in aBond.keys():
                        aType = aBond['_chem_comp_bond.value_order']
                    elif '_chem_comp_bond.type' in aBond.keys():
                        aType = aBond['_chem_comp_bond.type']
                    else:
                        print("No bond-order for the bond between atoms %s and %s"
                              %(aBond['_chem_comp_bond.atom_id_1'],
                                aBond['_chem_comp_bond.atom_id_2']))
                        sys.exit()
                    
                    bLen = "0.0"
                    dBLen = "0.01"
                    aMmCif.write("%s%s%s%s%s%s\n"
                                 % (self.monomRoot.ljust(8),
                                 aBond['_chem_comp_bond.atom_id_1'].ljust(10), 
                                 aBond['_chem_comp_bond.atom_id_2'].ljust(10),  
                                 aType.ljust(10), bLen.ljust(10), dBLen.ljust(10)))
            
            
                for aBond in self.addedHBs:
                    aType = ""
                    if '_chem_comp_bond.value_order' in aBond.keys():
                        aType = aBond['_chem_comp_bond.value_order']
                    elif '_chem_comp_bond.type' in aBond.keys():
                        aType = aBond['_chem_comp_bond.type']
                    else:
                        print("No bond-order for the bond between atoms %s and %s"
                              %(aBond['_chem_comp_bond.atom_id_1'],
                                aBond['_chem_comp_bond.atom_id_2']))
                        sys.exit()
                    
                    bLen = "0.0"
                    dBLen = "0.01"
                    aMmCif.write("%s%s%s%s%s%s\n"
                                 % (self.monomRoot.ljust(8),
                                 aBond['_chem_comp_bond.atom_id_1'].ljust(10), 
                                 aBond['_chem_comp_bond.atom_id_2'].ljust(10),  
                                 aType.ljust(10), bLen.ljust(10), dBLen.ljust(10)))
            
            aMmCif.close()
    
 
    
    def outInterMCif(self, tMmcifName, tAtoms, tBonds, tVersionInfo):
        
        try:
            aMmCif = open(tMmcifName, "w")
        except IOError:
            print(tMmcifName, " Could not be opened for writing")
        else:
            # Header section

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
            aMmCif.write(self.dataDiescriptor[0]+"\n")
            
            aLigId = ""
            if len(self.monomRoot) < 3:
                aLigId = self.monomRoot
            else:
                aLigId = self.monomRoot[:3]                
            aMmCif.write("data_comp_%s\n"%aLigId)
            # Atom section
            aMmCif.write("loop_\n")
            aMmCif.write("_chem_comp_atom.comp_id\n")
            aMmCif.write("_chem_comp_atom.atom_id\n")
            aMmCif.write("_chem_comp_atom.type_symbol\n")
            aMmCif.write("_chem_comp_atom.charge\n")
            aMmCif.write("_chem_comp_atom.model_Cartn_x\n")
            aMmCif.write("_chem_comp_atom.model_Cartn_y\n")
            aMmCif.write("_chem_comp_atom.model_Cartn_z\n")
            
            for aAtom in tAtoms:
                posX =0.0
                posY =0.0
                posZ =0.0
                
                if "_chem_comp_atom.model_Cartn_x" in aAtom.keys() and \
                    aAtom['_chem_comp_atom.model_Cartn_x'].find("?")==-1:
                    posX = float(aAtom['_chem_comp_atom.model_Cartn_x'])
                    posY = float(aAtom['_chem_comp_atom.model_Cartn_y'])
                    posZ = float(aAtom['_chem_comp_atom.model_Cartn_z'])
                elif '_chem_comp_atom.x' in aAtom.keys():
                    posX = float(aAtom['_chem_comp_atom.x'])
                    posY = float(aAtom['_chem_comp_atom.y'])
                    posZ = float(aAtom['_chem_comp_atom.z'])
                
                aMmCif.write("%s%s%s%10.2f%10.4f%10.4f%10.4f\n"
                             % (self.monomRoot.ljust(8), aAtom['_chem_comp_atom.atom_id'].ljust(10),
                                aAtom['_chem_comp_atom.type_symbol'].ljust(6), 
                                float(aAtom['_chem_comp_atom.charge']), posX, posY, posZ))
                
                
            # Bond section
            if len(tBonds):
                aMmCif.write("#\n")
                aMmCif.write("loop_\n")
                aMmCif.write("_chem_comp_bond.comp_id\n")
                aMmCif.write("_chem_comp_bond.atom_id_1\n")
                aMmCif.write("_chem_comp_bond.atom_id_2\n")
                aMmCif.write("_chem_comp_bond.type\n")
                aMmCif.write("_chem_comp_bond.value_dist_nucleus\n")
                aMmCif.write("_chem_comp_bond.value_dist_nucleus_esd\n")
                aMmCif.write("_chem_comp_bond.value_dist\n")
                aMmCif.write("_chem_comp_bond.value_dist_esd\n")
                
                for aBond in tBonds:
                    aType = ""
                    if '_chem_comp_bond.value_order' in aBond.keys():
                        aType = aBond['_chem_comp_bond.value_order']
                    elif '_chem_comp_bond.type' in aBond.keys():
                        aType = aBond['_chem_comp_bond.type']
                    else:
                        print("No bond-order for the bond between atoms %s and %s"
                              %(aBond['_chem_comp_bond.atom_id_1'],
                                aBond['_chem_comp_bond.atom_id_2']))
                        sys.exit()
                    
                    bLen = "%s"%aBond["_chem_comp_bond.value_dist"]
                    dBLen = "0.04"
                    aMmCif.write("%s%s%s%s%s%s%s%s\n"
                                 % (self.monomRoot.ljust(8),
                                 aBond['_chem_comp_bond.atom_id_1'].ljust(10), 
                                 aBond['_chem_comp_bond.atom_id_2'].ljust(10),  
                                 aType.ljust(10), bLen.ljust(10), dBLen.ljust(10), 
                                 bLen.ljust(10), dBLen.ljust(10)))
            if len(tVersionInfo):
                aMmCif.write("loop_\n")
                aMmCif.write("_acedrg_chem_comp_descriptor.comp_id\n")
                aMmCif.write("_acedrg_chem_comp_descriptor.program_name\n")
                aMmCif.write("_acedrg_chem_comp_descriptor.program_version\n")
                aMmCif.write("_acedrg_chem_comp_descriptor.type\n")
                if 'ACEDRG_VERSION' in tVersionInfo:
                    aL = "%s%s%s%s\n"%(self.monomRoot.ljust(8), "acedrg".ljust(21),
                                       tVersionInfo['ACEDRG_VERSION'].ljust(21),
                                       '\"dictionary generator\"'.ljust(40))
                    aMmCif.write(aL)
                if "DATABASE_VERSION" in tVersionInfo:
                    aL = "%s%s%s%s\n"%(self.monomRoot.ljust(8), "acedrg_database".ljust(21),
                                       tVersionInfo['DATABASE_VERSION'].ljust(21),
                                       '\"data source\"'.ljust(40))
                    aMmCif.write(aL)   
                if "RDKit_VERSION" in tVersionInfo:
                    aL = "%s%s%s%s\n"%(self.monomRoot.ljust(8), "rdkit".ljust(21),
                                       tVersionInfo['RDKit_VERSION'].ljust(21),
                                       '\"Chemoinformatics tool\"'.ljust(40))                     
                    aMmCif.write(aL)
                if "SERVALCAT_VERSION" in tVersionInfo:
                    aL = "%s%s%s%s\n"%(self.monomRoot.ljust(8), "servalcat".ljust(21),
                                       tVersionInfo['SERVALCAT_VERSION'].ljust(21),
                                       '\"optimization tool\"'.ljust(40))
            
            aMmCif.close()        
            
    def setMetalNBAngs(self, tAtoms):
        
        self.setHYAtoms(self, tAtoms)
        
        for aMA in self.metalConnAtomsMap.keys():
            for aMNB in self.metalConnAtomsMap[aMA]:
                pass
        
    def setHYAtoms(self, tAtoms):
        
        pass 
        
    def modiCif(self, tInCif, tOutCif):
        
        
        hIds = []
        for aHA in self.addedHs:
            hIds.append(aHA['_chem_comp_atom.atom_id'])
            
        
        for aA in self.metalConnFullAtoms:
            #numC   = -self.checkVal(aA)
            #aA['_chem_comp_atom.charge'] = numC
            #if numC !=0:
            print(aA['_chem_comp_atom.atom_id'], " carries ", aA['_chem_comp_atom.charge'], " charge")
        print("number of metal conn atoms ", len(self.metalConnFullAtoms))
        for aMA in self.metalAtoms:
            print("Metal Atom ", aMA['_chem_comp_atom.atom_id'])
            numC =0
            for aA in self.metalConnFullAtoms:
                #print( aA['_chem_comp_atom.atom_id'])
                #print(self.metalConnAtomsMap)
                if aMA['_chem_comp_atom.atom_id'] in self.metalConnAtomsMap and\
                    aA['_chem_comp_atom.atom_id'] in self.metalConnAtomsMap[aMA['_chem_comp_atom.atom_id']]:
                      numC-=aA['_chem_comp_atom.charge']
            aMA['_chem_comp_atom.charge'] = numC
            print("Its charge ", aMA['_chem_comp_atom.charge'])
            
        
        
        
        # Modify cif 
        try:
            aMmCif = open(tInCif, "r")
        except IOError:
            print(tInCif, " Could not be opened for writing")
        else:
            allLs = aMmCif.readlines()
            aMmCif.close()
            
            
            lExistAng = False
            for aL in allLs:
                if aL.find("_chem_comp_angle.") !=-1:
                    lExistAng = True
                    break
                
            lSP = False
            lSR = False
            
            for aL in allLs:
                if aL.find('loop_') != -1:
                    lSP = False
                    lSR = False
                elif lSP or lSR:
                    strGrp = aL.strip().split()
                    if lSP :
                        if len(strGrp)==4:
                            idP  = strGrp[1]
                            idPA = strGrp[2]
                            if not idP in self.simpP:
                                self.simpP[idP] = []
                            self.simpP[idP].append(idPA)
                    elif lSR :
                        if len(strGrp)==4:
                            idR  = strGrp[1]
                            idRA = strGrp[2]
                            idAR = strGrp[3]
                            if idAR.find("YES") !=-1:
                                if not idR in self.simpR:
                                    self.simpR[idR] = []
                                self.simpR[idR].append(idRA)
                if aL.find('_chem_comp_plane_atom.dist_esd') !=-1:
                    lSP = True
                    lSR = False
                elif aL.find('_chem_comp_ring_atom.is_aromatic_ring') !=-1:
                    lSR = True
                    lSP = False
                    
            #print(self.simpP)
            #print(self.simpR)
    
    
            
            aFSYS = FileTransformer()
            aFSYS.mmCifReader(tInCif)
            
            angSumMap = {}
            for aAng in aFSYS.angles: 
                idCen = aAng["_chem_comp_angle.atom_id_2"] 
                if not idCen in angSumMap.keys():
                    angSumMap[idCen]=float(aAng['_chem_comp_angle.value_angle'])
                else:
                    angSumMap[idCen]+=float(aAng['_chem_comp_angle.value_angle'])
            
            #for aA in angSumMap:
            #    print(" ang sum for ", aA, " is ", angSumMap[aA])
            
            for aMA in self.metalConnAtomsMap.keys():
                #print("For metal atom ", aMA)
                for aMN in self.metalConnAtomsMap[aMA]:
                    aSP = self.atmHybr[aMN]
                    #print("atom ", aMN, " hybr ", aSP)
                    #if aMN in self.atmNonHMap.keys():
                    for aNN in self.nonMAtmConnsMap[aMN]:
                        #print("NB atom ", aMN, " has the following angles: ")
                        #print("Angle among %s and %s and %s"%(aMA, aMN, aNN))
                        self.setASpeAng(aMA, aMN, aNN, aSP, angSumMap)
            
            self.setMetalPA()
            
            
            
            lA = False
            lB = False
            lANG = False
            lP   = False
            lDis = False
            newLines = []
            for aL in allLs:
                if lA and aL.find('_chem_comp_atom') == -1:
                    self.writeNewAtomCifLine(newLines)
                    lA=False
                elif lB and aL.find('_chem_comp_bond')==-1:
                    self.writeNewBondCifLine(newLines)
                    lB=False
                elif lANG and aL.find('_chem_comp_angle.')==-1:
                    if lExistAng:
                        self.writeNewAngCifLine(newLines)
                    lANG=False
                elif lP and aL.find('_chem_comp_plane_atom.')==-1:
                    self.writeNewPlCifLine(newLines)
                    lP=False
                elif not  lA and aL.find('_chem_comp_atom') != -1:
                    lA=True
                    lB=False
                    lANG = False
                    lP=False
                elif not lB and aL.find('_chem_comp_bond') != -1:
                    lB=True
                    lA=False
                    lANG=False
                    lP=False
                elif not lANG and aL.find('_chem_comp_angle') != -1:
                    lANG=True
                    lB=False
                    lA=False
                    lP=False
                elif not lP and aL.find('_chem_comp_plane_atom') !=-1:
                    lP = True
                    lA = False
                    lB = False
                    lANG = False
                elif not lDis and aL.find('_acedrg_chem_comp_descriptor.') !=-1:
                    if not lExistAng:
                        self.writeNewAngCifLine2(newLines)
                    lDis = True
                    lP = False
                    lA = False
                    lB = False
                    lANG = False
                    
                if self.noNewHLine(aL, hIds):
                    newLines.append(aL)
                else:
                    pass 
                #else:
                #    print(aL)
                
            """
            lA = False
            newLines2 = []

            for aL2 in newLines:
                if lA and aL2.find('loop_')!=-1:
                    lA = False
                    newLines2.append(aL2)
                elif lA and aL2.find('_chem_comp_atom') == -1:
                    [lC, aCA] =self.chargeL(aL2, cIds)
                    if lC:
                        self.replChargeLine(aL2, aCA, newLines2)
                    else:
                        newLines2.append(aL2)
                elif not lA and aL2.find('_chem_comp_atom.pdbx_model_Cartn_z_ideal') !=-1:
                    lA = True
                    newLines2.append(aL2)
                else:
                    newLines2.append(aL2)
                    
            """
            
            try:
                aOutCif = open(tOutCif, "w")
            except IOError:
                print(tOutCif, " Could not be opened for writing")
            else:
                for aL in newLines:
                    aOutCif.write(aL)
                aOutCif.close()
        
   
    def noNewHLine(self, tLine, tHIds):
        
        aRet = True
        
        for aHId in tHIds:
            if tLine.find(aHId) !=-1:
                aRet = False
                break
        
        return aRet
        
    def chargeL(self, tL, tCA):
        
        aRet = False
        aRetA = None
        
        for aA in tCA.keys():
            if tL.find(aA) !=-1:
                strs = tL.strip().split()
                if len(strs) > 2:
                    if strs[1].find(aA) !=-1:
                        aRet = True
                        aRetA = tCA[aA]
                        break
        return [aRet, aRetA]
    
    def replChargeLine(self, tL, tCA, tLines):
        
        strs = tL.strip().split()
        aNewL = "%s%s%s%s%s%s%s%s%s%s%s%s\n"%(strs[0].ljust(8), strs[1].ljust(8), strs[2].ljust(8),
                                            strs[3].ljust(8), strs[4].ljust(8), str(tCA['_chem_comp_atom.charge']).ljust(8),
                                            strs[6].ljust(12), strs[7].ljust(12), strs[8].ljust(12), 
                                            strs[9].ljust(12), strs[10].ljust(12),
                                            strs[11].ljust(12))
        
        tLines.append(aNewL)
        
        
    def writeNewAtomCifLine(self, tLines):
        
        
            
        for aMA in self.metalAtoms:
            posX =0.0
            posY =0.0
            posZ =0.0
            
            if "_chem_comp_atom.model_Cartn_x" in aMA.keys():
                if  aMA['_chem_comp_atom.model_Cartn_x'].find("?")==-1:
                     posX = float(aMA['_chem_comp_atom.model_Cartn_x'])
                     posY = float(aMA['_chem_comp_atom.model_Cartn_y'])
                     posZ = float(aMA['_chem_comp_atom.model_Cartn_z'])
                elif '_chem_comp_atom.pdbx_model_Cartn_x_ideal'  in aMA.keys()\
                    and aMA['_chem_comp_atom.pdbx_model_Cartn_x_ideal'].find("?")==-1:
                        posX = float(aMA['_chem_comp_atom.pdbx_model_Cartn_x_ideal'])
                        posY = float(aMA['_chem_comp_atom.pdbx_model_Cartn_y_ideal'])
                        posZ = float(aMA['_chem_comp_atom.pdbx_model_Cartn_z_ideal'])
            elif '_chem_comp_atom.x' in aMA.keys():
                posX = float(aMA['_chem_comp_atom.x'])
                posY = float(aMA['_chem_comp_atom.y'])
                posZ = float(aMA['_chem_comp_atom.z'])
            aLine = "%s%s%s%s%s%10.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n"\
                    % (self.monomRoot.ljust(8), aMA['_chem_comp_atom.atom_id'].ljust(10),
                       aMA['_chem_comp_atom.atom_id'].ljust(10),
                       aMA['_chem_comp_atom.type_symbol'].ljust(6), 
                       aMA['_chem_comp_atom.type_symbol'].upper().ljust(6),
                       float(aMA['_chem_comp_atom.charge']), posX, posY, posZ, posX, posY, posZ)
            tLines.append(aLine)
            
    def writeNewBondCifLine(self, tLines):
        
        for aB in self.metalBonds:
            aType = ""
            if '_chem_comp_bond.value_order' in aB.keys():
                aType = aB['_chem_comp_bond.value_order']
            elif '_chem_comp_bond.type' in aB.keys():
                aType = aB['_chem_comp_bond.type']
            aLine = "%s%s%s%s%s%s%s%s%s\n"\
                  % (self.monomRoot.ljust(8),
                     aB['_chem_comp_bond.atom_id_1'].ljust(10), 
                     aB['_chem_comp_bond.atom_id_2'].ljust(10),  
                     aType.ljust(10), "n".ljust(8), aB["_chem_comp_bond.value_dist"].ljust(10), 
                     "0.04".ljust(10), aB["_chem_comp_bond.value_dist"].ljust(10), "0.04".ljust(10))
            tLines.append(aLine)
    
    def writeNewAngCifLine(self, tLines):
        
        for aA in self.speAngs:
        
            aLine = "%s%s%s%s%s%s\n"\
                  % (self.monomRoot.ljust(8),
                     aA['_chem_comp_angle.atom_id_1'].ljust(10), 
                     aA['_chem_comp_angle.atom_id_2'].ljust(10),
                     aA['_chem_comp_angle.atom_id_3'].ljust(10),
                     aA["_chem_comp_angle.value_angle"].ljust(16), 
                     aA["_chem_comp_angle.value_angle_esd"].ljust(10))
            tLines.append(aLine)
            
    def writeNewAngCifLine2(self, tLines):
        
        tLines.append("_chem_comp_angle.comp_id\n")
        tLines.append("_chem_comp_angle.atom_id_1\n")
        tLines.append("_chem_comp_angle.atom_id_2\n")
        tLines.append("_chem_comp_angle.atom_id_3\n")
        tLines.append("_chem_comp_angle.value_angle\n")
        tLines.append("_chem_comp_angle.value_angle_esd\n")
        
        for aA in self.speAngs:
            aLine = "%s%s%s%s%s%s\n"\
                  % (self.monomRoot.ljust(8),
                     aA['_chem_comp_angle.atom_id_1'].ljust(10), 
                     aA['_chem_comp_angle.atom_id_2'].ljust(10),
                     aA['_chem_comp_angle.atom_id_3'].ljust(10),
                     aA["_chem_comp_angle.value_angle"].ljust(16), 
                     aA["_chem_comp_angle.value_angle_esd"].ljust(10))
            tLines.append(aLine)
        tLines.append("loop_\n")
    
    def writeNewPlCifLine(self, tLines):
        
        pl_esd = "0.020"
        for aMA in self.metalInPs:
            for aPl in self.metalInPs[aMA]:
                aLine = "%s%s%s%s\n"\
                  % (self.monomRoot.ljust(8),
                     aPl.ljust(10), 
                     aMA.ljust(10),
                     pl_esd.ljust(10))
                    
            tLines.append(aLine)
        
    def writeNewAtomPdbLines(self, tLines):
        
        idx = len(self.metalAtoms)
        
        for aA in self.metalAtoms:
            sIdx = str(idx)
            posX =0.0
            posY =0.0
            posZ =0.0
            
            if "_chem_comp_atom.model_Cartn_x" in aA.keys():
                if aA['_chem_comp_atom.model_Cartn_x'].find("?")==-1:
                    posX = aA['_chem_comp_atom.model_Cartn_x']
                    posY = aA['_chem_comp_atom.model_Cartn_y']
                    posZ = aA['_chem_comp_atom.model_Cartn_z']
                elif '_chem_comp_atom.pdbx_model_Cartn_x_ideal'  in aA.keys()\
                     and aA['_chem_comp_atom.pdbx_model_Cartn_x_ideal'].find("?")==-1:
                    posX = aA['_chem_comp_atom.pdbx_model_Cartn_x_ideal']
                    posY = aA['_chem_comp_atom.pdbx_model_Cartn_y_ideal']
                    posZ = aA['_chem_comp_atom.pdbx_model_Cartn_z_ideal']    
            elif '_chem_comp_atom.x' in aA.keys():
                posX = aA['_chem_comp_atom.x']
                posY = aA['_chem_comp_atom.y']
                posZ = aA['_chem_comp_atom.z']
            
            seqNum  = "X"
            empt0   = " "
            empt   = "         "
            charge = int(aA['_chem_comp_atom.charge'])
            sCharge = ""
            if charge > 0:
                sCharge = "+" + str(charge) 
            elif charge < 0:
                   sCharge = "0"+ str(abs(charge)) 
                
            aLine = "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%("HETATM".ljust(6), sIdx.rjust(5), 
                                                      aA['_chem_comp_atom.atom_id'].rjust(5),
                                                      self.monomRoot.rjust(4), seqNum.rjust(2), "1".rjust(4),empt0.rjust(4),
                                                            posX.rjust(8), posY.rjust(8),posZ.rjust(8),
                                                            "1.00".rjust(6), "0.50".rjust(6), empt.rjust(10), 
                                                            aA['_chem_comp_atom.type_symbol'].rjust(2), sCharge.rjust(2))
            tLines.append(aLine)
            idx+=1
            
    
    
    def runAcedrg(self, tInCif, tOutRoot):
        
        if not os.path.isdir(tOutRoot+"_TMP"):
            os.mkdir(tOutRoot+"_TMP")
        aMtConnFileName = os.path.join(tOutRoot+"_TMP", "mtconn.list")
        aMC = open(aMtConnFileName, "w")
        for aK in self.connMAMap.keys():
            for aM in self.connMAMap[aK]:
                aMC.write("BONDING  %s%s\n"%(aK.ljust(10), aM.ljust(10)))
        aMC.close()    
        
        self._cmdline   = "acedrg -c %s  -r %s -o %s --mtList %s "%(tInCif, self.monomRoot, tOutRoot, aMtConnFileName)
        self._log_name  = tOutRoot + "_acedrg.log"
        #print(self._cmdline)
        self.runExitCode = self.subExecute()
    
    def runRefmac(self, tInCif, tInPdb, tOutRoot):
        
        
        aOutPdb  = tOutRoot + ".pdb"
        self._log_name = tOutRoot + "_refmac.log"
        self._cmdline  = "refmac5 xyzin %s libin %s xyzout %s"%(tInPdb, tInCif,  aOutPdb)
        self._cmdline += " <<eof\n"
        self._cmdline += "ncyc 40   \n"
        self._cmdline += "make hout yes    \n"
        self._cmdline += "end       \n"
        self._cmdline += "eof       \n"
        
        self.runExitCode = self.subExecute()
        
   
    

        
    
    def executeOld(self):
        
        print("workMode is ", self.workMode)
        
        if self.workMode == 1:
            if os.path.isfile(self.inMmCifName) :
               print("inName : ", self.inMmCifName)
               self.fileConv.mmCifReader(self.inMmCifName)
               print("Number of atoms in the cif file ", len(self.fileConv.atoms))
               
               if len(self.fileConv.atoms) > 0:
               
                   if len(self.fileConv.dataDescriptor):
                       self.setMonoRoot(self.fileConv.dataDescriptor)
                   else:
                       self.monomRoot = "LIG"
                   print("Mon-root is ", self.monomRoot)
                   
                   self.getNewMolWithoutMetal2(self.fileConv.atoms, 
                                               self.fileConv.bonds,
                                               self.outRstCifName)
                   
                   
                   tmpMmcifName = self.monomRoot + "_tmpIn.cif"
                   
                   self.setCCP4MonDataDescritor()
                   self.outInitModifiedCif(tmpMmcifName)
                   
                   if os.path.isfile(tmpMmcifName):
                       print("acedrg is running")
                       outTmpRootName = self.outRoot + "_tmp"
                       self.runAcedrg(tmpMmcifName, outTmpRootName)
                       
                       outTmpCifName  =  outTmpRootName + ".cif"
                       outTmpPDBName  =  outTmpRootName + ".pdb"
                       
                       print("Stage 1 ")
                       for aA in self.metalAtoms:
                           print("Metal atom ", aA['_chem_comp_atom.atom_id'], 
                                 " carries ", aA['_chem_comp_atom.charge'], 
                                 " formal charge ")
                       if os.path.isfile(outTmpCifName) and os.path.isfile(outTmpPDBName):
                           
                           outTmp2RootName = self.outRoot   + "_tmp2"
                           outTmp2CifName  = outTmp2RootName + ".cif"
                           outTmp2PDBName  = outTmp2RootName + ".pdb"
                           
                           #self.modiCifAndPdb(outTmpCifName, outTmpPDBName, 
                           #                   outTmp2CifName, outTmp2PDBName)
                           self.modiCifAndPdb2(outTmpCifName, outTmpPDBName, 
                                              outTmp2CifName, outTmp2PDBName)
                           
                           print("Stage 3 ")
                           for aA in self.metalAtoms:
                               print("Metal atom ", aA['_chem_comp_atom.atom_id'], 
                                     " carries ", aA['_chem_comp_atom.charge'], 
                                     " formal charge ")
                           
                           if os.path.isfile(outTmp2CifName) and os.path.isfile(outTmp2PDBName):
                               print ("Refmac is running")
                               outFinRootName = self.outRoot   + "_final"
                               
                               self.runRefmac(outTmp2CifName, outTmp2PDBName, outFinRootName)
                               #os.system("cp %s %s"%(outTmp2CifName, outFinRootName+".cif"))
                               print("The final pdb is ", outFinRootName+".pdb") 
                               for aA in self.metalAtoms:
                                   print("Metal atom ", aA['_chem_comp_atom.atom_id'], 
                                         " carries ", aA['_chem_comp_atom.charge'], 
                                         " formal charge ")
                                         
                               
        
        elif self.workMode == 2:
            if os.path.isfile(self.inSmiName) :
                print("input file is %s"%self.inSmiName)
                self.initMols('smi', self.inSmiName, self.mols)
                idxMol = 0
                self.metalAtoms={}
                self.metalConnAtoms={}
                for aMol in self.mols:
                    self.metalAtoms[idxMol]      = []
                    self.metalConnAtoms[idxMol]  = []
                    self.modiMol(aMol, self.metalAtoms[idxMol],  self.metalConnAtoms[idxMol],
                            self.organicSec, self.defaultBo)
                    
                    
               

#def main():
#    hemTObj = hemT(sys.argv[1:])

#if __name__ == "__main__":
#    main()
