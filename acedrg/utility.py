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

from builtins import range
import os,os.path,sys
import platform
import glob,shutil
import re,string
from optparse import OptionParser 
import time
import math
import select
import random

import numpy    as np
import networkx as nx

from rdkit import Chem
from rdkit.Chem import AllChem

# Functions for numbers

def isInt(t_string):

    aRet = True

    t_string = t_string.strip()
    try:
        aInt = int(t_string)
    except ValueError:
        aRet = False

    return aRet    

def isFloat(t_string):

    aRet = True

    t_string = t_string.strip()
    try:
        aInt = float(t_string)
    except ValueError:
        aRet = False

    return aRet 

# Functions for list comparison

def listComp(a_list, b_list):

    if a_list[1] > b_list[1]:
        return 1
    elif a_list[1] < b_list[1]:
        return -1

    if len(a_list[0]) < len(b_list[0]):
        return 1
    elif len(a_list[0]) > len(b_list[0]):
        return -1

    if a_list[0] > b_list[0]:
        return 1
    elif a_list[0] == b_list[0]:
        return 0
    elif a_list[0] < b_list[0]:
        return -1

def listComp2(a_list, b_list):

    if a_list[0] > b_list[0]:
        return 1
    elif a_list[0] < b_list[0]:
        return -1
    else :
        return 0

def listCompDes(a_list, b_list):

    if int(a_list[1]) < int(b_list[1]):
        return 1
    elif int(a_list[1]) > int(b_list[1]):
        return -1
    else :
        return 0

def listCompAcd(a_list, b_list):

    if int(a_list[0]) > int(b_list[0]):
        return 1
    elif int(a_list[0]) < int(b_list[0]):
        return -1
    else :
        return 0

# Functions for set a bool dictionary
def setBoolDict(tKey, tVK, tDict, tVD):

    tDict[tKey] = tVK

    for aK in list(tDict.keys()):
        if aK != tKey:
            tDict[aK] = tVD

# Special line splitting function

def splitLineSpa(tLine):

    # split a line by empty spaces, but take account of "\"" and "\'"
    reStrs = []
    aL = tLine.strip()
    tS = "\'"
    tD = "\""
    nS = aL.find(tS)
    nD = aL.find(tD)
    if nS != -1 or nD !=-1:
        aSep = tD
        if nS != -1 and nD !=-1:
            if nS < nD:
                aSep = tS 
        elif nS !=-1:
            aSep = tS
        blocs = [] 
        aStr = ""
        lIn   = False
        for aC in aL:
            if aC == aSep:
                if len(aStr) !=0:
                    if lIn:
                        aStr += aC
                        aPair = [lIn, aStr]
                        blocs.append(aPair)
                        aStr = ""
                        lIn = False
                    else:
                        aPair = [lIn, aStr]
                        blocs.append(aPair)
                        aStr = ""
                        aStr += aC
                        lIn  = True 
                else:
                    aStr += aC
                    lIn  = True 
            else :
                aStr +=aC
        
        if len(aStr) !=0:
            aPair = [lIn, aStr]
            blocs.append(aPair)

        for aB in blocs:
            if not aB[0]:
                aStrs2 = aB[1].strip().split()
                for aStr in  aStrs2:
                    if len(aStr.strip()):
                        reStrs.append(aStr)
            else:
                reStrs.append(aB[1])           
                         
    else:
        allStrs = tLine.strip().split()
        for aStr in allStrs:
            reStrs.append(aStr)   

    return reStrs 

def countPrime(tL):

    lCount = False
    nCount = 0
    nCount = 0
    aSep = "\""

    for aC in tL.strip():
        if aC.find(aSep) !=-1:
            nCount +=1
   
    if nCount > 1 :
        lCount = True

    return lCount 

def aLineToAList2(tL, tList):

    tL = tL.strip()

    lP = False
     
    str1 = ""
    str2 = ""
    str3 = ""

    sep  = "\'"
    b   = 0
    e   = 0
    lP  = False
   
    pos =0
    for aS in tL:
        if lP:
           if aS.find(sep) !=-1:
              e = pos
        elif aS.find(sep) !=-1:
           lP = True
           b  = pos
           e  = pos    
        pos +=1
    if b==e:
        tList = aL.strip().split()
    else:
        strs1 = tL[0:b].split()
        str2  = tL[b:e+1]
        strs3 = tL[e+1:len(tL)].split()
    
        for aE in strs1:
            tList.append(aE)

        tList.append(str2)

        for aE in strs3:
            tList.append(aE)

        print(tList)


def aLineToAlist2(tL, tList):
       
            aTS = ""
            l0  = False
            l1  = False
            l2  = False   

            for aS in tL:
                #print(aS)
                if aS ==" ":
                    if not l1 and not l2 and len(aTS)!=0:
                        tList.append(aTS.strip()) 
                        #print("1 ",aTS)
                        l1 = False
                        l2 = False
                        aTS = ""
                    elif l1 and not l2 and len(aTS) !=0:
                        tList.append(aTS.strip()) 
                        #print("2 ",aTS)
                        aTS = ""
                        l1 = False
                        l2 = False
                    elif l2:
                        aTS = aTS + aS 
                elif aS.find("\"") !=-1:
                    if l2:
                        aTS = "\"" + aTS + "\""
                        tList.append(aTS.strip())
                        #print("3 ",aTS)
                        aTS = ""
                        l1  = False  
                        l2  = False
                    elif l1:
                        aTS = aTS + aS
                    else:
                        l2 = True
                        l1 = False
                elif aS.find("\'") !=-1:
                    if l1:
                        if not l2:
                            aTS = "\'" + aTS + "\'"
                            tList.append(aTS)
                            #print("4 ",aTS)
                            aTS = ""
                            l1  = False
                            l2  = False
                        else:
                            aTS = aTS + aS
                            l1  = True
                    else:
                        aTS = aTS + aS
                        l1 = True
                else:
                    aTS = aTS + aS
               
                    
            if aTS != "":
                tList.append(aTS.strip())    

def aLineToAlist(tL, tList):

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
                aTS = aTS + " "
                #print "\"" + aTS + "\""
            elif l0:
                if len(aTS):
                    tList.append(aTS.strip())
                    aTS = ""
                    l0  = False
            else:
                l0  = True
    if aTS != "":
        tList.append(aTS.strip())    

def splitLineSpa2(tLine):

    reStrs = []

    aSep1 = "\'"
    aSep2 = "\""
    if tLine.find(aSep1) !=-1 or tLine.find(aSep2) !=-1:
        set1 = tLine.strip().split(aSep1)
        n1   = len(set1)
        set2 = tLine.strip().split(aSep2)
        n2   = len(set2)
    
        if n1 >=n2:
            splitLineSpa3(set1, aSep1, reStrs)
        else:
            splitLineSpa3(set2, aSep2, reStrs)
    else:
        reStrs = tLine.strip().split()
         
    return reStrs

def splitLineSpa3(tSet, tSep, tStrs):

    n = len(tSet)
 
    if n > 3:
        aM = ""
        for i in range(1, n-1):
            aM +=tSet[i] + tSep  
        aM1 = aM[:-1]
        TB = tSet[0].strip().split()
        for aT in TB:
            tStrs.append(aT)
        tStrs.append(aM1)
        TE = tSet[-1].strip().split()
        for aT in TE:
            tStrs.append(aT)
    elif n ==3:
        TB = tSet[0].strip().split()
        for aT in TB:
            tStrs.append(aT)
        tStrs.append(tSet[1])
        TE = tSet[-1].strip().split()
        for aT in TE:
            tStrs.append(aT)

def setNameByNumPrime(tStr):

    aRet = tStr
    aL = len(tStr)
    pos1 =0
    pos2 =aL-1
    if len(tStr) > 0:
        if tStr.count("\'") ==1:
            if tStr[pos1] !="\"" and tStr[pos2] !="\"":
                aRet= "\"" + tStr + "\""
        elif tStr.count("\"") ==1:
            if tStr[pos1] !="\'" and tStr[pos2] !="\'":
                aRet= "\'" + tStr + "\'"
 
    return aRet

def BondOrderS2N(tBS):

    aBN = -1

    aBS = tBS.upper()
    if aBS.find("SING") !=-1:
        aBN = 1      
    elif aBS.find("DOUB") !=-1:
        aBN = 2
    elif aBS.find("TRIP") !=-1:
        aBN = 3
    elif aBS.find("AROM") !=-1:
        aBN =1.5
    else:
        aBN = 1   # e.g. metal coordination bonds

    return aBN


def checkRepAtomTypes(tFileName, tS):

    aRet = True
    aFile = open(tFileName, "r")
    allLs = aFile.readlines()
    aFile.close()
   
    nAllTP = 0
    allAtmTP = {}
    lAtmSec = False
    for aL in allLs:
        if aL.find("H FormType:") !=-1:
            lAtmSec = False
        elif aL.find("ATOMS:") != -1:
            lAtmSec = True
        elif lAtmSec:
            strgrp = aL.strip().split()
            if len(strgrp)==4:
                if strgrp[1].find("H")==-1:
                    nAllTP +=1
                    aTP = strgrp[3].strip()
                    if not aTP in allAtmTP.keys():
                        allAtmTP[aTP] = 1
                    else:
                        allAtmTP[aTP] += 1

    
    nMax = 0
    maxTP = ""
    if nAllTP > 0:
        for aTP in allAtmTP.keys():
            if allAtmTP[aTP] > nMax:
                nMax = allAtmTP[aTP]
                maxTP = aTP
        r = float(nMax)/nAllTP

        if r >= tS :
            aRet = False    
        print ("number of All non-H atoms %d "%nAllTP)
        print ("Largest repeated type %s with %d atoms"%(maxTP, nMax))
        print("RatioMax is  %f "%r) 
        print("All Atom types : ")
        for aK in allAtmTP.keys():
            print("%s    %d"%(aK, allAtmTP[aK]))

    else:
        aRet = False

    return aRet

# Transfer datas-tructures such as atoms and bonds from the rdkit style to dictionary style (gemmi style)

def DataStrTransferAtomsAndBonds(tMol, tDAtoms, tDBonds, tMonoRoot):
    
    allAtoms = tMol.GetAtoms()
    conf = tMol.GetConformer()
    
    aIdSet = {}
    for aAt in allAtoms:
        aDAt = {}
        aDAt["_chem_comp_atom.comp_id"] = tMonoRoot
        aDAt["_chem_comp_atom.type_symbol"] = aAt.GetSymbol().upper()
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
        aDAt["_chem_comp_atom.atom_conn"] = []
        pos = conf.GetAtomPosition(aAt.GetIdx())
        #print("x=%5.4f y=%5.4f z=%5.4f" % (pos.x, pos.y, pos.z))  
        aDAt["_chem_comp_atom.x"] = pos.x
        aDAt["_chem_comp_atom.y"] = pos.y
        aDAt["_chem_comp_atom.z"] = pos.z
        aDAt["_chem_comp_atom.model_Cartn_x"] = pos.x
        aDAt["_chem_comp_atom.model_Cartn_y"] = pos.y
        aDAt["_chem_comp_atom.model_Cartn_z"] = pos.z
        aDAt["_chem_comp_atom.pdbx_model_Cartn_x_ideal"] = pos.x
        aDAt["_chem_comp_atom.pdbx_model_Cartn_y_ideal"] = pos.y
        aDAt["_chem_comp_atom.pdbx_model_Cartn_z_ideal"] = pos.z
        
        tDAtoms.append(aDAt)
    
    for aBo in tMol.GetBonds():
        aDB = {}
        aDB["_chem_comp_bond.bond_serial_number"] = aBo.GetIdx()
        idx1 = aBo.GetBeginAtomIdx()
        id1  = tDAtoms[idx1]["_chem_comp_atom.atom_id"]
        idx2 = aBo.GetEndAtomIdx()
        id2  = tDAtoms[idx2]["_chem_comp_atom.atom_id"]
        aDB["_chem_comp_bond.atom_serial_number_1"] = idx1
        aDB["_chem_comp_bond.atom_serial_number_2"] = idx2
        aDB["_chem_comp_bond.atom_id_1"] = id1
        aDB["_chem_comp_bond.atom_id_2"] = id2
        
        aDB["_chem_comp_bond.value_order"] =str(aBo.GetBondType()) 
        tDAtoms[idx1]["_chem_comp_atom.atom_conn"].append(idx2)
        tDAtoms[idx2]["_chem_comp_atom.atom_conn"].append(idx1)
        tDBonds.append(aDB)
        
     
# Geometries releated

def getLengBetween2Pos(tPos1, tPos2):
    
    aVec = np.array([tPos1[0] - tPos2[0], tPos1[1] - tPos2[1], tPos1[2] - tPos2[2]])
    aLeng = np.linalg.norm(aVec)
    
    return aLeng

def getPlaneNormal(tPos1, tPos2, tPos3):
    # Return a unit normal vector of the plane by 3-points
    pts = np.array([[tPos1[0], tPos1[1], tPos1[2]],
                   [tPos2[0], tPos2[1], tPos2[2]],
                   [tPos3[0], tPos3[1], tPos3[2]]])
    
    vec1 = pts[0, :] - pts[1, :]
    v1l  = np.linalg.norm(vec1)
    vec1   = vec1/v1l
    vec2 = pts[0, :] - pts[2, :]
    v2l  = np.linalg.norm(vec2)
    vec2   = vec2/v2l
    vec3 = np.cross(vec1, vec2)
    vec3 = vec3/np.linalg.norm(vec3)
    return vec3

def getMassCenterPos(tSetPos):
    
    
    xSum=0.0
    ySum=0.0
    zSum=0.0
    
    for aP in tSetPos:
        xSum += aP[0]
        ySum += aP[1]
        zSum += aP[2]
        
    numPs = len(tSetPos)
    
    return [xSum/numPs, ySum/numPs, zSum/numPs]
    
def getAUnitDirVec2P(tPos1, tPos2):
    
    aDV = np.array([tPos1[0] - tPos2[0], tPos1[1] - tPos2[1], tPos1[2] - tPos2[2]])
    nL  = np.linalg.norm(aDV)
    
    return aDV/nL

def getAUnitDirFromCrossProd(tVec1, tVec2):
    
    aV1 = np.array(tVec1)
    aV2 = np.array(tVec1)

    aCP = np.linalg.cross(aV1, aV2)
    aCPL = np.linalg.norm(aCP)
    aRet = None 
    if aCPL > 0:
        aRet = aCP/aCPL
        
    return aRet
    
    

def getAUnitDirVec(tPos1, tSetPos):
    
    Pos2 = getMassCenterPos(tSetPos)
    
    return getAUnitDirVec2P(tPos1, Pos2)


def getAAngFrom3Ps(tPosC, tPos1, tPos2):
    
    a = np.array([tPos1[0],  tPos1[1], tPos1[2]])
    b = np.array([tPosC[0], tPosC[1], tPosC[2]])
    c = np.array([tPos2[0],  tPos2[1], tPos2[2]])

    ba = a - b
    l_ba = np.linalg.norm(ba)
    bc = c - b
    l_bc = np.linalg.norm(bc)
    
    cosine_angle = np.dot(ba/l_ba, bc/l_bc)
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)

def setOneAtomCoordsOutCB(tIdxCB, tIdxNonCB, tAtomsCB, tAtomsNCB, tLeng):
    print("here ", tAtomsCB[tIdxCB]["_chem_comp_atom.atom_id"])
    print(tAtomsCB[tIdxCB]["_chem_comp_atom.atom_conn"])
    if len(tAtomsCB[tIdxCB]["_chem_comp_atom.atom_conn"]) > 3:
        aSetPos = []
        for aIdx in tAtomsCB[tIdxCB]["_chem_comp_atom.atom_conn"]:
            x = float(tAtomsCB[aIdx]["_chem_comp_atom.model_Cartn_x"])
            y = float(tAtomsCB[aIdx]["_chem_comp_atom.model_Cartn_y"])
            z = float(tAtomsCB[aIdx]["_chem_comp_atom.model_Cartn_z"])
            aPos = [x,y,z]
            aSetPos.append(aPos)
            
        x1 = float(tAtomsCB[tIdxCB]["_chem_comp_atom.model_Cartn_x"])
        y1 = float(tAtomsCB[tIdxCB]["_chem_comp_atom.model_Cartn_y"])
        z1 = float(tAtomsCB[tIdxCB]["_chem_comp_atom.model_Cartn_z"])
        aPosNonH = [x1, y1, z1]
        aDirV= getAUnitDirVec(aPosNonH, aSetPos)
        
        oldX = float(tAtomsNCB[tIdxNonCB]["_chem_comp_atom.x"]) 
        oldY = float(tAtomsNCB[tIdxNonCB]["_chem_comp_atom.y"]) 
        oldZ = float(tAtomsNCB[tIdxNonCB]["_chem_comp_atom.z"])
        
        
        newX = x1 + tLeng*aDirV[0]
        newY = y1 + tLeng*aDirV[1]
        newZ = z1 + tLeng*aDirV[2]
        
        
        replaceAtomCoords(tAtomsNCB[tIdxNonCB], newX, newY, newZ)
        
        #print("new coords for ", tAtomsNCB[tIdxNonCB]["_chem_comp_atom.atom_id"])
        #print("X: ", newX)
        #print("Y: ", newY)
        #print("Z: ", newZ)
        #print("conn atom ", tAtomsCB[tIdxCB]["_chem_comp_atom.atom_id"])
        #print("Its x ", x1)
        #print("Its y ", y1)
        #print("Its z ", z1)
        #print("The leng between ",  getLengBetween2Pos([x1, y1,z1], [newX, newY, newZ]))
        
        #tAtomsNCB[tIdxNonCB]["_chem_comp_atom.model_Cartn_x"] = newX
        #tAtomsNCB[tIdxNonCB]["_chem_comp_atom.model_Cartn_y"] = newY
        #tAtomsNCB[tIdxNonCB]["_chem_comp_atom.model_Cartn_z"] = newZ
        
        deltaP = [newX-oldX, newY-oldY, newZ-oldZ]
        
        
        return deltaP
    
def replaceAtomCoords(tAtom, tX, tY, tZ):
    
    # changes in all labels
    #print("tX=", tX)
    
    tAtom["_chem_comp_atom.model_Cartn_x"] = tX
    tAtom["_chem_comp_atom.model_Cartn_y"] = tY
    tAtom["_chem_comp_atom.model_Cartn_z"] = tZ
    
    tAtom["_chem_comp_atom.x"] = tX
    tAtom["_chem_comp_atom.y"] = tY
    tAtom["_chem_comp_atom.z"] = tZ

    tAtom["_chem_comp_atom.pdbx_model_Cartn_x_ideal"] = tX
    tAtom["_chem_comp_atom.pdbx_model_Cartn_y_ideal"] = tY
    tAtom["_chem_comp_atom.pdbx_model_Cartn_z_ideal"] = tZ
    #print(tAtom)            

# Equiv-class related algorithms

def getLinkedGroups(tAtoms, tLinks, tMols):
    
    tGroups    = {}
    nAts       = len(tAtoms)
    tGroups[0] = 0
    for i in  range(1, nAts):
        id_i = tAtoms[i]['_chem_comp_atom.atom_id']
        tGroups[i] = i
        for j in range(i):
            id_j = tAtoms[j]['_chem_comp_atom.atom_id']
            tGroups[j] = tGroups[tGroups[j]]
            if id_j in tLinks[id_i]:
                tGroups[tGroups[tGroups[j]]] = i
    
    for i in range(nAts):
        tGroups[i] = tGroups[tGroups[i]];
    
    tMols = {}
    
    for iG in range(len(tGroups)):
        idxM = tGroups[iG]
        if not idxM in tMols:
            tMols[idxM] = []
        tMols[idxM].append(iG)

    print("number of fragments is ", len(tMols))
    aIdxF = 1
    for aIdxM in tMols.keys():
        print("Fragment ", aIdxF, " contains the following atoms: ")
        aIdxF+=1
        for aIdxA in tMols[aIdxM]:
            print("Atom : %s "%tAtoms[aIdxA]['_chem_comp_atom.atom_id'])
            
def getLinkedGroups2(tAtoms, tLinks):
    
    tGroups    = {}
    nAts       = len(tAtoms)
    tGroups[0] = 0
    for i in  range(1, nAts):
        id_i = tAtoms[i]['_chem_comp_atom.atom_id']
        tGroups[i] = i
        for j in range(i):
            id_j = tAtoms[j]['_chem_comp_atom.atom_id']
            tGroups[j] = tGroups[tGroups[j]]
            if id_j in tLinks[id_i]:
                tGroups[tGroups[tGroups[j]]] = i
    
    for i in range(nAts):
        tGroups[i] = tGroups[tGroups[i]];
    
    tMols = {}
    
    for iG in range(len(tGroups)):
        idxM = tGroups[iG]
        if not idxM in tMols:
            tMols[idxM] = []
        tMols[idxM].append(iG)

    print("number of fragments is ", len(tMols))
    aIdxF = 1
    for aIdxM in tMols.keys():
        print("Fragment ", aIdxF, " contains the following atoms: ")
        aIdxF+=1
        for aIdxA in tMols[aIdxM]:
            print("Atom : %s "%tAtoms[aIdxA]['_chem_comp_atom.atom_id']) 
    return tMols    

# Graph-related algorithms 

def aMolToAGraph(tAtoms, tBonds, tFileId):
    
    weightMap = {}
    # Using periodic tables 
    weightMap["B"] = 2
    weightMap["C"] = 3
    weightMap["N"] = 4
    weightMap["O"] = 5
    weightMap["S"] = 6
    weightMap["SE"]= 7
    weightMap["P"] = 8
    weightMap["AS"] = 9
    weightMap["SI"] = 9
    weightMap["GA"] = 9
    weightMap["GI"] = 9
    weightMap["IN"] = 9
    weightMap["OTHER"] =1
    
    
    
    aG = nx.Graph(name=tFileId)
    # Add atoms info to a node
    for aIdx, aAt in enumerate(tAtoms):
        
        aG.add_node(aIdx)
        
        for aK in aAt.keys():
            aG.nodes[aIdx][aK] = aAt[aK]
            elem = aAt["_chem_comp_atom.type_symbol"].upper()
            #if not elem in weightMap:
            #    aG.nodes[aIdx]["_chem_comp_atom.type_symbol"] = "FE"
            
    #print("=======================================")
    #print("The node properties of the graph are : ")
    #print(aG.nodes.data()) 
    
    for aB in tBonds:
        idxA1 = int(aB["_chem_comp_bond.atom_serial_number_1"])
        idxA2 = int(aB["_chem_comp_bond.atom_serial_number_2"])
        elem1 = tAtoms[idxA1]["_chem_comp_atom.type_symbol"].upper()
        elem2 = tAtoms[idxA2]["_chem_comp_atom.type_symbol"].upper()
        aWeight = 0 
        if elem1 in weightMap:
            aWeight += weightMap[elem1]
        else:
            aWeight +=1
        if elem2 in weightMap:
            aWeight +=  weightMap[elem2]
        else:
            aWeight +=1    
        
        aG.add_edge(idxA1, idxA2, weight=aWeight)
        for aK in aB.keys():
            aG[idxA1][idxA2][aK] = aB[aK]
    
    #print("=======================================")
    #print("The edge properties of the graph are : ")
    #print(aG.edges.data())

    return aG

def MatchTwoGraphs(tG1, tG2):
    
    aGM = nx.algorithms.isomorphism.GraphMatcher(tG1, tG2, node_match= lambda n1,n2:n1['_chem_comp_atom.type_symbol']==n2['_chem_comp_atom.type_symbol'], edge_match= lambda e1,e2: e1['weight'] == e2['weight'])
    
    return aGM.is_isomorphic(), aGM.mapping.items()


def cutRings(tMol, tBrokenIds):    
    
    origBonds = []
    for aB in tMol["bonds"]:
        origBonds.append(aB)
    tmpBonds = []
    extRings =[]
    selectedAtms    = []
    selectedRingMap = {}
    ringAtomMap     = {}
    numRingMap      = {}
    for aRIdx in tMol["rings"]:
        for aRA in tMol["rings"][aRIdx]:
            aId = aRA["_chem_comp_ring_atom.atom_id"]
            if not aId in numRingMap:
                numRingMap[aId] = []
            numRingMap[aId].append(aRIdx)
    for aRIdx in tMol["rings"]:
        #print(aRIdx)
        for aRA in tMol["rings"][aRIdx]:
            aId = aRA["_chem_comp_ring_atom.atom_id"]
            if not aRIdx in extRings and len(numRingMap[aId])==1:
                extRings.append(aRIdx)
                selectedAtms.append(aRA["_chem_comp_ring_atom.atom_id"])
                print("%s is selected"%aRA["_chem_comp_ring_atom.atom_id"])
                break                                    
        if not aRIdx in ringAtomMap:
            ringAtomMap[aRIdx] = []
        for aRA in tMol["rings"][aRIdx]:
            aId = aRA["_chem_comp_ring_atom.atom_id"]
            selectedRingMap[aId] = aRIdx
            ringAtomMap[aRIdx].append(aId)
            
            
    #print(selectedAtms)
    #print(selectedRingMap)
    #print(numRingMap)
    
    extRings =[]
    for aB in tMol["bonds"]:
        id1 = aB["_chem_comp_bond.atom_id_1"]
        id2 = aB["_chem_comp_bond.atom_id_2"]
        #print("Check bond between %s and %s "%(id1, id2))
        lInc = False
        if id1 in selectedAtms and id1 in selectedRingMap and id1 in numRingMap:
            if not selectedRingMap[id1] in extRings and not len(numRingMap[id1]) > 1:
                if id2 in ringAtomMap[selectedRingMap[id1]]:
                    lInc = True
                    extRings.append(selectedRingMap[id1])
        elif id2 in selectedAtms and id2 in selectedRingMap and id2 in numRingMap:
            if not selectedRingMap[id2] in extRings and not len(numRingMap[id2]) > 1:
                if id1 in ringAtomMap[selectedRingMap[id2]]:
                    lInc = True
                    extRings.append(selectedRingMap[id2])
         
        if not lInc:
            tmpBonds.append(aB)
        else:
            print("Bond between %s and %s is broken"%(id1, id2))
            if not id1 in tBrokenIds:
                tBrokenIds.append(id1)
            if not id2 in tBrokenIds:
                tBrokenIds.append(id2)
            
    tMol["bonds"] = [] 
    for aB in tmpBonds:
        tMol["bonds"].append(aB)
    

def setTreeNodeOrder(tConn, tStartId, tTree, tStartList):
    
    if not tStartId in tStartList: 
        tStartList.append(tStartId)
    #print("tStartList =", tStartList)
    if tStartId in tConn :
        #print("Start id ", tStartId)
        #print("Its children ")
        #print(tConn[tStartId])
    
        for aId in tConn[tStartId] :
            #print(aId, " is a ch of ", tStartId)
            if not aId in tStartList:
                #print(aId, " is assigned ")
                tTree[aId]["atom_back"] = tStartId
                if not "atom_forward" in tTree[aId]:
                    tTree[aId]["atom_forward"] =[]
                #print("its connection ", tConn[aId] )
                if len(tConn[aId]) > 1 :
                    for aNId in tConn[aId]:
                        if aNId != tStartId:
                            #print("Next Lev ", aNId)
                            tTree[aId]["atom_forward"].append(aNId)
                            setTreeNodeOrder(tConn, aId, tTree, tStartList)
                else:
                    pass
                    #print("No start")
                    #print("Back to ", aId, " next" )
       
        return             
          

        