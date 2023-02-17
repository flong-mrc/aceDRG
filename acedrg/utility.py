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

# Functions for numbers

def isInt(t_string):

    aRet = True

    t_string = t_string.strip()
    try:
        aInt = int(t_string)
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
