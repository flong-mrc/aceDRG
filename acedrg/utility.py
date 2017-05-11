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
import select
import random

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

    for aK in tDict.keys():
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
