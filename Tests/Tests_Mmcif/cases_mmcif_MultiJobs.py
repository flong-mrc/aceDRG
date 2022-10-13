#!/usr/bin/env ccp4-python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
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

numAllJobs =0
numJobsS   =0
numJobsF   =0

doneNames = []
failNames = []

lib1 = "/Users/flong/DB/PDB_related/PDB_Ligands/Cif"
lib2 = "/Users/flong/CCP4/ccp4-8.0/lib/data/monomers"
inRoot = os.getcwd()

inListFN = os.path.join(inRoot, "inMmcifCCD", "in.list")
inListF  = open(inListFN, "r")
a        = inListF.readlines()
inListF.close()
inList1 = []
for aL in a:
    strGrp=aL.strip().split()
    for aE in strGrp:
        inList1.append(aE.strip())

inListFN = os.path.join(inRoot, "inMmcif", "in.list")
inListF  = open(inListFN, "r")
b        = inListF.readlines()
inListF.close()
inList2 = []
for aL in b:
    strGrp=aL.strip().split()
    for aE in strGrp:
        inList2.append(aE.strip())


for aCif in inList1:
    ligName = aCif.strip().split(".")[0]
    subDir = ligName[0].lower()
    outRoot = "Test_" + ligName + "_fromMmcifCCD_usingCoords"    
    logName = outRoot + ".log"
    inF     = os.path.join(lib1, subDir, aCif)
    cmdLine = "acedrg -c %s -o %s -p > %s"\
              %(inF, outRoot, logName)
    print(cmdLine)
    numAllJobs += 1
    lRun=os.system(cmdLine)
    if lRun :
        print("%s runtime error "%outRoot)
    outCif = "%s.cif"%outRoot
    if os.path.isfile(outCif):
        numJobsS   +=1
        doneNames.append(outRoot)
    else:
        numJobsF   +=1
        failNames.append(outRoot)
    outRoot = "Test_" + ligName + "_fromMmcifCCD"    
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s > %s"\
              %(inF, outRoot, logName)
    print(cmdLine)
    numAllJobs += 1
    lRun=os.system(cmdLine)
    if lRun :
        print("%s runtime error "%outRoot)
    outCif = "%s.cif"%outRoot
    if os.path.isfile(outCif):
        numJobsS   +=1
        doneNames.append(outRoot)
    else:
        numJobsF   +=1
        failNames.append(outRoot)

for aCif in inList2:
    ligName = aCif.strip().split(".")[0]
    subDir = ligName[0].lower()
    outRoot = "Test_" + ligName + "_fromMmcifCCP4_usingCoords"    
    inF     = os.path.join(lib2, subDir, aCif)
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s -p > %s"\
              %(inF, outRoot, logName)
    print(cmdLine)
    numAllJobs += 1
    lRun=os.system(cmdLine)
    if lRun :
        print("%s runtime error "%outRoot)
    outCif = "%s.cif"%outRoot
    if os.path.isfile(outCif):
        numJobsS   +=1
        doneNames.append(outRoot)
    else:
        numJobsF   +=1
        failNames.append(outRoot)

    outRoot = "Test_" + ligName + "_fromMmcifCCP4"    
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s > %s"\
              %(inF, outRoot,  logName)
    print(cmdLine)
    numAllJobs += 1
    lRun = os.system(cmdLine)
    if lRun :
        print("%s runtime error "%outRoot)
    outCif = "%s.cif"%outRoot
    if os.path.isfile(outCif):
        numJobsS   +=1
        doneNames.append(outRoot)
    else:
        numJobsF   +=1
        failNames.append(outRoot)

print("Total Number of jobs running", numAllJobs)
print("Total Number of job successfully finished", numJobsS)
if numJobsF > 0:
    print("Total Number of job failed %d. They are: "%numJobsF)
    for aName in failNames:
        print(aName)


