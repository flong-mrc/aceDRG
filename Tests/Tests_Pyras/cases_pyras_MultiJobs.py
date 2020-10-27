#!/usr/bin/env python3

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

inRoot = os.getcwd()
indir1 = os.path.join(inRoot, "inPyraCCD", "*.cif")
indir2 = os.path.join(inRoot, "inPyraCCP4", "*.cif")

for aCif in glob.glob(indir1):
    ligName = os.path.basename(aCif).strip().split(".")[0]
    outRoot = "Test_" + ligName + "_fromPyraCCD_usingCoords"    
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s -p > %s"\
              %(aCif, outRoot, logName)
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
    outRoot = "Test_" + ligName + "_fromPyraCCD"    
    logName = outRoot + ".log"
    cmdLine = "acedrg  -c %s -o %s > %s"\
              %(aCif, outRoot, logName)
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

for aCif in glob.glob(indir2):
    ligName = os.path.basename(aCif).strip().split(".")[0]
    outRoot = "Test_" + ligName + "_fromPyraCCP4_usingCoords"    
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s -p > %s"\
              %(aCif, outRoot, logName)
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

    outRoot = "Test_" + ligName + "_fromPyraCCP4"    
    logName = outRoot + ".log"
    cmdLine = "acedrg -c %s -o %s > %s"\
              %(aCif, outRoot,  logName)
    print(cmdLine)
    numAllJobs += 1
    os.system(cmdLine)
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


