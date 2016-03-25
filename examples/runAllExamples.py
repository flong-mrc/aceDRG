import os, os.path, sys, glob, string

numAllJobs =0
numJobsS   =0
numJobsF   =0

if len(sys.argv) > 1:
    print "Usuage: python runAllExamples.py"
    sys.exit()

doneNames = []
failNames = []

print "========================================================="
print "|  run all examples with input SMILES files at ./inSmil |"
print "========================================================="
for aF in glob.glob("./inSmi/*.smiles"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_fromSmiles"%aFRoot
        cmdL = "acedrg -i %s  -o  %s "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)
            

print "========================================================="
print "|  run all examples with input MOL files at ./inMOL     |"
print "========================================================="
for aF in glob.glob("./inMol/*.mol"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_fromMol"%aFRoot
        cmdL = "acedrg -m %s  -o  %s_fromMol -p "%(aF, aFRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "========================================================="
print "|  run all examples with input SYBYL/MOL2 files at ./inMOL2     |"
print "========================================================="
for aF in glob.glob("./inMol2/*.mol2"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_SYBYL_Mol2"%aFRoot
        cmdL = "acedrg -g %s  -o  %s -p "%(aF, rRoot) 
        numAllJobs += 1
        print cmdL
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "============================================================"
print "|  run all examples with input mmCif files at ./inMmcifPDB |"
print "============================================================"
for aF in glob.glob("./inMmcifPDB/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_fromMmcifPDB"%aFRoot
        cmdL = "acedrg -c %s  -o  %s -p "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "========================================================="
print "|  run all examples with input mmCif files at ./inMmcif |"
print "========================================================="
for aF in glob.glob("./inMmcif/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_fromMmcif"%aFRoot
        cmdL = "acedrg -c %s  -o  %s -p "%(aF, rRoot) 
        numAllJobs += 1
        print cmdL
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "========================================================="
print "|  run all examples containing pyranose at ./inPyraCCD  |"
print "========================================================="
for aF in glob.glob("./inPyraCCD/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_inPyraCCD"%aFRoot
        cmdL = "acedrg -c %s  -o  %s -p "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)

print "========================================================="
print "|  run all examples containing pyranose at ./inPyraCCP4 |"
print "========================================================="
for aF in glob.glob("./inPyraCCP4/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_inPyraCCP4"%aFRoot
        cmdL = "acedrg -c %s  -o  %s -p "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "========================================================="
print "|  run all examples in inFuns                           |"
print "========================================================="
for aF in glob.glob("./inFuns/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_Funs"%aFRoot
        cmdL = "acedrg -c %s -o  %s -p "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)


print "=============================================================="
print "|  run those cif in CCP4 monomer lib which had some problems |"
print "=============================================================="
for aF in glob.glob("./inDictProb/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_inDictProb"%aFRoot
        cmdL = "acedrg -c %s  -o  %s -p "%(aF, rRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)

print "=============================================================="
print "|  run all examples with input mmCif files at ./inProbsCases |"
print "=============================================================="
for aF in glob.glob("./inProbsCases/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        rRoot = "%s_inProbsCases_P"%aFRoot
        cmdL = "acedrg -c %s  -o  %s_inProbsCases_P -p "%(aF, aFRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)
         
        """
        rRoot = "%s_inProbsCases_T"%aFRoot
        cmdL = "acedrg -c %s  -o  %s_inProbsCases_T  "%(aF, aFRoot) 
        print cmdL
        numAllJobs += 1
        lRun=os.system(cmdL)
        if lRun :
            print "%s runtime error "%rRoot
        outCif = "%s.cif"%rRoot
        if os.path.isfile(outCif):
            numJobsS   +=1
            doneNames.append(rRoot)
        else:
            numJobsF   +=1
            failNames.append(rRoot)
        """
"""
print "========================================================="
print "|  run all examples with input small molecule Cif files |"
print "|  at ./inStdcif to generate molecules and derive       |"
print "|  unique bond lengths and angles, and atom types       |"
print "========================================================="
for aF in glob.glob("./inStdcif/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) !=0:
        cmdL = "acedrg -e -b %s  -r  %s -o  %s"%(aF, aFRoot, aFRoot) 
        print cmdL
        os.system(cmdL)

print "========================================================="
print "|  run all small molecule Cif files at ./inStdcif and   |"
print "|  and generate a table containing unique bond lengths  |"
print "|  and angles, and atom types in those cif files        |"
print "========================================================="
if os.path.isdir("inStdcif"):
    cmdL = "acedrg -e -d %s -o %s " %("inStdcif", "Summary_inStdcif") 
    print cmdL
    os.system(cmdL)
"""

print "Total Number of jobs running for dictionary generation ", numAllJobs 
print "Total Number of job successfully finished", numJobsS 
print "Total Number of job failed %d. They are: "%numJobsF
for aName in failNames:
    print aName

