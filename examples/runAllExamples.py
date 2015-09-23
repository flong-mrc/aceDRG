import os, os.path, sys, glob, string

if len(sys.argv) > 1:
    print "Usuage: python runAllExamples.py"
    sys.exit()

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
        cmdL = "acedrg -i %s  -o  %s_fromSmiles "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

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
        cmdL = "acedrg -m %s  -o  %s_fromMol "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)


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
        cmdL = "acedrg -g %s  -o  %s_from_SYBYL_Mol2 "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

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
        cmdL = "acedrg -c %s  -o  %s_fromMmcifPDB "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

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
        cmdL = "acedrg -c %s  -o  %s_fromMmcif "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

print "========================================================="
print "|  run all examples containing pyranose at ./inPyro     |"
print "========================================================="
for aF in glob.glob("./inPyro/*.cif"):
    aFRoot = os.path.basename(aF).strip().split(".")[0].strip()
    if len(aFRoot) >= 3:
        MonoName = aFRoot[-3:]
    else:
        MonoName = aFRoot

    if len(aFRoot) !=0:
        cmdL = "acedrg -c %s  -o  %s_fromMmcif -k 10 "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

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
        cmdL = "acedrg -c %s -o  %s_Funs "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)

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
        cmdL = "acedrg -c %s  -o  %s_fromMmcif "%(aF, aFRoot) 
        print cmdL
        os.system(cmdL)
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

