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


