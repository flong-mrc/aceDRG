#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=/sw/bin/g++-4 -O3
CXX=/sw/bin/g++-4 -O3
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Debug
CND_DISTDIR=bin
CND_BUILDDIR=build

# Include 

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/kernel/angle.o \
	${OBJECTDIR}/src/kernel/chain.o \
	${OBJECTDIR}/src/kernel/coordPolyhedra.o \
	${OBJECTDIR}/src/forcefield/getDerivs.o \
	${OBJECTDIR}/src/kernel/plane.o \
	${OBJECTDIR}/src/kernel/utility.o \
	${OBJECTDIR}/src/kernel/molecule.o \
	${OBJECTDIR}/src/kernel/restraintList.o \
	${OBJECTDIR}/src/kernel/MolGenerator.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/kernel/libglink.o \
	${OBJECTDIR}/src/kernel/PDBFile.o \
	${OBJECTDIR}/src/kernel/CCP4AtomType.o \
	${OBJECTDIR}/src/kernel/codBondAndAngleGroups.o \
	${OBJECTDIR}/src/kernel/codClassify.o \
	${OBJECTDIR}/src/kernel/atomAssembly.o \
	${OBJECTDIR}/src/kernel/chiral.o \
	${OBJECTDIR}/src/kernel/file.o \
	${OBJECTDIR}/src/kernel/getExtraRestrs.o \
	${OBJECTDIR}/src/kernel/torsion.o \
	${OBJECTDIR}/src/kernel/periodicTable.o \
	${OBJECTDIR}/src/kernel/SmileTool.o \
	${OBJECTDIR}/src/kernel/inputParams.o \
	${OBJECTDIR}/src/kernel/DictCifFILE.o \
	${OBJECTDIR}/src/kernel/atom.o \
	${OBJECTDIR}/src/go/GlobOpt.o \
	${OBJECTDIR}/src/kernel/bond.o \
	${OBJECTDIR}/src/kernel/DnaRna.o \
	${OBJECTDIR}/src/kernel/libgmodel.o \
	${OBJECTDIR}/src/kernel/TransCoord.o \
	${OBJECTDIR}/src/kernel/AllSystem.o \
	${OBJECTDIR}/src/kernel/secondaryStructures.o \
	${OBJECTDIR}/src/kernel/ssbond.o \
	${OBJECTDIR}/src/kernel/crystInfo.o \
	${OBJECTDIR}/src/kernel/residue.o \
	${OBJECTDIR}/src/kernel/checkEnvAndGetMode.o \
	${OBJECTDIR}/src/kernel/atomsTree.o \
	${OBJECTDIR}/src/kernel/neighbList.o \
	${OBJECTDIR}/src/kernel/ring.o \
	${OBJECTDIR}/src/forcefield/getObjValue.o \
	${OBJECTDIR}/src/kernel/MolSdfFile.o \
	${OBJECTDIR}/src/go/LinAlg.o \
	${OBJECTDIR}/src/kernel/ExtraRestrDictFile.o \
	${OBJECTDIR}/src/kernel/MonomerLib.o 


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-lsqlite3
CXXFLAGS=-lsqlite3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libmol

${CND_DISTDIR}/libmol: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}
	${LINK.cc} -o ${CND_DISTDIR}/libmol ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/src/kernel/angle.o: src/kernel/angle.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/angle.o src/kernel/angle.cpp

${OBJECTDIR}/src/kernel/chain.o: src/kernel/chain.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/chain.o src/kernel/chain.cpp

${OBJECTDIR}/src/kernel/coordPolyhedra.o: src/kernel/coordPolyhedra.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/coordPolyhedra.o src/kernel/coordPolyhedra.cpp

${OBJECTDIR}/src/forcefield/getDerivs.o: src/forcefield/getDerivs.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/forcefield
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/forcefield/getDerivs.o src/forcefield/getDerivs.cpp

${OBJECTDIR}/src/kernel/plane.o: src/kernel/plane.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/plane.o src/kernel/plane.cpp

${OBJECTDIR}/src/kernel/utility.o: src/kernel/utility.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/utility.o src/kernel/utility.cpp

${OBJECTDIR}/src/kernel/molecule.o: src/kernel/molecule.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/molecule.o src/kernel/molecule.cpp

${OBJECTDIR}/src/kernel/restraintList.o: src/kernel/restraintList.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/restraintList.o src/kernel/restraintList.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/kernel/libglink.o: src/kernel/libglink.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/libglink.o src/kernel/libglink.cpp

${OBJECTDIR}/src/kernel/PDBFile.o: src/kernel/PDBFile.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/PDBFile.o src/kernel/PDBFile.cpp

${OBJECTDIR}/src/kernel/CCP4AtomType.o: src/kernel/CCP4AtomType.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/CCP4AtomType.o src/kernel/CCP4AtomType.cpp

${OBJECTDIR}/src/kernel/codBondAndAngleGroups.o: src/kernel/codBondAndAngleGroups.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/codBondAndAngleGroups.o src/kernel/codBondAndAngleGroups.cpp

${OBJECTDIR}/src/kernel/codClassify.o: src/kernel/codClassify.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/codClassify.o src/kernel/codClassify.cpp

${OBJECTDIR}/src/kernel/atomAssembly.o: src/kernel/atomAssembly.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/atomAssembly.o src/kernel/atomAssembly.cpp

${OBJECTDIR}/src/kernel/file.o: src/kernel/file.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/file.o src/kernel/file.cpp

${OBJECTDIR}/src/kernel/getExtraRestrs.o: src/kernel/getExtraRestrs.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/getExtraRestrs.o src/kernel/getExtraRestrs.cpp

${OBJECTDIR}/src/kernel/chiral.o: src/kernel/chiral.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/chiral.o src/kernel/chiral.cpp

${OBJECTDIR}/src/kernel/torsion.o: src/kernel/torsion.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/torsion.o src/kernel/torsion.cpp

${OBJECTDIR}/src/kernel/periodicTable.o: src/kernel/periodicTable.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/periodicTable.o src/kernel/periodicTable.cpp

${OBJECTDIR}/src/kernel/SmileTool.o: src/kernel/SmileTool.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/SmileTool.o src/kernel/SmileTool.cpp

${OBJECTDIR}/src/kernel/inputParams.o: src/kernel/inputParams.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/inputParams.o src/kernel/inputParams.cpp

${OBJECTDIR}/src/kernel/DictCifFILE.o: src/kernel/DictCifFILE.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/DictCifFILE.o src/kernel/DictCifFILE.cpp

${OBJECTDIR}/src/kernel/atom.o: src/kernel/atom.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/atom.o src/kernel/atom.cpp

${OBJECTDIR}/src/go/GlobOpt.o: src/go/GlobOpt.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/go
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/go/GlobOpt.o src/go/GlobOpt.cpp

${OBJECTDIR}/src/kernel/DnaRna.o: src/kernel/DnaRna.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/DnaRna.o src/kernel/DnaRna.cpp

${OBJECTDIR}/src/kernel/bond.o: src/kernel/bond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/bond.o src/kernel/bond.cpp

${OBJECTDIR}/src/kernel/libgmodel.o: src/kernel/libgmodel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/libgmodel.o src/kernel/libgmodel.cpp

${OBJECTDIR}/src/kernel/TransCoord.o: src/kernel/TransCoord.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/TransCoord.o src/kernel/TransCoord.cpp

${OBJECTDIR}/src/kernel/AllSystem.o: src/kernel/AllSystem.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/AllSystem.o src/kernel/AllSystem.cpp

${OBJECTDIR}/src/kernel/secondaryStructures.o: src/kernel/secondaryStructures.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/secondaryStructures.o src/kernel/secondaryStructures.cpp

${OBJECTDIR}/src/kernel/ssbond.o: src/kernel/ssbond.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/ssbond.o src/kernel/ssbond.cpp

${OBJECTDIR}/src/kernel/residue.o: src/kernel/residue.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/residue.o src/kernel/residue.cpp

${OBJECTDIR}/src/kernel/checkEnvAndGetMode.o: src/kernel/checkEnvAndGetMode.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/checkEnvAndGetMode.o src/kernel/checkEnvAndGetMode.cpp

${OBJECTDIR}/src/kernel/atomsTree.o: src/kernel/atomsTree.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/atomsTree.o src/kernel/atomsTree.cpp

${OBJECTDIR}/src/kernel/neighbList.o: src/kernel/neighbList.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/neighbList.o src/kernel/neighbList.cpp

${OBJECTDIR}/src/kernel/ring.o: src/kernel/ring.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/ring.o src/kernel/ring.cpp

${OBJECTDIR}/src/forcefield/getObjValue.o: src/forcefield/getObjValue.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/forcefield
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/forcefield/getObjValue.o src/forcefield/getObjValue.cpp

${OBJECTDIR}/src/kernel/MolSdfFile.o: src/kernel/MolSdfFile.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/MolSdfFile.o src/kernel/MolSdfFile.cpp

${OBJECTDIR}/src/go/LinAlg.o: src/go/LinAlg.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/go
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/go/LinAlg.o src/go/LinAlg.cpp

${OBJECTDIR}/src/kernel/MonomerLib.o: src/kernel/MonomerLib.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/MonomerLib.o src/kernel/MonomerLib.cpp

${OBJECTDIR}/src/kernel/ExtraRestrDictFile.o: src/kernel/ExtraRestrDictFile.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -g -Wall -Iinclude/kernel -Iinclude/math -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/ExtraRestrDictFile.o src/kernel/ExtraRestrDictFile.cpp

${OBJECTDIR}/src/kernel/crystInfo.o: src/kernel/crystInfo.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -I/Applications/ccp4-6.3.0/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/crystInfo.o src/kernel/crystInfo.cpp

${OBJECTDIR}/src/kernel/MolGenerator.o: src/kernel/MolGenerator.cpp 
	${MKDIR} -p ${OBJECTDIR}/src/kernel
	${RM} $@.d
	$(COMPILE.cc) -Wall -Iinclude/kernel -Iinclude/math -Iinclude/go -Iinclude/forcefield -I/Applications/ccp4-6.3.0/include -MMD -MP -MF $@.d -o ${OBJECTDIR}/src/kernel/MolGenerator.o src/kernel/MolGenerator.cpp


# Clean Targets
clean:
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/libmol

