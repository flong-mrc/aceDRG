cmake_minimum_required(VERSION 2.8)
project(acedrg CXX)

include_directories(${CMAKE_SOURCE_DIR}/include/kernel
                    ${CMAKE_SOURCE_DIR}/include/math
                    ${CMAKE_SOURCE_DIR}/include/go
                    ${CMAKE_SOURCE_DIR}/include/forcefield)

set(SRC
    src/kernel/angle.cpp
    src/kernel/chain.cpp
    src/kernel/coordPolyhedra.cpp
    src/forcefield/getDerivs.cpp
    src/kernel/plane.cpp
    src/kernel/utility.cpp
    src/kernel/molecule.cpp
    src/kernel/restraintList.cpp
    src/kernel/MolGenerator.cpp
    main.cpp
    src/kernel/libglink.cpp
    src/kernel/PDBFile.cpp
    src/kernel/CCP4AtomType.cpp
    src/kernel/codBondAndAngleGroups.cpp
    src/kernel/codClassify.cpp
    src/kernel/atomAssembly.cpp
    src/kernel/chiral.cpp
    src/kernel/file.cpp
    src/kernel/getExtraRestrs.cpp
    src/kernel/torsion.cpp
    src/kernel/periodicTable.cpp
    src/kernel/SmileTool.cpp
    src/kernel/inputParams.cpp
    src/kernel/DictCifFILE.cpp
    src/kernel/atom.cpp
    src/go/GlobOpt.cpp
    src/kernel/bond.cpp
    src/kernel/DnaRna.cpp
    src/kernel/libgmodel.cpp
    src/kernel/TransCoord.cpp
    src/kernel/AllSystem.cpp
    src/kernel/secondaryStructures.cpp
    src/kernel/ssbond.cpp
    src/kernel/crystInfo.cpp
    src/kernel/residue.cpp
    src/kernel/checkEnvAndGetMode.cpp
    src/kernel/atomsTree.cpp
    src/kernel/neighbList.cpp
    src/kernel/ring.cpp
    src/forcefield/getObjValue.cpp
    src/kernel/MolSdfFile.cpp
    src/go/LinAlg.cpp
    src/kernel/ExtraRestrDictFile.cpp
    src/kernel/MonomerLib.cpp
    src/kernel/chemPropSet.cpp
)

add_executable(libmol ${SRC})

set(CCP4I_DIR ACEDRG_ccp4i)
set(CCP4I_DEST share/ccp4i)

install(TARGETS libmol DESTINATION libexec)
install(PROGRAMS bin/acedrg DESTINATION bin)
install(DIRECTORY tables DESTINATION share/acedrg/)

install(FILES ${CCP4I_DIR}/scripts/acedrg.script
        DESTINATION ${CCP4I_DEST}/scripts/)
install(FILES ${CCP4I_DIR}/templates/acedrg.com
        DESTINATION ${CCP4I_DEST}/templates/)
install(FILES ${CCP4I_DIR}/tasks/acedrg.def
              ${CCP4I_DIR}/tasks/acedrg.tcl
        DESTINATION ${CCP4I_DEST}/tasks/)

