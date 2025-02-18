cmake_minimum_required(VERSION 3.12)
project(acedrg C CXX)

find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter)
if(WIN32)
        string(REPLACE "\\" "/" Python_SITELIB "${Python_SITELIB}")
endif(WIN32)

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
if (MSVC)
set(SRC ${SRC} src/kernel/win32/getopt.c)
endif (MSVC)

add_executable(libmol ${SRC})

install(TARGETS libmol      DESTINATION libexec)
install(DIRECTORY ../acedrg-tables-git/tables    DESTINATION share/acedrg)
install(DIRECTORY acedrg    DESTINATION ${Python_SITELIB})
install(DIRECTORY utils     DESTINATION ${Python_SITELIB}/acedrg)
install(PROGRAMS bin/acedrg DESTINATION bin)
if (WIN32)
    install(PROGRAMS bin/acedrg.bat DESTINATION bin)
endif()
install(DIRECTORY ACEDRG_ccp4i/scripts   DESTINATION share/ccp4i)
install(DIRECTORY ACEDRG_ccp4i/tasks     DESTINATION share/ccp4i)
install(DIRECTORY ACEDRG_ccp4i/templates DESTINATION share/ccp4i)

