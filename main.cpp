/* 
 * File:   main.cpp
 * Author: flong
 *
 * Created on January 24, 2012, 11:40 AM
 * last updated on Feb 12, 2012
 */

#ifndef CIFFILE_H
#include "include/kernel/DictCifFile.h"
#endif

#ifndef PDBFILE_H
#include "include/kernel/PDBFile.h"
#include "SmileTool.h"
#endif

#ifndef MOLSDFFILE_H
#include "include/kernel/MolSdfFile.h"
#endif

#ifndef EXTRARESTRDICTFILE_H
#include "include/kernel/ExtraRestrDictFile.h"
#endif

#ifndef MONOMERLIB_H
#include "include/kernel/MonomerLib.h"
#endif

#ifndef CODCLASSIFY_H
#include "include/kernel/codClassify.h"
#endif

#ifndef CODBONDANDANGLEGROUPS_H
#include "codBondAndAngleGroups.h"
#endif

#ifndef CHECKENVANDGETMODE_H
#include "CheckEnvAndGetMode.h"
#endif

#ifndef UTILITY_H
#include "include/kernel/utility.h"
#endif

#ifndef GLOBOPT_H
#include "include/go/GlobOpt.h"
#endif

#ifndef MOLGENERATOR_H
#include "MolGenerator.h"
#endif

#ifndef LINALG_H
#include "include/go/LinAlg.h"
#endif


//using namespace GO;
//using namespace FF;

/*
 *  A temporary version of main, some re-organization under way
 */
int main(int argc, char** argv) {
    
    
    LIBMOL::CheckEnvAndGetMode AJob(argc, argv);
    
    // std::cout << "The work mode is " << AJob.workMode << std::endl;
    
    /*
    for (std::map<ID,ID>::iterator iKW=checkMode.IOEntries.begin();
            iKW !=checkMode.IOEntries.end(); iKW++)
    {   
        std::cout << iKW->first << "\t" << iKW->second << std::endl;
    }
    */
    
    if (AJob.workMode == 11 
       || AJob.workMode ==12
       || AJob.workMode ==13
       || AJob.workMode ==14)
    {
        if (AJob.workMode == 11 || AJob.workMode ==12)
        {
            if (AJob.workMode ==12)
            {
                if (AJob.IOEntries.find("inSmiName") !=  AJob.IOEntries.end())
                {
                    if (LIBMOL::isFileExist(AJob.IOEntries["inSmiName"].c_str()))
                    {
                        LIBMOL::SmileTool aSmiConvertor;
                        aSmiConvertor.SmileToCifUsingLibcheck(AJob.IOEntries["inSmiName"].c_str(),
                                                              AJob.IOEntries["inCifName"]);
                        std::cout << "The initial Cif file generated base on input SMILE file is: " 
                                  << std::endl <<  AJob.IOEntries["inCifName"] << std::endl;
                    }
                    else
                    {
                        std::cout << "You should provide a correct smile file name and location "
                                  << std::endl;
                        exit(1);
                    }
                }
                else 
                {
                    std::cout << "Program Bug, this mode requires a SMILE file "
                              << std::endl;
                }
            }
            
            LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
            
            LIBMOL::AllSystem   aTargetSystem(dataFromCif); 
        
            if ( (int)aTargetSystem.allAtoms.size() > 0)
            {   
                aTargetSystem.setupAllTargetValuesFromCOD(AJob.IOEntries["userOutName"].c_str(), 
                                                          AJob.IOEntries["monoRootName"]);
               
                // aTargetSystem.chiralExch();
                if ( !aTargetSystem.containMetal()) 
                {
                    // coordinate generator works only for a system of "organic" atoms. Tempo!
                    GO::FindGlobMin  aGlobMinSystem(aTargetSystem, 
                                                    AJob.IOEntries["userOutName"].c_str(),
                                                    AJob.IOEntries["monoRootName"]);           
                    aGlobMinSystem.Driver();
            
                
                    LIBMOL::outMMCif(AJob.IOEntries["userOutName"].c_str(),
                                     AJob.IOEntries["monoRootName"], 
                                     aGlobMinSystem.allAtoms,
                                     aTargetSystem.allHAtomIdx,
                                     aGlobMinSystem.allBonds,
                                     aGlobMinSystem.allAngles,
                                     aGlobMinSystem.allTorsions,
                                     aGlobMinSystem.allRings,
                                     aGlobMinSystem.allPlanes,
                                     aGlobMinSystem.allChirals);
            
                    LIBMOL::outPDB(AJob.IOEntries["userOutName"].c_str(),
                                   AJob.IOEntries["monoRootName"], 
                                   aGlobMinSystem.allAtoms);
                }
                else
                {
                    LIBMOL::outMMCif(AJob.IOEntries["userOutName"].c_str(),
                                     AJob.IOEntries["monoRootName"], 
                                     aTargetSystem.allAtoms,
                                     aTargetSystem.allHAtomIdx,
                                     aTargetSystem.allBonds,
                                     aTargetSystem.allAngles,
                                     aTargetSystem.allTorsions,
                                     aTargetSystem.allRings,
                                     aTargetSystem.allPlanes,
                                     aTargetSystem.allChirals);
            
                    LIBMOL::outPDB(AJob.IOEntries["userOutName"].c_str(),
                                   AJob.IOEntries["monoRootName"], 
                                   aTargetSystem.allAtoms);
                    std::cout << "The system contain metal atoms." << std::endl
                              << "Restraint cif is the only output. "
                              << "No atom coordinates will be generated at the momemt. " << std::endl
                              << "They will be available soon. " << std::endl;
                }
            }
            else
            {
                std::cout << "The input Cif file " << AJob.IOEntries["inCifName"] 
                          << " contains NO atoms! check the file " << std::endl;
            }
        }
        else if (AJob.workMode == 13)
        {
            LIBMOL::MolSdfFile dataFromSdf(AJob.IOEntries["inSdfName"], std::ios::in);
            
            
            if (dataFromSdf.allMols.size() !=0)
            {
                int i =1;
                
                for (std::vector<LIBMOL::Molecule>::iterator iMol=dataFromSdf.allMols.begin();
                        iMol !=dataFromSdf.allMols.end(); iMol++)
                {
                    
                    LIBMOL::AllSystem   aTargetSystem(*iMol);
                    
                    if ( (int)aTargetSystem.allAtoms.size() > 0)
                    {
                        std::string tOutName(AJob.IOEntries["userOutName"]);
                        tOutName.append("_");
                        tOutName.append(LIBMOL::IntToStr(i));
                        aTargetSystem.setupAllTargetValuesFromCOD(tOutName.c_str(), 
                                                                  AJob.IOEntries["monoRootName"]);
                        // aTargetSystem.chiralExch();
                        if ( !aTargetSystem.containMetal()) 
                        {
                             
                            // coordinate generator works only for a system of "organic" atoms. Tempo!
                            GO::FindGlobMin  aGlobMinSystem(aTargetSystem, 
                                                            tOutName.c_str(),
                                                            AJob.IOEntries["monoRootName"]);           
                            aGlobMinSystem.Driver();
            
                            LIBMOL::outMMCif(tOutName.c_str(),
                                             AJob.IOEntries["monoRootName"], 
                                             aGlobMinSystem.allAtoms,
                                             aTargetSystem.allHAtomIdx,
                                             aGlobMinSystem.allBonds,
                                             aGlobMinSystem.allAngles,
                                             aGlobMinSystem.allTorsions,
                                             aGlobMinSystem.allRings,
                                             aGlobMinSystem.allPlanes,
                                             aGlobMinSystem.allChirals);
            
                            LIBMOL::outPDB(tOutName.c_str(),
                                           AJob.IOEntries["monoRootName"], 
                                           aGlobMinSystem.allAtoms);
                        }
                        else
                        {
                            LIBMOL::outMMCif(tOutName.c_str(),
                                             AJob.IOEntries["monoRootName"], 
                                             aTargetSystem.allAtoms,
                                             aTargetSystem.allHAtomIdx,
                                             aTargetSystem.allBonds,
                                             aTargetSystem.allAngles,
                                             aTargetSystem.allTorsions,
                                             aTargetSystem.allRings,
                                             aTargetSystem.allPlanes,
                                             aTargetSystem.allChirals);
                            
                            LIBMOL::outPDB(tOutName.c_str(),
                                           AJob.IOEntries["monoRootName"], 
                                           aTargetSystem.allAtoms);
                    
                            std::cout << "The system contain metal atoms." << std::endl
                                      << "Restraint cif is the only output. "
                                      << "No atom coordinates will be generated at the momemt. " << std::endl
                                      << "They will be available soon. " << std::endl;
                        }
                    }
                    i++;
                }
            }
            else
            {
                std::cout << "The input Sdf/Mol file " << AJob.IOEntries["inSdfName"] 
                          << " contains NO atoms! check the file " << std::endl;
            }
        }
    }
    else if (AJob.workMode == 2)
    {
        
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
        LIBMOL::CodClassify  aCodSystem(dataFromCif);
        aCodSystem.getAnglesFromPDB(AJob.IOEntries["inPdbName"]);
        // Temp for outRestraintCif2
        aCodSystem.outRestraintCif2(AJob.IOEntries["userOutName"].c_str(), AJob.IOEntries["monoRootName"]);
        aCodSystem.outPDB(AJob.IOEntries["userOutName"].c_str(), AJob.IOEntries["monoRootName"]);
       
    }
    else if (AJob.workMode==31 || AJob.workMode==32 || AJob.workMode==33)
    {
        std::cout << "WorkMode: Molecule generation" << std::endl;
        int aNBDepth = LIBMOL::StrToInt(AJob.IOEntries["NBDepth"]);
        if (AJob.workMode==31)
        {
            std::cout << "Input cif " << AJob.IOEntries["inCifNameB"] << std::endl;
            LIBMOL::GenCifFile  dataFromCif(AJob.IOEntries["inCifNameB"], std::ios::in);
            LIBMOL::MolGenerator  aMolCreator(dataFromCif, aNBDepth);
            aMolCreator.execute(AJob.IOEntries["userOutName"].c_str());
        }
        else if (AJob.workMode==32)
        {
            
            std::cout << "Input cif " << AJob.IOEntries["inCifName"] << std::endl;
            LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
            LIBMOL::MolGenerator  aMolCreator(dataFromCif, aNBDepth);
            aMolCreator.execute(AJob.IOEntries["userOutName"].c_str());
        }
        else if (AJob.workMode==33)
        {
            std::cout << "The directory of input cif files " 
                      << AJob.IOEntries["inCifNameB"] << std::endl;
            
        }
        
    }
    else if (AJob.workMode == 41)
    {
        LIBMOL::DictCifFile aTargetSystem(AJob.IOEntries["inCifName"], std::ios::in);
        
        if ( (int)aTargetSystem.allAtoms.size() > 0)
        {     
            // CodClassify  aCodSystem(aTargetSystem.allAtoms);
            LIBMOL::CodClassify  aCodSystem(aTargetSystem);
            aCodSystem.codAtomClassify(2);
            aCodSystem.outAtomTypes(AJob.IOEntries["monoRootName"]);
        }
    }
    else if (AJob.workMode == 42)
    {
        LIBMOL::MolSdfFile aMolFile(AJob.IOEntries["inSdfName"],  std::ios::in);
        LIBMOL::CodClassify  aCodSystem(aMolFile.allMols);
        aCodSystem.codAtomClassify(2);
        aCodSystem.outAtomTypes(AJob.IOEntries["monoRootName"]);
    }
    
    return 0;
}

