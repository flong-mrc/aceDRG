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

#ifndef SMILETOOL_H
#include "include/kernel/SmileTool.h"
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

#ifndef ALLSYSTEM_H
#include "AllSystem.h"
#endif

#ifndef CCP4ATOMTYPE_H
#include "CCP4AtomType.h"
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
    
    std::cout << "The work mode is " << AJob.workMode << std::endl;
   
    
    for (std::map<LIBMOL::ID,LIBMOL::ID>::iterator iKW=AJob.IOEntries.begin();
            iKW !=AJob.IOEntries.end(); iKW++)
    {   
        std::cout << iKW->first << "\t" << iKW->second << std::endl;
    }
    
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
                                           
            
            LIBMOL::AllSystem   aTargetSystem(dataFromCif, AJob.IOEntries["libMolTabDir"]); 
        
            if ( (int)aTargetSystem.allAtoms.size() > 0)
            {   
                aTargetSystem.setupAllTargetValuesFromCOD(AJob.IOEntries["userOutName"].c_str(), 
                                                          AJob.IOEntries["monoRootName"], 
                                                          AJob.IOEntries["libMolTabDir"]);
                //if(AJob.IOEntries["NoGeoOpt"].substr(0,1).compare("y")==0)
                //{
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
                //}
                // aTargetSystem.chiralExch();
                
                /*
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
                */
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
                    
                    LIBMOL::AllSystem   aTargetSystem(*iMol, AJob.IOEntries["libMolTabDir"]);
                    
                    if ( (int)aTargetSystem.allAtoms.size() > 0)
                    {
                        std::string tOutName(AJob.IOEntries["userOutName"]);
                        if(dataFromSdf.allMols.size()>1)
                        {
                            tOutName.append("_");
                            tOutName.append(LIBMOL::IntToStr(i));
                        }
                        
                        aTargetSystem.setupAllTargetValuesFromCOD(tOutName.c_str(), 
                                                   AJob.IOEntries["monoRootName"],
                                                   AJob.IOEntries["libMolTabDir"]);
                        
                        /*
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
                        */
                        // aTargetSystem.chiralExch();
                        
                        
                        if ( !aTargetSystem.containMetal()) 
                        {
                             
                            // coordinate generator works only for a system of "organic" atoms. Tempo!
                            GO::FindGlobMin  aGlobMinSystem(aTargetSystem, 
                                                            tOutName.c_str(),
                                                            AJob.IOEntries["monoRootName"]);  
            
                            aGlobMinSystem.Driver();
            
                            //aGlobMinSystem.PreIdealization();
                            
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
    else if (AJob.workMode == 21 || AJob.workMode == 22 )
    {
        
        std::cout <<" 1 " << std::endl;
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"].c_str(),
                                        AJob.IOEntries["inPdbName"].c_str());
        std::cout << " 2 " << std::endl;
        exit(1);
        std::string tOutName(AJob.IOEntries["userOutName"]);
        LIBMOL::outMMCif3Secs(tOutName.c_str(), AJob.IOEntries["monoRootName"],
                      dataFromCif.allAtoms, dataFromCif.allUnchangedBlocks);
                         
        // dataFromCif
        //LIBMOL::CodClassify  aCodSystem(dataFromCif);
        //aCodSystem.getAnglesFromPDB(AJob.IOEntries["inPdbName"]);
        //Temp for outRestraintCif2
        //aCodSystem.outRestraintCif2(AJob.IOEntries["userOutName"].c_str(), AJob.IOEntries["monoRootName"]);
        //aCodSystem.outPDB(AJob.IOEntries["userOutName"].c_str(), AJob.IOEntries["monoRootName"]);
        
       
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
            LIBMOL::CodClassify  aCodSystem(aTargetSystem, AJob.IOEntries["libMolTabDir"]);
            aCodSystem.codAtomClassify(2);
            outCodAndCcp4AtomTypes(AJob.IOEntries["userOutName"].c_str(),
                                   aCodSystem.allAtoms);
            
            LIBMOL::CCP4AtomType  aCPP4TypeTool(aCodSystem.allAtoms, aCodSystem.allRings);
            aCPP4TypeTool.setAllAtomsCCP4Type();
            for (int i=0; i < (int)aCPP4TypeTool.allAtoms.size(); i++)
            {
                if (aCPP4TypeTool.allAtoms[i].ccp4Type.compare(aCodSystem.allAtoms[i].ccp4Type) !=0)
                {
                    std::cout << aCodSystem.allAtoms[i].id << "    " 
                              << aCPP4TypeTool.allAtoms[i].ccp4Type
                              << "     " << aCodSystem.allAtoms[i].ccp4Type << std::endl;
                }
            }
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

