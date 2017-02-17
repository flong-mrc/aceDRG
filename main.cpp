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

#ifndef RING_H
#include "ring.h"
#endif

//using namespace GO;
//using namespace FF;

/*
 *  A temporary version of main, some re-organization under way
 */
int main(int argc, char** argv) {
    
    LIBMOL::CheckEnvAndGetMode AJob(argc, argv);
    
    std::cout << "workMode " << AJob.workMode << std::endl;
    std::cout << "user output name " << AJob.IOEntries["userOutName"] << std::endl;
    
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
                 aTargetSystem.setupAllTargetValuesFromCOD2(AJob.IOEntries["userOutName"].c_str(), 
                                                            AJob.IOEntries["monoRootName"], 
                                                            AJob.IOEntries["libMolTabDir"]);
                
                 LIBMOL::outMMCif(AJob.IOEntries["userOutName"].c_str(),
                                  AJob.IOEntries["monoRootName"], 
                                  aTargetSystem.propComp,
                                 aTargetSystem.allAtoms,
                                 // aTargetSystem.allHAtomIdx,
                                 aTargetSystem.allBonds,
                                 aTargetSystem.allAngles,
                                 aTargetSystem.miniTorsions,
                                 aTargetSystem.allRingsV,
                                 aTargetSystem.allPlanes,
                                 aTargetSystem.allChirals);
                 
                LIBMOL::outPDB(AJob.IOEntries["userOutName"].c_str(),
                               AJob.IOEntries["monoRootName"], 
                               aTargetSystem.allAtoms);
                
                LIBMOL::outB_and_A_Levels(AJob.IOEntries["userOutName"].c_str(),
                                          aTargetSystem.allAtoms,
                                          aTargetSystem.allBonds,
                                          aTargetSystem.allAngles);
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
                        // std::cout << "input name" << AJob.IOEntries["userOutName"] << std::endl;
                        
                        
                        if(dataFromSdf.allMols.size()>1)
                        {
                            tOutName.append("_");
                            tOutName.append(LIBMOL::IntToStr(i));
                        }
                        
                        
                        aTargetSystem.setupAllTargetValuesFromCOD2(tOutName.c_str(), 
                                                   AJob.IOEntries["monoRootName"],
                                                   AJob.IOEntries["libMolTabDir"]);
                        
                        
                        LIBMOL::outMMCif(tOutName.c_str(),
                                                 AJob.IOEntries["monoRootName"], 
                                                 aTargetSystem.propComp,
                                                 aTargetSystem.allAtoms,
                                                 // aTargetSystem.allHAtomIdx,
                                                 aTargetSystem.allBonds,
                                                 aTargetSystem.allAngles,
                                                 aTargetSystem.allTorsions,
                                                 aTargetSystem.allRingsV,
                                                 aTargetSystem.allPlanes,
                                                 aTargetSystem.allChirals);
                        
                        
                        LIBMOL::outPDB(tOutName.c_str(),
                                AJob.IOEntries["monoRootName"], 
                                aTargetSystem.allAtoms);
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
        else if (AJob.workMode == 14)
        {
            LIBMOL::SYBLMol2File dataFromMol2(AJob.IOEntries["inMol2Name"], std::ios::in);
            
            if (dataFromMol2.atoms.size() >0)
            {
                LIBMOL::AllSystem   aTargetSystem(dataFromMol2, AJob.IOEntries["libMolTabDir"]);
                
                LIBMOL::outMMCif2(AJob.IOEntries["userOutName"].c_str(),
                                     AJob.IOEntries["monoRootName"], 
                                     aTargetSystem.propComp,
                                     aTargetSystem.allAtoms,
                                     // aTargetSystem.allHAtomIdx,
                                     aTargetSystem.allBonds,
                                     aTargetSystem.allAngles,
                                     aTargetSystem.allTorsions,
                                     aTargetSystem.allRingsV,
                                     aTargetSystem.allPlanes,
                                     aTargetSystem.allChirals);
            }
                 
        }
    }
    else if (AJob.workMode == 21 || AJob.workMode == 22 )
    {
        
       
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"].c_str(),
                                        AJob.IOEntries["inPdbName"].c_str());
        
        std::string tOutName(AJob.IOEntries["userOutName"]);
        LIBMOL::outMMCif3Secs(tOutName.c_str(), AJob.IOEntries["monoRootName"],
                      dataFromCif.allAtoms, dataFromCif.allUnchangedBlocks);
                         
       
    }
    else if (AJob.workMode==31 || AJob.workMode==311 || AJob.workMode ==312
             || AJob.workMode==32 || AJob.workMode==33)
    {
        // std::cout << "WorkMode: Molecule generation" << std::endl;
        int aNBDepth = LIBMOL::StrToInt(AJob.IOEntries["NBDepth"]);
        
        if (AJob.workMode==31 || AJob.workMode==311 || AJob.workMode ==312)
        {
            std::cout << "Input cif " << AJob.IOEntries["inCifNameB"] << std::endl;
            LIBMOL::GenCifFile  dataFromCif(AJob.IOEntries["inCifNameB"], std::ios::in);
            std::cout << "Number of crystal in the system is " 
                      << dataFromCif.allCryst.size() << std::endl;
           
            // TEMP, CSD do not provide several parameters. Rely on CSD
            // search criteria. 
            if (dataFromCif.notPowder && dataFromCif.resolOK
                && dataFromCif.RFactorOK && dataFromCif.colidOK 
                && !dataFromCif.hasHeavyCalcAtoms && !dataFromCif.lErr)
            //if(!dataFromCif.lErr)
            {
                std::cout << "The structure is from single crystallographic x-ray "
                          << std::endl;
                std::cout << "R factor satisfies the requirement" << std::endl;
                LIBMOL::MolGenerator  aMolCreator(dataFromCif, aNBDepth);
                std::cout << "workMode " << AJob.workMode << std::endl;
               
                
                aMolCreator.aLibmolTabDir = AJob.IOEntries["libMolTabDir"];
                
                if (AJob.workMode==31 || AJob.workMode==311)
                {
                    aMolCreator.execute(AJob.IOEntries["userOutName"].c_str());
                }
                else if (AJob.workMode ==312)
                {
                    std::cout << "Studies related metal atoms " << std::endl;
                    //std::cout << "Input cif " << AJob.IOEntries["inCifNameB"] << std::endl;
                    // std::cout << "Contain Metal " << dataFromCif.hasMetal << std::endl;
                    
                    if (dataFromCif.hasMetal)
                    {
                        std::cout << "The system contain metal atoms " << std::endl;
                        std::cout << "Those metal atoms are : " << std::endl;
                        for (std::vector<LIBMOL::AtomDict>::iterator iAt=dataFromCif.allAtoms.begin();
                             iAt != dataFromCif.allAtoms.end(); iAt++)
                        {
                            if (iAt->isMetal)
                            {
                                std::cout << "Atom " << iAt->id << " of element " 
                                          << iAt->chemType << std::endl;
                            }
                        }
                
                        LIBMOL::outMetalAtomInfo(AJob.IOEntries["userOutName"].c_str(),
                                        dataFromCif);
                        //std::cout << "The structure is from single crystallographic x-ray "
                        //  << std::endl;
                        //std::cout << "R factor satisfies the requirement" << std::endl;
                
                        aMolCreator.executeMet(AJob.IOEntries["userOutName"].c_str());
                    }
                    else
                    {
                        dataFromCif.errMsg.push_back("The system contains no metal atoms \n");
                        LIBMOL::writeMsgFile(AJob.IOEntries["userOutName"],
                                         dataFromCif.errMsg);
                    }
                }
            }
            else
            {
                if(dataFromCif.errMsg.size() !=0)
                {
                    LIBMOL::writeMsgFile(AJob.IOEntries["userOutName"],
                                         dataFromCif.errMsg);
                }
                
                if (!dataFromCif.notPowder)
                {
                    std::cout << "REJECTED STRUCTURE: Powder Diffraction  " 
                              << std::endl; 
                }
                else if (!dataFromCif.resolOK)
                {
                    std::cout << "REJECTED STRUCTURE: Problems in Resolution or Low resolutions " 
                              << std::endl; 
                }
                else if ( !dataFromCif.RFactorOK)
                {
                    std::cout << "REJECTED STRUCTURE: R factor related " << std::endl; 
                    std::cout << "The data will not be converted to molecules because of : "
                              << std::endl << "(1) high R factors" << std::endl
                              << "or " << std::endl << "(2) no R factors in the data"
                              << std::endl;
                }
                else if (!dataFromCif.colidOK)
                {
                    std::cout << "REJECTED STRUCTURE: Too many atoms have less than 1.0 occp " 
                              << std::endl; 
                }
                else if (dataFromCif.hasHeavyCalcAtoms)
                {
                    std::cout << "REJECTED STRUCTURE: the system contains non-H atoms that"
                              << " are obtained from theoretical calculations " 
                              << std::endl;
                }
                else
                {
                    std::cout << "Check ! the structure has not been converted to " 
                              << "molecules because of unknown reasons !" << std::endl; 
                }
            }
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
       
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
                                           
        // LIBMOL::DictCifFile aTargetSystem(AJob.IOEntries["inCifName"], std::ios::in);
        
        if (dataFromCif.allAtoms.size() > 0)
        {     
            LIBMOL::AllSystem   aTargetSystem(dataFromCif, AJob.IOEntries["libMolTabDir"]); 
            LIBMOL::CodClassify  aCodSystem(aTargetSystem.allAtoms,
                                            aTargetSystem.allHAtomIdx, 
                                            aTargetSystem.allBonds, 
                                            aTargetSystem.allAngles, 
                                            aTargetSystem.allTorsions, 
                                            aTargetSystem.allChirals, 
                                            aTargetSystem.allPlanes, 
                                            aTargetSystem.allRings,
                                            AJob.IOEntries["libMolTabDir"], 2);
            // aCodSystem.codAtomClassifyNew2(2);
            LIBMOL::outAtomTypesAndConnections(AJob.IOEntries["userOutName"].c_str(),
                                               aCodSystem.allAtoms,
                                               aCodSystem.allBonds);
        }
    }
    else if (AJob.workMode == 43)
    {
        LIBMOL::MolSdfFile aMolFile(AJob.IOEntries["inSdfName"].c_str());
        std::cout << "Number of molecules in " << AJob.IOEntries["inSdfName"] 
                  << " is " << aMolFile.allMols.size() << std::endl;
        
        std::map<std::string, std::vector<std::string> >    allAtomTypes;
        std::vector<std::string>                            allBondLines;
        std::vector<std::string>                            allAngleLines;
        
        for (unsigned i=0; i < aMolFile.allMols.size(); i++)
        {
            std::cout << "DEAL WITH MOLECULE " << i+1 << " in " 
                      << AJob.IOEntries["inSdfName"] << std::endl;
            
            LIBMOL::AllSystem   aTargetSystem(aMolFile.allMols[i], AJob.IOEntries["libMolTabDir"]);
            
            LIBMOL::CodClassify  aCodSystem(aTargetSystem.allAtoms,
                                            aTargetSystem.allHAtomIdx, 
                                            aTargetSystem.allBonds, 
                                            aTargetSystem.allAngles, 
                                            aTargetSystem.allTorsions, 
                                            aTargetSystem.allChirals, 
                                            aTargetSystem.allPlanes, 
                                            aTargetSystem.allRings,
                                            AJob.IOEntries["libMolTabDir"], 2);
            
            
            LIBMOL::accumInfoMols(aMolFile.allMols[i].comments[0],
                                  aCodSystem.allAtoms,
                                  aCodSystem.allBonds,
                                  aCodSystem.allAngles,
                                  allAtomTypes, 
                                  allBondLines, 
                                  allAngleLines);
        }
        
        std::string  aAtName = AJob.IOEntries["userOutName"] + "_atomType.txt";
        std::ofstream outAtomTypes(aAtName.c_str());
        for (std::map<std::string,std::vector<std::string> >::iterator iAt=allAtomTypes.begin();
                    iAt != allAtomTypes.end(); iAt++)
        {
             outAtomTypes << iAt->first << "  :  "
                      << iAt->second.size() << std::endl;
        }
        
        std::string  aBaName = AJob.IOEntries["userOutName"] + "_bonds.txt";
        std::ofstream outBonds(aBaName.c_str());
        for (std::vector<std::string>::iterator iB=allBondLines.begin();
              iB != allBondLines.end(); iB++)
        {
            outBonds << *iB;
        }
        
    }
    else if (AJob.workMode == 900)
    {
        if (AJob.IOEntries.find("Type1") !=AJob.IOEntries.end()
            and AJob.IOEntries.find("Type2") !=AJob.IOEntries.end())
        {   
            LIBMOL::isomorGraph graphTool;
            LIBMOL::Graph  g1, g2;
            
            std::vector<std::map<int, int> > matchedList;
            
            graphTool.graphMatch(AJob.IOEntries["Type1"].c_str(), 
                                     g1,
                                     AJob.IOEntries["Type2"].c_str(),
                                     g2,
                                     matchedList,
                                     1);
            if (matchedList.size()>0)
            {
                graphTool.outputMatchedGraphs(g1, g2, matchedList, 1,
                                              AJob.IOEntries["userOutName"].c_str());
                std::cout << "Graph matching done. " << std::endl;
            }
            else
            {
                std::cout << "No matched graphs found. " << std::endl;
            }
        }
    }
    else if (AJob.workMode == 910)
    {
        // Hukel MO applications
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
        LIBMOL::CodClassify aClassifiedSys(dataFromCif, AJob.IOEntries["libMolTabDir"]);
        
        LIBMOL::HuckelMOSuite  aMoTool;
        aMoTool.setWorkMode(2);
        //aMoTool.execute(dataFromCif.allAtoms, dataFromCif.allBonds);
        aMoTool.execute2(aClassifiedSys.allAtoms, aClassifiedSys.allBonds, 
                         aClassifiedSys.allRingsV);   
        if(AJob.IOEntries.find("userOutName") 
             == AJob.IOEntries.end())
        {
            AJob.IOEntries["userOutName"] = "Libmol_atoms_bonds.txt";
        }
        aMoTool.outBoAndChList(AJob.IOEntries["userOutName"].c_str(),
                               aClassifiedSys.allAtoms, aClassifiedSys.allBonds);    
        
    }
    else if (AJob.workMode == 1001)
    {
        
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
                                               
        LIBMOL::AllSystem   aTargetSystem(dataFromCif, AJob.IOEntries["libMolTabDir"]); 
        
        std::cout << "number of rings " << aTargetSystem.allRingsV.size() << std::endl;
        
        if ( (int)aTargetSystem.allAtoms.size() > 0)
        {
            LIBMOL::KekulizeMol aKTool;
            aKTool.execute(aTargetSystem.allAtoms, 
                           aTargetSystem.allBonds,
                           aTargetSystem.allRingsV);
            
            LIBMOL::outBoAndChList(AJob.IOEntries["userOutName"].c_str(),
                                   aTargetSystem.allAtoms, 
                                   aTargetSystem.allBonds);
        }
    }
    
    
    return 0;
    
}

