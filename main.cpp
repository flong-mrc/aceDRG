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
                
                //aTargetSystem.setupAllTargetValuesFromCOD(AJob.IOEntries["userOutName"].c_str(), 
                //                                           AJob.IOEntries["monoRootName"], 
                //                                           AJob.IOEntries["libMolTabDir"]);
                
                LIBMOL::outMMCif(AJob.IOEntries["userOutName"].c_str(),
                                 AJob.IOEntries["monoRootName"], 
                                 aTargetSystem.allAtoms,
                                 // aTargetSystem.allHAtomIdx,
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
                                                 aTargetSystem.allAtoms,
                                                 // aTargetSystem.allHAtomIdx,
                                                 aTargetSystem.allBonds,
                                                 aTargetSystem.allAngles,
                                                 aTargetSystem.allTorsions,
                                                 aTargetSystem.allRings,
                                                 aTargetSystem.allPlanes,
                                                 aTargetSystem.allChirals);
                        
                        
                        LIBMOL::outPDB(tOutName.c_str(),
                                AJob.IOEntries["monoRootName"], 
                                aTargetSystem.allAtoms);
                        
                        
                        // aTargetSystem.chiralExch();
                        
                        /*
                        if ( !aTargetSystem.containMetal()) 
                        {
                             
                            // coordinate generator works only for a system of "organic" atoms. Tempo!
                            GO::FindGlobMin  aGlobMinSystem(aTargetSystem, 
                                                            tOutName.c_str(),
                                                            AJob.IOEntries["monoRootName"]);  
            
                            aGlobMinSystem.Driver();
            
                            if (aGlobMinSystem.lGlob)
                            {
                                //aGlobMinSystem.PreIdealization();
                            
                                LIBMOL::outMMCif(tOutName.c_str(),
                                                 AJob.IOEntries["monoRootName"], 
                                                 aGlobMinSystem.allAtoms,
                                                 // aTargetSystem.allHAtomIdx,
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
                                std::cout << "Using initial coords " << std::endl;
                                LIBMOL::outMMCif(tOutName.c_str(),
                                             AJob.IOEntries["monoRootName"], 
                                             aTargetSystem.allAtoms,
                                             // aTargetSystem.allHAtomIdx,
                                             aTargetSystem.allBonds,
                                             aTargetSystem.allAngles,
                                             aTargetSystem.allTorsions,
                                             aTargetSystem.allRings,
                                             aTargetSystem.allPlanes,
                                             aTargetSystem.allChirals);
                            
                                LIBMOL::outPDB(tOutName.c_str(),
                                               AJob.IOEntries["monoRootName"], 
                                               aTargetSystem.allAtoms);
                            }
                           
                        }
                        else
                        {
                            LIBMOL::outMMCif(tOutName.c_str(),
                                             AJob.IOEntries["monoRootName"], 
                                             aTargetSystem.allAtoms,
                                             // aTargetSystem.allHAtomIdx,
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
                         */
                        
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
                
                if ( (int)aTargetSystem.allAtoms.size() > 0)
                {   
                    aTargetSystem.setupAllTargetValuesFromCOD2(AJob.IOEntries["userOutName"].c_str(), 
                                                              AJob.IOEntries["monoRootName"], 
                                                              AJob.IOEntries["libMolTabDir"]);
               
                    LIBMOL::outMMCif(AJob.IOEntries["userOutName"].c_str(),
                                     AJob.IOEntries["monoRootName"], 
                                     aTargetSystem.allAtoms,
                                     // aTargetSystem.allHAtomIdx,
                                     aTargetSystem.allBonds,
                                     aTargetSystem.allAngles,
                                     aTargetSystem.allTorsions,
                                     aTargetSystem.allRings,
                                     aTargetSystem.allPlanes,
                                     aTargetSystem.allChirals);
            
                    LIBMOL::outPDB(AJob.IOEntries["userOutName"].c_str(),
                                   AJob.IOEntries["monoRootName"], 
                                   aTargetSystem.allAtoms);
                }
                else
                {
                    std::cout << "The input MOl2 file " << AJob.IOEntries["inCifName"] 
                              << " contains NO atoms! check the file " << std::endl;
                }
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
            if (dataFromCif.RFactorOK)
            {
                std::cout << "R factor satisfies the requirement" << std::endl;
                LIBMOL::MolGenerator  aMolCreator(dataFromCif, aNBDepth);
                aMolCreator.aLibmolTabDir = AJob.IOEntries["libMolTabDir"];
                aMolCreator.execute(AJob.IOEntries["userOutName"].c_str());
            }
            else
            {
                std::cout << "REJECTED STRUCTURE: R factor related " << std::endl; 
                std::cout << "The data will not be converted to molecules because of : "
                          << std::endl << "(1) high R factors" << std::endl
                          << "or " << std::endl << "(2) no R factors in the data"
                          << std::endl;
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
       
        
        LIBMOL::DictCifFile aTargetSystem(AJob.IOEntries["inCifName"], std::ios::in);
        
        if ( (int)aTargetSystem.allAtoms.size() > 0)
        {     
            // CodClassify  aCodSystem(aTargetSystem.allAtoms);
            std::vector<LIBMOL::AtomDict>  tmpAtomSet;
            
            LIBMOL::CodClassify  aCodSystem(aTargetSystem, AJob.IOEntries["libMolTabDir"]);
            aCodSystem.codAtomClassifyNew2(2);
            
            for (std::vector<LIBMOL::AtomDict>::iterator iAt=aCodSystem.allAtoms.begin();
                    iAt !=aCodSystem.allAtoms.end(); iAt++)
            {
                tmpAtomSet.push_back(*iAt);
            }
            /*
            std::cout <<  "There are " << aCodSystem.allRings.size() << " rings. They are: "
                      << std::endl;
            
            for (std::map<std::string, std::vector<LIBMOL::RingDict> > ::iterator iR1=aCodSystem.allRings.begin();
                    iR1 !=aCodSystem.allRings.end(); iR1++)
            {
                std::cout << "Ring representation " << iR1->first << std::endl;
                for (std::vector<LIBMOL::RingDict>::iterator iR11=iR1->second.begin();
                        iR11 !=iR1->second.end(); iR11++)
                {
                    std::cout << "The ring consists of atoms: " << std::endl;
                    for (std::vector<LIBMOL::AtomDict>::iterator iAt1=iR11->atoms.begin();
                            iAt1 !=iR11->atoms.end(); iAt1++)
                    {
                        std::cout << iAt1->id << std::endl;
                    }
                }
                
                std::cout << std::endl;
                
            }
            */
            LIBMOL::CodClassify  aCodSystem2(aTargetSystem.allAtoms);
            aCodSystem2.setAtomsBondingAndChiralCenter();
            aCodSystem2.codAtomClassifyNew2(2);
            
            std::cout << "Different symbol system now: " << std::endl;
            for (unsigned i=0; i < aCodSystem.allAtoms.size(); i++)
            {
                if (aCodSystem.allAtoms[i].codClass.compare(tmpAtomSet[i].codClass)
                        !=0)
                {
                    std::cout << "Atom " << aCodSystem.allAtoms[i].id 
                              << " new Cod type " << aCodSystem2.allAtoms[i].codClass
                              << " old Cod type " << tmpAtomSet[i].codClass 
                              << std::endl;
                }
            }
            
            
            
            /*
            LIBMOL::ringTools aRingTool;
            aTargetSystem.allRings.clear();
            int nMaxRing = 7, nDep=2;
            aRingTool.detectRingFromAtoms(aTargetSystem.allAtoms,
                                          aTargetSystem.allRings, nDep, nMaxRing);
            
            aTargetSystem.allRingsV.clear();
            
            
            std::cout <<  "(2) Different tool. There are " 
                      << aTargetSystem.allRings.size() << " rings. They are: "
                      << std::endl;
            
            for (std::map<std::string, std::vector<LIBMOL::RingDict> > ::iterator iR1=aTargetSystem.allRings.begin();
                    iR1 !=aTargetSystem.allRings.end(); iR1++)
            {
                
                std::cout << "(2)Ring representation " << iR1->first << std::endl;
                for (std::vector<LIBMOL::RingDict>::iterator iR11=iR1->second.begin();
                        iR11 !=iR1->second.end(); iR11++)
                {
                    aTargetSystem.allRingsV.push_back(*iR11);
                    std::cout << "The ring consists of atoms: " << std::endl;
                    for (std::vector<LIBMOL::AtomDict>::iterator iAt1=iR11->atoms.begin();
                            iAt1 !=iR11->atoms.end(); iAt1++)
                    {
                        std::cout << iAt1->id << std::endl;
                    }
                }
                
                std::cout << std::endl;
                
            }
            
            if (aTargetSystem.allRingsV.size())
            {
                LIBMOL::checkAndSetupPlanes(aTargetSystem.allRingsV, aTargetSystem.allPlanes, aTargetSystem.allAtoms);
            }
            
           
            
            aRingTool.setAtomsRingRepreS(aTargetSystem.allAtoms,
                                         aTargetSystem.allRingsV);
            
            
            */
            
            LIBMOL::outCodAndCcp4AtomTypes(AJob.IOEntries["userOutName"].c_str(),
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
            
          
            if (aCodSystem.allRings.size() != 0)
            {
                std::vector<LIBMOL::RingDict> aRingSys;
                for (std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> >::iterator iMR=
                        aCodSystem.allRings.begin(); iMR !=aCodSystem.allRings.end(); iMR++)
                {
                    for (std::vector<LIBMOL::RingDict>::iterator iR=iMR->second.begin();
                            iR !=iMR->second.end(); iR++)
                    {
                        aRingSys.push_back(*iR);
                    }
                }
                
                // LIBMOL::TestAromaticity(aRingSys, aCodSystem.allAtoms);
                LIBMOL::checkAndSetupPlanes(aRingSys, aCodSystem.allPlanes, aCodSystem.allAtoms);
                
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
    else if (AJob.workMode == 900)
    {
        if (AJob.IOEntries.find("codAtomStr") !=AJob.IOEntries.end())
        {
            LIBMOL::CodClassify aTool;
            
            
            
            std::ifstream    aInFile(AJob.IOEntries["codAtomStr"].c_str());
            if (aInFile.is_open())
            {
                int i=0;
                while(!aInFile.eof())
                {
                    std::string aRecord;
                    std::getline(aInFile, aRecord);
                    aRecord = LIBMOL::TrimSpaces(aRecord);
                    std::cout << aRecord << std::endl;
                    if (aRecord.size())
                    {
                        std::vector<std::string> aSetStrs;
                        LIBMOL::StrTokenize(aRecord, aSetStrs);
                        /*
                        if (aSetStrs.size() ==2)
                        {
                            LIBMOL::BondDict aBond;
                            
                            for (unsigned j=0; j < aSetStrs.size(); j++)
                            {
                                LIBMOL::AtomDict aAtom;
                                aAtom.id = "Test_" + LIBMOL::IntToStr(int(i));
                                aAtom.seriNum = i;
                                std::cout << aAtom.id << std::endl;
                                aTool.codClassToAtom2(aSetStrs[j], aAtom);
                                std::vector<LIBMOL::ID> tConns;
                                LIBMOL::StrTokenize(aAtom.codNB2Symb, tConns, ':');
                                for (int k=0; k < (int)tConns.size(); k++)
                                {
                                    aAtom.connAtoms.push_back(k);
                                }
                                aTool.allAtoms.push_back(aAtom);
                                aBond.fullAtoms[aAtom.id] = i;
                                aBond.atomsIdx.push_back(i);
                                
                                i++;
                                
                            }
                            
                            aTool.allBonds.push_back(aBond);
                        }
                         */
                        if (aSetStrs.size() ==3)
                        {
                            LIBMOL::AngleDict aAng;
                            for (unsigned j=0; j < aSetStrs.size(); j++)
                            {
                                LIBMOL::AtomDict aAtom;
                                aAtom.id = "Test_" + LIBMOL::IntToStr(int(i));
                                aAtom.seriNum = i;
                                std::cout << aAtom.id << std::endl;
                                aTool.codClassToAtom2(aSetStrs[j], aAtom);
                                std::vector<LIBMOL::ID> tConns;
                                LIBMOL::StrTokenize(aAtom.codNB2Symb, tConns, ':');
                                for (int k=0; k < (int)tConns.size(); k++)
                                {
                                    aAtom.connAtoms.push_back(k);
                                }
                                aTool.allAtoms.push_back(aAtom);
                                aAng.atoms.push_back(i);
                                i++;
                                
                            }
                            aTool.allAngles.push_back(aAng);
                        }
                        
                    }
                         
                    
                }
                aInFile.close();
                
                std::cout << "number of angles " << aTool.allAngles.size() << std::endl;
                std::cout << "number of atoms " << aTool.allAtoms.size() << std::endl;
                
                
                if (aTool.allAtoms.size() && aTool.allAngles.size())
                {
                    
                    aTool.hashingAtoms2();
                    //aTool.setOrgBondHeadHashList2();
                    //aTool.groupCodOrgBonds2();
                    //aTool.setOrgAngleHeadHashList22();
                    
                    
                    aTool.groupCodOrgAngles22();
                    aTool.searchCodAngles2();
                    
                    /*
                    for (std::vector<LIBMOL::BondDict>::iterator iBo=aTool.allBonds.begin();
                            iBo !=aTool.allBonds.end(); iBo++)
                    {
                        aTool.searchCodOrgBonds2(iBo);
                    }
                     */            
                }
            }
            
        }
    }
    
    return 0;
    
}

