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

#ifndef PERIODICTABLE_H
#include "include/kernel/periodicTable.h"
#endif

#ifndef GLOBOPT_H
#include "include/go/GlobOpt.h"
#endif

#ifndef MOLGENERATOR_H
#include "MolGenerator.h"
#endif

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
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

#ifdef _MSC_VER
#include <ciso646>
#endif


//using namespace GO;
//using namespace FF;

/*
 *  A temporary version of main, some re-organization under way
 */
int main(int argc, char** argv) {

    LIBMOL::CheckEnvAndGetMode AJob(argc, argv);

    std::cout << "workMode " << AJob.workMode << std::endl;
    //std::cout << "user output name " << AJob.IOEntries["userOutName"]
    //          << std::endl;


    //for (std::map<LIBMOL::ID,LIBMOL::ID>::iterator iKW=AJob.IOEntries.begin();
    //        iKW !=AJob.IOEntries.end(); iKW++)
    //{
    //    std::cout << iKW->first << "\t" << iKW->second << std::endl;
    //}

    if (AJob.workMode == 11
        || AJob.workMode ==12
        || AJob.workMode ==13
        || AJob.workMode ==14)
    {
        if (AJob.workMode == 11 || AJob.workMode ==12 ||AJob.workMode==1116)
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

            LIBMOL::AllSystem   aTargetSystem(dataFromCif,
                                              AJob.IOEntries["libMolTabDir"],
                                              AJob.upperBondSig,
                                              AJob.lowBondSig,
                                              AJob.upperAngleSig,
                                              AJob.lowAngleSig);

            if (AJob.IOEntries.find("inMC") !=  AJob.IOEntries.end())
            {
                aTargetSystem.getMetaLConn(AJob.IOEntries["inMC"].c_str());
                LIBMOL::setAtomsBondingAndChiralCenter(aTargetSystem.allAtoms);
            }

            if (AJob.IOEntries.find("modifiedPl") !=  AJob.IOEntries.end())
            {
                aTargetSystem.lMdPls = true;
            }

            if (AJob.IOEntries.find("PeptidesOnly") !=  AJob.IOEntries.end())
            {
                aTargetSystem.isPeptide = true;
            }

            if ( (int)aTargetSystem.allAtoms.size() > 0 && AJob.workMode != 11162)
            {

                 aTargetSystem.setupAllTargetValuesFromCOD2(
                                AJob.IOEntries["userOutName"].c_str(),
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
                                  aTargetSystem.allChirals,
                                  aTargetSystem.upperBondSig,
                                  aTargetSystem.lowBondSig,
                                  aTargetSystem.upperAngleSig,
                                  aTargetSystem.lowAngleSig,
                                  aTargetSystem.HydrDistTable);


                LIBMOL::outPDB(AJob.IOEntries["userOutName"].c_str(),
                               AJob.IOEntries["monoRootName"],
                               aTargetSystem.allAtoms, 1);

                LIBMOL::outB_and_A_Levels(AJob.IOEntries["userOutName"].c_str(),
                                          aTargetSystem.allAtoms,
                                          aTargetSystem.allBonds,
                                          aTargetSystem.allAngles);
            }
            else if (AJob.workMode == 11162)
            {

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
                                  aTargetSystem.allChirals,
                                  aTargetSystem.upperBondSig,
                                  aTargetSystem.lowBondSig,
                                  aTargetSystem.upperAngleSig,
                                  aTargetSystem.lowAngleSig,
                                  aTargetSystem.HydrDistTable);
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
                                                 aTargetSystem.allChirals,
                                                 aTargetSystem.upperBondSig,
                                                 aTargetSystem.lowBondSig,
                                                 aTargetSystem.upperAngleSig,
                                                 aTargetSystem.lowAngleSig,
                                                 aTargetSystem.HydrDistTable);


                        LIBMOL::outPDB(tOutName.c_str(),
                                AJob.IOEntries["monoRootName"],
                                aTargetSystem.allAtoms, 1);
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
    else if (AJob.workMode==31
             || AJob.workMode==311
             || AJob.workMode==3111
             || AJob.workMode ==312
             || AJob.workMode ==313
             || AJob.workMode ==314
             || AJob.workMode ==315
             || AJob.workMode==32
             || AJob.workMode==33)
    {
        //std::cout << "WorkMode: Molecule generation" << std::endl;
        int aNBDepth = LIBMOL::StrToInt(AJob.IOEntries["NBDepth"]);

        if (   AJob.workMode==31
            || AJob.workMode==311
            || AJob.workMode==3111
            || AJob.workMode ==312
            || AJob.workMode ==313
            || AJob.workMode ==314
            || AJob.workMode ==315
            || AJob.workMode ==316)
        {
            /*
            std::cout << "Input cif " << AJob.IOEntries["inCifNameB"] << std::endl;
            if (AJob.IOEntries.find("MolGenParaFile")==
                AJob.IOEntries.end())
            {
                LIBMOL::GenCifFile
                dataFromCif(AJob.IOEntries["inCifNameB"], std::ios::in);
                std::cout << "Number of crystal in the system is "
                      << dataFromCif.allCryst.size() << std::endl;
            }
             */

            LIBMOL::GenCifFile
            dataFromCif(AJob.IOEntries["inCifNameB"],
                        AJob.IOEntries["MolGenParaFile"],
                        std::ios::in);
            std::cout << "Number of crystal in the system is "
                      << dataFromCif.allCryst.size() << std::endl;

            // TEMP, CSD do not provide several parameters. Rely on CSD
            // search criteria.

            if (dataFromCif.allAtoms.size() > 0)
            {
                dataFromCif.outAtomElems(AJob.IOEntries["userOutName"]);
            }


            //if (dataFromCif.notPowder && dataFromCif.resolOK
            //    && dataFromCif.RFactorOK && dataFromCif.colidOK
            //    && !dataFromCif.hasHeavyCalcAtoms && !dataFromCif.lErr)
            //if(!dataFromCif.lErr)

            if (dataFromCif.checkOverAll(AJob.workMode))
            {
                if (AJob.workMode ==3111)
                {
                    std::cout << "The structure is from neutron diffractions "
                              << std::endl;
                }
                else
                {
                    //std::cout << "The structure is from single crystallographic x-ray "
                    //          << std::endl;
                    //std::cout << "R factor satisfies the requirement" << std::endl;
                    std::cout << "Pass temporarily" << std::endl;
                }

                //std::cout << "workMode " << AJob.workMode << std::endl;
                //std::cout << "LibmolTabDir " << AJob.IOEntries["libMolTabDir"]
                //          << std::endl;

                dataFromCif.outCrystInfo(AJob.IOEntries["userOutName"]);

                LIBMOL::MolGenerator  aMolCreator(dataFromCif, aNBDepth);

                if (AJob.IOEntries.find("libMolTabDir")
                      != AJob.IOEntries.end())
                {
                    aMolCreator.aLibmolTabDir = AJob.IOEntries["libMolTabDir"];
                }


                if (AJob.workMode==31 || AJob.workMode==311)
                {
                    aMolCreator.execute1(AJob.IOEntries["UserParaFile"].c_str(),
                                         AJob.IOEntries["userOutName"].c_str());
                }
                else if (AJob.workMode==3111)
                {
                    aMolCreator.executeNeuD(AJob.IOEntries["userOutName"].c_str());
                }
                else if (AJob.workMode ==312 || AJob.workMode==313)
                {
                    std::cout << "Studies on metal atoms " << std::endl;

                    if (dataFromCif.hasMetal)
                    {
                        std::cout << "The system contain metal atoms " << std::endl;
                        std::cout << "Those metal atoms are : " << std::endl;
                        for (std::vector<LIBMOL::AtomDict>::iterator
                             iAt=dataFromCif.allAtoms.begin();
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
                        if (AJob.workMode==312)
                        {
                            aMolCreator.executeMet(AJob.IOEntries["UserParaFile"].c_str(),
                                                   AJob.IOEntries["userOutName"].c_str());

                            if(dataFromCif.errMsg.size() !=0)
                            {
                                dataFromCif.reOrdErrMsg();
                                if (aMolCreator.allMsg.size() !=0)
                                {
                                    for (std::vector<std::string>::iterator
                                        iErr=aMolCreator.allMsg.begin();
                                        iErr!=aMolCreator.allMsg.end(); iErr++)
                                        {
                                            dataFromCif.errMsg.push_back(*iErr);
                                        }

                                }
                                LIBMOL::writeMsgFile(AJob.IOEntries["userOutName"],
                                                   dataFromCif.errMsg);
                                for (std::vector<std::string>::iterator
                                     iErr=dataFromCif.errMsg.begin();
                                     iErr!=dataFromCif.errMsg.end(); iErr++)
                                {
                                    std::cout <<(*iErr) << std::endl;
                                }
                            }
                        }
                        else if (AJob.workMode==313)
                        {
                            aMolCreator.executeMetRange(
                                        AJob.IOEntries["UserParaFile"].c_str(),
                                        AJob.IOEntries["userOutName"].c_str());
                        }
                    }
                    else
                    {
                        dataFromCif.errMsg.push_back("The system contains no metal atoms \n");
                        LIBMOL::writeMsgFile(AJob.IOEntries["userOutName"],
                                             dataFromCif.errMsg);
                    }
                }
                else if (AJob.workMode==314)
                {

                    std::cout << "Studies of neighbor distribution of "
                                "certain non-metal atoms " << std::endl;
                    /*
                    LIBMOL::PeriodicTable aPTab;
                    double aDDelta = 0.3;


                    if (AJob.IOEntries.find("distDelta")
                            !=AJob.IOEntries.end())
                    {
                        aDDelta =
                        LIBMOL::StrToReal(AJob.IOEntries["distDelta"]);
                    }
                    */

                    aMolCreator.executeSelectedAtomRange(
                                        AJob.IOEntries["UserParaFile"].c_str(),
                                        AJob.IOEntries["userOutName"].c_str());
                }
                else if (AJob.workMode==315)
                {
                    std::cout << "Detect possible H-bond candidates "
                              << std::endl;


                }
            }
            else
            {
                if(dataFromCif.errMsg.size() !=0)
                {
                    dataFromCif.reOrdErrMsg();
                    LIBMOL::writeMsgFile(AJob.IOEntries["userOutName"],
                                         dataFromCif.errMsg);
                    for (std::vector<std::string>::iterator
                         iErr=dataFromCif.errMsg.begin();
                         iErr!=dataFromCif.errMsg.end(); iErr++)
                    {
                        std::cout <<(*iErr) << std::endl;
                    }
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
            aMolCreator.execute(AJob.IOEntries["UserParaFile"].c_str(),
                                AJob.IOEntries["userOutName"].c_str());
        }
        else if (AJob.workMode==33)
        {
            std::cout << "The directory of input cif files "
                      << AJob.IOEntries["inCifNameB"] << std::endl;
        }

    }
    else if (AJob.workMode == 350 || AJob.workMode == 351)
    {
        if (AJob.IOEntries.find("inPdbName") !=AJob.IOEntries.end())
        {

            std::cout << "Input coordinate file "
                      << AJob.IOEntries["inPdbName"] << std::endl;
            std::cout << "Output dict file is " << AJob.IOEntries["userOutName"] << std::endl;
            LIBMOL::PDBFile dataFromPdb(AJob.IOEntries["inPdbName"].c_str(), std::ios::in);
            LIBMOL::MolGenerator  aMolCreator(dataFromPdb);
            aMolCreator.executePdb(AJob.IOEntries["userOutName"].c_str(), AJob.workMode);
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

            //aCodSystem.codAtomClassifyNew2(2);
            LIBMOL::setAtomFormTypes(aCodSystem.allAtoms);
            LIBMOL::outAtomTypesAndConnections(AJob.IOEntries["userOutName"].c_str(),
                                               AJob.IOEntries["monoRootName"],
                                               aCodSystem.allAtoms,
                                               aCodSystem.allBonds,
                                               aCodSystem.allRingsV);

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
    else if (AJob.workMode == 800)
    {


        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);

        LIBMOL::AllSystem   aTargetSystem(dataFromCif,
                                          AJob.IOEntries["libMolTabDir"],
                                          AJob.upperBondSig,
                                          AJob.lowBondSig,
                                          AJob.upperAngleSig,
                                          AJob.lowAngleSig);


        LIBMOL::outProElecDistances(AJob.IOEntries["userOutName"].c_str(),
                                    aTargetSystem);
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
    else if (AJob.workMode == 920 || AJob.workMode == 930)
    {
        LIBMOL::DictCifFile dataFromCif(AJob.IOEntries["inCifName"], std::ios::in);
        LIBMOL::setAtomRingProps(dataFromCif.allAtoms, dataFromCif.allRingsV);
        std::cout << "number of rings " << dataFromCif.allRingsV.size() << std::endl;

        if(dataFromCif.allAtoms.size()> 0 && dataFromCif.allBonds.size() > 0)
        {
            if (AJob.workMode == 920 )
            {
                std::cout << "Kekulize the molecule " << std::endl;
                LIBMOL::KekulizeMol aKTool;
                std::map<std::string, int>               aHMap;
                aKTool.execute(dataFromCif.allAtoms,
                               dataFromCif.allBonds,
                               dataFromCif.allRingsV,
                               aHMap);
                aKTool.outBondsAndHAtms(dataFromCif.allBonds, aHMap, AJob.IOEntries["userOutName"]);
                //std::cout << "Kekulize done " << std::endl;
            }
            else if (AJob.workMode == 930 )
            {
                std::cout << "Kekulize the molecule and add charges. " << std::endl;
                LIBMOL::KekulizeMol aKTool;
                aKTool.executeBC(dataFromCif.allAtoms,
                               dataFromCif.allBonds,
                               dataFromCif.allRingsV);
                std::string aName = AJob.IOEntries["userOutName"] + ".cif";
                outMMCif3(aName.c_str(),
                         AJob.IOEntries["monoRootName"],
                         dataFromCif.allAtoms,
                         dataFromCif.allBonds,
                         dataFromCif.allRingsV);

                aKTool.outBandC(AJob.IOEntries["userOutName"].c_str(),
                         dataFromCif.allAtoms,
                         dataFromCif.allBonds);
            }
        }

    }
    else if (AJob.workMode == 940)
    {

    }
    /*
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
    */
    else if(AJob.workMode==1002)
    {
        std::cout << "Input cif "
                  << AJob.IOEntries["inCifNameB"] << std::endl;
        LIBMOL::GenCifFile  dataFromCif(AJob.IOEntries["inCifNameB"], std::ios::in);
        std::cout << "Number of H atoms " << dataFromCif.allHAtomIdx.size()
                  << std::endl;

        if (dataFromCif.allHAtomIdx.size() >0)
        {

            AJob.IOEntries["userOutName"] = "HasH_" +
                                            AJob.IOEntries["userOutName"]
                                          + ".msg";
            std::cout << AJob.IOEntries["userOutName"] << std::endl;
            std::ofstream aHasHF(AJob.IOEntries["userOutName"].c_str());

            if (aHasHF.is_open())
            {
                for (std::vector<int>::iterator
                     iH=dataFromCif.allHAtomIdx.begin();
                     iH !=dataFromCif.allHAtomIdx.end(); iH++)
                {
                    aHasHF << dataFromCif.allAtoms[*iH].id << std::endl;
                }

                aHasHF.close();
            }
        }
        else
        {
            AJob.IOEntries["userOutName"] = "NoH_" +
                                            AJob.IOEntries["userOutName"]
                                          + ".msg";
            std::cout << AJob.IOEntries["userOutName"] << std::endl;

            std::ofstream aNoHF(AJob.IOEntries["userOutName"].c_str());

        }
    }

    return 0;

}
