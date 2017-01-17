/* 
 * File:   MolGenerator.cpp
 * Author: flong
 *
 * Generator a molecule from input sets of atoms and symmetry operators 
 * 
 * Created on August 13, 2013, 5:56 PM
 */

#include "MolGenerator.h"
#include "codClassify.h"

namespace LIBMOL {

    MolGenerator::MolGenerator() : myNBDepth(1),
    lColid(false),
    lHasMetal(false) {
    }

    MolGenerator::MolGenerator(const GenCifFile& tCifObj, int tNBDepth) {

        myNBDepth = tNBDepth;
        lColid = false;
        lHasMetal = false;

        // std::cout << "myNBDepth " << myNBDepth << std::endl;

        // std::cout << "Number of crystals " << tCifObj.allCryst.size() << std::endl;


        for (std::vector<CrystInfo>::const_iterator iC = tCifObj.allCryst.begin();
                iC != tCifObj.allCryst.end(); iC++) {
            allCryst.push_back(*iC);
        }

        for (std::vector<AtomDict>::const_iterator iA = tCifObj.allAtoms.begin();
                iA != tCifObj.allAtoms.end(); iA++) {
            initAtoms.push_back(*iA);
        }
        // std::cout << "number of initial atoms " << initAtoms.size() << std::endl;

        for (std::vector<std::string>::const_iterator iM = tCifObj.errMsg.begin();
                iM != tCifObj.errMsg.end(); iM++) {
            allMsg.push_back(*iM);
        }

    }

    MolGenerator::MolGenerator(const DictCifFile & tCifObj, int tNBDepth) {
        myNBDepth = tNBDepth;

        for (std::vector<AtomDict>::const_iterator iA = tCifObj.allAtoms.begin();
                iA != tCifObj.allAtoms.end(); iA++) {
            initAtoms.push_back(*iA);
        }

        // setAllAtomEXcessElectrons(initAtoms);

    }

    MolGenerator::~MolGenerator() {
    }

    void MolGenerator::checkMetal(std::vector<ID>& tMeTab) {
        lHasMetal = false;

        for (std::vector<AtomDict>::iterator iAt = initAtoms.begin();
                iAt != initAtoms.end(); iAt++) {
            if (isMetal(tMeTab, iAt->chemType)) {
                lHasMetal = true;
                break;
            }
        }
    }

    void MolGenerator::execute(FileName tOutName) {
        CCP4DictParas tCCP4EnergParas;
        ccp4DictParas.push_back(tCCP4EnergParas);
        PeriodicTable aPTable;


        // std::cout << "number of crystal " << allCryst.size() << std::endl;
        std::cout << "number of atom read in " << std::endl;

        if (initAtoms.size() > 1) {
            for (std::vector<CrystInfo>::iterator iCryst = allCryst.begin();
                    iCryst != allCryst.end(); iCryst++) {
                for (std::vector<AtomDict>::iterator iA = initAtoms.begin();
                        iA != initAtoms.end(); iA++) {
                    //FractToOrtho(iA->fracCoords, iA->coords, iCryst->itsCell->a,
                    //            iCryst->itsCell->b, iCryst->itsCell->c, iCryst->itsCell->alpha,
                    //            iCryst->itsCell->beta, iCryst->itsCell->gamma);
                    if (iA->ocp < 1.0000001) {
                        //packAtomIntoCell((*iA));
                        iA->sId = "555";
                        allAtoms.push_back(*iA);
                        refAtoms.push_back(*iA);
                        // std::cout << "Is in preCell " << iA->isInPreCell << std::endl;
                    }
                }


                //outPDB("initAtoms.pdb", "UNL", initAtoms);

                std::cout << "Number of atoms read from the input file "
                        << initAtoms.size() << std::endl;

                symmAtomGen(iCryst, aPTable);

                std::cout << "number of atoms in a center unit cell "
                        << allAtoms.size() << std::endl;


                if (!lColid) 
                {
                    buildRefAtoms(iCryst);
                    // getUniqueBonds(aPTable);
                    getUniqueAtomLinks(aPTable, iCryst);

                    getMolByEqClassInCell();

                    buildAndValidMols(aPTable, iCryst);

                    // getAtomTypeMols();

                    std::vector<Molecule> aSetOfFiniteMols, aSetOfInfMols;
                    checkInfMols(aSetOfInfMols, aSetOfFiniteMols);

                    allMolecules.clear();
                    for (std::vector<Molecule>::iterator iMol
                            = aSetOfFiniteMols.begin();
                            iMol != aSetOfFiniteMols.end(); iMol++) 
                    {
                        allMolecules.push_back(*iMol);
                    }

                    HuckelMOSuite aMoTool;
                    aMoTool.setWorkMode(2);
                    for (std::vector<Molecule>::iterator iMol
                            = allMolecules.begin();
                            iMol != allMolecules.end(); iMol++) 
                    {
                        std::cout << std::endl
                                << "Huckel MO related section for molecules "
                                << iMol->seriNum << std::endl;


                        aMoTool.execute2(iMol->atoms, iMol->allBonds, iMol->rings);
                        /*
                        for (std::vector<BondDict>::iterator iBo=iMol->allBonds.begin();
                                iBo !=iMol->allBonds.end(); iBo++)
                        {
                            std::cout << "Now bond " << iBo->seriNum 
                                      << " has order " << iBo->orderN << std::endl;
                        }
                         */
                    }

                    getOverallBondAndAnglesNew();
                    outTables(tOutName, allMolecules, aSetOfInfMols);
                } else {
                    outMsg(tOutName);
                }
            }
        }
    }

    void MolGenerator::executeMet(FileName tOutName) 
    {
        if (initAtoms.size() > 0) 
        {
            CCP4DictParas tCCP4EnergParas;
            ccp4DictParas.push_back(tCCP4EnergParas);
            PeriodicTable aPTable;

            std::vector<ID> allMetals;
            initMetalTab(allMetals);
            checkMetal(allMetals);

            if (lHasMetal) 
            {
                for (std::vector<CrystInfo>::iterator iCryst = allCryst.begin();
                        iCryst != allCryst.end(); iCryst++) 
                {
                    for (std::vector<AtomDict>::iterator iA = initAtoms.begin();
                            iA != initAtoms.end(); iA++) 
                    {
                        iA->sId = "555";
                        allAtoms.push_back(*iA);
                        refAtoms.push_back(*iA);
                        // std::cout << "Is in preCell " << iA->isInPreCell << std::endl;
                    }

                    //outPDB("initAtoms.pdb", "UNL", initAtoms);

                    std::cout << "Number of atoms read from the input file "
                            << initAtoms.size() << std::endl;
                    std::cout << "Number of atoms before symmetrical ops is "
                            << allAtoms.size() << std::endl;
                    symmAtomGen(iCryst, aPTable);

                    std::cout << "number of atoms in a center unit cell "
                            << allAtoms.size() << std::endl;
                   
                    if (!lColid) 
                    {
                        
                        buildRefAtoms(iCryst);

                        // getUniqueBonds(aPTable);
                        getUniqueAtomLinks(aPTable, iCryst);
                        // No molecules will be generated
                        /*
                        buildMetalAtomCoordMap(iCryst);
                        outMetalAtomCoordInfo(tOutName);
                        */
                        // Try the new method related the new class "Metalcluster"
                        buildMetalClusters(iCryst);
                        outMetalClusterInfo(tOutName);

                    }
                }
            } 
            else 
            {
                std::cout << "No metal atoms in the assym unit "
                        << std::endl;
            }

        }




    }

    void MolGenerator::symmAtomGen(std::vector<CrystInfo>::iterator tCrys,
            PeriodicTable & tPTable) {
        for (std::vector<AtomDict>::iterator iA = initAtoms.begin();
                iA != initAtoms.end(); iA++) {
            // int tDictMult = iA->symmMult;
            // iA->symmMult =1;
            unsigned nPre = allAtoms.size();

            if (iA->ocp < 1.000001) {
                for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator
                    iOp = tCrys->itsSpaceGroup->sgOp.begin();
                        iOp != tCrys->itsSpaceGroup->sgOp.end(); iOp++) {
                    getOneSymmAtom(iA, iOp, tCrys, tPTable);
                }
            }
            unsigned nPro = allAtoms.size();

            unsigned nDiff = nPro - nPre;
            // std::cout << nDiff << " symm-generated atoms are added " << std::endl;
        }
        std::cout << "number of atoms in Assym " << initAtoms.size() << std::endl;
        std::cout << "number of ops " << tCrys->itsSpaceGroup->sgOp.size() << std::endl;
        std::cout << "number of atoms in the unit cell " << allAtoms.size() << std::endl;
        
    }

    void MolGenerator::getOneSymmAtom(std::vector<AtomDict>::iterator tCurAtom,
            std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
            std::vector<CrystInfo>::iterator tCryst,
            PeriodicTable & tPTab) {
        std::vector<REAL> startFracCoords, endFracCoords;
        for (unsigned i = 0; i < 3; i++) {
            startFracCoords.push_back(tCurAtom->fracCoords[i]);
            endFracCoords.push_back(0.0);
        }
        startFracCoords.push_back(1.0);
        endFracCoords.push_back(0.0);

        matMultVec(tOp->second, startFracCoords, endFracCoords);
        //TranslateIntoUnitCell(endFracCoords, tFracCoords);
        /*
        std::cout << "effects of symm operators" << std::endl;
        std::cout << "initial fract coords " << std::endl;
        for (std::vector<REAL>::iterator iX=startFracCoords.begin();
                iX != startFracCoords.end(); iX++)
        {
            std::cout << *iX << std::endl;
        }
        
        std::cout << "Opt matrix: " << std::endl;
        for (std::vector<std::vector<REAL> >::iterator iRow=tOp->second.begin();
                iRow !=tOp->second.end(); iRow++)
        {
            for (std::vector<REAL>::iterator iEl=iRow->begin();
                    iEl !=iRow->end(); iEl++)
            {
                std::cout << *iEl << "\t";
            }
            std::cout << std::endl;
        }
       
            
        std::cout << "symm-site fract coords" << std::endl;
        for (std::vector<REAL>::iterator iX=endFracCoords.begin();
                iX != endFracCoords.end(); iX++)
        {
            std::cout << *iX << std::endl;
        }
         */

        if (!colidAtom(endFracCoords, allAtoms, 0)) {
            //std::cout << "The symmetrically generated atom should be included in the unit cell "
            //          << std::endl;

            AtomDict tAtom(*tCurAtom);
            tAtom.isInPreCell = false;
            tAtom.symmOp = tOp->first;
            tAtom.fromOrig = tCurAtom->seriNum;
            tAtom.coords.clear();
            tAtom.fracCoords.clear();
            //std::cout << tAtom.id << std::endl;
            //std::cout << "symm atom ocp " << tAtom.ocp << std::endl;

            
            for (unsigned ix = 0; ix < endFracCoords.size() - 1; ix++) {
                tAtom.fracCoords.push_back(endFracCoords[ix]);
            }

            //FractToOrtho(endFracCoords, tAtom.coords, tCryst->itsCell->a,
            //         tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
            //         tCryst->itsCell->beta, tCryst->itsCell->gamma);


            tAtom.seriNum = (int) allAtoms.size();
            packAtomIntoCell(tAtom);
            if (!colidAtom(tAtom, allAtoms, 0)) {
                tAtom.sId = "555";
                allAtoms.push_back(tAtom);
                refAtoms.push_back(tAtom);
            }
        }


    }

    void MolGenerator::packAtomIntoCell(AtomDict & tAtm) {



        for (std::vector<REAL>::iterator iX = tAtm.fracCoords.begin();
                iX != tAtm.fracCoords.end(); iX++) {
            if ((*iX) > 1.00000001) {
                *iX -= 1.0;
            } else if ((*iX) < -0.00000001) {
                *iX += 1;
            }
        }
    }

    void MolGenerator::buildRefAtoms(std::vector<CrystInfo>::iterator iCryst) {
        REAL fraCell = 1.0, tCM = 1000000.0;
        REAL aShell = 7.0;
        if (iCryst->itsCell->a < tCM) {
            tCM = iCryst->itsCell->a;
        }

        if (iCryst->itsCell->b < tCM) {
            tCM = iCryst->itsCell->b;
        }

        if (iCryst->itsCell->c < tCM) {
            tCM = iCryst->itsCell->c;
        }

        if (fabs(tCM - 1000000.0) < 0.1) {
            std::cout << "No cell lengths or huge values of cell lengths. Program stops "
                    << std::endl;
            exit(1);
        }

        getFracReal(tCM, fraCell, aShell);

        // std::cout << "fraCell " << fraCell << std::endl;


        for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                iA != allAtoms.end(); iA++) {
            for (int i = -myNBDepth; i < myNBDepth + 1; i++) {
                for (int j = -myNBDepth; j < myNBDepth + 1; j++) {
                    for (int k = -myNBDepth; k < myNBDepth + 1; k++) {
                        if (!(i == 0 && j == 0 && k == 0)) {

                            std::vector<REAL> tFx;
                            tFx.push_back(iA->fracCoords[0] + i);
                            tFx.push_back(iA->fracCoords[1] + j);
                            tFx.push_back(iA->fracCoords[2] + k);

                            // Take a zone surround the center cell 
                            //REAL tRange1 = -fraCell, tRange2 = 1 + fraCell;
                            //REAL tRange1 = -0.8, tRange2 = 1.8;
                            //if ((tFx[0] > tRange1 && tFx[0] <= tRange2)
                            //        && (tFx[1] > tRange1 && tFx[1] <= tRange2)
                            //        && (tFx[2] > tRange1 && tFx[2] <= tRange2)
                            //        && !colidAtom(tFx, refAtoms, 0)) {
                            if (!colidAtom(tFx, refAtoms, 0)) {
                                AtomDict aAtom(*iA);
                                aAtom.isInPreCell = false;

                                aAtom.sId = IntToStr(5 + i) + IntToStr(5 + j) + IntToStr(5 + k);
                                aAtom.symmOp = iA->symmOp;
                                aAtom.seriNum = (int) refAtoms.size();
                                aAtom.fromOrig = iA->fromOrig;

                                aAtom.fracCoords[0] = tFx[0];
                                aAtom.fracCoords[1] = tFx[1];
                                aAtom.fracCoords[2] = tFx[2];
                                //FractToOrtho(aAtom.fracCoords, aAtom.coords, 
                                //             iCryst->itsCell->a, iCryst->itsCell->b, 
                                //             iCryst->itsCell->c, iCryst->itsCell->alpha,
                                //             iCryst->itsCell->beta, iCryst->itsCell->gamma);
                                refAtoms.push_back(aAtom);
                            }
                        }
                        // std::cout << "size of refAtoms " << refAtoms.size() << std::endl;
                    }
                }
            }
        }

        swithAtoms(iCryst);

        std::cout << "total number of atoms in system  "
                << allAtoms.size() << std::endl;

        /*
        for (std::vector<AtomDict>::iterator iRA1=allAtoms.begin();
                iRA1 != allAtoms.end(); iRA1++)
        {
            if (iRA1->id == "Mo1" && iRA1->seriNum==0 )
            {
                for (std::vector<AtomDict>::iterator iRA2=allAtoms.begin();
                     iRA2 != allAtoms.end(); iRA2++)
                {
                    if (iRA2->id == "Mo1")
                    {
                        // packAtomIntoCell(*iRA2);
                        REAL rD = getBondLenFromFracCoords(iRA1->fracCoords, 
                                                           iRA2->fracCoords,
                        iCryst->itsCell->a, iCryst->itsCell->b,
                        iCryst->itsCell->c, iCryst->itsCell->alpha,
                        iCryst->itsCell->beta, iCryst->itsCell->gamma);
                        //std::cout << "Dist is " << rD << std::endl; 
                    }
                }
            }
        }
         */
    }

    void MolGenerator::addOneSetRefAtoms(AtomDict & tCurAtom,
            std::vector<CrystInfo>::iterator tCryst) {
        for (int i = -myNBDepth; i < myNBDepth + 1; i++) {
            for (int j = -myNBDepth; j < myNBDepth + 1; j++) {
                for (int k = -myNBDepth; k < myNBDepth + 1; k++) {
                    if (!(i == 0 && j == 0 && k == 0)) {
                        AtomDict aAtom(tCurAtom);
                        //std::string label = "_"+ IntToStr(5-i) + IntToStr(5-j) + IntToStr(5-k);
                        //aAtom.sId = label;
                        aAtom.seriNum = (int) refAtoms.size();
                        aAtom.fracCoords.clear();
                        aAtom.fracCoords.push_back(tCurAtom.fracCoords[0] + i);
                        aAtom.fracCoords.push_back(tCurAtom.fracCoords[1] + j);
                        aAtom.fracCoords.push_back(tCurAtom.fracCoords[2] + k);
                        //FractToOrtho(aAtom.fracCoords, aAtom.coords, tCryst->itsCell->a,
                        //             tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                        //             tCryst->itsCell->beta, tCryst->itsCell->gamma);
                        refAtoms.push_back(aAtom);
                        /*
                        std::cout << "Added ref atom: ID-> " << aAtom.id 
                                  << "  symm ID  " << aAtom.sId << std::endl
                                  << "Coords : " << std::endl;
                       
                        for (std::vector<REAL>::iterator iC=aAtom.coords.begin();
                                iC != aAtom.coords.end(); iC++)
                        {
                            std::cout << *iC << std::endl;
                        }
                         */
                    }
                }
            }
        }
    }

    void MolGenerator::swithAtoms(std::vector<CrystInfo>::iterator tCryst) {
        allAtoms.clear();
        for (std::vector<AtomDict>::iterator iAt = refAtoms.begin();
                iAt != refAtoms.end(); iAt++) {
            allAtoms.push_back(*iAt);
        }
        refAtoms.clear();


        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) {
            iAt->coords.clear();
            FractToOrtho(iAt->fracCoords, iAt->coords, tCryst->itsCell->a,
                    tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                    tCryst->itsCell->beta, tCryst->itsCell->gamma);
        }

        //std::cout << "Number of atoms in allAtoms " 
        //          << allAtoms.size() << std::endl;


        int iC = 0;
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) {
            if (iAt->sId == "555") {
                iC++;
            }
        }
        //std::cout << "number of atoms in the center unit cell is "
        //        << iC << std::endl;


    }

    void MolGenerator::cleanUnconnAtoms() {
        std::vector<AtomDict> tAtoms;
        std::map<int, int> tAtMap;
        int j = 0;
        for (int i = 0; i < (int) allAtoms.size(); i++) {
            if (allAtoms[i].connAtoms.size() != 0) {
                tAtoms.push_back(allAtoms[i]);
                tAtMap[i] = j;
                j++;
            }
        }



        //std::cout << "Before clean, number of atoms is " << allAtoms.size() << std::endl;
        //std::cout << "After clean, number of atoms is " << tAtoms.size() << std::endl;

        allAtoms.clear();

        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++) {
            for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                    iCo != iAt->connAtoms.end(); iCo++) {
                (*iCo) = tAtMap[(*iCo)];
            }
            iAt->seriNum = (int) allAtoms.size();
            allAtoms.push_back(*iAt);
        }

        /*
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            std::cout << "Atom " << iAt->id << "(cell idx " << iAt->sId
                      << ") bonds to " 
                      << iAt->connAtoms.size() << " atoms. They are: " << std::endl;
            for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                    iCo !=iAt->connAtoms.end(); iCo++)
            {
                std::cout << allAtoms[*iCo].id << "(cell idx "
                          << allAtoms[*iCo].sId << ")" << std::endl;
            }
                       
        }
         */

    }



    // get atom connections within an unit cell

    void MolGenerator::setUniqueAtomLinks(PeriodicTable & tPTab) {
        REAL covalent_sensitivity = 0.20;

        for (unsigned i = 0; i < allAtoms.size(); i++) {
            for (unsigned j = i + 1; j < allAtoms.size(); j++) {
                // REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                std::vector<REAL> linkRange;
                getBondingRangePairAtoms(allAtoms[i], allAtoms[j],
                        covalent_sensitivity, tPTab,
                        linkRange);
                if (linkRange[0] > 0.20 && linkRange[1] > 0.20) {
                    if (rD > linkRange[0] && rD < linkRange[1]) {
                        setOneUniqueBondCell(i, j, rD);
                        allAtoms[i].connAtoms.push_back(j);
                        allAtoms[j].connAtoms.push_back(i);
                    }
                }
            }
        }

        //std::cout << "Number of atoms in cell " << allAtoms.size() << std::endl;
        //std::cout << "Number of bonds (not unique bonds) : " << bonds.size() << std::endl;

        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) {
            std::cout << "atom " << iAt->seriNum << " is connected to the following atoms :"
                    << std::endl;
            for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                    iNB != iAt->connAtoms.end(); iNB++) {
                std::cout << "atom " << *iNB << std::endl;
            }
        }

    }

    void MolGenerator::setUniqueAtomLinks(PeriodicTable & tPTab,
            std::vector<CrystInfo>::iterator tCryst) {
        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 = 0.25;
        REAL covalent_sensitivity2 = 0.50;
        for (unsigned i = 0; i < allAtoms.size(); i++) {
            std::cout << std::endl << "================================================================" << std::endl;
            std::cout << "Look for bonds to atom " << allAtoms[i].id
                    // << " its sID " << allAtoms[i].sId 
                    << " its serial number " << allAtoms[i].seriNum << std::endl;

            for (unsigned j = i + 1; j < allAtoms.size(); j++) {
                // REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);

                REAL rD = getBondLenFromFracCoords(allAtoms[i].fracCoords, allAtoms[j].fracCoords,
                        tCryst->itsCell->a, tCryst->itsCell->b,
                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                        tCryst->itsCell->beta, tCryst->itsCell->gamma);
                std::vector<REAL> linkRange;
                if ((!allAtoms[i].isMetal) && (!allAtoms[j].isMetal)) {
                    covalent_sensitivity = covalent_sensitivity1;
                } else {
                    covalent_sensitivity = covalent_sensitivity2;
                }
                getBondingRangePairAtoms2(allAtoms[i], allAtoms[j],
                        covalent_sensitivity, tPTab,
                        linkRange);


                if (linkRange[0] > 0.20 && linkRange[1] > 0.20) {

                    if (rD < linkRange[1]) {
                        setOneUniqueBondCell(i, j, rD);
                        allAtoms[i].connAtoms.push_back(j);
                        allAtoms[j].connAtoms.push_back(i);
                    }
                }
            }
        }

        std::cout << "Number of atoms in cell " << allAtoms.size() << std::endl;
        std::cout << "Number of bonds (not unique bonds) : " << bonds.size() << std::endl;

    }

    bool MolGenerator::inBonds(int tIdx1, int tIdx2,
            std::vector<BondDict>& tBonds) {
        //std::cout << "input atom " << tIdx1 << " and " << tIdx2 << std::endl;

        for (std::vector<BondDict>::iterator iBo = tBonds.begin();
                iBo != tBonds.end(); iBo++) {
            //std::cout << "bond atom index " << iBo->atomsIdx[0]
            //          << " and " << iBo->atomsIdx[1] << std::endl;

            if ((iBo->atomsIdx[0] == tIdx1 && iBo->atomsIdx[1] == tIdx2)
                    ||
                    (iBo->atomsIdx[0] == tIdx2 && iBo->atomsIdx[1] == tIdx1)) {
                return true;
            }
        }
        return false;
    }

    void MolGenerator::getUniqueBonds(PeriodicTable & tPTab) {
        /*
        
        REAL covalent_sensitivity;
        
        
        NeighbListDict  tNBListOfSystem;
        
        int aDim  = 3;
        int aMode = 0;
        LIBMOL::REAL tCellL     = 3.5;
        LIBMOL::REAL tCellShell = 0.5;
        
        tNBListOfSystem.building(allAtoms, aDim, tCellL, tCellShell, aMode);
        
        std::cout << "NB list for refAtoms set " << std::endl;   
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->sId=="555")
            {
                std::cout << "atom "<< iAt->id << "(serial number " 
                          << iAt->seriNum  <<") has " << iAt->neighbAtoms.size() 
                          << " NB atoms " << std::endl;
                
            }
        }
        
        exit(1);
        
        
        // std::vector<std::string>   existBondID;
        int j=0;
        for (unsigned i=0; i <refAtoms.size(); i++)
        {
            if (refAtoms[i].sId=="")
            {
                j++;
                //std::cout << "Look for bonds to atom " << refAtoms[i].id 
                //           << " its sID " << refAtoms[i].sId << std::endl;
                //std::cout << "Its has " << refAtoms[i].neighbAtoms.size() 
                //          << " neighbor atoms. " <<std::endl;
                for (std::vector<int>::iterator iNB=refAtoms[i].neighbAtoms.begin();
                        iNB !=refAtoms[i].neighbAtoms.end(); iNB++)
                {
                    REAL rD = distanceV(refAtoms[i].coords, refAtoms[(*iNB)].coords);
                    std::vector<REAL> bondRange;
                    getBondingRangePairAtoms(refAtoms[i], refAtoms[(*iNB)],
                                         covalent_sensitivity, tPTab, 
                                         bondRange);
                    
                    //std::cout << "NB Atom " << refAtoms[(*iNB)].id << " its sID " 
                    //          << refAtoms[(*iNB)].sId << std::endl;
                    //std::cout << "Distance " << rD << std::endl;
                    //std::cout << "Range between " << bondRange[0]
                    //          << " and " << bondRange[1] << std::endl;
                              
                
                    if (bondRange[0] >0.20 && bondRange[1] >0.20)
                    {
                        if (rD > bondRange[0] && rD < bondRange[1])
                        {
                            // setOneUniqueBondCrys(i, *iNB, rD);
                            if (std::find(refAtoms[i].connAtoms.begin(), refAtoms[i].connAtoms.end(), *iNB)
                                    ==refAtoms[i].connAtoms.end())
                            {
                                refAtoms[i].connAtoms.push_back(*iNB);
                            }
                        }
                    }
                }
            }
        }
        // std::cout << "Number of atoms in the unit cell considered " << j << std::endl;
        
        setUniqueAtomLinks(tPTab);
        
         */
    }

    void MolGenerator::getUniqueAtomLinks(PeriodicTable & tPTab,
                         std::vector<CrystInfo>::iterator tCryst) 
    {

        /*  
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            std::cout << "Coordinates for atom " << iAt->id << " (cell "
                      << iAt->sId << ")" << std::endl;
                
            for (std::vector<REAL>::iterator iX=iAt->coords.begin();
                    iX !=iAt->coords.end(); iX++)
            {
                std::cout << *iX << std::endl;
            }
        }
         */

        REAL covalent_sensitivity = 0.15;
        REAL covalent_sensitivity1 = 0.22;
        REAL covalent_sensitivity2 = 0.25;
        REAL covalent_sensitivity3 = 0.30;
        REAL covalent_sensitivity4 = 0.60;

        NeighbListDict tNBListOfSystem;

        int aDim = 3;
        int aMode = 0;
        LIBMOL::REAL tCellL = 3.5;
        LIBMOL::REAL tCellShell = 0.5;

        tNBListOfSystem.building(allAtoms, aDim, tCellL, tCellShell, aMode);

        /*
        std::cout << "NB list for allAtoms set " << std::endl;   


        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                std::cout << "atom "<< iAt->id << "(serial number " 
                          << iAt->seriNum  <<") has " << iAt->neighbAtoms.size() 
                          << " NB atoms " << std::endl;
                
                
                    for (std::vector<int>::iterator iNB=iAt->neighbAtoms.begin();
                            iNB != iAt->neighbAtoms.end(); iNB++)
                    {
                        std::cout << "NB atom " << allAtoms[*iNB].id << std::endl;
                    }
                
            }
        }
        */

        // std::vector<std::string>   existBondID;
        // int j=0;
        for (unsigned i = 0; i < allAtoms.size(); i++) 
        {
            //if (allAtoms[i].sId=="555")
            //{
            //j++;
            std::cout << "Look for bonds to atom " << allAtoms[i].id
                      << "(serial number  " << allAtoms[i].seriNum
                      << ") " << std::endl;
            
            for (std::vector<int>::iterator iNB = allAtoms[i].neighbAtoms.begin();
                    iNB != allAtoms[i].neighbAtoms.end(); iNB++) 
            {

                REAL rD = getBondLenFromFracCoords(allAtoms[i].fracCoords, allAtoms[(*iNB)].fracCoords,
                        tCryst->itsCell->a, tCryst->itsCell->b,
                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                        tCryst->itsCell->beta, tCryst->itsCell->gamma);

                std::vector<REAL> bondRange;

                if ((allAtoms[i].chemType == "H" && allAtoms[*iNB].chemType == "O")
                        || (allAtoms[i].chemType == "O" && allAtoms[*iNB].chemType == "H")) {

                    covalent_sensitivity = covalent_sensitivity4;

                } else if (allAtoms[*iNB].chemType == "O" || allAtoms[i].chemType == "O") {

                    covalent_sensitivity = covalent_sensitivity3;
                } else if (allAtoms[*iNB].chemType == "H"
                        || allAtoms[i].chemType == "H") {
                    covalent_sensitivity = covalent_sensitivity3;
                } else if ((!allAtoms[i].isMetal) && (!allAtoms[*iNB].isMetal)) {
                    covalent_sensitivity = covalent_sensitivity1;
                } else {
                    covalent_sensitivity = covalent_sensitivity2;
                }

                //std::cout << "covalent_sensitivity=" << covalent_sensitivity << std::endl; 

                getBondingRangePairAtoms2(allAtoms[i], allAtoms[(*iNB)],
                        covalent_sensitivity, tPTab,
                        bondRange);
                
                std::cout << "Distance between: " << std::endl
                                      << "Atom 1 " << allAtoms[i].id 
                                      << " of serial number " 
                                      << allAtoms[i].seriNum
                                      << " from original atom " 
                                      << allAtoms[allAtoms[i].fromOrig].id 
                                      << std::endl
                                      << "Atom 2 " << allAtoms[*iNB].id 
                                      << " of serial number " << allAtoms[*iNB].seriNum
                                      << " from original atom " 
                                      << allAtoms[allAtoms[*iNB].fromOrig].id 
                                      << " is " << rD << std::endl;
                std::cout << "Range between " << bondRange[0]
                                      << " and " << bondRange[1] << std::endl;
                
                
                if (bondRange[0] > 0.20 && bondRange[1] > 0.20) 
                {
                    if (rD > bondRange[0] && rD < bondRange[1]) 
                    {
                        
                        // setOneUniqueBondCell(i, *iNB, rD);
                        if (std::find(allAtoms[i].connAtoms.begin(), allAtoms[i].connAtoms.end(), *iNB)
                                == allAtoms[i].connAtoms.end()) {
                            allAtoms[i].connAtoms.push_back(*iNB);
                        }
                        if (std::find(allAtoms[*iNB].connAtoms.begin(), 
                                      allAtoms[*iNB].connAtoms.end(), i)
                                      == allAtoms[*iNB].connAtoms.end()) 
                        {
                            allAtoms[*iNB].connAtoms.push_back(i);
                        }
                       
                        if (allAtoms[i].isMetal || allAtoms[*iNB].isMetal)
                        {
                            
                            //std::cout << "Its has " << allAtoms[i].neighbAtoms.size()
                            //          << " neighbor atoms. " << std::endl;
                           /*
                            std::cout << "Distance between: " << std::endl
                                      << "Atom 1 " << allAtoms[i].id 
                                      << " of serial number " 
                                      << allAtoms[i].seriNum
                                      << " from original atom " 
                                      << allAtoms[allAtoms[i].fromOrig].id 
                                      << std::endl
                                      << "Atom 2 " << allAtoms[*iNB].id 
                                      << " of serial number " << allAtoms[*iNB].seriNum
                                      << " from original atom " 
                                      << allAtoms[allAtoms[*iNB].fromOrig].id 
                                      << " is " << rD << std::endl;
                            std::cout << "Range between " << bondRange[0]
                                      << " and " << bondRange[1] << std::endl;
                            
                            */
                            std::cout << "a bond between " 
                                      << allAtoms[i].id << " and "
                                      << allAtoms[*iNB].id 
                                      << " is added to the bond_list_cell "
                                      << std::endl << "Its bond length is " << rD
                                      << std::endl;
                        }
                    }
                }
            }
        }

        // std::cout << "Number of atoms in the unit cell considered " << j << std::endl;
    
        // cleanUnconnAtoms();
        /*
        std::cout << "*********************************" << std::endl;
        std::cout << "Before check " << std::endl
                  << "*********************************" << std::endl;
                
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) {
            if (iAt->isInPreCell) {
                std::cout << "Atom " << iAt->id << " (serial number " << iAt->seriNum
                        << " ) bonds to " << iAt->connAtoms.size() << " atoms"
                        << std::endl;
                // std::cout << "The NB atoms are :" << std::endl;
                for (std::vector<int>::iterator iC = iAt->connAtoms.begin();
                        iC != iAt->connAtoms.end(); iC++) {
                    std::cout << "A 1st NB atoms is :" << std::endl;
                    std::cout << "Atom " << allAtoms[*iC].id
                            << " of serial number " << allAtoms[*iC].seriNum
                            << std::endl;
                    //std::cout << "it connects to "
                    //        << allAtoms[allAtoms[*iC].fromOrig].connAtoms.size()
                    //        << " 2nd NB atoms " << std::endl << std::endl;
                }
            }
        }
        */
        checkAtomLinks();

        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) {
            if (iAt->isInPreCell) {
                std::cout << "Atom " << iAt->id << " (serial number " << iAt->seriNum
                        << " ) bonds to " << iAt->connAtoms.size() << " atoms"
                        << std::endl;
                // std::cout << "The NB atoms are :" << std::endl;
                for (std::vector<int>::iterator iC = iAt->connAtoms.begin();
                        iC != iAt->connAtoms.end(); iC++) {
                    std::cout << "A 1st NB atoms is :" << std::endl;
                    std::cout << "Atom " << allAtoms[*iC].id
                            << " of serial number " << allAtoms[*iC].seriNum
                            << std::endl;
                    std::cout << "it connects to "
                            << allAtoms[allAtoms[*iC].fromOrig].connAtoms.size()
                            << " 2nd NB atoms " << std::endl << std::endl;
                }
            }
        }
        
        // setUniqueAtomLinks(tPTab, tCryst);
    }
    
    void MolGenerator::checkAtomLinks()
    {
        if (allAtoms.size() >0)
        {
            for (std::vector<AtomDict>::iterator iAt =allAtoms.begin();
                    iAt !=allAtoms.end(); iAt++)
            {
                if (iAt->isInPreCell)
                {
                    std::vector<int> confirmedLinks;
                    // std::vector<int> falseLinks;
                    
                    std::cout << "************************************" << std::endl;
                    std::cout << "Check atom " << iAt->id << " now " << std::endl;
                    std::cout << "************************************" << std::endl;
                    
                    for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                            iCo != iAt->connAtoms.end(); iCo++)
                    {   
                        int tOrigAt = allAtoms[*iCo].fromOrig;
                        //std::cout << "One connected atom is " << allAtoms[*iCo].id
                        //          << " of serial number " << allAtoms[*iCo].seriNum 
                        //          << " which is from atom " <<  allAtoms[tOrigAt].id
                        //          << " of serial number " 
                        //          << allAtoms[tOrigAt].seriNum <<  std::endl;
                       
                        std::vector<std::string> tAtConns;
                        //std::cout << "The assym cell atom " << allAtoms[tOrigAt].id
                        //          << " connected : " << std::endl;
                        for (std::vector<int>::iterator 
                             iNA=allAtoms[tOrigAt].connAtoms.begin();
                             iNA!=allAtoms[tOrigAt].connAtoms.end(); iNA++)
                        {
                            //std::cout << "atom " << allAtoms[*iNA].id << std::endl;
                            std::string tID = TrimSpaces(allAtoms[*iNA].id);
                            tAtConns.push_back(tID);
                        }
                        
                        if (std::find(tAtConns.begin(), tAtConns.end(), iAt->id)
                                  !=tAtConns.end())
                        {
                            std::cout << iAt->id << " is confirmed bonding to " 
                                      << "atom " <<  allAtoms[*iCo].id << std::endl;
                            confirmedLinks.push_back(*iCo);
                        }
                        
                        /*
                        else
                        {
                            falseLinks.push_back(*iCo);
                        }
                         */
                    }
                    
                    /*
                    if (falseLinks.size()==0)
                    {
                        std::cout << "All bonds to atom " << iAt->id 
                                  << " is confirmed " << std::endl;
                    }
                    else
                    {
                         std::cout << falseLinks.size() << " atoms "
                                   << " are linked to atom " 
                                   << iAt->id << " by wrong symm ops " << std::endl;
                        std::cout << "They are atoms : " << std::endl;
                        for (std::vector<int>::iterator iNB=falseLinks.begin();
                                iNB != falseLinks.end(); iNB++)
                        {
                            std::cout << allAtoms[*iNB].id << std::endl;
                        }         
                    } 
                    */
                    
                    if (confirmedLinks.size() !=0)
                    {
                        iAt->connAtoms.clear();
                        for (std::vector<int>::iterator iNC=confirmedLinks.begin();
                                iNC !=confirmedLinks.end(); iNC++)
                        {
                            iAt->connAtoms.push_back(*iNC);
                        }
                    }
                }
            }
        }
    }

    void MolGenerator::getUniqueBondsMols(Molecule& tMol,
            std::vector<CrystInfo>::iterator tCryst) {
        std::map<std::string, int> aMA;
        std::vector<REAL> aV;

        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            if (iAt->isInPreCell) {
                for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                        iCo != iAt->connAtoms.end(); iCo++) {
                    std::string str1 = iAt->id + "_" + tMol.atoms[(*iCo)].id;
                    std::string str2 = tMol.atoms[(*iCo)].id + "_" + iAt->id;
                    if ((*iCo) > iAt->seriNum &&
                            (aMA.find(str1) == aMA.end() && aMA.find(str2) == aMA.end())) {
                        aMA[str1] = 1;
                        REAL rD = getBondLenFromFracCoords(iAt->fracCoords, tMol.atoms[(*iCo)].fracCoords,
                                tCryst->itsCell->a, tCryst->itsCell->b,
                                tCryst->itsCell->c, tCryst->itsCell->alpha,
                                tCryst->itsCell->beta, tCryst->itsCell->gamma);
                        BondDict aBond;
                        aBond.atomsIdx.push_back(iAt->seriNum);
                        aBond.atomsIdx.push_back(*iCo);
                        aBond.value = rD;
                        tMol.bonds.push_back(aBond);
                        aV.push_back(aBond.value);
                        //std::cout << "a bond between " << iAt->id << " and "
                        //          << tMol.atoms[(*iCo)].id << " is added to the bond_list_cell " 
                        //          << std::endl << "Its bond length is " << aBond.value 
                        //          << std::endl;
                    } else {
                        REAL rD = getBondLenFromFracCoords(iAt->fracCoords, tMol.atoms[(*iCo)].fracCoords,
                                tCryst->itsCell->a, tCryst->itsCell->b,
                                tCryst->itsCell->c, tCryst->itsCell->alpha,
                                tCryst->itsCell->beta, tCryst->itsCell->gamma);
                        if (!inVectABS(aV, rD, 0.00001)) {
                            BondDict aBond;
                            aBond.atomsIdx.push_back(iAt->seriNum);
                            aBond.atomsIdx.push_back(*iCo);
                            aBond.value = rD;
                            tMol.bonds.push_back(aBond);
                            aV.push_back(aBond.value);
                            //std::cout << "a bond between " << iAt->id << " and "
                            //          << tMol.atoms[(*iCo)].id << " is added to the bond_list_cell " 
                            //          << std::endl << "Its bond length is(unique) " << aBond.value 
                            //          << std::endl;
                        }

                    }
                }
            }
        }
    }

    void MolGenerator::getUniqueBondsMols2(Molecule& tMol) {

        std::vector<REAL> aV;

        // std::cout << "number of pre-bonds " << tMol.allBonds.size() << std::endl;

        for (std::vector<BondDict>::iterator iB = tMol.allBonds.begin();
                iB != tMol.allBonds.end(); iB++) {

            if (aV.size() == 0) {
                tMol.bonds.push_back(*iB);
            } else if (!inVectABS(aV, iB->value, 0.00000001)) {
                tMol.bonds.push_back(*iB);
            }
            /* symm generated bonds should be the same
            else 
            {
                std::cout << "bond length duplicated "
                          << iB->value << std::endl
                          << "Between atoms "
                          << tMol.atoms[iB->atomsIdx[0]].id << " and "
                          << tMol.atoms[iB->atomsIdx[1]].id << std::endl;
            }
             */
            aV.push_back(iB->value);
        }

        //std::cout << "Number of bonds in the mol " << tMol.bonds.size() << std::endl;

    }

    void MolGenerator::getAllBondsInOneMol(Molecule & tMol,
            std::vector<CrystInfo>::iterator tCryst) {
        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                    iCo != iAt->connAtoms.end(); iCo++) {
                if ((*iCo) > iAt->seriNum) {
                    REAL rD = getBondLenFromFracCoords(iAt->fracCoords, tMol.atoms[(*iCo)].fracCoords,
                            tCryst->itsCell->a, tCryst->itsCell->b,
                            tCryst->itsCell->c, tCryst->itsCell->alpha,
                            tCryst->itsCell->beta, tCryst->itsCell->gamma);

                    BondDict aBond;
                    aBond.seriNum = tMol.allBonds.size();
                    aBond.atomsIdx.push_back(iAt->seriNum);
                    aBond.atomsIdx.push_back(*iCo);
                    aBond.atoms.push_back(iAt->id);
                    aBond.atoms.push_back(tMol.atoms[*iCo].id);
                    aBond.value = rD;
                    aBond.molIdx = tMol.seriNum;
                    aBond.isInSameRing = checkIfBondInSameRing(tMol.atoms, iAt->seriNum, *iCo);
                    tMol.allBonds.push_back(aBond);

                    /*
                    std::cout << "a bond between " << iAt->id << " and "
                              << tMol.atoms[(*iCo)].id << " is added to the bond_list_cell " 
                              << std::endl << "Its bond length is " << aBond.value 
                              << std::endl;
                    if (aBond.isInSameRing)
                    {
                        std::cout << "Two atoms in the bond are in the same ring " << std::endl;
                    }
                    else
                    {
                        std::cout << "Two atoms in the bond are not in the same ring " << std::endl;
                    }
                     */
                }
            }
        }

        // std::cout << "Number of bonds is " << tMol.allBonds.size() << std::endl;
        // Modify and setup atom's sp features


        modAtomsBondingAndChiralCenter(tMol.atoms, tMol.allBonds, tMol.angles, tMol.rings, 1);

        setAtomsNB1NB2_SP(tMol.atoms);

        // Tests on the extra-electrons scheme
        setAtomNFormalCharge(tMol);
        setAllAtomEXcessElectrons(tMol.atoms);
        setAtomsNB1NB2_exElectrons(tMol.atoms);


        setBondsAndAngles_NB1NB2_SP(tMol.atoms, tMol.allBonds, tMol.angles);
        setBondsAndAngles_NB1NB2_EE(tMol.atoms, tMol.allBonds, tMol.angles);

        getUniqueBondsMols2(tMol);

    }

    void MolGenerator::getBondingRangePairAtoms(AtomDict & tAtm1,
            AtomDict & tAtm2,
            REAL tSens,
            PeriodicTable & tPTab,
            std::vector<REAL> & tRange) {
        tRange.clear();
        tRange.push_back(-1.0);
        tRange.push_back(-1.0);


        if (tPTab.elemProps.find(tAtm1.chemType) != tPTab.elemProps.end()
                && tPTab.elemProps.find(tAtm2.chemType) != tPTab.elemProps.end()) {
            if (tPTab.elemProps[tAtm1.chemType].find("cova")
                    != tPTab.elemProps[tAtm1.chemType].end()
                    && tPTab.elemProps[tAtm2.chemType].find("cova")
                    != tPTab.elemProps[tAtm2.chemType].end()) {
                if (tPTab.elemProps[tAtm1.chemType]["cova"] > 0.2
                        && tPTab.elemProps[tAtm2.chemType]["cova"] > 0.2) {
                    REAL tRad, tExtraD;
                    tRad = tPTab.elemProps[tAtm1.chemType]["cova"]
                            + tPTab.elemProps[tAtm2.chemType]["cova"];
                    tExtraD = tSens*tRad;
                    tRange[0] = tRad - tExtraD;
                    if (tRange[0] < 0.2) {
                        tRange[0] = 0.2;
                    }
                    tRange[1] = tRad + tExtraD;
                } else {
                    std::cout << "Bug in at least one of atom covalent radius! "
                            << std::endl << "atom " << tAtm1.id << " coval: "
                            << tPTab.elemProps[tAtm1.chemType]["cova"] << std::endl
                            << "atom " << tAtm2.id << " coval: "
                            << tPTab.elemProps[tAtm2.chemType]["cova"] << std::endl;
                }
            } else {
                std::cout << "Bug! At least one of covalent radius for "
                        << tAtm1.id << " or " << tAtm2.id
                        << " is not defined the internal periodic table"
                        << std::endl;
            }
        } else {
            std::cout << "Bug! At least one element im "
                    << tAtm1.id << " or " << tAtm2.id
                    << " is not defined the internal periodic table "
                    << std::endl;
        }
    }

    void MolGenerator::getBondingRangePairAtoms2(AtomDict & tAtm1,
            AtomDict & tAtm2,
            REAL tSens,
            PeriodicTable & tPTab,
            std::vector<REAL> & tRange) 
    {
        tRange.clear();
        tRange.push_back(-1.0);
        tRange.push_back(-1.0);

        REAL comp1 = 0.0, comp2 = 0.0;

        if (tPTab.elemProps.find(tAtm1.chemType) != tPTab.elemProps.end()
                && tPTab.elemProps.find(tAtm2.chemType) != tPTab.elemProps.end()) {
            if (tAtm1.formalCharge != 0 && tAtm2.formalCharge != 0) {

                if (tAtm1.formalCharge > 0 && tAtm2.formalCharge < 0) {
                    if (tPTab.elemProps[tAtm1.chemType].find("ionM+")
                            != tPTab.elemProps[tAtm1.chemType].end()
                            && tPTab.elemProps[tAtm2.chemType].find("ionM-")
                            != tPTab.elemProps[tAtm2.chemType].end()) {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["ionM+"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["ionM-"];
                    } else {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                    }
                } else if (tAtm1.formalCharge < 0 && tAtm2.formalCharge > 0) {
                    if (tPTab.elemProps[tAtm1.chemType].find("ionM-")
                            != tPTab.elemProps[tAtm1.chemType].end()
                            && tPTab.elemProps[tAtm2.chemType].find("ionM+")
                            != tPTab.elemProps[tAtm2.chemType].end()) {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["ionM-"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["ionM+"];
                    } else {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                    }
                } else {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                }
            }
            else if (tAtm1.isMetal && tAtm2.isMetal) {
                if (tPTab.elemProps[tAtm1.chemType].find("ionM+")
                        != tPTab.elemProps[tAtm1.chemType].end()
                        && tPTab.elemProps[tAtm2.chemType].find("ionM-")
                        != tPTab.elemProps[tAtm2.chemType].end()) {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["ionM+"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["ionM-"];
                }
                else if (tPTab.elemProps[tAtm2.chemType].find("ionM+")
                        != tPTab.elemProps[tAtm2.chemType].end()
                        && tPTab.elemProps[tAtm1.chemType].find("ionM-")
                        != tPTab.elemProps[tAtm1.chemType].end()) {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["ionM-"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["ionM+"];
                }
                else {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                }
            }
            else {
                comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                //std::cout << "atom 1 " << tAtm1.chemType << std::endl;
                //std::cout << "inside comp1 " << comp1 << std::endl;
                comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                //std::cout << "atom 2 " << tAtm2.chemType << std::endl;
                //std::cout << "comp2 " << comp2 << std::endl;
            }

            REAL tRad, tExtraD;
            tRad = comp1 + comp2;
            tExtraD = tSens*tRad;
            if (tAtm2.chemType == "H" || tAtm1.chemType == "H" 
                || tAtm1.chemType=="D" || tAtm2.chemType =="D") {
                tRange[0] = tRad - 0.5 * tExtraD;
            } else {
                tRange[0] = tRad - 1.55 * tExtraD;
            }
            if (tRange[0] < 0.2) {
                tRange[0] = 0.2;
            }
            tRange[1] = tRad + 0.8 * tExtraD;

            //std::cout <<  "comp1 " << comp1 << std::endl;
            //std::cout <<  "comp2 " << comp2 << std::endl;
            //std::cout << "tRad " << tRad << std::endl;
            //std::cout << "tSen " << tSens << std::endl;
            //std::cout << "tExtraD " << std::endl;

        } else {
            std::cout << "Bug! At least one of covalent radius for "
                    << tAtm1.id << " or " << tAtm2.id
                    << " is not defined the internal periodic table"
                    << std::endl;
        }
    }

    void MolGenerator::setOneUniqueBondCrys(int tIdxAtm1, int tIdxAtm2,
            REAL rD) {

        ID id1 = refAtoms[tIdxAtm1].id, id2 = refAtoms[tIdxAtm2].id;
        bool tFound = false;

        for (std::vector<BondDict>::iterator iB = bonds.begin();
                iB != bonds.end(); iB++) {
            ID id3 = refAtoms[iB->atomsIdx[0]].id,
                    id4 = refAtoms[iB->atomsIdx[1]].id;
            if ((id1.compare(id3) == 0 && id2.compare(id4) == 0)
                    || (id1.compare(id4) == 0 && id2.compare(id3) == 0)) {
                if (fabs(iB->value - rD) <= 0.00001) {
                    tFound = true;
                    break;
                }
            }
        }

        if (!tFound) {
            BondDict aBond;
            aBond.atomsIdx.push_back(tIdxAtm1);
            aBond.atomsIdx.push_back(tIdxAtm2);
            aBond.value = rD;
            bonds.push_back(aBond);
            //std::cout << "a bond between " << refAtoms[tIdxAtm1].id << " and "
            //          << refAtoms[tIdxAtm2].id << " is added to the b_list " 
            //          << std::endl << "Its bond length is " << aBond.value 
            //          << std::endl;
        }
    }

    void MolGenerator::setOneUniqueBondCell(int tIdxAtm1, int tIdxAtm2,
            REAL rD) {
        /*
        ID id1=allAtoms[tIdxAtm1].id, id2 = allAtoms[tIdxAtm2].id;
        bool tFound=false;
        
        for (std::vector<BondDict>::iterator iB=bonds.begin();
                iB !=bonds.end(); iB++)
        {
            ID id3=allAtoms[iB->atomsIdx[0]].id, 
               id4=allAtoms[iB->atomsIdx[1]].id;
            if ((id1.compare(id3)==0 && id2.compare(id4)==0)
                 || (id1.compare(id4)==0 && id2.compare(id3)==0))
            {
                if (fabs(iB->value-rD) <=0.00001)
                {
                    tFound=true;
                    break;
                }
            }
        }
         */
        //if (!tFound)
        //{

        BondDict aBond;
        aBond.atomsIdx.push_back(tIdxAtm1);
        aBond.atomsIdx.push_back(tIdxAtm2);
        aBond.value = rD;
        bonds.push_back(aBond);


        std::cout << "a bond between " << allAtoms[tIdxAtm1].id << " and "
                << allAtoms[tIdxAtm2].id << " is added to the bond_list_cell "
                << std::endl << "Its bond length is " << aBond.value
                << std::endl;

        //}
    }

    REAL MolGenerator::getBondLenFromFracCoords(std::vector<REAL>& tCoord1,
            std::vector<REAL>& tCoord2,
            REAL a, REAL b, REAL c,
            REAL alpha, REAL beta, REAL gamma) {
        if (tCoord1.size() == tCoord2.size() && tCoord1.size() == 3) {
            std::vector<REAL> deltX;
            for (unsigned i = 0; i < 3; i++) {
                deltX.push_back((tCoord2[i] - tCoord1[i]));
            }
            
            return sqrt(a * a * deltX[0] * deltX[0] + b * b * deltX[1] * deltX[1] 
                        + c * c * deltX[2] * deltX[2]
                        + 2 * b * c * deltX[1] * deltX[2] * cos(alpha * PI180)
                        + 2 * c * a * deltX[2] * deltX[0] * cos(beta * PI180)
                        + 2 * a * b * deltX[0] * deltX[1] * cos(gamma * PI180));
            
        } 
        else 
        {
            std::cout << "Error: dim of atom 1 coords " << tCoord1.size() << std::endl;
            std::cout << "Error: dim of atom 2 coords " << tCoord2.size() << std::endl;
            return -1.0;
        }
    }

    void MolGenerator::getUniqAngles() 
    {
        for (int i = 0; i < (int) allAtoms.size(); i++) {
            for (int j = 0; j < (int) allAtoms[i].connAtoms.size(); j++) {
                for (int k = j + 1; k < (int) allAtoms[i].connAtoms.size(); k++) {
                    int i1 = allAtoms[i].connAtoms[j];
                    int i2 = allAtoms[i].connAtoms[k];
                    //std::cout << "Angle between " << allAtoms[i].id 
                    //          << "(center) and " << allAtoms[i1].id
                    //          << " and " << allAtoms[i2].id << std::endl;

                    AngleDict aAng;
                    aAng.anchorID = allAtoms[i].id;
                    aAng.anchorPos = i;
                    aAng.atoms.push_back(i);

                    if ((int) allAtoms[i1].connAtoms.size() >=
                            (int) allAtoms[i1].connAtoms.size()) {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                    } else {
                        aAng.atoms.push_back(i2);
                        aAng.atoms.push_back(i1);
                    }
                    std::vector<REAL> tV1, tV2;
                    for (int iC = 0; iC < (int) allAtoms[i].coords.size(); iC++) {
                        tV1.push_back(allAtoms[i1].coords[iC] - allAtoms[i].coords[iC]);
                        tV2.push_back(allAtoms[i2].coords[iC] - allAtoms[i].coords[iC]);
                    }
                    aAng.value = getAngle2V(tV1, tV2);
                    aAng.sigValue = 3.0;
                    aAng.numCodValues = 0;
                    angles.push_back(aAng);
                }
            }
        }
        /*
        for (std::vector<AngleDict>::iterator iA=angles.begin();
                iA !=angles.end(); iA++)
        {
            std::cout << "angle between " << allAtoms[iA->atoms[0]].id 
                      << "(center) and " <<  allAtoms[iA->atoms[1]].id
                      << " and " << allAtoms[iA->atoms[2]].id 
                      << " is " << iA->value*PID180 << std::endl;
        }
         */
    }

    void MolGenerator::getUniqAngles(std::vector<CrystInfo>::iterator tCryst) {
        for (int i = 0; i < (int) allAtoms.size(); i++) {
            for (int j = 0; j < (int) allAtoms[i].connAtoms.size(); j++) {
                for (int k = j + 1; k < (int) allAtoms[i].connAtoms.size(); k++) {
                    int i1 = allAtoms[i].connAtoms[j];
                    int i2 = allAtoms[i].connAtoms[k];
                    //std::cout << "Angle between " << allAtoms[i].id 
                    //          << "(center) and " << allAtoms[i1].id
                    //          << " and " << allAtoms[i2].id << std::endl;

                    AngleDict aAng;
                    aAng.anchorID = allAtoms[i].id;
                    aAng.anchorPos = i;
                    aAng.atoms.push_back(i);

                    if ((int) allAtoms[i1].connAtoms.size() >=
                            (int) allAtoms[i2].connAtoms.size()) {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                    } else {
                        aAng.atoms.push_back(i2);
                        aAng.atoms.push_back(i1);
                    }
                    /*
                    std::vector<REAL> tV1, tV2;    
                    for (int iC=0; iC < (int)allAtoms[i].coords.size(); iC++)
                    {
                        tV1.push_back(allAtoms[i1].coords[iC]-allAtoms[i].coords[iC]);
                        tV2.push_back(allAtoms[i2].coords[iC]-allAtoms[i].coords[iC]);
                    }
                     */
                    // aAng.value        = getAngle2V(tV1, tV2);
                    aAng.value = getAngleValueFromFracCoords(allAtoms[i],
                            allAtoms[i1], allAtoms[i2],
                            tCryst->itsCell->a, tCryst->itsCell->b,
                            tCryst->itsCell->c, tCryst->itsCell->alpha,
                            tCryst->itsCell->beta, tCryst->itsCell->gamma);

                    aAng.sigValue = 3.0;
                    aAng.numCodValues = 0;
                    angles.push_back(aAng);
                }
            }
        }

        /*
        for (std::vector<AngleDict>::iterator iA=angles.begin();
                iA !=angles.end(); iA++)
        {
            std::cout << "angle between " << allAtoms[iA->atoms[0]].id 
                      << "(center) and " <<  allAtoms[iA->atoms[1]].id
                      << " and " << allAtoms[iA->atoms[2]].id 
                      << " is " << iA->value*PID180 << std::endl;
        }
         */
    }

    void MolGenerator::getUniqAngleMols(Molecule & tMol,
            std::vector<CrystInfo>::iterator tCryst) {
        std::map<std::string, int> aMM;
        std::vector<REAL> aV;

        for (int i = 0; i < (int) tMol.atoms.size(); i++) {
            if (tMol.atoms[i].isInPreCell) {
                for (int j = 0; j < (int) tMol.atoms[i].connAtoms.size(); j++) {
                    for (int k = j + 1; k < (int) tMol.atoms[i].connAtoms.size(); k++) {

                        int i1 = tMol.atoms[i].connAtoms[j];
                        int i2 = tMol.atoms[i].connAtoms[k];
                        if (inBonds(i, i1, tMol.bonds) || inBonds(i, i2, tMol.bonds)) {
                            //std::cout << "Angle between atom " << tMol.atoms[i].id 
                            //          << "(" << tMol.atoms[i].seriNum << ") "
                            //          << "(center) coords " << std::endl;
                            //for (std::vector<REAL>::iterator iX=tMol.atoms[i].fracCoords.begin();
                            //        iX!=tMol.atoms[i].fracCoords.end(); iX++)
                            //{
                            //    std::cout << *iX << std::endl;
                            //}


                            //std::cout  << " atom " << tMol.atoms[i1].id
                            //          << "(" << tMol.atoms[i1].seriNum << ") " << std::endl;
                            //for (std::vector<REAL>::iterator iX=tMol.atoms[i1].fracCoords.begin();
                            //        iX!=tMol.atoms[i1].fracCoords.end(); iX++)
                            //{
                            //    std::cout << *iX << std::endl;
                            //}
                            //std::cout << " and atom " << tMol.atoms[i2].id 
                            //          << "(" << tMol.atoms[i2].seriNum << ") "
                            //          << std::endl;
                            //for (std::vector<REAL>::iterator iX=tMol.atoms[i2].fracCoords.begin();
                            //        iX!=tMol.atoms[i2].fracCoords.end(); iX++)
                            //{
                            //    std::cout << *iX << std::endl;
                            //}

                            std::string str1, str2;
                            str1 = tMol.atoms[i1].id + "_" + tMol.atoms[i2].id;
                            str2 = tMol.atoms[i2].id + "_" + tMol.atoms[i1].id;
                            if (aMM.find(str1) == aMM.end() && aMM.find(str2) == aMM.end()) {
                                aMM[str1] = 1;
                                AngleDict aAng;
                                aAng.anchorID = tMol.atoms[i].id;
                                aAng.anchorPos = i;
                                aAng.atoms.push_back(i); // this is for individual molecule output
                                aAng.atomChemTypes.push_back(tMol.atoms[i].chemType);
                                aAng.atomsCodClasses.push_back(tMol.atoms[i].codClass); // this for the overall file

                                if ((int) tMol.atoms[i1].connAtoms.size() >=
                                        (int) tMol.atoms[i2].connAtoms.size()) {
                                    aAng.atoms.push_back(i1);
                                    aAng.atoms.push_back(i2);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                } else {
                                    aAng.atoms.push_back(i2);
                                    aAng.atoms.push_back(i1);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                }

                                aAng.value = getAngleValueFromFracCoords(tMol.atoms[i],
                                        tMol.atoms[i1], tMol.atoms[i2],
                                        tCryst->itsCell->a, tCryst->itsCell->b,
                                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                                        tCryst->itsCell->beta, tCryst->itsCell->gamma);
                                // std::cout << "new key " << aAng.value << std::endl;
                                aAng.sigValue = 3.0;
                                aAng.numCodValues = 0;
                                tMol.angles.push_back(aAng);
                                aV.push_back(aAng.value);
                            } else {
                                REAL rD = getAngleValueFromFracCoords(tMol.atoms[i],
                                        tMol.atoms[i1], tMol.atoms[i2],
                                        tCryst->itsCell->a, tCryst->itsCell->b,
                                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                                        tCryst->itsCell->beta, tCryst->itsCell->gamma);
                                // std::cout << "rD " << rD << std::endl;
                                if (!inVectABS(aV, rD, 0.0001)) {
                                    aMM[str1] = 1;
                                    AngleDict aAng;
                                    aAng.anchorID = tMol.atoms[i].id;
                                    aAng.anchorPos = i;
                                    aAng.atoms.push_back(i); // this is for individual molecule output
                                    aAng.atomChemTypes.push_back(tMol.atoms[i].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i].codClass); // this for the overall file

                                    if ((int) tMol.atoms[i1].connAtoms.size() >=
                                            (int) tMol.atoms[i2].connAtoms.size()) {
                                        aAng.atoms.push_back(i1);
                                        aAng.atoms.push_back(i2);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                    } else {
                                        aAng.atoms.push_back(i2);
                                        aAng.atoms.push_back(i1);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                    }
                                    aAng.value = rD;
                                    // std::cout << "new value " << aAng.value << std::endl;
                                    aAng.sigValue = 3.0;
                                    aAng.numCodValues = 0;
                                    tMol.angles.push_back(aAng);
                                    aV.push_back(rD);
                                }
                            }
                        }
                    }
                }
            }
        }

        /*
        for (std::vector<RingDict>::iterator iR=tMol.rings.begin();
                iR !=tMol.rings.end(); iR++)
        {
            std::cout << "For ring " << iR->rep << std::endl;
            std::cout << "Its representation  " << iR->sRep << std::endl;
            std::cout << "Its atoms are : " << std::endl;
            for (std::vector<AtomDict>::iterator iA=iR->atoms.begin();
                    iA !=iR->atoms.end(); iA++)
            {
                std::cout << "atom " << iA->id << std::endl;
            }
            std::cout << std::endl;
        }
         */

        for (std::vector<AngleDict>::iterator iAn = tMol.angles.begin();
                iAn != tMol.angles.end(); iAn++) {
            iAn->isInSameRing = checkIfAngleInSameRing(tMol.atoms,
                    tMol.rings,
                    iAn->atoms[0],
                    iAn->atoms[1], iAn->atoms[2]);
        }

        /*
        for (std::vector<AngleDict>::iterator iA=angles.begin();
                iA !=angles.end(); iA++)
        {
            std::cout << "angle between " << tMol.atoms[iA->atoms[0]].id 
                      << "(center) and " <<  tMol.atoms[iA->atoms[1]].id
                      << " and " << tMol.atoms[iA->atoms[2]].id 
                      << " is " << iA->value*PID180 << std::endl;
        }
         */
    }

    REAL MolGenerator::getAngleValueFromFracCoords(AtomDict & tAtCen,
            AtomDict & tAt1,
            AtomDict & tAt2,
            REAL a, REAL b, REAL c,
            REAL alpha, REAL beta, REAL gamma) {


        if (tAtCen.fracCoords.size() == 3 && tAtCen.fracCoords.size() == tAt1.fracCoords.size()
                && tAtCen.fracCoords.size() == tAt2.fracCoords.size()) {
            std::vector<REAL> tV1, tV2;
            for (unsigned i = 0; i < 3; i++) {
                tV1.push_back(tAt1.fracCoords[i] - tAtCen.fracCoords[i]);
                tV2.push_back(tAt2.fracCoords[i] - tAtCen.fracCoords[i]);
            }

            REAL lenTv1 = getBondLenFromFracCoords(tAtCen.fracCoords, tAt1.fracCoords,
                    a, b, c, alpha, beta, gamma);
            REAL lenTv2 = getBondLenFromFracCoords(tAtCen.fracCoords, tAt2.fracCoords,
                    a, b, c, alpha, beta, gamma);
            //std::cout << "lenTv1 " << lenTv1 << std::endl;
            //std::cout << "lenTv2 " << lenTv2 << std::endl;

            if (lenTv1 > 0.000001 && lenTv2 > 0.000001) {
                REAL coF = (a * a * tV1[0] * tV2[0] + b * b * tV1[1] * tV2[1] + c * c * tV1[2] * tV2[2]
                        + b * c * (tV1[1] * tV2[2] + tV1[2] * tV2[1]) * cos(alpha * PI180)
                        + c * a * (tV1[2] * tV2[0] + tV1[0] * tV2[2]) * cos(beta * PI180)
                        + a * b * (tV1[0] * tV2[1] + tV1[1] * tV2[0]) * cos(gamma * PI180)) / (lenTv1 * lenTv2);

                //std::cout << "coF " << coF << std::endl;
                // The float point error : ABS might > 1.0 or < -1.0
                if (coF >= 1.0) {
                    coF -= 0.0000000001;
                } else if (coF <= -1.0) {
                    coF += 0.0000000001;
                }
                return acos(coF);
            } else {
                std::cout << "There is at least one pair of atoms overlapped " << std::endl;
                return 0.0;
            }
        } else {
            std::cout << "Error: check atom coordinate dimensions " << std::endl;
            return 0.0;
        }

    }

    void MolGenerator::getMolByEqClassInCell() {

        std::map<int, int> classNum;
        classNum[0] = 0;
        for (unsigned i = 1; i < allAtoms.size(); i++) {
            classNum[i] = i;
            for (unsigned j = 0; j <= i - 1; j++) {
                classNum[j] = classNum[classNum[j]];
                if (std::find(allAtoms[i].connAtoms.begin(),
                        allAtoms[i].connAtoms.end(), j)
                        != allAtoms[i].connAtoms.end()) {
                    classNum[classNum[classNum[j]]] = i;
                }
            }
        }

        // final sweeping 
        for (unsigned i = 0; i < allAtoms.size(); i++) {
            classNum[i] = classNum[classNum[i]];
        }


        std::map<unsigned, std::vector<int> > aMoleculesInCell;
        // get molecules in a unit cell
        for (unsigned i = 0; i < classNum.size(); i++) {
            aMoleculesInCell[classNum[i]].push_back(i);
        }

        //std::cout << "Number of molecules after EQ classification is " 
        //          << aMoleculesInCell.size() << std::endl;

        //deleteNonCenCellMols(aMoleculesInCell);
        deleteNonASUAtomCellMols(aMoleculesInCell);

        std::cout << "Number of molecules containing atoms in ASU is "
                << moleculesInCell.size() << std::endl;
        /*
        for (std::map<unsigned, std::vector<int> >::iterator iMol=moleculesInCell.begin(); 
                  iMol !=moleculesInCell.end(); iMol++)
        {
            std::cout << "Molecule " << iMol->first << " contains " 
                      << iMol->second.size() << " atoms " << std::endl 
                      << "The are: " << std::endl;
            for (std::vector<int>::iterator iAt=iMol->second.begin();
                    iAt !=iMol->second.end(); iAt++)
            {
                std::cout << "Atom " << allAtoms[*iAt].seriNum << " "
                          << allAtoms[*iAt].id  << " with connections "
                          << allAtoms[*iAt].connAtoms.size() <<  std::endl;
            }                   
        }
       
         */

    }

    void MolGenerator::buildAndValidMols(PeriodicTable & tPTab,
            std::vector<CrystInfo>::iterator tCryst) {

        std::cout << "Number of molecules after linked equiv class: "
                << moleculesInCell.size() << std::endl;

        // std::vector<Molecule> tmpMols;

        int nMol = 0;
        int allAsy = 0;
        for (std::map<unsigned, std::vector<int> >::iterator iMol = moleculesInCell.begin();
                iMol != moleculesInCell.end(); iMol++) {
            // Check if these molecule with atom index only contain atoms in ASU.
            // If yes, materialize the molecule.

            // check if there is any atom with occupancy less than 0.95 
            // and then set individual molecules
            if (isASUAtomInMol(iMol)) {

                std::cout << "\n---------------------------------------" << std::endl;

                Molecule aMol;
                aMol.seriNum = nMol;
                nMol++;
                //std::cout << "Its serial number " 
                //          << aMol.seriNum << std::endl;
                //std::cout << "The molecule contains "  
                //          << iMol->second.size() << std::endl;

                std::map<int, int> tAtmMap;
                int i = 0;
                for (std::vector<int>::iterator iAt = iMol->second.begin();
                        iAt != iMol->second.end(); iAt++) {
                    //std::cout << "atom serial number " << *iAt << std::endl;
                    tAtmMap[*iAt] = i;
                    aMol.atoms.push_back(allAtoms[*iAt]);
                    i++;
                }

                aMol.setFormula();
                //std::cout << "The formula for the current molecule is "
                //          << aMol.formula << std::endl;

                int nAsy = 0;
                for (std::vector<AtomDict>::iterator iAt = aMol.atoms.begin();
                        iAt != aMol.atoms.end(); iAt++) {
                    if (iAt->isInPreCell) {
                        nAsy++;
                    }
                }
                allAsy += nAsy;
                //std::cout << "For this molecule, number of atoms in the assym cell is "
                //          << nAsy << std::endl;

                // Change the connection index for atoms in that molecule.
                i = 0;
                for (std::vector<AtomDict>::iterator iAt = aMol.atoms.begin();
                        iAt != aMol.atoms.end(); iAt++) {
                    iAt->seriNum = i;
                    i++;
                    for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                            iCo != iAt->connAtoms.end(); iCo++) {
                        (*iCo) = tAtmMap[(*iCo)];
                    }
                }

                std::cout << "The number of atoms in molecule " << aMol.seriNum
                        << " is " << aMol.atoms.size() << std::endl;



                //   << "They are : " << std::endl;
                /*
                for (std::vector<AtomDict>::iterator iAt=aMol.atoms.begin();
                        iAt !=aMol.atoms.end(); iAt++)
                {
                    if (iAt->isInPreCell)
                    {
                        std::cout <<"Atoms " << iAt->id << "(" << iAt->seriNum
                                  << ") in cell " << iAt->sId << std::endl;
                        if (iAt->connAtoms.size() !=0)
                        {
                            std::cout << "======================" << std::endl;
                            std::cout << "It connects: " << std::endl; 
                            for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                                  iCo !=iAt->connAtoms.end(); iCo++)
                            {
                                std::cout << "Atom " << aMol.atoms[*iCo].id << "("
                                          << aMol.atoms[*iCo].seriNum << ") in cell " 
                                          << aMol.atoms[*iCo].sId << std::endl;
                            }
                            std::cout << "======================" << std::endl;
                        }
                    }
                }
               
                */

                //getAtomTypeOneMol(aMol);
                getAtomTypeOneMolNew(aMol);

                // getUniqueBondsMols(aMol, tCryst);
                getAllBondsInOneMol(aMol, tCryst);


                std::string aMolErrInfo;
                bool aPass = validateMolecule(aMol, tPTab, aMolErrInfo);

                if (aPass && aMol.bonds.size() > 0) {
                    //getAtomTypeOneMol(aMol);
                    getUniqAngleMols(aMol, tCryst);
                    setAnglesSPSigns(aMol.atoms, aMol.angles);
                    // Change the molecule serial number 
                    allMolecules.push_back(aMol);
                    std::string aMsg;
                    aMsg = "Molecule " + IntToStr(aMol.seriNum) + " is included.\n"
                            + "It contains " + IntToStr(aMol.atoms.size())
                            + " atoms\n";
                    validedMolMsg[aMol.seriNum] = aMsg;

                    std::cout << "Molecule " << aMol.seriNum << " is included "
                            << std::endl << "It contains " << aMol.atoms.size()
                            << " atoms" << std::endl;
                }
                else {
                    std::string aMsg = aMolErrInfo + "\n"
                            + "Molecule " + IntToStr(aMol.seriNum)
                            + " is rejected\n";
                    errMolMsg[aMol.seriNum] = aMsg;

                    std::cout << aMolErrInfo << std::endl;
                    std::cout << "Molecule " << aMol.seriNum
                            << " is rejected " << std::endl;
                }
            }
        }

        std::cout << "Number of molecules in allMolecules after validation "
                << allMolecules.size() << std::endl;
        std::cout << "Total number of atoms in Asy cell covered "
                << allAsy << std::endl;
        std::cout << "Final check, whole structure " << std::endl;

        std::string errInfoStr;
        bool tOkStruct = validateBondValueDiffStruct(allMolecules, errInfoStr);
        if (!tOkStruct) {
            allMolecules.clear();
            std::cout << errInfoStr << std::endl;
        }

    }

    bool MolGenerator::checEquiMoles(std::vector<Molecule> & tSetMols,
            Molecule& tMol) {
        bool tAccept = true;

        if (tMol.atoms.size() == 0 || tMol.formula.size() == 0) {
            tAccept = false;
        } else {
            std::vector<std::string> atmIDs;
            for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                    iAt != tMol.atoms.end(); iAt++) {
                atmIDs.push_back(iAt->id);
            }

            for (std::vector<Molecule>::iterator iMol = tSetMols.begin();
                    iMol != tSetMols.end(); iMol++) {
                if (iMol->seriNum != tMol.seriNum) {
                    if ((iMol->atoms.size() == tMol.atoms.size() &&
                            iMol->formula.compare(tMol.formula) == 0)) {
                        for (std::vector<AtomDict>::iterator iA = iMol->atoms.begin();
                                iA != iMol->atoms.end(); iA++) {
                            if (iA->chemType.find("H") == std::string::npos) {
                                if ((std::find(atmIDs.begin(), atmIDs.end(), iA->id)
                                        != atmIDs.end())) {
                                    tAccept = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        return tAccept;
    }

    bool MolGenerator::checkAsmAtomsInMol(Molecule& tMol,
            std::string& tErrInfo) {
        bool aReturn = true;

        int nAsym = 0;

        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            if (iAt->isInPreCell) {
                nAsym++;
            }
        }

        if (nAsym == 0) {
            tErrInfo = "Mol " + IntToStr(tMol.seriNum)
                    + " contains no Asym cell atoms \n";
            aReturn = false;
        }
        std::cout << "Mol " << tMol.seriNum << " contains "
                << nAsym << " Asym cell atoms " << std::endl;

        return aReturn;
    }

    bool MolGenerator::checkAtomOcp(Molecule & tMol,
            std::string & tErrInfo) {

        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            //std::cout << "Atom " << iAt->id << " occupancy "
            //           << iAt->ocp << std::endl;

            if (iAt->ocp < 0.99) {
                tErrInfo = "Atom " + iAt->id + " has an occupancy "
                        + RealToStr(iAt->ocp) + ", small than 1.0!";

                return false;
            }

        }

        return true;

    }

    void MolGenerator::getMolByEqClassInCrys() {

        std::map<int, int> classNum;
        classNum[0] = 0;
        for (unsigned i = 1; i < refAtoms.size(); i++) {
            classNum[i] = i;
            for (unsigned j = 0; j <= i - 1; j++) {
                classNum[j] = classNum[classNum[j]];
                if (std::find(refAtoms[i].connAtoms.begin(),
                        refAtoms[i].connAtoms.end(), j)
                        != refAtoms[i].connAtoms.end()) {
                    classNum[classNum[classNum[j]]] = i;
                }
            }
        }

        // final sweeping 
        for (unsigned i = 0; i < refAtoms.size(); i++) {
            classNum[i] = classNum[classNum[i]];
        }


        // get molecules in a finite crystal 
        for (unsigned i = 0; i < classNum.size(); i++) {
            moleculesInCryst[classNum[i]].push_back(i);
        }
    }

    void MolGenerator::deleteNonCenCellMols(std::map<unsigned, std::vector<int> >
            & tMoleculesInCell) {
        std::map<unsigned, bool> tMolMap;
        for (std::map<unsigned, std::vector<int> >::iterator iM = tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++) {
            bool tIn = false;
            for (std::vector<int>::iterator iAt = iM->second.begin();
                    iAt != iM->second.end(); iAt++) {
                if (allAtoms[*iAt].sId == "555") {
                    tIn = true;
                    break;
                }
            }

            tMolMap[iM->first] = tIn;
            /*
            std::cout << "Mol " << iM->first;
            if (tIn)
            { 
                std::cout << " contains  atoms in the center cell " << std::endl;
            }
            else
            {
                std::cout << " does not contains  atoms in the center cell " << std::endl;
            }
             */
        }

        // Delete those molecules which contain no atoms in the center cell

        for (std::map<unsigned, bool>::iterator iM = tMolMap.begin();
                iM != tMolMap.end(); iM++) {
            if (!iM->second) {
                tMoleculesInCell.erase(iM->first);
            }
        }


        unsigned iMol = 1;
        for (std::map<unsigned, std::vector<int> >::iterator iM = tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++) {
            for (std::vector<int>::iterator iAt = iM->second.begin();
                    iAt != iM->second.end(); iAt++) {
                moleculesInCell[iMol].push_back(*iAt);
            }

            iMol++;
        }

        std::cout << "Only " << moleculesInCell.size()
                << " need to be considered. " << std::endl;
        // "They are: " << std::endl;
        /*
        for (std::map<unsigned, std::vector<int> >::iterator iM=moleculesInCell.begin();
                iM !=moleculesInCell.end(); iM++)
        {
            std::cout << "Molecule " << iM->first << " which contains: " << std::endl;
            for (std::vector<int>::iterator iAt=iM->second.begin();
                    iAt !=iM->second.end(); iAt++)
            {
                std::cout << "Atom " << allAtoms[*iAt].id 
                          <<"(" << allAtoms[*iAt].seriNum << ")" 
                          << " in cell " << allAtoms[*iAt].sId << std::endl;
                std::cout << "It connects to: " << std::endl;
                for (std::vector<int>::iterator iNB=allAtoms[*iAt].connAtoms.begin();
                        iNB != allAtoms[*iAt].connAtoms.end(); iNB++)
                {
                    std::cout << "Atom " << allAtoms[*iNB].id << "("
                              << allAtoms[*iNB].seriNum << ")"
                              << " in cell " << allAtoms[*iNB].sId << std::endl;
                }
                          
                        
            }
                      
        }
         */

    }

    void MolGenerator::deleteNonASUAtomCellMols(std::map<unsigned, std::vector<int> >& tMoleculesInCell) {


        std::map<unsigned, bool> tMolMap;
        for (std::map<unsigned, std::vector<int> >::iterator iM = tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++) {


            tMolMap[iM->first] = isASUAtomInMol(iM);
            /*
            std::cout << "Mol " << iM->first;
            if (tMolMap[iM->first])
            { 
                std::cout << " contains  atoms in the ASU" << std::endl;
            }
            else
            {
                std::cout << " does not contains  atoms in the ASU" << std::endl;
            }
             **/
        }

        // Delete those molecules which contain no atoms in the center cell

        for (std::map<unsigned, bool>::iterator iM = tMolMap.begin();
                iM != tMolMap.end(); iM++) {
            if (!iM->second) {
                tMoleculesInCell.erase(iM->first);
            }
        }


        unsigned iMol = 1;
        for (std::map<unsigned, std::vector<int> >::iterator iM = tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++) {
            for (std::vector<int>::iterator iAt = iM->second.begin();
                    iAt != iM->second.end(); iAt++) {
                moleculesInCell[iMol].push_back(*iAt);
            }

            iMol++;
        }

        //std::cout << "Only " << moleculesInCell.size() 
        //          << " need to be considered. They are: " << std::endl;
        /*
        for (std::map<unsigned, std::vector<int> >::iterator iM=moleculesInCell.begin();
                iM !=moleculesInCell.end(); iM++)
        {
            std::cout << "Molecule " << iM->first << " which contains: " << std::endl;
            for (std::vector<int>::iterator iAt=iM->second.begin();
                    iAt !=iM->second.end(); iAt++)
            {
                std::cout << "Atom " << allAtoms[*iAt].id 
                          <<"(" << allAtoms[*iAt].seriNum << ")" 
                          << " in cell " << allAtoms[*iAt].sId << std::endl;
                std::cout << "It connects to: " << std::endl;
                for (std::vector<int>::iterator iNB=allAtoms[*iAt].connAtoms.begin();
                        iNB != allAtoms[*iAt].connAtoms.end(); iNB++)
                {
                    std::cout << "Atom " << allAtoms[*iNB].id << "("
                              << allAtoms[*iNB].seriNum << ")"
                              << " in cell " << allAtoms[*iNB].sId << std::endl;
                }
                          
                        
            }
                      
        }
         */

    }

    void MolGenerator::checkAtomElementID(std::vector<AtomDict> & tAtoms) {
        PeriodicTable aPTab;

        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++) {
            if (iAt->chemType.empty() ||
                    aPTab.elements.find(iAt->chemType) == aPTab.elements.end()) {
                if (!iAt->id.empty()) {
                    if (iAt->id.size() >= 2) {
                        if (aPTab.elements.find(iAt->id.substr(0, 2))
                                != aPTab.elements.end()) {
                            iAt->chemType = iAt->id.substr(0, 2);
                        } else if (aPTab.elements.find(iAt->id.substr(0, 1))
                                != aPTab.elements.end()) {
                            iAt->chemType = iAt->id.substr(0, 1);
                        }
                    } else if (iAt->id.size() == 1) {
                        if (aPTab.elements.find(iAt->id)
                                != aPTab.elements.end()) {
                            iAt->chemType = iAt->id;
                        }
                    }
                } else {
                    std::cout << "What is atom " << iAt->id << "'s element type "
                            << std::endl;
                    exit(1);
                }
            }
        }

    }

    bool MolGenerator::colidAtom(AtomDict& tAtom,
            std::vector<AtomDict>& tRefAtoms, int tMode) {
        for (std::vector<AtomDict>::iterator iAt = tRefAtoms.begin();
                iAt != tRefAtoms.end(); iAt++) {
            REAL da = fabs(tAtom.fracCoords[0] - iAt->fracCoords[0]);
            REAL db = fabs(tAtom.fracCoords[1] - iAt->fracCoords[1]);
            REAL dc = fabs(tAtom.fracCoords[2] - iAt->fracCoords[2]);
            if (da < 0.01 && db < 0.01 && dc < 0.01) {
                if (tMode == 1) {
                    allMsg.push_back("REJECTED STRUCTURE:  atom collisions.");
                    allMsg.push_back("Atom " + tAtom.id + " has collision during symmgen \n");
                    std::cout << "Atom " << tAtom.id << " has collision during symmgen \n";
                    lColid = true;
                }
                return true;
            }
        }
        return false;
    }

    bool MolGenerator::colidAtom(std::vector<REAL> & tFracX,
            std::vector<AtomDict>& tRefAtoms, int tMode) {
        for (std::vector<AtomDict>::iterator iAt = tRefAtoms.begin();
                iAt != tRefAtoms.end(); iAt++) {
            REAL da = fabs(tFracX[0] - iAt->fracCoords[0]);
            REAL db = fabs(tFracX[1] - iAt->fracCoords[1]);
            REAL dc = fabs(tFracX[2] - iAt->fracCoords[2]);
            if (da < 0.005 && db < 0.005 && dc < 0.005) {
                if (tMode == 1) {
                    allMsg.push_back("REJECTED STRUCTURE:  atom collisions.");
                    allMsg.push_back("Atom " + iAt->id + " has collision during symmgen \n");
                    std::cout << "Atom " << iAt->id << " has collision during symmgen \n";
                    lColid = true;
                }
                return true;
            }
        }

        return false;
    }

    bool MolGenerator::isASUAtomInMol(std::map<unsigned, std::vector<int> >::iterator tMol) {
        for (std::vector<int>::iterator iAt = tMol->second.begin();
                iAt != tMol->second.end(); iAt++) {
            if (allAtoms[*iAt].isInPreCell) {
                return true;
            }
        }

        return false;
    }

    bool MolGenerator::connMetal(std::vector<int> & tIdxs,
            std::vector<AtomDict>& tAtoms) {
        for (std::vector<int>::iterator iIdx = tIdxs.begin();
                iIdx != tIdxs.end(); iIdx++) {
            if (tAtoms[*iIdx].isMetal) {
                return true;
            }
        }

        return false;
    }

    bool MolGenerator::connMetal2ndNB(std::vector<int>& tIdxs,
            std::vector<AtomDict>& tAtoms) {
        //std::cout << "Check if two atoms in a bond connect to a metal atom"
        //          << std::endl;

        for (std::vector<int>::iterator iIdx = tIdxs.begin();
                iIdx != tIdxs.end(); iIdx++) {
            // std::cout << "Atom " << tAtoms[*iIdx].id << " : " << std::endl;
            if (tAtoms[*iIdx].isMetal) {
                return true;
            }
            for (std::vector<int>::iterator iNB = tAtoms[*iIdx].connAtoms.begin();
                    iNB != tAtoms[*iIdx].connAtoms.end(); iNB++) {
                if (tAtoms[*iNB].isMetal) {
                    std::cout << "its 1stNB " << tAtoms[*iNB].id
                            << " is metal" << std::endl;
                    return true;
                }
                for (std::vector<int>::iterator i2NB = tAtoms[*iNB].connAtoms.begin();
                        i2NB != tAtoms[*iNB].connAtoms.end(); i2NB++) {
                    if (tAtoms[*i2NB].isMetal) {
                        std::cout << "its 2ndNB " << tAtoms[*i2NB].id
                                << " is metal" << std::endl;
                        return true;
                    }
                }
            }
        }

        return false;
    }

    bool MolGenerator::validateBonds(std::vector<BondDict>::iterator tBo,
            std::string & tErrInfo,
            PeriodicTable & tPTab) {
        // Bond properties (number of connection) check
        std::cout << "Bond between atoms " << allAtoms[tBo->atomsIdx[0]].id
                << " and " << allAtoms[tBo->atomsIdx[1]].id << std::endl;

        if ((allAtoms[tBo->atomsIdx[0]].connAtoms.size() == 1 &&
                allAtoms[tBo->atomsIdx[1]].connAtoms.size() == 1)) {
            tErrInfo = "Both atoms connect to only one atom for atom"
                    + allAtoms[tBo->atomsIdx[0]].id + " and "
                    + allAtoms[tBo->atomsIdx[1]].id;
            return false;
        } else if ((allAtoms[tBo->atomsIdx[0]].connAtoms.size() > 4 &&
                !allAtoms[tBo->atomsIdx[0]].isMetal)
                || (allAtoms[tBo->atomsIdx[1]].connAtoms.size() > 4 &&
                !allAtoms[tBo->atomsIdx[1]].isMetal)) {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[0]].connAtoms.size())
                    + "Atom " + allAtoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[1]].connAtoms.size());

            return false;
        }// C
        else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("C") == 0
                && allAtoms[tBo->atomsIdx[0]].connAtoms.size() < 3)
                ||
                (allAtoms[tBo->atomsIdx[1]].chemType.compare("C") == 0
                && allAtoms[tBo->atomsIdx[1]].connAtoms.size() < 3)) {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[0]].connAtoms.size())
                    + "Atom " + allAtoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[1]].connAtoms.size());
            return false;
        }// B
        else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("B") == 0
                && allAtoms[tBo->atomsIdx[0]].connAtoms.size() < 2)
                ||
                (allAtoms[tBo->atomsIdx[1]].chemType.compare("B") == 0
                && allAtoms[tBo->atomsIdx[1]].connAtoms.size() < 2)) {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[0]].connAtoms.size())
                    + "Atom " + allAtoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[1]].connAtoms.size());
            return false;
        } else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("H") == 0
                && allAtoms[tBo->atomsIdx[0]].connAtoms.size() > 1)
                ||
                (allAtoms[tBo->atomsIdx[1]].chemType.compare("H") == 0
                && allAtoms[tBo->atomsIdx[1]].connAtoms.size() > 1)) {
            tErrInfo = "Reject the molecule! Atom " + allAtoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[0]].connAtoms.size())
                    + "Atom " + allAtoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                    + " has connections " + IntToStr((int) allAtoms[tBo->atomsIdx[1]].connAtoms.size());
            return false;
        }

        // bond length check again, this time we use both high and 
        // low bounds of a bond length

        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 = 0.20;
        REAL covalent_sensitivity2 = 0.50;
        REAL covalent_sensitivity3 = 0.30;
        REAL covalent_sensitivity4 = 0.80;

        std::vector<REAL> linkRange;
        if ((allAtoms[tBo->atomsIdx[0]].chemType == "H" && allAtoms[tBo->atomsIdx[1]].chemType == "O")
                || (allAtoms[tBo->atomsIdx[0]].chemType == "O" && allAtoms[tBo->atomsIdx[1]].chemType == "H")) {
            covalent_sensitivity = covalent_sensitivity4;
        } else if (allAtoms[tBo->atomsIdx[1]].chemType == "O"
                || allAtoms[tBo->atomsIdx[0]].chemType == "O") {
            covalent_sensitivity = covalent_sensitivity3;
        } else if (allAtoms[tBo->atomsIdx[1]].chemType == "H"
                || allAtoms[tBo->atomsIdx[0]].chemType == "H") {
            covalent_sensitivity = covalent_sensitivity2;
        } else if ((!allAtoms[tBo->atomsIdx[0]].isMetal) && (!allAtoms[tBo->atomsIdx[1]].isMetal)) {
            covalent_sensitivity = covalent_sensitivity1;
        } else {
            covalent_sensitivity = covalent_sensitivity2;
        }
        getBondingRangePairAtoms2(allAtoms[tBo->atomsIdx[0]], allAtoms[tBo->atomsIdx[1]],
                covalent_sensitivity, tPTab, linkRange);

        if (tBo->value < linkRange[0] && tBo->value > linkRange[1]) {
            tErrInfo = "Bond between " + allAtoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                    + " and " + allAtoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                    + " is " + RealToStr(tBo->value)
                    + "It should be between " + RealToStr(linkRange[0]) + " and "
                    + RealToStr(linkRange[1]);
            return false;
        }

        return true;

    }

    bool MolGenerator::validateBonds(std::vector<BondDict>::iterator tBo,
            Molecule & tMol,
            std::string & tErrInfo,
            PeriodicTable & tPTab) {
        // Bond properties (number of connection) check
        /*
        std::cout << "Bond between atoms " << tMol.atoms[tBo->atomsIdx[0]].id
                  << " and " << tMol.atoms[tBo->atomsIdx[1]].id << std::endl;
        std::cout << "Atom " << tMol.atoms[tBo->atomsIdx[0]].id << " connects "
                  << tMol.atoms[tBo->atomsIdx[0]].connAtoms.size() << std::endl
                  << "Atom " << tMol.atoms[tBo->atomsIdx[1]].id << " connects "
                  << tMol.atoms[tBo->atomsIdx[1]].connAtoms.size() << std::endl;
        std::cout << "It's element ID " << tMol.atoms[tBo->atomsIdx[0]].chemType << std::endl;
        std::cout << "It's element ID " << tMol.atoms[tBo->atomsIdx[1]].chemType << std::endl;
        std::cout << "comp H " << tMol.atoms[tBo->atomsIdx[1]].chemType.compare("H") << std::endl;
        std::cout << "comp H " << tMol.atoms[tBo->atomsIdx[0]].chemType.compare("H") << std::endl;
        std::cout << "comp H " << tMol.atoms[tBo->atomsIdx[1]].chemType.compare("H") << std::endl;
         */
        if ((tMol.atoms[tBo->atomsIdx[0]].connAtoms.size() == 1 &&
                tMol.atoms[tBo->atomsIdx[1]].connAtoms.size() == 1)) {
            tErrInfo = "Both atoms connect to only one atom for atom"
                    + tMol.atoms[tBo->atomsIdx[0]].id + " and "
                    + tMol.atoms[tBo->atomsIdx[1]].id;
            return false;
        }

        // bond length check again, this time we use BOTH(!) the high and 
        // low bounds of a bond length

        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 = 0.20;
        REAL covalent_sensitivity2 = 0.50;
        REAL covalent_sensitivity3 = 0.30;
        REAL covalent_sensitivity4 = 0.80;

        std::vector<REAL> linkRange;
        if ((tMol.atoms[tBo->atomsIdx[0]].chemType == "H" && tMol.atoms[tBo->atomsIdx[1]].chemType == "O")
                || (tMol.atoms[tBo->atomsIdx[0]].chemType == "O" && tMol.atoms[tBo->atomsIdx[1]].chemType == "H")) {
            covalent_sensitivity = covalent_sensitivity4;
        } else if (tMol.atoms[tBo->atomsIdx[1]].chemType == "O"
                || tMol.atoms[tBo->atomsIdx[0]].chemType == "O") {
            covalent_sensitivity = covalent_sensitivity3;
        } else if (allAtoms[tBo->atomsIdx[1]].chemType == "H"
                || allAtoms[tBo->atomsIdx[0]].chemType == "H") {
            covalent_sensitivity = covalent_sensitivity2;
        } else if ((!tMol.atoms[tBo->atomsIdx[0]].isMetal) && (!tMol.atoms[tBo->atomsIdx[1]].isMetal)) {
            covalent_sensitivity = covalent_sensitivity1;
        } else {
            covalent_sensitivity = covalent_sensitivity2;
        }
        getBondingRangePairAtoms2(tMol.atoms[tBo->atomsIdx[0]], tMol.atoms[tBo->atomsIdx[1]],
                covalent_sensitivity, tPTab, linkRange);
        //std::cout << "atom1 " << tMol.atoms[tBo->atomsIdx[0]].id << std::endl 
        //          << "atom2 " << tMol.atoms[tBo->atomsIdx[1]].id << std::endl;
        //std::cout << "bond: low boundary  " << linkRange[0] <<  std::endl;
        //std::cout << "high boundary " << linkRange[1] << std::endl;
        //std::cout << "bond value  " << tBo->value << std::endl; 

        if (tBo->value < linkRange[0] || tBo->value > linkRange[1]) {
            tErrInfo = "Bond between " + tMol.atoms[tBo->atomsIdx[0]].id
                    + " serial number " + IntToStr(tMol.atoms[tBo->atomsIdx[0]].seriNum)
                    + " and " + tMol.atoms[tBo->atomsIdx[1]].id
                    + " serial number " + IntToStr(tMol.atoms[tBo->atomsIdx[1]].seriNum)
                    + " is " + RealToStr(tBo->value)
                    + "It should be between " + RealToStr(linkRange[0]) + " and "
                    + RealToStr(linkRange[1]);
            std::cout << tErrInfo << std::endl;
            return false;
        }
        //else
        //{
        //    std::cout << "The bond is validated " << std::endl;
        //}

        return true;

    }

    bool MolGenerator::validateBondValueSame(Molecule& tMol,
            std::string& tErrInfo) {
        bool aReturn = true;

        std::map<ID, std::vector<BondDict> > aBondMap;

        for (std::vector<BondDict>::iterator iBo = tMol.bonds.begin();
                iBo != tMol.bonds.end(); iBo++) {
            // Check bonds with atoms in asymmetric unit only.
            // symmetry generated atoms allowed to duplicate bond values
            // with the original one. But we only include bonds with at least 
            // one atoms in an asymmetrical unit. 
            // This step moved here from Python scripts. Remember cancel the 
            // step there !!!!

            std::string sValue = RealToStr(iBo->value);
            /*
            std::cout << "Check bond between " << tMol.atoms[iBo->atomsIdx[0]].id
                      << " and " << tMol.atoms[iBo->atomsIdx[1]].id
                      << ". The value is " << sValue << std::endl;
            std::cout << "Is atom " << "in asymm cell ? " 
                      << tMol.atoms[iBo->atomsIdx[0]].isInPreCell 
                      << std::endl;
            std::cout << "Is atom " << tMol.atoms[iBo->atomsIdx[1]].id
                      << "in asymm cell ? " << tMol.atoms[iBo->atomsIdx[1]].isInPreCell 
                      << std::endl;
             */

            if (tMol.atoms[iBo->atomsIdx[0]].isInPreCell
                    || tMol.atoms[iBo->atomsIdx[1]].isInPreCell) {

                if (aBondMap.find(sValue) == aBondMap.end()) {

                    aBondMap[sValue].push_back(*iBo);
                } else {
                    tErrInfo = "Bond lengths " + sValue
                            + " duplicated among the bond between ";
                    tErrInfo += (" atoms " + tMol.atoms[iBo->atomsIdx[0]].id
                            + " and " + tMol.atoms[iBo->atomsIdx[1]].id
                            + ", ");

                    tErrInfo += ("and the bond between atoms " +
                            aBondMap[sValue][0].atoms[0]
                            + " and " + aBondMap[sValue][0].atoms[1]
                            + " \n");
                    aReturn = false;
                    break;
                }
            }
        }

        return aReturn;
    }

    bool MolGenerator::validateBondValueDiff(Molecule& tMol,
            std::string& tErrInfo) {
        std::map<ID, std::vector<REAL> > aBondMap;
        REAL diffThreshold = 0.04;

        for (std::vector<BondDict>::iterator iBo = tMol.bonds.begin();
                iBo != tMol.bonds.end(); iBo++) {
            if (tMol.atoms[iBo->atomsIdx[0]].chemType.compare("H") != 0
                    && tMol.atoms[iBo->atomsIdx[1]].chemType.compare("H")) {
                std::list<std::string> tCids;
                tCids.push_back(tMol.atoms[iBo->atomsIdx[0]].codClass);
                tCids.push_back(tMol.atoms[iBo->atomsIdx[1]].codClass);
                tCids.sort(compareNoCase2);

                std::string aCombID;
                int nRS = 0;
                for (std::list<std::string>::iterator iId = tCids.begin();
                        iId != tCids.end(); iId++) {
                    if (nRS == 0) {
                        aCombID.append(*iId);
                    } else {
                        aCombID.append("_" + *iId);
                    }
                    nRS++;
                }
                // std::cout << "bond_atom_IDS = " <<  aCombID
                //          << std::endl;

                if (aBondMap.find(aCombID) == aBondMap.end()) {
                    // std::cout << iBo->value << " added" << std::endl;
                    aBondMap[aCombID].push_back(iBo->value);
                } else {
                    // std::cout << iBo->value << " under check " << std::endl;
                    if (!tMol.atoms[iBo->atomsIdx[0]].isMetal
                            && !tMol.atoms[iBo->atomsIdx[1]].isMetal) {
                        if (outVectAbsDiff(aBondMap[aCombID], iBo->value, diffThreshold)) {
                            tErrInfo.append("bond between atom " + tMol.atoms[iBo->atomsIdx[0]].id +
                                    " and " + tMol.atoms[iBo->atomsIdx[1]].id +
                                    " is very different with another bond of the same atom types."
                                    + " Molecule rejected \n");
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    bool MolGenerator::validateBondValueDiffStruct(std::vector<Molecule> & tMols,
            std::string& tErrInfo) {
        std::map<ID, std::vector<REAL> > aBondMap;
        REAL diffThreshold = 0.05;
        for (std::vector<Molecule>::iterator iMol = tMols.begin();
                iMol != tMols.end(); iMol++) {
            for (std::vector<BondDict>::iterator iBo = iMol->bonds.begin();
                    iBo != iMol->bonds.end(); iBo++) {
                if (iMol->atoms[iBo->atomsIdx[0]].chemType.compare("H") != 0
                        && iMol->atoms[iBo->atomsIdx[1]].chemType.compare("H")) {
                    std::list<std::string> tCids;
                    tCids.push_back(iMol->atoms[iBo->atomsIdx[0]].codClass);
                    tCids.push_back(iMol->atoms[iBo->atomsIdx[1]].codClass);
                    tCids.sort(compareNoCase2);

                    std::string aCombID;
                    int nRS = 0;
                    for (std::list<std::string>::iterator iId = tCids.begin();
                            iId != tCids.end(); iId++) {
                        if (nRS == 0) {
                            aCombID.append(*iId);
                        } else {
                            aCombID.append("_" + *iId);
                        }
                        nRS++;
                    }
                    // std::cout << "bond_atom_IDS = " <<  aCombID
                    //          << std::endl;

                    if (aBondMap.find(aCombID) == aBondMap.end()) {
                        // std::cout << iBo->value << " added" << std::endl;
                        aBondMap[aCombID].push_back(iBo->value);
                    } else {
                        // std::cout << iBo->value << " under check " << std::endl;

                        if (outVectAbsDiff(aBondMap[aCombID], iBo->value, diffThreshold)) {
                            tErrInfo.append("REJECTED STRCUTRUE: DIFF OF BOND VALUES TOO BIG\nbond between atom "
                                    + iMol->atoms[iBo->atomsIdx[0]].id +
                                    " and " + iMol->atoms[iBo->atomsIdx[1]].id +
                                    " is very different with another bond of the same atom types."
                                    + " Molecule rejected \n");
                            return false;
                        }
                    }
                }
            }
        }

        return true;
    }

    bool MolGenerator::validateAtomLinks(Molecule & tMol,
            PeriodicTable & tPTab,
            std::string & tErrInfo) {
        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            if (iAt->isInPreCell) {
                int nMax = 4;
                if (iAt->chemType.compare("B") == 0) {
                    nMax = 5;
                }
                if (iAt->connAtoms.size() > nMax && !iAt->isMetal
                        && !connMetal(iAt->connAtoms, tMol.atoms)) {
                    tErrInfo = "Reject Molecule! Atom " + iAt->id
                            + " serial number " + IntToStr(iAt->seriNum)
                            + " has connections " + IntToStr((int) iAt->connAtoms.size());
                    return false;
                }// C
                else if (iAt->chemType.compare("C") == 0 && iAt->connAtoms.size() < 2) {
                    tErrInfo = "Reject Molecule! Atom " + iAt->id
                            + " serial number " + IntToStr(iAt->seriNum)
                            + " has connections " + IntToStr((int) iAt->connAtoms.size());
                    return false;
                }// O
                else if (iAt->chemType.compare("O") == 0 && iAt->connAtoms.size() > 2) {
                    bool tFind = false;

                    for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                            iCo != iAt->connAtoms.end(); iCo++) {
                        if (tMol.atoms[*iCo].isMetal) {
                            tFind = true;
                            break;
                        }
                    }

                    if (!tFind) {
                        tErrInfo = "Reject Molecule: covalent bonded O atom has more"
                                " than 2 bonds!   Atom " + iAt->id
                                + " serial number " + IntToStr(iAt->seriNum)
                                + " has connections " + IntToStr((int) iAt->connAtoms.size());

                        return false;
                    }
                }// B
                else if (iAt->chemType.compare("B") == 0
                        && iAt->connAtoms.size() < 2) {
                    tErrInfo = "Reject molecule ! B Atom " + iAt->id
                            + " serial number " + IntToStr(iAt->seriNum)
                            + " has connections " + IntToStr((int) iAt->connAtoms.size());
                    return false;
                }// H 
                else if (iAt->chemType.compare("H") == 0 && iAt->connAtoms.size() > 1) {
                    bool findO = false;
                    for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                            iNB != iAt->connAtoms.end(); iNB++) {
                        if (tMol.atoms[*iNB].chemType == "O") {
                            findO = true;
                            break;
                        }
                    }
                    if ((findO && iAt->connAtoms.size() > 2)
                            || (!findO && iAt->connAtoms.size() > 1)) {
                        tErrInfo = "Reject the molecule! H Atom " + iAt->id
                                + " serial number " + IntToStr(iAt->seriNum)
                                + " has connections " + IntToStr((int) iAt->connAtoms.size());
                        return false;
                    }
                }// Halogen group
                else if (tPTab.elements[iAt->chemType]["group"] == 17
                        && iAt->connAtoms.size() > 1
                        && !connMetal(iAt->connAtoms, tMol.atoms)) {
                    tErrInfo = "Reject the molecule! Connection > 1 for Halogen Atom "
                            + iAt->id
                            + " serial number " + IntToStr(iAt->seriNum)
                            + " has connections " + IntToStr((int) iAt->connAtoms.size());

                    return false;
                }
            }
        }

        return true;

    }

    bool MolGenerator::validateMolecule(Molecule & tMol,
            PeriodicTable & tPTab,
            std::string & tErrInfo) {

        if (!checkAsmAtomsInMol(tMol, tErrInfo)) {
            return false;
        }
        if (!checkAtomOcp(tMol, tErrInfo)) {
            return false;
        }
        if (!validateAtomLinks(tMol, tPTab, tErrInfo)) {
            return false;
        }

        std::cout << "Bond validations: number of bonds in the molecule "
                << tMol.bonds.size() << std::endl;

        for (std::vector<BondDict>::iterator iBo = tMol.bonds.begin();
                iBo != tMol.bonds.end(); iBo++) {
            if (!validateBonds(iBo, tMol, tErrInfo, tPTab)) {
                return false;
            }
        }


        if (!validateBondValueSame(tMol, tErrInfo)) {
            return false;
        }


        // addition value difference check
        if (!validateBondValueDiff(tMol, tErrInfo)) {
            return false;
        }

        return true;

    }

    void MolGenerator::getAtomTypeOneMol(Molecule& tMol) {
        //std::cout << "Number of atoms in this molecule is "
        //        << tMol.atoms.size() << std::endl;

        CodClassify aCodSys(tMol.atoms);


        aCodSys.codAtomClassify2(2);


        tMol.atoms.clear();
        for (std::vector<AtomDict>::iterator iAt = aCodSys.allAtoms.begin();
                iAt != aCodSys.allAtoms.end(); iAt++) {
            tMol.atoms.push_back(*iAt);
        }

        /*      
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                 iAt !=tMol.atoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                std::cout << "Atom " << iAt->id << " has COD class id "
                                     << iAt->codClass << std::endl;
            }
        }
         */
    }

    void MolGenerator::getAtomTypeOneMolNew(Molecule& tMol) 
    {
        //std::cout << "Number of atoms in this molecule is "
        //        << tMol.atoms.size() << std::endl;

        CodClassify aCodSys(tMol.atoms);

        setAtomsBondingAndChiralCenter(aCodSys.allAtoms);


        aCodSys.codAtomClassifyNew3(2);

        /*
        for (std::vector<AtomDict>::iterator iAt=aCodSys.allAtoms.begin();
                iAt !=aCodSys.allAtoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                std::cout << "atom " << iAt->id << " is read-in atom " 
                          << std::endl;
            }
            else
            {
                std::cout << "atom " << iAt->id << " is symm-gen atom " 
                          << std::endl;
            }
        }
         */



        tMol.atoms.clear();
        for (std::vector<AtomDict>::iterator iAt = aCodSys.allAtoms.begin();
                iAt != aCodSys.allAtoms.end(); iAt++) 
        {
            tMol.atoms.push_back(*iAt);
        }

        /*
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                 iAt !=tMol.atoms.end(); iAt++)
        {
           
            std::cout << "Atom " << iAt->id << " has COD class id "
                                 << iAt->codClass << std::endl;
            if (iAt->inRings.size() !=0)
            {
                std::cout << "it is in " << iAt->inRings.size() 
                          << " ring(s) " << std::endl;
            }
                
        }
         */
        tMol.rings.clear();
        for (std::vector<RingDict>::iterator iR = aCodSys.allRingsV.begin();
                iR != aCodSys.allRingsV.end(); iR++) {
            tMol.rings.push_back(*iR);
        }

        // re-index atom's inRing idx 
        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            iAt->inRings.clear();
        }

        int idxR = 0;
        for (std::vector<RingDict>::iterator iR = tMol.rings.begin();
                iR != tMol.rings.end(); iR++) {
            std::vector<ID> tAtIds;
            for (std::vector<AtomDict>::iterator iAt = iR->atoms.begin();
                    iAt != iR->atoms.end(); iAt++) {
                tAtIds.push_back(iAt->id);
            }
            for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                    iAt != tMol.atoms.end(); iAt++) {
                if (std::find(tAtIds.begin(), tAtIds.end(), iAt->id)
                        != tAtIds.end()) {
                    iAt->inRings.push_back(idxR);
                }
            }

            idxR++;
        }



        // check 
        /*
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++)
        {
            if(iAt->inRings.size() !=0)
            {
                std::cout << "Atom " << iAt->id 
                          << " is in the following ring(s): " << std::endl;
                for (std::vector<int>::iterator idxR = iAt->inRings.begin();
                        idxR !=iAt->inRings.end(); idxR++)
                {
                    std::cout << tMol.rings[*idxR].rep << std::endl;
                }
            }
        }
         */


    }

    void MolGenerator::setAtomNFormalCharge(Molecule& tMol) {
        // check if the input cif file does not put formal charges in some N atoms
        // see 1000006 for example 
        for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++) {
            if (iAt->chemType.compare("N") == 0) {
                if (iAt->connAtoms.size() == 4 && iAt->formalCharge != 1) {
                    iAt->formalCharge = 1.0;
                }
            }
        }
    }

    void MolGenerator::getAtomTypeMols() {
        for (unsigned i = 0; i < allMolecules.size(); i++) {
            getAtomTypeOneMol(allMolecules[i]);
        }
    }

    void MolGenerator::getOverallBondAndAngles() {
        std::map<std::string, std::vector<REAL> > aBM, aAM;

        for (std::vector<Molecule>::iterator iMol = allMolecules.begin();
                iMol != allMolecules.end(); iMol++) {
            //std::cout << "Mol " << iMol->seriNum << " has " << iMol->bonds.size() 
            //          << " bonds " << std::endl;

            for (std::vector<BondDict>::iterator iBo = iMol->bonds.begin();
                    iBo != iMol->bonds.end(); iBo++) {
                std::vector<sortMap3> tVec;
                struct sortMap3 tMap;
                tMap.key = iMol->atoms[iBo->atomsIdx[0]].codClass;
                tMap.val = iMol->atoms[iBo->atomsIdx[0]].id;
                tVec.push_back(tMap);
                tMap.key = iMol->atoms[iBo->atomsIdx[1]].codClass;
                tMap.val = iMol->atoms[iBo->atomsIdx[1]].id;
                tVec.push_back(tMap);
                std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                std::string tKey = tVec[0].key + "_" + tVec[1].key;

                if (aBM.find(tKey) == aBM.end()) {
                    // std::cout << "new key " << tKey << std::endl
                    //          << "value is " << iBo->value << std::endl;
                    iBo->atoms.clear();
                    iBo->atomsCodClasses.clear();
                    for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                            iAt != tVec.end(); iAt++) {
                        iBo->atoms.push_back(iAt->val);
                        iBo->atomsCodClasses.push_back(iAt->key);
                    }
                    bonds.push_back(*iBo);
                    aBM[tKey].push_back(iBo->value);
                } else if (!inVectABS(aBM[tKey], iBo->value, 0.00001)) {
                    // std::cout << "new value " << iBo->value << std::endl
                    //          << "Key is " << tKey << std::endl;
                    iBo->atoms.clear();
                    iBo->atomsCodClasses.clear();
                    for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                            iAt != tVec.end(); iAt++) {
                        iBo->atoms.push_back(iAt->val);
                        iBo->atomsCodClasses.push_back(iAt->key);
                    }
                    bonds.push_back(*iBo);
                    aBM[tKey].push_back(iBo->value);
                }


            }
            //std::cout << "Total Number of unique bonds are " << bonds.size() << std::endl;
            //std::cout << "this molecule has " << iMol->angles.size() 
            //          << " angles " << std::endl;

            for (std::vector<AngleDict>::iterator iAn = iMol->angles.begin();
                    iAn != iMol->angles.end(); iAn++) {

                std::vector<sortMap3> tVec;
                struct sortMap3 tMap;
                tMap.key = iMol->atoms[iAn->atoms[1]].codClass;
                tMap.val = iMol->atoms[iAn->atoms[1]].id;
                tVec.push_back(tMap);
                tMap.key = iMol->atoms[iAn->atoms[2]].codClass;
                tMap.val = iMol->atoms[iAn->atoms[2]].id;
                tVec.push_back(tMap);
                std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                std::string tKey = iMol->atoms[iAn->atoms[0]].codClass
                        + "_" + tVec[0].key + "_" + tVec[1].key;
                if (aAM.find(tKey) == aAM.end()) {
                    iAn->atomChemTypes.clear();
                    iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].id);
                    iAn->atomsCodClasses.clear();
                    iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                    for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                            iAt != tVec.end(); iAt++) {
                        iAn->atomChemTypes.push_back(iAt->val);
                        iAn->atomsCodClasses.push_back(iAt->key);
                    }
                    angles.push_back(*iAn);
                    aAM[tKey].push_back(iAn->value);
                } else if (!inVectABS(aAM[tKey], iAn->value, 0.0001)) {
                    iAn->atomChemTypes.clear();
                    iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].id);
                    iAn->atomsCodClasses.clear();
                    iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                    for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                            iAt != tVec.end(); iAt++) {
                        iAn->atomChemTypes.push_back(iAt->val);
                        iAn->atomsCodClasses.push_back(iAt->key);
                    }
                    angles.push_back(*iAn);
                    aAM[tKey].push_back(iAn->value);
                }
            }
            // std::cout << "Total number of unique angles are " << angles.size() << std::endl;
        }

        std::cout << "Total Number of unique bonds are " << bonds.size() << std::endl;
        std::cout << "Total number of unique angles are " << angles.size() << std::endl;

    }

    void MolGenerator::getOverallBondAndAnglesNew() {
        std::map<std::string, std::vector<REAL> > aBM, aAM;

        for (std::vector<Molecule>::iterator iMol = allMolecules.begin();
                iMol != allMolecules.end(); iMol++) {
            //std::cout << "Mol " << iMol->seriNum << " has " << iMol->bonds.size()
            //          << " bonds " << std::endl;

            int nVB = 0, nPB = 0;
            for (std::vector<BondDict>::iterator iBo = iMol->allBonds.begin();
                    iBo != iMol->allBonds.end(); iBo++) {
                if (iMol->atoms[iBo->atomsIdx[0]].isInPreCell
                        || iMol->atoms[iBo->atomsIdx[1]].isInPreCell) {
                    nVB++;
                    std::vector<int> aIdxSet;
                    aIdxSet.push_back(iBo->atomsIdx[0]);
                    aIdxSet.push_back(iBo->atomsIdx[1]);
                    /*
                    std::cout << "Check a bond " << std::endl;
                    std::cout << "atom 1 "       << iMol->atoms[iBo->atomsIdx[0]].id
                              << " with type "   << iMol->atoms[iBo->atomsIdx[0]].codClass
                              << std::endl 
                              << "atom 2 " << iMol->atoms[iBo->atomsIdx[1]].id
                              << " with type " << iMol->atoms[iBo->atomsIdx[1]].codClass
                              << std::endl << "OrderN " << iBo->orderN << std::endl;
                     */
                    if (!connMetal2ndNB(aIdxSet, iMol->atoms)) {
                        //std::cout << "bond contains no metal up to the 2NB " << std::endl;
                        std::vector<sortMap3> tVec;
                        struct sortMap3 tMap;
                        tMap.key = iMol->atoms[iBo->atomsIdx[0]].codClass;
                        tMap.val = iMol->atoms[iBo->atomsIdx[0]].id;
                        tVec.push_back(tMap);
                        tMap.key = iMol->atoms[iBo->atomsIdx[1]].codClass;
                        tMap.val = iMol->atoms[iBo->atomsIdx[1]].id;
                        tVec.push_back(tMap);
                        std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                        std::string tKey = tVec[0].key + "_" + tVec[1].key;

                        if (aBM.find(tKey) == aBM.end()) {
                            //std::cout << "new key " << tKey << std::endl
                            //          << "value is " << iBo->value << std::endl;
                            iBo->atoms.clear();
                            iBo->atomsCodClasses.clear();
                            for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                                    iAt != tVec.end(); iAt++) {
                                iBo->atoms.push_back(iAt->val);
                                iBo->atomsCodClasses.push_back(iAt->key);
                            }
                            bonds.push_back(*iBo);
                            aBM[tKey].push_back(iBo->value);

                        }
                        else if (!inVectABS(aBM[tKey], iBo->value, 0.00001)) {
                            // std::cout << "new value " << iBo->value << std::endl
                            //          << "Key is " << tKey << std::endl;
                            iBo->atoms.clear();
                            iBo->atomsCodClasses.clear();
                            for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                                    iAt != tVec.end(); iAt++) {
                                iBo->atoms.push_back(iAt->val);
                                iBo->atomsCodClasses.push_back(iAt->key);
                            }
                            bonds.push_back(*iBo);
                            aBM[tKey].push_back(iBo->value);
                        }

                    } else {
                        nPB++;
                        bondsMetal2NB[iMol->seriNum].push_back(iBo->seriNum);
                    }

                }
            }
            if (nPB == nVB) {
                noContriMols.push_back(iMol->seriNum);
            }

            //std::cout << "Total Number of unique bonds are " << bonds.size() << std::endl;
            //std::cout << "this molecule has " << iMol->angles.size() 
            //          << " angles " << std::endl;

            for (std::vector<AngleDict>::iterator iAn = iMol->angles.begin();
                    iAn != iMol->angles.end(); iAn++) {
                if (iMol->atoms[iAn->atoms[0]].isInPreCell
                        || iMol->atoms[iAn->atoms[1]].isInPreCell
                        || iMol->atoms[iAn->atoms[2]].isInPreCell) {
                    std::vector<int> aIdxSet;
                    aIdxSet.push_back(iAn->atoms[0]);

                    if (!connMetal2ndNB(aIdxSet, iMol->atoms)) {
                        std::vector<sortMap3> tVec;
                        struct sortMap3 tMap;
                        tMap.key = iMol->atoms[iAn->atoms[1]].codClass;
                        tMap.val = iMol->atoms[iAn->atoms[1]].id;
                        tVec.push_back(tMap);
                        tMap.key = iMol->atoms[iAn->atoms[2]].codClass;
                        tMap.val = iMol->atoms[iAn->atoms[2]].id;
                        tVec.push_back(tMap);
                        std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                        std::string tKey = iMol->atoms[iAn->atoms[0]].codClass
                                + "_" + tVec[0].key + "_" + tVec[1].key;
                        if (aAM.find(tKey) == aAM.end()) {
                            iAn->atomChemTypes.clear();
                            iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].id);
                            iAn->atomsCodClasses.clear();
                            iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                            for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                                    iAt != tVec.end(); iAt++) {
                                iAn->atomChemTypes.push_back(iAt->val);
                                iAn->atomsCodClasses.push_back(iAt->key);
                            }
                            angles.push_back(*iAn);
                            aAM[tKey].push_back(iAn->value);
                        } else if (!inVectABS(aAM[tKey], iAn->value, 0.0001)) {
                            iAn->atomChemTypes.clear();
                            iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].id);
                            iAn->atomsCodClasses.clear();
                            iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                            for (std::vector<sortMap3>::iterator iAt = tVec.begin();
                                    iAt != tVec.end(); iAt++) {
                                iAn->atomChemTypes.push_back(iAt->val);
                                iAn->atomsCodClasses.push_back(iAt->key);
                            }
                            angles.push_back(*iAn);
                            aAM[tKey].push_back(iAn->value);
                        }
                    } else {
                        anglesMetal2NB[iMol->seriNum].push_back(iAn->seriNum);
                    }
                }
            }

            // std::cout << "Total number of unique angles are " << angles.size() << std::endl;
        }

        std::cout << "Total Number of unique bonds are " << bonds.size() << std::endl;
        /*
        for (std::vector<BondDict>::iterator iBo=bonds.begin();
                iBo !=bonds.end(); iBo++)
        {
            std::cout << "Order N is " << iBo->orderN << std::endl;
        }
         */
        std::cout << "Total number of unique angles are " << angles.size() << std::endl;
    }

    void MolGenerator::outTableMols(std::ofstream & tMolTabs,
            Molecule & tMol) {
        if (tMol.atoms.size() > 0) {
            tMol.calcSumExcessElecs();
            tMol.calcSumCharges();

            tMolTabs << "data_mol_" << IntToStr(tMol.seriNum)
                    << std::endl;
            tMolTabs << "loop_" << std::endl
                    << "is_infinite_molecule " << "\t";
            if (tMol.isInf) {
                tMolTabs << "yes";
            } else {
                tMolTabs << "no";
            }
            tMolTabs << std::endl;
            tMolTabs << "total_number_of_atoms\t"
                    << tMol.atoms.size() << std::endl;
            tMolTabs << "total_number_of_excess_electrons\t" << tMol.sumExcessElecs
                    << std::endl;
            tMolTabs << "total_number_of_charges\t" << tMol.sumCharges
                    << std::endl;

            tMolTabs << "loop_" << std::endl
                    << "_chem_comp_atom.serial_number" << std::endl
                    << "_chem_comp_atom.atom_id " << std::endl
                    << "_chem_comp_atom.element_symbol" << std::endl
                    << "_chem_comp_atom.formal_charge" << std::endl
                    << "_chem_comp_atom.hybr" << std::endl
                    << "_chem_comp_atom.nb1_nb2_sp" << std::endl
                    << "_chem_comp_atom.extraEls" << std::endl
                    << "_chem_comp_atom.nb1_nb2_extra_elec" << std::endl
                    << "_chem_comp_atom.cod_type" << std::endl;

            for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                    iAt != tMol.atoms.end(); iAt++) {

                tMolTabs << std::setw(6) << iAt->seriNum + 1
                        << std::setw(6) << iAt->id
                        << std::setw(4) << iAt->chemType
                        << std::setw(6) << iAt->formalCharge
                        << std::setw(10) << iAt->hybrid << "    "
                        << std::setw(iAt->codNB1NB2_SP.size() + 4) << iAt->codNB1NB2_SP << "    "
                        << std::setw(5) << iAt->excessElec << "    "
                        << std::setw(iAt->codNB1NB2_ExElec.size() + 4)
                        << iAt->codNB1NB2_ExElec << "    "
                        << iAt->codClass << std::endl;
            }
            tMolTabs << std::endl;
        } else {
            std::cout << "There is no atoms in the molecule" << std::endl;
        }

        if (tMol.allBonds.size() > 0) {
            std::cout << std::endl;
            tMolTabs << "loop_" << std::endl
                    << "_chem_comp_bond.bond_serial_number" << std::endl
                    << "_chem_comp_bond.atom1_serial_number" << std::endl
                    << "_chem_comp_bond.atom2_serial_number" << std::endl
                    << "_chem_comp_bond.atom1_id" << std::endl
                    << "_chem_comp_bond.atom2_id" << std::endl
                    << "_chem_comp_bond.value_dist" << std::endl
                    << "_chem_comp_bond.is_in_same_ring" << std::endl
                    << "_chem_comp_bond.order" << std::endl;
            int nBo = 1;
            for (std::vector<BondDict>::iterator iBo = tMol.allBonds.begin();
                    iBo != tMol.allBonds.end(); iBo++) {
                std::vector<int> aIdxSet;
                aIdxSet.push_back(iBo->atomsIdx[0]);
                aIdxSet.push_back(iBo->atomsIdx[1]);
                if (!connMetal2ndNB(aIdxSet, tMol.atoms)) {
                    std::string tStr;
                    if (iBo->isInSameRing) {
                        tStr = "Y";
                    } else {
                        tStr = "N";
                    }

                    tMolTabs << std::setw(6) << nBo
                            << std::setw(6) << iBo->atomsIdx[0] + 1
                            << std::setw(6) << iBo->atomsIdx[1] + 1
                            << std::setw(4) << tMol.atoms[iBo->atomsIdx[0]].chemType
                            << std::setw(4) << tMol.atoms[iBo->atomsIdx[1]].chemType
                            << std::setw(10) << iBo->value
                            << std::setw(6) << tStr
                            << std::setw(10) << iBo->orderN << std::endl;
                    nBo++;
                }

            }

            std::cout << std::endl;
        } else {
            std::cout << "There is no bonds in the molecule" << std::endl;
        }

        if (tMol.angles.size() > 0) {
            tMolTabs << "loop_" << std::endl
                    << "_chem_comp_angle.angle_serial_number" << std::endl
                    << "_chem_comp_angle.atom1_serial_number" << std::endl
                    << "_chem_comp_angle.atom2_serial_number" << std::endl
                    << "_chem_comp_angle.atom3_serial_number" << std::endl
                    << "_chem_comp_angle.atom1_id" << std::endl
                    << "_chem_comp_angle.atom2_id" << std::endl
                    << "_chem_comp_angle.atom3_id" << std::endl
                    << "_chem_comp_angle.value_angle" << std::endl
                    << "_chem_comp_angle.ring_size_of_angle" << std::endl;
            int nAn = 1;
            for (std::vector<AngleDict>::iterator iAn = tMol.angles.begin();
                    iAn != tMol.angles.end(); iAn++) {


                tMolTabs << std::setw(6) << nAn
                        << std::setw(6) << iAn->atoms[0] + 1
                        << std::setw(6) << iAn->atoms[1] + 1
                        << std::setw(6) << iAn->atoms[2] + 1
                        << std::setw(4) << tMol.atoms[iAn->atoms[0]].chemType
                        << std::setw(4) << tMol.atoms[iAn->atoms[1]].chemType
                        << std::setw(4) << tMol.atoms[iAn->atoms[2]].chemType
                        << std::setw(10) << iAn->value * PID180
                        << std::setw(6) << iAn->isInSameRing << std::endl;
                nAn++;
            }

            std::cout << std::endl;
        }

        /*
        if (tMol.atoms.size() >0)
        {
            int sumEx=sumExElectrons(tMol.atoms);
            
            tMolTabs << "loop_" << std::endl
                     << "_chem_comp_exElectrons.sum" << std::endl;
            tMolTabs << std::setw(6)  << sumEx << std::endl << std::endl;
        }
         */


    }

    void MolGenerator::outTableBAndA(FileName tBAndAFName) {
        // std::cout << tBAndAFName << std::endl;


        std::ofstream aBAndAF(tBAndAFName);
        if (aBAndAF.is_open()) {
            if (bonds.size() > 0) {
                aBAndAF << "Atom1_COD_type        Atom2_COD_type      Atom1_id       Atom2_id       "
                        << "Atom1_sp       Atom2_sp    atom1_NB1NB2_sp    atom2_NB1NB2_sp   "
                        << "atom1_NB2_extraEls       atom2_NB2_extraEls    "
                        << "Bond_length       BondInRing    "
                        << "BondOrder         Idx_of_Mol" << std::endl;
                for (std::vector<BondDict>::iterator iBo = bonds.begin();
                        iBo != bonds.end(); iBo++) {

                    std::string tStr;
                    if (iBo->isInSameRing) {
                        tStr = "Y";
                    } else {
                        tStr = "N";
                    }

                    if (iBo->atomsCodClasses.size() == 2) {
                        aBAndAF << iBo->atomsCodClasses[0] << "\t"
                                << iBo->atomsCodClasses[1] << "\t"
                                << iBo->atoms[0] << "\t" << iBo->atoms[1] << "\t"
                                << iBo->atomSPs[iBo->atoms[0]] << "\t"
                                << iBo->atomSPs[iBo->atoms[1]] << "\t"
                                << iBo->atomNB1NB2SPs[iBo->atoms[0]] << "\t"
                                << iBo->atomNB1NB2SPs[iBo->atoms[1]] << "\t"
                                << iBo->atomNB2ExtraEls[iBo->atoms[0]] << "\t"
                                << iBo->atomNB2ExtraEls[iBo->atoms[1]] << "\t"
                                << iBo->value << "\t"
                                << tStr << "\t"
                                << iBo->orderN << "\t"
                                << iBo->molIdx << std::endl;
                    }
                }

                aBAndAF << std::endl;
            }


            if (angles.size() != 0) {
                aBAndAF << "Center_Atom_COD_type        Atom1_COD_type          Atom2_COD_type     "
                        << "Center_Atom_id      Atom1_id     Atom2_id     "
                        << "Center_Atom_sp     Atom1_sp    Atom2_sp    "
                        << "Center_Atom_NB1NB2_sp   Atom1_NB1NB2_sp    Atom2_NB1NB2_sp"
                        << "Angle_value      RingSize"
                        << std::endl;
                for (std::vector<AngleDict>::iterator iAn = angles.begin();
                        iAn != angles.end(); iAn++) {
                    aBAndAF << iAn->atomsCodClasses[0] << "\t"
                            << iAn->atomsCodClasses[1] << "\t"
                            << iAn->atomsCodClasses[2] << "\t"
                            << iAn->atomChemTypes[0] << "\t"
                            << iAn->atomChemTypes[1] << "\t"
                            << iAn->atomChemTypes[2] << "\t"
                            << iAn->atomsSPStats[iAn->atomChemTypes[0]] << "\t"
                            << iAn->atomsSPStats[iAn->atomChemTypes[1]] << "\t"
                            << iAn->atomsSPStats[iAn->atomChemTypes[2]] << "\t"
                            << iAn->atomsNB1NB2SPStats[iAn->atomChemTypes[0]] << "\t"
                            << iAn->atomsNB1NB2SPStats[iAn->atomChemTypes[1]] << "\t"
                            << iAn->atomsNB1NB2SPStats[iAn->atomChemTypes[2]] << "\t"
                            << iAn->value * PID180 << "\t"
                            << iAn->isInSameRing << std::endl;
                }
                aBAndAF << std::endl;
            }
            aBAndAF.close();
        }

    }

    void MolGenerator::setTableSpAndChirals(Molecule & tMol,
            std::map<std::string, std::map<std::string,
            std::map<std::string, std::map<std::string,
            REAL> > > >& tAtomSpAndChMap) {
        if (tMol.atoms.size() > 0) {
            for (std::vector<AtomDict>::iterator iAt = tMol.atoms.begin();
                    iAt != tMol.atoms.end(); iAt++) {
                StrUpper(iAt->hybrid);
                if ((iAt->chemType.find("N") != std::string::npos
                        || iAt->chemType.find("B") != std::string::npos)
                        && (iAt->hybrid.find("SP2") != std::string::npos
                        || iAt->hybrid.find("SP3") != std::string::npos)
                        && iAt->isInPreCell) {
                    bool lCal = false;
                    if (tAtomSpAndChMap.find(iAt->chemType) == tAtomSpAndChMap.end()) {
                        lCal = true;
                    }
                    if (tAtomSpAndChMap[iAt->chemType].find(iAt->id)
                            == tAtomSpAndChMap[iAt->chemType].end()) {
                        lCal = true;
                    } else if (tAtomSpAndChMap[iAt->chemType][iAt->id].find(iAt->hybrid)
                            == tAtomSpAndChMap[iAt->chemType][iAt->id].end()) {
                        lCal = true;
                    }

                    if (lCal && iAt->connAtoms.size() == 3) {
                        REAL nH = 0.0;
                        for (std::vector<int>::iterator iCo = iAt->connAtoms.begin();
                                iCo != iAt->connAtoms.end(); iCo++) {
                            if (tMol.atoms[*iCo].chemType.find("H") != std::string::npos) {
                                nH++;
                            }
                        }

                        tAtomSpAndChMap[iAt->chemType][iAt->id]
                                [iAt->hybrid]["numH"] = nH;

                        std::vector<REAL> tVect1, tVect2, tVect3;
                        std::cout << "Atom " << iAt->id << " is connected to " << std::endl
                                << "Atom " << tMol.atoms[iAt->connAtoms[0]].id << std::endl
                                << "Atom " << tMol.atoms[iAt->connAtoms[1]].id << std::endl
                                << "Atom " << tMol.atoms[iAt->connAtoms[2]].id << std::endl;

                        diffVects(iAt->coords, tMol.atoms[iAt->connAtoms[0]].coords, tVect1);
                        std::cout << "bond 1 " << lengthV(tVect1) << std::endl;
                        diffVects(iAt->coords, tMol.atoms[iAt->connAtoms[1]].coords, tVect2);
                        std::cout << "bond 2 " << lengthV(tVect2) << std::endl;
                        diffVects(iAt->coords, tMol.atoms[iAt->connAtoms[2]].coords, tVect3);
                        std::cout << "bond 3 " << lengthV(tVect3) << std::endl;
                        if (tVect1.size() == tVect2.size()
                                && tVect1.size() == tVect3.size()
                                && tVect1.size() != 0) {
                            // REAL tCV=calNormalizedChiralVol(tVect1, tVect2, tVect3);
                            //std::string sCri = RealToStr(aCri);
                            //std::string aSup;

                            std::vector<REAL> tV12, tV13, tV23;
                            crossP2V(tVect1, tVect2, tV12);
                            crossP2V(tVect1, tVect3, tV13);
                            crossP2V(tVect2, tVect3, tV23);

                            REAL ang1 = getAngle2V(tV12, tV13) * PID180;
                            transAng(ang1);
                            REAL ang2 = getAngle2V(tV13, tV23) * PID180;
                            transAng(ang2);
                            REAL ang3 = getAngle2V(tV12, tV23) * PID180;
                            transAng(ang3);

                            REAL ang = (ang1 + ang2 + ang3) / 3.0;

                            /*
                            REAL ang  = ang1;
                            if (ang2 > ang)
                            {
                                ang = ang2;
                            }
                            if (ang3 > ang)
                            {
                                ang = ang3;
                            }
                             */

                            tAtomSpAndChMap[iAt->chemType][iAt->id]
                                    [iAt->hybrid]["ang"] = ang;
                        }

                    }
                }
            }
        }
    }

    void MolGenerator::contMetal2NB(int & tNB, int & tNA) {
        for (std::map<int, std::vector<int> >::iterator iM = bondsMetal2NB.begin();
                iM != bondsMetal2NB.end(); iM++) {
            tNB += iM->second.size();
        }
        for (std::map<int, std::vector<int> >::iterator iM = anglesMetal2NB.begin();
                iM != anglesMetal2NB.end(); iM++) {
            tNA += iM->second.size();
        }

    }

    void MolGenerator::outTables(FileName tOutName,
            std::vector<Molecule> & tFinMols,
            std::vector<Molecule> & tInfMols) {
        int aNB = 0, aNA = 0;
        contMetal2NB(aNB, aNA);

        std::cout << "Bonds with metals at the secondary NB " << aNB << std::endl;
        std::cout << "Angles with metals at the secondary NB " << aNB << std::endl;
        Name aFName(tOutName);
        // std::cout << "Output root is " << aFName << std::endl;
        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');
        Name rootFName;


        for (unsigned jF = 0; jF < nameComps.size(); jF++) {
            rootFName.append(nameComps[jF]);
        }
        if (rootFName.size() == 0) {
            rootFName.append("Current");
        }
        unsigned nAll = validedMolMsg.size() + errMolMsg.size();
        if (nAll != 0) {

            Name msgFName(rootFName);
            msgFName.append("_msg.txt");
            std::ofstream msg(msgFName.c_str());
            msg << "Total number of molecules : " << nAll << std::endl;
            msg << "Total number of validated molecules "
                    << validedMolMsg.size() << std::endl;
            msg << "Number of molecules without contributions "
                    << noContriMols.size() << std::endl;
            msg << "Number of validated bonds " << bonds.size() << std::endl;
            msg << "Number of excluded bonds because of metal NB "
                    << aNB << std::endl;
            msg << "Number of excluded angles because of metal NB "
                    << aNA << std::endl;
            for (std::map<int, std::string>::iterator iVM = validedMolMsg.begin();
                    iVM != validedMolMsg.end(); iVM++) {
                msg << "VM Begin: " << std::endl
                        << iVM->second
                        << "VM End:" << std::endl;
                //std::cout << "For molecule " << iVM->first 
                //          << ", the message is " << std::endl
                //          << iVM->second << std::endl;
            }

            msg << "Total number of deleted molecules "
                    << errMolMsg.size() << std::endl;
            for (std::map<int, std::string>::iterator iEM = errMolMsg.begin();
                    iEM != errMolMsg.end(); iEM++) {
                msg << "DM Begin: " << std::endl
                        << "molecule " << iEM->first
                        << ", the message is " << std::endl
                        << iEM->second << "DM End:" << std::endl;


                //std::cout << "For molecule " << iEM->first 
                //          << ", the message is " << std::endl
                //          << iEM->second << std::endl;
            }
            msg.close();
        }

        if (bonds.size() != 0) {
            Name bondAndAngleFName(rootFName);
            bondAndAngleFName.append("_unique_bond_and_angles.txt");
            outTableBAndA(bondAndAngleFName.c_str());
        }

        if (allMolecules.size() != 0) {
            Name allMolsFName(rootFName);
            allMolsFName.append("_all_mols.txt");
            std::ofstream aMolTable(allMolsFName.c_str());
            if (aMolTable.is_open()) {
                outMolsInfo(aMolTable, tFinMols, tInfMols);
            }
            else {
                std::cout << allMolsFName << " can not be open for writing "
                        << std::endl;
            }

            std::map<std::string, std::map< std::string, std::map<std::string,
                    std::map<std::string, REAL> > > >atomSpAndChMap;
            Name atomSpAndChMapFName(rootFName);
            atomSpAndChMapFName.append("_atom_sp_ang.list");
            std::ofstream atomSpAndChMapF(atomSpAndChMapFName.c_str());
            if (atomSpAndChMapF.is_open()) {
                for (unsigned iMol = 0; iMol < allMolecules.size(); iMol++) {
                    setTableSpAndChirals(allMolecules[iMol], atomSpAndChMap);
                }

                for (std::map<std::string, std::map<std::string, std::map<std::string,
                        std::map<std::string, REAL> > > >::iterator
                        iMap = atomSpAndChMap.begin(); iMap != atomSpAndChMap.end();
                        iMap++) {
                    for (std::map<std::string, std::map<std::string,
                            std::map<std::string, REAL> > >::iterator
                            iSC = iMap->second.begin(); iSC != iMap->second.end();
                            iSC++) {
                        for (std::map<std::string, std::map<std::string,
                                REAL> >::iterator iP = iSC->second.begin();
                                iP != iSC->second.end(); iP++) {
                            atomSpAndChMapF << iMap->first << "\t"
                                    << iSC->first << "\t"
                                    << iP->first << "\t"
                                    << iP->second["numH"] << "\t"
                                    << iP->second["ang"] << std::endl;
                        }
                    }
                }
                atomSpAndChMapF.close();
            }
            else {
                std::cout << allMolsFName << " can not be open for writing "
                        << std::endl;
            }
        }
    }

    void MolGenerator::outMsg(FileName tOutName) {
        Name aFName(tOutName);

        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');
        Name rootFName;

        for (unsigned jF = 0; jF < nameComps.size(); jF++) {
            rootFName.append(nameComps[jF]);
        }
        if (rootFName.size() == 0) {
            rootFName.append("Current");
        }

        Name msgFName(rootFName);
        msgFName.append("_msg.txt");
        std::ofstream msg(msgFName.c_str());

        if (allMsg.size()) {
            for (std::vector<std::string>::iterator iS = allMsg.begin();
                    iS != allMsg.end(); iS++) {
                msg << *iS;
            }
        }
        msg.close();
    }

    void MolGenerator::getOutFileRoot(FileName tOutName, Name & tRootName) {
        Name aFName(tOutName);

        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');

        for (unsigned jF = 0; jF < nameComps.size() - 1; jF++) {
            tRootName.append(nameComps[jF]);
        }
        if (tRootName.size() == 0) {
            tRootName.append("Current");
        }
    }

    void MolGenerator::outMolsInfo(std::ofstream & tMolTabs,
            std::vector<Molecule>& tFinMols,
            std::vector<Molecule>& tInfMols) {


        if (tMolTabs.is_open()) {

            tMolTabs << "loop_" << std::endl
                    << "num_finite_mols\t" << tFinMols.size() << std::endl
                    << "num_infinite_mols\t" << tInfMols.size() << std::endl;


            if (tFinMols.size() > 0) {
                for (unsigned i = 0; i < tFinMols.size(); i++) {
                    outTableMols(tMolTabs, tFinMols[i]);
                }
            }

            if (tInfMols.size() > 0) {
                for (unsigned i = 0; i < tInfMols.size(); i++) {
                    outTableMols(tMolTabs, tInfMols[i]);
                }
            }
            tMolTabs.close();
        }
    }

    void MolGenerator::checkInfMols(std::vector<Molecule> & aSetInfMols,
            std::vector<Molecule> & aSetFinMols) {
        std::vector<Molecule> aSetUniqMols;

        for (unsigned iMol = 0; iMol < allMolecules.size(); iMol++) {
            if (checEquiMoles(aSetUniqMols, allMolecules[iMol])) {
                aSetUniqMols.push_back(allMolecules[iMol]);

            }
        }

        if (aSetUniqMols.size() != 0) {
            for (std::vector<Molecule>::iterator iMol = aSetUniqMols.begin();
                    iMol != aSetUniqMols.end(); iMol++) {
                if (!checkOneMolInf(iMol)) {
                    aSetFinMols.push_back(*iMol);
                } else {
                    iMol->isInf = true;
                    aSetInfMols.push_back(*iMol);
                }
            }
        }

    }

    bool MolGenerator::checkOneMolInf(std::vector<Molecule>::iterator tMol) {
        bool aInfMol = false;
        // Two criteria are used. Both should be satisfied
        if (allCryst.size() != 0) {
            if (allCryst[0].itsSpaceGroup != NullPoint) {
                int maxMolAtms = initAtoms.size() *
                        allCryst[0].itsSpaceGroup->sgOp.size();
                if (tMol->atoms.size() > maxMolAtms) {
                    aInfMol = true;
                }
            }
        }

        if (!aInfMol) {
            for (unsigned i = 0; i < tMol->atoms.size(); i++) {
                // Check if a molecule containing both an atom and its 
                // translation copy.  
                // 3 conditions:
                // 1) One atom's seriNum equals to the other's fromOrig
                // 2) They belong to the same symm operator
                // 3) They are from different unit cell (sId different)
                for (unsigned j = i + 1; j < tMol->atoms.size(); j++) {
                    if (tMol->atoms[i].seriNum == tMol->atoms[j].fromOrig
                            && tMol->atoms[i].symmOp == tMol->atoms[j].symmOp
                            && tMol->atoms[i].sId != tMol->atoms[i].sId) {
                        aInfMol = true;
                        break;
                    }
                }
            }
        }
        return aInfMol;
    }

    void MolGenerator::buildMetalClusters(
                              std::vector<CrystInfo>::iterator tCryst)
    {
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) 
        {
            if (iAt->isInPreCell && iAt->isMetal && checkNBAtomOccp(iAt)) 
            {
                bool lMC=true;
                metalCluster aMC;
                
                aMC.metSeril = iAt->seriNum;
                
                std::cout << "A Metal cluster is found at atom " 
                          << iAt->id << " of serial number " 
                          << iAt->seriNum << std::endl;
                
                std::cout << "It connects to " << iAt->connAtoms.size()
                          << " atoms " << std::endl;
                
                for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                        iNB != iAt->connAtoms.end(); iNB++) 
                {
                    if (checkNBAtomOccp(allAtoms[*iNB])
                        && allAtoms[*iNB].chemType.compare("H") !=0
                        && allAtoms[*iNB].chemType.compare("D") !=0)
                    {
                        aMC.ligandSerilSet.push_back(*iNB);
                        
                        for (std::vector<int>::iterator 
                             iNB2=allAtoms[*iNB].connAtoms.begin();
                             iNB2 !=allAtoms[*iNB].connAtoms.end(); iNB2++)
                        {
                            if(*iNB2 !=iAt->seriNum)
                            {
                                aMC.ligandNBs[*iNB].push_back(*iNB2);
                            }
                        }
                    }
                    else
                    {
                        lMC = false;
                        break;
                    }
                   
                }
                
                if (lMC)
                {
                    aMC.setMetClusterFormu(allAtoms);
                    aMC.buildBondAndAngleMap(allAtoms, tCryst);
                    allMetalClusters.push_back(aMC);
                }
            }
        }
    }
    
    void MolGenerator::buildMetalAtomCoordMap(
            std::vector<CrystInfo>::iterator tCryst) 
    {
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++) 
        {
            if (iAt->isInPreCell && iAt->isMetal && checkNBAtomOccp(iAt)) 
            {
                std::cout << "Atom " << iAt->id << " connects to "
                        << iAt->connAtoms.size()
                        << " atoms " << std::endl;
                for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                        iNB != iAt->connAtoms.end(); iNB++) 
                {
                    if (allAtoms[*iNB].chemType.compare("H") !=0
                        && allAtoms[*iNB].chemType.compare("D") !=0)
                    {
                        connMapMetalAtoms[iAt->seriNum].push_back(*iNB);
                        std::cout << "neighbor atom " << allAtoms[*iNB].id
                                  << " of " << allAtoms[*iNB].seriNum << std::endl;
                        REAL rD = getBondLenFromFracCoords(iAt->fracCoords,
                                                 allAtoms[*iNB].fracCoords,
                                    tCryst->itsCell->a, tCryst->itsCell->b,
                                    tCryst->itsCell->c, tCryst->itsCell->alpha,
                                    tCryst->itsCell->beta, tCryst->itsCell->gamma);
                        std::cout << "Distance " << rD << std::endl;
                        metalRelatedBonds[iAt->seriNum][*iNB] = rD;
                        for (std::vector<int>::iterator iNB2 = iAt->connAtoms.begin();
                                iNB2 != iAt->connAtoms.end(); iNB2++) 
                        {
                            if (*iNB2 > *iNB && allAtoms[*iNB2].chemType.compare("H") !=0
                                 && allAtoms[*iNB2].chemType.compare("D") !=0) 
                            {
                                std::cout << " the other NB " << allAtoms[*iNB2].id
                                          << " of serial " << allAtoms[*iNB2].seriNum
                                          << std::endl;

                                REAL rA = getAngleValueFromFracCoords(*iAt,
                                          allAtoms[*iNB], allAtoms[*iNB2],
                                          tCryst->itsCell->a, tCryst->itsCell->b,
                                          tCryst->itsCell->c, tCryst->itsCell->alpha,
                                          tCryst->itsCell->beta, tCryst->itsCell->gamma);
                                std::cout << "Angle " << rA * PID180 << std::endl;
                                metalRelatedAngles[iAt->seriNum][*iNB][*iNB2] = rA*PID180;
                            }
                        }
                    }
                }
            }
        }
    }

    bool MolGenerator::checkNBAtomOccp(std::vector<AtomDict>::iterator tAtm) 
    {
        // For metal atoms, need to be sure upto the second NB atoms have 
        // full occp
        bool aRe = true;
        for (std::vector<int>::iterator lCo = tAtm->connAtoms.begin();
                lCo != tAtm->connAtoms.end(); lCo++) 
        {
            if (allAtoms[*lCo].ocp <= 0.99) 
            {
                aRe = false;
                break;
            }
            
            
        }

        return aRe;

    }

    bool MolGenerator::checkNBAtomOccp(AtomDict & tAtm) 
    {
        // For metal atoms, need to be sure upto the second NB atoms have 
        // full occp
        bool aRe = true;
        for (std::vector<int>::iterator lCo = tAtm.connAtoms.begin();
                lCo != tAtm.connAtoms.end(); lCo++) 
        {
            if (allAtoms[*lCo].ocp <= 0.99) 
            {
                aRe = false;
                break;
            }
        }

        return aRe;

    }
    
    int MolGenerator::getNumOrgNB(std::vector<AtomDict>& tAtoms,
            int tIdx, std::vector<std::string> & tOrgTab) 
    {
        int nOrg = 0;
        for (std::vector<int>::iterator iCo = tAtoms[tIdx].connAtoms.begin();
                iCo != tAtoms[tIdx].connAtoms.end(); iCo++) {
            if (isOrganc(tOrgTab, tAtoms[*iCo].chemType) 
                && tAtoms[*iCo].chemType.compare("H")!=0
                && tAtoms[*iCo].chemType.compare("D")!=0) 
            {
                nOrg++;
            }
        }

        return nOrg;
    }

    void MolGenerator::outMetalAtomCoordInfo(FileName tOutName) 
    {
        std::vector<std::string> aOrgTab;
        initOrgTable(aOrgTab);

        Name aFName(tOutName);
        // std::cout << "Output root is " << aFName << std::endl;
        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');
        Name rootFName;


        for (unsigned jF = 0; jF < nameComps.size(); jF++) {
            rootFName.append(nameComps[jF]);
        }

        if (rootFName.size() == 0) {
            rootFName.append("Current");
        }

        if (allMsg.size()) {

        }

        if (connMapMetalAtoms.size() != 0) {
            std::vector<REAL> uniqBV, uniqAV;

            Name bondAndAngleFName(rootFName);
            bondAndAngleFName.append("_unique_bond_and_angles.txt");
            std::ofstream aBAndAF(bondAndAngleFName.c_str());
            if (aBAndAF.is_open()) {
                aBAndAF << "Element1\tElement2\t"
                        << "AtomName1\tAtomName2\t"
                        << "CoordinationNum1\tCoordinationNum1OrganicOnly\t"
                        << "CoordinationNumber2\t"
                        << "BondLength" << std::endl;
                for (std::map<int, std::vector<int> >::iterator
                    iM = connMapMetalAtoms.begin();
                        iM != connMapMetalAtoms.end(); iM++) {
                    const int nOrgConns = getNumOrgNB(allAtoms, iM->first, aOrgTab);

                    for (std::vector<int>::iterator iP = iM->second.begin();
                            iP != iM->second.end(); iP++) {
                        bool lDo = false;
                        

                        if (!inVectABS(uniqBV, metalRelatedBonds[iM->first][*iP],
                                0.00001)) {
                            lDo = true;
                        }
                        
                        //
                        bool l2NB
                             =checkNBAtomOccp(allAtoms[allAtoms[*iP].fromOrig]);
                        
                        if (lDo && l2NB) 
                        {
                            std::cout << "Bond : atom " << allAtoms[iM->first].id
                                      << " and " << allAtoms[*iP].id
                                      << " value " << metalRelatedBonds[iM->first][*iP]
                                      << std::endl;
                            
                            int nNBConn = 
                             allAtoms[allAtoms[*iP].fromOrig].connAtoms.size();
                            aBAndAF << allAtoms[iM->first].chemType << "\t"
                                    << allAtoms[*iP].chemType << "\t"
                                    << allAtoms[iM->first].id << "\t"
                                    << allAtoms[*iP].id << "\t"
                                    << iM->second.size() << "\t"
                                    << nOrgConns << "\t"
                                    << nNBConn << "\t"
                                    << metalRelatedBonds[iM->first][*iP] << std::endl;
                            uniqBV.push_back(metalRelatedBonds[iM->first][*iP]);
                        }
                    }
                }

                if (metalRelatedAngles.size() > 0) {
                    aBAndAF << "ElementCenter\tElement1\tElement2\t"
                            << "AtomNameCenter\tAtomName1\t"
                            << "AtomName2\tCoordinationNumber\t"
                            << "BondAngle" << std::endl;

                    for (std::map<int, std::map<int, std::map<int, REAL> > >::iterator
                        iMA = metalRelatedAngles.begin();
                            iMA != metalRelatedAngles.end(); iMA++) {
                        for (std::map<int, std::map<int, REAL> >::iterator
                            iNB1 = iMA->second.begin(); iNB1 != iMA->second.end();
                                iNB1++) {
                            for (std::map<int, REAL>::iterator
                                iNB2 = iNB1->second.begin();
                                    iNB2 != iNB1->second.end(); iNB2++) {
                                if (!inVectABS(uniqAV,
                                        metalRelatedAngles[iMA->first][iNB1->first][iNB2->first],
                                        0.00001)) {
                                    aBAndAF << allAtoms[iMA->first].chemType << "\t"
                                            << allAtoms[iNB1->first].chemType << "\t"
                                            << allAtoms[iNB2->first].chemType << "\t"
                                            << allAtoms[iMA->first].id << "\t"
                                            << allAtoms[iNB1->first].id << "\t"
                                            << allAtoms[iNB2->first].id << "\t"
                                            << connMapMetalAtoms[iMA->first].size() << "\t"
                                            << metalRelatedAngles[iMA->first]
                                            [iNB1->first][iNB2->first]
                                            << std::endl;

                                    uniqAV.push_back(metalRelatedAngles[iMA->first]
                                            [iNB1->first][iNB2->first]);
                                }
                            }
                        }
                    }
                }

                aBAndAF.close();
            }
        }

    }
    
    void MolGenerator::outMetalClusterInfo(FileName tOutName)
    {
        std::vector<std::string> aOrgTab;
        initOrgTable(aOrgTab);

        Name aFName(tOutName);
        // std::cout << "Output root is " << aFName << std::endl;
        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');
        Name rootFName;

        for (unsigned jF = 0; jF < nameComps.size(); jF++) {
            rootFName.append(nameComps[jF]);
        }

        if (rootFName.size() == 0) {
            rootFName.append("Current");
        }

        if (allMsg.size()) 
        {
        }
        
        if (allMetalClusters.size() !=0)
        {
            Name bondAndAngleFName(rootFName);
            bondAndAngleFName.append("_unique_bond_and_angles.txt");
            std::ofstream aBAndAF(bondAndAngleFName.c_str());
            if (aBAndAF.is_open()) 
            {
                aBAndAF << "MetalElement\tLigandElement2\t"
                        << "AtomName1\tAtomName2\t"
                        << "CoordinationNum1\tCoordinationNum1OrganicOnly\t"
                        << "CoordinationNumber2\t"
                        << "BondLength" << std::endl;
                
            }
        }
        
    }
    

}
