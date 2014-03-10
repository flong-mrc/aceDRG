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

namespace LIBMOL
{
    MolGenerator::MolGenerator():myNBDepth(1)
    {
    }
    
    
    MolGenerator::MolGenerator(const GenCifFile& tCifObj, int tNBDepth)
    {
        
        myNBDepth = tNBDepth;
        
        // std::cout << "myNBDepth " << myNBDepth << std::endl;
        
        std::cout << "Number of crystals " << tCifObj.allCryst.size() << std::endl;
        
        
        for (std::vector<CrystInfo>::const_iterator iC=tCifObj.allCryst.begin();
                iC!=tCifObj.allCryst.end(); iC++)
        {
            allCryst.push_back(*iC);
        }
        
        for (std::vector<AtomDict>::const_iterator iA=tCifObj.allAtoms.begin();
                iA!=tCifObj.allAtoms.end(); iA++)
        {
            initAtoms.push_back(*iA);
        }
        std::cout << "number of initial atoms " << initAtoms.size() << std::endl;
        
    }
    
    MolGenerator::MolGenerator(const DictCifFile & tCifObj, int tNBDepth)
    {
        
        myNBDepth = tNBDepth;
        
        for (std::vector<AtomDict>::const_iterator iA=tCifObj.allAtoms.begin();
                iA!=tCifObj.allAtoms.end(); iA++)
        {
            initAtoms.push_back(*iA);
        }
    }
    
    
    MolGenerator::~MolGenerator()
    {
    }
    
    
    void MolGenerator::execute(FileName tOutName)
    {
        
        CCP4DictParas tCCP4EnergParas;
        ccp4DictParas.push_back(tCCP4EnergParas);
        PeriodicTable aPTable;
        
        
        if (initAtoms.size() !=0)
        {
            for (std::vector<CrystInfo>::iterator iCryst=allCryst.begin();
                    iCryst !=allCryst.end(); iCryst++)
            {
                for (std::vector<AtomDict>::iterator iA=initAtoms.begin();
                    iA !=initAtoms.end(); iA++)
                {
                    FractToOrtho(iA->fracCoords, iA->coords, iCryst->itsCell->a,
                                iCryst->itsCell->b, iCryst->itsCell->c, iCryst->itsCell->alpha,
                                iCryst->itsCell->beta, iCryst->itsCell->gamma);
                    allAtoms.push_back(*iA);
                    refAtoms.push_back(*iA);
                }
                //outPDB("initAtoms.pdb", "UNL", initAtoms);
                   
                buildRefAtoms(iCryst);
                
                std::cout << "Number of atoms read from the input file "
                          << initAtoms.size() << std::endl; 
                std:: cout << "number of ref atoms associated with initial atoms " 
                           << refAtoms.size() << std::endl;
                symmAtomGen(iCryst, aPTable);
                
                std::cout << "number of atoms in a unit cell " 
                          << allAtoms.size() << std::endl;
                std::cout << "Number of atoms in refAtoms " 
                          << refAtoms.size() << std::endl;
                
                getUniqueBonds(aPTable);
                
                std::cout << "The following are the bonds we have found :" << std::endl;
                for (std::vector<BondDict>::iterator iB=bonds.begin();
                        iB !=bonds.end(); iB++)
                {
                    std::cout << "Bond between atom " << refAtoms[iB->atomsIdx[0]].id 
                              << " and " << refAtoms[iB->atomsIdx[1]].id << std::endl 
                              << "Its bond length is " << iB->value << std::endl;
                }
                
                getUniqAngles();
                
                CodClassify aCodSys(allAtoms);
                aCodSys.codAtomClassify(2);
                
                allAtoms.clear();
                for (std::vector<AtomDict>::iterator iAt=aCodSys.allAtoms.begin();
                        iAt!=aCodSys.allAtoms.end(); iAt++)
                {
                    allAtoms.push_back(*iAt);
                }
                
                getMolByEqClassInCell();
                
                if (allMolecules.size() > 0)
                {
                    for (unsigned  i=0; i < allMolecules.size(); i++)
                    {
                        if (allMolecules[i].atoms.size() >0 &&
                            allMolecules[i].bonds.size() >0 &&
                            allMolecules[i].angles.size() >0)
                        {
                            Name aFName(tOutName);
                            std::vector<std::string> nameComps;
                            StrTokenize(aFName, nameComps, '.');
                            Name bFName;
                            for (unsigned jF=0; jF < nameComps.size()-1; jF++)
                            {
                                bFName.append(nameComps[jF]);
                            }
                            bFName.append("_mol" + IntToStr(i+1) + "."+ nameComps[nameComps.size()-1]);
                         
                            outTables(bFName.c_str(), allMolecules[i]);
                        }
                    }
                }
            }
        }
    }
    
    void MolGenerator::buildRefAtoms(std::vector<CrystInfo>::iterator  iCryst)
    {
        for (std::vector<AtomDict>::iterator iA=initAtoms.begin();
                            iA !=initAtoms.end(); iA++)
        {
            for (int i=-myNBDepth; i < myNBDepth+1; i++)
            {
                for (int j=-myNBDepth; j <myNBDepth+1; j++)
                {
                    for (int k=-myNBDepth; k <myNBDepth+1; k++)
                    {
                    
                        std::string label = "_"+ IntToStr(i) + IntToStr(j) + IntToStr(k);
                        AtomDict aAtom(*iA);
                        if(!(i==0 && j==0 && k==0))
                        {
                            aAtom.sId = label;
                            aAtom.symmOp = iA->symmOp;
                            aAtom.seriNum = (int)refAtoms.size();
                            aAtom.fracCoords[0] = iA->fracCoords[0] + i;
                            aAtom.fracCoords[1] = iA->fracCoords[1] + j;
                            aAtom.fracCoords[2] = iA->fracCoords[2] + k;
                        
                            FractToOrtho(aAtom.fracCoords, aAtom.coords, iCryst->itsCell->a,
                                iCryst->itsCell->b, iCryst->itsCell->c, iCryst->itsCell->alpha,
                                iCryst->itsCell->beta, iCryst->itsCell->gamma);
                            refAtoms.push_back(aAtom);
                        }
                        
                        // std::cout << "size of refAtoms " << refAtoms.size() << std::endl;
                    }
                }
            }
        }
        
        /*
        for (std::vector<AtomDict>::iterator iRA=refAtoms.begin();
                iRA != refAtoms.end(); iRA++)
        {
          
                std::cout << "For atom " << iRA->id << " and ref ID " << iRA->sId  
                          << " : " << std::endl;
                std::cout << "x = " << iRA->fracCoords[0] << std::endl
                          << "y = " << iRA->fracCoords[1] << std::endl
                          << "z = " << iRA->fracCoords[2] << std::endl;  
        }
        */
       
        
    }
    
    void MolGenerator::symmAtomGen(std::vector<CrystInfo>::iterator  tCrys,
                                   PeriodicTable  & tPTable)
    {   
        for (std::vector<AtomDict>::iterator iA=initAtoms.begin();
                    iA !=initAtoms.end(); iA++)
        {
            // int tDictMult = iA->symmMult;
            // iA->symmMult =1;
            for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator
                        iOp=tCrys->itsSpaceGroup->sgOp.begin();
                        iOp!=tCrys->itsSpaceGroup->sgOp.end(); iOp++)
            {
                getOneSymmAtom(iA, iOp, tCrys, tPTable);
            }
        }
        
    }
    
    void MolGenerator::getOneSymmAtom(std::vector<AtomDict>::iterator   tCurAtom,
                                      std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
                                      std::vector<CrystInfo>::iterator  tCryst,
                                      PeriodicTable & tPTab)
    {
        
        std::vector<REAL>  startFracCoords, endFracCoords;
        for (unsigned i=0; i <3; i++)
        {
            startFracCoords.push_back(tCurAtom->fracCoords[i]);
            endFracCoords.push_back(0.0);
        }
        startFracCoords.push_back(1.0);
        endFracCoords.push_back(0.0);
        
        matMultVec(tOp->second, startFracCoords, endFracCoords);
        //TranslateIntoUnitCell(endFracCoords, tFracCoords);
        
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
        
        bool tUnique = checkUniqueAtom(endFracCoords);
        
        if (tUnique)
        {   
            std::cout << "The symmetrically generated atom should be included in the unit cell "
                      << std::endl;
            
            AtomDict tAtom(*tCurAtom);
            tAtom.symmOp=tOp->first;
            tAtom.coords.clear();
            //std::cout << tAtom.id << std::endl;
            //std::cout << "symm atom ocp " << tAtom.ocp << std::endl;
            for (std::vector<REAL>::iterator iFX=endFracCoords.begin();
                    iFX !=endFracCoords.end(); iFX++)
            {
                tAtom.fracCoords.push_back(*iFX);
            }
            
            FractToOrtho(endFracCoords, tAtom.coords, tCryst->itsCell->a,
                     tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                     tCryst->itsCell->beta, tCryst->itsCell->gamma);
            tAtom.seriNum = (int)allAtoms.size();
            allAtoms.push_back(tAtom);
            tAtom.seriNum = (int)refAtoms.size();
            refAtoms.push_back(tAtom);
            //std::cout << "the last atom in allAtoms is now atom " << allAtoms[(int)allAtoms.size()-1].seriNum << std::endl;
            //std::cout << "the last atom in refAtoms is now atom " << refAtoms[(int)refAtoms.size()-1].seriNum << std::endl;
            /*
            std::cout << "symm-site coords " << std::endl;
            for (std::vector<REAL>::iterator iX=tAtom.coords.begin();
                    iX != tAtom.coords.end(); iX++)
            {
                std::cout << *iX << std::endl;
            }
            */
            // std::cout << "number of atoms the unit cell " << allAtoms.size() << std::endl;
            addOneSetRefAtoms(tAtom, tCryst);
            // std::cout << "number of ref atoms " << refAtoms.size() << std::endl;
            //tCurAtom->symmMult++;
            
        }
        
        /*
        // check symmMult
        if (tCurAtom->symmMult !=tDictMult)
        {
            std::cout << "Calculated atom multiply is  " << tCurAtom->symmMult 
                      << "  Value read from the cif is " << tDictMult << std::endl;
            exit(1);
        } 
        */  
    }
    
    bool MolGenerator::checkUniqueAtom(std::vector<REAL> & tFracX)
    {
        bool tU = true;
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            
            if (fabs(iA->fracCoords[0]-tFracX[0]) <0.0005 &&
                fabs(iA->fracCoords[1]-tFracX[1]) <0.0005 &&
                fabs(iA->fracCoords[2]-tFracX[2]) <0.0005)
            {
                tU = false;
                break;
            }
                    
        }
        
        return tU;
        
    }
    
    void MolGenerator::addOneSetRefAtoms(AtomDict & tCurAtom,
                                         std::vector<CrystInfo>::iterator  tCryst)
    {
        for (int i=-myNBDepth; i <myNBDepth+1; i++)
        {
            for (int j=-myNBDepth; j < myNBDepth+1; j++)
            {
                for (int k=-myNBDepth; k <myNBDepth+1; k++)
                {
                    
                    
                    if(!(i==0 && j==0 && k==0))
                    {
                        AtomDict aAtom(tCurAtom);
                        std::string label = "_"+ IntToStr(i) + IntToStr(j) + IntToStr(k);
                        aAtom.sId = label;
                        aAtom.seriNum = (int)refAtoms.size();
                        aAtom.fracCoords[0] = tCurAtom.fracCoords[0] + i;
                        aAtom.fracCoords[1] = tCurAtom.fracCoords[1] + j;
                        aAtom.fracCoords[2] = tCurAtom.fracCoords[2] + k;
                        FractToOrtho(aAtom.fracCoords, aAtom.coords, tCryst->itsCell->a,
                                     tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                                     tCryst->itsCell->beta, tCryst->itsCell->gamma);
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
    
    // get atom connections within an unit cell
    void MolGenerator::setUniqueAtomLinks(PeriodicTable & tPTab)
    {
        REAL covalent_sensitivity =0.35;
        
        for (unsigned i=0; i < allAtoms.size(); i++)
        {
            for (unsigned j=i+1; j < allAtoms.size(); j++)
            {
                REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                std::vector<REAL> linkRange;
                getBondingRangePairAtoms(allAtoms[i], allAtoms[j],
                                         covalent_sensitivity, tPTab,
                                         linkRange);
                if (linkRange[0] >0.20 && linkRange[1] >0.20)
                {
                    if (rD > linkRange[0] && rD < linkRange[1])
                    {
                        setOneUniqueBondCell(i, j, rD);
                        allAtoms[i].connAtoms.push_back(j);
                        allAtoms[j].connAtoms.push_back(i);
                    }
                }
            }
        }
        
        std::cout << "Number of atoms in cell " << allAtoms.size() << std::endl;
        std::cout << "Number of bonds (not unique bonds) : " << bonds.size() << std::endl;
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin(); 
                iAt !=allAtoms.end(); iAt++)
        {
            std::cout << "atom " << iAt->seriNum << " is connected to the following atoms :"
                      << std::endl;
            for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                    iNB !=iAt->connAtoms.end(); iNB++)
            {
                std::cout << "atom " << *iNB << std::endl;
            }
        }
       
        
    }
    
    void MolGenerator::getUniqueBonds(PeriodicTable & tPTab)
    {
        REAL covalent_sensitivity =0.35;
        NeighbListDict  tNBListOfSystem;
        
        int aDim  = 3;
        int aMode = 0;
        LIBMOL::REAL tCellL     = 3.5;
        LIBMOL::REAL tCellShell = 0.5;
        
        tNBListOfSystem.building(refAtoms, aDim, tCellL, tCellShell, aMode);
        
        // std::cout << "NB list for refAtoms set " << std::endl;   
        
        
        // std::vector<std::string>   existBondID;
        int j=0;
        for (unsigned i=0; i <refAtoms.size(); i++)
        {
            if (refAtoms[i].sId=="")
            {
                j++;
                std::cout << "Look for bonds to atom " << refAtoms[i].id 
                           << " its sID " << refAtoms[i].sId << std::endl;
                std::cout << "Its has " << refAtoms[i].neighbAtoms.size() 
                          << " neighbor atoms. " <<std::endl;
                for (std::vector<int>::iterator iNB=refAtoms[i].neighbAtoms.begin();
                        iNB !=refAtoms[i].neighbAtoms.end(); iNB++)
                {
                    REAL rD = distanceV(refAtoms[i].coords, refAtoms[(*iNB)].coords);
                    std::vector<REAL> bondRange;
                    getBondingRangePairAtoms(refAtoms[i], refAtoms[(*iNB)],
                                         covalent_sensitivity, tPTab, 
                                         bondRange);
                    
                    std::cout << "NB Atom " << refAtoms[(*iNB)].id << " its sID " 
                              << refAtoms[(*iNB)].sId << std::endl;
                    std::cout << "Distance " << rD << std::endl;
                    std::cout << "Range between " << bondRange[0]
                              << " and " << bondRange[1] << std::endl;
                              
                
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
        
        
    }
    
    void MolGenerator::getBondingRangePairAtoms(AtomDict          & tAtm1, 
                                                AtomDict          & tAtm2, 
                                                REAL              tExtraD,
                                                PeriodicTable     & tPTab,
                                                std::vector<REAL> & tRange)
    {
        tRange.clear();
        tRange.push_back(-1.0);
        tRange.push_back(-1.0);
        
        //std::cout << "atom 1 " << tAtm1.chemType << std::endl;
        //std::cout << "atom 2 " << tAtm2.chemType << std::endl;
        
        if (tPTab.elements.find(tAtm1.chemType) !=tPTab.elements.end()
            && tPTab.elements.find(tAtm2.chemType) !=tPTab.elements.end())
        {
            if (tPTab.elemProps[tAtm1.chemType].find("cova") 
                     !=tPTab.elemProps[tAtm1.chemType].end()
                && tPTab.elemProps[tAtm2.chemType].find("cova") 
                     !=tPTab.elemProps[tAtm2.chemType].end() )
            {
                if(tPTab.elemProps[tAtm1.chemType]["cova"] > 0.2
                   && tPTab.elemProps[tAtm2.chemType]["cova"] > 0.2)
                {
                    tRange[0] =  tPTab.elemProps[tAtm1.chemType]["cova"]
                                +tPTab.elemProps[tAtm2.chemType]["cova"]
                                -tExtraD;
                    tRange[1] =  tPTab.elemProps[tAtm1.chemType]["cova"]
                                +tPTab.elemProps[tAtm2.chemType]["cova"]
                                +tExtraD;
                }
                else
                {
                    std::cout << "Bug in at least one of atom covalent radius! "
                              << std::endl << "atom " << tAtm1.id << " coval: "
                              << tPTab.elemProps[tAtm1.chemType]["cova"] << std::endl
                              << "atom " << tAtm2.id  << " coval: " 
                              << tPTab.elemProps[tAtm2.chemType]["cova"] << std::endl;
                }
            }
            else 
            {
                std::cout << "Bug! At least one of covalent radius for "
                          << tAtm1.id << " or " << tAtm2.id 
                          << " is not defined the internal periodic table" 
                          << std::endl;
            }
        }
        else
        {
            std::cout << "Bug! At least one element im "
                      << tAtm1.id << " or " << tAtm2.id 
                      << " is not defined the internal periodic table " 
                      << std::endl;
        }
    }
    
    void MolGenerator::setOneUniqueBondCrys(int tIdxAtm1, int tIdxAtm2, 
                                            REAL rD)
    {
        
        ID id1=refAtoms[tIdxAtm1].id, id2 = refAtoms[tIdxAtm2].id;
        bool tFound=false;
        
        for (std::vector<BondDict>::iterator iB=bonds.begin();
                iB !=bonds.end(); iB++)
        {
            ID id3=refAtoms[iB->atomsIdx[0]].id, 
               id4=refAtoms[iB->atomsIdx[1]].id;
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
        
        if (!tFound)
        {
            BondDict aBond;
            aBond.atomsIdx.push_back(tIdxAtm1);
            aBond.atomsIdx.push_back(tIdxAtm2);
            aBond.value = rD;
            bonds.push_back(aBond);
            std::cout << "a bond between " << refAtoms[tIdxAtm1].id << " and "
                      << refAtoms[tIdxAtm2].id << " is added to the b_list " 
                      << std::endl << "Its bond length is " << aBond.value 
                      << std::endl;
        }
    }
    
    void MolGenerator::setOneUniqueBondCell(int tIdxAtm1, int tIdxAtm2, 
                                            REAL rD)
    {
        
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
        
        if (!tFound)
        {
            BondDict aBond;
            aBond.atomsIdx.push_back(tIdxAtm1);
            aBond.atomsIdx.push_back(tIdxAtm2);
            aBond.value = rD;
            bonds.push_back(aBond);
            std::cout << "a bond between " << allAtoms[tIdxAtm1].id << " and "
                      << allAtoms[tIdxAtm2].id << " is added to the bond_list_cell " 
                      << std::endl << "Its bond length is " << aBond.value 
                      << std::endl;
        }
    }
    
    void MolGenerator::getUniqAngles()
    {
        for(int i=0; i < (int)allAtoms.size(); i++)
        {   
            for (int j=0; j < (int)allAtoms[i].connAtoms.size(); j++)
            {
                for (int k=j+1; k < (int)allAtoms[i].connAtoms.size(); k++)
                {
                    int i1 = allAtoms[i].connAtoms[j];
                    int i2 = allAtoms[i].connAtoms[k];
                    //std::cout << "Angle between " << allAtoms[i].id 
                    //          << "(center) and " << allAtoms[i1].id
                    //          << " and " << allAtoms[i2].id << std::endl;
                        
                    AngleDict aAng;
                    aAng.anchorID  = allAtoms[i].id;
                    aAng.anchorPos = i;
                    aAng.atoms.push_back(i);
                       
                    if ((int) allAtoms[i1].connAtoms.size() >=
                        (int) allAtoms[i1].connAtoms.size())
                    {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                    }
                    else
                    {
                        aAng.atoms.push_back(i2);
                        aAng.atoms.push_back(i1);
                    }
                    std::vector<REAL> tV1, tV2;    
                    for (int iC=0; iC < (int)allAtoms[i].coords.size(); iC++)
                    {
                        tV1.push_back(allAtoms[i1].coords[iC]-allAtoms[i].coords[iC]);
                        tV2.push_back(allAtoms[i2].coords[iC]-allAtoms[i].coords[iC]);
                    }
                    aAng.value        = getAngle2V(tV1, tV2);
                    aAng.sigValue     = 3.0;
                    aAng.numCodValues = 0;
                    angles.push_back(aAng);
                }
            }
        }   
        
        for (std::vector<AngleDict>::iterator iA=angles.begin();
                iA !=angles.end(); iA++)
        {
            std::cout << "angle between " << allAtoms[iA->atoms[0]].id 
                      << "(center) and " <<  allAtoms[iA->atoms[1]].id
                      << " and " << allAtoms[iA->atoms[2]].id 
                      << " is " << iA->value*PID180 << std::endl;
        }
    }
    
    void MolGenerator::getMolByEqClassInCell()
    {
        std::map<int, int> classNum;
        classNum[0] = 0;
        for (unsigned i=1; i < allAtoms.size(); i++)
        {
            classNum[i] = i;
            for (unsigned j=0; j <=i-1;j++)
            {
                classNum[j]=classNum[classNum[j]];
                if (std::find(allAtoms[i].connAtoms.begin(),
                              allAtoms[i].connAtoms.end(), j)
                          !=allAtoms[i].connAtoms.end())
                {
                    classNum[classNum[classNum[j]]]=i;
                }
            }
        }
        
        // final sweeping 
        for (unsigned i=0; i < allAtoms.size(); i++ )
        {
            classNum[i]=classNum[classNum[i]];     
        }
        
        
        // get molecules in a unit cell
        for (unsigned i=0; i < classNum.size(); i++)
        {
            moleculesInCell[classNum[i]].push_back(i);
        }
        
        std::cout << "Number of molecules is " << moleculesInCell.size() << std::endl;
        
        int i_mol =0;
        for (std::map<unsigned, std::vector<int> >::iterator iMol=moleculesInCell.begin(); 
                  iMol !=moleculesInCell.end(); iMol++)
        {
            i_mol++;
            std::cout << "Molecule " << i_mol << " contains " 
                      << iMol->second.size() << " atoms " << std::endl 
                      << "The are: " << std::endl;
            for (std::vector<int>::iterator iAt=iMol->second.begin();
                    iAt !=iMol->second.end(); iAt++)
            {
                std::cout << "Atom " << allAtoms[*iAt].seriNum << " "
                          << allAtoms[*iAt].id << std::endl;
            }                   
        }
        
        getMolsInCell();
     
    }
    
    void MolGenerator::getMolsInCell()
    {
        int i=0;
        for (std::map<unsigned, std::vector<int> >::iterator iMol=moleculesInCell.begin(); 
                  iMol !=moleculesInCell.end(); iMol++)
        {
            // check if there is any atom with occupancy less than 0.95 
            // and then set individual molecules
            if (checkAtomOcp(iMol->second))
            {
                Molecule aMol;
                aMol.seriNum = i;
                i++;
            
                for (std::vector<int>::iterator iAt=iMol->second.begin(); 
                        iAt !=iMol->second.end(); iAt++)
                {
                    aMol.atoms.push_back(allAtoms[*iAt]);
                }
                
                for (std::vector<BondDict>::iterator iB=bonds.begin();
                      iB !=bonds.end(); iB++)
                {
                    
                    if (iB->atomsIdx.size()==2)
                    {
                        
                        if ((std::find(iMol->second.begin(), iMol->second.end(), iB->atomsIdx[0])
                               != iMol->second.end()) 
                             && (std::find(iMol->second.begin(), iMol->second.end(), iB->atomsIdx[1])
                               != iMol->second.end()) )
                        {
                            aMol.bonds.push_back(*iB);
                        }
                    }
                }
                
                for (std::vector<AngleDict>::iterator iAn=angles.begin();
                        iAn !=angles.end(); iAn++)
                {
                    if (iAn->atoms.size() == 3)
                    {
                        if ( (std::find(iMol->second.begin(), iMol->second.end(), iAn->atoms[0])
                                !=iMol->second.end())  &&
                             (std::find(iMol->second.begin(), iMol->second.end(), iAn->atoms[1])
                                !=iMol->second.end())  && 
                             (std::find(iMol->second.begin(), iMol->second.end(), iAn->atoms[2])
                                !=iMol->second.end()) )
                        {
                            aMol.angles.push_back(*iAn);
                        }
                    }
                }
                allMolecules.push_back(aMol);
            }
        }
       
    }
    
    bool MolGenerator::checkAtomOcp(std::vector<int>& tMol)
    {
        bool tOcpOne = true;
        for (std::vector<int>::iterator iAt=tMol.begin();
                iAt !=tMol.end(); iAt++)
        {
            //std::cout << "atom " << *iAt << " id: " << allAtoms[*iAt].id
            //          << " ocp " << allAtoms[*iAt].ocp << std::endl;
            if (allAtoms[*iAt].ocp < 0.99)
            {
                tOcpOne = false;
                break;
            }
        }
        
        return tOcpOne;
        
    }
    
    void MolGenerator::getMolByEqClassInCrys()
    {
        
        std::map<int, int> classNum;
        classNum[0] = 0;
        for (unsigned i=1; i < refAtoms.size(); i++)
        {
            classNum[i] = i;
            for (unsigned j=0; j <=i-1;j++)
            {
                classNum[j]=classNum[classNum[j]];
                if (std::find(refAtoms[i].connAtoms.begin(),
                              refAtoms[i].connAtoms.end(), j)
                          !=refAtoms[i].connAtoms.end())
                {
                    classNum[classNum[classNum[j]]]=i;
                }
            }
        }
        
        // final sweeping 
        for (unsigned i=0; i < refAtoms.size(); i++ )
        {
            classNum[i]=classNum[classNum[i]];     
        }
        
        
        // get molecules in a finite crystal 
        for (unsigned i=0; i < classNum.size(); i++)
        {
            moleculesInCryst[classNum[i]].push_back(i);
        }
    }
    
    void MolGenerator::outTables(FileName tOutName,
                                 Molecule & tMol)
    {
        std::ofstream outTableF(tOutName);
        
        if (outTableF.is_open())
        {
            if (tMol.atoms.size() >0)
            {
                outTableF << "loop_" << std::endl
                          << "_chem_comp_atom.serial_number" << std::endl
                          << "_chem_comp_atom.atom_id " << std::endl
                          << "_chem_comp_atom.element_symbol" << std::endl
                          << "_chem_comp_atom.cod_type" << std::endl;
                for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin(); 
                        iAt !=tMol.atoms.end(); iAt++)
                {
                    outTableF << std::setw(6) << iAt->seriNum+1   
                              << std::setw(6) << iAt->id    
                              << std::setw(5) << iAt->chemType << "    "
                              << iAt->codClass << std::endl;
                }
                
                outTableF << std::endl;
            }
            else
            {
                std::cout << "There is no atoms in the molecule" << std::endl;
            }
            
            if (tMol.bonds.size() >0)
            {
                std::cout << std::endl;
                outTableF << "loop_" << std::endl  
                      << "_chem_comp_bond.bond_serial_number"  << std::endl
                      << "_chem_comp_bond.atom1_serial_number" << std::endl
                      << "_chem_comp_bond.atom2_serial_number" << std::endl
                      << "_chem_comp_bond.atom1_element_symbol" << std::endl
                      << "_chem_comp_bond.atom2_element_symbol" << std::endl
                      << "_chem_comp_bond.value_dist"<< std::endl;
                int nBo = 1;
                for (std::vector<BondDict>::iterator iBo=tMol.bonds.begin();
                        iBo !=tMol.bonds.end(); iBo++)
                {
                    outTableF << std::setw(6) << nBo
                              << std::setw(6) << iBo->atomsIdx[0] + 1 
                              << std::setw(6) << iBo->atomsIdx[1] + 1
                              << std::setw(4) << allAtoms[iBo->atomsIdx[0]].chemType
                              << std::setw(4) << allAtoms[iBo->atomsIdx[1]].chemType
                              << std::setw(10)<< iBo->value << std::endl;
                    nBo++;
                }
            
                std::cout << std::endl;
            }
            else
            {
                std::cout << "There is no bonds in the molecule" << std::endl;
            }
            
            if (tMol.angles.size() > 0)
            {
                outTableF << "loop_" << std::endl
                          << "_chem_comp_angle.angle_serial_number" << std::endl
                          << "_chem_comp_angle.atom1_serial_number" << std::endl
                          << "_chem_comp_angle.atom2_serial_number" << std::endl
                          << "_chem_comp_angle.atom3_serial_number" << std::endl
                          << "_chem_comp_angle.atom1_element_symbol" << std::endl
                          << "_chem_comp_angle.atom2_element_symbol" << std::endl
                          << "_chem_comp_angle.atom3_element_symbol" << std::endl
                          << "_chem_comp_angle.value_angle"          << std::endl;
                int nAn = 1;
                for (std::vector<AngleDict>::iterator iAn=tMol.angles.begin();
                        iAn !=tMol.angles.end(); iAn++)
                {
                    outTableF << std::setw(6)  << nAn
                              << std::setw(6)  << iAn->atoms[0] +1 
                              << std::setw(6)  << iAn->atoms[1] +1
                              << std::setw(6)  << iAn->atoms[2] +1
                              << std::setw(4)  << allAtoms[iAn->atoms[0]].chemType 
                              << std::setw(4)  << allAtoms[iAn->atoms[1]].chemType 
                              << std::setw(4)  << allAtoms[iAn->atoms[2]].chemType 
                              << std::setw(10) << iAn->value*PID180 << std::endl;
                    nAn++;
                }
                
                std::cout << std::endl;
            }
            
            outTableF.close();
            
        }
        else
        {
            std::cout << tOutName << " can not been opened for writing ! "
                      << std::endl;
        }
    }
}