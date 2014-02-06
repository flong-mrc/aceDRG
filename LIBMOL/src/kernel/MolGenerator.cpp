/* 
 * File:   MolGenerator.cpp
 * Author: flong
 *
 * Generator a molecule from input sets of atoms and symmetry operators 
 * 
 * Created on August 13, 2013, 5:56 PM
 */

#include "MolGenerator.h"

namespace LIBMOL
{
    MolGenerator::MolGenerator():myNBDepth(1)
    {
    }
    
    
    MolGenerator::MolGenerator(const GenCifFile& tCifObj, int tNBDepth)
    {
        
        myNBDepth = tNBDepth;
        
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
        std::cout << "number of unique atoms " << initAtoms.size() << std::endl;
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
    
    
    void MolGenerator::execute()
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
                }
                
                buildRefAtoms(iCryst);
                 
                symmAtomGen(iCryst, aPTable);
                
                getUniqueBonds(aPTable);
                
                getUniqAngles();
                
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
                            aAtom.id += label;
                        
                            aAtom.fracCoords[0] = iA->fracCoords[0] + i;
                            aAtom.fracCoords[1] = iA->fracCoords[1] + j;
                            aAtom.fracCoords[2] = iA->fracCoords[2] + k;
                        
                            FractToOrtho(aAtom.fracCoords, aAtom.coords, iCryst->itsCell->a,
                                iCryst->itsCell->b, iCryst->itsCell->c, iCryst->itsCell->alpha,
                                iCryst->itsCell->beta, iCryst->itsCell->gamma);
                        }
                        refAtoms.push_back(aAtom);
                    }
                }
            }
        }
        
        /*
        for (std::vector<AtomDict>::iterator iRA=refAtoms.begin();
                iRA != refAtoms.end(); iRA++)
        {
            std::cout << "For atom " << iRA->id << " : " << std::endl;
            std::cout << "x = " << iRA->coords[0] << std::endl
                      << "y = " << iRA->coords[1] << std::endl
                      << "z = " << iRA->coords[2] << std::endl;
        }
        */
        
        
    }
    
    void MolGenerator::symmAtomGen(std::vector<CrystInfo>::iterator  tCrys,
                                   PeriodicTable  & tPTable)
    {   
        for (std::vector<AtomDict>::iterator iA=initAtoms.begin();
                    iA !=initAtoms.end(); iA++)
        {
            int tDictMult = iA->symmMult;
            iA->symmMult =1;
            for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator
                        iOp=tCrys->itsSpaceGroup->sgOp.begin();
                        iOp!=tCrys->itsSpaceGroup->sgOp.end(); iOp++)
            {
                getOneSymmAtom(iA, iOp, tCrys, tPTable);
            }
            std::cout << "After symm Opts: " << std::endl
                      << "atom site-symm-mult in the input file " 
                      << tDictMult << std::endl
                      << "atom site-symm-mult from symm opts " 
                      << iA->symmMult << std::endl;
        
        }
        
    }
    
    void MolGenerator::getOneSymmAtom(std::vector<AtomDict>::iterator   tCurAtom,
                                      std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
                                      std::vector<CrystInfo>::iterator  tCryst,
                                      PeriodicTable & tPTab)
    {
        
        std::vector<REAL>   tFracCoords, tCoords, startFracCoords, endFracCoords;
        for (unsigned i=0; i <3; i++)
        {
            startFracCoords.push_back(tCurAtom->fracCoords[i]);
            endFracCoords.push_back(0.0);
        }
        startFracCoords.push_back(1.0);
        endFracCoords.push_back(0.0);
        
        matMultVec(tOp->second, startFracCoords, endFracCoords);
        
        TranslateIntoUnitCell(endFracCoords, tFracCoords);
        
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
                
        std::cout << "symm-site fract coords " << std::endl;
        for (std::vector<REAL>::iterator iX=tFracCoords.begin();
                iX != tFracCoords.end(); iX++)
        {
            std::cout << *iX << std::endl;
        }
       
        bool tUnique = checkUniqueAtom(tFracCoords);
        
        
        std::cout << "initial  coords " << std::endl;
        for (std::vector<REAL>::iterator iX=tCurAtom->coords.begin();
                iX != tCurAtom->coords.end(); iX++)
        {
            std::cout << *iX << std::endl;
        }
        
        std::cout << "symm-site coords " << std::endl;
        for (std::vector<REAL>::iterator iX=tCoords.begin();
                iX != tCoords.end(); iX++)
        {
            std::cout << *iX << std::endl;
        }
        
        if (tUnique)
        {   
            std::cout << "The symmetrically generated atom should be included in the unit cell "
                      << std::endl;
            
            AtomDict tAtom(*tCurAtom);
            tAtom.sId = tOp->first;
            tAtom.coords.clear();
            
            for (std::vector<REAL>::iterator iFX=tFracCoords.begin();
                    iFX !=tFracCoords.end(); iFX++)
            {
                tAtom.fracCoords.push_back(*iFX);
            }
            
            FractToOrtho(tFracCoords, tAtom.coords, tCryst->itsCell->a,
                     tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                     tCryst->itsCell->beta, tCryst->itsCell->gamma);
            
            allAtoms.push_back(tAtom);
       
            std::cout << "number of atoms the unit cell " << allAtoms.size() << std::endl;
            addOneSetRefAtoms(tAtom, tCryst);
            std::cout << "number of ref atoms " << refAtoms.size() << std::endl;
            tCurAtom->symmMult++;
            
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
            
            if (fabs(iA->fracCoords[0]-tFracX[0]) <0.0001 &&
                fabs(iA->fracCoords[1]-tFracX[1]) <0.0001 &&
                fabs(iA->fracCoords[2]-tFracX[2]) <0.0001)
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
                    std::string label = "_"+ IntToStr(i) + IntToStr(j) + IntToStr(k);
                    AtomDict aAtom(tCurAtom);
                    if(!(i==0 && j==0 && k==0))
                    {
                        aAtom.sId = label;
                    
                        aAtom.fracCoords[0] = tCurAtom.fracCoords[0] + i;
                        aAtom.fracCoords[1] = tCurAtom.fracCoords[1] + j;
                        aAtom.fracCoords[2] = tCurAtom.fracCoords[2] + k;
                        
                        FractToOrtho(aAtom.fracCoords, aAtom.coords, tCryst->itsCell->a,
                                     tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                                     tCryst->itsCell->beta, tCryst->itsCell->gamma);
                        refAtoms.push_back(aAtom);
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
                        allAtoms[i].connAtoms.push_back(j);
                        allAtoms[j].connAtoms.push_back(i);
                    }
                }
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
        
        std::vector<std::string>   existBondID;
        
        for (unsigned i=0; i <refAtoms.size(); i++)
        {
            if (refAtoms[i].fracCoords[0] >=0.0 && refAtoms[i].fracCoords[0] <1.0
                 && refAtoms[i].fracCoords[1] >=0.0 && refAtoms[i].fracCoords[1] <1.0
                 && refAtoms[i].fracCoords[2] >=0.0 && refAtoms[i].fracCoords[2] <1.0 )
            {
                for (std::vector<int>::iterator iNB=refAtoms[i].neighbAtoms.begin();
                        iNB !=refAtoms[i].neighbAtoms.end(); iNB++)
                {
                    REAL rD = distanceV(refAtoms[i].coords, refAtoms[(*iNB)].coords);
                    std::vector<REAL> bondRange;
                    getBondingRangePairAtoms(refAtoms[i], refAtoms[(*iNB)],
                                         covalent_sensitivity, tPTab, 
                                         bondRange);
                
                    if (bondRange[0] >0.20 && bondRange[1] >0.20)
                    {
                        if (rD > bondRange[0] && rD < bondRange[1])
                        {
                            setOneUniqueBond(i, *iNB, rD);
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
        
        
        if (tPTab.elements.find(tAtm1.chemType) !=tPTab.elements.end()
            && tPTab.elements.find(tAtm2.chemType) !=tPTab.elements.end())
        {
            if (tPTab.elements[tAtm1.chemType].find("cova") 
                     !=tPTab.elements[tAtm1.chemType].end()
                && tPTab.elements[tAtm2.chemType].find("cova") 
                     !=tPTab.elements[tAtm2.chemType].end() )
            {
                if(tPTab.elements[tAtm1.chemType]["cova"] > 0.2
                   && tPTab.elements[tAtm2.chemType]["cova"] > 0.2)
                {
                    tRange[0] =  tPTab.elements[tAtm1.chemType]["cova"]
                                +tPTab.elements[tAtm2.chemType]["cova"]
                                -tExtraD;
                    tRange[1] =  tPTab.elements[tAtm1.chemType]["cova"]
                                +tPTab.elements[tAtm2.chemType]["cova"]
                                +tExtraD;
                }
                else
                {
                    std::cout << "Bug! At least one of covalent radius for "
                              << tAtm1.id << " or " << tAtm2.id 
                              << " is wrong " << std::endl;
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
    
    void MolGenerator::setOneUniqueBond(int tIdxAtm1, int tIdxAtm2, 
                                        REAL rD)
    {
        BondDict aBond;
        aBond.atomsIdx.push_back(tIdxAtm1);
        aBond.atomsIdx.push_back(tIdxAtm2);
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
            aBond.value = rD;
            bonds.push_back(aBond);
        }
    }
    
    void MolGenerator::getUniqAngles()
    {
        
    }
    
    void MolGenerator::getMolsInCell()
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
}