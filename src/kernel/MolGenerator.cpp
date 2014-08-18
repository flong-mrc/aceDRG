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
        // std::cout << "number of initial atoms " << initAtoms.size() << std::endl;
        
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
                    //FractToOrtho(iA->fracCoords, iA->coords, iCryst->itsCell->a,
                    //            iCryst->itsCell->b, iCryst->itsCell->c, iCryst->itsCell->alpha,
                    //            iCryst->itsCell->beta, iCryst->itsCell->gamma);
                    if (iA->ocp < 1.0000001)
                    {
                        packAtomIntoCell((*iA));
                        iA->sId ="555";
                        allAtoms.push_back(*iA);
                        refAtoms.push_back(*iA);
                    }
                }
                
                //outPDB("initAtoms.pdb", "UNL", initAtoms);
                
                std::cout << "Number of atoms read from the input file "
                          << initAtoms.size() << std::endl; 
                
                symmAtomGen(iCryst, aPTable);
                
                std::cout << "number of atoms in a center unit cell " 
                          << allAtoms.size() << std::endl;
                
                buildRefAtoms(iCryst);     
                
                
                // getUniqueBonds(aPTable);
                getUniqueAtomLinks(aPTable, iCryst); 
                
                getMolByEqClassInCell();
                
                buildAndValidMols(aPTable, iCryst);
                
                // getAtomTypeMols();
                
                getOverallBondAndAngles();
                
                outTables(tOutName);
                
                /*
                if (allMolecules.size() > 0)
                {
               
                    for (unsigned  i=0; i < allMolecules.size(); i++)
                    {
                        if (allMolecules[i].atoms.size() >0 &&
                            allMolecules[i].bonds.size() >0 )
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
                 */
            }
        }
    }
    

    void MolGenerator::symmAtomGen(std::vector<CrystInfo>::iterator  tCrys,
                                   PeriodicTable  & tPTable)
    {   
        for (std::vector<AtomDict>::iterator iA=initAtoms.begin();
                    iA !=initAtoms.end(); iA++)
        {
            // int tDictMult = iA->symmMult;
            // iA->symmMult =1;
            if (iA->ocp < 1.000001)
            {
                for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator
                        iOp=tCrys->itsSpaceGroup->sgOp.begin();
                        iOp!=tCrys->itsSpaceGroup->sgOp.end(); iOp++)
                {
                    getOneSymmAtom(iA, iOp, tCrys, tPTable);
                }
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
        
        if (!colidAtom(endFracCoords,allAtoms))
        {   
            //std::cout << "The symmetrically generated atom should be included in the unit cell "
            //          << std::endl;
            
            AtomDict tAtom(*tCurAtom);
            tAtom.isInPreCell = false;
            tAtom.symmOp=tOp->first;
            tAtom.coords.clear();
            tAtom.fracCoords.clear();
            //std::cout << tAtom.id << std::endl;
            //std::cout << "symm atom ocp " << tAtom.ocp << std::endl;
            
            for (unsigned ix=0; ix < endFracCoords.size()-1; ix++)
            {
                tAtom.fracCoords.push_back(endFracCoords[ix]);
            }
            
            //FractToOrtho(endFracCoords, tAtom.coords, tCryst->itsCell->a,
            //         tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
            //         tCryst->itsCell->beta, tCryst->itsCell->gamma);
            
            
            tAtom.seriNum = (int)allAtoms.size();
            packAtomIntoCell(tAtom);
            if (!colidAtom(tAtom, allAtoms))
            {
                tAtom.sId = "555";
                allAtoms.push_back(tAtom);
                refAtoms.push_back(tAtom);
            }
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
    
    void MolGenerator::packAtomIntoCell(AtomDict & tAtm)
    {   
        
        
        
        for (std::vector<REAL>::iterator iX=tAtm.fracCoords.begin();
                iX != tAtm.fracCoords.end(); iX++)
        {
            if((*iX) >1.00000001)
            {
                *iX-=1.0;
            }
            else if ((*iX) < -0.00000001)
            {
                *iX+=1;
            }
        }
    }
    
    void MolGenerator::buildRefAtoms(std::vector<CrystInfo>::iterator  iCryst)
    {
        REAL fraCell=1.0, tCM=1000000.0;
        REAL aShell=7.0;
        if(iCryst->itsCell->a < tCM)
        {
            tCM = iCryst->itsCell->a;
        }
        
        if (iCryst->itsCell->b < tCM) 
        {
            tCM = iCryst->itsCell->b;
        }
        
        if (iCryst->itsCell->c < tCM) 
        {
            tCM = iCryst->itsCell->c;
        }
        
        if (fabs(tCM-1000000.0) < 0.1)
        {
            std::cout << "No cell lengths or huge values of cell lengths. Program stops " 
                      << std::endl;
            exit(1);
        }
        
        getFracReal(tCM, fraCell, aShell);
        
        std::cout << "fraCell " << fraCell << std::endl;
        
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                            iA !=allAtoms.end(); iA++)
        {
            for (int i=-myNBDepth; i < myNBDepth+1; i++)
            {
                for (int j=-myNBDepth; j <myNBDepth+1; j++)
                {
                    for (int k=-myNBDepth; k <myNBDepth+1; k++)
                    {   
                        if(!(i==0 && j==0 && k==0))
                        {
                            
                            std::vector<REAL>  tFx;
                            tFx.push_back(iA->fracCoords[0] + i);
                            tFx.push_back(iA->fracCoords[1] + j);
                            tFx.push_back(iA->fracCoords[2] + k);
                           
                            // Take a zone surround the center cell 
                            REAL tRange1 = -fraCell, tRange2 = 1+fraCell;
                            //REAL tRange1 = -0.8, tRange2 = 1.8;
                            if ((tFx[0] > tRange1 && tFx[0] <=tRange2)
                                && (tFx[1] > tRange1 && tFx[1] <=tRange2)
                                && (tFx[2] > tRange1 && tFx[2] <=tRange2)
                                && !colidAtom(tFx, refAtoms))
                            {
                                AtomDict aAtom(*iA);
                                aAtom.isInPreCell = false;
                                
                                aAtom.sId = IntToStr(5+i) + IntToStr(5+j) + IntToStr(5+k);;
                                aAtom.symmOp = iA->symmOp;
                                aAtom.seriNum = (int)refAtoms.size();
                        
                                //FractToOrtho(aAtom.fracCoords, aAtom.coords, 
                                //             iCryst->itsCell->a, iCryst->itsCell->b, 
                                //             iCryst->itsCell->c, iCryst->itsCell->alpha,
                                //             iCryst->itsCell->beta, iCryst->itsCell->gamma);
                                aAtom.fracCoords[0] = tFx[0];
                                aAtom.fracCoords[1] = tFx[1];
                                aAtom.fracCoords[2] = tFx[2];
                                refAtoms.push_back(aAtom);
                            }
                        }   
                        // std::cout << "size of refAtoms " << refAtoms.size() << std::endl;
                    }
                }
            }
        }
        
        swithAtoms(iCryst);
        
        /*
        for (std::vector<AtomDict>::iterator iRA=allAtoms.begin();
                iRA != allAtoms.end(); iRA++)
        {
            if (iRA->sId=="555" )
            {
                std::cout << "For atom " << iRA->id << " and ref ID " << iRA->sId  
                          << " : " << std::endl;
                std::cout << "x = " << iRA->coords[0] << std::endl
                          << "y = " << iRA->coords[1] << std::endl
                          << "z = " << iRA->coords[2] << std::endl;
            }
        }
        */
       
        
        
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
                        //std::string label = "_"+ IntToStr(5-i) + IntToStr(5-j) + IntToStr(5-k);
                        //aAtom.sId = label;
                        aAtom.seriNum = (int)refAtoms.size();
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
    
    void MolGenerator::swithAtoms(std::vector<CrystInfo>::iterator tCryst)
    {
        allAtoms.clear();
        for (std::vector<AtomDict>::iterator iAt = refAtoms.begin();
                  iAt !=refAtoms.end(); iAt++)
        {
            allAtoms.push_back(*iAt);
        }
        refAtoms.clear();
                
        
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            iAt->coords.clear();
            FractToOrtho(iAt->fracCoords, iAt->coords, tCryst->itsCell->a,
                         tCryst->itsCell->b, tCryst->itsCell->c, tCryst->itsCell->alpha,
                         tCryst->itsCell->beta, tCryst->itsCell->gamma);
        }
        
        std::cout << "Number of atoms in allAtoms " 
                  << allAtoms.size() << std::endl;
        
        
        int iC=0;
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->sId=="555")
            {
                iC++;
            }
        }
        std::cout << "number of atoms in the center unit cell is " 
                  << iC << std::endl;
      
    }
    
    void MolGenerator::cleanUnconnAtoms()
    {
        std::vector<AtomDict> tAtoms;
        std::map<int, int> tAtMap;
        int j=0;
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            if (allAtoms[i].connAtoms.size() !=0)
            {
                tAtoms.push_back(allAtoms[i]);
                tAtMap[i]=j;
                j++;
            }
        }
        
        
        
        //std::cout << "Before clean, number of atoms is " << allAtoms.size() << std::endl;
        //std::cout << "After clean, number of atoms is " << tAtoms.size() << std::endl;
        
        allAtoms.clear();
        
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                    iCo != iAt->connAtoms.end(); iCo++)
            {
                (*iCo)=tAtMap[(*iCo)];
            }
            iAt->seriNum = (int)allAtoms.size();
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
    void MolGenerator::setUniqueAtomLinks(PeriodicTable & tPTab)
    {
        REAL covalent_sensitivity =0.20;
        
        for (unsigned i=0; i < allAtoms.size(); i++)
        {
            for (unsigned j=i+1; j < allAtoms.size(); j++)
            {
               // REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);
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
        
        //std::cout << "Number of atoms in cell " << allAtoms.size() << std::endl;
        //std::cout << "Number of bonds (not unique bonds) : " << bonds.size() << std::endl;
        /*
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
        */
        
    }
    
    void MolGenerator::setUniqueAtomLinks(PeriodicTable & tPTab,
                                          std::vector<CrystInfo>::iterator tCryst)
    {
        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 =0.25;
        REAL covalent_sensitivity2 =0.50;
        for (unsigned i=0; i < allAtoms.size(); i++)
        {
            std::cout << std::endl << "================================================================" << std::endl;
            std::cout << "Look for bonds to atom " << allAtoms[i].id 
                      // << " its sID " << allAtoms[i].sId 
                      << " its serial number " << allAtoms[i].seriNum << std::endl;
             
            for (unsigned j=i+1; j < allAtoms.size(); j++)
            {
                // REAL rD = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                
                REAL rD = getBondLenFromFracCoords(allAtoms[i].fracCoords, allAtoms[j].fracCoords,
                                                  tCryst->itsCell->a, tCryst->itsCell->b, 
                                                  tCryst->itsCell->c, tCryst->itsCell->alpha,
                                                  tCryst->itsCell->beta, tCryst->itsCell->gamma);
                std::vector<REAL> linkRange;
                if ((!allAtoms[i].isMetal) && (!allAtoms[j].isMetal))
                {
                    covalent_sensitivity = covalent_sensitivity1;
                }
                else
                {
                    covalent_sensitivity = covalent_sensitivity2;
                }
                getBondingRangePairAtoms2(allAtoms[i], allAtoms[j],
                                         covalent_sensitivity, tPTab,
                                         linkRange);
                
                
                if (linkRange[0] >0.20 && linkRange[1] >0.20)
                {
                    //std::cout << "NB Atom " << allAtoms[j].id << " its sID " 
                    //      << allAtoms[j].sId << " serial number " 
                    //      << allAtoms[j].seriNum << std::endl;
                    //std::cout << "Distance " << rD << std::endl;
                    //std::cout << "Range between " << linkRange[0]
                    //      << " and " << linkRange[1] << std::endl;
                    // if (rD > linkRange[0] && rD < linkRange[1])
                    if (rD < linkRange[1])
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
        /*
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin(); 
                iAt !=allAtoms.end(); iAt++)
        {
            //std::cout << "atom " << iAt->seriNum << " is connected to the following atoms :"
            //          << std::endl;
            for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                    iNB !=iAt->connAtoms.end(); iNB++)
            {
                std::cout << "atom " << *iNB << std::endl;
            }
        }
         */
       
        
    }
    
    bool MolGenerator::inBonds(int tIdx1, int tIdx2, 
                               std::vector<BondDict>& tBonds)
    {
        //std::cout << "input atom " << tIdx1 << " and " << tIdx2 << std::endl;
        
        for (std::vector<BondDict>::iterator iBo=tBonds.begin();
                iBo!=tBonds.end(); iBo++)
        {
            //std::cout << "bond atom index " << iBo->atomsIdx[0]
            //          << " and " << iBo->atomsIdx[1] << std::endl;
            
            if ((iBo->atomsIdx[0]==tIdx1 && iBo->atomsIdx[1]==tIdx2)
                 ||
                 (iBo->atomsIdx[0]==tIdx2 && iBo->atomsIdx[1]==tIdx1))
            {
                return true;
            }
        }
        return false;
    }
    
    void MolGenerator::getUniqueBonds(PeriodicTable & tPTab)
    {
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
    
    
    void MolGenerator::getUniqueAtomLinks(PeriodicTable                 & tPTab, 
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
        
        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 =0.15;
        REAL covalent_sensitivity2 =0.22;
        REAL covalent_sensitivity3 =0.20;
        
        NeighbListDict  tNBListOfSystem;
        
        int aDim  = 3;
        int aMode = 0;
        LIBMOL::REAL tCellL     = 3.5;
        LIBMOL::REAL tCellShell = 0.5;
        
        tNBListOfSystem.building(allAtoms, aDim, tCellL, tCellShell, aMode);
        
       
        // std::cout << "NB list for allAtoms set " << std::endl;   
       
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                //std::cout << "atom "<< iAt->id << "(serial number " 
                //          << iAt->seriNum  <<") has " << iAt->neighbAtoms.size() 
                //          << " NB atoms " << std::endl;
                /*
                if (iAt->id =="C22")
                {
                    for (std::vector<int>::iterator iNB=iAt->neighbAtoms.begin();
                            iNB != iAt->neighbAtoms.end(); iNB++)
                    {
                        std::cout << "NB atom " << allAtoms[*iNB].id << std::endl;
                    }
                }
                */
            }
        }
       
       
        
        // std::vector<std::string>   existBondID;
        // int j=0;
        for (unsigned i=0; i <allAtoms.size(); i++)
        {
            //if (allAtoms[i].sId=="555")
            //{
                //j++;
                std::cout << "Look for bonds to atom " << allAtoms[i].id 
                          << "(serial number  " << allAtoms[i].seriNum
                          << ") " << std::endl;
                std::cout << "Its has " << allAtoms[i].neighbAtoms.size() 
                          << " neighbor atoms. " <<std::endl;
                for (std::vector<int>::iterator iNB=allAtoms[i].neighbAtoms.begin();
                        iNB !=allAtoms[i].neighbAtoms.end(); iNB++)
                {
                    // REAL rD = distanceV(refAtoms[i].coords, refAtoms[(*iNB)].coords);
                    std::cout << "NB Atom " << refAtoms[(*iNB)].id << " its sID " 
                              << refAtoms[(*iNB)].sId << std::endl;
                    REAL rD = getBondLenFromFracCoords(allAtoms[i].fracCoords, allAtoms[(*iNB)].fracCoords,
                                                       tCryst->itsCell->a, tCryst->itsCell->b, 
                                                       tCryst->itsCell->c, tCryst->itsCell->alpha,
                                                       tCryst->itsCell->beta, tCryst->itsCell->gamma);
                    
                    std::vector<REAL> bondRange;
                    if (allAtoms[i].chemType=="H" || allAtoms[*iNB].chemType=="H")
                    {
                        covalent_sensitivity = covalent_sensitivity3;
                    }
                    else if ((!allAtoms[i].isMetal) && (!allAtoms[*iNB].isMetal))
                    {
                        covalent_sensitivity = covalent_sensitivity1;
                    }
                    
                    else
                    {
                        covalent_sensitivity = covalent_sensitivity2;
                    }
                    
                    std::cout << "covalent_sensitivity=" << covalent_sensitivity << std::endl; 

                    getBondingRangePairAtoms2(allAtoms[i], allAtoms[(*iNB)],
                                         covalent_sensitivity, tPTab, 
                                         bondRange);
                    
                    
                    
                        std::cout << "Range between " << bondRange[0]
                                      << " and " << bondRange[1] << std::endl;
                        std::cout << "Distance between " << allAtoms[i].id << " and "
                                  << allAtoms[*iNB].id << " is " << rD << std::endl;
                   
                   
                    
                    
                    if (bondRange[0] >0.20 && bondRange[1] >0.20)
                    {
                        if (rD < bondRange[1])
                        {
                            
                              
                            // setOneUniqueBondCell(i, *iNB, rD);
                            if (std::find(allAtoms[i].connAtoms.begin(), allAtoms[i].connAtoms.end(), *iNB)
                                    ==allAtoms[i].connAtoms.end())
                            {
                                allAtoms[i].connAtoms.push_back(*iNB);
                            }
                            if (std::find(allAtoms[*iNB].connAtoms.begin(), allAtoms[*iNB].connAtoms.end(), i)
                                    ==allAtoms[*iNB].connAtoms.end())
                            {
                                allAtoms[*iNB].connAtoms.push_back(i);
                            }
                            
                            // if (allAtoms[i].id=="C22" && allAtoms[i].isInPreCell)
                            
                            std::cout << "a bond between " << allAtoms[i].id << " and "
                             << allAtoms[*iNB].id << " is added to the bond_list_cell " 
                             << std::endl << "Its bond length is " << rD 
                             << std::endl;
                            
                        }
                    }
                }
            //}
        }
        
        // std::cout << "Number of atoms in the unit cell considered " << j << std::endl;
        
        cleanUnconnAtoms();
     
        
        
        //int iUA=0;
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                std::cout << "Atom " << iAt->id << " (" << iAt->seriNum 
                          << ") bonds to " << iAt->connAtoms.size() << " atoms"
                          << std::endl;
                //if (iAt->isInPreCell)
                //{
                //    iUA++;
                //}
            }
        }
        
      
        //std::cout << "there are " << iUA << " atoms in ASU " << std::endl;
                
        // setUniqueAtomLinks(tPTab, tCryst);
        
        
    }
    
    void MolGenerator::getUniqueBondsMols(Molecule& tMol, 
                                          std::vector<CrystInfo>::iterator tCryst)
    {   
        std::map<std::string, int> aMA;
        std::vector<REAL> aV;
        
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                      iCo !=iAt->connAtoms.end(); iCo++)
                {
                    std::string str1=iAt->id + "_" + tMol.atoms[(*iCo)].id;
                    std::string str2=tMol.atoms[(*iCo)].id + "_" + iAt->id;
                    if ((*iCo) > iAt->seriNum &&
                         (aMA.find(str1)==aMA.end() && aMA.find(str2)==aMA.end()))
                    {
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
                        std::cout << "a bond between " << iAt->id << " and "
                                  << tMol.atoms[(*iCo)].id << " is added to the bond_list_cell " 
                                  << std::endl << "Its bond length is " << aBond.value 
                                  << std::endl;
                    }
                    else
                    {
                       REAL rD = getBondLenFromFracCoords(iAt->fracCoords, tMol.atoms[(*iCo)].fracCoords,
                                                       tCryst->itsCell->a, tCryst->itsCell->b, 
                                                       tCryst->itsCell->c, tCryst->itsCell->alpha,
                                                       tCryst->itsCell->beta, tCryst->itsCell->gamma); 
                       if (!inVectABS(aV, rD, 0.00001))
                       {
                           BondDict aBond;
                           aBond.atomsIdx.push_back(iAt->seriNum);
                           aBond.atomsIdx.push_back(*iCo);
                           aBond.value = rD;
                           tMol.bonds.push_back(aBond);
                           aV.push_back(aBond.value);            
                           std::cout << "a bond between " << iAt->id << " and "
                                     << tMol.atoms[(*iCo)].id << " is added to the bond_list_cell " 
                                     << std::endl << "Its bond length is(unique) " << aBond.value 
                                     << std::endl;
                       }
                       
                    }
                }
            }
        }
    }
    
    void MolGenerator::getUniqueBondsMols2(Molecule& tMol)
    {
        
        std::vector<REAL> aV;
        
        for (std::vector<BondDict>::iterator iB=tMol.allBonds.begin();
                iB !=tMol.allBonds.end(); iB++)
        {
            // Consider bonds involve atoms in an assym-cell only 
            if (tMol.atoms[iB->atomsIdx[0]].isInPreCell 
                || tMol.atoms[iB->atomsIdx[1]].isInPreCell)
            {
                if (!inVectABS(aV, iB->value, 0.00001))
                {
                    tMol.bonds.push_back(*iB);
                }
                
                aV.push_back(iB->value);
            }
        }
    }
    
    void MolGenerator::getAllBondsInOneMol(Molecule & tMol, 
                                           std::vector<CrystInfo>::iterator tCryst)
    {
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                iAt != tMol.atoms.end(); iAt++)
        {
            for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                      iCo !=iAt->connAtoms.end(); iCo++)
            {
                if ((*iCo) > iAt->seriNum)         
                {
                    REAL rD = getBondLenFromFracCoords(iAt->fracCoords, tMol.atoms[(*iCo)].fracCoords,
                                                       tCryst->itsCell->a, tCryst->itsCell->b, 
                                                       tCryst->itsCell->c, tCryst->itsCell->alpha,
                                                       tCryst->itsCell->beta, tCryst->itsCell->gamma);
                    BondDict aBond;
                    aBond.atomsIdx.push_back(iAt->seriNum);
                    aBond.atomsIdx.push_back(*iCo);
                    aBond.value = rD;
                    aBond.isInSameRing = aBond.checkIfInSameRing(tMol.atoms, iAt->seriNum, *iCo);
                    tMol.allBonds.push_back(aBond);            
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
                }
            }
        }
        
        getUniqueBondsMols2(tMol);
   
    }
    
    
    void MolGenerator::getBondingRangePairAtoms(AtomDict          & tAtm1, 
                                                AtomDict          & tAtm2, 
                                                REAL              tSens,
                                                PeriodicTable     & tPTab,
                                                std::vector<REAL> & tRange)
    {
        tRange.clear();
        tRange.push_back(-1.0);
        tRange.push_back(-1.0);
       
        
        if (tPTab.elemProps.find(tAtm1.chemType) !=tPTab.elemProps.end()
            && tPTab.elemProps.find(tAtm2.chemType) !=tPTab.elemProps.end())
        {
            if (tPTab.elemProps[tAtm1.chemType].find("cova") 
                     !=tPTab.elemProps[tAtm1.chemType].end()
                && tPTab.elemProps[tAtm2.chemType].find("cova") 
                     !=tPTab.elemProps[tAtm2.chemType].end() )
            {
                if(tPTab.elemProps[tAtm1.chemType]["cova"] > 0.2
                   && tPTab.elemProps[tAtm2.chemType]["cova"] > 0.2)
                {
                    REAL tRad, tExtraD;
                    tRad  =  tPTab.elemProps[tAtm1.chemType]["cova"]
                            +tPTab.elemProps[tAtm2.chemType]["cova"];
                    tExtraD= tSens*tRad;
                    tRange[0] =   tRad -tExtraD;
                    if (tRange[0] < 0.2)
                    {
                        tRange[0] =0.2;
                    }
                    tRange[1] =  tRad + tExtraD;
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
    
    
    void MolGenerator::getBondingRangePairAtoms2(AtomDict          & tAtm1, 
                                                 AtomDict          & tAtm2, 
                                                 REAL              tSens,
                                                 PeriodicTable     & tPTab,
                                                 std::vector<REAL> & tRange)
    {
        tRange.clear();
        tRange.push_back(-1.0);
        tRange.push_back(-1.0);
        
        REAL comp1=0.0, comp2=0.0;
        
        if (tPTab.elemProps.find(tAtm1.chemType) !=tPTab.elemProps.end()
            && tPTab.elemProps.find(tAtm2.chemType) !=tPTab.elemProps.end())
        {
            if (tAtm1.formalCharge !=0 && tAtm2.formalCharge !=0)
            {
                
                if (tAtm1.formalCharge >0 && tAtm2.formalCharge <0)
                {
                    if (tPTab.elemProps[tAtm1.chemType].find("ionM+")
                            !=tPTab.elemProps[tAtm1.chemType].end()
                        && tPTab.elemProps[tAtm2.chemType].find("ionM-")
                            !=tPTab.elemProps[tAtm2.chemType].end())
                    {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["ionM+"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["ionM-"];
                    }
                    else
                    {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                    }
                }
                else if (tAtm1.formalCharge < 0 && tAtm2.formalCharge > 0)
                {
                    if (tPTab.elemProps[tAtm1.chemType].find("ionM-")
                            !=tPTab.elemProps[tAtm1.chemType].end()
                         && tPTab.elemProps[tAtm2.chemType].find("ionM+")
                            !=tPTab.elemProps[tAtm2.chemType].end() )
                    {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["ionM-"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["ionM+"];
                    }
                    else
                    {
                        comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                        comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                    }
                }
                else
                {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                }
            }
            else if (tAtm1.isMetal && tAtm2.isMetal)
            {
                if (tPTab.elemProps[tAtm1.chemType].find("ionM+")
                            !=tPTab.elemProps[tAtm1.chemType].end()
                         && tPTab.elemProps[tAtm2.chemType].find("ionM+")
                            !=tPTab.elemProps[tAtm2.chemType].end())
                {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["ionM+"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["ionM+"];
                }
                else
                {
                    comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                    comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                }
            }
            else
            {
                comp1 = tPTab.elemProps[tAtm1.chemType]["cova"];
                //std::cout << "atom 1 " << tAtm1.chemType << std::endl;
                //std::cout << "inside comp1 " << comp1 << std::endl;
                comp2 = tPTab.elemProps[tAtm2.chemType]["cova"];
                //std::cout << "atom 2 " << tAtm2.chemType << std::endl;
                //std::cout << "comp2 " << comp2 << std::endl;
            }
            
            REAL tRad, tExtraD;
            tRad  =  comp1 + comp2;      
            tExtraD= tSens*tRad;
            tRange[0] =   tRad -tExtraD;
            if (tRange[0] < 0.2)
            {
                tRange[0] =0.2;
            }
            tRange[1] =  tRad + tExtraD;
            
            std::cout <<  "comp1 " << comp1 << std::endl;
            std::cout <<  "comp2 " << comp2 << std::endl;
            std::cout << "tRad " << tRad << std::endl;
            std::cout << "tSen " << tSens << std::endl;
            std::cout << "tExtraD " << std::endl;
           
        }
        else 
        {
           std::cout << "Bug! At least one of covalent radius for "
                     << tAtm1.id << " or " << tAtm2.id 
                     << " is not defined the internal periodic table" 
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
            //std::cout << "a bond between " << refAtoms[tIdxAtm1].id << " and "
            //          << refAtoms[tIdxAtm2].id << " is added to the b_list " 
            //          << std::endl << "Its bond length is " << aBond.value 
            //          << std::endl;
        }
    }
    
    void MolGenerator::setOneUniqueBondCell(int tIdxAtm1, int tIdxAtm2, 
                                            REAL rD)
    {
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
                                                REAL alpha, REAL beta, REAL gamma)
    {   
        if (tCoord1.size() == tCoord2.size() && tCoord1.size()==3)
        {
            std::vector<REAL> deltX;
            for (unsigned i=0; i <3; i++ )
            {
                deltX.push_back((tCoord2[i]-tCoord1[i]));
            }
            /*
            std::cout << "a2 " << a*a << std::endl;
            std::cout << "beta " << beta << std::endl;
            std::cout << "c " << c << std::endl;
            std::cout << "a " << a << std::endl;
            std::cout << "PI180 " << PI180 << std::endl;
            std::cout << "3.141592653589/180.0 " <<  3.141592653589/180.0 << std::endl;
            std::cout << "cos(beta*PI180) " << cos(beta*PI180) << std::endl;
            std::cout << "c*a*cos(beta*PI180) " << c*a*cos(beta*PI180) << std::endl;
             */
            return sqrt(a*a*deltX[0]*deltX[0] + b*b*deltX[1]*deltX[1] + c*c*deltX[2]*deltX[2]
                    + 2*b*c*deltX[1]*deltX[2]*cos(alpha*PI180)
                    + 2*c*a*deltX[2]*deltX[0]*cos(beta*PI180)
                    + 2*a*b*deltX[0]*deltX[1]*cos(gamma*PI180));
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
    
    void MolGenerator::getUniqAngles(std::vector<CrystInfo>::iterator tCryst)
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
                        (int) allAtoms[i2].connAtoms.size())
                    {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                    }
                    else
                    {
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
                    aAng.value  = getAngleValueFromFracCoords(allAtoms[i],
                                    allAtoms[i1], allAtoms[i2],
                                    tCryst->itsCell->a, tCryst->itsCell->b, 
                                    tCryst->itsCell->c, tCryst->itsCell->alpha,
                                    tCryst->itsCell->beta, tCryst->itsCell->gamma);
                    
                    aAng.sigValue     = 3.0;
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
    
        
    
    void MolGenerator::getUniqAngleMols(Molecule    & tMol,
                                        std::vector<CrystInfo>::iterator tCryst)
    {
        std::map<std::string, int> aMM;
        std::vector<REAL> aV;
        
        for(int i=0; i < (int)tMol.atoms.size(); i++)
        {   
            if (tMol.atoms[i].isInPreCell)
            {
                for (int j=0; j < (int)tMol.atoms[i].connAtoms.size(); j++)
                {
                    for (int k=j+1; k < (int)tMol.atoms[i].connAtoms.size(); k++)
                    {
                        
                        int i1 = tMol.atoms[i].connAtoms[j];
                        int i2 = tMol.atoms[i].connAtoms[k];
                        if (inBonds(i, i1, tMol.bonds) || inBonds(i, i2, tMol.bonds))
                        {
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
                            if (aMM.find(str1) ==aMM.end() && aMM.find(str2) == aMM.end())
                            {
                                aMM[str1] =1;
                                AngleDict aAng;
                                aAng.anchorID  = tMol.atoms[i].id;
                                aAng.anchorPos = i;
                                aAng.atoms.push_back(i);  // this is for individual molecule output
                                aAng.atomChemTypes.push_back(tMol.atoms[i].chemType);
                                aAng.atomsCodClasses.push_back(tMol.atoms[i].codClass); // this for the overall file
                       
                                if ((int) tMol.atoms[i1].connAtoms.size() >=
                                    (int) tMol.atoms[i2].connAtoms.size())
                                {
                                    aAng.atoms.push_back(i1);
                                    aAng.atoms.push_back(i2);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                }
                                else
                                {
                                    aAng.atoms.push_back(i2);
                                    aAng.atoms.push_back(i1);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                    aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                }
                    
                                aAng.value  = getAngleValueFromFracCoords(tMol.atoms[i],
                                              tMol.atoms[i1], tMol.atoms[i2],
                                              tCryst->itsCell->a, tCryst->itsCell->b, 
                                              tCryst->itsCell->c, tCryst->itsCell->alpha,
                                              tCryst->itsCell->beta, tCryst->itsCell->gamma);
                                // std::cout << "new key " << aAng.value << std::endl;
                               aAng.sigValue     = 3.0;
                               aAng.numCodValues = 0;
                               tMol.angles.push_back(aAng);
                               aV.push_back(aAng.value);
                            }
                            else
                            {
                                REAL rD= getAngleValueFromFracCoords(tMol.atoms[i],
                                            tMol.atoms[i1], tMol.atoms[i2],
                                            tCryst->itsCell->a, tCryst->itsCell->b, 
                                            tCryst->itsCell->c, tCryst->itsCell->alpha,
                                            tCryst->itsCell->beta, tCryst->itsCell->gamma);
                                // std::cout << "rD " << rD << std::endl;
                                if (!inVectABS(aV, rD, 0.0001))
                                {
                                    aMM[str1] =1;
                                    AngleDict aAng;
                                    aAng.anchorID  = tMol.atoms[i].id;
                                    aAng.anchorPos = i;
                                    aAng.atoms.push_back(i);  // this is for individual molecule output
                                    aAng.atomChemTypes.push_back(tMol.atoms[i].chemType);
                                    aAng.atomsCodClasses.push_back(tMol.atoms[i].codClass); // this for the overall file
                       
                                    if ((int) tMol.atoms[i1].connAtoms.size() >=
                                        (int) tMol.atoms[i2].connAtoms.size())
                                    {
                                        aAng.atoms.push_back(i1);
                                        aAng.atoms.push_back(i2);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                    }
                                    else
                                    {
                                        aAng.atoms.push_back(i2);
                                        aAng.atoms.push_back(i1);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i2].chemType);
                                        aAng.atomChemTypes.push_back(tMol.atoms[i1].chemType);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i2].codClass);
                                        aAng.atomsCodClasses.push_back(tMol.atoms[i1].codClass);
                                    }
                                    aAng.value  = rD;
                                    // std::cout << "new value " << aAng.value << std::endl;
                                    aAng.sigValue     = 3.0;
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
    
    REAL MolGenerator::getAngleValueFromFracCoords(AtomDict  & tAtCen,
                                                   AtomDict  & tAt1, 
                                                   AtomDict  & tAt2,
                                                   REAL a, REAL b, REAL c, 
                                                   REAL alpha, REAL beta, REAL gamma)
    {
        
        
        if (tAtCen.fracCoords.size()==3 && tAtCen.fracCoords.size()== tAt1.fracCoords.size()
              && tAtCen.fracCoords.size()== tAt2.fracCoords.size())
        {
            std::vector<REAL> tV1, tV2;
            for (unsigned i=0; i < 3; i++)
            {
                tV1.push_back(tAt1.fracCoords[i] - tAtCen.fracCoords[i]);
                tV2.push_back(tAt2.fracCoords[i] - tAtCen.fracCoords[i]);
            }
            
            REAL lenTv1=getBondLenFromFracCoords(tAtCen.fracCoords, tAt1.fracCoords,
                                                 a, b, c, alpha, beta, gamma);
            REAL lenTv2=getBondLenFromFracCoords(tAtCen.fracCoords, tAt2.fracCoords,
                                                 a, b, c, alpha, beta, gamma);
            //std::cout << "lenTv1 " << lenTv1 << std::endl;
            //std::cout << "lenTv2 " << lenTv2 << std::endl;
            
            if (lenTv1 > 0.000001 && lenTv2 > 0.000001)
            {
                REAL coF=  (a*a*tV1[0]*tV2[0] + b*b*tV1[1]*tV2[1] + c*c*tV1[2]*tV2[2]
                           + b*c*(tV1[1]*tV2[2]+tV1[2]*tV2[1])*cos(alpha*PI180)
                           + c*a*(tV1[2]*tV2[0]+tV1[0]*tV2[2])*cos(beta*PI180)
                           + a*b*(tV1[0]*tV2[1]+tV1[1]*tV2[0])*cos(gamma*PI180))/(lenTv1*lenTv2);
            
                //std::cout << "coF " << coF << std::endl;
                // The float point error : ABS might > 1.0 or < -1.0
                if (coF >=1.0)
                {
                    coF-=0.0000000001;
                }
                else if (coF <=-1.0) 
                {
                    coF +=0.0000000001;
                }
                return acos(coF);
            }
            else
            {
                std::cout << "There is at least one pair of atoms overlapped " << std::endl;
                return 0.0;
            }
        }
        else
        {
            std::cout << "Error: check atom coordinate dimensions " << std::endl;
            return 0.0;
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
        
        
        std::map<unsigned, std::vector<int> >  aMoleculesInCell;
        // get molecules in a unit cell
        for (unsigned i=0; i < classNum.size(); i++)
        {
            aMoleculesInCell[classNum[i]].push_back(i);
        }
        
        std::cout << "Number of molecules is after EQ classification " 
                  << aMoleculesInCell.size() << std::endl;
        
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
                                     std::vector<CrystInfo>::iterator  tCryst)
    {
        
        //std::cout << "Number of molecules after linked equiv class: " 
        //          << moleculesInCell.size() << std::endl;
        int nMol =0;
        for (std::map<unsigned, std::vector<int> >::iterator iMol=moleculesInCell.begin(); 
                  iMol !=moleculesInCell.end(); iMol++)
        {
            // Check if these molecule with atom index only contain atoms in ASU.
            // If yes, materialize the molecule.
            
            // check if there is any atom with occupancy less than 0.95 
            // and then set individual molecules
            //if (isASUAtomInMol(iMol))
            //{
                
            
                Molecule aMol;
                aMol.seriNum = nMol;
                nMol++;
                // std::cout << "Its serial number " 
                //           << aMol.seriNum << std::endl;
                std::map<int, int> tAtmMap;
                int i=0;
                for (std::vector<int>::iterator iAt=iMol->second.begin(); 
                       iAt !=iMol->second.end(); iAt++)
                {
                    tAtmMap[*iAt]=i;
                    aMol.atoms.push_back(allAtoms[*iAt]);
                    i++;
                }
            
                // Change the connection index for atoms in that molecule.
                i=0;
                for (std::vector<AtomDict>::iterator iAt=aMol.atoms.begin();
                       iAt !=aMol.atoms.end(); iAt++)
                {
                    iAt->seriNum=i;
                    i++;
                    for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                           iCo !=iAt->connAtoms.end(); iCo++)
                    {
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
               getAtomTypeOneMol(aMol);
               
               // getUniqueBondsMols(aMol, tCryst);
               getAllBondsInOneMol(aMol, tCryst);
               
               /*
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
               */    
            
               std::string  aMolErrInfo;    
               bool aPass = validateMolecule(aMol, tPTab, aMolErrInfo);
               if ( aPass && aMol.bonds.size() >0 )
               {
                   //getAtomTypeOneMol(aMol);
                   getUniqAngleMols(aMol, tCryst);
                   allMolecules.push_back(aMol);
                   std::cout << "Molecule " << aMol.seriNum << " is included "
                             << std::endl << "It contains " << aMol.atoms.size() 
                             << " atoms" << std::endl; 
               }
               else
               {
                   std::cout << aMolErrInfo << std::endl;
                   std::cout << "Molecule " << aMol.seriNum 
                             << " is rejected " << std::endl;
               }
            //}
        }
        std::cout << "Number of molecules in allMolecules after validation "
                  << allMolecules.size() << std::endl;
      
    }
    
    bool MolGenerator::checkAtomOcp(Molecule    & tMol,
                                    std::string & tErrInfo)
    {
        
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin(); 
                        iAt !=tMol.atoms.end(); iAt++)
        {
           //std::cout << "Atom " << iAt->id << " occupancy "
           //           << iAt->ocp << std::endl;
            
            if (iAt->ocp < 0.99)
            {
                tErrInfo = "Atom " + iAt->id + " has an occupancy "
                          + RealToStr(iAt->ocp) + ", small than 1.0!";
                
                return false;
            }
            
        }
        
        return true;
        
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
    
    void MolGenerator::deleteNonCenCellMols(std::map<unsigned, std::vector<int> >  
                                   & tMoleculesInCell)
    {
        std::map<unsigned, bool> tMolMap;
        for (std::map<unsigned, std::vector<int> >::iterator iM=tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++)
        {
            bool tIn=false;
            for (std::vector<int>::iterator iAt=iM->second.begin();
                    iAt !=iM->second.end(); iAt++)
            {
                if (allAtoms[*iAt].sId=="555")
                {
                    tIn=true;
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
        
        for (std::map<unsigned, bool>::iterator iM=tMolMap.begin();
                iM !=tMolMap.end(); iM++)
        {
            if (!iM->second)
            {
                tMoleculesInCell.erase(iM->first);
            }
        }
       
        
        unsigned iMol=1;
        for (std::map<unsigned, std::vector<int> >::iterator iM=tMoleculesInCell.begin();
                iM !=tMoleculesInCell.end(); iM++)
        {
            for (std::vector<int>::iterator iAt=iM->second.begin();
                    iAt !=iM->second.end(); iAt++)
            {
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
    
    void MolGenerator::deleteNonASUAtomCellMols(std::map<unsigned,std::vector<int> >& tMoleculesInCell)
    {
        
        
        std::map<unsigned, bool> tMolMap;
        for (std::map<unsigned, std::vector<int> >::iterator iM=tMoleculesInCell.begin();
                iM != tMoleculesInCell.end(); iM++)
        {
            
            
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
        
        for (std::map<unsigned, bool>::iterator iM=tMolMap.begin();
                iM !=tMolMap.end(); iM++)
        {
            if (!iM->second)
            {
                tMoleculesInCell.erase(iM->first);
            }
        }
       
        
        unsigned iMol=1;
        for (std::map<unsigned, std::vector<int> >::iterator iM=tMoleculesInCell.begin();
                iM !=tMoleculesInCell.end(); iM++)
        {
            for (std::vector<int>::iterator iAt=iM->second.begin();
                    iAt !=iM->second.end(); iAt++)
            {
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
    
    void MolGenerator::checkAtomElementID(std::vector<AtomDict> & tAtoms)
    {
        PeriodicTable aPTab;
        
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            if (iAt->chemType.empty() || 
                aPTab.elements.find(iAt->chemType) == aPTab.elements.end())
                     
            {
                if (!iAt->id.empty())
                {
                    if (iAt->id.size() >=2)
                    {
                        if(aPTab.elements.find(iAt->id.substr(0,2)) 
                                !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id.substr(0,2);
                        }
                        else if (aPTab.elements.find(iAt->id.substr(0,1))
                                  !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id.substr(0,1);
                        }
                    }
                    else if(iAt->id.size()==1)
                    {
                        if (aPTab.elements.find(iAt->id)
                              !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id;
                        }
                    }
                }
                else
                {
                    std::cout << "What is atom " << iAt->id << "'s element type "
                              << std::endl;
                    exit(1);
                }
            }
        }
        
    }
    
    bool MolGenerator::colidAtom(AtomDict& tAtom, 
                                 std::vector<AtomDict>& tRefAtoms)
    {
        for (std::vector<AtomDict>::iterator iAt=tRefAtoms.begin();
                iAt !=tRefAtoms.end(); iAt++)
        {
            REAL da = fabs(tAtom.fracCoords[0]-iAt->fracCoords[0]);
            REAL db = fabs(tAtom.fracCoords[1]-iAt->fracCoords[1]);
            REAL dc = fabs(tAtom.fracCoords[2]-iAt->fracCoords[2]);
            if (da < 0.01 && db < 0.01 && dc < 0.01)
            {
                return true;
            }
        }
        
        return false;
    }
    
    bool MolGenerator::colidAtom(std::vector<REAL>    & tFracX, 
                                 std::vector<AtomDict>& tRefAtoms)
    {
        for (std::vector<AtomDict>::iterator iAt=tRefAtoms.begin();
                iAt !=tRefAtoms.end(); iAt++)
        {
            REAL da = fabs(tFracX[0]-iAt->fracCoords[0]);
            REAL db = fabs(tFracX[1]-iAt->fracCoords[1]);
            REAL dc = fabs(tFracX[2]-iAt->fracCoords[2]);
            if (da < 0.005 && db < 0.005 && dc < 0.005)
            {
                return true;
            }
        }
        
        return false;
    }
    
    
    
    bool MolGenerator::isASUAtomInMol(std::map<unsigned,std::vector<int> >::iterator tMol)
    {
        for (std::vector<int>::iterator iAt=tMol->second.begin(); 
                     iAt !=tMol->second.end(); iAt++)
        {
            if (allAtoms[*iAt].isInPreCell)
            {
                return true;
            }
        }
        
        return false;             
    }
    
    bool MolGenerator::connMetal(std::vector<int>     &  tIdxs, 
                                 std::vector<AtomDict>& tAtoms)
    {
        for (std::vector<int>::iterator iIdx=tIdxs.begin();
                iIdx !=tIdxs.end(); iIdx++)
        {
            if (tAtoms[*iIdx].isMetal)
            {
                return true;
            }
        }
        
        return false;
    }
    
    bool MolGenerator::validateBonds(std::vector<BondDict>::iterator tBo, 
                                     std::string & tErrInfo,
                                     PeriodicTable & tPTab)
    {
        // Bond properties (number of connection) check
        std::cout << "Bond between atoms " << allAtoms[tBo->atomsIdx[0]].id
                  << " and " << allAtoms[tBo->atomsIdx[1]].id << std::endl;
        
        if ((allAtoms[tBo->atomsIdx[0]].connAtoms.size()==1 && 
              allAtoms[tBo->atomsIdx[1]].connAtoms.size()==1))
        {
            tErrInfo = "Both atoms connect to only one atom for atom" 
                       + allAtoms[tBo->atomsIdx[0]].id  + " and "
                       + allAtoms[tBo->atomsIdx[1]].id ;
            return false;
        }
        else if ((allAtoms[tBo->atomsIdx[0]].connAtoms.size() >4 && 
                  !allAtoms[tBo->atomsIdx[0]].isMetal)
                || (allAtoms[tBo->atomsIdx[1]].connAtoms.size() > 4 && 
                  !allAtoms[tBo->atomsIdx[1]].isMetal))
        {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id 
                      + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[0]].connAtoms.size()) 
                      + "Atom " + allAtoms[tBo->atomsIdx[1]].id 
                      + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[1]].connAtoms.size()); 
            
            return false;
        }
        // C
        else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("C")==0 
                  && allAtoms[tBo->atomsIdx[0]].connAtoms.size() <3) 
                  ||
                  (allAtoms[tBo->atomsIdx[1]].chemType.compare("C")==0 
                  && allAtoms[tBo->atomsIdx[1]].connAtoms.size() <3))
        {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id 
                      + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[0]].connAtoms.size()) 
                      + "Atom " + allAtoms[tBo->atomsIdx[1]].id 
                      + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[1]].connAtoms.size()); 
            return false;
        }
        // B
        else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("B")==0 
                  && allAtoms[tBo->atomsIdx[0]].connAtoms.size() <2) 
                  ||
                  (allAtoms[tBo->atomsIdx[1]].chemType.compare("B")==0 
                  && allAtoms[tBo->atomsIdx[1]].connAtoms.size() <2))
        {
            tErrInfo = "Atom " + allAtoms[tBo->atomsIdx[0]].id 
                      + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[0]].connAtoms.size()) 
                      + "Atom " + allAtoms[tBo->atomsIdx[1]].id 
                      + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[1]].connAtoms.size());  
            return false;
        }
        else if ((allAtoms[tBo->atomsIdx[0]].chemType.compare("H")==0 
                  && allAtoms[tBo->atomsIdx[0]].connAtoms.size() >1) 
                  ||
                  (allAtoms[tBo->atomsIdx[1]].chemType.compare("H")==0 
                  && allAtoms[tBo->atomsIdx[1]].connAtoms.size() > 1))
        {
            tErrInfo ="Reject the molecule! Atom " + allAtoms[tBo->atomsIdx[0]].id 
                      + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[0]].connAtoms.size()) 
                      + "Atom " + allAtoms[tBo->atomsIdx[1]].id 
                      + " serial number " + IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                      + " has connections " + IntToStr((int)allAtoms[tBo->atomsIdx[1]].connAtoms.size());  
            return false;
        } 
       
        // bond length check again, this time we use both high and 
        // low bounds of a bond length
        
        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 =0.25;
        REAL covalent_sensitivity2 =0.50;
        
        std::vector<REAL> linkRange;
                
        
        if ((!allAtoms[tBo->atomsIdx[0]].isMetal) && (!allAtoms[tBo->atomsIdx[1]].isMetal))
        {
            covalent_sensitivity = covalent_sensitivity1;
        }
        else
        {
            covalent_sensitivity = covalent_sensitivity2;
        }
        getBondingRangePairAtoms2(allAtoms[tBo->atomsIdx[0]], allAtoms[tBo->atomsIdx[1]],
                                  covalent_sensitivity, tPTab, linkRange);
        
        if (tBo->value <linkRange[0] && tBo->value > linkRange[1])
        {
            tErrInfo = "Bond between " + allAtoms[tBo->atomsIdx[0]].id 
                     + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[0]].seriNum)
                     + " and " + allAtoms[tBo->atomsIdx[1]].id 
                     + " serial number " +  IntToStr(allAtoms[tBo->atomsIdx[1]].seriNum)
                     + " is " + RealToStr(tBo->value)
                     + "It should be between " + RealToStr(linkRange[0]) + " and " 
                     + RealToStr(linkRange[1]);
            return false;
        }
        
        return true;
       
    }
    
    bool MolGenerator::validateBonds(std::vector<BondDict>::iterator tBo, 
                                     Molecule      & tMol,
                                     std::string & tErrInfo,
                                     PeriodicTable & tPTab)
    {
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
        if ((tMol.atoms[tBo->atomsIdx[0]].connAtoms.size()==1 && 
              tMol.atoms[tBo->atomsIdx[1]].connAtoms.size()==1))
        {
            tErrInfo = "Both atoms connect to only one atom for atom" 
                       + tMol.atoms[tBo->atomsIdx[0]].id  + " and "
                       + tMol.atoms[tBo->atomsIdx[1]].id ;
            return false;
        }
        
        // bond length check again, this time we use BOTH(!) the high and 
        // low bounds of a bond length
        
        REAL covalent_sensitivity;
        REAL covalent_sensitivity1 =0.25;
        REAL covalent_sensitivity2 =0.50;
        
        std::vector<REAL> linkRange;
                
        
        if ((!tMol.atoms[tBo->atomsIdx[0]].isMetal) && (!tMol.atoms[tBo->atomsIdx[1]].isMetal))
        {
            covalent_sensitivity = covalent_sensitivity1;
        }
        else
        {
            covalent_sensitivity = covalent_sensitivity2;
        }
        getBondingRangePairAtoms2(tMol.atoms[tBo->atomsIdx[0]], tMol.atoms[tBo->atomsIdx[1]],
                                  covalent_sensitivity, tPTab, linkRange);
        
        if (tBo->value <linkRange[0] && tBo->value > linkRange[1])
        {
            tErrInfo = "Bond between " + tMol.atoms[tBo->atomsIdx[0]].id 
                     + " serial number " +  IntToStr(tMol.atoms[tBo->atomsIdx[0]].seriNum)
                     + " and " + tMol.atoms[tBo->atomsIdx[1]].id 
                     + " serial number " +  IntToStr(tMol.atoms[tBo->atomsIdx[1]].seriNum)
                     + " is " + RealToStr(tBo->value)
                     + "It should be between " + RealToStr(linkRange[0]) + " and " 
                     + RealToStr(linkRange[1]);
            return false;
        }
        
        return true;
       
    }
    
    bool MolGenerator::validateAtomLinks(Molecule & tMol, 
                                         PeriodicTable & tPTab,
                                         std::string   & tErrInfo) 
                                         
    {
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                iAt !=tMol.atoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
            if (iAt->connAtoms.size() >4 && !iAt->isMetal
                && !connMetal(iAt->connAtoms, tMol.atoms))
            {
                tErrInfo = "Reject Molecule! Atom " + iAt->id 
                           + " serial number " +  IntToStr(iAt->seriNum)
                           + " has connections " + IntToStr((int)iAt->connAtoms.size()); 
                return false;
            }
            // C
            else if (iAt->chemType.compare("C")==0 && iAt->connAtoms.size() <3)
            {
                tErrInfo =   "Reject Molecule! Atom " + iAt->id 
                              + " serial number " +  IntToStr(iAt->seriNum)
                              + " has connections " + IntToStr((int)iAt->connAtoms.size());
                return false;
            }
            // O
            else if (iAt->chemType.compare("O")==0 && iAt->connAtoms.size() >2)     
            {
                bool tFind=false;
            
                for (std::vector<int>::iterator iCo=iAt->connAtoms.begin();
                        iCo !=iAt->connAtoms.end(); iCo++)
                {
                    if (tMol.atoms[*iCo].isMetal)
                    {
                        tFind = true;
                        break;
                    }
                }
            
                if (!tFind)
                {
                    tErrInfo = "Reject Molecule: covalent bonded O atom has more"
                               " than 2 bonds!   Atom " + iAt->id 
                             + " serial number " +  IntToStr(iAt->seriNum)
                             + " has connections " + IntToStr((int)iAt->connAtoms.size()); 
                             
                    return false;
                }
            }
            // B
            else if (iAt->chemType.compare("B")==0 
                  && iAt->connAtoms.size() <2) 
                  
            {
                tErrInfo = "Reject molecule ! B Atom " + iAt->id 
                          + " serial number " +  IntToStr(iAt->seriNum)
                          + " has connections " + IntToStr((int)iAt->connAtoms.size());
                return false;
                          
            }
            // H 
            else if (iAt->chemType.compare("H")==0 && iAt->connAtoms.size() >1) 
                  
            {
                tErrInfo =  "Reject the molecule! H Atom " + iAt->id 
                          + " serial number "   +  IntToStr(iAt->seriNum)
                          + " has connections " +  IntToStr((int)iAt->connAtoms.size()); 
                return false;
            }
            // Halogen group
            else if (tPTab.elements[iAt->chemType]["group"] == 17  
                  && iAt->connAtoms.size() >1 
                  && !connMetal(iAt->connAtoms, tMol.atoms))     
            {
                tErrInfo ="Reject the molecule! Connection > 1 for Halogen Atom " 
                          + iAt->id 
                          + " serial number " +  IntToStr(iAt->seriNum)
                          + " has connections " + IntToStr((int)iAt->connAtoms.size()); 
                        
                return false;
            }
        }
        }
        
        return true;
        
    }
    
    bool MolGenerator::validateMolecule(Molecule      & tMol, 
                                        PeriodicTable & tPTab, 
                                        std::string   & tErrInfo)
    {
        
        if (!checkAtomOcp(tMol, tErrInfo))
        {
            return false;
        }
        if (!validateAtomLinks(tMol, tPTab, tErrInfo))
        {
            return false;
        }
        for (std::vector<BondDict>::iterator iBo=tMol.bonds.begin();
                        iBo !=tMol.bonds.end(); iBo++)
        {
            if (!validateBonds(iBo, tMol, tErrInfo, tPTab))
            {
                return false;
            }
        }
        
        return true;
        
    }
    
    void MolGenerator::getAtomTypeOneMol(Molecule& tMol)
    {
        //std::cout << "Number of atoms in this molecule is " 
        //          << tMol.atoms.size() << std::endl;
        
        CodClassify aCodSys(tMol.atoms);
        
        aCodSys.codAtomClassify2(2);
        
        tMol.atoms.clear();         
        for (std::vector<AtomDict>::iterator iAt=aCodSys.allAtoms.begin();
                iAt!=aCodSys.allAtoms.end(); iAt++)
        {
            tMol.atoms.push_back(*iAt);
        }
                
                
        for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                 iAt !=tMol.atoms.end(); iAt++)
        {
            if (iAt->isInPreCell)
            {
                std::cout << "Atom " << iAt->id << " has COD class id "
                                     << iAt->codClass << std::endl;
            }
        }
                        
    }
   
    void MolGenerator::getAtomTypeMols()
    {
        for(unsigned i=0; i < allMolecules.size(); i++)
        {
            getAtomTypeOneMol(allMolecules[i]);
        }
    }
    
    void MolGenerator::getOverallBondAndAngles()
    {
        std::map<std::string, std::vector<REAL> > aBM, aAM;
        
        for (std::vector<Molecule>::iterator iMol=allMolecules.begin();
                iMol !=allMolecules.end(); iMol++)
        {
            //std::cout << "Mol " << iMol->seriNum << " has " << iMol->bonds.size() 
            //          << " bonds " << std::endl;
            
            for (std::vector<BondDict>::iterator iBo=iMol->bonds.begin();
                    iBo !=iMol->bonds.end(); iBo++)
            {
                std::vector<sortMap3> tVec;
                struct sortMap3 tMap;
                tMap.key = iMol->atoms[iBo->atomsIdx[0]].codClass;
                tMap.val = iMol->atoms[iBo->atomsIdx[0]].chemType;
                tVec.push_back(tMap);
                tMap.key = iMol->atoms[iBo->atomsIdx[1]].codClass;
                tMap.val = iMol->atoms[iBo->atomsIdx[1]].chemType;
                tVec.push_back(tMap);
                std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                std::string tKey = tVec[0].key + "_" + tVec[1].key;
                
                if (aBM.find(tKey) == aBM.end())
                {
                    // std::cout << "new key " << tKey << std::endl
                    //          << "value is " << iBo->value << std::endl;
                    iBo->atoms.clear();
                    iBo->atomsCodClasses.clear();
                    for (std::vector<sortMap3>::iterator iAt=tVec.begin();
                           iAt !=tVec.end(); iAt++)
                    {
                        iBo->atoms.push_back(iAt->val);
                        iBo->atomsCodClasses.push_back(iAt->key);
                    }
                    bonds.push_back(*iBo);
                    aBM[tKey].push_back(iBo->value);
                }
                else if (!inVectABS(aBM[tKey], iBo->value, 0.00001))
                {
                    // std::cout << "new value " << iBo->value << std::endl
                    //          << "Key is " << tKey << std::endl;
                    iBo->atoms.clear();
                    iBo->atomsCodClasses.clear();
                    for (std::vector<sortMap3>::iterator iAt=tVec.begin();
                          iAt !=tVec.end(); iAt++)
                    {
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
            
            for(std::vector<AngleDict>::iterator iAn=iMol->angles.begin();
                    iAn !=iMol->angles.end(); iAn++)
            {
               
                std::vector<sortMap3> tVec;
                struct sortMap3 tMap;
                tMap.key = iMol->atoms[iAn->atoms[1]].codClass;
                tMap.val = iMol->atoms[iAn->atoms[1]].chemType;
                tVec.push_back(tMap);
                tMap.key = iMol->atoms[iAn->atoms[2]].codClass;
                tMap.val = iMol->atoms[iAn->atoms[2]].chemType;
                tVec.push_back(tMap);
                std::sort(tVec.begin(), tVec.end(), desSortMapKey3);
                std::string tKey = iMol->atoms[iAn->atoms[0]].codClass 
                                   + "_" + tVec[0].key + "_" + tVec[1].key;
                if(aAM.find(tKey)==aAM.end())
                {
                    iAn->atomChemTypes.clear();
                    iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].chemType);
                    iAn->atomsCodClasses.clear();
                    iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                    for (std::vector<sortMap3>::iterator iAt=tVec.begin();
                            iAt !=tVec.end(); iAt++)
                    {
                        iAn->atomChemTypes.push_back(iAt->val);
                        iAn->atomsCodClasses.push_back(iAt->key);
                    }
                    angles.push_back(*iAn);
                    aAM[tKey].push_back(iAn->value);
                }
                else if (!inVectABS(aAM[tKey], iAn->value, 0.0001))
                {
                    iAn->atomChemTypes.clear();
                    iAn->atomChemTypes.push_back(iMol->atoms[iAn->atoms[0]].chemType);
                    iAn->atomsCodClasses.clear();
                    iAn->atomsCodClasses.push_back(iMol->atoms[iAn->atoms[0]].codClass);
                    for (std::vector<sortMap3>::iterator iAt=tVec.begin();
                            iAt !=tVec.end(); iAt++)
                    {
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
    
    void MolGenerator::outTableMols(std::ofstream & tMolTabs,
                                 Molecule & tMol)
    {   
            if (tMol.atoms.size() >0)
            {
                tMolTabs << "data_mol_" <<IntToStr(tMol.seriNum +1) 
                         << std::endl << std::endl;
                        
                tMolTabs << "loop_" << std::endl
                          << "_chem_comp_atom.serial_number"  << std::endl
                          << "_chem_comp_atom.atom_id "       << std::endl
                          << "_chem_comp_atom.element_symbol" << std::endl
                          << "_chem_comp_atom.cod_type"       << std::endl;
                
                for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin(); 
                        iAt !=tMol.atoms.end(); iAt++)
                {
                    
                    tMolTabs << std::setw(6) << iAt->seriNum+1   
                              << std::setw(10) << iAt->id    
                              << std::setw(5) << iAt->chemType << "    "
                              << iAt->codClass << std::endl;
                }
                
                tMolTabs << std::endl;
            }
            else
            {
                std::cout << "There is no atoms in the molecule" << std::endl;
            }
            
            if (tMol.bonds.size() >0)
            {
                std::cout << std::endl;
                tMolTabs << "loop_" << std::endl  
                      << "_chem_comp_bond.bond_serial_number"  << std::endl
                      << "_chem_comp_bond.atom1_serial_number" << std::endl
                      << "_chem_comp_bond.atom2_serial_number" << std::endl
                      << "_chem_comp_bond.atom1_element_symbol" << std::endl
                      << "_chem_comp_bond.atom2_element_symbol" << std::endl
                      << "_chem_comp_bond.value_dist" << std::endl
                      << "_chem_comp_bond.is_in_same_ring" << std::endl;
                int nBo = 1;
                for (std::vector<BondDict>::iterator iBo=tMol.bonds.begin();
                        iBo !=tMol.bonds.end(); iBo++)
                {
                    std::string tStr;
                    if (iBo->isInSameRing)
                    {
                        tStr="Y";
                    }
                    else
                    {
                        tStr="N";
                    }
                   
                    tMolTabs << std::setw(6) << nBo
                             << std::setw(6) << iBo->atomsIdx[0] + 1 
                             << std::setw(6) << iBo->atomsIdx[1] + 1
                             << std::setw(4) << tMol.atoms[iBo->atomsIdx[0]].chemType
                             << std::setw(4) << tMol.atoms[iBo->atomsIdx[1]].chemType
                             << std::setw(10)<< iBo->value 
                             << std::setw(6) << tStr << std::endl;
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
                tMolTabs << "loop_" << std::endl
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
                    tMolTabs << std::setw(6)  << nAn
                              << std::setw(6)  << iAn->atoms[0] +1 
                              << std::setw(6)  << iAn->atoms[1] +1
                              << std::setw(6)  << iAn->atoms[2] +1
                              << std::setw(4)  << tMol.atoms[iAn->atoms[0]].chemType 
                              << std::setw(4)  << tMol.atoms[iAn->atoms[1]].chemType 
                              << std::setw(4)  << tMol.atoms[iAn->atoms[2]].chemType 
                              << std::setw(10) << iAn->value*PID180 << std::endl;
                    nAn++;
                }
                
                std::cout << std::endl;
            }
            
            std::cout << std::endl;  
    }
    
    void MolGenerator::outTableBAndA(FileName tBAndAFName)
    {
        std::cout << tBAndAFName << std::endl;

        
        std::ofstream aBAndAF(tBAndAFName);
        if (aBAndAF.is_open())
        {
            if (bonds.size() >0)
            {
                aBAndAF << "Atom1_COD_type        Atom2_COD_type          Atom1_chem_type     Atom2_chem_type     Bond_length       BondInRing"
                        << std::endl;
                for (std::vector<BondDict>::iterator iBo=bonds.begin();
                        iBo !=bonds.end(); iBo++)
                {
                    std::string tStr;
                    if (iBo->isInSameRing)
                    {
                        tStr="Y";
                    }
                    else
                    {
                        tStr="N";
                    }
                    
                    if (iBo->atomsCodClasses.size()==2)
                    {
                          aBAndAF <<  iBo->atomsCodClasses[0] << "\t"
                                  <<  iBo->atomsCodClasses[1] << "\t"
                                  <<  iBo->atoms[0] << "\t" << iBo->atoms[1] 
                                  << "\t"  << iBo->value << "\t"  
                                  <<  tStr << std::endl;
                    }
                }
                
                aBAndAF << std::endl;
            }
           
            
            if (angles.size() !=0)
            {
                aBAndAF << "Center_Atom_COD_type        Atom1_COD_type          Atom2_COD_type     " 
                        << "Center_Atom_chem_type      Atom1_chem_type     Atom2_chem_type     Angle_length"
                        << std::endl;
                for (std::vector<AngleDict>::iterator iAn=angles.begin();
                        iAn !=angles.end(); iAn++)
                {
                    aBAndAF << iAn->atomsCodClasses[0] << "\t" 
                            << iAn->atomsCodClasses[1] << "\t"
                            << iAn->atomsCodClasses[2] << "\t"
                            << iAn->atomChemTypes[0] << "\t"
                            << iAn->atomChemTypes[1] << "\t"
                            << iAn->atomChemTypes[2] << "\t"
                            << iAn->value*PID180 << std::endl;
                }
                aBAndAF << std::endl;
                
            }
           
            aBAndAF.close();
        }
    }
    
    void MolGenerator::outTables(FileName tOutName)
    {
        Name aFName(tOutName);
            
        std::vector<std::string> nameComps;
        StrTokenize(aFName, nameComps, '.');
        Name rootFName;
        for (unsigned jF=0; jF < nameComps.size()-1; jF++)
        {
            rootFName.append(nameComps[jF]);
        }
            
        if (rootFName.size() ==0)
        {
            rootFName.append("Current");
        }
        
        if (bonds.size() !=0)
        {   
            Name bondAndAngleFName(rootFName);
            bondAndAngleFName.append("_unique_bond_and_angles.txt");   
            outTableBAndA(bondAndAngleFName.c_str());
        }
        
        if (allMolecules.size() !=0)
        {
            Name allMolsFName(rootFName);
            allMolsFName.append("_all_mols.txt");
            std::ofstream aMolTable(allMolsFName.c_str());
            if (aMolTable.is_open())
            {
                for (unsigned iMol=0; iMol < allMolecules.size(); iMol++)
                {
                    outTableMols( aMolTable, allMolecules[iMol]);
                }
            }
            else
            {
                std::cout << allMolsFName << " can not be open for writing " 
                          << std::endl;
            }
        }
    }
}
