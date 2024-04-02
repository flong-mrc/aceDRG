
/* 
 * File:   GlobOpt.h
 * Author: flong
 *
 * last updated on February 13, 2013, 17:28 PM
 */

#include "GlobOpt.h"
#include "AllSystem.h"
#include "atomsTree.h"

namespace GO
{
    
    OptimSet::OptimSet():objValue(1.0e8)
    {
    }
    
    OptimSet::~OptimSet()
    {    
    }
    
    extern bool sortOptSetByValues(const OptimSet & tSetA, const OptimSet & tSetB)
    {
        return tSetA.objValue < tSetB.objValue;
    }
    
    FindGlobMin::FindGlobMin():lGlob(true),
                               gWorkMode(1),
                               lWorkMode(1),
                               workSpace(1),
                               lComp(100),
                               curObjValue(1.0e8),
                               preObjValue(1.0e8),
                               numOptimSets(50),
                               hasIniCoords(false),
                               itsFileNameRoot(NullString),
                               itsMonoRoot(NullString)                           
    {
    }
    
    FindGlobMin::FindGlobMin(const LIBMOL::AllSystem& aSystem, 
                             std::string    tFileNameeRoot,
                             std::string    tMonoRoot):lGlob(true),
                                                       gWorkMode(1),
                                                       lWorkMode(1),
                                                       workSpace(1),
                                                       lComp(100),
                                                       curObjValue(1.0e8),
                                                       preObjValue(1.0e8),
                                                       numOptimSets((int)aSystem.allAtoms.size())
    {
        for (std::vector<LIBMOL::AtomDict>::const_iterator iAt=aSystem.allAtoms.begin();
                iAt != aSystem.allAtoms.end(); iAt++)
        {
            
            allAtoms.push_back(*iAt);
        }
        
        for (std::vector<int>::const_iterator iHA=aSystem.allHAtomIdx.begin();
                iHA != aSystem.allHAtomIdx.end(); iHA++)
        {
            allHAtomIdx.push_back(*iHA);
        }
        
        for (std::vector<LIBMOL::BondDict>::const_iterator iBo=aSystem.allBonds.begin();
                iBo != aSystem.allBonds.end(); iBo++)
        {
            allBonds.push_back(*iBo);
        }
        
        for (std::vector<LIBMOL::AngleDict>::const_iterator iAn=aSystem.allAngles.begin();
                iAn !=aSystem.allAngles.end(); iAn++)
        {
            allAngles.push_back(*iAn);
        }
        
        for (std::vector<LIBMOL::TorsionDict>::const_iterator iTo=aSystem.allTorsions.begin();
                iTo !=aSystem.allTorsions.end(); iTo++)
        {
            allTorsions.push_back(*iTo);
        }
        
        for (std::vector<LIBMOL::ChiralDict>::const_iterator iCh =aSystem.allChirals.begin();
                iCh != aSystem.allChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<LIBMOL::PlaneDict>::const_iterator iPl =aSystem.allPlanes.begin();
                iPl !=aSystem.allPlanes.end(); iPl++)
        {
            allPlanes.push_back(*iPl);
        }
        
        for (std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> >::const_iterator iRm=aSystem.allRings.begin();
                iRm !=aSystem.allRings.end(); iRm++)
        {
            for (std::vector<LIBMOL::RingDict>::const_iterator aR = iRm->second.begin();
                    aR != iRm->second.end(); aR++)
            {
                allRings[iRm->first].push_back(*aR);
            }
        }
        
        for (std::vector<LIBMOL::RingDict>::const_iterator iRv=aSystem.allRingsV.begin();
                iRv !=aSystem.allRingsV.end(); iRv++)
        {
            allRingsV.push_back(*iRv);
        }
        
        hasIniCoords    = aSystem.hasCoords;
        
        
        itsFileNameRoot = tFileNameeRoot;
        itsMonoRoot     = tMonoRoot;
        
        if (aSystem.hasCoords)
        {
            numOptimSets = 5;
        }
        else
        {
            numOptimSets = (int)allAtoms.size();
        }
        
        std::cout << "The system for global optimization is initiated" << std::endl;
    }
    
    FindGlobMin::FindGlobMin(const LIBMOL::AllSystem& aSystem):gWorkMode(1),
                                                       lWorkMode(1),
                                                       workSpace(1),
                                                       lComp(100),
                                                       curObjValue(1.0e8),
                                                       preObjValue(1.0e8),
                                                       numOptimSets((int)aSystem.allAtoms.size())
    {
        for (std::vector<LIBMOL::AtomDict>::const_iterator iAt=aSystem.allAtoms.begin();
                iAt != aSystem.allAtoms.end(); iAt++)
        {
            allAtoms.push_back(*iAt);
        }
        
        for (std::vector<int>::const_iterator iHA=aSystem.allHAtomIdx.begin();
                iHA != aSystem.allHAtomIdx.end(); iHA++)
        {
            allHAtomIdx.push_back(*iHA);
        }
        
        for (std::vector<LIBMOL::BondDict>::const_iterator iBo=aSystem.allBonds.begin();
                iBo != aSystem.allBonds.end(); iBo++)
        {
            allBonds.push_back(*iBo);
        }
        
        for (std::vector<LIBMOL::AngleDict>::const_iterator iAn=aSystem.allAngles.begin();
                iAn !=aSystem.allAngles.end(); iAn++)
        {
            allAngles.push_back(*iAn);
        }
        
        for (std::vector<LIBMOL::TorsionDict>::const_iterator iTo=aSystem.allTorsions.begin();
                iTo !=aSystem.allTorsions.end(); iTo++)
        {
            allTorsions.push_back(*iTo);
        }
        
        for (std::vector<LIBMOL::ChiralDict>::const_iterator iCh =aSystem.allChirals.begin();
                iCh != aSystem.allChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<LIBMOL::PlaneDict>::const_iterator iPl =aSystem.allPlanes.begin();
                iPl !=aSystem.allPlanes.end(); iPl++)
        {
            allPlanes.push_back(*iPl);
        }
        
        for (std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> >::const_iterator iRm=aSystem.allRings.begin();
                iRm !=aSystem.allRings.end(); iRm++)
        {
            for (std::vector<LIBMOL::RingDict>::const_iterator aR = iRm->second.begin();
                    aR != iRm->second.end(); aR++)
            {
                allRings[iRm->first].push_back(*aR);
            }
        }
        
        for (std::vector<LIBMOL::RingDict>::const_iterator iRv=aSystem.allRingsV.begin();
                iRv !=aSystem.allRingsV.end(); iRv++)
        {
            allRingsV.push_back(*iRv);
        }
        
        if (aSystem.hasCoords)
        {
            numOptimSets = 5;
        }
        else
        {
            numOptimSets = (int)allAtoms.size();
        }
            
        // std::cout << "The system for global optimization is initiated" << std::endl;
    }
    
    FindGlobMin::FindGlobMin(const LIBMOL::AllSystem& aSystem,
                             const int     nOpt):gWorkMode(1),
                                                 lWorkMode(1),
                                                 workSpace(1),
                                                 lComp(100),
                                                 curObjValue(1.0e8),
                                                 preObjValue(1.0e8),
                                                 numOptimSets(nOpt)
    {
        for (std::vector<LIBMOL::AtomDict>::const_iterator iAt=aSystem.allAtoms.begin();
                iAt != aSystem.allAtoms.end(); iAt++)
        {
            
            allAtoms.push_back(*iAt);
        }
        
        for (std::vector<int>::const_iterator iHA=aSystem.allHAtomIdx.begin();
                iHA != aSystem.allHAtomIdx.end(); iHA++)
        {
            allHAtomIdx.push_back(*iHA);
        }
        
        for (std::vector<LIBMOL::BondDict>::const_iterator iBo=aSystem.allBonds.begin();
                iBo != aSystem.allBonds.end(); iBo++)
        {
            allBonds.push_back(*iBo);
        }
        
        for (std::vector<LIBMOL::AngleDict>::const_iterator iAn=aSystem.allAngles.begin();
                iAn !=aSystem.allAngles.end(); iAn++)
        {
            allAngles.push_back(*iAn);
        }
        
        for (std::vector<LIBMOL::TorsionDict>::const_iterator iTo=aSystem.allTorsions.begin();
                iTo !=aSystem.allTorsions.end(); iTo++)
        {
            allTorsions.push_back(*iTo);
        }
        
        for (std::vector<LIBMOL::ChiralDict>::const_iterator iCh =aSystem.allChirals.begin();
                iCh != aSystem.allChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<LIBMOL::PlaneDict>::const_iterator iPl =aSystem.allPlanes.begin();
                iPl !=aSystem.allPlanes.end(); iPl++)
        {
            allPlanes.push_back(*iPl);
        }
        
        for (std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> >::const_iterator iRm=aSystem.allRings.begin();
                iRm !=aSystem.allRings.end(); iRm++)
        {
            for (std::vector<LIBMOL::RingDict>::const_iterator aR = iRm->second.begin();
                    aR != iRm->second.end(); aR++)
            {
                allRings[iRm->first].push_back(*aR);
            }
        }
        
        for (std::vector<LIBMOL::RingDict>::const_iterator iRv=aSystem.allRingsV.begin();
                iRv !=aSystem.allRingsV.end(); iRv++)
        {
            allRingsV.push_back(*iRv);
        }
        hasIniCoords    = aSystem.hasCoords;    
        // std::cout << "The system for global optimization is initiated" << std::endl;
    }
    
    
    FindGlobMin::~FindGlobMin()
    {
    }
    
    
    // Initial atom positions 
    bool FindGlobMin::SetDefinedInitPositions()
    {
        // The coordinates linked with tree structures
        
        // LIBMOL::buildAtomTree tTreeTool;
        
        // tTreeTool.setAtomsMST(allAtoms, allBonds);
        
        //std::cout << " a tree-like structure is built " << std::endl;
         // convert all angle value and sigma to radius
        
        
        LIBMOL::TransCoords   tTransTool;
        
        // tTransTool.generateCoordTorsToCart(allAtoms,allBonds, allAngles, allTorsions, 
        // allRings, allPlanes, allChirals);
       
        bool tt=tTransTool.generateCoordTorsToCart3(allAtoms,allBonds, allAngles, allTorsions, 
                                            allRingsV, allPlanes, allChirals);
        
        
        std::cout << "Done a set of coordinates based on the tree structure " << std::endl;
        // std::vector<std::string> parts;
        //LIBMOL::StrTokenize(itsFileNameRoot, parts, '.');
        
        //std::string tFileName = parts[0] + "_tree";
        // std::cout << tFileName << std::endl;
        int aMode =1;
        LIBMOL::outPDB("tt.pdb", "UNK", allAtoms, aMode);
       
        
        /*
        for (std::vector<LIBMOL::AngleDict>::iterator iA=allAngles.begin();
                iA !=allAngles.end(); iA++)
        {
            iA->value      = iA->value*PI180;
            iA->valueST    = iA->valueST*PI180;
            iA->sigValue   = iA->sigValue*PI180;
            iA->sigValueST = iA->sigValueST*PI180;
            // std::cout << "Converted to " << iA->value << std::endl;
        } 
         */   
        // initial forces on atoms to zero;
        /*
        for (std::vector<LIBMOL::AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            int tDim = (int)iA->coords.size();
            for (int i=0; i < tDim; i++)
            {
                iA->forces.push_back(0.0);
            }
        }
        */
       
        // Initial optimization of the atom coordinates
        for (int i=0; i <(int)allAtoms.size(); i++)
        { 
            std::cout << "Atom coordinates for atom " << allAtoms[i].id << std::endl
                      << "X: " << allAtoms[i].coords[0] << "  Y:  " 
                      << allAtoms[i].coords[1] << "  Z: " 
                      << allAtoms[i].coords[2] << std::endl;
        }
        
        return tt;
        /*
        
        // Check 
        std::cout << "enter here " << std::endl;
        
        for (int i=0; i <(int)allAtoms.size(); i++)
        {   
            for (std::vector<int>::iterator iNB=allAtoms[i].connAtoms.begin();
                                iNB !=allAtoms[i].connAtoms.end(); iNB++)
            {
                //std::cout << "Atom coordinates for atom " << allAtoms[i].id << std::endl
                //          << "X: " << allAtoms[i].coords[0] << "  Y:  " 
                //          << allAtoms[i].coords[1] << "  Z: " 
                //          << allAtoms[i].coords[2] << std::endl;
                //std::cout << "Atom coordinates for NB atom " << allAtoms[*iNB].id << std::endl
                //          << "X: " << allAtoms[*iNB].coords[0] << "  Y:  " 
                //          << allAtoms[*iNB].coords[1] << "  Z: " 
                //          << allAtoms[*iNB].coords[2] << std::endl;
                LIBMOL::REAL tDist = LIBMOL::distanceV(allAtoms[i].coords, allAtoms[*iNB].coords);
                std::cout << "Bond length "  << tDistMat[i][*iNB] 
                          << " :  initial distance " << tDist << std::endl;
        
            }
        
        }
        
        for (std::vector<LIBMOL::AngleDict>::iterator iAN = allAngles.begin();
                iAN != allAngles.end(); iAN++)
        {
            std::vector<LIBMOL::REAL> v21, v31;
            for (int j=0; j < (int)allAtoms[iAN->atoms[0]].coords.size(); j++)
            {
                v21.push_back((allAtoms[iAN->atoms[1]].coords[j]-allAtoms[iAN->atoms[0]].coords[j]));
                v31.push_back((allAtoms[iAN->atoms[2]].coords[j]-allAtoms[iAN->atoms[0]].coords[j]));
            }
            LIBMOL::REAL tAngle = LIBMOL::getAngle2V(v21, v31);
            std::cout << "The angle between " << allAtoms[iAN->atoms[1]].id << " "
                               << allAtoms[iAN->atoms[0]].id << "  " << allAtoms[iAN->atoms[2]].id 
                               << " : "  << std::endl
                               << " Dictionary value " << iAN->value 
                               << " : initial value "    << tAngle*PID180 << std::endl;
        }
        
        LIBMOL::outPDB("2.pdb", "DDI", allAtoms);
       
        */
        
        
    }
    
    void FindGlobMin::SetRanDomInitPositions()
    {
        for(std::vector<LIBMOL::AtomDict>::iterator iAt=allAtoms.begin();
                iAt!=allAtoms.end(); iAt++)
        {
            for(std::vector<LIBMOL::REAL>::iterator iCo=iAt->coords.begin();
                    iCo !=iAt->coords.end(); iCo++)
            {
                LIBMOL::REAL r1 = LIBMOL::GetRand();
                
                if (r1 <=0.50)
                {
                    *iCo = LIBMOL::GetRand();
                }
                else
                {
                    *iCo = LIBMOL::GetRand();
                }           
            }
        }
    }
    
    void FindGlobMin::UpdatePositionsWithRandomShifts()
    {
        LIBMOL::REAL r_delta=0.01; 
        
        int tSize = (int)allOptimSets.size();
        std::cout << "We have " << tSize << " optimal sets " << std::endl;
        
        if (tSize%10==0)
        {
            r_delta=0.5*LIBMOL::GetRand();
        }
        else if (tSize%5==0)
        {
            r_delta=0.3*LIBMOL::GetRand();
        }
        else if (tSize%2==0)
        {
            r_delta=0.05*LIBMOL::GetRand();
        }
        else
        {
            r_delta=0.1*LIBMOL::GetRand();
        }
        
        // std::cout << "r_delta " << r_delta << std::endl;
        
        if(tSize <= 1)
        {
            
            for (std::vector<LIBMOL::AtomDict>::iterator iAt=allAtoms.begin();
                    iAt !=allAtoms.end(); iAt++)
            {
                for (std::vector<LIBMOL::REAL>::iterator iCo=iAt->coords.begin();
                         iCo !=iAt->coords.end(); iCo++)
                {
                    LIBMOL::REAL tShift=0.0, tRan=0.0;
                    tRan = LIBMOL::GetRand();
                    tShift = r_delta*tRan;
                    *iCo += tShift;
                    //std::cout << "tRan " << tRan << std::endl;
                    //std::cout << "shift " << tShift << std::endl;
                }
            }
        }
        else 
        {
            for (int i=0; i < (int)allAtoms.size(); i++)
            {
                for (int j=0; j <(int)allAtoms[i].coords.size(); j++)
                {
                    LIBMOL::REAL tShift = r_delta*LIBMOL::GetRand();
                    allAtoms[i].coords[j] = allOptimSets[0].atoms[i].coords[j] + tShift;
                    // std::cout << "shift " << tShift << std::endl;
                }
            }
        }
        
        
        
        
    }
    
    void FindGlobMin::PreIdealization()
    {
         for (std::vector<LIBMOL::AngleDict>::iterator iA=allAngles.begin();
                iA !=allAngles.end(); iA++)
         {
            //iA->value      = iA->value*PI180;
            iA->valueST    = iA->valueST*PI180;
            // iA->sigValue   = iA->sigValue*PI180;
            iA->sigValueST = iA->sigValueST*PI180;
            //std::cout << "Dict value Converted to " << iA->valueST << std::endl;
            //std::cout << "Dict sigma value converted to " << iA->sigValueST << std::endl;
         }   
        
         // A set of initial coordinates is done
         hasIniCoords = false;
         if (!hasIniCoords)
         {
             lGlob=SetDefinedInitPositions();
         }
        
        //else
        //{
        //    LIBMOL::outPDB("Init.pdb", "XXX", allAtoms);
        //}
   
        /*
        for (std::vector<LIBMOL::AtomDict>::iterator iA=allAtoms.begin();
              iA !=allAtoms.end();iA++)
        {
            std::cout << "For atom " << iA->id 
                      << ", its initial coords are: " <<std::endl;
            for (std::vector<LIBMOL::REAL>::iterator iCo=iA->coords.begin();
                    iCo !=iA->coords.end(); iCo++)
            {
                std::cout << *iCo << std::endl;
            }
            
        }
        
         */
        
       
        // do a Cartesian local optimization 
        lWorkMode = 1;
        workSpace = 1;
        // vdw only
        //std::cout << "VDW only \n";
        //lComp =101;
        //LocalMin();
        
        // std::cout << "Pre-stage " << std::endl;
        // lComp =100;
        // LocalMin();
        
        
        //std::cout << "Pre-stage 2 \n";
        //lComp =100;
        //LocalMin();
        //AddToOptim();
        // LIBMOL::outPDB("Pre-Idealized.pdb", "XXX", allAtoms);
        //std::cout << "Pre-Idealizations finished" << std::endl;
        // SelectBestOpt();
    }
    
    void FindGlobMin::ProIdealization()
    {
        std::cout << "Pro-idealization " << std::endl;
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            for (int j=0; j < (int)allAtoms[i].coords.size(); j++)
            {
                allAtoms[i].coords[j] = allOptimSets[0].atoms[i].coords[j];
            }
        }
        std::cout << "The starting objective value " 
                  << allOptimSets[0].objValue << std::endl;
        
        singleCompsIdealization();
        
        std::cout << " The final full idealization \n";
        lComp =100;
        LocalMin();
        int aMode =1;
        if (curObjValue < 1.5*allOptimSets[0].objValue)
        {
            AddToOptim();
            LIBMOL::outPDB("Pro-Idealized.pdb", "DDI", allAtoms, aMode);
        }
        else
        {
            LIBMOL::outPDB("Pro-Idealized.pdb", "DDI", allOptimSets[0].atoms, aMode);
        }
        std::cout << "Pro-Idealizations finished" << std::endl;
        
    }
    
    void FindGlobMin::singleCompsIdealization()
    {
        // correct bonds
        std::cout << "Bond only \n";
        lComp =1;
        LocalMin();
        // angle only
        std::cout << "Angle only \n";
        lComp =2;
        LocalMin();
        
        // Plane only
        std::cout << "Plane only \n";
        lComp =4;
        LocalMin();
        lComp=100;
        // vdw only
        //std::cout << "VDW only \n";
        //lComp =101;
        //LocalMin();
    }
    
    void FindGlobMin::LocalMin()
    {
        FindLocalMin  toolLocalMin;
        
        toolLocalMin.workMode  = lWorkMode;
        toolLocalMin.workSpace = workSpace;
        toolLocalMin.lComp     = lComp;
        
        toolLocalMin.Driver(allAtoms, allBonds, allAngles, allTorsions,
                            allRingsV, allPlanes, allChirals);
       
        preObjValue = curObjValue;
        curObjValue = toolLocalMin.curObjValue;
    }
    
    void FindGlobMin::AddToOptim()
    {
        OptimSet aOptSet;
        aOptSet.objValue = curObjValue;
        for (std::vector<LIBMOL::AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            aOptSet.atoms.push_back(*iA);
        }
        
        allOptimSets.push_back(aOptSet);
        if ((int)allOptimSets.size() !=0)
        {
            if (curObjValue < allOptimSets[0].objValue )
            {
                
                allOptimSets.push_back(aOptSet);
                std::sort(allOptimSets.begin(), allOptimSets.end(), sortOptSetByValues);
                std::cout << "The " << (int)allOptimSets.size() 
                          << "the OptimSet is stored " << std::endl;
                std::cout << "Its curObjValue " << aOptSet.objValue << std::endl;
            }
        }
        else
        {
          allOptimSets.push_back(aOptSet); 
          std::cout << "The " << (int)allOptimSets.size() 
                    << "th OptimSet is stored " << std::endl;
          std::cout << "Its curObjValue " << aOptSet.objValue << std::endl;
        }
        /*
        std::cout << "the lowest objective value is " 
                  << allOptimSets[0].objValue << std::endl;
        std::cout << "The first atom coords in this atom set are: " << std::endl;
        std::cout << "X = " << allOptimSets[0].atoms[0].coords[0] << std::endl
                  << "Y = " << allOptimSets[0].atoms[0].coords[1] << std::endl
                  << "Z = " << allOptimSets[0].atoms[0].coords[2] << std::endl;
        */
        
    }
    
    void FindGlobMin::SelectBestOpt()
    {
        if (allOptimSets.size() > 0)
        {
            allAtoms.clear();
            for (std::vector<LIBMOL::AtomDict>::iterator iA=allOptimSets[0].atoms.begin();
                    iA !=allOptimSets[0].atoms.end(); iA++)
            {
                allAtoms.push_back(*iA);
            }
        }
    }
    // 1 The global optimization scheme combing barrier tunneling 
    // and local minimization techniques
    
    void FindGlobMin::TunnellingAndMinDriver()
    { 
        
        lComp =100;
        std::cout << "Global min stage starts " << std::endl;
        // numOptimSets = (int)allAtoms.size();
        numOptimSets = 10;
        std::cout << "number of local min to be found " 
                  << numOptimSets << std::endl;
        
        
        do 
        {
            Tunnelling_Newton_EXP();
            LocalMin();
            
            AddToOptim();
            std::cout << "number of local min found " 
                      << (int)allOptimSets.size() 
                      << std::endl;
            if((int)allOptimSets.size()%20==0)
            {
                singleCompsIdealization();
            }
        }while((int)allOptimSets.size() < numOptimSets);
        
        
    }
    
    void FindGlobMin::TunnellingAndMinDriver(int nOpt)
    { 
        
        lComp =100;
        std::cout << "Global min stage starts " << std::endl;
        // numOptimSets = (int)allAtoms.size();
        numOptimSets = nOpt;
        std::cout << "number of local min to be found " 
                  << numOptimSets << std::endl;
        
        
        do 
        {
            Tunnelling_Newton_EXP();
            LocalMin();
            
            AddToOptim();
            std::cout << "number of local min found " 
                      << (int)allOptimSets.size() 
                      << std::endl;
            if((int)allOptimSets.size()%20==0)
            {
                singleCompsIdealization();
            }
        }while((int)allOptimSets.size() < numOptimSets);
        
        
    }
    
    //        Search postion x so that F(x) <= 0 along the Newton's
    //        direction, where F(x) is the objective function. 
    //        The solution found here will be used as an  
    //        initial point for a new local minimization.

    void FindGlobMin::Tunnelling_Newton_EXP()
    {
        
        int numAtoms = (int)allAtoms.size(); 
        int dim      = (int)allAtoms[0].coords.size();
        int numVars;
        
        if(workSpace == 1)
        {
            numVars = numAtoms*dim;
        }
        else if (workSpace == 2)
        {
            numVars = (int)allTorsions.size();
        }
        else
        {
            std::cout << "in which space are work  ? " << std::endl;
            exit(1);
        }
        
        
        if(!numVars)
        {
            std::cout << "The number of the variables is zero!"
                      << std::endl;
            exit(1);
        }
        
        int idx_error = 0;
        LIBMOL::REAL bal_value, bal_deriv_sq;

        //std::cout << "numVars is " << numVars << std::endl;
        //std::cout << "work space is " << workSpace << std::endl;
        
        int i, j;

        if(workSpace == 1)
        {
            // Newton's method in Cartesian space
            
            LIBMOL::REAL ** bal_deriv = new LIBMOL::REAL * [numAtoms];
            for (i = 0; i < numAtoms; i++)
            {
                bal_deriv [i] = new LIBMOL::REAL [dim];
                for (j = 0; j < dim; j++)
                {
                    bal_deriv[i][j] = 1.0;
                }
            }  
            
            LIBMOL::REAL ** position_m;
            position_m = new LIBMOL::REAL * [numAtoms];
            for (i = 0; i < numAtoms; i++)
            {
                position_m [i] = new LIBMOL::REAL  [dim];
                for (j = 0; j < dim; j++)
                {
                    position_m[i][j] = 0.0;
                }
            }
           
            int      N_max1 = 10;
            int      N_max2 = 200;
            int      i_loop_inner, i_loop_outer;

            LIBMOL::REAL belta_adapt = 0.0;
         
            i_loop_outer  = 1;
            
            do   
            {  
                // Initiate a random starting point
                UpdatePositionsWithRandomShifts();                  
                
                
                SetBalanceValueAndDeriv_EXP(idx_error, belta_adapt, bal_value, 
                                            bal_deriv, position_m);

                // std::cout << "init balance function value " << bal_value  << std::endl;
                // std::cout << "idx_error is " <<    idx_error << endl;
  
                // Starting tunnelling process along the selected direction
                // according to the Newton method 
           
                //belta_adapt = 0.2; 
                if(bal_value >= 1.0e-6 && !idx_error )
                {
               
                    i_loop_inner = 0;           
                    do 
                    {
                        bal_deriv_sq = 0.0;
                        for (i = 0; i < numAtoms; i++)
                        {
                            for (j = 0; j < dim; j++)
                            {
                                bal_deriv_sq+= 
                                bal_deriv[i][j]*bal_deriv[i][j];
                            }
                        }
                        
                        if (bal_deriv_sq < 1.0 )
                        {
                            bal_deriv_sq =  sqrt(bal_deriv_sq);
                        }
     
                        for (i = 0; i < numAtoms; i++)
                        {
                            //  std::cout << " For atom " << i+1 << endl; 
                            for (j = 0; j < dim; j++)
                            {
                                position_m[i][j]       = allAtoms[i].coords[j];
                                allAtoms[i].coords[j]  = allAtoms[i].coords[j]
                                       -bal_value*bal_deriv[i][j]/bal_deriv_sq;

                                // cout << "Deriv  in " << j+1 
                                //     << " is " << bal_deriv[i][j] << endl; 
                                //cout << "Deriv_sq  in " << bal_deriv_sq << endl;
                           
                                //cout << "Previous coord in " << j+1 
                                //     << " is " << position_m[i][j] << endl;
                                //cout << "Updated coord in " << j+1 
                                //     << " is " << atoms[i].coords[j] << endl;
                            }
                        }
                        
                        SetBalanceValueAndDeriv_EXP(idx_error, belta_adapt, bal_value, 
                                                    bal_deriv, position_m);
                  
                        i_loop_inner++;
                        //  cout << " at " <<  i_loop_inner << "th inner loop " << endl;
                        //  cout << "balance function value " << bal_value  << endl;
                        // cout << "idx_error is " <<    idx_error << endl;  
                        // cout << "Continue ? " << endl;
                        // cin.get();
                    
                    }while(!idx_error && bal_value >= 1.0e-4 && i_loop_inner < N_max2);
                }
                // std::cout << "Local eq point found, Continue ? " << std::endl;
                // cin.get();
                i_loop_outer++;
            
            }while(idx_error && bal_value >= 1.0e-4 && i_loop_outer < N_max1  );
            
            
            
            
            for (i = 0; i < numAtoms; i++)
            {
                delete [] bal_deriv [i];
                bal_deriv[i] = NULL;

                delete [] position_m [i];
                position_m [i] = NULL;
            }
            delete [] bal_deriv;
            bal_deriv  = NULL;
            delete [] position_m;
            position_m = NULL;
            
        }       
            
    }
    
    
    void FindGlobMin::SetBalanceValueAndDeriv_EXP(int& idx_error, 
                                                  LIBMOL::REAL coeffi_ad, 
                                                  LIBMOL::REAL& bal_value, 
                                                  LIBMOL::REAL** bal_deriv, 
                                                  LIBMOL::REAL** pos_m)
    {
        // std::cout << "New function in test " << std::endl;
        
        // int numOptims= (int)allOptimSets[0].atomGroups.size();
        int numAtoms = (int)allAtoms.size();
        int dim      = (int)allAtoms[0].coords.size();
        
        int i, j;

        idx_error = 0;
        
        FF::GetObjValue toolObjVandD;

        if(workSpace == 1)      // test Cartesian space first
        {
            
            toolObjVandD.workSpace = 1;
            LIBMOL::REAL tObjValue = toolObjVandD.getAll(allAtoms, allBonds, allAngles, 
                                                         allTorsions, allPlanes, allChirals);
   
            bal_value = tObjValue - curObjValue;
      
            // std::cout << "bal_value " << bal_value << std::endl;
            
            
            if (bal_value <= 0)
            { 
                return;
            }
            else 
            {
                int tOptSize = (int)allOptimSets.size();
                
                if (tOptSize > 0)
                {
                    LIBMOL::REAL bal_value_sqrt;
                    bal_value_sqrt = sqrt(bal_value); // test here

                    LIBMOL::REAL r_diff       = 0.0;
                    LIBMOL::REAL r_diff_sq    = 0.0;
                    LIBMOL::REAL r_diff_sqrt  = 0.0;

                    LIBMOL::REAL ** r_diff_comp;
                    r_diff_comp = new LIBMOL::REAL * [numAtoms];
                    for (i = 0; i < numAtoms; i++)
                    {
                        r_diff_comp[i] = new LIBMOL::REAL [dim];
                        for (j=0; j < dim; j++)
                        {
                            r_diff_comp[i][j] = 0.0;
                        }
                    }
                
                    for (i = 0; i < numAtoms; i++)
                    {
                        for (j = 0; j < dim; j++)
                        {                   
                            r_diff_comp[i][j]  = allAtoms[i].coords[j] - 
		            allOptimSets[0].atoms[i].coords[j];
                            r_diff_sq += pow(r_diff_comp[i][j], 2.0); 
                        }
                    }
                
                    r_diff_sqrt = sqrt(r_diff_sq);

                    if(r_diff_sqrt > 1.0e-8)
                    {
                    
                        // a. calculate the value of the balance function
                        r_diff   = (1+r_diff_sqrt)/r_diff_sqrt;
                        bal_value = bal_value_sqrt*r_diff;

                        // b. calculate the derivs of the balance function
                        // SetInForces();
                        LIBMOL::REAL bal_inv = 1.0/bal_value_sqrt;
                        LIBMOL::REAL AA, BB;
                        AA = 0.50*bal_inv*(r_diff);
                        BB = -bal_value_sqrt/(r_diff_sq*r_diff_sqrt);
                        for(i = 0; i < numAtoms; i++)
                        {
                            for (j = 0; j < dim; j++)
                            {
                                bal_deriv[i][j]  
                                    =  AA*allAtoms[i].forces[j]
			            + BB*r_diff_comp[i][j];
		                //   cout << "bal_deriv["<<i<<"]["<<j<<"]= " 
                                //        << bal_deriv[i][j] 
                                //        << endl; 
                            }
                        }

                        for (i = 0; i < numAtoms; i++)
                        {
                            delete [] r_diff_comp[i];
                            r_diff_comp[i] = NULL;
                        }              
                        delete [] r_diff_comp;
                        r_diff_comp = NULL; 
                        return; 
                    }
	            else 
                    {
                        for (i = 0; i < numAtoms; i++)
                        {
                            delete [] r_diff_comp[i];
                            r_diff_comp[i] = NULL;
                        }              
                        delete [] r_diff_comp;
                        r_diff_comp = NULL;
  
                        idx_error = 1;
                        return;
                    }
                }
            }
        }
    }
    
    // Overall controller of the global optimization section
    
    void FindGlobMin::Driver()
    {
        time_t tStart, tEnd;
        std::time (&tStart);
        std::cout << "Geometrical optimization started at " << std::ctime(&tStart);
        
        // hasIniCoords = false;
        
        PreIdealization(); 
        
        if (lGlob)
        {
            if (gWorkMode == 1)
            {
                TunnellingAndMinDriver();
            }
        
            // ProIdealization();
            std::cout << "The following are optimized sets " << std::endl;
            // std::sort(allOptimSets.begin(), allOptimSets.end(), sortOptSetByValues);
            for (std::vector<OptimSet>::iterator iO=allOptimSets.begin(); 
                    iO != allOptimSets.end(); iO++)
            {
                std::cout << "Set objective value " << iO->objValue << std::endl;
            }
        
            // SelectBestOpt();
        
            std::time(&tEnd);
            std::cout << " Geometrical optimization finished at " << std::ctime(&tEnd);
      
            std::cout  << "it takes " << std::setprecision(3) << std::difftime(tEnd,tStart) 
                       << " seconds" << std::endl;
        }
    }
    
    
    void FindGlobMin::Driver(bool useCoords=false, int nOpt=0)
    {
        if (useCoords)
        {
           UpdatePositionsWithRandomShifts(); 
        }
        else
        {
            ProIdealization();
        }
        
        if (nOpt==0)
        {
           numOptimSets = (int)allAtoms.size(); 
        }
        else
        {
            numOptimSets = nOpt;
        }
        
        if (gWorkMode == 1)
        {
            TunnellingAndMinDriver(nOpt);
        }
        
        ProIdealization();
        
        SelectBestOpt();
       
    }
        
}
