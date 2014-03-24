
/* 
 * File:   AllSystem.cpp
 * Author: flong
 *
 * Created on July 5, 2012, 10:28 AM
 */

#include "AllSystem.h"

namespace LIBMOL
{
    AllSystem::AllSystem():hasCoords(false),
                           hasCCP4Type(false),
                           itsContainMetal(false),
                           itsCurAngleSeriNum(ZeroInt),
                           itsCurAngle(NullPoint),
                           itsCurTorsionSeriNum(ZeroInt),
                           itsCurTorsion(NullPoint),
                           itsCurChiralSeriNum(ZeroInt),
                           itsCurChiral(NullPoint)
    {      
    }
 
    AllSystem::AllSystem(const AllSystem& tAllSys):hasCoords(tAllSys.hasCoords),
                                                   hasCCP4Type(tAllSys.hasCCP4Type),
                                                   itsContainMetal(tAllSys.itsContainMetal),
                                                   itsCurAngleSeriNum(ZeroInt),
                                                   itsCurAngle(NullPoint),
                                                   itsCurTorsionSeriNum(ZeroInt),
                                                   itsCurTorsion(NullPoint),
                                                   itsCurChiralSeriNum(ZeroInt),
                                                   itsCurChiral(NullPoint)
    {      
        AddAtoms(tAllSys.allAtoms);
        AddBonds(tAllSys.allBonds);
        AddAngles(tAllSys.allAngles);
        AddTorsions(tAllSys.allTorsions);
        AddChirals(tAllSys.allChirals);
        AddPlanes(tAllSys.allPlanes);
        AddRings(tAllSys.allRings);
        setAtomsMetalType();
    }
    
    
    AllSystem::AllSystem(const std::vector<AtomDict>& tAllAtoms,
                         const std::vector<int> &  tAllHAtomIdx,
                         const std::vector<BondDict>& tAllBonds, 
                         const std::vector<AngleDict>& tAllAngles, 
                         const std::vector<TorsionDict>& tAllTorsions, 
                         const std::vector<ChiralDict>& tAllChirals, 
                         const std::vector<PlaneDict>& tAllPlanes, 
                         const std::map<ID, std::vector<RingDict> >& tAllRings)
                          :hasCoords(false),
                           hasCCP4Type(false),
                           itsContainMetal(false),
                           itsCurAngleSeriNum(ZeroInt),
                           itsCurAngle(NullPoint),
                           itsCurTorsionSeriNum(ZeroInt),
                           itsCurTorsion(NullPoint),
                           itsCurChiralSeriNum(ZeroInt),
                           itsCurChiral(NullPoint)
    {
        
        AddAtoms(tAllAtoms);
        for (std::vector<int>::const_iterator iH = tAllHAtomIdx.begin();
                iH != tAllHAtomIdx.end(); iH++)
        {
            allHAtomIdx.push_back(*iH);
        }
        
        AddBonds(tAllBonds);
        AddAngles(tAllAngles);
        AddTorsions(tAllTorsions);
        AddChirals(tAllChirals);
        AddPlanes(tAllPlanes);
        AddRings(tAllRings);
        setAtomsMetalType();
        
    }
    
    
    AllSystem::AllSystem(DictCifFile& tCifObj):hasCoords(tCifObj.hasCoords),
                                               hasCCP4Type(tCifObj.hasCCP4Type),
                                               itsContainMetal(false),
                                               itsCurAngleSeriNum(ZeroInt),
                                               itsCurAngle(NullPoint),
                                               itsCurTorsionSeriNum(ZeroInt),
                                               itsCurTorsion(NullPoint),
                                               itsCurChiralSeriNum(ZeroInt),
                                               itsCurChiral(NullPoint)
    {
        AddAtoms(tCifObj.allAtoms);
        AddBonds(tCifObj.allBonds);
        // AddTorsions(tCifObj.allTorsions);
        AddChirals(tCifObj.allChirals);
        setSysProps();
        
    }
    
    AllSystem::AllSystem(Molecule& tMol):hasCoords(tMol.hasCoords),
                                               itsContainMetal(false),
                                               itsCurAngleSeriNum(ZeroInt),
                                               itsCurAngle(NullPoint),
                                               itsCurTorsionSeriNum(ZeroInt),
                                               itsCurTorsion(NullPoint),
                                               itsCurChiralSeriNum(ZeroInt),
                                               itsCurChiral(NullPoint)
    {
        AddAtoms(tMol.atoms);
        AddBonds(tMol.bonds);
        AddChirals(tMol.chirals);
        setSysProps();
    }
    
    AllSystem::AllSystem(const CodClassify& tProCodSys)
                          :itsCurAngleSeriNum(ZeroInt),
                           itsCurAngle(NullPoint),
                           itsCurTorsionSeriNum(ZeroInt),
                           itsCurTorsion(NullPoint),
                           itsCurChiralSeriNum(ZeroInt),
                           itsCurChiral(NullPoint)
    {
        AddAtoms(tProCodSys.allAtoms);
        AddBonds(tProCodSys.allBonds);
        AddAngles(tProCodSys.allAngles);
        AddTorsions(tProCodSys.allTorsions);
        AddChirals(tProCodSys.allChirals);
        AddPlanes(tProCodSys.allPlanes);
        AddRings(tProCodSys.allRings);
    }
    
    AllSystem::~AllSystem()
    {        
        if(itsCurAngle)
        {
            delete itsCurAngle;
            itsCurAngle = NULL;
        }
        if(itsCurTorsion)
        {
            delete itsCurTorsion;
            itsCurTorsion = NULL;
        }
        if(itsCurChiral)
        {
            delete itsCurChiral;
            itsCurChiral = NULL;
        }
    }
   
    void AllSystem::resetSystem(CodClassify& tCodSys)
    {
        
        allAtoms.clear();
        allBonds.clear();
        allAngles.clear();
        allTorsions.clear();
        allChirals.clear();
        allPlanes.clear();
        allRings.clear();
        
        AddAtoms(tCodSys.allAtoms);
        AddBonds(tCodSys.allBonds);
        AddAngles(tCodSys.allAngles);
        AddTorsions(tCodSys.allTorsions);
        AddChirals(tCodSys.allChirals);
        AddPlanes(tCodSys.allPlanes);
        AddRings(tCodSys.allRings);
    }
    
    void AllSystem::AddAtom(const AtomDict & tAtom)
    {
        allAtoms.push_back(tAtom);
    }
    
    void AllSystem::AddAtoms(const std::vector<AtomDict>& tAllAtoms)
    {
        for (std::vector<AtomDict>::const_iterator iA=tAllAtoms.begin();
                iA !=tAllAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
    }
    
    void AllSystem::AddBonds(const std::vector<BondDict>& tAllBonds)
    {
        for (std::vector<BondDict>::const_iterator iB=tAllBonds.begin();
                iB != tAllBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
    }
    
    void AllSystem::AddBond(const BondDict& tBond)
    {
        allBonds.push_back(tBond);
    }
    
    void AllSystem::AddAngles(const std::vector<AngleDict>& tAllAngles)
    {
        for (std::vector<AngleDict>::const_iterator iA=tAllAngles.begin();
                iA != tAllAngles.end(); iA++)
        {
            allAngles.push_back(*iA);
        }
    }
    
    void AllSystem::AddAngle(const AngleDict& tAngle)
    {
        allAngles.push_back(tAngle);
    }
    
    
    void AllSystem::AddTorsions(const std::vector<TorsionDict>& tAllTorsions)
    {
        for (std::vector<TorsionDict>::const_iterator iT=tAllTorsions.begin();
                iT != tAllTorsions.end(); iT++)
        {
            allTorsions.push_back(*iT);        
        }
    }
    
    void AllSystem::AddTorsion(const TorsionDict& tTorsion)
    {
        allTorsions.push_back(tTorsion);
    }
    
    void AllSystem::AddChirals(const std::vector<ChiralDict>& tAllChirals)
    {
        for (std::vector<ChiralDict>::const_iterator iC=tAllChirals.begin();
                iC != tAllChirals.end(); iC++)
        {
            allChirals.push_back(*iC);
        }
    }
    
    void AllSystem::AddChiral(const ChiralDict& tChiral)
    {
        allChirals.push_back(tChiral);
    }
    
    void AllSystem::AddPlanes(const std::vector<PlaneDict>& tAllPlanes)
    {
        for(std::vector<PlaneDict>::const_iterator iP=tAllPlanes.begin();
                iP!=tAllPlanes.end(); iP++)
        {
            allPlanes.push_back(*iP);
        }
    }
    
    void AllSystem::AddPlane(const PlaneDict& tPlane)
    {
        allPlanes.push_back(tPlane);
    }
    
    void AllSystem::AddRings(const std::map<ID, std::vector<RingDict> > & tAllRings)
    {
        
        for (std::map<ID, std::vector<RingDict> >::const_iterator iM=tAllRings.begin();
                iM!=tAllRings.end(); iM++)
        {
            for (std::vector<RingDict>::const_iterator iR=iM->second.begin();
                    iR != iM->second.end(); iR++)
            {
                allRings[iM->first].push_back(*iR);
            }
        }
    }
    
    void AllSystem::AddRing(const RingDict& tRing)
    {
        allRingsV.push_back(tRing);
    }
    
    int AllSystem::atomPosition(ID tID)
    {
        for (int i=0; i<(int)allAtoms.size(); i++)
        {
            if(allAtoms[i].id.compare(tID)==0)
            {
                return i;
            }
        }
        
        return -1;
        
    }
    
    
    void AllSystem::setSysProps()
    {
        setHydroAtomConnect();
            
        setAtomsBondingAndChiralCenter();
      
        setAtomsCChemType();
   
            
        setAtomsMetalType();
        
        setAtomsPartialCharges();
          
        setAllAngles();
           
        ringDetecting();
        
        for (std::map<ID, std::vector<RingDict> >::iterator iMR=allRings.begin();
                    iMR != allRings.end(); iMR++)
        {
            for (std::vector<RingDict>::iterator iR=iMR->second.begin();
                        iR != iMR->second.end(); iR++)
            {
                
                iR->setAtmsLink(allAtoms);
                allRingsV.push_back(*iR);
            }
        }
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
               iA !=allAtoms.end(); iA++)
        {
            iA->inRings.clear();
        }
        
        for (int i=0; i < (int)allRingsV.size(); i++)
        {
            for (std::vector<AtomDict>::iterator iA=allRingsV[i].atoms.begin();
                    iA != allRingsV[i].atoms.end(); iA++)
            {
                int tIdx = atomPosition(iA->id);
                if (tIdx >=0 && tIdx < (int)allAtoms.size())
                {
                    allAtoms[tIdx].inRings.push_back(i);
                }
                else
                {
                    std::cout << "The atom " << iA->id 
                              << " in ring " << allRingsV[i].rep 
                              << " does not exists, bug !"
                              << std::endl;
                }
            }
        }
        /*
        std::cout << "There are " << allRingsV.size() << " rings " << std::endl
                  << "These are: " << std::endl;
        for (std::map<ID, std::vector<RingDict> >::iterator iR=allRings.begin();
                iR !=allRings.end(); iR++)
        {
            std::cout << "Ring " << iR->first << "." << std::endl;
            for (std::vector<RingDict>::iterator iRV=iR->second.begin();
                    iRV !=iR->second.end(); iRV++)
            {
                std::cout << "Ring size:  " << iRV->atoms.size() << std::endl;;
                for (std::vector<AtomDict>::iterator iA=iRV->atoms.begin(); 
                        iA != iRV->atoms.end(); iA++)
                {
                    std::cout<< "Atom " << iA->id << std::endl;
                }
            }
        }
        */
        // Double check 
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if ((int)iA->inRings.size() !=0)
            {
                std::cout << "Atom " << iA->id << " is in " 
                          << (int)iA->inRings.size() 
                          << " ring(s) " << std::endl;
                for (std::vector<int>::iterator iR=iA->inRings.begin();
                        iR != iA->inRings.end(); iR++)
                {
                    std::cout << "Ring " << allRingsV[*iR].rep << std::endl;
                }
            }
        }
         */
              
        detectPlaneGroups(); 
       
        setAllTorsions2();
        
        if (!hasCCP4Type)
        {
            setAtomsCCP4Type();
        }
        
        for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
        {
            std::cout << "\nAtom " << iA->seriNum << " : " << std::endl
                      << "Its ID : " << iA->id << std::endl
                      << "Its Chemical Type : " << iA->chemType << std::endl
                      << "Its bonding index : "   << iA->bondingIdx << std::endl
                      // << "Its COD chemical type : " << iA->cChemType << std::endl
                      << "Its CCP4 atom type : "  << iA->ccp4Type << std::endl 
                      << "Its residue Name: " << iA->resName<< std::endl;
            std::cout << "Its connected atoms are : " << std::endl;
            for (std::vector<int>::iterator iSer= iA->connAtoms.begin();
                        iSer != iA->connAtoms.end(); iSer++)
            {
                std::cout << allAtoms[*iSer].id << std::endl;
            }
        
        }
    }
    
    int AllSystem::getNumSpecAtomConnect(int idxAtm, ID tChemType)
    {
        int nA=0;
        for (std::vector<int>::iterator iC=allAtoms[idxAtm].connAtoms.begin();
                iC !=allAtoms[idxAtm].connAtoms.end(); iC++)
        {
            if (allAtoms[*iC].chemType.compare(tChemType)==0)
            {
                nA++;
            }
        }
        return nA;
    }
    
    int AllSystem::getNumSpecAtomConnect(std::vector<AtomDict>::iterator iA,
                                    ID tChemType)
    {
        int nO=0;
        for (std::vector<int>::iterator iC=iA->connAtoms.begin();
                iC !=iA->connAtoms.end(); iC++)
        {
            if (allAtoms[*iC].chemType.compare(tChemType)==0)
            {
                nO++;
            }
        }
        return nO;
    }
    
    int AllSystem::getNumOxyConnect(std::vector<AtomDict>::iterator iA)
    {
        int nO=0;
        for (std::vector<int>::iterator iC=iA->connAtoms.begin();
                iC !=iA->connAtoms.end(); iC++)
        {
            if (allAtoms[*iC].chemType.compare("O")==0)
            {
                nO++;
            }
        }
        return nO;
    }
    
    void AllSystem::setHydroAtomConnect()
    {
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA!=allAtoms.end(); iA++)
        {
            if(iA->chemType.compare("H")==0)
            {
                allHAtomIdx.push_back(iA->seriNum);
            }
            else
            {
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                    iNB !=iA->connAtoms.end(); iNB++)
                {
                    if (allAtoms[*iNB].chemType.compare("H")==0)
                    {
                        iA->connHAtoms.push_back(*iNB);
                    }
                }
            }
        }
        
        // Check
        /*
        for (std::vector<AtomDict>::iterator iHA=allHAtomIdx.begin();
                iHA !=allHAtomIdx.end(); iHA++)
        {
            std::cout << "Atom " << iHA->seriNum << "its element ID " 
                    << iHA->chemType << std::endl;
        }
        
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA!=allAtoms.end(); iA++)
        {
            std::cout << "Atom " << iA->seriNum << " connects to " 
                    << (int)iA->connHAtoms.size() << " H atoms " << std::endl;
            if ((int)iA->connHAtoms.size() !=0)
            {
                std::cout << "These atoms are : " << std::endl;
                for (std::vector<int>::iterator iH=iA->connHAtoms.begin();
                        iH !=iA->connHAtoms.end(); iH++)
                {
                    std::cout << "Atom " << allAtoms[*iH].id << std::endl;
                }
            }
        }
        */
        
    }
    
    
    void AllSystem::setAtomsBondingAndChiralCenter()
    {
        
        // First round
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("C")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                }
                else if(t_len==2)
                {
                    iAt->chiralIdx  = 0;
                    if (getNumOxyConnect(iAt)==1)
                    {
                        // water is removed 
                        iAt->bondingIdx=2;
                    }
                    else
                    {
                        iAt->bondingIdx=1;
                    }
                }
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    } 
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp
                //{
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                //}
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                } 
                else if (t_len==1)
                {
                    // triple bond 
                    iAt->chiralIdx = 0;
                    iAt->bondingIdx= 1;
                }
            }
            else if (iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                    iAt->formalCharge = -1;
                }
            }
            else if (iAt->chemType.compare("O")==0)
            {
                if ((int)iAt->connAtoms.size()==2)
                {
                    iAt->bondingIdx = 2;
                }
            }
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    
                    if (iAt->chiralIdx ==0)
                    {
                        std::vector<ID> atps;
                        for (std::vector<int>::iterator iNA=iAt->connAtoms.begin();
                                iNA !=iAt->connAtoms.end(); iNA++)
                        {
                            if (std::find(atps.begin(), atps.end(), allAtoms[*iNA].chemType)==atps.end())
                            {
                                atps.push_back(allAtoms[*iNA].chemType);
                            }
                        }
                        if ((int)atps.size() >2)
                        {
                            iAt->chiralIdx  = 2;
                        }
                        else
                        {
                            iAt->chiralIdx =0;
                        }
                        // iAt->chiralIdx  = 2;
                    }
                   
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==3)
                {
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(allAtoms[*iCA].bondingIdx == 2)
                        {
                            l_sp2 = true;
                        }
                    }
                    if (l_sp2)
                    {
                        // Now we can say this atom is in sp2 orbits 
                        iAt->chiralIdx  =  0;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        if (iAt->chiralIdx ==0)
                        {
                            iAt->chiralIdx  = 2;
                        }
           
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
        
        // Further check if a chiral center is a real one
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
        }
        
        /*
                // First round
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }
            if (iAt->chemType.compare("C")==0)
            {
                //int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp 
                //{   // should do on the next round when all NB atoms are set
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                // }
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx =  2;
                } 
            }
            else if (iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
            }
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {            
            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                int t_len =0;
                for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                     iConn !=iAt->connAtoms.end(); iConn++)
                {
                    if(!allAtoms[*iConn].isMetal)
                    {
                        t_len++;
                    }
                }
                if(t_len==3)
                {
                    
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(allAtoms[*iCA].bondingIdx == 2)
                        {
                            l_sp2 = true;
                        }
                    }
                    if (l_sp2)
                    {
                        // Now we can say this atom is in sp2 orbits 
                        iAt->chiralIdx  = -1;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        iAt->chiralIdx  =  1;
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
 
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx ==1)
            {
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
        }
        */
        // No need for the third round, those could be defined in 
        // the first round
        // Check
        /*
        std::cout << "Chiral and plane feather for atoms in the system" 
                  << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            if (iAt->chiralIdx == -1)
            {
                std::cout << "Atom " << iAt->id << " may be in planes " 
                        << std::endl;
            }
            else if (iAt->chiralIdx == 1)
            {
                std::cout << "Atom " << iAt->id << " may be in a chiral center "
                        << std::endl;
            }
            else if (iAt->chiralIdx==2)
            {
                std::cout << "Atom " << iAt->id 
                        << " may be in a chiral center but the volume sign undefined "
                        << std::endl;
            }
            else 
            {
                std::cout << "Atom " << iAt->id << " is not a chiral center "
                        << std::endl;
            }
        } 
         */  
        
    }
    
    void AllSystem::setAtomsCChemType()
    {
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            // transfer the characters from 2nd position to low case
            if((int)iA->chemType.size()!=1)
            {
                ID tS = iA->chemType.substr(1);
                StrLower(tS);
                iA->chemType = iA->chemType[0] + tS;
            }
            
            // Embed planar info into the chemType to match COD definition 
            iA->cChemType = iA->chemType;
            if (iA->bondingIdx ==2)
            {
                 StrLower(iA->cChemType);
            }
            
            //std::cout << "Atom ID: " << iA->id << std::endl
            //        << "Atom chemType: " << iA->chemType << std::endl
            //        << "Atom chemType of COD classes: " << iA->cChemType << std::endl;
        }
    }
    
    void AllSystem::setAtomsCCP4Type()
    {
        // should be done after rings and chiral centers are detected
        
        CCP4AtomType  aCPP4TypeTool(allAtoms, allRings);
        aCPP4TypeTool.setAllAtomsCCP4Type();
        for (int i=0; i < (int)aCPP4TypeTool.allAtoms.size(); i++)
        {
            allAtoms[i].ccp4Type = aCPP4TypeTool.allAtoms[i].ccp4Type;
            //std::cout << "Atom " << allAtoms[i].id 
            //          << " CCP4 atom type is " << allAtoms[i].ccp4Type 
            //          << std::endl;
        }
    }
    void AllSystem::setDefaultCoordGeos()
    {
        DefaultCoordGeos[2]  = "LINEAR";
        DefaultCoordGeos[3]  = "TRIGONAL-PLANAR";
        DefaultCoordGeos[4]  = "TETRAHEDRAL";
        DefaultCoordGeos[5]  = "TRIGONAL-BIPYRAMID";
        DefaultCoordGeos[6]  = "OCTAHEDRAL";
        DefaultCoordGeos[7]  = "CAPPED-OCTAHEDRAL";
        DefaultCoordGeos[8]  = "CUBIC";
        DefaultCoordGeos[9]  = "TRICAPPED-TRIGONAL-PRISMATIC";
        DefaultCoordGeos[10] = "BICAPPED-SQUARE-ANTIPRISMATIC";
        DefaultCoordGeos[11] = "ALL-FACE-CAPPED-TRIGONAL-PRISMATIC";
        DefaultCoordGeos[12] = "CUBOCTAHEDRON";
        
        std::string clibMonDir(std::getenv("CLIBD_MON"));
        std::string metDefCoordGeoFileName = clibMonDir + "/allMetalDefCoordGeos.table";
        std::ifstream metDefCoordGeoFile(metDefCoordGeoFileName.c_str());
        
        if(metDefCoordGeoFile.is_open())
        {
            std::string tRecord="";
            
            while(!metDefCoordGeoFile.eof())
            {
                std::getline(metDefCoordGeoFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                
                if ((int)tBuf.size() ==3)
                {
                    int cn = StrToInt(TrimSpaces(tBuf[1]));
                    
                    DefaultCoordGeos2[TrimSpaces(tBuf[0])][cn]=TrimSpaces(tBuf[2]);
                }
            }
            metDefCoordGeoFile.close();
        }
        
    }   
    
    const bool AllSystem::containMetal()
    {
        return itsContainMetal;
    }
    
    void AllSystem::setAtomsMetalType()
    {
        setDefaultCoordGeos();
        
        ID metals[] = {"Li", "li", "Na", "na", "K",  "k",  "Rb", "rb", "Cs", "cs", "Fr", "fr",
                     "Be", "be", "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba", "Ra", "ra",
                     "Sc", "sc", "Y",  "y",
                     "B",  "b",  "Si", "si", "Ge", "ge", "As", "as", "Sb", "sb", "Te", "te", "Po", "po",
                     "Ti", "ti", "Zr", "zr", "Hf", "hf", "Rf", "rf",
                     "V",  "v"   "Nb", "nb", "Ta", "ta", "Db", "db", 
                     "Cr", "cr", "Mo", "mo", "W",  "w",  "Sg", "sg", 
                     "Mn", "mn", "Tc", "tc", "Re", "re", "Bh", "bh",  
                     "Fe", "fe", "Ru", "ru", "Os", "os", "Hs", "hs",   
                     "Co", "co", "Rh", "rh", "Ir", "ir", "Mt", "mt",  
                     "Ni", "ni", "Pd", "pd", "Pt", "pt", "Ds", "ds",  
                     "Cu", "cu", "Ag", "ag", "Au", "au", "Rg", "rg",   
                     "Zn", "zn", "Cd", "cd", "Hg", "hg",   
                     "Al", "al", "Ga", "ga", "In", "in", "Ti", "ti", 
                     "Sn", "sn", "Pb", "pb", "Bi", "bi"};
        
        MetalTable.assign(metals, metals+122);
        
                
        //std::cout << "Metal Elements :" << std::endl;
        //for (std::vector<ID>::iterator iM =MetalTable.begin();
        //       iM !=MetalTable.end(); iM++)
        //{
        //    std::cout << *iM << std::endl;
        //}
        //std::cout << "Metal Elements :" << std::endl;
        //for (std::vector<ID>::iterator iM =MetalTable.begin();
        //       iM !=MetalTable.end(); iM++)
        //{
        //    std::cout << *iM << std::endl;
        //}
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        { 
            std::vector<ID>::iterator iFind = std::find(MetalTable.begin(), 
                                                        MetalTable.end(), iA->chemType);
            if (iFind !=MetalTable.end())
            {
                iA->isMetal = true;
                std::cout <<iA->id << " is a metal atom." << std::endl;
                itsContainMetal = true;
                
                int cn = (int)iA->connAtoms.size();
                ID  ct = iA->chemType;
       
                if (iA->metalGeo.empty())
                {
                    std::map<ID, std::map<int,ID> >::iterator iFind1 
                                        =DefaultCoordGeos2.find(ct);
                    if (iFind1 !=DefaultCoordGeos2.end())
                    {
                        std::map<int,ID>::iterator iFind2
                        = DefaultCoordGeos2[ct].find(cn);
                        if (iFind2 !=DefaultCoordGeos2[ct].end())
                        {
                            iA->metalGeo=DefaultCoordGeos2[ct][cn];
                        }
                        else
                        {
                            iA->metalGeo = DefaultCoordGeos[cn];
                        }
                    }
                    else
                    {
                        iA->metalGeo = DefaultCoordGeos[cn];
                    }
                }
                
                std::cout << "Its coordination number is set to " << cn << std::endl
                          << "its default coordination geometry is set to " 
                          << iA->metalGeo << std::endl;
            }
        }        
        
    }
    
    void AllSystem::setAtomsPartialCharges()
    {
        // Need to use bond-order or valence in p-table later on 
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if (iA->parCharge == -1.0)
            {
                // make the charges dislocated 
                if (iA->chemType =="O" && (int)iA->connAtoms.size()==1)
                {
                    int tRoot = iA->connAtoms[0];
                    std::vector<int> iOs;
                    for (std::vector<int>::iterator iNA=allAtoms[tRoot].connAtoms.begin();
                            iNA !=allAtoms[tRoot].connAtoms.end(); iNA++)
                    {
                        if (allAtoms[*iNA].id !=iA->id && allAtoms[*iNA].chemType=="O"
                             && (int)allAtoms[*iNA].connAtoms.size()==1)
                        {
                            iOs.push_back(*iNA);
                        }
                    }
                    
                    if ((int)iOs.size()==1)
                    {
                        iA->parCharge              = 0.5;
                        allAtoms[iOs[0]].parCharge = 0.5;
                    }
                }
                
            }
        }
    }
    
    // setAllAngles() may not needed in the future
    void AllSystem::setAllAngles()
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
                    std::vector<int> tVec;
                    aAng.anchorID  = allAtoms[i].id;
                    aAng.anchorPos = i;
                    aAng.atoms.push_back(i);
                       
                    if ((int) allAtoms[i1].connAtoms.size() >=
                        (int) allAtoms[i1].connAtoms.size())
                    {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                        tVec.push_back(i1);
                        tVec.push_back(i2);
                    }
                    else
                    {
                        aAng.atoms.push_back(i2);
                        aAng.atoms.push_back(i1);
                        tVec.push_back(i1);
                        tVec.push_back(i2);
                    }
                        
                    aAng.value        = 0.0;
                    aAng.sigValue     = 3.0;
                    aAng.numCodValues = 0;
                    allAngles.push_back(aAng);
                    allAnglesIdxs[i].push_back(tVec);
                }
            }
        }       

        // Check 
        std::cout << "There are " << (int)allAngles.size() 
                  << " angles in the system " << std::endl;
        
        std::cout << "These angles are : " << std::endl;
        
        for (std::map<int, std::vector<std::vector<int> > >::iterator tAG 
               =allAnglesIdxs.begin(); tAG != allAnglesIdxs.end(); tAG++)
        {
            std::cout << "There are " << (int)tAG->second.size() 
            << " Angles centered on atom " << allAtoms[tAG->first].id << std::endl
                    << "They are " << std::endl;
                    
            for(std::vector<std::vector<int> >::iterator iAN =tAG->second.begin();
                    iAN != tAG->second.end(); iAN++)
            {
                std::cout << "atoms: " << allAtoms[tAG->first].id <<",  ";
                for (std::vector<int>::iterator iAt = iAN->begin(); 
                        iAt !=iAN->end(); iAt++)
                {
                  std::cout << allAtoms[*iAt].id << ", ";
                }
                std::cout<<std::endl;
            }
        }
    }
   
    void AllSystem::checkAngConstraints()
    {
        std::map<int, std::vector<int> >  tSP2Angs, tSP3Angs;
        
        // Get all sp2 and sp3 angles centered at different atoms
        
        for (unsigned i=0;  i < allAngles.size(); i++)
        {
            int iCenAtom= allAngles[i].atoms[0];
            if (allAtoms[iCenAtom].bondingIdx==2)
            {
                tSP2Angs[iCenAtom].push_back(i);
            }
            else if (allAtoms[iCenAtom].bondingIdx==3)
            {
                tSP3Angs[iCenAtom].push_back(i);
            }
        }
        
        // check all sp2 atoms
        for (std::map<int, std::vector<int> >::iterator iSetAngs=tSP2Angs.begin();
                iSetAngs !=tSP2Angs.end(); iSetAngs++)
        {
            checkSP2Constraints(iSetAngs->second);
        }
        
        // check all sp3 atoms
        for (std::map<int, std::vector<int> >::iterator iSetAngs=tSP3Angs.begin();
                iSetAngs !=tSP3Angs.end(); iSetAngs++)
        {
            checkSP3Constraints(iSetAngs->second);
        }
    }
    
    void AllSystem::checkSP2Constraints(std::vector<int> tAngIdxs)
    {
        REAL sp2Sum  =360.00;
        REAL sp2Diff =0.0;
        REAL tAngSum = 0.0;
        if (tAngIdxs.size()==3)
        {
            for (std::vector<int>::iterator iA=tAngIdxs.begin();
                    iA !=tAngIdxs.end(); iA++)
            {
                tAngSum+=(allAngles[*iA].value);
            }
            
            sp2Diff = sp2Sum - tAngSum;
            
            if (fabs(sp2Diff) >0.01)
            {
                sp2Diff= sp2Diff/3.0;
            
                // add the diff to individual angles
                tAngSum=0.0;
                for (std::vector<int>::iterator iA=tAngIdxs.begin();
                        iA !=tAngIdxs.end(); iA++)
                {
                    allAngles[*iA].value +=sp2Diff;
                    tAngSum+=(allAngles[*iA].value);
                }
                
                // Further check.  
                sp2Diff = sp2Sum - tAngSum;
                //  Now sp2Diff should be negligible small, add it to an angle
                allAngles[tAngIdxs[0]].value += sp2Diff;
            }
            else 
            {
                // very small, just add it to an angle 
                allAngles[tAngIdxs[0]].value += sp2Diff;
            }
            
            //Debug, output  the modified sum, should be 360.00
            tAngSum=0.0;
            for (std::vector<int>::iterator iA=tAngIdxs.begin();
                        iA !=tAngIdxs.end(); iA++)
            {
                tAngSum+=(allAngles[*iA].value);
            }
            
            std::cout << "Sum of values for angles centered at atom " 
                      << allAtoms[allAngles[tAngIdxs[0]].atoms[0]].id 
                      << " is " << tAngSum << std::endl;
        }
    }
    
    void AllSystem::checkSP3Constraints(std::vector<int> tAngIdxs)
    {
        
    }
    
    /*
    void AllSystem::setOneTorsion(std::vector<int> tAV, REAL tValue,
                                    int tPeriod)
    {
        
        TorsionDict tTor;
      
        for (int i=0; i <(int)tAV.size(); i++) 
        {
            tTor.atoms.push_back(tAV[i]);
        }
        
        tTor.value  = tValue;
        tTor.period = tPeriod;
        
        allTorsions.push_back(tTor);
        
    }
    
    */
    
    void AllSystem::setOneTorsion(std::vector<int> tAV, REAL tValue,
                                    int tPeriod)
    {
        
        TorsionDict tTor;
      
        tTor.seriNum = itsCurTorsionSeriNum;
        for (int i=0; i <(int)tAV.size(); i++) 
        {
            tTor.atoms.push_back(tAV[i]);
        }
        
        tTor.value  = tValue;
        tTor.valueST = tValue;
        tTor.period = tPeriod;
        
        allTorsions.push_back(tTor);
        itsCurTorsionSeriNum++;
        
    }
    
    void AllSystem::SetOneSP1SP1Bond(int tIdx1, int tIdx2)
    {
        std::vector<int> aTS;
        int i1=-1, i2=-1;
        for(std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                iA1 !=allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if(*iA1 !=tIdx2)
            {
                i1= *iA1;
                break;
            }
        }
            
        if (i1 !=-1)
        {
            aTS.push_back(i1);
        }
        else
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx1].id << " has one connection"
                      << std::endl;
            return;
        }
             
        aTS.push_back(tIdx1);
        aTS.push_back(tIdx2);
            
        for(std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                iA2 !=allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if(*iA2 !=tIdx1)
            {
                i2= *iA2;
                break;
            }
        }    
          
        if (i2 !=-1)
        {
            aTS.push_back(i2);
        }
        else
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx2].id << " has one connection"
                      << std::endl;
            return;
        }
        
        setOneTorsion(aTS, 180.0, 1);
         
        std::cout << "Torsion on  " << allAtoms[allTorsions[itsCurTorsionSeriNum-1].atoms[0]].id
                << " --> " << allAtoms[allTorsions[itsCurTorsionSeriNum-1].atoms[1]].id
                << " --> " << allAtoms[allTorsions[itsCurTorsionSeriNum-1].atoms[2]].id
                << " --> " << allAtoms[allTorsions[itsCurTorsionSeriNum-1].atoms[3]].id
                << " is " << allTorsions[itsCurTorsionSeriNum-1].valueST << std::endl;
    }
    
    void AllSystem::SetOneSP1SP2Bond(int tIdx1, int tIdx2)
    {
        int i1=-1;
        for(std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                iA1 !=allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if(*iA1 !=tIdx2)
            {
                i1= *iA1;
                break;
            }
        }
            
        if (i1 ==-1)
        {
         
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx1].id << " has one connection"
                      << std::endl;
            return;
        }
             
        
        std::vector<int> tV2; 
        // Ring atom first 
        for(std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                iA2 !=allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if(*iA2 !=tIdx1 && AtomsInSameRing(allAtoms[i1], allAtoms[*iA2], allRingsV))
            {
                tV2.push_back(*iA2);
                break;
            }
        }
        
        for(std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                iA2 !=allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if(*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }    
          
        if ((int)tV2.size() < 1)
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx2].id << " has one connection"
                      << std::endl;
            return;
        }
        else if ((int)tV2.size() > 2)
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx2].id << " has mor than 3 connections. Not sp2 "
                      << std::endl;
            return;
        }
        
        REAL va[1][2];
        int per = 1;
        if(AtomsInSameRing(allAtoms[i1], allAtoms[tV2[0]], allRingsV))
        {   
            va[0][0] =        0.0;
            va[0][1] =      180.0;
        }
        else
        {
            va[0][0] =        90.0;
            va[0][1] =       -90.0;
        }
        
        int i=0;        
        for (std::vector<int>::iterator iA2=tV2.begin();
                iA2 !=tV2.end(); iA2++)
        {
            std::vector<int> aTS;
            aTS.push_back(i1);
            aTS.push_back(tIdx1);
            aTS.push_back(tIdx2);
            aTS.push_back(*iA2);
            setOneTorsion(aTS, va[0][i], per);
            i++;
        }
         
    }
   
    void AllSystem::SetOneSP1SP3Bond(int tIdx1, int tIdx2)
    {
        int i1=-1;
        
        for(std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                iA1 !=allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            
            if(*iA1 !=tIdx2)
            {
                i1= *iA1;
                break;
            }
        }
            
        if (i1 ==-1)
        {
         
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx1].id << " has one connection"
                      << std::endl;
            return;
        }
             
        
        std::vector<int> tV2; 
       
        // Ring atom first 
        for(std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                iA2 !=allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if(AtomsInSameRing(allAtoms[i1], allAtoms[*iA2], allRingsV))
            {
                tV2.push_back(*iA2);
                break;
            }
        }
        
        for(std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                iA2 !=allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if(*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }   
        
        
        if ((int)tV2.size() == 0)
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx2].id << " has one connection"
                      << std::endl;
            return;
        }
        else if ((int)tV2.size() > 3)
        {
            std::cout << "No torsion angles for the bond of atom "
                      << allAtoms[tIdx1].id << " and " << allAtoms[tIdx2].id
                      << " because " << allAtoms[tIdx2].id << " has more than 4 connections. Not sp3 "
                      << std::endl;
            return;
        }
        
        REAL va[1][3];
        int per = 3;
        
        if(AtomsInSameRing(allAtoms[i1], allAtoms[tV2[0]], allRingsV))
        {
            va[0][0] =       0.0;
            va[0][1] =     120.0;
            va[0][2] =    -120.0;  
        }
        else
        {
            va[0][0] =     180.0;
            va[0][1] =     -60.0;
            va[0][2] =      60.0;
        }
        
        int i=0; 
        //need to use chiral properties
        if ((int)allAtoms[tIdx2].chiralIdx !=0)
        {
            
        }
        else
        {
            for (std::vector<int>::iterator iA2=tV2.begin();
                iA2 !=tV2.end(); iA2++)
            {
                std::vector<int> aTS;
                aTS.push_back(i1);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(*iA2);
                setOneTorsion(aTS, va[0][i], per);
                i++;
            }    
        }
    }
    
    
    void AllSystem::SetOneSP2SP2Bond(int tIdx1, int tIdx2)
    {
            // two sp2 atoms
            std::vector<int> tV1, tV2;
            
            int tS1=-1, tS2=-1;
        
            for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
            {
                tS1 =-1;
                for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
                {
                    tS2 =-1;
                    if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                    {
                        if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                        {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            // std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //         << std::endl;
                            break;
                        }
                    }
                }
                
                if(tS1 !=-1 && tS2 !=-1)
                {
                    break;
                }
            }
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 != tIdx2 && *iA1 !=tS1)
                {
                    tV1.push_back(*iA1);
                }
            }
            
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && *iA2 !=tS2)
                {
                    tV2.push_back(*iA2);
                }
            }
            
            if (tS1 ==-1 && tS2 ==-1 
                && (int)tV1.size()==2 && (int)tV2.size()==2)
            {
                int tm;
                // check if the ending atom connects H atoms only
                bool tH0=true, tH1=true;
                for (std::vector<int>::iterator iConn=allAtoms[tV1[0]].connAtoms.begin();
                        iConn != allAtoms[tV1[0]].connAtoms.end(); iConn++)
                {
                    if (allAtoms[*iConn].chemType !="H" && *iConn !=tIdx1)
                    {
                        tH0=false;
                        break;
                    }
                }
                for (std::vector<int>::iterator iConn=allAtoms[tV1[1]].connAtoms.begin();
                        iConn != allAtoms[tV1[1]].connAtoms.end(); iConn++)
                {
                    if (allAtoms[*iConn].chemType !="H" && *iConn !=tIdx1 )
                    {
                        tH1=false;
                        break;
                    }
                }
                
                /*
                if (tH0)
                {
                    if (!tH1 ||((int)allAtoms[tV1[0]].connAtoms.size() < (int)allAtoms[tV1[1]].connAtoms.size()))
                    {
                        tm     = tV1[0];
                        tV1[0] = tV1[1];
                        tV1[1] = tm;
                    }
                }
                else if((int)allAtoms[tV1[0]].connAtoms.size() < (int)allAtoms[tV1[1]].connAtoms.size() && !tH1)
                {
                    tm  =tV1[0];
                    tV1[0] = tV1[1];
                    tV1[1] = tm;
                }
                */
                
                if (tH0 && !tH1)
                {
                    tm     = tV1[0];
                    tV1[0] = tV1[1];
                    tV1[1] = tm;
                }
                else if((int)allAtoms[tV1[0]].connAtoms.size() < (int)allAtoms[tV1[1]].connAtoms.size())
                {
                    tm  =tV1[0];
                    tV1[0] = tV1[1];
                    tV1[1] = tm;
                }
                
                
                tH0=true, tH1=true;
                for (std::vector<int>::iterator iConn=allAtoms[tV2[0]].connAtoms.begin();
                        iConn != allAtoms[tV2[0]].connAtoms.end(); iConn++)
                {
                    if (allAtoms[*iConn].chemType !="H"  && *iConn !=tIdx2)
                    {
                        tH0=false;
                        break;
                    }
                }
                for (std::vector<int>::iterator iConn=allAtoms[tV2[1]].connAtoms.begin();
                        iConn != allAtoms[tV2[1]].connAtoms.end(); iConn++)
                {
                    if (allAtoms[*iConn].chemType !="H" && *iConn !=tIdx2)
                    {
                        tH1=false;
                        break;
                    }
                }
                
                // std::cout << "tV2[0]: " << allAtoms[tV2[0]].id << " tV2[1] " << allAtoms[tV2[1]].id << std::endl;
                if (tH0)
                {
                    if (!tH1 ||((int)allAtoms[tV2[0]].connAtoms.size() < (int)allAtoms[tV2[1]].connAtoms.size()))
                    {
                        tm     = tV2[0];
                        tV2[0] = tV2[1];
                        tV2[1] = tm;
                    }
                }
                else if((int)allAtoms[tV2[0]].connAtoms.size() < (int)allAtoms[tV2[1]].connAtoms.size() && !tH1)
                {
                    tm  =tV2[0];
                    tV2[0] = tV2[1];
                    tV2[1] = tm;
                }
                // std::cout << "tV2: " << " tH0-> " << tH0 << " tH1 " << tH1 << std::endl;
            }
            
            
            
            REAL va[2][2];
            int per = 2;
            if (tS1 !=-1 && tS2 !=-1)
            {
                va[0][0] =     0.0;
                va[0][1] =   180.0;
                va[1][0] =   180.0;
                va[1][1] =     0.0;
            }
            else
            {
                va[0][0] =  180.0;
                va[0][1] =    0.0;
                va[1][0] =    0.0;
                va[1][1] =  180.0;
            }
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            } 
            
            /*
            std::cout << "id1 " << allAtoms[tIdx1].id 
                  << " and id2 " << allAtoms[tIdx2].id << std::endl;
            
            if ((allAtoms[tIdx1].id=="C14" && allAtoms[tIdx2].id=="C23")
                 || (allAtoms[tIdx2].id=="C23" && allAtoms[tIdx1].id=="C14"))
            {
                std::cout << "tS1 = " << tS1 << std::endl
                          << "tS2 = " << tS2 << std::endl;
                
            }
            */
        /*   
        std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
         */
    }
    
    void AllSystem::SetOneSP2SP3Bond(int tIdx1, int tIdx2, std::string tF)
    {
                
        
        // one sp2 atom and one sp3 atom
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                            
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                
                break;
                
            }
        }
            
        REAL va1  =   0.0;
        REAL va2  =  60.0;
        int  per  =     6;
        
        REAL va[2][3];
        if (tS1 !=-1 && tS2 !=-1)
        {
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
            /*
            if (tF=="even")
            {
                va[0][0] =    va2;
                va[0][1] =  3*va2;
                va[0][2] =   -va2;
        
                va[1][0] =  -2*va2;
                va[1][1] =     va1;
                va[1][2] =   2*va2;
            }
            else
            {
                va[0][0] =   -va2;
                va[0][1] =    va2;
                va[0][2] =   3*va2;
        
                va[1][0] =   2*va2;
                va[1][1] =  -2*va2;
                va[1][2] =     va1;
            }
             */
        }
        else
        {  /*
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
           */  
            
            va[0][0] =  3*va2;
            va[0][1] =  -va2;
            va[0][2] =   va2;
            
            va[1][0] =    va1;
            va[1][1] =  2*va2;
            va[1][2] = -2*va2;
        }
        
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 != tIdx2 && *iA1 != tS1)
            {
                tV1.push_back(*iA1);
            }
        }
        
        // for sp3 atoms
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            
            /*
            if(tS1 !=-1 && tS2 !=-1)
            {
                std::cout << "Chiral center " << allAtoms[tIdx2].id << std::endl;
                std::cout << "Chiral atoms\n";
                for (int i=0; i < (int)allChirals[iCh].atoms.size(); i++)
                {
                    std::cout << allAtoms[allChirals[iCh].atoms[i]].id << std::endl;
                }
                exit(1);
                        
            }
            */
            // buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
           /*
            //if(allAtoms[tIdx2].id =="C15" || allAtoms[tIdx2].id =="C9")
            //{
             std::cout << "Chiral serial number is " << iCh << std::endl;   
             std::cout << "root atom is  " << allAtoms[tIdx1].id << std::endl;
             std::cout << "atom  " << allAtoms[tIdx2].id << " mutable size is "
                      << (int)allChirals[iCh].mutTable.size() << std::endl
                      << std::endl << allChirals[iCh].mutTable[tIdx1][0]
                      << " and " << allChirals[iCh].mutTable[tIdx1][1] << std::endl
                      << " chiral atom seq is " << std::endl; 
             for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
             {
                 std::cout << "atom " << allAtoms[*iA].id << std::endl;
             }
           
            //}
            */
        }
        else
        {
            int nH = getNumSpecAtomConnect(tIdx2, "H");
            if (allAtoms[tIdx1].chemType.compare("H") !=0 && nH==2)
            {   
  
                for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                        iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
                {
                    if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                          allAtoms[*iA2].chemType=="H")
                    {
                        tV2.push_back(*iA2);
                        break;
                    }
                }
                
                for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                        iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
                {
                    if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                          allAtoms[*iA2].chemType!="H")
                    {
                        tV2.push_back(*iA2);
                        break;
                    }
                }
            }
        }
        
        // put those atoms which are not, such as H atoms, in chiral 
        // list into the output list
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        if ((int)tV1.size() <=(int)tV2.size())
        {
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }
        }
        else
        {
            for (int i =0; i < (int)tV2.size(); i++)
            {
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV2[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV1[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }   
        }
        
        /*
        std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        */
       
    }
    
    
    void AllSystem::SetOneSP2SP3Bond(int tIdx1, int tIdx2)
    {
                
        
        
        // one sp2 atom and one sp3 atom
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        std::vector<int> tS3;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            if(*iAt1 !=tIdx2)
            {
                tS3.push_back(*iAt1);
                tS1 =-1;
                for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                           iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
                {
                    tS2 =-1;
                    if (*iAt2 !=tIdx1)
                    {
                        if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                        {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;   
                        }
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                
                break;
                
            }
        }
            
        REAL va1;
        REAL va2;
        int  per  =     6;
        
        REAL va[2][3];
        if (tS1 !=-1 && tS2 !=-1)
        {
            va1  =   0.0;
            va2  =  60.0;
            
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
            /*
            va[0][0] =    va2;
            va[0][1] =  3*va2;
            va[0][2] =   -va2;
        
            va[1][0] =  -2*va2;
            va[1][1] =     va1;
            va[1][2] =   2*va2;
             */
        }
        else if ((int)tS3.size()==2 && (AtomsInSameRing(allAtoms[tS3[0]], allAtoms[tS3[1]], allRingsV)))
        {
            va1      =   90.0;
            va2      =   30.0;
            
            va[0][0] =  2*va1-va2;
            va[0][1] =  -va1;
            va[0][2] =   va2;
        
            va[1][0] =  -va2;
            va[1][1] =   va1;
            va[1][2] = -(2*va1-va2);
        }
        else
        {
            va1  =   0.0;
            va2  =  60.0;
            
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
            
            /*
            
           
            va[0][0] =    va1;
            va[0][1] =  -(2*va1-va2);
            va[0][2] =   -va2;
        
            va[1][0] =  -va1;
            va[1][1] =   va2;
            va[1][2] =   2*va1-va2;
            
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
             ////////////////
            va[0][0] =  3*va2;
            va[0][1] =  -va2;
            va[0][2] =   va2;
            
            va[1][0] =    va1;
            va[1][1] =  2*va2;
            va[1][2] = -2*va2;
            */
        }
        
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 != tIdx2 && *iA1 != tS1)
            {
                tV1.push_back(*iA1);
            }
        }
        
        // for sp3 atoms
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            
            /*
            if(tS1 !=-1 && tS2 !=-1)
            {
                std::cout << "Chiral center " << allAtoms[tIdx2].id << std::endl;
                std::cout << "Chiral atoms\n";
                for (int i=0; i < (int)allChirals[iCh].atoms.size(); i++)
                {
                    std::cout << allAtoms[allChirals[iCh].atoms[i]].id << std::endl;
                }
                exit(1);
                        
            }
            */
            // buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
           /*
            //if(allAtoms[tIdx2].id =="C15" || allAtoms[tIdx2].id =="C9")
            //{
             std::cout << "Chiral serial number is " << iCh << std::endl;   
             std::cout << "root atom is  " << allAtoms[tIdx1].id << std::endl;
             std::cout << "atom  " << allAtoms[tIdx2].id << " mutable size is "
                      << (int)allChirals[iCh].mutTable.size() << std::endl
                      << std::endl << allChirals[iCh].mutTable[tIdx1][0]
                      << " and " << allChirals[iCh].mutTable[tIdx1][1] << std::endl
                      << " chiral atom seq is " << std::endl; 
             for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
             {
                 std::cout << "atom " << allAtoms[*iA].id << std::endl;
             }
           
            //}
            */
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType=="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
        
        // put those atoms which are not, such as H atoms, in chiral 
        // list into the output list
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        if ((int)tV1.size() <=(int)tV2.size())
        {
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }
        }
        else
        {
            for (int i =0; i < (int)tV2.size(); i++)
            {
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV2[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV1[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }   
        }
        /*
        std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
       */
       
    }
    
    
    
    void AllSystem::SetOneSP3SP3Bond(int tIdx1, int tIdx2)
    {
        
        REAL va[3][3];
        
        
        
        
        /*
        va[0][0] = 60.0;
        va[0][1] = 180.0;
        va[0][2] = -60.0;
        
        va[1][0] = -60.0;
        va[1][1] =  60.0;
        va[1][2] = 180.0;
        
        va[2][0] = 180.0;
        va[2][1] = -60.0;
        va[2][2] =  60.0;
        */
        
        
        va[0][0] = 180.0;
        va[0][1] = -60.0;
        va[0][2] =  60.0;
        
        va[1][0] =  60.0;
        va[1][1] = 180.0;
        va[1][2] = -60.0;
        
        va[2][0] = -60.0;
        va[2][1] =  60.0;
        va[2][2] = 180.0;
       
        
        int  per =     3;
        
        // One sp3 atom and one sp3 atom
        // Several procedures
        // 1. if there are two atoms in the same ring
        
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        // 2. take consider of chirals 
        // std::cout << "tS1 " << tS1 << " tS2 " << tS2 << std::endl;
        
        if ((int)allAtoms[tIdx1].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx1].inChirals[0];
            
            // buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
        }
        else
        {
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType !="H")
                {
                    tV1.push_back(*iA1);
                    break;
                }
            }
        }
        
       
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            //buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType !="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
       
        
       
           
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        /*
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        */
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
    }
    
    void AllSystem::SetOneSP3SP3Bond(int tIdx1, int tIdx2, std::string tF)
    {
        
        REAL va[3][3];
        if (tF=="even")
        {
            va[0][0] =  60.0;
            va[0][1] = 180.0;
            va[0][2] = -60.0;
        
            va[1][0] = -60.0;
            va[1][1] =  60.0;
            va[1][2] = 180.0;
        
            va[2][0] = 180.0;
            va[2][1] = -60.0;
            va[2][2] =  60.0;
        }
        else if(tF=="odd")
        {
            va[0][0] = -60.0;
            va[0][1] =  60.0;
            va[0][2] = 180.0;
        
            va[1][0] = 180.0;
            va[1][1] = -60.0;
            va[1][2] =  60.0;
        
            va[2][0] =  60.0;
            va[2][1] = 180.0;
            va[2][2] = -60.0;
        }
        else
        {
            std::cout << "what is the sequence idx of this torsion in a ring, even or odd "
                    << std::endl;
            exit(1);
        }
        
        int  per =     3;
        
        // One sp3 atom and one sp3 atom
        // Several procedures
        // 1. if there are two atoms in the same ring
        
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            // std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        // 2. take consider of chirals 
        // std::cout << "tS1 " << tS1 << " tS2 " << tS2 << std::endl;
        
        if ((int)allAtoms[tIdx1].inChirals.size() !=0)
        {
            
            int iCh = allAtoms[tIdx1].inChirals[0];
            
            // buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
        }
        else
        {
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType !="H")
                {
                    tV1.push_back(*iA1);
                    break;
                }
            }
        }
        
       
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            //buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType !="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
        
           
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        /*
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        */
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
    }
    
    void AllSystem::SetOneSP3SP3Bond4H(int tIdx1, int tIdx2)
    {
        // Each sp3 center atoms connected two H atoms
        // Two sp3 center atoms are not in the same ring
        
        REAL va[3][3];
        va[0][0] = 180.0;
        va[0][1] = -60.0;
        va[0][2] =  60.0;
        
        va[1][0] =  60.0;
        va[1][1] =  180.0;
        va[1][2] = -60.0;
        
        va[2][0] = -60.0;
        va[2][1] =  60.0;
        va[2][2] = 180.0;
        
        int  per =     3;
        
        std::vector<int> tV1, tV2;
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType !="H")
            {
                tV1.push_back(*iA1);
                break;
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType !="H")
            {
                tV2.push_back(*iA2);
                break;
            }
        }
        
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        
        /*
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
         */
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
   
    }
    
    bool AllSystem::checkSP3SP34H(int tIdx1, int tIdx2)
    {
        int  n1  = 0;
        int  n2  = 0;
        bool h4  = false;
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && allAtoms[*iA1].chemType =="H")
            {
                n1++;
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2!=tIdx1 && allAtoms[*iA2].chemType =="H")
            {
                n2++;
            }
        }
        
        if (n1==2 && n2==2)
        {
            h4 = true;
        }
        
        return h4;
    }
    
    void AllSystem::SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                             int tPer, REAL tIniValue)
    {
        // atoms connected to atoms tIdx1, tIdx2
        int iS=-1;
        std::vector<int> tV1, tV2;
        int tS1=-1, tS2=-1;
        //std::cout << "For atom " << allAtoms[tIdx1].id << " and atom " 
        //          << allAtoms[tIdx2].id << "connected atoms are " << std::endl;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            //std::cout << "1 linked " << allAtoms[*iAt1].id << std::endl;
            
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                // std::cout << "2 linked " << allAtoms[*iAt2].id << std::endl;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        if(tS1 ==-1 && tS2 ==-1)
        {
            //First find a atom which is not H as a start atom
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if(allAtoms[*iA1].chemType !="H" && *iA1 !=tIdx2)
                {
                    iS=*iA1;
                    break;
                }
            }
            
            // if in some cases, only H has involved, use them.
            if (iS==-1)
            {
                iS=allAtoms[tIdx1].connAtoms[0];
            }
            
            tV1.push_back(iS);
        }
        
            if ((int)allAtoms[tIdx1].inChirals.size() !=0)
            {
                int iCh = allAtoms[tIdx1].inChirals[0];
                //buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
                //std::cout << "buildChiralCluster2 " << std::endl;
                buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
            }
            
        // std::cout << "tV1 " << (int)tV1.size() << std::endl;
        /*
        for(std::vector<int>::iterator iV1=tV1.begin(); iV1 != tV1.end();
                iV1++)
        {
            std::cout << allAtoms[*iV1].id  << std::endl;
        }
         */ 
                
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
                {
                    tV1.push_back(*iA1);
                }
            }
           
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
                {
                    tV2.push_back(*iA2);
                }
            }
            
           // std::cout << "number of tV2 " << (int)tV2.size() << std::endl;
            
            REAL perValue = 360.0/tPer;
            REAL iniValue;
            
            int n = 0;
            
            for (int i=0; i <(int)tV1.size(); i++)
            {
                if(tS1 !=-1 && tS2 !=-1)
                {
                   iniValue= 60.0 +  n*perValue;
                }
                else
                {
                    iniValue=tIniValue + n*perValue;
                }
                // std::cout << "iniValue " << iniValue << std::endl;
                
                int m =0;
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    // std::cout << "tV1 " << i << " value " << allAtoms[tV1[i]].id << std::endl;
                    // std::cout << "tV2 " << j << " value " << allAtoms[tV2[j]].id << std::endl;
                    std::vector<int> TS;
                    TS.push_back(tV1[i]);
                    TS.push_back(tIdx1);
                    TS.push_back(tIdx2);
                    TS.push_back(tV2[j]);
                    
                 
                    REAL curValue = iniValue +m*perValue;
                  
                    if (curValue > 360.0)
                    {
                        curValue = curValue -360.0;
                    }
                    else if (curValue > 180.0)
                    {
                        curValue = curValue -360.0;
                    }
                    setOneTorsion(TS, curValue, tPer);
                    //std::cout << "atom 1 " << allAtoms[tV1[i]].id
                    //          << " atom 4 " << allAtoms[tV2[j]].id
                    //          << " tor " << curValue << std::endl;
                    m++;
                }
                n++;
            }
            
            /*
            std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
            for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
            {
                std::cout << "atom " << allAtoms[*iA].id << std::endl;
            }
                
            std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                      << allAtoms[tIdx2].id << " form torsion using " << std::endl;
            for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
            {
                std::cout << "atom " << allAtoms[*iA].id << std::endl;
            }
             */
            
    }
    
    void AllSystem::SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                             int tPer, REAL tIniValue, 
                                             std::string tF)
    {
        // atoms connected to atoms tIdx1, tIdx2
        int iS=-1;
        std::vector<int> tV1, tV2;
        int tS1=-1, tS2=-1;
        //std::cout << "For atom " << allAtoms[tIdx1].id << " and atom " 
        //          << allAtoms[tIdx2].id << "connected atoms are " << std::endl;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            //std::cout << "1 linked " << allAtoms[*iAt1].id << std::endl;
            
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                //std::cout << "2 linked " << allAtoms[*iAt2].id << std::endl;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        if(tS1 ==-1 && tS2 ==-1)
        {
            //First find a atom which is not H as a start atom
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if(allAtoms[*iA1].chemType !="H" && *iA1 !=tIdx2)
                {
                    iS=*iA1;
                    break;
                }
            }
            
            // if in some cases, only H has involved, use them.
            if (iS==-1)
            {
                iS=allAtoms[tIdx1].connAtoms[0];
            }
            
            tV1.push_back(iS);
        }
        
            if ((int)allAtoms[tIdx1].inChirals.size() !=0)
            {
                int iCh = allAtoms[tIdx1].inChirals[0];
                //buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
                //std::cout << "buildChiralCluster2 " << std::endl;
                buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
            }
            
        // std::cout << "tV1 " << (int)tV1.size() << std::endl;
        /*
        for(std::vector<int>::iterator iV1=tV1.begin(); iV1 != tV1.end();
                iV1++)
        {
            std::cout << allAtoms[*iV1].id  << std::endl;
        }
         */
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
                {
                    tV1.push_back(*iA1);
                }
            }
           
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
                {
                    tV2.push_back(*iA2);
                }
            }
            
           // std::cout << "number of tV2 " << (int)tV2.size() << std::endl;
            
            REAL perValue = 360.0/tPer;
            REAL iniValue;
            
            int n = 0;
            
            for (int i=0; i <(int)tV1.size(); i++)
            {
                if(tS1 !=-1 && tS2 !=-1)
                {
                    if (tF=="even")
                    {
                        iniValue= 60.0 +  n*perValue;
                    }
                    else if (tF=="old")
                    {
                        iniValue= -60.0 +  n*perValue;
                    }
                }
                else
                {
                    iniValue=tIniValue + n*perValue;
                }
                // std::cout << "iniValue " << iniValue << std::endl;
                
                int m =0;
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    // std::cout << "tV1 " << i << " value " << allAtoms[tV1[i]].id << std::endl;
                    // std::cout << "tV2 " << j << " value " << allAtoms[tV2[j]].id << std::endl;
                    std::vector<int> TS;
                    TS.push_back(tV1[i]);
                    TS.push_back(tIdx1);
                    TS.push_back(tIdx2);
                    TS.push_back(tV2[j]);
                    
                 
                    REAL curValue = iniValue +m*perValue;
                  
                    if (curValue > 360.0)
                    {
                        curValue = curValue -360.0;
                    }
                    else if (curValue > 180.0)
                    {
                        curValue = curValue -360.0;
                    }
                    setOneTorsion(TS, curValue, tPer);
                    //std::cout << "atom 1 " << allAtoms[tV1[i]].id
                     //         << " atom 4 " << allAtoms[tV2[j]].id
                     //         << " tor " << curValue << std::endl;
                    m++;
                }
                n++;
            }
            
            /*
            std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
             */
    }
    
    
    void AllSystem::setTorsionFromOneBond(int tIdx1, int tIdx2)
    {
        int bIdx1 =  allAtoms[tIdx1].bondingIdx;
        ID  id1   =  allAtoms[tIdx1].chemType;
        
        int bIdx2 =  allAtoms[tIdx2].bondingIdx;
        ID  id2   =  allAtoms[tIdx2].chemType;
        /*
        std::cout << "id1 "    << allAtoms[tIdx1].id << " serial num " << tIdx1 
                  << " bIdx1 " << bIdx1
                  << " and id2 " << allAtoms[tIdx2].id  << " serial num " << tIdx2
                  << " bIdx2 "   << bIdx2 << std::endl;
        */
        std::vector<ID> OxyCol;
        
        OxyCol.push_back("O");
        OxyCol.push_back("S");
        OxyCol.push_back("Se");
        OxyCol.push_back("Te");
        OxyCol.push_back("Po");
        
        std::vector<ID>::iterator iFind1 = std::find(OxyCol.begin(), OxyCol.end(), id1);
        std::vector<ID>::iterator iFind2 = std::find(OxyCol.begin(), OxyCol.end(), id2);
        
        if ((iFind1 != OxyCol.end() && bIdx1==2) && 
                (iFind2 != OxyCol.end() && bIdx2==3))
        {
            // two Oxy column atoms with sp3 orbiting
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 2, 90.0)" << std::endl;
            //std::cout << "atom 1 "  << allAtoms[tIdx1].id 
            //          << " atom 2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 90.0);   
        }
        else if ((iFind1 != OxyCol.end() && bIdx1==2) && bIdx2==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 12, 180.0)" << std::endl;
            //std::cout << "atom 1 sp2 " << allAtoms[tIdx1].id 
            //        << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 180.0);
            
        }
        else if ((iFind2 != OxyCol.end() && bIdx2==2) && bIdx1==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond" << std::endl;
            // std::cout << "atom 1 sp3 " << allAtoms[tIdx1].id 
            //          << "atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx1, tIdx2, 3, 180.0);
        }
        else if (bIdx1==2 && bIdx2==2)
        {
            //std::cout << "SetOneSP2SP2Bond" << std::endl;
            SetOneSP2SP2Bond(tIdx1, tIdx2);
            
        }
        else if ((bIdx1==2 && bIdx2==3) || (bIdx1==3 && bIdx2==2))
        {
            if (bIdx1==2 && bIdx2==3)
            {
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                // std::cout << " atom 1 sp2 " << allAtoms[tIdx1].id 
                //          << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx1, tIdx2);
            }
            else
            {  
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                // std::cout << " atom 1 sp3 " << allAtoms[tIdx1].id 
                //           << " atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx2, tIdx1);
            }
        }
        else if (bIdx1==3 && bIdx2==3)
        {
            
            if (!AtomsInSameRing(allAtoms[tIdx1], allAtoms[tIdx2], allRingsV)
                 && checkSP3SP34H(tIdx1, tIdx2))
            {
                //std::cout << "SetOneSP3SP3Bond4H" << std::endl;
                SetOneSP3SP3Bond4H(tIdx1, tIdx2);
            }
            else
            {
                //std::cout << "SetOneSP3SP3Bond" << std::endl;
                SetOneSP3SP3Bond(tIdx1, tIdx2);
            }
        }
        else if (bIdx1==1 && bIdx2==1)
        {
            SetOneSP1SP1Bond(tIdx1, tIdx2);
        }
        else if(bIdx1==0 && bIdx2==0)
        {
            SetOneSP1SP1Bond(tIdx1, tIdx2);
        }
        else if((bIdx1==0 || bIdx1==1) && bIdx2==2)
        {
            SetOneSP1SP2Bond(tIdx1, tIdx2);
        }
        else if((bIdx2==0||bIdx2==1) && bIdx1==2)
        {
            SetOneSP1SP2Bond(tIdx2, tIdx1);
        }
       
        else if((bIdx1==0 ||bIdx1==1) && bIdx2==3)
        {
            SetOneSP1SP3Bond(tIdx1, tIdx2);
        }
        else if((bIdx2==0||bIdx2==1) && bIdx1==3)
        {
            SetOneSP1SP3Bond(tIdx2, tIdx1);
        }
        
        
    }
     
    void AllSystem::setTorsionFromOneBond(int tIdx1, int tIdx2, std::string tF)
    {
        int bIdx1 =  allAtoms[tIdx1].bondingIdx;
        ID  id1   =  allAtoms[tIdx1].chemType;
        
        int bIdx2 =  allAtoms[tIdx2].bondingIdx;
        ID  id2   =  allAtoms[tIdx2].chemType;
        /*
        std::cout << "id1 "    << allAtoms[tIdx1].id << " serial num " << tIdx1 
                  << " bIdx1 " << bIdx1
                  << " and id2 " << allAtoms[tIdx2].id  << " serial num " << tIdx2
                  << " bIdx2 "   << bIdx2 << std::endl;
        */
        std::vector<ID> OxyCol;
        
        OxyCol.push_back("O");
        OxyCol.push_back("S");
        OxyCol.push_back("Se");
        OxyCol.push_back("Te");
        OxyCol.push_back("Po");
        
        std::vector<ID>::iterator iFind1 = std::find(OxyCol.begin(), OxyCol.end(), id1);
        std::vector<ID>::iterator iFind2 = std::find(OxyCol.begin(), OxyCol.end(), id2);
        
        if ((iFind1 != OxyCol.end() && bIdx1==2) && 
                (iFind2 != OxyCol.end() && bIdx2==3))
        {
            // two Oxy column atoms with sp3 orbiting
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 2, 90.0)" << std::endl;
            //std::cout << "atom 1 "  << allAtoms[tIdx1].id 
            //          << " atom 2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 90.0);   
        }
        else if ((iFind1 != OxyCol.end() && bIdx1==2) && bIdx2==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 12, 180.0)" << std::endl;
            // std::cout << "atom 1 sp2 " << allAtoms[tIdx1].id 
            //         << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 180.0, tF);
            
        }
        else if ((iFind2 != OxyCol.end() && bIdx2==2) && bIdx1==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond" << std::endl;
            // std::cout << "atom 1 sp3 " << allAtoms[tIdx1].id 
            //           << "atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx1, tIdx2, 3, 180.0, tF);
        }
        else if (bIdx1==2 && bIdx2==2)
        {
            //std::cout << "SetOneSP2SP2Bond" << std::endl;
            SetOneSP2SP2Bond(tIdx1, tIdx2);
            
        }
        else if ((bIdx1==2 && bIdx2==3) || (bIdx1==3 && bIdx2==2))
        {
            if (bIdx1==2 && bIdx2==3)
            {
                // std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                // std::cout << " atom 1 sp2 " << allAtoms[tIdx1].id 
                //           << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx1, tIdx2);
            }
            else
            {  
                // std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                // std::cout << " atom 1 sp3 " << allAtoms[tIdx1].id 
                //           << " atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx2, tIdx1);
            }
        }
        else if (bIdx1==3 && bIdx2==3)
        {
            // std::cout << "SetOneSP3SP3Bond" << std::endl;
            SetOneSP3SP3Bond(tIdx1, tIdx2, tF);
        }
        else if(bIdx1==1 && bIdx2==1)
        {
            SetOneSP1SP1Bond(tIdx1, tIdx2);
        }
        else if(bIdx1==0 && bIdx2==0)
        {
            SetOneSP1SP1Bond(tIdx1, tIdx2);
        }
        else if((bIdx1==0 || bIdx1==1) && bIdx2==2)
        {
            SetOneSP1SP2Bond(tIdx1, tIdx2);
        }
        else if((bIdx2==0 ||bIdx2==1) && bIdx1==2)
        {
            SetOneSP1SP2Bond(tIdx2, tIdx1);
        }  
        else if(bIdx1==0  && bIdx2==3)
        {
            SetOneSP1SP3Bond(tIdx1, tIdx2);
        }
        else if(bIdx2==0 && bIdx1==3)
        {
            SetOneSP1SP3Bond(tIdx2, tIdx1);
        }
    }
    
    void AllSystem::setAllTorsions()
    {
        // loop over all bonds to get the torsion angles
        for (std::vector<BondDict>::iterator iABo= allBonds.begin();
                iABo != allBonds.end(); iABo++)
        {
            if ((int)iABo->fullAtoms.size() ==2)
            {
                std::vector<int> tPos;
                for(std::map<ID, int>::iterator iAM=iABo->fullAtoms.begin();
                        iAM !=iABo->fullAtoms.end(); iAM++)
                {
                    //std::cout << "atom " << iAM->first << std::endl;
                    // std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                    //        << " atoms " << std::endl;
                    if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                    {
                        tPos.push_back(iAM->second);
                    }
                }
                // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
                if((int)tPos.size() ==2)
                {
                    setTorsionFromOneBond(tPos[0], tPos[1]);
                }
            }
        }
        
        std::cout << "All torsions have been setup " << std::endl;
       
    }
        
    void AllSystem::setAllTorsions2()
    {
        
        std::vector<int>  tDone;
        
        // First set all torsion within all rings
        for (std::vector<RingDict>::iterator iR=allRingsV.begin();
                iR != allRingsV.end(); iR++)
        {
            setAllTorsionsInOneRing(tDone, *iR);
        }   
        
        // find all torsion not involved rings
        for (std::vector<BondDict>::iterator iABo= allBonds.begin();
                iABo != allBonds.end(); iABo++)
        {
            if (std::find(tDone.begin(), tDone.end(), iABo->seriNum)==tDone.end())
            {
                std::vector<int> tPos;
                for(std::map<ID, int>::iterator iAM=iABo->fullAtoms.begin();
                        iAM !=iABo->fullAtoms.end(); iAM++)
                {
                    //std::cout << "atom " << iAM->first << std::endl;
                    //std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                    //          << " atoms " << std::endl;
                    if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                    {
                        tPos.push_back(iAM->second);
                    }
                }
                // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
                if((int)tPos.size() ==2)
                {
                    //std::cout << "Set torsion angles around the bond of atom "
                    //          << allAtoms[tPos[0]].id << " and "
                    //          << allAtoms[tPos[1]].id << std::endl;
                    
                    setTorsionFromOneBond(tPos[0], tPos[1]);
                }
            }
        }
        
        std::cout << "All torsions have been setup " << std::endl;
       
    }
    
    void AllSystem::setAllTorsionsInOneRing(std::vector<int> & tBs, 
                                              RingDict & tR)
    {
        
        std::vector<int> tAs, tLinkA, tBos;
        
        // A list of idx of atoms in the ring
        // std::cout << " atoms in the ring are: " << std::endl;
        for (int i=0; i < (int)tR.atoms.size(); i++)
        {
            tAs.push_back(tR.atoms[i].seriNum);
            // std::cout << tR.atoms[i].seriNum << std::endl;
        }
        
        tLinkA.push_back(tR.atoms[0].seriNum);
        
        int iCur   =tR.atoms[0].seriNum;
        int iLoop  =1;
        // std::cout << "ring rep " << tR.rep << " size " << (int)tAs.size() << std::endl;
       
        
        while ((int)tLinkA.size() < (int)tAs.size()
               && iLoop < (int)tAs.size())
        {   
            // std::cout << "atom " << allAtoms[iCur].seriNum << std::endl;
            for (std::vector<int>::iterator iC=allAtoms[iCur].connAtoms.begin();  
                    iC!=allAtoms[iCur].connAtoms.end(); iC++)
            {
                // std::cout << "connection " << *iC << std::endl;
                if(std::find(tAs.begin(), tAs.end(), *iC) !=tAs.end()
                   && std::find(tLinkA.begin(), tLinkA.end(), *iC) ==tLinkA.end())
                {
                    
                    int iB=getBond(allBonds, allAtoms[iCur].seriNum, allAtoms[*iC].seriNum);
                    if (iB >=0)
                    {
                        // std::cout << "find " << *iC << std::endl;
                        tBos.push_back(iB);
                        tLinkA.push_back(*iC);
                        iCur=*iC;
                        break;
                    }
                }
            }
            
            //double protection 
            iLoop++;
        }
        // Last bond
        /*
        int idxL = (int)tLinkA.size()-1;
        int iB=getBond(allBonds, allAtoms[tLinkA[0]].seriNum, 
                       allAtoms[tLinkA[idxL]].seriNum);
        if (iB >=0 && std::find(tBos.begin(), tBos.end(), iB)==tBos.end())
        {
            tBos.push_back(iB);
        }
        */
        if ((int)tLinkA.size() != (int)tR.atoms.size())
        {
            /*
            std::cout << "could not all linked atoms in ring " << tR.rep << std::endl;
            std::cout << "atoms in the ring are : " << std::endl;
            for (std::vector<AtomDict>::iterator iA=tR.atoms.begin(); 
                    iA !=tR.atoms.end(); iA++)
            {
                std::cout << iA->id << ",\t";
            }
            std::cout << std::endl;
            std::cout << "the links (bonds) found: " << std::endl;
             
            for (std::vector<int>::iterator iL=tLinkA.begin(); 
                    iL != tLinkA.end(); iL++)
            {
                std::cout << allAtoms[*iL].id << std::endl;
            }
             */
        }
        //std::cout << "ring reps : " << tR.rep << std::endl;
        //std::cout << "it has " << (int)tBos.size() << std::endl;
        // Now we have bonds in the ring in sequence
        for (int i=0; i < (int)tBos.size(); i++)
        {
            std::string flip;
            if (i%2==0)
            {
                flip = "even";
            }
            else
            {
                flip = "odd";
            }
            
            std::vector<int> tPos;
            for(std::map<ID, int>::iterator iAM=allBonds[tBos[i]].fullAtoms.begin();
                        iAM !=allBonds[tBos[i]].fullAtoms.end(); iAM++)
            {
                //std::cout << "atom " << iAM->first << std::endl;
                // std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                //        << " atoms " << std::endl;
                if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                {
                    tPos.push_back(iAM->second);
                }
            }
            // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
            if((int)tPos.size() ==2)
            {    
                //std::cout << "set torsion angles for the bond of atoms "
                //          << allAtoms[tPos[0]].id << " and " 
                //          << allAtoms[tPos[1]].id 
                //          << " with " << flip << std::endl;
                       
                setTorsionFromOneBond(tPos[0], tPos[1], flip);
                tBs.push_back(tBos[i]);
            }
        }
        
    }
    
    void AllSystem::setTorsionIdxFromOneBond(int tIdx1, int tIdx2)
    {
                
        std::cout << "For the bond consisting of atoms  " << allAtoms[tIdx1].id 
                  << " and " << allAtoms[tIdx2].id << std::endl
                  << "It has following torsion angles: " << std::endl; 
                
                
       for (std::vector<int>::iterator iAt1= allAtoms[tIdx1].connAtoms.begin();
            iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
       {
           for (std::vector<int>::iterator iAt2 = allAtoms[tIdx2].connAtoms.begin();
                iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
           {
               if (*iAt1 != tIdx2 && *iAt2 != tIdx1)
               {
                   TorsionDict aTorsion;         
                   aTorsion.atoms.push_back(*iAt1);
                   aTorsion.atoms.push_back(tIdx1);
                   aTorsion.atoms.push_back(tIdx2);
                   aTorsion.atoms.push_back(*iAt2);
                   std::cout << "Torsion: " << allAtoms[*iAt1].id 
                             << ", " << allAtoms[tIdx1].id 
                             << ", " << allAtoms[tIdx2].id
                             << ", " << allAtoms[*iAt2].id << std::endl; 
                                    
                   allTorsions.push_back(aTorsion);
               }
           }
        }
        
    }
    
        // Ring related 
    void AllSystem::ringDetecting()
    {
        // tempo, set ring size to 6
        int maxSize = 7;
        // std::vector<AtomDict> atomsInPath;
        std::map<int, ID>  atomsInPath;
        std::map<int, ID>  atomsSeen;
       
        // 1. loops beginning from all atoms in the system
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
   
            //std::cout << "-----------------" << std::endl;
            //std::cout << "starting atom  " << iA->id << std::endl;
            //std::cout << "-----------------" << std::endl;
            
            //if((int)tempIDs.size() !=0)
            //{
            //    tempIDs.clear();
            //}
            
            int preSeriNum = -999;
            int startLev   = 1;
            atomsInPath.clear();
            atomsSeen.clear();
            // atomsInPath.push_back((*iA));
            // atomsInPath.insert(std::pair<int, ID>(iA->seriNum, iA->chemType))
            // tempIDs.insert(std::pair<int, ID>(iA->seriNum, iA->id));
            
            // 2. loop from its bonded atoms recursively
            // checkOnePathSec(atomsInPath, tempIDs, *iA, maxSize, iA);
            // checkOnePathSec(atomsInPath, *iA, maxSize, iA);
            checkOnePathSec(*iA, maxSize, iA, preSeriNum,  startLev, atomsSeen, atomsInPath);
            
            /*
            if ((int)iA->ringRep.size() >0)
            {
                
                std::cout << "\nAtom " << iA->id << " is in "
                        << (int)iA->ringRep.size() << " rings "
                        << std::endl;
                for (std::map<std::string, int>::iterator iMa = iA->ringRep.begin();
                       iMa != iA->ringRep.end(); iMa++ )
                {
                    std::cout << "Ring: " << iMa->first << "; Size " << iMa->second 
                            <<std::endl;
                } 
            } 
            else
            {
                std::cout << "\nAtom " << iA->id << " is in no ring" <<std::endl;
            }
            */
           
            // std::cout << "finish atom  " << iA->seriNum << std::endl;
        }
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if ((int)iA->ringRep.size() >0)
            {
                std::cout << "\nAtom " << iA->id << " is in  " 
                        << (int)iA->ringRep.size() << " ring(s)" << std::endl;
                for (std::map<std::string, int>::iterator iRR=iA->ringRep.begin();
                        iRR != iA->ringRep.end(); iRR++)
                {
                    std::cout << "Ring: " << iRR->first << ", size : "
                            << iRR->second << std::endl;
                }
            }
        }
    }
    
 
    // Version 2
    void AllSystem::checkOnePathSec(AtomDict                & curAto,
                                      int                       iMax,
                                      std::vector<AtomDict>::iterator iOriAto,
                                      int                       SeriNumPreAto,  
                                      int                       curLev,
                                      std::map<int, ID>       & seenAtomIDs,
                                      std::map<int, ID>       & atomIDsInPath)
    {
       
        if ( curLev <iMax )
        {  
            int NachbarpunkteDetected = 0;
            
            // Check Nachbarpunkte
            for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
            {
                int tSeriNum = allAtoms[*tNBA].seriNum;
                
                // Find  Nachbarpunkte of the current path if the current atom
                // (1) should not be the the atom beginning the path. (it is a ring)
                // (2) should not be the one just walked through in the last step
                // (3) is in a list of atoms we have seen
                if (tSeriNum != iOriAto->seriNum && tSeriNum != SeriNumPreAto 
                        && seenAtomIDs.find(tSeriNum) !=seenAtomIDs.end())
                {
                    
                    NachbarpunkteDetected = 1;
                    
                    // found Nachbarpunkte, stop this path
                    /*
                    
                    std::cout << "atom : " <<  allAtoms[*tNBA].id << std::endl;
                    
                    std::cout << "Nachbarpunkte found, stop this path " << std::endl;
                    for (std::map<int, ID>::iterator iS=seenAtomIDs.begin();
                            iS != seenAtomIDs.end(); iS++)
                    {
                        std::cout << "atom : " << iS->first 
                                << " : " << iS ->second << std::endl;
                    }
                    */              
                }  
            }
            
            // check if a ring is formed 
            if (!NachbarpunkteDetected)
            {
                for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
                {
                    int tSeriNum = allAtoms[*tNBA].seriNum;
                    
                    if (tSeriNum == iOriAto->seriNum && tSeriNum != SeriNumPreAto
                        && curLev > 2)
                    {
                        //std::cout << iOriAto->id << " : "
                        //          << curAto.id << " : " << allAtoms[*tNBA].id
                        //          << " : " << SeriNumPreAto 
                        //          << " Find a ring." << std::endl;    
                        //   sort the atoms in the seenAtom vector and 
                        //   check if this ring is already found (the same 
                        //   ring of the same atom should be found twice at
                        //   least because of the walking algorithm.
                        // FIND A RING !
                        atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                        std::list<std::string> tAllIds;
                        std::vector <AtomDict> ttAtoms;
                        for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                        {
                            //if (iSee->first != iOriAto->seriNum)
                            //{
                            tAllIds.push_back(iSee->second);
                            int posIdx = atomPosition(iSee->second);
                            ttAtoms.push_back(allAtoms[posIdx]);
                                // tRing.atoms.push_back(*tAt);
                                //std::cout << iSee->second << std::endl;   
                            //}
                        }
                        RingDict aRingDict(ttAtoms);                 
                        
                        tAllIds.sort(compareNoCase);
                        //std::string tRepStr(iOriAto->id);
                        std::string tRepStr;
                        for (std::list<std::string>::iterator iAI =tAllIds.begin();
                                    iAI != tAllIds.end(); iAI++)
                        {
                            tRepStr.append(*iAI);
                        }
                            
                        iOriAto->ringRep[tRepStr] = (int)atomIDsInPath.size();
                        
                        aRingDict.rep = tRepStr;
                        std::map<ID, std::vector<RingDict> >::iterator iFindRing=allRings.find(tRepStr);
                        if (iFindRing == allRings.end())
                        {
                            allRings[tRepStr].push_back(aRingDict);
                        }
                        
                        atomIDsInPath.erase(curAto.seriNum);  
                        NachbarpunkteDetected = 1;
                        
                    }
                    
                }
            }
              
            if (! NachbarpunkteDetected)
            {
                // found no Nachbarpunkte and no ring
                // descend the new path in the neighborhood graph:
                
                
                int tNewLev = curLev + 1;
                seenAtomIDs.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                if (tNewLev < iMax)
                {
                   /*
                        std::cout << "atom " << curAto.id 
                                  << " finds no Nachbarpunkte in it neighbor  " 
                                  << std::endl << "Descend into the next atom "
                                  << std::endl;
                    
                    */
                    
                    atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                    for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                         tNBA != curAto.connAtoms.end(); tNBA++)
                    {
                        if(curLev==1)
                        {
                            // tempo list of atoms in a path
                            if((int)seenAtomIDs.size() !=0)
                            {
                                seenAtomIDs.clear();
                            }
                            if((int)atomIDsInPath.size() !=0)
                            {
                                atomIDsInPath.clear();
                            }
                            //std::cout << "after clear, the size is " 
                            //          << (int)seenAtomIDs.size() << std::endl;
                            seenAtomIDs.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                            atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                        }
                        if (SeriNumPreAto != allAtoms[*tNBA].seriNum)
                        {
                            /*
                                std::cout << std::endl << "Current size " << curLev << std::endl;
                                std::cout << "Orig atom : " << iOriAto->id
                                          << " Prev atom : " << SeriNumPreAto
                                          << " Curent atom :  " << curAto.id 
                                          << std::endl << std::endl;
                                std::cout << "NB atom : " << allAtoms[*tNBA].id << std::endl;  
                            */
                            // This is a new atom and append the atom to seen-atom list
                            // and call function checkOnePathSec() recursively
                            int tPreSeriNum = curAto.seriNum; 
                            checkOnePathSec(allAtoms[*tNBA], iMax, iOriAto, tPreSeriNum, 
                                            tNewLev, seenAtomIDs, atomIDsInPath);
                        }
                    }
                    atomIDsInPath.erase(curAto.seriNum);
                    seenAtomIDs.erase(curAto.seriNum);
                }
                atomIDsInPath.erase(curAto.seriNum);
                seenAtomIDs.erase(curAto.seriNum);
            }
        }
        else
        {
            atomIDsInPath.erase(curAto.seriNum);
            seenAtomIDs.erase(curAto.seriNum);
        }
    }
    
    
    void AllSystem::outRingSec(AtomDict &tAtom)
    {
        int numRings = (int)tAtom.ringRep.size();
       
        if (numRings)
        {
            
            std::map<int, int> sizeMap;
            
            
            for (std::map<std::string, int>::iterator iMR=tAtom.ringRep.begin();
                    iMR != tAtom.ringRep.end(); iMR++)
            {
                if (sizeMap.find(iMR->second) ==sizeMap.end())
                {
                    sizeMap[iMR->second] =1;
                }
                else
                {
                    sizeMap[iMR->second]++;
                }   
            }
            
            tAtom.codClass.append("[");
            int i =0;
            int j = (int)sizeMap.size();
            for (std::map<int, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                std::string tSize = IntToStr(iSMa->first);
                std::string tNum  = IntToStr(iSMa->second);
                
                if(iSMa->second >= 3)
                {
                    tAtom.codClass.append(tNum + "x" + tSize);
                }
                else if (iSMa->second==2)
                {
                    tAtom.codClass.append( tSize + "," + tSize);
                }   
                else if (iSMa->second==1)
                {
                    tAtom.codClass.append(tSize);
                }
                       
                
                if(i != j-1)
                {
                    tAtom.codClass.append(",");
                }
                else
                {
                    tAtom.codClass.append("]");
                }
                    
                i++;
            }
        }
    }
    
    std::string  AllSystem::outRingSecStr(AtomDict &tAtom)
    {
        std::string tS1 = "";
        int numRings = (int)tAtom.ringRep.size();
        
        if (numRings)
        {
            
            std::map<int, int> sizeMap;
            
            
            for (std::map<std::string, int>::iterator iMR=tAtom.ringRep.begin();
                    iMR != tAtom.ringRep.end(); iMR++)
            {
                if (sizeMap.find(iMR->second) ==sizeMap.end())
                {
                    sizeMap[iMR->second] =1;
                }
                else
                {
                    sizeMap[iMR->second]++;
                }   
            }
            
            tS1.append("[");
            int i =0;
            int j = (int)sizeMap.size();
            for (std::map<int, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                std::string tSize = IntToStr(iSMa->first);
                std::string tNum  = IntToStr(iSMa->second);
               
                if(iSMa->second !=1)
                {
                    tS1.append(tNum + "x" + tSize);
                }
                else
                {
                    tS1.append(tSize);
                }
                
                if(i != j-1)
                {
                    tS1.append(",");
                }
                else
                {
                    tS1.append("]");
                }
                    
                i++;
            }
        }
        
        return tS1;
        
    }
    
    void AllSystem::chiralExch()
    {
        
    }
    
    
    // Plane related 
    void AllSystem::detectPlaneGroups()
    {
        
        groupOrgAtomsToPlanes();
        //groupMetAndLigandAtomsToPlanes();
   
        //Check
        std::cout<< "There are " << (int)allPlanes.size() 
                <<" Planes in the system" << std::endl;
        for (int i=0; i < (int)allPlanes.size(); i++)
        {
            std::cout<<"Plane " <<i+1 << " contains "
                    << (int)allPlanes[i].atoms.size() 
                    << " atoms. They are: "<< std::endl;
            for (std::map<ID, int>::iterator iAt=allPlanes[i].atoms.begin();
                    iAt!=allPlanes[i].atoms.end(); iAt++)
            {
                std::cout<<iAt->first << "\t";
            }
            std::cout<<std::endl;
        }
        
    }
    
    void AllSystem::groupOrgAtomsToPlanes()
    {
        std::vector<PlaneDict> smalPls;
        
        // Find the smallest planes
        setSmallestPLs(smalPls);
        
        // Merge the above planes to the large plane groups
        mergeLargePLGroups(smalPls);
        
    }
    
    void AllSystem::setSmallestPLs(std::vector<PlaneDict>& tSmaPls)
    {
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            if (allAtoms[i].bondingIdx == 2)
            {
                PlaneDict tSmaPl;
                ID        aID     = allAtoms[i].id;
                tSmaPl.archID  = aID;
                tSmaPl.archPos = i;
                tSmaPl.atoms[aID] = i;
                
                for (int j=0; j< (int)allAtoms[i].connAtoms.size(); j++)
                {
                    int tPos  = allAtoms[i].connAtoms[j];
                    if (!allAtoms[tPos].isMetal)
                    {
                        ID  tId   = allAtoms[tPos].id;
                        tSmaPl.atoms[tId] = tPos;
                    }
                }
                tSmaPls.push_back(tSmaPl);
            }
        }
        /*
        std::cout <<"There are " << (int)tSmaPls.size() 
                  << " sp2 planes " << std::endl;
        for (std::vector<PlaneDict>::iterator iP=tSmaPls.begin();
                iP !=tSmaPls.end(); iP++)
        {
            std::cout <<"small plane: " << std::endl;
            std::cout << "Center atom " << iP->archID << std::endl;
            for (std::map<ID, int>::iterator iA=iP->atoms.begin();
                    iA !=iP->atoms.end(); iA++)
            {
                if (iA->first != iP->archID )
                {
                    std::cout << "Other atom " << iA->first << std::endl;
                }
            }
        }
        
        */
    }
    
    void AllSystem::mergeLargePLGroups(std::vector<PlaneDict>& tSmaPls)
    {
        //std::cout <<"There are " << (int)tSmaPls.size() 
        //          << " sp2 planes " << std::endl;
        std::map<int, std::vector<int> > allPlsClasses;
        
        // assign a e-id to each plane set
        // allPlsClasses[e-classId] = vector of small plane indexes 
        for (int i=0; i < (int)tSmaPls.size(); i++)
        {
            for (std::map<ID,int>::iterator iSA=tSmaPls[i].atoms.begin();
                   iSA!=tSmaPls[i].atoms.end(); iSA++ )
            {
                allPlsClasses[i].push_back(iSA->second);
            }
        }
             
        // get equiv-class by ring-relation
        
        for (int i=0; i < (int)tSmaPls.size(); i++)
        {
            bool iMerge = false;
            /*
            std::cout << "small plane " << i << std::endl;
            std::cout << "Archor atoms: " <<tSmaPls[i].archID << std::endl;
            for (std::map<ID, int>::iterator iSA1=tSmaPls[i].atoms.begin();
                   iSA1 != tSmaPls[i].atoms.end(); iSA1++)
            {
                std::cout << iSA1->first << "\t";
            }
            std::cout<<std::endl;
            */
            
            for (int j=i+1; j <(int)tSmaPls.size(); j++)
            {
                /*
                std::cout << "linked small plane " << j << std::endl;
                std::cout << "Archor atoms: " << tSmaPls[j].archID << std::endl;
                for (std::map<ID, int>::iterator iSA2=tSmaPls[j].atoms.begin();
                   iSA2 != tSmaPls[j].atoms.end(); iSA2++)
                {
                    std::cout << iSA2->first << "\t";
                }
                std::cout<<std::endl;
                
                 */
                
                std::map<int, std::vector<int> >::iterator tFindM;
                tFindM = allPlsClasses.find(i);
                
                if (tFindM != allPlsClasses.end())
                {
                    //std::cout << "test ring-related " << std::endl;
                    //std::cout << isInSameRing(tSmaPls[i],tSmaPls[j]) << std::endl;
                   
                    if(isInSameRing(tSmaPls[i],tSmaPls[j]))
                    {   
                        // std::cout << "Merged planes " << i <<"  "<< j << std::endl;
                        for (std::vector<int>::iterator iI = allPlsClasses[i].begin();
                                iI !=allPlsClasses[i].end(); iI++ )
                        {
                            std::vector<int>::iterator tFindV;
                            tFindV = std::find(allPlsClasses[j].begin(), 
                            allPlsClasses[j].end(), *iI);
                            if (tFindV ==allPlsClasses[j].end())
                            {
                                allPlsClasses[j].push_back(*iI); 
                            }
                        }
                        //std::cout << "Plane " << j << "now has atoms " 
                        //          << (int)allPlsClasses[j].size() << std::endl;
                        iMerge = true;        
                    }    
                }
            }
            
            if (iMerge)
            {
                allPlsClasses.erase(i);
            }
        }  
        
        // Further
        
        std::vector<int> delKeys;
        for(std::map<int, std::vector<int> >::iterator iPl1=allPlsClasses.begin();
                iPl1!=allPlsClasses.end(); iPl1++)
        {
            for(std::map<int, std::vector<int> >::iterator iPl2=allPlsClasses.begin();
                iPl2!=allPlsClasses.end(); iPl2++)
            {
                if (iPl1->first < iPl2->first)
                {   
                    if(furtherM(iPl1->second, iPl2->second))
                    {    
                        
                        if((int)iPl1->second.size() >= (int)iPl2->second.size())
                        {
                            for (std::vector<int>::iterator iV=iPl2->second.begin();
                                    iV!=iPl2->second.end(); iV++)
                            {
                                
                                std::vector<int>::iterator tFindV;
                                tFindV = std::find(iPl1->second.begin(), 
                                         iPl1->second.end(), *iV);
                                if (tFindV ==iPl1->second.end())
                                {
                                    iPl1->second.push_back(*iV); 
                                }
                            }
                            delKeys.push_back(iPl2->first);   
                        }
                        else
                        {
                            for (std::vector<int>::iterator iV=iPl1->second.begin();
                                    iV!=iPl1->second.end(); iV++)
                            {
                                std::vector<int>::iterator tFindV;
                                tFindV = std::find(iPl2->second.begin(), 
                                         iPl2->second.end(), *iV);
                                if (tFindV ==iPl2->second.end())
                                {
                                    iPl2->second.push_back(*iV); 
                                }
                            }
                            delKeys.push_back(iPl1->first);
                        }
                    }
                }
            }
        }
        
        for(std::vector<int>::iterator iDel=delKeys.begin();
                iDel!=delKeys.end(); iDel++)
        {
            if(allPlsClasses.find(*iDel) != allPlsClasses.end())
            {
                allPlsClasses.erase(*iDel);
            }
        }
       
        
        // Put equiv-classes in a vector of planeDicts 
        for (std::map<int, std::vector<int> >::iterator iPC=allPlsClasses.begin();
                iPC!=allPlsClasses.end(); iPC++ )
        {
            if (iPC->second.size() >3 )
            {
                PlaneDict aPl;
            
                for (std::vector<int>::iterator iId=iPC->second.begin();
                       iId !=iPC->second.end(); iId++ )
                {
                    aPl.atoms[allAtoms[*iId].id]=*iId;
                }
                allPlanes.push_back(aPl);
            }
        } 
    }

    bool AllSystem::isInSameRing(PlaneDict & tP1, PlaneDict & tP2)
    {
        if (std::find(allAtoms[tP1.archPos].connAtoms.begin(),
                allAtoms[tP1.archPos].connAtoms.end(), tP2.archPos)
                ==allAtoms[tP1.archPos].connAtoms.end())
        {
            // std::cout << "Not in " << tP1.archID << "Neighbor " << std::endl; 
            return false;
        } 
        
        if ((int)allAtoms[tP1.archPos].ringRep.size())
        {
                for (std::map<std::string, int>::iterator iM1=allAtoms[tP1.archPos].ringRep.begin();
                        iM1 !=allAtoms[tP1.archPos].ringRep.end(); iM1++)
                {
                    // std::cout<< "1: " << iM1->first << std::endl;
                    if ((int)allAtoms[tP2.archPos].ringRep.size())
                    {
                       /* 
                        for (std::map<std::string, int>::iterator iM2=allAtoms[tP2.archPos].ringRep.begin();
                                 iM2 !=allAtoms[tP2.archPos].ringRep.end(); iM2++)
                        {
                            std::cout<< "2: " << iM2->first << std::endl;
                        }
                        */
                            
                        if(allAtoms[tP2.archPos].ringRep.find(iM1->first) !=
                                   allAtoms[tP2.archPos].ringRep.end())
                        {
                            if (allRings[iM1->first][0].isPlanar)
                            {
                                return true;
                            }
                        }
                    }
                }
            
        }
        
        return false;   
    }
    
    bool AllSystem::furtherM(std::vector<int> &tV1, std::vector<int> &tV2)
    {
        int nFind = 0;
        for(std::vector<int>::iterator iV=tV1.begin();
                iV !=tV1.end(); iV++)
        {
           std::vector<int>::iterator tFindV;
           tFindV = std::find(tV2.begin(), tV2.end(), *iV);
           if (tFindV !=tV2.end())
           {
               nFind++; 
           } 
        }
        
        if (nFind >=3)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
        
    
    
    void AllSystem::setupAllTargetValuesFromCOD(ID tOutName, ID tMonoName)
    {
       CodClassify  aCodSystem(allAtoms, allHAtomIdx, allBonds, allAngles, 
                               allTorsions, allChirals, allPlanes, allRings);
        
     
       aCodSystem.setupAllTargetValues();
       
      
      /* 
       outMMCif(tOutName.c_str(),
                tMonoName, 
                aCodSystem.allAtoms,
                aCodSystem.allHAtomIdx,
                aCodSystem.allBonds,
                aCodSystem.allAngles,
                aCodSystem.allTorsions,
                aCodSystem.allRings,
                aCodSystem.allPlanes,
                aCodSystem.allChirals);
       */
       resetSystem(aCodSystem);
    
       /*
       for(std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            std::cout << "Atom " << iAt->id << " is in " 
                      << (int)iAt->ringRep.size() << std::endl;
        }
        
        */ 
       
    }
    
    void AllSystem::OutputRestraintCif(FileName tCifName)
    {
        /*
        std::ofstream outCif(tCifName);
        
        if (outCif.is_open())
        {
            OutputRestraintCifHead(outCif);
            
            if ((int)allBonds.size())
            {
                OutputRestraintCifBonds(outCif);
            }
            
            if ((int)allAngles.size())
            {
                OutputRestraintCifAngles(outCif);
            }
            
            if ((int)allTorsions.size())
            {
                OutputRestraintCifTosions(outCif);
            }
            
            if((int)allPlanes.size())
            {
                OutputRestraintCifPlanes(outCif);
            }
            
            if((int)allChirals.size())
            {
                OutputRestraintCifChirals(outCif);
            }
            
            if ((int)allRings.size())
            {
                OutputRestraintCifRings(outCif);
            }
        }
         */
        
    }
    /*
    void AllSystem::OutputRestraintCifHead(std::ofstream& tCifFile)
    {
        // Tempo 
        tCifFile << "global_  " << std::endl;
        tCifFile << "_entry.id  XXXXXX" << std::endl;
        tCifFile << "_struct.keywords   '----' " << std::endl;
        tCifFile << "_autid.creation_date   XX-XXX-XX" << std::endl;
        tCifFile << "_struct.title  coordinates of XXXXX  from program: XXXXXX" << std::endl;
        tCifFile << "_lib.name     mon_lib   " << std::endl;
        tCifFile << "_lib.version  XX.XX.XX  " << std::endl;
        tCifFile << "_lib.update   XX/XX/XX  " << std::endl;
        tCifFile << "# ------------------------------------------------" << std::endl;
        tCifFile << "data_restraints " << std::endl << "# " << std::endl;
        
        tCifFile << "#" << std::endl << "loop_" << std::endl <<"#" << std::endl;
        tCifFile << "_restr.record"  << std::endl;
        tCifFile << "_restr.number"  << std::endl;
        tCifFile << "_restr.label"   << std::endl;
        tCifFile << "_restr.period"  << std::endl;
        tCifFile << "_restr.atom_id_1"  << std::endl;
        tCifFile << "_restr.atom_id_2"  << std::endl;
        tCifFile << "_restr.atom_id_3"  << std::endl;
        tCifFile << "_restr.atom_id_4"  << std::endl;
        tCifFile << "_restr.value"      << std::endl;
        tCifFile << "_restr.dev"        << std::endl;
        tCifFile << "_restr.val_obs"    << std::endl;
        tCifFile << "_restr.dist"       << std::endl;
        tCifFile << "_restr.dist_dev"   << std::endl;
        tCifFile << "_restr.econst"     << std::endl;
        
    }
    
    void AllSystem::OutputRestraintCifBonds(std::ofstream  & tCifFile)
    {

        
        int i =1;
        tCifFile.precision(3);
        
        for (std::vector<BondDict>::iterator iB=allBonds.begin();
                iB !=allBonds.end(); iB++)
        {
            int j = iB->fullAtoms[iB->atoms[0]];
            int k = iB->fullAtoms[iB->atoms[1]];
            
            tCifFile  << std::fixed << "BOND"  << "\t" << i  <<  "\t"
                      << iB->order  << "\t" << "."  << "\t" << j  
                      <<  "\t" << k << "\t" << "."  << "\t"
                      <<  "."  << iB->value   << "\t" << iB->sigValue
                      <<  "."  << "\t" << "."  << "\t" << "." << "."
                      <<  "#    " << iB->atoms[0] << "   " << iB->atoms[1]
                      << std::endl;
            i++;
        }
        
    }
    
    void AllSystem::OutputRestraintCifAngles(std::ofstream & tCifFile)
    {
        int i =1;
        tCifFile.precision(3);
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA!=allAngles.end(); iA++)
        {
            tCifFile  << std::fixed << "ANGL"  << "\t" << i  <<  "\t"
                      << "."  << "\t" << "."   << "\t" << iA->atoms[0] 
                      <<  "\t" << iA->atoms[1] << "\t" << iA->atoms[2]  
                      << "\t"  <<  "."  << iA->value   << "\t" << iA->sigValue
                      <<  "."  << "\t" << "."  << "\t" << "."  << "\t" << "."
                      <<  "#    " << allAtoms[iA->atoms[0]].id 
                      << "   " << allAtoms[iA->atoms[1]].id
                      << "   " << allAtoms[iA->atoms[2]].id << std::endl;
            i++;
        }
        
    }
    
    void AllSystem::OutputRestraintCifTosions(std::ofstream & tCifFile)
    {
        int i =1;
        tCifFile.precision(3);
        
        for (std::vector<TorsionDict>::iterator iT=allTorsions.begin();
                iT != allTorsions.end(); iT++)
        {
            i++;
        }
    }
    
    void AllSystem::OutputRestraintCifChirals(std::ofstream & tCifFile)
    {
        int i =1;
        tCifFile.precision(3);
        
        for (std::vector<ChiralDict>::iterator iCh = allChirals.begin();
                iCh != allChirals.end(); iCh++)
        {
            
            tCifFile  << std::fixed << "CHIR"  << "\t" << i  <<  "\t"
                      << iCh->value << "\t" << "." << "\t" << iCh->atoms[0]  
                      << "\t" << iCh->atoms[1]  << "\t" << iCh->atoms[2]
                      << iCh->value   << "\t" << "."
                      <<  "."  << "\t" << "."  << "\t" << "."  << "\t" << "."
                      <<  "#    " << allAtoms[iCh->atoms[0]].id 
                      << "   " << allAtoms[iCh->atoms[1]].id
                      << "   " << allAtoms[iCh->atoms[2]].id 
                      << std::endl;
            i++;
        }
        
        
    }
        
    void AllSystem::OutputRestraintCifPlanes(std::ofstream & tCifFile)
    {
        int i =1;
        tCifFile.precision(3);
        
        for (std::vector<PlaneDict>::iterator iP=allPlanes.begin();
                iP != allPlanes.end(); iP++)
        {
            for (std::map<ID, int>::iterator iA = iP->atoms.begin();
                    iA != iP->atoms.end(); iA++)
            {
                tCifFile  << std::fixed << "PLAN"  << "\t" << i  <<  "\t"
                          << "." << "\t" << "." << "\t" << iA->second  
                          << "\t" << "."  << "\t" << "." << "\t"  
                          << "\t" << "0.020"   << "\t" << "." << "\t" 
                          <<  "."  << "\t" << "."  << "\t" << "."  << "\t" << "."
                          <<  "#    " << iA->first << std::endl;
            }
        }
        
    }
    
    void AllSystem::OutputRestraintCifRings(std::ofstream  & tCifFile)
    {   
    }
     */
}
