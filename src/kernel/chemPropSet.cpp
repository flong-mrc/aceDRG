
/* 
 * File:   chemPropSet.cpp
 * Author: flong
 *
 * Created on July 3, 2014, 1:39 PM
 */

#include "chemPropSet.h"

namespace LIBMOL
{
    
    extern bool assignElementType(PeriodicTable & tP, std::string tStr, 
                                  std::vector<AtomDict>::iterator tAtom)
    {
        bool tFind = false;
        if (tP.elements.find(tStr) !=tP.elements.end())
        {
            tAtom->chemType = tStr;
            tFind = true;
        }
        else
        {
            if (tStr.size() > 1)
            {
                std::string tSubStr = tStr.substr(0, tStr.size()-1);
                if(tP.elements.find(tSubStr) != tP.elements.end())
                {
                    tAtom->chemType = tSubStr;
                    tFind = true;
                }
            }
        }
        
        return tFind;
    }

    
    extern int getNumOxyConnect(std::vector<AtomDict>  &  tAtoms,
                                std::vector<AtomDict>::iterator iA)
    {
        int nO=0;
        for (std::vector<int>::iterator iC=iA->connAtoms.begin();
                iC !=iA->connAtoms.end(); iC++)
        {
            if (tAtoms[*iC].chemType.compare("O")==0)
            {
                nO++;
            }
        }
        return nO;
    }
    
    extern void getHydroAtomConnect(std::vector<AtomDict>  &  tAtoms)
    {
        for(std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA!=tAtoms.end(); iA++)
        {
            if(iA->chemType.compare("H") !=0)
            {
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                    iNB !=iA->connAtoms.end(); iNB++)
                {
                    if (tAtoms[*iNB].chemType.compare("H")==0)
                    {
                        iA->connHAtoms.push_back(*iNB);
                    }
                }
            }
        }
        
        // Check
        /*
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA!=allAtoms.end(); iA++)
        {
            std::cout << "Atom " << iA->id << " connects to " 
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
    
    // supplement function 
    extern void mdChiralByClasses(std::vector<AtomDict>::iterator tAt,
                                  std::vector<AtomDict>        &  tAtoms)
    {
        std::vector<ID> atps;
                        
        for (std::vector<int>::iterator iNA=tAt->connAtoms.begin();
                                iNA !=tAt->connAtoms.end(); iNA++)
        {
            if (std::find(atps.begin(), atps.end(), tAtoms[*iNA].codClass)==atps.end())
            {
                atps.push_back(tAtoms[*iNA].chemType);
            }
        }
        
        if ((int)atps.size() >2)
        {
            tAt->chiralIdx  = 2;
        }
        else
        {
            tAt->chiralIdx =0;
        }
    }
    
    // Set atom's bonding features (sp, sp2, sp3 and chiral center) based on 
    // the atom's connections.
   
    extern void setAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms)
    {
  
        
        // First round
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            /*
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!tAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }
            */
            int t_len = (int)iAt->connAtoms.size();
            //std::cout << "Atom " << iAt->id << std::endl
            //        <<  " connect to  " << t_len << std::endl;
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
                    if (getNumOxyConnect(tAtoms, iAt)==1)
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
            else if (iAt->chemType.compare("N")==0 
                    || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    iAt->chiralIdx  = 2;
                    iAt->bondingIdx = 3;  
                }
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("O")==0)
            {
                if ((int)iAt->connAtoms.size()==2)
                {
                    iAt->bondingIdx = 3;
                }
                else if (iAt->connAtoms.size()==1)
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
                        mdChiralByClasses(iAt, tAtoms);
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
            else if (iAt->chemType.compare("S")==0 || iAt->chemType.compare("SE")==0 )
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3 || t_len==2)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3; 
                }
                /*
                else if (t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx =4; // --> spd 
                }
                 */
            }
            // std::cout << "its chiralIdx " << iAt->chiralIdx << std::endl;
        }
        
        // more conditions 
        // Do oxygen atom first, to see if an Oxygen atom of two connections 
        // is sp2. The default one in the above step is sp3 
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!tAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("O")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                
                if(t_len==2)
                {
                    if (iAt->parCharge ==0.0)
                    {
                        int nH = 0;
                        for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                                 iCA != iAt->connAtoms.end(); iCA++)
                        {
                            if(tAtoms[*iCA].chemType.compare("H")==0)
                            {
                                nH++;
                            }
                        }
                        
                        if (nH !=1)
                        {
                            bool l_sp2 = false;
                            for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                                 iCA != iAt->connAtoms.end(); iCA++)
                            {
                                if(tAtoms[*iCA].bondingIdx == 2)
                                {
                                    l_sp2 = true;
                                    break;
                                }
                            }
                            
                            if (l_sp2)
                            {
                                // Now we can say this atom is in sp2 orbits 
                                iAt->chiralIdx  =  0;
                                iAt->bondingIdx =  2;
                            }
                        }
                    }
                } 
            }
            //std::cout << "Again atom " << iAt->id << " its chiralIdx " 
            //          << iAt->chiralIdx << std::endl;
        }
        
        // Then N and B atoms
        std::map<int, int> preBondingIdx;
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            preBondingIdx[iAt->seriNum] = iAt->bondingIdx;
        }
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!tAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }
            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                
                if(t_len==3)
                {
                    if (iAt->parCharge ==0.0)
                    {
                        
                        bool l_sp2 = false;
                        for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                                 iCA != iAt->connAtoms.end(); iCA++)
                        {
                            //if(tAtoms[*iCA].bondingIdx == 2)
                            if (preBondingIdx[*iCA]==2)
                            {
                                l_sp2 = true;
                                break;
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
                    else if (iAt->parCharge ==1.0)
                    {
                        iAt->chiralIdx  =  0;
                        iAt->bondingIdx =  2;
                    }
                } 
            }
            //std::cout << "Again atom " << iAt->id << " its chiralIdx " 
            //          << iAt->chiralIdx << std::endl;
        }
        
        // Further check if a chiral center is a real one
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = tAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(tAtoms[*iNB].id);
                    }
                }
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
            
            
        }
  
        // Check
        /*
        std::cout << "Chiral and plane feather for atoms in the system" 
                  << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            if (iAt->chiralIdx == -1)
            {
                std::cout << "Atom " << iAt->id << " may be  in a negative chiral center " 
                        << std::endl;
            }
            else if (iAt->chiralIdx == 1)
            {
                std::cout << "Atom " << iAt->id << " may be in a positive chiral center "
                        << std::endl;
            }
            else if (iAt->chiralIdx==2)
            {
                std::cout << "Atom " << iAt->id 
                        << " may be in a chiral center but the volume sign undefined "
                        << std::endl;
            }
            else if (iAt->chiralIdx==0)
            {
                std::cout << "Atom " << iAt->id 
                        << " is not a chiral center" << std::endl;
            }
        }
         */
    }
    
    extern void modAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds, 
                                               std::vector<AngleDict> & tAngles,
                                               std::vector<RingDict>  & tRings,
                                               int                      tMode)
    {
        REAL angCri = 10.0;
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            if((iA->chemType.compare("N")==0 || iA->chemType.compare("B")==0) 
                  && (iA->connAtoms.size() == 3))
            {
                std::cout << "Check " << iA->id << " now " << std::endl;
                std::cout << "Its initial sp is " << iA->bondingIdx << std::endl;
                bool lAromRs = false;
                for (std::vector<int>::iterator iCo=iA->connAtoms.begin();
                        iCo !=iA->connAtoms.end(); iCo++)
                {
                    if (tAtoms[*iCo].inRings.size() !=0)
                    {  
                        //std::cout << "connected atom " << tAtoms[*iCo].id 
                        //          << "is in rings: " << std::endl;
                        
                        for (std::vector<int>::iterator iR=tAtoms[*iCo].inRings.begin();
                               iR !=tAtoms[*iCo].inRings.end(); iR++)
                        {
                            // std::cout << tRings[*iR].rep << std::endl;
                            if(tRings[*iR].isAromatic || tAtoms[*iCo].bondingIdx==2)
                            {
                                // std::cout << "It is sp2 related " << std::endl;
                                lAromRs=true;
                                break;
                            }
                        }
                    }
                    //else if (tAtoms[*iCo].bondingIdx==2)
                    //{
                        //std::cout << "Check NB atom " << tAtoms[*iCo].id << std::endl;
                        //std::cout << "It is sp2 related " << std::endl;
                        //lAromRs=true;
                        //break;
                    //}
                }
                
                // Use nAromRs, =1 will make decision temporarily 
                // adjust in future. 
                // std::cout << "lAromRs " << lAromRs << std::endl;
                if (lAromRs)
                {
                    
                    //if (checkBridgeStruct(tAtoms, tRings, iA->seriNum))
                    //{
                        // std::cout << "Inside 1" << std::endl;
                    //    iA->chiralIdx  = 5;
                    //    iA->bondingIdx = 3;
                    //}
                    if (tMode==1)
                    {
                        if (iA->isInPreCell)
                        {
                           // if (confirmPlaneByChiralVol(tAtoms, iA))
                           if (confirmPlaneByAngle(tAtoms, iA, angCri)) 
                           {
                               iA->chiralIdx  = 5; // New value 
                               iA->bondingIdx = 2; // still keep it sp2, will use together
                                            // with chiralIdx
                               //std::cout << "inside 2 " << std::endl;
                           }
                           else
                           {
                               iA->chiralIdx  = 5;
                               iA->bondingIdx = 3;
                           }
                        }      
                    }
                    else
                    {
                       if (confirmPlaneByAngle(tAtoms, iA, angCri))
                       {
                            iA->chiralIdx  = 5; // New value 
                            iA->bondingIdx = 2; // still keep it sp2, will use together
                                            // with chiralIdx
                               //std::cout << "inside 2 " << std::endl;
                        }
                        else
                        {
                            iA->chiralIdx  = 5;
                            iA->bondingIdx = 3;
                        }
                    }
                }
                else
                {
                    //std::cout << "workMode : " << tMode << std::endl;
                    if (tMode==1)
                    {
                        if (iA->isInPreCell)
                        {
                            if (confirmPlaneByAngle(tAtoms, iA, angCri))
                            {
                                iA->chiralIdx  = 5; // New value 
                                iA->bondingIdx = 2; // still keep it sp2, will use together
                                // with chiralIdx
                                //std::cout << "inside 2 " << std::endl;
                            }
                            else
                            {
                                iA->chiralIdx  = 5;
                                iA->bondingIdx = 3;
                            }
                        }
                    }
                    else
                    {
                        if (confirmPlaneByAngle(tAtoms, iA, angCri))
                        {
                            iA->chiralIdx  = 5; // New value 
                            iA->bondingIdx = 2; // still keep it sp2, will use together
                                // with chiralIdx
                                //std::cout << "inside 2 " << std::endl;
                        }
                        else
                        {
                            iA->chiralIdx  = 5;
                            iA->bondingIdx = 3;
                        }
                    }
                }
                /*
                 * else will keep its initial setting
                else
                {
                    iA->chiralIdx  = 2;
                    iA->bondingIdx = 3;
                }
                 */
                std::cout << "Here Its hybridization is sp" << iA->bondingIdx << std::endl; 
            }
            
            iA->hybrid = strTransSP(iA->bondingIdx);
        }
        
        // Remember atom hybridizations in bonds and Angles 
        // ( will be used for DB analysis)  
        for (std::vector<BondDict>::iterator iB=tBonds.begin();
                iB !=tBonds.end(); iB++)
        {
            iB->atomSPs.clear();
            for (std::vector<int>::iterator iA=iB->atomsIdx.begin();
                    iA !=iB->atomsIdx.end(); iA++)
            {
                //std::cout << "Atom " << tAtoms[*iA].id  
                //          << "\t" << " hybrid " 
                //          << tAtoms[*iA].hybrid << std::endl;
                iB->atomSPs[tAtoms[*iA].id] = tAtoms[*iA].hybrid;
            }
        }   
    }
    
    extern void setAtomsNB1NB2_SP(std::vector<AtomDict> & tAtoms)
    {
        // make sure codAtmRoot are there
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            std::vector<ID>  cIDs;
            StrTokenize(iAt->codClass, cIDs, '(');
            if (cIDs.size() >0)
            {
                iAt->codAtmRoot = TrimSpaces(cIDs[0]);
            }
            else
            {
                std::cout << "can not find root symbol for atom class "
                          << iAt->codClass << std::endl;
                exit(1);
            }
        }
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            int aAt_seri = iAt->seriNum;
            std::vector<ID> aNB1_NB2SP_Set;
            for (std::vector<int>::iterator iNB1=iAt->connAtoms.begin();
                    iNB1 != iAt->connAtoms.end(); iNB1++)
            {
                ID aNB1_main = tAtoms[*iNB1].codAtmRoot;
                ID aNB2SpStr="";
                std::vector<int> aNB2SpSet;
                for (std::vector<int>::iterator iNB2=tAtoms[*iNB1].connAtoms.begin();
                        iNB2 != tAtoms[*iNB1].connAtoms.end(); iNB2++)
                {   
                    aNB2SpSet.push_back(tAtoms[*iNB2].bondingIdx);
                }
                
                std::sort(aNB2SpSet.begin(), aNB2SpSet.end(), std::greater<int>());
                for (unsigned i=0; i < aNB2SpSet.size(); i++)
                {
                    aNB2SpStr.append(IntToStr(aNB2SpSet[i]));
                    if (i != aNB2SpSet.size()-1)
                    {
                        aNB2SpStr.append("_");
                    }
                }
                aNB1_NB2SP_Set.push_back(aNB1_main + "-"+aNB2SpStr);
            }
            
            std::sort(aNB1_NB2SP_Set.begin(), aNB1_NB2SP_Set.end(), compareNoCase);
            
            iAt->codNB1NB2_SP.clear();
            for (unsigned i=0; i < aNB1_NB2SP_Set.size(); i++)
            {
                iAt->codNB1NB2_SP.append(aNB1_NB2SP_Set[i]);
                if (i != aNB1_NB2SP_Set.size()-1)
                {
                    iAt->codNB1NB2_SP.append(":");
                }
            }  
        }   
    }
    
    extern void setAtomsNB1NB2_exElectrons(std::vector<AtomDict> & tAtoms)
    {
                // make sure codAtmRoot are there
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            std::vector<ID>  cIDs;
            StrTokenize(iAt->codClass, cIDs, '(');
            if (cIDs.size() >0)
            {
                iAt->codAtmRoot = TrimSpaces(cIDs[0]);
            }
            else
            {
                std::cout << "can not find root symbol for atom class "
                          << iAt->codClass << std::endl;
                exit(1);
            }
        }
        
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            int aAt_seri = iAt->seriNum;
            std::vector<ID> aNB1_NB2SP_Set;
            for (std::vector<int>::iterator iNB1=iAt->connAtoms.begin();
                    iNB1 != iAt->connAtoms.end(); iNB1++)
            {
                ID aNB1_main = tAtoms[*iNB1].codAtmRoot;
                ID aNB2SpStr="";
                std::vector<int> aNB2SpSet;
                for (std::vector<int>::iterator iNB2=tAtoms[*iNB1].connAtoms.begin();
                        iNB2 != tAtoms[*iNB1].connAtoms.end(); iNB2++)
                {   
                    aNB2SpSet.push_back(tAtoms[*iNB2].excessElec);
                }
                
                std::sort(aNB2SpSet.begin(), aNB2SpSet.end(), std::greater<int>());
                
                for (unsigned i=0; i < aNB2SpSet.size(); i++)
                {
                    aNB2SpStr.append(IntToStr(aNB2SpSet[i]));
                    if (i != aNB2SpSet.size()-1)
                    {
                        aNB2SpStr.append("_");
                    }
                }
                aNB1_NB2SP_Set.push_back(aNB1_main + "-"+aNB2SpStr);
            }
            
            std::sort(aNB1_NB2SP_Set.begin(), aNB1_NB2SP_Set.end(), compareNoCase);
            
            iAt->codNB1NB2_ExElec.clear();
            for (unsigned i=0; i < aNB1_NB2SP_Set.size(); i++)
            {
                iAt->codNB1NB2_ExElec.append(aNB1_NB2SP_Set[i]);
                if (i != aNB1_NB2SP_Set.size()-1)
                {
                    iAt->codNB1NB2_ExElec.append(":");
                }
            }  
        }
    }
    
    extern void setBondsAndAngles_NB1NB2_SP(std::vector<AtomDict> & tAtoms,
                                            std::vector<BondDict> & tBonds,
                                            std::vector<AngleDict> & tAngles)
    {
        for (std::vector<BondDict>::iterator iB=tBonds.begin();
                iB !=tBonds.end(); iB++)
        {
            iB->atomNB1NB2SPs.clear();
            for (std::vector<int>::iterator iA=iB->atomsIdx.begin();
                    iA !=iB->atomsIdx.end(); iA++)
            {
                //std::cout << "Atom " << tAtoms[*iA].id  
                //          << "\t" << " hybrid " 
                //          << tAtoms[*iA].hybrid << std::endl;
                iB->atomNB1NB2SPs[tAtoms[*iA].id] = tAtoms[*iA].codNB1NB2_SP;
            }
        }   
        
        for (std::vector<AngleDict>::iterator iAn=tAngles.begin();
                iAn !=tAngles.end(); iAn++)
        {
            iAn->atomsNB1NB2SPStats.clear();
            // std::cout << "XXXXX angle " << std::endl;
            for (std::vector<int>::iterator iAt=iAn->atoms.begin();
                    iAt !=iAn->atoms.end(); iAt++)
            {
                //std::cout << "Atom " << tAtoms[*iAt].id  
                //          << "\t" << " hybrid " << tAtoms[*iAt].hybrid << std::endl;
                iAn->atomsNB1NB2SPStats[tAtoms[*iAt].id] =  tAtoms[*iAt].codNB1NB2_SP;
            }
        }
    }
    
    
    extern void setBondsAndAngles_NB1NB2_EE(std::vector<AtomDict> & tAtoms,
                                            std::vector<BondDict> & tBonds,
                                            std::vector<AngleDict> & tAngles)
    {
        for (std::vector<BondDict>::iterator iB=tBonds.begin();
                iB !=tBonds.end(); iB++)
        {
            iB->atomNB2ExtraEls.clear();
            for (std::vector<int>::iterator iA=iB->atomsIdx.begin();
                    iA !=iB->atomsIdx.end(); iA++)
            {
                //std::cout << "Atom " << tAtoms[*iA].id  
                //          << "\t" << " hybrid " 
                //          << tAtoms[*iA].hybrid << std::endl;
                iB->atomNB2ExtraEls[tAtoms[*iA].id] = tAtoms[*iA].codNB1NB2_ExElec;
            }
        }   
        
        /*
        for (std::vector<AngleDict>::iterator iAn=tAngles.begin();
                iAn !=tAngles.end(); iAn++)
        {
            iAn->atomsNB1NB2SPStats.clear();
            // std::cout << "XXXXX angle " << std::endl;
            for (std::vector<int>::iterator iAt=iAn->atoms.begin();
                    iAt !=iAn->atoms.end(); iAt++)
            {
                //std::cout << "Atom " << tAtoms[*iAt].id  
                //          << "\t" << " hybrid " << tAtoms[*iAt].hybrid << std::endl;
                iAn->atomsNB1NB2SPStats[tAtoms[*iAt].id] =  tAtoms[*iAt].codNB1NB2_SP;
            }
        }
         */
    }
    
    
    extern void reIndexAtomInRing(std::vector<AtomDict> & tAtoms,
                                  std::vector<RingDict> & tRings)
    {
                // re-index atom's inRing idx 
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            iAt->inRings.clear();
        }
        
        int idxR = 0;
        for (std::vector<RingDict>::iterator iR=tRings.begin();
                iR!=tRings.end(); iR++)
        {
            std::vector<ID> tAtIds;
            for(std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
            {
                tAtIds.push_back(iAt->id);
            }
            for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
            {
                if(std::find(tAtIds.begin(), tAtIds.end(), iAt->id)
                    !=tAtIds.end())
                {
                    iAt->inRings.push_back(idxR);
                }
            }
             
            idxR++;
        }
    }
    
    extern void setAnglesSPSigns(std::vector<AtomDict>  & tAtoms,
                                 std::vector<AngleDict> & tAngles)
    {
        for (std::vector<AngleDict>::iterator iAn=tAngles.begin();
                iAn !=tAngles.end(); iAn++)
        {
            iAn->atomsSPStats.clear();
            // std::cout << "XXXXX angle " << std::endl;
            for (std::vector<int>::iterator iAt=iAn->atoms.begin();
                    iAt !=iAn->atoms.end(); iAt++)
            {
                //std::cout << "Atom " << tAtoms[*iAt].id  
                //          << "\t" << " hybrid " << tAtoms[*iAt].hybrid << std::endl;
                iAn->atomsSPStats[tAtoms[*iAt].id] =  tAtoms[*iAt].hybrid;
            }
        }
    }
    
    extern bool confirmPlaneByChiralVol(std::vector<AtomDict> & tAtoms,
                                        std::vector<AtomDict>::iterator tA)
    {
        bool tP = false;
        // Ignore whatever sign the chiral volume is, just check if a plane formed
        std::vector<REAL> vec1, vec2, vec3;
        std::cout << "Atom name " << tA->id << std::endl;
        std::cout << "number of bonding " << tA->connAtoms.size() << std::endl;
        // std::cout << " coordinates ? " << tA->coordExist << std::endl;
        
        if (tA->connAtoms.size() >=3)
        {
            int nH=0;
            for (unsigned i=0; i < 3; i++)
            {
                if (tAtoms[tA->connAtoms[i]].chemType.find("H") !=std::string::npos)
                {
                    nH+=1;
                }
            }
            for (unsigned i=0; i < 3; i++)
            {
                vec1.push_back(tAtoms[tA->connAtoms[0]].coords[i]-tA->coords[i]);
                vec2.push_back(tAtoms[tA->connAtoms[1]].coords[i]-tA->coords[i]);
                vec3.push_back(tAtoms[tA->connAtoms[2]].coords[i]-tA->coords[i]);
            }
            
            REAL aVol=calNormalizedChiralVol(vec1, vec2, vec3);
            
            std::cout << "Chi vol " << aVol << std::endl;
            std::cout << "number of H connected " << nH << std::endl;
            
            REAL tB;
            if (nH>1)
            {
                tB=0.25;
            }
            else
            {
                tB=0.15;
            }
            if (fabs(aVol) < tB)
            {
                tP = true;
            }
        }
        return tP;
    }
    
    extern bool confirmPlaneByAngle(std::vector<AtomDict> & tAtoms,
                                    std::vector<AtomDict>::iterator tA,
                                    REAL                    tCri)
    
    {
        bool lP = false;
        // Ignore whatever sign the chiral volume is, just check if a plane formed
        std::vector<REAL> vec1, vec2, vec3;
        std::cout << "Atom name " << tA->id << std::endl;
        std::cout << "number of bonding " << tA->connAtoms.size() << std::endl;
        // std::cout << " coordinates ? " << tA->coordExist << std::endl;
        
        if (tA->connAtoms.size() >=3)
        {
            int nH=0;
            for (unsigned i=0; i < 3; i++)
            {
                if (tAtoms[tA->connAtoms[i]].chemType.find("H") !=std::string::npos)
                {
                    nH+=1;
                }
            }
            for (unsigned i=0; i < 3; i++)
            {
                vec1.push_back(tAtoms[tA->connAtoms[0]].coords[i]-tA->coords[i]);
                vec2.push_back(tAtoms[tA->connAtoms[1]].coords[i]-tA->coords[i]);
                vec3.push_back(tAtoms[tA->connAtoms[2]].coords[i]-tA->coords[i]);
            }
            
            lP= checkPlaneAng3V(vec1, vec2, vec3, tCri);
        }
        
        return lP;
    }
    
    
    
    extern std::string strTransSP(int tSP)
    {
        std::string sSP("SP-NON");
        if (tSP==1)
        {
            sSP="SP1";
        }
        else if (tSP==2)
        {
            sSP="SP2";
        }
        else if (tSP==3)
        {
            sSP="SP3";
        }
        else if (tSP==4)
        {
            sSP="SPD";
        }
        return sSP;
    }
    
    extern bool checkBridgeStruct(std::vector<AtomDict> & tAtoms,
                                  std::vector<RingDict> & tRings,
                                  int                     anchorIdx)
    {
        bool inB = false;
        std::cout << "Check bridge now: " << std::endl; 
        std::cout << "Check atom " << tAtoms[anchorIdx].id << std::endl;
        std::cout << "It is in rings " << tAtoms[anchorIdx].inRings.size() << std::endl;
        
        if (tAtoms[anchorIdx].inRings.size() >1)
        {   
            std::vector<int> nShare1, nShare2;
            for (std::vector<int>::iterator iConn=tAtoms[anchorIdx].connAtoms.begin();
                     iConn!=tAtoms[anchorIdx].connAtoms.end(); iConn++)
            {
                nShare1.clear();
                std::cout << "connected atom " << tAtoms[*iConn].id << std::endl;
                if (tAtoms[*iConn].inRings.size() > 1)
                {
                    for (std::vector<int>::iterator iR1=tAtoms[*iConn].inRings.begin();
                                    iR1 !=tAtoms[*iConn].inRings.end(); iR1++)
                    {
                        if (std::find(tAtoms[anchorIdx].inRings.begin(),tAtoms[anchorIdx].inRings.end(),
                                     *iR1) !=tAtoms[anchorIdx].inRings.end())
                        {
                            nShare1.push_back(*iR1);
                        }
                    }
                }
                
                if (nShare1.size() > 1 && tAtoms[*iConn].connAtoms.size() > 1)
                {
                    for (std::vector<int>::iterator iCo2=tAtoms[*iConn].connAtoms.begin();
                            iCo2!=tAtoms[*iConn].connAtoms.end(); iCo2++)
                    {
                        if (*iCo2 != anchorIdx && tAtoms[*iCo2].inRings.size() > 1)
                        {
                            nShare2.clear();
                            for (std::vector<int>::iterator iR=tAtoms[*iCo2].inRings.begin();
                                    iR !=tAtoms[*iCo2].inRings.end(); iR++)
                            {
                                if (std::find(tAtoms[anchorIdx].inRings.begin(),tAtoms[anchorIdx].inRings.end(),
                                     *iR) !=tAtoms[anchorIdx].inRings.end())
                                {
                                    nShare2.push_back(*iR);
                                }
                            }
                            if(nShare2.size() > 1)
                            {
                                std::cout << " atom " << tAtoms[anchorIdx].id 
                                          << " is in the bridge structure " << std::endl;
                                inB=true;
                                break;
                            }
                        }
                    }
                }
                
                if (inB)
                {
                    break;
                }
            }
        }
        return inB;
    }
    
    // Check protonation state of an atom and decide how many H atoms will added 
    // to bond it.
    
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  std::vector<AtomDict>   & tAtoms,
                                  std::vector<BondDict>   & tBonds,
                                  PeriodicTable & tTab)
    {
        // REAL aNumH=0.0;
        
        REAL Order  = getTotalBondOrder(tBonds, tAtoms, tIA);
        REAL Diff1  = Order -tIA->formalCharge;
        std::cout << "For atom " << tIA->id << std::endl;
        std::cout << "Total bond order is " << Order << std::endl;
        std::cout << "formal charge is " << tIA->formalCharge << std::endl;
        
        std::map<int, int> Vals;
        
        Vals[1] = tTab.elements[tIA->chemType]["val"];
        if (tTab.extraValences.find(tIA->chemType) !=tTab.extraValences.end())
        {
            int i=2;
            for (std::vector<int>::iterator iM=tTab.extraValences[tIA->chemType].begin();
                    iM !=tTab.extraValences[tIA->chemType].end(); iM++)
            {
                 Vals[i] = *iM;
                 i++;
            }
        }
        
        // Now check all possibilities of element valences 
        REAL minD=100.0;
        
        for (std::map<int, int>::iterator iM=Vals.begin();
                iM !=Vals.end(); iM++)
        {
            REAL Diff2  = (REAL)iM->second
                           -Diff1;
            std::cout << "Val " << iM->second << " Diff1 " << Diff1 << std::endl;
            
            std::cout <<" diff2 is " << Diff2 << std::endl;
            if (fabs(Diff2) <0.000001)
            {
                return 0.0;
            }
            else if (Diff2 >0 && Diff2 < minD )
            {
                minD = Diff2;
                if (minD >8)
                {
                    std::cout << "Bond order or valance error " << std::endl;
                    exit(1);
                }
            }
            
        }
        
        if (minD >8)
        {
            std::cout << "Bond order or valance error " << std::endl;
            exit(1);
        }
        
        std::cout << minD <<" H atom should be added to bind atom " 
                  << tIA->id << std::endl;
        
        return minD;
   
    }
    
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  Molecule   & tMol, 
                                  PeriodicTable & tTab)
    {
        // REAL aNumH=0.0;
        
        REAL Order  = getTotalBondOrder(tMol, tIA);
        REAL Diff1  = Order -tIA->formalCharge;
        std::cout << "For atom " << tIA->id << std::endl;
        std::cout << "Total bond order is " << Order << std::endl;
        std::cout << "formal charge is " << tIA->formalCharge << std::endl;
        
        std::map<int, int> Vals;
        
        Vals[1] = tTab.elements[tIA->chemType]["val"];
        if (tTab.extraValences.find(tIA->chemType) !=tTab.extraValences.end())
        {
            int i=2;
            for (std::vector<int>::iterator iM=tTab.extraValences[tIA->chemType].begin();
                    iM !=tTab.extraValences[tIA->chemType].end(); iM++)
            {
                 Vals[i] = *iM;
                 i++;
            }
        }
        
        // Now check all possibilities of element valences 
        REAL minD=100.0;
        
        for (std::map<int, int>::iterator iM=Vals.begin();
                iM !=Vals.end(); iM++)
        {
            REAL Diff2  = (REAL)iM->second
                           -Diff1;
      
            //std::cout <<" diff2 is " << Diff2 << std::endl;
            if (fabs(Diff2) <0.000001)
            {
                return 0.0;
            }
            else if (Diff2 >0 && Diff2 < minD )
            {
                minD = Diff2;
            }
            
            if (minD >8)
            {
                std::cout << "Bond order or valance error " << std::endl;
                exit(1);
            }
            
            
        }
        
        std::cout << minD <<" H atom should be added to bind atom " 
                  << tIA->id << std::endl;
        
        return minD;
   
    }
    
    
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder)
    {
        REAL aNumH = 0.0;
        StrUpper(tIA->resName);
        // Currently using PH=7 
        // 1. O in Carboxy-Terminus, PKa=2, in ASP and GLU pKa=4
        // No H need to be added
        
        // 2. O in Tyrosine, Pka=10
        //std::cout << "Protonated stat " << tIA->id << std::endl;
        //std::cout << "total order << " << tTolBondOrder << std::endl;
        if (tTolBondOrder==1.0 && tIA->charge==0 )
        {
            aNumH=1.0;
          
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        StrUpper(tIA->resName);
        REAL tVal = getTotalBondOrder(tMol, tIA);
        // Currently using PH=7 
        // 1. O in Carboxy-Terminus, PKa=2, in ASP and GLU pKa=4
        // No H need to be added
        
        // 2. O in Tyrosine, Pka=10
        if (tVal==1.0 && tIA->charge==0 )
        {
            aNumH=1.0;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder, 
                                std::vector<AtomDict> & tAllAtoms)
    {
        REAL aNumH = 0.0;
        StrUpper(tIA->resName);
        // REAL tVal = getTotalBondOrder(tMol, tIA);
        // std::cout << "tVal for N " << tVal << std::endl;
        
        if ( tTolBondOrder < 4)
        {
            if(tIA->resName.compare("LYS")==0 
               || tIA->resName.compare("ARG")==0)
            {
                REAL tD = 4 - tTolBondOrder;  
                aNumH=tD;
                tIA->formalCharge = tD;
            }
            else
            {
                // check Amino-Terminus, pKa=10 
                int tC=0, tH=0, tOther=0;
                for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                        iNB != tIA->connAtoms.end(); iNB++)
                {
                    if(tAllAtoms[*iNB].chemType.compare("C")==0)
                    {
                        tC++;
                    }
                    else if(tAllAtoms[*iNB].chemType.compare("H")==0)
                    {
                        tH++;
                    }
                    else
                    {
                        tOther++;
                    }
                }
                
                if (tC==1 && tOther==0)
                {
                    //It is Amino-Terminus, need to put a charge here
                    aNumH =4 - tTolBondOrder;
                    
                }
                
            }
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        //StrUpper(tIA->resName);
        REAL tVal = getTotalBondOrder(tMol, tIA);
        
        
        // std::cout << "tVal for N " << tVal << std::endl;
        
        
        if (tVal < 4)
        {
            if(tIA->resName.compare("LYS")==0 
               || tIA->resName.compare("ARG")==0)
            {
                REAL tD = 4 - tVal;  
                aNumH=tD;
                tIA->formalCharge = tD;
            }
            else
            {
                // check Amino-Terminus, pKa=10 
                int tC=0, tH=0, tOther=0;
                for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                        iNB != tIA->connAtoms.end(); iNB++)
                {
                    if(tMol.atoms[*iNB].chemType.compare("C")==0)
                    {
                        tC++;
                    }
                    else if(tMol.atoms[*iNB].chemType.compare("H")==0)
                    {
                        tH++;
                    }
                    else
                    {
                        tOther++;
                    }
                }
                
                if (tC==1 && tOther==0)
                {
                    //It is Amino-Terminus, need to put a charge here
                    aNumH =4 - tVal;
                    
                }
                
            }
        }
        
        
        return aNumH;
    }
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder)
    {
        REAL aNumH = 0.0;
                
        if (tTolBondOrder ==1 && tIA->formalCharge==0 )
        {
            aNumH=1.0;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        
        REAL tVal = getTotalBondOrder(tMol, tIA);
        
        if (tVal==1 && tIA->formalCharge==0 )
        {
            aNumH=1.0;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder)
    {
        REAL aNumH = 0.0;
        
        if ( tTolBondOrder < 4)
        {
            aNumH = 4-tTolBondOrder;
        }
        else if (tTolBondOrder > 4)
        {
            std::cout << "Warning: Atom " << tIA->id << " has bond order "
                      << tTolBondOrder << std::endl;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        
        REAL tVal = getTotalBondOrder(tMol, tIA);
        
        if (tVal < 4)
        {
            aNumH = 4-tVal;
        }
        else if (tVal > 4)
        {
            std::cout << "Warning: Atom " << tIA->id << " has bond order "
                      << tVal << std::endl;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
        
    }
    
    // These two functions obtain bond-order of one bond and total bond-order of 
    // around one atom when all individual bond-order values exist
    
    extern REAL getTotalBondOrder(Molecule   & tMol, 
                                  std::vector<AtomDict>::iterator tIA)
    {
        REAL tVal = 0.0;
        for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                    iNB !=tIA->connAtoms.end(); iNB++)
        {
            REAL aOrd = getBondOrder(tMol, tIA->seriNum, *iNB);
            //std::cout << "bond order between atom " << tIA->seriNum+1 
            //          << " and " << tMol.atoms[*iNB].seriNum+1
            //          << " is " << aOrd << std::endl;
            if (aOrd >0)
            {
                tVal +=aOrd;
                // std::cout << "total order now " << tVal << std::endl;
            }
            else
            {
                std::cout << "Can not find the bond between atoms " << tIA->id 
                          << " serial number " << tIA->seriNum + 1
                          << " and " << tMol.atoms[*iNB].id 
                          << " serial number " << tMol.atoms[*iNB].seriNum+1
                          << std::endl;
                std::cout << "Some thing is wrong in the Bond list " << std::endl;
                exit(1);
            }
        }
        
        return tVal;
    }
    
    
    extern REAL getTotalBondOrder(std::vector<BondDict>   & tBonds, 
                                  std::vector<AtomDict>   & tAtoms,
                                  std::vector<AtomDict>::iterator tIA)
    {
        REAL tVal = 0.0;
        //std::cout << "see atom " << tIA->id << std::endl;
        
        //std::cout << "It connected " << tIA->connAtoms.size() << " atoms " << std::endl;
        
        for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                    iNB !=tIA->connAtoms.end(); iNB++)
        {  
            //std::cout << "connected atom " << *iNB << std::endl;
            REAL aOrd = getBondOrder(tBonds, tIA->seriNum, *iNB);
            //std::cout << "bond order between atom " << tIA->seriNum+1 
            //          << " and " << tAtoms[*iNB].seriNum+1
            //          << " is " << aOrd << std::endl;
            if (aOrd >0)
            {
                tVal +=aOrd;
                //std::cout << "total order now " << tVal << std::endl;
            }
            else
            {
                std::cout << "Can not find the bond between atoms " << tIA->id 
                          << " serial number " << tIA->seriNum + 1
                          << " and " << tAtoms[*iNB].id 
                          << " serial number " << tAtoms[*iNB].seriNum+1
                          << std::endl;
                std::cout << "Some thing is wrong in the Bond list " << std::endl;
                exit(1);
            }
        }
        
        return tVal;
    }
    
    
    
    extern REAL getBondOrder(Molecule & tMol, int tIdx1, int tIdx2)
    {
        REAL tOrd = -1.0;
        
        for (std::vector<BondDict>::iterator iB=tMol.bonds.begin();
                iB !=tMol.bonds.end(); iB++)
        {
            if ((iB->atomsIdx[0]==tIdx1 && iB->atomsIdx[1]==tIdx2)
                || (iB->atomsIdx[0]==tIdx2 && iB->atomsIdx[1]==tIdx1))
            {
                tOrd = StrToReal(iB->order);
                if (tOrd ==4.0)
                {
                    tOrd = 1.5;
                }
                
                //std::cout << "Bond order " << iB->order << std::endl
                //          << tOrd << std::endl;
                break;
            }
        }
        
        return tOrd;
        
    }
    
    
    extern REAL getBondOrder(std::vector<BondDict> & tBonds,
                             int tIdx1, int tIdx2)
    {
        REAL tOrd = -1.0;
        
        for (std::vector<BondDict>::iterator iB=tBonds.begin();
                iB !=tBonds.end(); iB++)
        {
            if ((iB->atomsIdx[0]==tIdx1 && iB->atomsIdx[1]==tIdx2)
                || (iB->atomsIdx[0]==tIdx2 && iB->atomsIdx[1]==tIdx1))
            {
                tOrd = StrToOrder2(iB->order);
                
                if (tOrd ==4.0)
                {
                    tOrd = 1.5;
                }
                
                //std::cout << "Bond order " << iB->order << std::endl
                //          << tOrd << std::endl;
                break;
            }
        }
        
        return tOrd;
        
    }
    
    extern REAL getBondOrderOneAtom(std::vector<BondDict> tBonds,
                                    std::vector<AtomDict> tAtoms,
                                    int tIdx)
    {
        
        REAL tOrd = -1.0;
        
        for (std::vector<int>::iterator iNB=tAtoms[tIdx].connAtoms.begin();
                iNB !=tAtoms[tIdx].connAtoms.end(); iNB++)
        {
            
        }
            
        
        return tOrd;
        
    }
    
    // Automated assignment of bond orders to a set of bonds
    extern void setAllBondOrders(std::vector<AtomDict> & tAtoms,
                                 std::vector<BondDict> & tBonds)
    {
        std::vector<int> unAssigned;
        
        for (unsigned int i=0; i < tBonds.size(); i++)
        {
            if (tAtoms[tBonds[i].atomsIdx[0]].chemType=="H" 
               || tAtoms[tBonds[i].atomsIdx[1]].chemType=="H"
               || tAtoms[tBonds[i].atomsIdx[0]].bondingIdx==3
               || tAtoms[tBonds[i].atomsIdx[1]].bondingIdx==3)
            {
                tBonds[i].order = "single";
            }
            else
            {
                unAssigned.push_back(i);
            }
            
        }
    }
    
    
    
    // Set Coordinates for added H atoms
    
    extern void  setOneHAtomCoordsSP3(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA)
    {   
        
        std::cout << "set H atom " << tIA->id  << " coords " << std::endl;
        
        TransCoords   aTransTool;
        
        int           aCen=tIA->connAtoms[0];
        
        int           root1=-1, root2=-1;
       
        for (std::vector<int>::iterator iR1=tAtoms[aCen].connAtoms.begin();
                iR1 !=tAtoms[aCen].connAtoms.end(); iR1++)
        {
            root1 = -1;
            root2 = -1;
            if(tAtoms[*iR1].chemType !="H")
            {
                root1 = tAtoms[*iR1].seriNum;
                for (std::vector<int>::iterator iR2=tAtoms[root1].connAtoms.begin();
                        iR2 !=tAtoms[root1].connAtoms.end(); iR2++)
                {
                    if (tAtoms[*iR2].chemType !="H" 
                        && tAtoms[*iR2].id !=tAtoms[aCen].id)
                    {
                        std::vector<REAL> v1, v2, v3;
                        for (unsigned iX=0; iX < tAtoms[*iR2].coords.size();
                                iX++)
                        {
                            v1.push_back(tAtoms[*iR1].coords[iX]-tAtoms[*iR2].coords[iX]);
                            v2.push_back(tAtoms[aCen].coords[iX]-tAtoms[*iR1].coords[iX]);
                        }
                        crossP2V(v1, v2, v3);
                        if (lengthV(v3) >0.00000001)
                        {
                            root2=tAtoms[*iR2].seriNum;
                            break;
                        }
                    }
                }
            }
            
            if (root1 !=-1 && root2 !=-1)
            {
                break;
            }
            
            
        }
        
        if (root1 !=-1 && root2 !=-1)
        {
            //std::cout << "root1 atom " << tAtoms[root1].id << std::endl;
            //std::cout << "root2 atom " << tAtoms[root2].id << std::endl;
            if (tIA->coords.size() ==0)
            {
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
            }
            
            std::vector<int> refAtmIdx;
            
            for (std::vector<int>::iterator iR1=tAtoms[aCen].connAtoms.begin();
                iR1 !=tAtoms[aCen].connAtoms.end(); iR1++)
            {
                if (tAtoms[*iR1].id != tAtoms[root1].id 
                     && tAtoms[*iR1].id != tIA->id
                     && tAtoms[*iR1].coordExist)
                {
                    refAtmIdx.push_back(tAtoms[*iR1].seriNum);
                }
            }
            
            REAL d = 0.95, alpha= 109.5*PI180;  // the guide values
            REAL aTor;
            
            //std::cout << "number of ref atom is " << refAtmIdx.size() << std::endl;
            
            if  (refAtmIdx.size() ==1)
            {
                aTor  = getTorsion(tAtoms[root2], tAtoms[root1], tAtoms[aCen],tAtoms[refAtmIdx[0]])*PI180
                             + 120.0*PI180;
            }
            else if  (refAtmIdx.size() ==2)
            {
                int tIdx;
                REAL tTor, tTor1, tTor2;
                /*
                for (unsigned i=0; i < 3; i++)
                {
                    std::cout << tAtoms[refAtmIdx[0]].id << "  "  << i << "  "
                               << tAtoms[refAtmIdx[0]].coords[i] << std::endl;
                    std::cout << tAtoms[refAtmIdx[1]].id << "  "  << i << "  "
                               << tAtoms[refAtmIdx[1]].coords[i] << std::endl;
                }
                */
                tTor1 = getTorsion(tAtoms[root2], tAtoms[root1], tAtoms[aCen],tAtoms[refAtmIdx[0]]);
                tTor2 = getTorsion(tAtoms[root2], tAtoms[root1], tAtoms[aCen],tAtoms[refAtmIdx[1]]);
                
                if (tTor1 > tTor2)
                {
                    tTor  = tTor1;
                    tIdx  = refAtmIdx[0];
                    tTor1 = tTor2;
                    refAtmIdx[0] = refAtmIdx[1];
                    tTor2 = tTor;
                    refAtmIdx[1] = tIdx;
                    
                }
                
                if (tTor2-tTor1 >=180)
                {
                    aTor  = (tTor1 + 120.0)*PI180;
                }
                else
                {
                    aTor  =  (tTor2 + 120.0)*PI180;
                }
                if (aTor > PI)
                {
                    aTor = aTor-2*PI;
                }
            }
            else
            {
                // only two atoms connect to tAtoms[aCen]
                aTor  = 120*PI180;
            }
            
            std::cout << "Tor is " << aTor*PID180 << std::endl;
            std::cout << "atom1 " << tAtoms[root2].id << " atom2 " 
                      << tAtoms[root1].id << " atom3 " << tAtoms[aCen].id
                      << " H atom added " << tIA->id << std::endl;
            aTransTool.growOneAtom(tAtoms[root2], tAtoms[root1], tAtoms[aCen], 
                                       tIA, d, alpha, aTor);
            std::cout << "check distance " << distanceV(tAtoms[aCen].coords,
                                                        tIA->coords) 
                                           << std::endl;
            std::vector<REAL> v1, v2;
            for (unsigned i=0; i < tIA->coords.size(); i++)
            {
                v1.push_back(tAtoms[aCen].coords[i]-tAtoms[root1].coords[i]);
                v2.push_back(tIA->coords[i]-tAtoms[aCen].coords[i]);
            }
            
            std::cout << "check angle: " << (PI-getAngle2V(v1, v2))*PID180
                                         << std::endl;
            
            //for (unsigned i=0; i < 3; i++)
            //{
            //    std::cout << i << " " << tIA->coords[i] << std::endl;
            //}
            //std::cout << "coords Tor = " << getTorsion(tAtoms[root2], tAtoms[root1], tAtoms[aCen], *tIA)
            //          << std::endl << std::endl;;
            tIA->coordExist = true;
        }
        
    }
    
    extern void  setOneHAtomCoordsSP2(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA)
    {
        std::cout << "set H atom " << tIA->id  << " coords " << std::endl;
        TransCoords   aTransTool;
        
        int           aCen=tIA->connAtoms[0] , aRef=-1;
        
        int           root1=-1, root2=-1;
       
        for (std::vector<int>::iterator iR1=tAtoms[aCen].connAtoms.begin();
                iR1 !=tAtoms[aCen].connAtoms.end(); iR1++)
        {
            root1 = -1;
            root2 = -1;
            if(tAtoms[*iR1].chemType !="H")
            {
                root1 = tAtoms[*iR1].seriNum;
                for (std::vector<int>::iterator iR2=tAtoms[root1].connAtoms.begin();
                        iR2 !=tAtoms[root1].connAtoms.end(); iR2++)
                {
                    if (tAtoms[*iR2].chemType !="H" 
                        && tAtoms[*iR2].id !=tAtoms[aCen].id)
                    {
                        root2=tAtoms[*iR2].seriNum;
                        break;
                    }
                }
            }
            if (root1 !=-1 && root2 !=-1)
            {
                break;
            }
        }
        
        if (root1 !=-1 && root2 !=-1)
        {
            if (tIA->coords.size() ==0)
            {
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
            }
            for (std::vector<int>::iterator iR1=tAtoms[aCen].connAtoms.begin();
                iR1 !=tAtoms[aCen].connAtoms.end(); iR1++)
            {
                if (tAtoms[*iR1].id != tAtoms[root1].id 
                    && tAtoms[*iR1].id != tIA->id)
                {
                    aRef = tAtoms[*iR1].seriNum;
                    break;
                }
                
            }
            
            REAL d = 0.95, alpha= 120*PI180;  // the guide values
            REAL aTor;
            if (aRef==-1)
            {
                // only two atoms connect to tAtoms[aCen]
                aTor  = PI;
            }
            else
            {
                aTor  = getTorsion(tAtoms[root2], tAtoms[root1], tAtoms[aCen],tAtoms[aRef])*PI180
                             + PI;
            }
            
            aTransTool.growOneAtom(tAtoms[root2], tAtoms[root1], tAtoms[aCen], 
                                       tIA, d, alpha, aTor);
            std::cout << "Tor is " << aTor*PID180 << std::endl;
            std::cout << "atom1 " << tAtoms[root2].id << " atom2 " 
                      << tAtoms[root1].id << " atom3 " << tAtoms[aCen].id
                      << " H atom added " << tIA->id << std::endl;
            std::cout << "check distance: " << distanceV(tAtoms[aCen].coords,
                                                        tIA->coords) 
                                           << std::endl;
            std::vector<REAL> v1, v2;
            for (unsigned i=0; i < tIA->coords.size(); i++)
            {
                v1.push_back(tAtoms[aCen].coords[i]-tAtoms[root1].coords[i]);
                v2.push_back(tIA->coords[i]-tAtoms[aCen].coords[i]);
            }
            
            std::cout << "check angle: " << (PI-getAngle2V(v1, v2))*PID180
                                         << std::endl;
            tIA->coordExist = true;
        }
        else
        {
            std::cout <<"The molecule is an isolated 4-atom cluster with "
                      << tAtoms[aCen].id << " at the center " << std::endl;
            exit(1);
        }
        
    }
    
    extern void  setOneHAtomCoordsSP(std::vector<AtomDict> & tAtoms,
                                     std::vector<AtomDict>::iterator tIA)
    {
        TransCoords   aTransTool;
        
        int           aCen=tIA->connAtoms[0];
        
        int           root1=-1, root2=-1;
       
        for (std::vector<int>::iterator iR1=tAtoms[aCen].connAtoms.begin();
                iR1 !=tAtoms[aCen].connAtoms.end(); iR1++)
        {
            root1 = -1;
            root2 = -1;
            if(tAtoms[*iR1].chemType !="H")
            {
                root1 = tAtoms[*iR1].seriNum;
                for (std::vector<int>::iterator iR2=tAtoms[root1].connAtoms.begin();
                        iR2 !=tAtoms[root1].connAtoms.end(); iR2++)
                {
                    if (tAtoms[*iR2].chemType !="H" 
                        && tAtoms[*iR2].id !=tAtoms[aCen].id)
                    {
                        root2=tAtoms[*iR2].seriNum;
                        break;
                    }
                }
            }
            if (root1 !=-1 && root2 !=-1)
            {
                break;
            }
        }
        
        if (root1 !=-1 && root2 !=-1)
        {
            if (tIA->coords.size() ==0)
            {
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
                tIA->coords.push_back(0.0);
            }
            REAL d = 0.95, alpha= 120*PI180;  // the guide values
            REAL aTor =  PI;
            aTransTool.growOneAtom(tAtoms[root2], tAtoms[root1], tAtoms[aCen], 
                                       tIA, d, alpha, aTor);
            tIA->coordExist = true;
        }
        else
        {
            std::cout <<"The molecule is an isolated 4-atom cluster with "
                      << tAtoms[aCen].id << " at the center " << std::endl;
            exit(1);
        }
    }
    
    extern void checkAllStereo(FileName tMdlIn, FileName tPdbIn,
                               FileName tPdbOut)
    {
        MolSdfFile aMolFile(tMdlIn,  std::ios::in);
        
        for (std::vector<Molecule>::iterator iMol=aMolFile.allMols.begin();
                iMol !=aMolFile.allMols.end(); iMol++)
        {
            checkStereoOneMol(iMol, tPdbIn);
        }
    }
    
    extern void   setAllAtomEXcessElectrons(std::vector<AtomDict> & tAtoms)
    {
        // Temporarily. should use the Periodic table object created before.
        
        std::vector<std::string> orgTab;
        initOrgTable(orgTab);
        
        std::map<ID, std::vector<int> > orgElemValMap;
        orgElemValMap["C"].push_back(4);
        orgElemValMap["N"].push_back(3);
        orgElemValMap["N"].push_back(5);
        orgElemValMap["O"].push_back(2);
        orgElemValMap["N"].push_back(3);
        orgElemValMap["S"].push_back(2);
        orgElemValMap["S"].push_back(4);
        orgElemValMap["S"].push_back(6);
        orgElemValMap["P"].push_back(5);
        orgElemValMap["SE"].push_back(2);
        orgElemValMap["SE"].push_back(4);
        orgElemValMap["SE"].push_back(6);
        orgElemValMap["B"].push_back(3);
        
        orgElemValMap["H"].push_back(1);
        orgElemValMap["F"].push_back(1);
        orgElemValMap["CL"].push_back(1);
        orgElemValMap["BR"].push_back(1);
        orgElemValMap["I"].push_back(1);
        orgElemValMap["AT"].push_back(1);
        
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            ID aElm = iAt->chemType;
            StrUpper(aElm);
            
            if (std::find(orgTab.begin(), orgTab.end(), aElm) != orgTab.end())
            {
                int valSize = (int)orgElemValMap[aElm].size();
                int orgNB =0;
                for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                            iNB != iAt->connAtoms.end(); iNB++)
                {   
                    ID aNBElm = tAtoms[*iNB].chemType;
                    StrUpper(aNBElm);
                    if(std::find(orgTab.begin(), orgTab.end(), aNBElm) != orgTab.end())
                    {
                        orgNB++;
                    }
                }
                int nExEls = orgElemValMap[aElm][0] + iAt->formalCharge - orgNB;
                if (nExEls < 0)
                {
                    int i = 1;
                
                    while (i < valSize)
                    {
                        nExEls = orgElemValMap[aElm][i] + iAt->formalCharge - iAt->connAtoms.size();
                        if (nExEls >=0)
                        {
                            break;
                        }
                        i++;
                    }
                }
            
            
                if (nExEls < 0)
                {
                    std::cout << "Error : the number of connections to Atom "
                              << iAt->id << ": "
                              << iAt->connAtoms.size() << " is larger than the valence "
                              << orgElemValMap[aElm][valSize-1] << " permits." 
                              << std::endl << "The formal charge is  "
                              << iAt->formalCharge
                              << std::endl;
                    std::cout << "Atom " << iAt->id << " has following connections: "
                              << std::endl;
                    for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                            iNB != iAt->connAtoms.end(); iNB++)
                    {
                        std::cout << tAtoms[*iNB].id << std::endl;
                    }
                }
                else
                {
                    iAt->excessElec = nExEls;
                }
            
                //std::cout << "For atom " << iAt->id << " : " << std::endl
                //          << "it connects " << iAt->connAtoms.size() 
                //          << " atom(s) " << std::endl
                //          << "its formal charge is " << iAt->formalCharge << std::endl 
                //          << "its number of EX electrons is " << iAt->excessElec 
                //          << std::endl;
            }   
        }   
    }
    
    extern int  sumExElectrons(std::vector<AtomDict> & tAtoms)
    {
        int sumExElec =0;
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            sumExElec += iAt->excessElec;
        }
        
        return sumExElec; 
    }
    extern void checkStereoOneMol(std::vector<Molecule>::iterator tMol,
                                  FileName tPdbIn)
    {
        
    }
    
    // New section for Protonated - de-protonated forms of functional groups
    
    extern void checkProtonatedCarBoxylicTerminus(std::vector<AtomDict>::iterator tIA,
                                               std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               REAL tPh, REAL tPka, 
                                               std::vector<int> tDoneAtomsIdx)
    {
        
        // This is for Carboxy-Terminus
        
        if (tIA->connAtoms.size() == 3)
        {
            //tDoneAtomsIdx.push_back(tIA->seriNum);
            // (1) Check if it is a carboxylic acid 
            std::vector<int> tAtIdx;   // remember which atoms are done
            
            std::map<int, int> Oxys, Others;    // key: serial number, value: number of connection 
            std::vector<int> doubleOs, singleOs;
            tAtIdx.push_back(tIA->seriNum);
            for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                   iNB !=tIA->connAtoms.end(); iNB++)
            {
                if (tAtoms[*iNB].chemType =="O")
                {
                    
                    int nS = (int)tAtoms[*iNB].connAtoms.size();
                    Oxys[*iNB]   = nS;
                    if (nS==2)
                    {
                        doubleOs.push_back(*iNB);
                    }
                    else if (nS==1)
                    {
                        singleOs.push_back(*iNB);
                    }
                }
                else 
                {
                    Others[*iNB] = tAtoms[*iNB].connAtoms.size();
                }
                tAtIdx.push_back(*iNB);
            }
             
            if (doubleOs.size() !=0)
            {
                for (std::vector<int>::iterator iOa=doubleOs.begin();
                        iOa !=doubleOs.end(); iOa++)
                {
                    for (std::vector<int>::iterator iOb=tAtoms[*iOa].connAtoms.begin();
                            iOb!=tAtoms[*iOa].connAtoms.end(); iOb++)
                    {
                        if (*iOb !=tIA->seriNum)
                        {
                            tAtIdx.push_back(*iOb);
                        }
                    }
                }
            }
                
            // Possible carboxylic acid Requirement:
            // C connects to two O and one other atom
            if (Oxys.size() ==2 && Others.size()==1)
            {
                std::cout << "It is a possible c-Ter " << std::endl;
                std::cout << "It connects  " << doubleOs.size() << " O atoms of double linked " << std::endl;
                std::cout << "It connects  " << singleOs.size() << " O atoms of single linked " << std::endl;
                // (2) Further check if it is a terminus
                // Two O atoms must connected with either one atom 
                // or two atoms with one atom is a H atom 
                bool isTer = true;
                
                if((doubleOs.size()==1 && singleOs.size()==1)
                   || singleOs.size()==2)
                {
                    if (doubleOs.size()==1)
                    {
                        std::cout << "O atom " << tAtoms[doubleOs[0]].id << std::endl;
                        for (unsigned iOc=0; iOc < tAtoms[doubleOs[0]].connAtoms.size(); 
                             iOc++)
                        {
                            
                            if (tAtoms[tAtoms[doubleOs[0]].connAtoms[iOc]].seriNum !=tIA->seriNum)
                            {
                                // Not carbon atom itself
                                std::cout << "atom " << tAtoms[tAtoms[doubleOs[0]].connAtoms[iOc]].id << std::endl;
                                if(tAtoms[tAtoms[doubleOs[0]].connAtoms[iOc]].chemType !="H")
                                {
                                    // Not H atom either, then not terminus 
                                    isTer = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    isTer = false;
                }
                
                
                if (isTer)
                {
                    // It is a C-terminus 
                    std::cout << "It is a C-ter " << std::endl;
                    // Default values: tPH = 7.0 tPka =4.0 (3-5) 
                    if (tPh > tPka) 
                    {
                        
                        // In this case, de-protonation should happen 
                        if (doubleOs.size()==1 && singleOs.size()==1)
                        {
                            // de-protonate the double connected O atom
                            int idxH =-1;
                            int posH =-1;
                            for (unsigned iOc=0; iOc !=tAtoms[doubleOs[0]].connAtoms.size();
                                 iOc++)
                            {
                                if (tAtoms[tAtoms[doubleOs[0]].connAtoms[iOc]].chemType=="H")
                                {
                                    idxH = tAtoms[doubleOs[0]].connAtoms[iOc];
                                    posH = iOc;
                                    break;
                                }
                            }
                            
                            if (idxH !=-1)
                            {
                                
                                // erase H atom and set formal charge on O to -1.0
                                ID aH =tAtoms[idxH].id;
                                tAtoms.erase(tAtoms.begin() + idxH);
                                std::cout << "Atom " << aH << " is deleted " << std::endl;
                                tAtoms[doubleOs[0]].formalCharge = -1.0;
                                tAtoms[doubleOs[0]].connAtoms.erase(tAtoms[doubleOs[0]].connAtoms.begin()+posH);
                                
                                // delete the bond between O and H
                                int aBIdx=getBond(tBonds, doubleOs[0], idxH);
                                if (aBIdx != -1)
                                {
                                    tBonds.erase(tBonds.begin() + aBIdx);
                                }
                            }
                        }  
                    }
                }
                for (std::vector<int>::iterator iDone = tAtIdx.begin();
                        iDone != tAtIdx.end(); iDone++)
                {
                    tDoneAtomsIdx.push_back(*iDone);
                }
            }
        }
    }
    
    
    extern void checkProtonatedSulfuricAcids(std::vector<AtomDict>::iterator tIA,
                                               std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               REAL tPH, std::vector<int> tDoneAtomsIdx)
    {
        // According to the empirical rules  from the examples (see the note,
        // assuming default pH value 7.0, pka value varies depending on connection)
        // (1) if S atom has two connections will not be de-protonated, except 
        // the non-H connected atom is on a aromatic ring.
        // (2) if S atom connects 3 atoms one of which is double bonded. 
        //     do nothing at the moment.
        // (3) if S atom connects to 4 atoms. Situation varies 
        tDoneAtomsIdx.push_back(tIA->seriNum);
        REAL tPKa;
        if (tIA->connAtoms.size()==2)
        {
            //De-protonated only when a user gives a very high PH value.
            tPKa = 10.0;
            
            if (tPH > tPKa)
            {
                
                for (std::vector<int>::iterator iNA=tIA->connAtoms.begin();
                        iNA !=tIA->connAtoms.end(); iNA++)
                {
                    
                    tDoneAtomsIdx.push_back(*iNA);
                    if (tAtoms[*iNA].chemType !="H")
                    {
                        tAtoms[*iNA].formalCharge = -1.0;
                    }
                }
            }
        }
        else if (tIA->connAtoms.size()==4)
        {
            std::vector<int> Oxys;
            
            for (std::vector<int>::iterator iOB=tIA->connAtoms.begin();
                    iOB !=tIA->connAtoms.end(); iOB++)
            {
                
                tDoneAtomsIdx.push_back(*iOB);
                if (tAtoms[*iOB].chemType=="O")
                {
                    Oxys.push_back(*iOB);
                }
            }
            
            if (Oxys.size()==3 || Oxys.size()==4)
            {
                // de-protonate  O atoms with a single bond.
                for (unsigned i=0; i < Oxys.size(); i++)
                {
                    int aBIdx=getBond(tBonds, tIA->seriNum, Oxys[i]);
                    if (aBIdx !=-1)
                    {
                        std::string tStr(tBonds[aBIdx].order);
                        StrUpper(tStr);
                        if (tStr.compare("1")==0 
                            || tStr.substr(0,4).compare("SING")==0)
                        {
                            tAtoms[Oxys[i]].formalCharge = -1.0;
                        }
                    }
                }
            }
        }
    }
    
    extern void checkProtonatedAminoTerminus(std::vector<AtomDict>::iterator tIA,
                                      std::vector<AtomDict> & tAtoms,
                                      std::vector<BondDict> & tBonds,
                                      REAL tPH, REAL tPka,
                                      std::vector<int> tDoneAtomsIdx)
    {
        
        tDoneAtomsIdx.push_back(tIA->seriNum);
       
        if (tIA->connAtoms.size() ==3)
        {
            std::cout << "It is a possible N-ter " << std::endl;
            
            std::vector<int> tAtIdx;   // remember which atoms are done
            tAtIdx.push_back(tIA->seriNum);
            
            std::vector<int> Hs, others;
            for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                    iNB !=tIA->connAtoms.end(); iNB++)
            {               
                if (tAtoms[*iNB].chemType == "H")
                {
                    Hs.push_back(*iNB);
                }
                else
                {
                    others.push_back(*iNB);
                }
                tAtIdx.push_back(*iNB);
            }
            std::cout << "It connects to " << Hs.size() << " H atoms " << std::endl
                      << " and " << others.size() << " other atoms " << std::endl;
            
            if (Hs.size()==2 && others.size()==1)
            {
                // It is Amino-Terminus
                if (tPH < tPka)
                {
                    tIA->formalCharge =1.0;
                    std::cout << "Atom " << tIA->id 
                              << " has formal charge " << tIA->formalCharge 
                              << std::endl;
                            
                    for (std::vector<int>::iterator iDone=tAtIdx.begin();
                           iDone !=tAtIdx.end(); iDone++)
                    {
                        tDoneAtomsIdx.push_back(*iDone);
                    }
                }
            }
        }
    }
    
    extern void checkProtonatedPAcids(std::vector<AtomDict>::iterator tIA,
                                      std::vector<AtomDict> & tAtoms,
                                      std::vector<BondDict> & tBonds,
                                      REAL tPH, std::vector<int> tDoneAtomsIdx)
    {
        if (tIA->connAtoms.size()==4)
        {
            
            tDoneAtomsIdx.push_back(tIA->seriNum);
            
            std::vector<int> Oxys, OxysS, bondOs, bondOtherSingle;
            
            for (std::vector<int>::iterator iOB=tIA->connAtoms.begin();
                    iOB!=tIA->connAtoms.end(); iOB++)
            {
                
                tDoneAtomsIdx.push_back(*iOB);
                
                std::string tStr;
                int aIdxB=getBond(tBonds, tIA->seriNum, *iOB);
                if (aIdxB !=-1)
                {
                    tStr.append(tBonds[aIdxB].order);
                    StrUpper(tStr);
                }
                if (tAtoms[*iOB].chemType=="O")
                {
                    Oxys.push_back(*iOB);
                    
                    if (tStr.compare("1")==0) 
                    {
                        bondOs.push_back(aIdxB);
                    }
                    else if (tStr.size() >=4)
                    {
                        if (tStr.substr(0,4).compare("SING")==0)
                        {
                            bondOs.push_back(aIdxB);
                        }
                    }
                }
                else
                {
                    if (tStr.compare("1")==0) 
                    {
                        bondOtherSingle.push_back(aIdxB);
                    }
                    else if (tStr.size() >=4)
                    {
                        if (tStr.substr(0,4).compare("SING")==0)
                        {
                            bondOtherSingle.push_back(aIdxB);
                        }
                    }
                }
                
                
            }
            
            
            if (Oxys.size()==2 && OxysS.size()==2 && bondOtherSingle.size()==2 )
            {
                // PO2 + 2 other single bonds
                tAtoms[OxysS[0]].formalCharge=-1.0;
                tBonds[bondOs[1]].order = "2";
            }
            else if (Oxys.size()==2 && OxysS.size()==1 && bondOtherSingle.size()==2)
            {
                
                // PO2 with one single O bond and one double O bond
                // the user assigned the valence correctly, just
                // make sure the formal charges are right.
                tAtoms[OxysS[0]].formalCharge=-1.0;
                
            }
            else if (Oxys.size()==3 && OxysS.size()==3 && bondOtherSingle.size()==1)
            {
                // PO3 + one other single bonds
                tAtoms[OxysS[0]].formalCharge=-1.0;
                tAtoms[OxysS[1]].formalCharge=-1.0;
                tBonds[bondOs[2]].order = "2";
            }
            else if (Oxys.size()==3 && OxysS.size()==2 && bondOtherSingle.size()==1)
            {
                // PO3 with two single O bonds and one double O bond
                // the user assigns the valence correctly, just
                // make sure the formal charges are right.
                tAtoms[OxysS[0]].formalCharge=-1.0;
                tAtoms[OxysS[1]].formalCharge=-1.0;
            }
            else if (Oxys.size()==4 && OxysS.size()==4)
            {
                // PO3 + one other single bonds
                tAtoms[OxysS[0]].formalCharge=-1.0;
                tAtoms[OxysS[1]].formalCharge=-1.0;
                tAtoms[OxysS[2]].formalCharge=-1.0;
                tBonds[bondOs[3]].order = "2";
            }
            else if (Oxys.size()==4 && OxysS.size()==3)
            {
                // PO4 with three single O bonds and one double O bond
                // the user assigns the valence correctly, just
                // make sure the formal charges are right.
                tAtoms[OxysS[0]].formalCharge=-1.0;
                tAtoms[OxysS[1]].formalCharge=-1.0;
                tAtoms[OxysS[2]].formalCharge=-1.0;
            }
            
        }
        
    }
    
    // The following function replaces the previous function addHAtomToMols().
    extern void ProtonateFunctionGroupInOneMol(std::vector<AtomDict>  & tAtoms,
                                               std::vector<BondDict>  & tBonds, 
                                               std::vector<int>       & tAddHIdx)
    {
        
        
        REAL PH=7.0;
        
        std::vector<int> doneAtomList;
        
        std::vector<std::string> orgElems;
        orgElems.push_back("C");
        orgElems.push_back("N");
        orgElems.push_back("S");
        orgElems.push_back("P");
        
        
        PeriodicTable aPTab;
        
        // Correct the formal charges for all atoms in the molecule.
        
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            
            if (std::find(orgElems.begin(), orgElems.end(), iA->chemType)
                 !=orgElems.end() 
                 && std::find(doneAtomList.begin(), doneAtomList.end(), iA->seriNum) 
                    ==doneAtomList.end())
            {
                // need to check protonated form
                if (iA->chemType =="C")
                {
                    REAL aPka = 4.0;
                    std::cout << "Get a C atom " << iA->id << std::endl;
                    checkProtonatedCarBoxylicTerminus(iA, tAtoms, tBonds, PH, aPka, doneAtomList);         
                }
                else if (iA->chemType =="N")
                {
                    REAL aPka = 10.0;
                    std::cout << "Get a N atom " << iA->id << std::endl;
                    checkProtonatedAminoTerminus(iA, tAtoms, tBonds, PH, aPka, doneAtomList);
                }   
                else if (iA->chemType =="S")
                {
                    checkProtonatedSulfuricAcids(iA, tAtoms, tBonds, PH, doneAtomList);
                }   
                else if (iA->chemType =="P")
                {
                    checkProtonatedPAcids(iA, tAtoms, tBonds, PH, doneAtomList);
                }  
                
            }
        }
        
        
        // Now add H atoms based on previous calculations
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            
            if (iA->chemType !="H")
            {
                int addH =(int)checkProtonateAll(iA, tAtoms, tBonds,  aPTab);
                
                // REAL addH =checkProtonated(iA, tIdxMol);
                if (addH != 0)
                {
                    adjustHAtoms(tAtoms, tBonds, iA->seriNum, addH, tAddHIdx);
                }
                
            }
        } 
        
        
         
    }
    
    // The following function overwrites addHAtoms in MolSdfFile class.
    // Add or delete H atoms 
    extern void adjustHAtoms(std::vector<AtomDict> & tAtoms,
                             std::vector<BondDict> & tBonds,
                             int                     tIdxAtm,
                             int                     tAddH,
                             std::vector<int>      & tHAtmIdxs)
    {
        // For setup H atom's name. How many H atoms it already connects
        
        int nH=0;
        std::vector<int> Hs;   // index of H atoms already connected. For deprotonized mainly
        
        for (std::vector<int>::iterator iH=tAtoms[tIdxAtm].connAtoms.begin();
                iH !=tAtoms[tIdxAtm].connAtoms.end(); iH++)
        {
            if (tAtoms[*iH].chemType =="H")
            {
                Hs.push_back(tAtoms[*iH].seriNum);
                nH++;
            }
        }
        
        std::string tS;
        getDigitSec(tAtoms[tIdxAtm].id, tS);
        
        if (tAddH > 0)
        {
            int j= 1;
            for (int i=0; i < (int)tAddH; i++)
            {
                AtomDict aH;
                aH.chemType = "H";
            
                if (tAtoms[tIdxAtm].chemType=="C")
                {
                    aH.id  = aH.chemType + tS + IntToStr(i+1);
                }
                else
                {
                    if ((nH+j)==1)
                    {
                        aH.id  = aH.chemType + tAtoms[tIdxAtm].id;
                    }
                    else
                    {
                        aH.id  = aH.chemType + tAtoms[tIdxAtm].id +  IntToStr(nH+j);
                    }
                }
                
                aH.seriNum =  (int)tAtoms.size();
                aH.connAtoms.push_back(tAtoms[tIdxAtm].seriNum);
                tAtoms[tIdxAtm].connAtoms.push_back(aH.seriNum);
                tAtoms[tIdxAtm].connHAtoms.push_back(aH.seriNum);
                tHAtmIdxs.push_back((int)tAtoms.size());
                tAtoms.push_back(aH);
                
                BondDict   aB;
                aB.seriNum = (int)tBonds.size();
                aB.order   = "single";
                aB.orderN  = 1.0;
                aB.atoms.push_back(tAtoms[tIdxAtm].id);
                aB.atoms.push_back(aH.id);
                aB.atomsIdx.push_back(tAtoms[tIdxAtm].seriNum);
                aB.atomsIdx.push_back(aH.seriNum);
                aB.fullAtoms[aH.id] = aH.seriNum;
                aB.fullAtoms[tAtoms[tIdxAtm].id] = tIdxAtm;
                tBonds.push_back(aB);
                
            }
        }
        else if (tAddH <0)
        {
            for (unsigned i=0; i < abs(tAddH); i++)
            {
                if (i < Hs.size())
                {
                    // Delete that H atom
                    tAtoms.erase(tAtoms.begin()+Hs[i]);
                    std::cout << "H atom " << tAtoms[Hs[i]].id 
                              << " has been deleted " << std::endl;
                }
                else
                {
                    std::cout << "Bug in de-protonation: not enough H atoms to remove " 
                              << std::endl;
                    break;
                }
                
                // delete the bond associated with that H atom
                int aHB=getBond(tBonds, tAtoms[tIdxAtm].seriNum, Hs[i]);
                if (aHB !=-1)
                {
                    tBonds.erase(tBonds.begin()+aHB);
                    std::cout << "Bond of atom " << tAtoms[tIdxAtm].id 
                              << " and H atom " <<  tAtoms[Hs[i]].id 
                              << " has been deleted " << std::endl;
                }
                else
                {
                    std::cout << "Bug: there is no bond between atom " 
                              << tAtoms[tIdxAtm].id 
                              << " and H atom " <<  tAtoms[Hs[i]].id << std::endl;
                }
            }
        }      
    }
    
    extern void setAllAddedHAtomCoords(std::vector<AtomDict> & tAtoms,
                                       std::vector<int>      & tHAtmIdxs)
    {
        
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            if (std::find(tHAtmIdxs.begin(), tHAtmIdxs.end(), iA->seriNum) !=tHAtmIdxs.end())
            {
                if (iA->connAtoms.size() ==1 )
                {
                    if (tAtoms[iA->connAtoms[0]].bondingIdx==3)
                    {
                        
                        setOneHAtomCoordsSP3(tAtoms, iA);
                    }
                    else if (tAtoms[iA->connAtoms[0]].bondingIdx==2)
                    {
                        setOneHAtomCoordsSP2(tAtoms, iA);
                    }
                    else if (tAtoms[iA->connAtoms[0]].bondingIdx==1)
                    {
                        setOneHAtomCoordsSP(tAtoms, iA);
                    }
                    else
                    {
                        std::cout << "The bondingIdx for atom " 
                                  <<  tAtoms[iA->connAtoms[0]].id
                                  << " is " << tAtoms[iA->connAtoms[0]].bondingIdx
                                  << std::endl;
                    }
                }
                else
                {
                    std::cout << "Bug. H atom " << iA->id << " connects to "
                              << iA->connAtoms.size() << " atom(s)" << std::endl;
                    exit(1);
                }
            }
        }
    }
    
    
    
    
    void extern setAromPlanes(std::vector<AtomDict>  & tAtoms,
                              std::vector<RingDict>  & tRings, 
                              std::vector<PlaneDict> & tPlans)
    {
        // 1. 
        std::vector<std::vector<int> > sP21Rs;
        for (std::vector<RingDict>::iterator iRi=tRings.begin();
                iRi !=tRings.end(); iRi++)
        {
            if (iRi->isPlanar)
            {
                std::vector<int> tAtms;
                for (std::vector<AtomDict>::iterator iAt=iRi->atoms.begin();
                        iAt !=iRi->atoms.end(); iAt++)
                {
                    tAtms.push_back(iAt->seriNum);
                }
                sP21Rs.push_back(tAtms);
            }
        }
        
        std::cout << "Number of rings containing sp2 and sp bonding atoms only is: "
                  << sP21Rs.size() << std::endl;
     
    }
    
    HuckelMOSuite::HuckelMOSuite()
    {
        // Default values, allow an input table to change them 
        orgAlphas["C"]  = -11.2;
        orgAlphas["N"]  = -11.2;
        orgAlphas["O"]  = -11.2;
        orgAlphas["S"]  = -11.2;
        orgAlphas["P"]  = -11.2;
        orgAlphas["B"]  = -11.2;
        orgAlphas["SE"] = -11.2;
        
        // Halogen atoms and H atom should not use orgAlphas
        
        orgBetas["C"]   = -0.7;
        orgBetas["N"]   = -0.7;
        orgBetas["O"]   = -0.7;
        orgBetas["S"]   = -0.7;
        orgBetas["P"]   = -0.7;
        orgBetas["B"]   = -0.7;
        orgAlphas["SE"] = -0.7;
        
        
    }
        
    HuckelMOSuite::~HuckelMOSuite()
    {
    }
     
    void HuckelMOSuite::setWorkMode(int tMode)
    {
        workMode = tMode;
    }
    void HuckelMOSuite::execute(std::vector<AtomDict>& tAtoms,
                                std::vector<BondDict>& tBonds)
    {
        if (workMode==1)
        {
            //Pick up pi electrons and fragment the system 
            initiaExElecs(tAtoms);
            PickOddAtoms(tAtoms);
            if (withExAtomIdxs.size() > 0)
            {
                partitionSysToSubGraphs(tAtoms);
                MOSolver(tAtoms);
            }
            else
            {
                std::cout << "No atoms are with free pi electrons, "
                      << std::endl << "and all  bonds are set " << std::endl;
            }
        
        }
        else if (workMode==2)
        {
            setBondOrderInSys(tAtoms, tBonds);
        }
        
        if (withExAtomIdxs.size() > 0)
        {
            partitionSysToSubGraphs(tAtoms);
            MOSolver(tAtoms);
            BondTrans(tBonds);
        }
        else
        {
            std::cout << "No atoms are with free pi electrons, "
                      << std::endl << "and all  bonds are set " << std::endl;
        }
    }
    
    void HuckelMOSuite::initiaExElecs(std::vector<AtomDict>& tAtoms)
    {
        // Initialization 
        setAllAtomEXcessElectrons(tAtoms);
        
        //Pick up pi electrons at the first stage 
        for (std::vector<AtomDict>::iterator iAt= tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            std::cout << "Atom " << iAt->id << "   " << iAt->seriNum << std::endl;
            if (iAt->excessElec !=0)
            {
                withExAtomIdxs.push_back(iAt->seriNum);
            }
            else
            {
                zeroExAtomIdxs.push_back(iAt->seriNum);
            }
        }
        
        // Check 
        /*
        std::cout << "Now those atoms are considered to be with pi electrons " 
                  << std::endl;
        for (std::vector<int>::iterator iAt=withExAtomIdxs.begin();
                iAt != withExAtomIdxs.end(); iAt++)
        {
            std::cout << "Atom " << tAtoms[*iAt].id << std::endl;
        }
         */
    }
    
    void HuckelMOSuite::PickPiElectrons(std::vector<AtomDict>& tAtoms)
    {
        // for the rounds after initialization 
        std::vector<int> tmpIdxList;
        for(std::vector<int>::iterator iAt=withExAtomIdxs.begin();
                iAt != withExAtomIdxs.end(); iAt++)
        {
            if (tAtoms[*iAt].excessElec == 0)
            {
                zeroExAtomIdxs.push_back(*iAt);
            }
            else
            {
                tmpIdxList.push_back(*iAt);
            }
        }
        
        withExAtomIdxs.clear();
        
        for(std::vector<int>::iterator iAt=tmpIdxList.begin();
                iAt != tmpIdxList.end(); iAt++)
        {
            withExAtomIdxs.push_back(*iAt);
        }
        
        
    } 
    
    void HuckelMOSuite::PickOddAtoms(std::vector<AtomDict>& tAtoms)
    {
        // Check those atoms appearing 1 ex-electron, but really not
        // Singly bonded O atoms
        std::vector<int> tmpIdxList1, tmpIdxList2;
        for (std::vector<int>::iterator aIdx=withExAtomIdxs.begin();
                aIdx !=withExAtomIdxs.end(); aIdx++)
        {
            if (tAtoms[*aIdx].connAtoms.size()==1 ) 
                // && tAtoms[*aIdx].chemType.compare("H") !=0) 
                // H should already be deleted
                // Remember both atoms at two ends of the connection
            {
                oddAtomIdxs[*aIdx]=tAtoms[*aIdx].connAtoms[0];
                tmpIdxList2.push_back(*aIdx);
                tmpIdxList2.push_back(tAtoms[*aIdx].connAtoms[0]);
            }
            else
            {
                tmpIdxList1.push_back(*aIdx);
            }
        }
        
        withExAtomIdxs.clear();
        
        // second round
        for (std::vector<int>::iterator aIdx=tmpIdxList1.begin();
                aIdx !=tmpIdxList1.end(); aIdx++)
        {
            if (std::find(tmpIdxList2.begin(), tmpIdxList2.end(), *aIdx)
                ==tmpIdxList2.end())
            {
                withExAtomIdxs.push_back(*aIdx);
            }
        }
        
        for (std::vector<int>::iterator iT=tmpIdxList2.begin();
                iT!=tmpIdxList2.end(); iT++)
        {
            zeroExAtomIdxs.push_back(*iT);
        }
        
        
        // Check 
        /*
        if (withExAtomIdxs.size()==0)
        {
            std::cout << "No atoms are considered to have free pi electrons"
                      << std::endl;
        }
        else
        {
            std::cout << "After singly connected non-H atoms are excluded. " << std::endl
                      << " Those atoms considered to be with pi electrons are: " 
                      << std::endl;
            for (std::vector<int>::iterator iAt=withExAtomIdxs.begin();
                     iAt != withExAtomIdxs.end(); iAt++)
            {
                std::cout << "Atom " << tAtoms[*iAt].id << std::endl;
            }
        }
        */
        std::cout << "Those are bonds excluded. " 
                  << std::endl;
        for (std::map<int, int>::iterator iAt=oddAtomIdxs.begin();
                iAt != oddAtomIdxs.end(); iAt++)
        {
            std::cout << "Bond between atom " << tAtoms[iAt->first].id 
                      << " and atom " << tAtoms[iAt->second].id << std::endl;
        }   
    }
    
    void HuckelMOSuite::setInitBondOrder(std::vector<AtomDict>             & tAtoms, 
                                         std::vector<BondDict>             & tBonds, 
                                         std::vector<int>                  & tCBondIdx, 
                                         std::map<int, std::vector<int> >  & tDelConn,
                                         std::map<int, int>                & tRemainVal)
    {
        
        tCBondIdx.clear();
        tDelConn.clear();
        
        // Dealt with the singly connected atoms first
        for (std::vector<int>::iterator iZA=zeroExAtomIdxs.begin();
                iZA !=zeroExAtomIdxs.end(); iZA++)
        {
            
            if (tAtoms[*iZA].connAtoms.size()==1)
            {
                int idxB = getBond(tBonds, *iZA, tAtoms[*iZA].connAtoms[0]);
                if (tAtoms[*iZA].chemType.compare("H")==0
                    || tAtoms[*iZA].chemType.compare("F")==0
                    || tAtoms[*iZA].chemType.compare("CL")==0
                    || tAtoms[*iZA].chemType.compare("BR")==0
                    || tAtoms[*iZA].chemType.compare("I")==0
                    || tAtoms[*iZA].chemType.compare("AR")==0)
                {
                    tCBondIdx.push_back(idxB);
                    tBonds[idxB].orderN = 1.0;
                    tRemainVal[tAtoms[*iZA].seriNum]--;
                    tRemainVal[tAtoms[*iZA].connAtoms[0]]--;
                    tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                    tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                }
                else if (tAtoms[*iZA].chemType.compare("O")==0)
                {
                    // assume O has charge -1, or 0
                    if (tAtoms[*iZA].charge==-1.0)
                    {
                        tCBondIdx.push_back(idxB);
                        tBonds[idxB].orderN = 1.0;
                        tRemainVal[tAtoms[*iZA].seriNum]--;
                        tRemainVal[tAtoms[*iZA].connAtoms[0]]--;
                        tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                        tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                    }
                    else if (tAtoms[*iZA].charge==0.0)
                    {
                        tCBondIdx.push_back(idxB);
                        tBonds[idxB].orderN = 2.0;
                        tRemainVal[tAtoms[*iZA].seriNum]-=2;
                        tRemainVal[tAtoms[*iZA].connAtoms[0]]-=2;
                        tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                        tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                    }
                }
                else if (tAtoms[*iZA].chemType.compare("N")==0 
                        || tAtoms[*iZA].chemType.compare("B")==0)
                {
                    if (tAtoms[*iZA].charge==0.0)
                    {
                        tCBondIdx.push_back(idxB);
                        tBonds[idxB].orderN = 3.0;
                        tRemainVal[tAtoms[*iZA].seriNum]-=3;
                        tRemainVal[tAtoms[*iZA].connAtoms[0]]-=3;
                        tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                        tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                    }
                    else
                    {
                        std::cout << "Singly connected N atom "
                                  << tAtoms[*iZA].id << " has charge "
                                  << tAtoms[*iZA].charge 
                                  << " check! " << std::endl;
                        exit(1);
                    }
                }
                else if (tAtoms[*iZA].chemType.compare("S")==0)
                {
                    if (tAtoms[*iZA].charge==-1.0)
                    {
                        tCBondIdx.push_back(idxB);
                        tBonds[idxB].orderN = 1.0;
                        tRemainVal[tAtoms[*iZA].seriNum]--;
                        tRemainVal[tAtoms[*iZA].connAtoms[0]]--;
                        tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                        tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                    }
                    else if (tAtoms[*iZA].charge==0.0
                             || tAtoms[*iZA].charge==-2.0)
                    {
                        tCBondIdx.push_back(idxB);
                        tBonds[idxB].orderN = 2.0;
                        tRemainVal[tAtoms[*iZA].seriNum]-=2;
                        tRemainVal[tAtoms[*iZA].connAtoms[0]]-=2;
                        tDelConn[tAtoms[*iZA].connAtoms[0]].push_back(tAtoms[*iZA].seriNum);
                        tDelConn[tAtoms[*iZA].seriNum].push_back(tAtoms[*iZA].connAtoms[0]);
                    }
                    else
                    {
                        std::cout << "Singly connected S atom "
                                  << tAtoms[*iZA].id << " has charge "
                                  << tAtoms[*iZA].charge 
                                  << " check! " << std::endl;
                        exit(1);
                    }
                }
                else if(tAtoms[*iZA].chemType.compare("C")==0 
                        || tAtoms[*iZA].chemType.compare("P")==0)
                {
                    std::cout << "Atom "
                              << tAtoms[*iZA].id << " connects to one atom, check!"
                              << std::endl;
                    exit(1);
                }   
                std::cout << "Atom " << tAtoms[*iZA].id 
                          << " has remain val " 
                          << tRemainVal[tAtoms[*iZA].seriNum] << std::endl;
            }
        }
       
        std::cout << "First step " << std::endl;
        std::cout << tCBondIdx.size() << "are decided bonds: " << std::endl;
        for (std::vector<int>::iterator iB=tCBondIdx.begin(); 
                iB !=tCBondIdx.end(); iB++)
        {
            std::cout << "Bond idx " << *iB << " which is between atom " 
                      << tBonds[*iB].atoms[0] << " and "
                      << tBonds[*iB].atoms[1] << std::endl;
            
        }
        
        
        // 2. Deal with definite single bonds
        for (std::vector<int>::iterator iZA=zeroExAtomIdxs.begin();
                iZA !=zeroExAtomIdxs.end(); iZA++)
        {
            if (tRemainVal[tAtoms[*iZA].seriNum] !=0)
            {
                //std::cout << "Atom " << tAtoms[*iZA].id 
                //          << " has remain val " 
                //          << tRemainVal[tAtoms[*iZA].seriNum] << std::endl;
                
                int nReConns;
                if (tDelConn.find(tAtoms[*iZA].seriNum) !=tDelConn.end())
                {
                    nReConns= tAtoms[*iZA].connAtoms.size() 
                              - tDelConn[tAtoms[*iZA].seriNum].size();
                }
                else
                {
                    nReConns= tAtoms[*iZA].connAtoms.size();
                }
                
                if (tRemainVal[tAtoms[*iZA].seriNum]==nReConns)
                {
                    // all bonds around this atom are single
                    for (std::vector<int>::iterator iCo=tAtoms[*iZA].connAtoms.begin();
                            iCo != tAtoms[*iZA].connAtoms.end(); iCo++)
                    {
                        int idxB = getBond(tBonds, *iZA, *iCo);
                        //std::cout << "Bond idx " << idxB << std::endl;
                        if (std::find(tCBondIdx.begin(), tCBondIdx.end(), idxB)
                                      ==tCBondIdx.end())
                        {
                            //std::cout << "Bond with atom " << tBonds[idxB].atoms[0]
                            //          << " and " << tBonds[idxB].atoms[1] << std::endl;
                            
                            tCBondIdx.push_back(idxB);
                            tBonds[idxB].orderN = 1.0;
                            tRemainVal[tAtoms[*iZA].seriNum]--;
                            tRemainVal[*iCo]--;
                            tDelConn[tAtoms[*iZA].seriNum].push_back(*iCo);
                            tDelConn[*iCo].push_back(tAtoms[*iZA].seriNum);
                        }
                    }
                }
            }
        }
        
        
        std::cout << "After second step " << std::endl;
        for (std::map<int, int>::iterator iZA=tRemainVal.begin();
                iZA !=tRemainVal.end(); iZA++)
        {
            if (iZA->second > 0)
            {
                std::cout << "Atom " << tAtoms[iZA->first].id 
                          << " has remain val " 
                          << iZA->second << std::endl;
            }
        }
        std::cout << "Total number of bonds " << tBonds.size() << std::endl;
        std::cout << tCBondIdx.size() << "are decided bonds: " << std::endl;
        for (std::vector<int>::iterator iB=tCBondIdx.begin(); 
                iB !=tCBondIdx.end(); iB++)
        {
            std::cout << "Bond idx " << *iB << " which is between atom " 
                      << tBonds[*iB].atoms[0] << " and "
                      << tBonds[*iB].atoms[1] << " with bond-order " 
                      << tBonds[*iB].orderN << std::endl;
            
        }
        std::cout << "Bonds remained unset are " << std::endl;
        for (unsigned i=0; i < tBonds.size(); i++)
        {
            if (std::find(tCBondIdx.begin(), tCBondIdx.end(), i)
                    ==tCBondIdx.end())
            {
                std::cout << "Bond between " << tBonds[i].atoms[0]
                          << " and " << tBonds[i].atoms[1] << std::endl;
            }
        }
        
        
        // 3. do while loop to set possible single, double and triple bonds 
        //    until nothing can be done
        
        std::map<int, std::vector<int> >  remainConns;
        
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            
            for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                    iNB !=iAt->connAtoms.end(); iNB++)
            {
                if (tDelConn.find(iAt->seriNum) != tDelConn.end())
                {
                    if (std::find(tDelConn[iAt->seriNum].begin(), 
                                  tDelConn[iAt->seriNum].end(), *iNB)
                             ==tDelConn[iAt->seriNum].end())
                    {
                        remainConns[iAt->seriNum].push_back(*iNB);
                    }
                }
                else
                {
                    remainConns[iAt->seriNum].push_back(*iNB);
                }
            }
        }
        
        for (std::map<int, std::vector<int> >::iterator iZA=remainConns.begin();
                iZA !=remainConns.end(); iZA++)
        {
            std::cout << "Atom " << tAtoms[iZA->first].id 
                                 << " has remain val " 
                                 << tRemainVal[iZA->first] 
                                 << " and connections " 
                                 << iZA->second.size() << std::endl;
        }
        
        
        int nDone;
        do 
        {
            nDone =0;
            setProBondOrdersOneLoop(nDone, tAtoms, tBonds, tCBondIdx, 
                                    remainConns, tDelConn, tRemainVal);
            std::cout << "nDone in this round " << nDone << std::endl;
        }while (tCBondIdx.size() != tBonds.size()
                && nDone !=0);
        
        withExAtomIdxs.clear();
        
        for (std::map<int, int>::iterator iVal=tRemainVal.begin();
                iVal != tRemainVal.end(); iVal++)
        {
            if (iVal->second > 0)
            {
                withExAtomIdxs.push_back(iVal->first);
            }
        }
        
        // Check
        std::cout << "Total number of bonds " << tBonds.size() << std::endl;
        
        std::cout << "Total number of bonds set " << tCBondIdx.size()
                  << std::endl;
        std::cout << "total number of atoms with free pi electrons is "
                  << withExAtomIdxs.size() << std::endl;
        
        if (tBonds.size() != tCBondIdx.size())
        {
            unsigned n = tBonds.size() - tCBondIdx.size(); 
            std::cout << n << " bonds remain to be decided their order " 
                      << std::endl;
            
            std::cout << "Those bonds are: " << std::endl;
            for (unsigned i=0; i < tBonds.size(); i++)
            {
                if (std::find(tCBondIdx.begin(), tCBondIdx.end(), i)
                      ==tCBondIdx.end())
                {
                    std::cout << "Bond between atom " 
                              << tBonds[i].atoms[0] << " and "
                              << tBonds[i].atoms[1] << std::endl;
                }
            }
        }
        else
        {
            std::cout << "All bonds have definite bond order. No need to go further "
                      << std::endl;
            std::cout << "Those bonds are: " << std::endl;
            for (std::vector<BondDict>::iterator iBo = tBonds.begin();
                    iBo != tBonds.end(); iBo++)
            {
                std::cout << "Bond between atom " << iBo->atoms[0]
                          << " and " << iBo->atoms[1] 
                          << " with bond order " << iBo->orderN << std::endl;
            }
        }
    }
    
    
    void HuckelMOSuite::setProBondOrdersOneLoop(int & nDone, std::vector<AtomDict>& tAtoms, 
                                                std::vector<BondDict> & tBonds, 
                                                std::vector<int>      & tCBondIdx, 
                                                std::map<int,std::vector<int> >& tRemainConns, 
                                                std::map<int,std::vector<int> >& tDelConn, 
                                                std::map<int,int>& tRemainVal)
    {
        std::vector<std::string> elems;
        elems.push_back("C");
        elems.push_back("N");
        elems.push_back("O");
        elems.push_back("S");
        elems.push_back("P");
        elems.push_back("B");
        
        std::map<int, std::vector<int> > tmpRC;
        for (std::map<int, std::vector<int> >::iterator iZA=tRemainConns.begin();
                iZA !=tRemainConns.end(); iZA++)
        { 
            for (std::vector<int>::iterator iCo=iZA->second.begin();
                    iCo != iZA->second.end(); iCo++)
            {
                tmpRC[iZA->first].push_back(*iCo);
            }
        }
        
        // Deal with the atoms with one remained connection first
        for (std::map<int, std::vector<int> >::iterator iZA=tmpRC.begin();
                iZA !=tmpRC.end(); iZA++)
        {  
            if (tRemainVal[iZA->first] > 0)
            {
                if (iZA->second.size()==1)
                {
                    int idxB = getBond(tBonds, iZA->first, iZA->second[0]);
                    if (std::find(tCBondIdx.begin(), tCBondIdx.end(), idxB)
                            ==tCBondIdx.end())
                    {
                        if (std::find(elems.begin(), elems.end(), tAtoms[iZA->first].chemType)
                                 !=elems.end())
                        {
                            if (tRemainVal[iZA->first]==1)
                            {
                                tCBondIdx.push_back(idxB);
                                tBonds[idxB].orderN = 1.0;
                                tRemainVal[iZA->first]--;
                                tRemainVal[iZA->second[0]]--;
                                tDelConn[iZA->second[0]].push_back(iZA->first);
                                tDelConn[iZA->first].push_back(iZA->second[0]);
                                nDone++;
                            }
                            else if (tRemainVal[iZA->first]==2)
                            {
                                tCBondIdx.push_back(idxB);
                                tBonds[idxB].orderN = 2.0;
                                tRemainVal[iZA->first]-=2;
                                tRemainVal[iZA->second[0]]-=2;
                                tDelConn[iZA->second[0]].push_back(iZA->first);
                                tDelConn[iZA->first].push_back(iZA->second[0]);
                                nDone++;
                            }
                            else if (tRemainVal[iZA->first]==3)
                            {
                                tCBondIdx.push_back(idxB);
                                tBonds[idxB].orderN = 3.0;
                                tRemainVal[iZA->first]-=3;
                                tRemainVal[iZA->second[0]]-=3;
                                tDelConn[iZA->second[0]].push_back(iZA->first);
                                tDelConn[iZA->first].push_back(iZA->second[0]);
                                nDone++;
                            }
                        }
                    }
                }
            }
        }
        
        tRemainConns.clear();
        
        for (std::map<int, std::vector<int> >::iterator iZA=tmpRC.begin();
                iZA !=tmpRC.end(); iZA++)
        {  
            for (std::vector<int>::iterator iCo=iZA->second.begin();
                    iCo != iZA->second.end(); iCo++)
            {
                if (tDelConn.find(iZA->first) !=tDelConn.end())
                {
                    if(std::find(tDelConn[iZA->first].begin(),
                                 tDelConn[iZA->first].end(), *iCo)
                            ==tDelConn[iZA->first].end())
                    {
                        tRemainConns[iZA->first].push_back(*iCo);
                    }
                }
                else
                {
                    tRemainConns[iZA->first].push_back(*iCo);
                }
            }
        }
        
    }
    
    void HuckelMOSuite::setBondOrderInSys(std::vector<AtomDict> & tAtoms,
                                          std::vector<BondDict> & tBonds)
    {
        
        initiaExElecs(tAtoms);
        PickOddAtoms(tAtoms);
        
        std::map<int, int> atomCurVals;
        PeriodicTable aPT;
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            atomCurVals[iAt->seriNum] = aPT.elements[iAt->chemType]["val"]
                                        + iAt->charge;
            std::cout << "Atom " << iAt->id << " has valence " 
                      << atomCurVals[iAt->seriNum] << std::endl; 
        }
       
        std::vector<int>    cBonds;
         
        std::map<int, std::vector<int> > delConns;
        setInitBondOrder(tAtoms, tBonds, cBonds, delConns, atomCurVals);
       
    }
    
    void HuckelMOSuite::partitionSysToSubGraphs(std::vector<AtomDict>& tAtoms)
    {
        std::cout << "input number of atoms to the graph partition is " 
                  << withExAtomIdxs.size() << std::endl;
        
        
        std::map<int, int> classNum;
        classNum[0] = 0;
        for (unsigned i=1; i < withExAtomIdxs.size(); i++)
        {
            classNum[i] = i;
            for (unsigned j=0; j <=i-1;j++)
            {
                classNum[j]=classNum[classNum[j]];
                if (std::find(tAtoms[withExAtomIdxs[i]].connAtoms.begin(),
                              tAtoms[withExAtomIdxs[i]].connAtoms.end(), withExAtomIdxs[j])
                          !=tAtoms[withExAtomIdxs[i]].connAtoms.end())
                {
                    classNum[classNum[classNum[j]]]=i;
                }
            }
        }
        
        // final sweeping 
        for (unsigned i=0; i < withExAtomIdxs.size(); i++ )
        {
            classNum[i]=classNum[classNum[i]];    
        }    
        
        std::map<int, std::vector<int> > tAllSubSys;
        
        tAllSubSys.clear();
        
        for (unsigned i=0; i < classNum.size(); i++)
        {
            tAllSubSys[withExAtomIdxs[classNum[i]]].push_back(withExAtomIdxs[i]);
        }
        
        
        int idx=1;
        for (std::map<int, std::vector<int> >::iterator iT=tAllSubSys.begin();
                iT != tAllSubSys.end(); iT++)
        {
            for (std::vector<int>::iterator iV=iT->second.begin(); 
                    iV !=iT->second.end(); iV++)
            {
                allSubGraphs[idx].push_back(*iV);
            }
            idx ++;
        }
        
        // Check
        std::cout << "There are " << allSubGraphs.size() 
                  << " subgraphs." << std::endl;
        for (std::map<int, std::vector<int> >::iterator iCla=allSubGraphs.begin();
                iCla !=allSubGraphs.end(); iCla++)
        {
            std::cout << "For subgraph " << iCla->first 
                      << ", it contains the following atoms "
                      << std::endl;
            for (std::vector<int>::iterator iAt=iCla->second.begin();
                    iAt !=iCla->second.end(); iAt++)
            {
                std::cout << "Atom " << tAtoms[*iAt].id << std::endl;
            }
        }
          
    }
    
    void HuckelMOSuite::setHMatrix(std::vector<AtomDict> & tAtoms,
                                   REAL ** tH,   
                                   std::vector<int> & tSubGraph)
    {
        
        for (unsigned i=0; i < tSubGraph.size(); i++)
        {   
            int curAtmIdx = tSubGraph[i];
            std::string tElem(tAtoms[curAtmIdx].chemType);
            StrUpper(tElem);
            if (orgAlphas.find(tElem) !=  orgAlphas.end())
            {
                tH[i][i] =  orgAlphas[tElem];
                for (unsigned j=0; j < tSubGraph.size(); j++)
                {
                    int curNBAtmidx =  tSubGraph[j];   
                    std::string tNBElem(tAtoms[curNBAtmidx].chemType);
                    StrUpper(tNBElem);
                    if (orgBetas.find(tNBElem)!=orgBetas.end())
                    {
                        if (std::find(tAtoms[curAtmIdx].connAtoms.begin(),
                                      tAtoms[curAtmIdx].connAtoms.end(), curNBAtmidx)
                                  !=  tAtoms[curAtmIdx].connAtoms.end())
                        {
                            tH[i][j]  = orgBetas[tElem];
                        }   
                    }
                    else
                    {
                        std::cout << "Non-organic element or Halogen element  "
                                  << "or H enters calculations. STOP "
                                  << std::endl
                                  << "Atom element symbol " << tNBElem 
                                  << std::endl;
                        exit(1);
                    }
                }
            }
            else
            {
                std::cout << "Non-organic element or Halogen element"
                          << "or H enters calculations. STOP"
                          << std::endl
                          << "Atom element symbol " << tElem 
                          << std::endl;
                exit(1);
            }   
        }
        
        // Check
        std::cout << "Check " << std::endl;
        for (unsigned i=0; i < tSubGraph.size(); i++)
        {
            for (unsigned j=0; j < tSubGraph.size(); j++ )
            {
                if (tH[i][j] !=0.0)
                {
                    std::cout << "H[" << i << "][" << j << "]= " 
                              << tH[i][j] << std::endl;
                    if (i==j)
                    {
                        std::cout << "atom " <<  tAtoms[tSubGraph[i]].id 
                                  << " self energy " << std::endl;
                    }
                    else
                    {
                        if (i < j)
                        {
                            std::cout << "which means atoms " << tAtoms[tSubGraph[i]].id 
                                      << " and " << tAtoms[tSubGraph[j]].id 
                                      << " are are bonded " << std::endl;
                        }
                    }
                }
                    
            }
        }
        
    }
    
    void HuckelMOSuite::getBondOrderFromOrb(int    tNOrbs,
                                            REAL** tEigenVect,
                                            std::vector<AtomDict>  & tAtoms,
                                            std::vector<int> & tSubGraph)
    {
        int nMaxOccp = tNOrbs/2;
        
        for (unsigned i=0; i < tNOrbs; i++)
        {
            for (unsigned j=0; j < tNOrbs; j++)
            {
                int idx1 = tSubGraph[i];
                int idx2 = tSubGraph[j];
                
                if (idx1 < idx2)
                {
                    if (std::find(tAtoms[idx1].connAtoms.begin(),
                                  tAtoms[idx1].connAtoms.end(), idx2) 
                        != tAtoms[idx1].connAtoms.end())
                    {
                        BondOrderFromMO[idx1][idx2] = 0.0;
                        for (unsigned mu=0; mu < nMaxOccp; mu++)
                        {
                            BondOrderFromMO[idx1][idx2] 
                                     +=(tEigenVect[mu][i]*tEigenVect[mu][j]);
                            
                        }
                        // double occupied orbitals 
                        BondOrderFromMO[idx1][idx2] = BondOrderFromMO[idx1][idx2]*2.0;
                        
                        std::cout << "Bond order between atom " << tAtoms[idx1].id
                                  << " and " << tAtoms[idx2].id <<  " is "
                                  << BondOrderFromMO[idx1][idx2] << std::endl;
                       
                    }
                }
            }
        } 
    }
    
    
    void HuckelMOSuite::MOSolver(std::vector<AtomDict>& tAtoms)
    {
        
        
        if (allSubGraphs.size() !=0)
        {
            for (std::map<int, std::vector<int> >::iterator iMa=allSubGraphs.begin();
                    iMa !=allSubGraphs.end(); iMa++)
            {
                std::cout << "Solve the Huckel MO for subgraph " << iMa->first
                          << std::endl;
           
                // Points to REAL are Need for using Jama 
                // for eigenvalue problem. stupid.
                
                REAL ** aH = new REAL * [iMa->second.size()];
                
                for (int i=0; i <iMa->second.size(); i++)
                {
                    aH[i] = new REAL [iMa->second.size()];
                    for (int j=0; j < iMa->second.size(); j++)
                    {
                        aH[i][j] = 0.0;
                    }
                }
                
                setHMatrix(tAtoms, aH, iMa->second); 
                
                REAL *  eigenValue_T = new REAL [iMa->second.size()];
                REAL ** eigenVect_T  = new REAL * [iMa->second.size()];
                for(unsigned i =0; i < iMa->second.size(); i++)
                {
                    eigenVect_T[i] = new REAL [iMa->second.size()];
                }
                
                EigenSolve(iMa->second.size(), aH, eigenValue_T, eigenVect_T);
                std::cout << "The eigenvalues obtained are: " << std::endl;
                
                for (unsigned i=0; i < iMa->second.size(); i++)
                {
                    std::cout <<"##############################################"
                              << std::endl;
                    std::cout << "Eigenvalue : " << eigenValue_T[i] << std::endl;
                    std::cout << "----------------------------------------------"
                              << std::endl;
                    std::cout << "The eigen-vector associated with it is: " 
                              << std::endl;
                    for (unsigned j=0; j < iMa->second.size(); j++)
                    {
                        std::cout << eigenVect_T[i][j] << std::endl;
                    }         
                }
                
                getBondOrderFromOrb(iMa->second.size(), eigenVect_T, 
                                    tAtoms, iMa->second);
                      
            }
        }
    }
    
    void HuckelMOSuite::BondTrans(std::vector<BondDict>& tBonds)
    {
        if (BondOrderFromMO.size() !=0)
        {
            for (std::map<int, std::map<int, REAL> >::iterator iM1=BondOrderFromMO.begin();
                    iM1 != BondOrderFromMO.end(); iM1++)
            {
                for (std::map<int, REAL>::iterator iM2=iM1->second.begin();
                        iM2 != iM1->second.end(); iM2++)
                {
                    int idxB = getBond(tBonds, iM1->first, iM2->first);
                    if (idxB >0 && idxB < tBonds.size())
                    {
                        tBonds[idxB].orderN = iM2->second;
                    }
                    else
                    {
                        std::cout << "Bug: no bond exists between atom "
                                  << iM1->first << " and " << iM2->first
                                  << std::endl;
                    }
                }
            }
        }
    }
    void HuckelMOSuite::outBoAndChList(FileName tFName, 
                                       std::vector<AtomDict>  & tAtoms,
                                       std::vector<BondDict>  & tBonds)
    {
        if (tAtoms.size() !=0 && tBonds.size() !=0)
        {
            std::ofstream outFBA(tFName);
            if(outFBA.is_open())
            { 
                // 1. Atom section
                outFBA << "loop_" << std::endl
                       << "_chem_comp_atom.serial_num" << std::endl
                       << "_chem_comp_atom.atom_id" << std::endl
                       << "_chem_comp_atom.charge" << std::endl;
                for (std::vector<AtomDict>::iterator iA = tAtoms.begin();
                        iA != tAtoms.end(); iA++)
                {
                    outFBA << std::setw(12) << iA->seriNum
                           << std::setw(6) << iA->id 
                           << std::setw(6) << iA->formalCharge << std::endl; 
                }       
                
                outFBA << "loop_" << std::endl
                       << "_chem_comp_bond.atom_id_1" << std::endl
                       << "_chem_comp_bond.atom_id_2" << std::endl
                       << "_chem_comp_bond.type" << std::endl;
                for (std::vector<BondDict>::iterator iB=tBonds.begin();
                          iB !=tBonds.end(); iB++)
                {
                    outFBA << std::setw(12)  << iB->atoms[0]  
                           << std::setw(12)  << iB->atoms[1]  
                           << std::setw(12)  << iB->orderN << std::endl; 
                }
                
            }
        }
        
    }
    
}