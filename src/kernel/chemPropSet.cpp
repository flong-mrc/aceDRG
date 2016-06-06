
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
    
    // Set atom's bonding features (sp, sp2, sp3 and chiral center) based on 
    // the atom's connections.
   
    extern void setAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms)
    {
  
        
        // First round
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
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
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
                        std::vector<ID> atps;
                        for (std::vector<int>::iterator iNA=iAt->connAtoms.begin();
                                iNA !=iAt->connAtoms.end(); iNA++)
                        {
                            if (std::find(atps.begin(), atps.end(), tAtoms[*iNA].chemType)==atps.end())
                            {
                                atps.push_back(tAtoms[*iNA].chemType);
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
    }
    
    extern void modAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds, 
                                               std::vector<AngleDict> & tAngles,
                                               std::vector<RingDict>  & tRings,
                                               int                      tMode)
    {
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
                        std::cout << "connected atom " << tAtoms[*iCo].id 
                                  << "is in rings: " << std::endl;
                        
                        for (std::vector<int>::iterator iR=tAtoms[*iCo].inRings.begin();
                               iR !=tAtoms[*iCo].inRings.end(); iR++)
                        {
                            std::cout << tRings[*iR].rep << std::endl;
                            if(tRings[*iR].isAromatic || tAtoms[*iCo].bondingIdx==2)
                            {
                                std::cout << "It is sp2 related " << std::endl;
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
                std::cout << "lAromRs " << lAromRs << std::endl;
                if (lAromRs)
                {
                    
                    if (checkBridgeStruct(tAtoms, tRings, iA->seriNum))
                    {
                        // std::cout << "Inside 1" << std::endl;
                        iA->chiralIdx  = 5;
                        iA->bondingIdx = 3;
                    }
                    else if (confirmPlaneByChiralVol(tAtoms, iA))
                    {
                        iA->chiralIdx  = 5; // New value 
                        iA->bondingIdx = 2; // still keep it sp2, will use together
                                            // with chiralIdx
                        //std::cout << "inside 2 " << std::endl;
                    }
                    else
                    {
                       // iA->chiralIdx  = 0;
                       // iA->bondingIdx = 2;
                    }
                }
                else
                {
                    
                    if (tMode==1)
                    {
                        if (confirmPlaneByChiralVol(tAtoms, iA))
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
                std::cout << "Its hybridization is sp" << iA->bondingIdx << std::endl; 
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
        
        if (tA->connAtoms.size() >=3 && tA->coordExist)
        {
            for (unsigned i=0; i < 3; i++)
            {
                vec1.push_back(tAtoms[tA->connAtoms[0]].coords[i]-tA->coords[i]);
                vec2.push_back(tAtoms[tA->connAtoms[1]].coords[i]-tA->coords[i]);
                vec3.push_back(tAtoms[tA->connAtoms[2]].coords[i]-tA->coords[i]);
            }
            REAL aVol=calChiralVol(vec1, vec2, vec3);
            if (fabs(aVol) < 0.001)
            {
                tP = true;
            }
        }
        return tP;
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
        //std::cout << "Check atom " << tAtoms[anchorIdx].id << std::endl;
        //std::cout << "It is in rings " << tAtoms[anchorIdx].inRings.size() << std::endl;
        
        if (tAtoms[anchorIdx].inRings.size() >1)
        {   
            std::vector<int> nShare1, nShare2;
            for (std::vector<int>::iterator iConn=tAtoms[anchorIdx].connAtoms.begin();
                     iConn!=tAtoms[anchorIdx].connAtoms.end(); iConn++)
            {
                nShare1.clear();
                // std::cout << "connected atom " << tAtoms[*iConn].id << std::endl;
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
            std::cout << "atom " << iA->id << std::endl;
            if (iA->chemType !="H")
            {
                int addH =(int)checkProtonateAll(iA, tAtoms, tBonds,  aPTab);
                std::cout << "Here 1" << std::endl;
                // REAL addH =checkProtonated(iA, tIdxMol);
                if (addH != 0)
                {
                    adjustHAtoms(tAtoms, tBonds, iA->seriNum, addH, tAddHIdx);
                }
                std::cout << "Here 2 " << std::endl;
            }
        } 
        
        std::cout << "Here 3" << std::endl;
        
        exit(1);
         
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
        
        // 2.
        
        
        
    }
    
}