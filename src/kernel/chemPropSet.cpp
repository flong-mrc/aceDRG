
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
            if(!iA->chemType.compare("H")==0)
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
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(tAtoms[*iCA].bondingIdx == 2)
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
    }
    
    // Check protonation state of an atom and decide how many H atoms will added 
    // to bond it.
    
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
        
        //std::cout << "set H atom " << tIA->id  << " coords " << std::endl;
        
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
    
    extern void checkProtonatedCarBoxylicAcids(std::vector<AtomDict>::iterator tIA,
                                               std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               REAL tPH)
    {
        if (tIA->chemType=="C")
        {
            std::map<int, int> Oxys, Others;    // key: serial number, value: number of connection 
            for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                   iNB !=tIA->connAtoms.end(); iNB++)
            {
                if (tAtoms[*iNB].chemType =="O")
                {
                    Oxys[*iNB] = tAtoms[*iNB].connAtoms.size();
                }
                else
                {
                    Others[*iNB] = tAtoms[*iNB].connAtoms.size();
                }
                    
                // Now check possible de-protonated form:
                // Required: C connects to two O and one other atom
                if (Oxys.size() ==2 && Others.size()==1)
                {
                    REAL tPKa = 4.0;
                    
                    bool tOD=true;
                    bool lDeloc =false;
                    for (std::map<int, int>::iterator iNo=Oxys.begin();
                            iNo !=Oxys.end(); iNo++)
                    {
                        // Further check Oxygen atoms connect only to C
                        if (iNo->second !=1)
                        {
                            tOD = false;
                            break;
                        }
                        else
                        {
                            int aBIdx=getBond(tBonds, tIA->seriNum, iNo->first);
                            if (aBIdx !=-1)
                            {
                                if (tBonds[aBIdx].order.substr(0,4).compare("deloc")==0)
                                {
                                    lDeloc = true;
                                    break;
                                }
                            }
                        }
                    }
                    
                    std::map<int, int>::iterator iOther= Others.begin();
                    if(tAtoms[iOther->first].connAtoms.size()==1)
                    {
                        tOD=false;
                    }
                        
                    if (tOD)
                    {
                        if (lDeloc || (tPH > tPKa))
                        {
                            // If no pKa and PH values are given.
                            // we take the default action, by 
                            // assuming tPKa between 3-5 and tPH =7
                            
                        
                            unsigned i=0;
                            for (std::map<int, int>::iterator iNo=Oxys.begin();
                                iNo !=Oxys.end(); iNo++)
                            {
                                if (i==0)
                                {
                                    // change this one to a double bond O atom
                                    int aBIdx=getBond(tBonds, tIA->seriNum, iNo->first);
                                    if (aBIdx !=-1)
                                    {
                                        tBonds[aBIdx].order="2";
                                    }
                                    else
                                    {
                                        std::cout << "Can not find the bond between atom "
                                                  << tIA->id << " and " << tAtoms[iNo->first].id
                                                  << std::endl;
                                        exit(1);
                                    }
                                }
                                else if (i==1)
                                {
                                    // add a charge to this O atoms and set it as single bond.
                                    tAtoms[iNo->first].formalCharge = -1.0;
                                    int aBIdx=getBond(tBonds, tIA->seriNum, iNo->first);
                                    if (aBIdx !=-1)
                                    {
                                        tBonds[aBIdx].order="1";
                                    }
                                    else
                                    {
                                        std::cout << "Can not find the bond between atom "
                                                  << tIA->id << " and " << tAtoms[iNo->first].id
                                                  << std::endl;
                                        exit(1);
                                    }
                                }
                                i++;
                            }
                        }    
                    }
                    
                }
            }
            
        }
    }
    
    
    extern void checkProtonatedSulfuricAcids(std::vector<AtomDict>::iterator tIA,
                                               std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               REAL tPH)
    {
        // According to the empirical rules  from the examples (see the note,
        // assuming default pH value 7.0, pka value varies depending on connection)
        // (1) if S atom has two connections will not be de-protonated, except 
        // the non-H connected atom is on a aromatic ring.
        // (2) if S atom connects 3 atoms one of which is double bonded. 
        //     do nothing at the moment.
        // (3) if S atom connects to 4 atoms. Situation varies 
        
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
    
    extern void checkProtonatedNAcids(std::vector<AtomDict>::iterator tIA,
                                      std::vector<AtomDict> & tAtoms,
                                      std::vector<BondDict> & tBonds,
                                      REAL tPH)
    {
        if (tIA->connAtoms.size()==2)
        {
            // Make sure the formal charge of N atom is 0.0 so that
            // later on a H atom will added (protonation)
            tIA->formalCharge = 0.0;
        }
        else if (tIA->connAtoms.size() ==3)
        {
            std::vector<int> bondOtherSingle;
            for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                    iNB !=tIA->connAtoms.end(); iNB++)
            {
                std::string tStr;
                int aIdxB=getBond(tBonds, tIA->seriNum, *iNB);
                if (aIdxB !=-1)
                {
                    tStr.append(tBonds[aIdxB].order);
                    StrUpper(tStr);
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
                
                if (bondOtherSingle.size()==3)
                {
                    tIA->formalCharge =1.0;
                }
            }
        }
    }
    
    extern void checkProtonatedPAcids(std::vector<AtomDict>::iterator tIA,
                                      std::vector<AtomDict> & tAtoms,
                                      std::vector<BondDict> & tBonds,
                                      REAL tPH)
    {
        if (tIA->connAtoms.size()==4)
        {
            std::vector<int> Oxys, OxysS, bondOs, bondOtherSingle;
            
            for (std::vector<int>::iterator iOB=tIA->connAtoms.begin();
                    iOB!=tIA->connAtoms.end(); iOB++)
            {
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
    
    
    void extern setAromPlanes(std::vector<AtomDict> & tAtoms,
                              std::vector<RingDict> & tRings, 
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