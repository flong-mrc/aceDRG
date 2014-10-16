
/* 
 * File:   ring.h
 * Author: flong
 *
 * Created on February 3, 2012, 10:26 AM
 */

#include "ring.h"
#include "molecule.h"

namespace LIBMOL
{
    
    Ring::Ring()
    {
    }
    
    Ring::Ring(const Ring& tRing)
    {
        
        for (std::list<AtomDict>::const_iterator iA=tRing.atoms.begin();
                iA != tRing.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }
        
        for (std::list<BondDict>::const_iterator iB=tRing.bonds.begin();
                iB !=tRing.bonds.end(); iB++)
        {
            bonds.push_back((*iB));
        }
        
        for (std::vector<Torsion>::const_iterator iT=tRing.torsions.begin();
                iT !=tRing.torsions.end(); iT++)
        {
            torsions.push_back((*iT));
        }
    }
    
    Ring::~Ring()
    {
        // may not need those
        if(!atoms.empty())
        {
            atoms.clear();
        }
        if (!bonds.empty())
        {
            bonds.clear();
        }
        if (!torsions.empty())
        {
            torsions.clear();
        }
    }
    
    //Another ring class for the dictionary
    RingDict::RingDict():isPlanar(false),
            isAromatic(false),
            rep(NullString)
    {
    }
    
    RingDict::RingDict(const RingDict& tR):isPlanar(tR.isPlanar),
            isAromatic(tR.isAromatic),
            rep(tR.rep),
            sRep(tR.sRep)
    {
        for (std::vector<AtomDict>::const_iterator iA=tR.atoms.begin();
                iA != tR.atoms.end(); iA++)
        {
            atoms.push_back(*iA);
        }
        
        for (std::map<int, std::map<ID, int> >::const_iterator iL=tR.atomsLink.begin();
                iL !=tR.atomsLink.end(); iL++)
        {
            for (std::map<ID, int>::const_iterator iI=iL->second.begin();
                    iI != iL->second.end(); iI++)
            {
                atomsLink[iL->first][iI->first]=iI->second;
            }
        }
    }
    
    RingDict::RingDict(const std::vector<AtomDict> tAtoms)
    {
        isPlanar = true;
        
        for (std::vector<AtomDict>::const_iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            atoms.push_back(*iA);
            if (iA->bondingIdx !=2)
            {
                isPlanar = false;
            }
        }
    }
    
    RingDict::~RingDict()
    {
    }
    
    void RingDict::setAtmsLink(std::vector<AtomDict> tAllAtoms)
    {
        atomsLink.clear();
        
        std::vector<int> tIdxs;
        int tSize = (int)tAllAtoms.size();
        int tRSize = (int)atoms.size();
        // std::cout << "Ring Size " << tRSize << std::endl;
        
        for (int i=0; i < tRSize; i++)
        {
            int iPos = atoms[i].atomPosition(tAllAtoms);
            if (iPos >=0 && iPos < tSize)
            {
                tIdxs.push_back(iPos);
            }
            else
            {
                std::cout << "Can not find the serial number of atom " 
                          << atoms[i].id << std::endl;
                exit(1);
            }
        }
        
        // std::cout << "Index set has " << (int)tIdxs.size() << " elements" << std::endl;
        
        // loop over atoms
        int iCur =atoms[0].seriNum;
        int iPre = -1;
        int iLoop = 0;
        while (iLoop < tRSize)
        {
            int i1=0;
            for (std::vector<int>::iterator iC=tAllAtoms[iCur].connAtoms.begin();
                    iC != tAllAtoms[iCur].connAtoms.end(); iC++)
            { 
                if (iLoop==0)
                {
                    if (std::find(tIdxs.begin(), tIdxs.end(), *iC) !=tIdxs.end())
                    {
                        if (i1==0)
                        {
                            // clockwise
                            atomsLink[iCur]["c"]  =*iC;
                            atomsLink[*iC]["ac"]  = iCur;
                            i1++;
                        }
                        else if (i1==1)
                        {
                            // anti-clockwise
                            atomsLink[iCur]["ac"]=*iC;
                            atomsLink[*iC]["c"] = iCur;
                            break;
                        }
                    }
                }
                else 
                {
                    if (std::find(tIdxs.begin(), tIdxs.end(), *iC) !=tIdxs.end()
                        && *iC != iPre)
                    {
                        atomsLink[iCur]["c"]  = *iC;
                        atomsLink[iCur]["ac"] = iPre;
                        break;
                    }
                        
                }
            }
            
            iPre = iCur;
            iCur = atomsLink[iCur]["c"];
            iLoop++;
        }
        
        // Check
        // std::cout << "The ring connection are " << std::endl;
        
        // std::cout << "Clockwise : " << std::endl;
        /*
        int i =0;
        std::map<int, std::map<ID, int> >::iterator iA=atomsLink.begin();
        iCur = iA->first;
        int iLink = iA->second["c"];
        do
        {
            // std::cout << "Atom " << tAllAtoms[iCur].id << " is linked to " 
            //           <<  tAllAtoms[iLink].id << std::endl;
            iCur  = iLink;
            iLink = atomsLink[iCur]["c"];
            i++;        
        }while(i < tRSize);
        */
    }
    
    void RingDict::setPlaneProp()
    {
        bool lSP2=true;
        
        for (std::vector<AtomDict>::iterator iA=atoms.begin();
                iA !=atoms.end(); iA++)
        {
            if (iA->bondingIdx !=2)
            {
                lSP2 =false;
                break;
            }
        }
        
        if (lSP2)
        {
            isPlanar = lSP2;
        }
    }
    
    
    extern bool AtomsInSameRing(AtomDict & tAt1, AtomDict tAt2, 
                std::vector<RingDict> & tRings)
    {
        bool tIn = false;
        
        for (std::vector<RingDict>::iterator iRi=tRings.begin();
                iRi !=tRings.end(); iRi++)
        {
            bool tF1= false, tF2=false;
            for (std::vector<AtomDict>::iterator iAt=iRi->atoms.begin();
                    iAt != iRi->atoms.end(); iAt++)
            {
                if (iAt->id==tAt1.id)
                {
                    tF1 = true;
                }
                else if (iAt->id == tAt2.id)
                {
                    tF2 = true;
                }
            }
            if (tF1 && tF2)
            {
                tIn = true;
                break;
            }
        }
        
        return tIn;
    }
    
    
    
    extern void buildOneRing(std::vector<AtomDict> & tAtoms, 
                             std::vector<AtomDict> sSet, 
                             std::vector<AtomDict> & doneSet)
    {
        
        
        
    }
    
    extern void mergePlaneRings(std::vector<RingDict> & tAllRings,
                                std::vector<std::vector<int> > & tRingAtoms)
    {
        
        std::cout << "Number of rings " << tAllRings.size() << std::endl;
      
        
        // Connected planar rings
        std::vector<int> DoneList;
        std::vector<std::vector<int> > mergedRings;       // The inside vector represents a set of merged rings 
        
        for (unsigned i=0; i < tAllRings.size(); i++)
        {
            if (tAllRings[i].isPlanar && std::find(DoneList.begin(), DoneList.end(), i) ==DoneList.end())
            {
                // std::cout << "Planar ring " << i << std::endl;
                
                DoneList.push_back(i);
                std::vector<int> curLinkedRing;
                curLinkedRing.push_back(i);
                findAllRingsConnectedOneRing(i, tAllRings,DoneList, curLinkedRing); 
                if (curLinkedRing.size() > 1)
                {
                    mergedRings.push_back(curLinkedRing);
                }
            }
            
        }
        
        // 
        for (std::vector<std::vector<int> >::iterator iMr =mergedRings.begin();
                iMr !=mergedRings.end(); iMr++)
        {
            std::vector<int> atmIdxs;
            for (std::vector<int>::iterator iR=iMr->begin();
                    iR !=iMr->end(); iR++)
            {
                
                for (std::vector<AtomDict>::iterator iAt=tAllRings[*iR].atoms.begin();
                        iAt !=tAllRings[*iR].atoms.end(); iAt++)
                {
                    if(std::find(atmIdxs.begin(), atmIdxs.end(), iAt->seriNum)
                            ==atmIdxs.end())
                    {
                        atmIdxs.push_back(iAt->seriNum);
                    }
                }
            }
            
            tRingAtoms.push_back(atmIdxs);
        }
        
        //
        
        
        /*
        
        if (mergedRings.size() !=0)
        {
            std::cout << "Merged planar rings : " << std::endl;
            for (std::vector<std::vector<int> >::iterator iMr=mergedRings.begin();
                    iMr != mergedRings.end(); iMr++)
            {
                std::cout << "A merged system contains " << iMr->size() 
                          << " rings. They are:  " << std::endl;
                for (std::vector<int>::iterator iR=iMr->begin();
                        iR != iMr->end(); iR++)
                {
                    std::cout << "The Ring with atoms " << std::endl;
                    for (std::vector<AtomDict>::iterator iA=tAllRings[*iR].atoms.begin();
                            iA != tAllRings[*iR].atoms.end(); iA++)
                    {
                        std::cout << iA->id << std::endl;
                    }
                }
            }
            
        }
        
        */
        
        
    }
    
    extern void mergePlaneRings(std::vector<RingDict>          & tAllRings,
                                std::vector<std::vector<int> > & tRingAtoms,
                                std::vector<AtomDict >         & tAtoms)
    {
        
        // std::cout << "Number of rings " << tAllRings.size() << std::endl;
      
        
        // Connected planar rings
        
        std::vector<int> DoneList;
        std::vector<std::vector<int> > mergedRings;       // The inside vector represents a set of merged rings 
        
        for (unsigned i=0; i < tAllRings.size(); i++)
        {
            if (tAllRings[i].isPlanar && std::find(DoneList.begin(), DoneList.end(), i) ==DoneList.end())
            {
                
                // std::cout << "Planar ring " << i << std::endl;
                
                DoneList.push_back(i);
                std::vector<int> curLinkedRing;
                curLinkedRing.push_back(i);
                findAllRingsConnectedOneRing(i, tAllRings,DoneList, curLinkedRing); 
                if (curLinkedRing.size() > 1)
                {
                    mergedRings.push_back(curLinkedRing);
                }
            }
            
        }
        
        
        
        if (mergedRings.size() !=0)
        {
            std::cout << "Merged planar rings : " << std::endl;
            for (std::vector<std::vector<int> >::iterator iMr=mergedRings.begin();
                    iMr != mergedRings.end(); iMr++)
            {
                std::cout << "A merged system contains " << iMr->size() 
                          << " rings. They are:  " << std::endl;
                for (std::vector<int>::iterator iR=iMr->begin();
                        iR != iMr->end(); iR++)
                {
                    std::cout << "The Ring with atoms " << std::endl;
                    for (std::vector<AtomDict>::iterator iA=tAllRings[*iR].atoms.begin();
                            iA != tAllRings[*iR].atoms.end(); iA++)
                    {
                        std::cout << iA->id << std::endl;
                    }
                }
            }
            
        }
        
        
        // 
        for (std::vector<std::vector<int> >::iterator iMr =mergedRings.begin();
                iMr !=mergedRings.end(); iMr++)
        {
            
            //std::cout << "A merged system contains " << iMr->size() 
            //          << " rings.  " << std::endl;
            
            std::vector<int> atmIdxs;
            for (std::vector<int>::iterator iR=iMr->begin();
                    iR !=iMr->end(); iR++)
            {
                for (std::vector<AtomDict>::iterator iAt=tAllRings[*iR].atoms.begin();
                        iAt !=tAllRings[*iR].atoms.end(); iAt++)
                {
                    
                    if(std::find(atmIdxs.begin(), atmIdxs.end(), iAt->seriNum)
                            ==atmIdxs.end())
                    {
                        atmIdxs.push_back(iAt->seriNum);
                    }
                }
            }
            
            std::cout<< "It contains " << atmIdxs.size() << " atoms " << std::endl;
                   
            if(checkAromaSys(atmIdxs, tAtoms))
            {
                std::cout << "It is an aromatic ring" << std::endl;
                tRingAtoms.push_back(atmIdxs);
            }
            else
            {
                std::cout << "It is not aromatic ring " << std::endl;
            }
        }
        
        
        /*
        std::cout << "Number of merged atom rings " 
                  << tRingAtoms.size() << std::endl;
        
        if (tRingAtoms.size() >0)
        {
            for (std::vector<std::vector<int> >::iterator iR=tRingAtoms.begin();
                   iR !=tRingAtoms.end(); iR++)
            {
                std::cout << "Atoms in a merged ring : " << std::endl;
                for (std::vector<int>::iterator iAt=iR->begin();
                        iAt !=iR->end(); iAt++)
                {
                    std::cout << tAtoms[*iAt].id << std::endl;
                }
            }
        }
        */
      
        
        
        
    }
    
    
    
    // A recursive function that get a set of planar rings sharing edges 
    // with each other.
    
    extern void findAllRingsConnectedOneRing(int tCurIdx, 
                                             std::vector<RingDict> & tRings,
                                             std::vector<int>      & tDoneList, 
                                             std::vector<int>      & tCurLinkedRing)
    
    {
        
        std::vector<int> tSeri;
        for (std::vector<AtomDict>::iterator iA=tRings[tCurIdx].atoms.begin();
                iA != tRings[tCurIdx].atoms.end(); iA++)
        {
            tSeri.push_back(iA->seriNum);
        }
        
        for (unsigned tNext=tCurIdx+1; tNext < tRings.size(); tNext++)
        {
            if (std::find(tDoneList.begin(), tDoneList.end(), tNext) == tDoneList.end()
                    && tRings[tNext].isPlanar)
            {
                int aA=0;
            
                for (std::vector<AtomDict>::iterator iNA=tRings[tNext].atoms.begin();
                        iNA != tRings[tNext].atoms.end(); iNA++)
                {
                    if (std::find(tSeri.begin(), tSeri.end(), iNA->seriNum) 
                          !=tSeri.end())
                    {
                        aA++;
                    }
                }
                
                if (aA > 1)
                {
                    tCurLinkedRing.push_back(tNext);
                    tDoneList.push_back(tNext);
                    findAllRingsConnectedOneRing(tNext, tRings, tDoneList, tCurLinkedRing);
                }
            }
        }
        
    }
    
    
    extern bool checkAromaSys(std::vector<int>      & tSubAtoms,
                              std::vector<AtomDict> & tAtoms)
    {
        bool lAr=false;
        REAL numPiAll=0.0;
        
        for (std::vector<int>::iterator iAt=tSubAtoms.begin();
                iAt !=tSubAtoms.end(); iAt++)
        {
            REAL numOneAtm = setPiForOneAtom(*iAt, tAtoms);
            std::cout << "Atom " << tAtoms[*iAt].id 
                      << " add " << numOneAtm << " pi atoms" << std::endl;
            if (numOneAtm > 0.0)
            {
                numPiAll +=numOneAtm;
            }
        }       
        
        std::cout << "number of pi electron in the system is " 
                  <<  numPiAll << std::endl;
        
        if (numPiAll > 0.0 && fabs(fmod(numPiAll, 4.0)-2.0) < 0.001)
        {
            lAr = true;
        }
        
        return lAr;
        
    }
    
    extern REAL setPiForOneAtom(int tIdx, std::vector<AtomDict> & tAtoms)
    {
        REAL aN=0.0;
        
        if (tAtoms[tIdx].bondingIdx ==2 && tAtoms[tIdx].formalCharge==0.0)
        {
            if (tAtoms[tIdx].chemType.compare("C") ==0)
            {           
                if (tAtoms[tIdx].connAtoms.size() ==3)
                {
                    bool aD=false;
                    for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size();
                           i++)
                    {
                        // check if there is a double bonded atom outside the ring 
                        int aNbIdx = tAtoms[tIdx].connAtoms[i];
                        if (tAtoms[aNbIdx].chemType.compare("O")==0 
                              && tAtoms[aNbIdx].connAtoms.size()==1)
                        {
                            // double bonded Oxy 
                            aD = true;
                        }
                        else if (tAtoms[aNbIdx].chemType.compare("C")==0
                                 && tAtoms[aNbIdx].connAtoms.size()==3)
                        {
                            int nH=0;
                            for (unsigned j=0; j < tAtoms[aNbIdx].connAtoms.size(); j++)
                            {
                                int bNbIdx = tAtoms[aNbIdx].connAtoms[j];
                                if (tAtoms[bNbIdx].chemType.compare("H")==0)
                                {
                                    nH++;
                                }
                            }
                            if (nH >=2)
                            {
                                aD=true;
                            }
                        }
                    }
                    if (!aD)
                    {
                        aN=1.0;
                    }
                }
                else if (tAtoms[tIdx].connAtoms.size() ==2)
                {
                    // Place holder in case for future. 
                    aN=0.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("N") ==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==2)
                {
                    // one electron becomes a pi electron.
                    aN=1.0;
                }
                else if (tAtoms[tIdx].connAtoms.size()==3)
                {
                    // the lone pair contributes two
                    aN=2.0;
                }
            }
            else if(tAtoms[tIdx].chemType.compare("O") ==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==2)
                {
                    // one lone pair contributes two
                    aN=2.0;
                }
            }
            else if(tAtoms[tIdx].chemType.compare("S") ==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==2)
                {
                    // one lone pair contributes two
                    aN=2.0;
                }
            }
            else if(tAtoms[tIdx].chemType.compare("P") ==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==3)
                {
                    // one lone pair contributes two
                    aN=2.0;
                }
                else if (tAtoms[tIdx].connAtoms.size()==2)
                {
                    // The lone pair contributes none and occupies a sp2 orbit. 
                    // Two other valence electrons occupy other two sp2 orbits.  
                    aN=1.0;
                }
            }
        }
        
        return aN;
    }
        
    extern void checkAndSetupPlanes(std::vector<RingDict>  & tAllRings,
                                    std::vector<PlaneDict> & tPlanes,
                                    std::vector<AtomDict>  & tAtoms)
    {
        
        // Check aromaticity for individual rings
        std::vector<RingDict> tAroRings;
        std::cout << "\nCheck aromaticity for individual ring " << std::endl;
                for (std::vector<RingDict>::iterator iR=tAllRings.begin();
                iR !=tAllRings.end(); iR++)
        {
            if (iR->isPlanar)
            {
                std::cout << "One Planar Ring! It contains " << iR->atoms.size()
                          << " atoms. They are : " << std::endl;
                std::vector<int> atmIdx;
                for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
                {
                    atmIdx.push_back(iAt->seriNum);
                    std::cout << iAt->id << std::endl;
                }        
                
                bool lAro = checkAromaSys(atmIdx, tAtoms);
                if (lAro)
                {
                    iR->isAromatic = true;
                    tAroRings.push_back(*iR);
                    std::cout << "It is an aromatic ring " << std::endl;
                }
                else
                {
                    iR->isAromatic = false;
                    std::cout << "It is not an aromatic ring " << std::endl;
                }
                
                std::cout << std::endl;
            }
            else
            {
                iR->isAromatic =false;
            }
            
        }
        
        
        std::vector<std::vector<int> >  ringAtomSets;
        mergePlaneRings(tAllRings, ringAtomSets, tAtoms);
        
        tPlanes.clear();
        
        setAllRingPlanes(tAllRings, tAtoms, tPlanes);
        
        // std::cout << "Check other planes " << std::endl;
        setAllOtherPlanes(tAllRings, tAtoms, tPlanes);
        
        std::cout << "There are " << tPlanes.size() 
                  << " planes. They are " << std::endl;
        
        for (std::vector<PlaneDict>::iterator iPl= tPlanes.begin();
                iPl !=tPlanes.end(); iPl++)
        {
            std::cout << "One Plane with " << iPl->atoms.size() 
                      << " atoms " << std::endl;
            for (std::map<ID, int>::iterator iAt=iPl->atoms.begin();
                    iAt !=iPl->atoms.end(); iAt++)
            {
                std::cout << iAt->first << std::endl;
            }
            std::cout << std::endl;
        }
        
        
    }
    
    extern void TestAromaticity(std::vector<RingDict>  & tAllRings,
                                std::vector<AtomDict>  & tAtoms)
    {
        
        std::cout << "\nCheck aromaticity for individual ring " << std::endl;
        
        for (std::vector<RingDict>::iterator iR=tAllRings.begin();
                iR !=tAllRings.end(); iR++)
        {
            if (iR->isPlanar)
            {
                std::cout << "One Planar Ring! It contains " << iR->atoms.size()
                          << " atoms. They are : " << std::endl;
                std::vector<int> atmIdx;
                for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
                {
                    atmIdx.push_back(iAt->seriNum);
                    std::cout << iAt->id << std::endl;
                }        
                
                bool lAro = checkAromaSys(atmIdx, tAtoms);
                if (lAro)
                {
                    std::cout << "It is an aromatic ring " << std::endl;
                }
                else
                {
                    std::cout << "It is not an aromatic ring " << std::endl;
                }
                
                std::cout << std::endl;
            }
            else
            {
                iR->isAromatic = false;
            }
            
        }
        
        
        // Now, test merged ring systems
        // std::cout << "Now Merged ring systems. " << std::endl;
        
        // std::vector<PlaneDict> mergedPlanes;
        
        // checkAndSetupPlanes(tAllRings, mergedPlanes, tAtoms);
        
    }
    
    extern void setAromaticBonds(std::vector<RingDict>  & tRings,
                                 std::vector<BondDict>  & tBonds)
    {
        for (std::vector<RingDict>::iterator iR=tRings.begin();
                iR !=tRings.end(); iR++)
        {
            if (iR->isAromatic)
            {
                for (unsigned i=0; i < iR->atoms.size(); i++)
                {
                    for (unsigned j=i+1; j < iR->atoms.size(); j++)
                    {
                        if (std::find(iR->atoms[i].connAtoms.begin(), 
                                  iR->atoms[i].connAtoms.end(), iR->atoms[j].seriNum)
                             !=iR->atoms[i].connAtoms.end())
                        {
                            int idxB=getBond(tBonds, iR->atoms[i].seriNum, iR->atoms[j].seriNum);
                            if (idxB >=0)
                            {
                                tBonds[idxB].isAromatic = true;
                                tBonds[idxB].order      = "aromatic";
                                std::cout << "Bond between atom "  
                                          << iR->atoms[i].id << " and "
                                          << iR->atoms[j].id << " is aromatic " 
                                          << std::endl;
                            }
                            else
                            {
                                std::cout << "Bug. can not find the bond between " 
                                          << iR->atoms[i].id << " and "
                                          << iR->atoms[j].id << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    
}