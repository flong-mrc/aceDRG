
/* 
 * File:   ring.h
 * Author: flong
 *
 * Created on February 3, 2012, 10:26 AM
 */

#include "ring.h"

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
            rep(NullString)
    {
    }
    
    RingDict::RingDict(const RingDict& tR):isPlanar(tR.isPlanar),
            rep(tR.rep)
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
    
}