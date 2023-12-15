
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
            isAromaticP(false),
            isAntiAroma(false),
            isSugar(false),
            sugarType(NullString),
            rep(NullString),
            sRep(NullString),
            extId(NullString)
    {
    }

    RingDict::RingDict(const RingDict& tR):isPlanar(tR.isPlanar),
            isAromatic(tR.isAromatic),
            isAromaticP(tR.isAromaticP),
            isAntiAroma(tR.isAntiAroma),
            isSugar(tR.isSugar),
            sugarType(tR.sugarType),
            rep(tR.rep),
            sRep(tR.sRep),
            extId(tR.extId)
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

        for (std::map<int, std::vector<int> >::const_iterator iL=tR.ringAtomLink.begin();
                iL !=tR.ringAtomLink.end(); iL++)
        {
            for (std::vector<int>::const_iterator iAt=iL->second.begin();
                    iAt !=iL->second.end(); iAt++)
            {
                ringAtomLink[iL->first].push_back(*iAt);
            }
        }

        for (std::map<ID, REAL>::const_iterator iSR=tR.sugarTors.begin();
                iSR !=tR.sugarTors.end(); iSR++)
        {
            sugarTors[iSR->first] = iSR->second;
        }

        for (std::vector<int>::const_iterator iBo=tR.bondIdxs.begin();
                iBo != tR.bondIdxs.end(); iBo++)
        {
            bondIdxs.push_back(*iBo);
        }

        for (std::vector<AtomDict>::const_iterator
             iSeqAtm=tR.seqedAtoms.begin(); iSeqAtm!=tR.seqedAtoms.end();
             iSeqAtm++)
        {
            seqedAtoms.push_back(*iSeqAtm);
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
        isAromatic  = false;
        isAromaticP = false;

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
        //

        //std::cout << "The ring connection are " << std::endl;

        //std::cout << "Clockwise : " << std::endl;

        int i =0;
        std::map<int, std::map<ID, int> >::iterator iA=atomsLink.begin();
        iCur = iA->first;
        int iLink = iA->second["c"];
        do
        {
            //std::cout << "Atom " << tAllAtoms[iCur].id << " is linked to "
            //           <<  tAllAtoms[iLink].id << std::endl;
            iCur  = iLink;
            iLink = atomsLink[iCur]["c"];
            i++;
        }while(i < tRSize);

    }

    void RingDict::setRingAtmsLinks()
    {

        if (atoms.size() > 0)
        {
            std::vector<int> tSerNums;
            ringAtomLink.clear();
            for (std::vector<AtomDict>::iterator iAt =atoms.begin();
                    iAt != atoms.end(); iAt++)
            {
                tSerNums.push_back(iAt->seriNum);

            }

            for (std::vector<AtomDict>::iterator iAt =atoms.begin();
                    iAt != atoms.end(); iAt++)
            {
                // std::cout << "Atom " << iAt->id << " serN " << iAt->seriNum << std::endl;
                for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                        iNB !=iAt->connAtoms.end(); iNB++)
                {
                    if (std::find(tSerNums.begin(), tSerNums.end(), *iNB) !=tSerNums.end())
                    {
                        // std::cout << "Atom " << *iNB
                        //          << " linked to " << iAt->seriNum << std::endl;
                        ringAtomLink[iAt->seriNum].push_back(*iNB);
                    }
                }
            }

            // Check the ring links
            for (std::map<int, std::vector<int> >::iterator iIdx=ringAtomLink.begin();
                   iIdx !=ringAtomLink.end(); iIdx++)
            {
                if (iIdx->second.size() !=2)
                {
                    std::cout << "ERROR: Atom " << iIdx->first
                              << " does not have 2 connection in ring "
                              << rep << std::endl;
                    exit(1);
                }
            }

        }

    }

    void RingDict::setPlaneProp()
    {
        bool lSP2=true;

        for (std::vector<AtomDict>::iterator iA=atoms.begin();
                iA !=atoms.end(); iA++)
        {
            if (iA->bondingIdx !=2 && iA->chemType.compare("N") !=0)
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

    void RingDict::setBondIdxs(std::vector<BondDict> & tBonds, int & tStartIdx)
    {
        bondIdxs.clear();
        int nAll =0;
        for (std::map<int, std::vector<int> >::iterator iL=ringAtomLink.begin();
                iL !=ringAtomLink.end(); iL++)
        {
            int nSet = 0;
            for (std::vector<int>::iterator iB=iL->second.begin();
                    iB != iL->second.end(); iB++)
            {
                int idxB = getBond(tBonds, iL->first, *iB);

                if (idxB !=-1)
                {
                    if (tBonds[idxB].order.find("AR")==tBonds[idxB].order.npos)
                    {
                        nSet+=1;
                    }
                    if (std::find(bondIdxs.begin(), bondIdxs.end(), idxB)
                            ==bondIdxs.end())
                    {
                        std::cout << "Bond " << idxB << std::endl;
                        std::cout << "Two atoms: " << tBonds[idxB].atoms[0]
                                  << " and " << tBonds[idxB].atoms[1] << std::endl;
                        std::cout << "Bond order " << tBonds[idxB].order << std::endl;
                        // std::cout << "Bond " << idxB << " is included " << std::endl;
                        bondIdxs.push_back(idxB);
                    }
                }
                if (nSet ==1 && tStartIdx ==0)
                {
                    tStartIdx = iL->first;
                }
                nAll+=nSet;
            }

        }
        if (nAll==bondIdxs.size())
        {
            tStartIdx = -1;
        }

        // check
        if (bondIdxs.size() != atoms.size())
        {

            std::cout << "Bug: number of bonds in the ring does not equal to "
                      << "number of atoms in the ring. " << std::endl;
            exit(1);
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

    extern int  checkIfAngleInSameRing(std::vector<AtomDict> & tAtoms,
                                       std::vector<RingDict> & tRings,
                                       int tIdxCen, int tIdx1, int tIdx2)
    {


        int nInRing = 0;
        std::vector<std::string> r1, r2;
        for (std::map<std::string, int>::iterator iM1=tAtoms[tIdx1].ringRep.begin();
                iM1 !=tAtoms[tIdx1].ringRep.end(); iM1++)
        {
            r1.push_back(iM1->first);
        }
        for (std::map<std::string, int>::iterator iM2=tAtoms[tIdx2].ringRep.begin();
                iM2 !=tAtoms[tIdx2].ringRep.end(); iM2++)
        {
            r2.push_back(iM2->first);
        }
        for (std::map<std::string, int>::iterator iR=tAtoms[tIdxCen].ringRep.begin();
                iR !=tAtoms[tIdxCen].ringRep.end(); iR++)
        {
            if (std::find(r1.begin(), r1.end(), iR->first) != r1.end()
                &&
                std::find(r2.begin(), r2.end(), iR->first) !=r2.end())
            {
                nInRing = iR->second;
                //std::cout << "Check atoms in an angle : " << std::endl;
                //std::cout << "center atom " << tAtoms[tIdxCen].id <<std::endl;
                //std::cout << "atom2 " << tAtoms[tIdx1].id <<std::endl;
                //std::cout << "atom3 " << tAtoms[tIdx2].id <<std::endl;
                //std::cout << "3 atoms are in ring " << iR->first
                //          << " of size " << iR->second <<  std::endl;
                break;
            }
        }

        return nInRing;
    }

    extern void buildOneRing(std::vector<AtomDict> & tAtoms,
                             std::vector<AtomDict> sSet,
                             std::vector<AtomDict> & doneSet)
    {



    }

    extern void mergePlaneRings(std::vector<RingDict> & tAllRings,
                                std::vector<std::vector<int> > & tRingAtoms)
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

        /*
        for (unsigned i=0; i < tAllRings.size(); i++)
        {
            std::cout << "Is ring " << i+1 << " planar ring ? "
                      << tAllRings[i].isPlanar << std::endl;
        }
        */
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

            std::cout << "number of Merged planar rings : "
                      << mergedRings.size() << std::endl;
            for (std::vector<std::vector<int> >::iterator iMr=mergedRings.begin();
                    iMr != mergedRings.end(); iMr++)
            {
                std::vector<int> atmIdxs;
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
                        if(std::find(atmIdxs.begin(), atmIdxs.end(), iA->seriNum)
                            ==atmIdxs.end())
                        {
                            atmIdxs.push_back(iA->seriNum);
                        }
                    }
                }

                //std::cout << "Overall, The merged system contains "
                //          << atmIdxs.size() << " atoms " << std::endl;


                if(checkAromaSys(atmIdxs, tAtoms))
                {
                    std::cout << "It is an aromatic system" << std::endl;
                    tRingAtoms.push_back(*iMr);
                    for (std::vector<int>::iterator iR=iMr->begin();
                        iR != iMr->end(); iR++)
                    {
                        tAllRings[*iR].isAromatic = true;
                    }
                }
                else
                {
                    std::cout << "It is not aromatic system for the merged system"
                              << std::endl;
                    // Further check
                    std::vector<int> atmIdxsOR;
                    for (std::vector<int>::iterator iR=iMr->begin();
                        iR != iMr->end(); iR++)
                    {
                        if (!tAllRings[*iR].isAromatic)
                        {
                           for (std::vector<AtomDict>::iterator iRAt=tAllRings[*iR].atoms.begin();
                                   iRAt !=tAllRings[*iR].atoms.end(); iRAt++)
                           {
                               atmIdxsOR.push_back(iRAt->seriNum);
                           }
                           if (checkAromaSys(atmIdxsOR, tAtoms))
                           {
                               tAllRings[*iR].isAromatic = true;
                           }
                        }
                    }


                }
            }

        }

        //

        /*
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
            REAL numOneAtm =0.0;
            //if(tAtoms[*iAt].chemType=="S"
            //   && tAtoms[*iAt].connAtoms.size()==4)
            //{
            //    numOneAtm = setPiForOne_S_Sp3_Atom(*iAt, tSubAtoms, tAtoms);
            //}
            //else
            //{


            numOneAtm = setPiForOneAtom(*iAt, tAtoms);
            //}


            std::cout << "Atom " << tAtoms[*iAt].id
                      << " its charge " << tAtoms[*iAt].charge << std::endl
                      << " add " << numOneAtm << " pi atoms" << std::endl;
            if (numOneAtm > 0.0)
            {
                numPiAll +=numOneAtm;
            }
        }

        //std::cout << "number of pi electron in the system is "
        //          <<  numPiAll << std::endl;

        if (numPiAll > 0.0 && fabs(fmod(numPiAll, 4.0)-2.0) < 0.001)
        {
            lAr = true;
        }

        return lAr;

    }

    extern bool checkAromaSys(std::vector<int>      & tSubAtoms,
                              std::vector<AtomDict> & tAtoms,
                              int                     tMode)
    {
        bool lAr=false;
        REAL numPiAll1=0.0;
        REAL numPiAll2=0.0;
        for (std::vector<int>::iterator iAt=tSubAtoms.begin();
                iAt !=tSubAtoms.end(); iAt++)
        {
            REAL numOneAtm1 =0.0;
            REAL numOneAtm2 =0.0;
            //if(tAtoms[*iAt].chemType=="S"
            //   && tAtoms[*iAt].connAtoms.size()==4)
            //{
            //    numOneAtm = setPiForOne_S_Sp3_Atom(*iAt, tSubAtoms, tAtoms);
            //}
            //else
            //{

            std::cout << "tMode " << tMode << std::endl;
            numOneAtm1 = setPiForOneAtomNoMetal(*iAt, tAtoms, tMode);
            numOneAtm2 = setPiForOneAtomAll(*iAt, tAtoms, tMode);
            //}

            std::cout << "XXX here Atom " << tAtoms[*iAt].id
                  << " its charge " << tAtoms[*iAt].charge << std::endl
                  << " add " << numOneAtm1 << " or "
                  <<  numOneAtm2 << " pi atoms" << std::endl;

            if (numOneAtm1 > 0.0)
            {
                numPiAll1 +=numOneAtm1;
            }
            if (numOneAtm2 > 0.0)
            {
                numPiAll2 +=numOneAtm2;
            }
        }

        //std::cout << "number of pi electron in the system is "
        //          <<  numPiAll << std::endl;

        if (numPiAll1 > 0.0 && fabs(fmod(numPiAll1, 4.0)-2.0) < 0.001)
        {
            lAr = true;
        }
        else if (numPiAll2 > 0.0 && fabs(fmod(numPiAll2, 4.0)-2.0) < 0.001)
        {
            lAr = true;
        }

        return lAr;

    }

    extern REAL getTotalPiElec(std::vector<int>      & tSubAtoms,
                              std::vector<AtomDict> & tAtoms,
                              int                     tMode)
    {

        REAL numPiAll1=0.0;
        REAL numPiAll2=0.0;
        REAL numPiAll =0.0;
        for (std::vector<int>::iterator iAt=tSubAtoms.begin();
                iAt !=tSubAtoms.end(); iAt++)
        {
            REAL numOneAtm1 =0.0;
            REAL numOneAtm2 =0.0;
            //if(tAtoms[*iAt].chemType=="S"
            //   && tAtoms[*iAt].connAtoms.size()==4)
            //{
            //    numOneAtm = setPiForOne_S_Sp3_Atom(*iAt, tSubAtoms, tAtoms);
            //}
            //else
            //{


            numOneAtm1= setPiForOneAtomNoMetal(*iAt, tAtoms, tMode);
            numOneAtm2 = setPiForOneAtomAll(*iAt, tAtoms, tMode);
            //}

            std::cout << "XXX here Atom " << tAtoms[*iAt].id
                  << " its charge " << tAtoms[*iAt].charge << std::endl
                  << " add " << numOneAtm1 << " or "
                  <<  numOneAtm2 << " pi atoms" << std::endl;

            if (numOneAtm1 > 0.0)
            {
                numPiAll1 +=numOneAtm1;
            }
            if (numOneAtm2 > 0.0)
            {
                numPiAll2 +=numOneAtm2;
            }
        }


        if (numPiAll2 > 0.0 && fabs(fmod(numPiAll2, 4.0)-2.0) < 0.001)
        {
            numPiAll= numPiAll2;
        }
        else
        {
            numPiAll= numPiAll1;
        }
        std::cout << "number of pi electron in the system is "
                  <<  numPiAll << std::endl;

        return numPiAll;

    }


    extern bool checkAromaSys(std::vector<int>      & tSubAtoms,
                              std::vector<AtomDict> & tAtoms,
                              std::vector<BondDict> & tBonds)
    {
        bool lAr=false;
        REAL numPiAll=0.0;

        for (std::vector<int>::iterator iAt=tSubAtoms.begin();
                iAt !=tSubAtoms.end(); iAt++)
        {
            REAL numOneAtm =0.0;
            //if(tAtoms[*iAt].chemType=="S"
            //   && tAtoms[*iAt].connAtoms.size()==4)
            //{
            //    numOneAtm = setPiForOne_S_Sp3_Atom(*iAt, tSubAtoms, tAtoms);
            //}
            //else
            //{
            std::cout << "2 Atom " << tAtoms[*iAt].id
                      << " its charge " << tAtoms[*iAt].charge << std::endl
                      << " add " << numOneAtm << " pi atoms" << std::endl;

            numOneAtm = setPiForOneAtom(*iAt, tAtoms, tBonds);
            //}


            std::cout << "3 Atom " << tAtoms[*iAt].id
                      << " its charge " << tAtoms[*iAt].charge << std::endl
                      << " add " << numOneAtm << " pi atoms" << std::endl;
            if (numOneAtm > 0.0)
            {
                numPiAll +=numOneAtm;
            }
        }

        //std::cout << "number of pi electron in the system is "
        //          <<  numPiAll << std::endl;

        if (numPiAll > 0.0 && fabs(fmod(numPiAll, 4.0)-2.0) < 0.001)
        {
            lAr = true;
        }

        return lAr;

    }


    extern REAL setPiForOneAtom(int tIdx, std::vector<AtomDict> & tAtoms)
    {
        REAL aN=0.0;
        std::cout << "Here charge is " << tAtoms[tIdx].charge << std::endl;
        std::cout << "formalCharge is " << tAtoms[tIdx].formalCharge
                  << std::endl;
        if (tAtoms[tIdx].bondingIdx ==2)
        {

            if (tAtoms[tIdx].charge==0.0)
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==3)
                    {
                        bool aD=false;
                        for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size();
                           i++)
                        {
                            // check if there is a double bonded atom outside
                            // the ring
                            int aNbIdx = tAtoms[tIdx].connAtoms[i];
                            if (tAtoms[aNbIdx].chemType.compare("O")==0
                                && tAtoms[aNbIdx].connAtoms.size()==1
                                && tAtoms[aNbIdx].charge ==0)
                            {
                                // double bonded Oxy
                                aD = true;
                            }
                            else if (tAtoms[aNbIdx].chemType.compare("C")==0
                                     && tAtoms[aNbIdx].connAtoms.size()==3)
                            {
                                int nH=0;
                                for (unsigned j=0;
                                     j < tAtoms[aNbIdx].connAtoms.size(); j++)
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
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].connAtoms.size()==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (tAtoms[tIdx].connAtoms.size()==3)
                    {
                        // the lone pair contributes two
                        aN=0.0;
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
                }
            }
            else
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {
                            aN =2.0;
                        }
                        else if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                            // cancel formal charge effect
                            //
                            // aN =2.0;
                        }
                    }
                    else if (tAtoms[tIdx].charge==1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                            // cancel formal charge effect
                            //
                            //aN = 2.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (tAtoms[tIdx].charge==1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            aN =1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {
                            aN =1.0;
                        }
                    }
                }
            }
        }
        else if (tAtoms[tIdx].bondingIdx ==3
                 && (tAtoms[tIdx].chemType.compare("N")==0
                     || tAtoms[tIdx].chemType.compare("B")==0))
                     // || tAtoms[tIdx].chemType.compare("P")==0
                     // || tAtoms[tIdx].chemType.compare("S")==0))
        {
            if (tAtoms[tIdx].chemType.compare("N")==0)
            {
                if (tAtoms[tIdx].charge==-1.0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==2)
                    {
                            aN =2.0;
                        }
                    }
                else if (tAtoms[tIdx].formalCharge==1.0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==3)
                    {                        // Place holder in case for future.
                        aN=1.0;
                    }
                }
                else
                {
                    aN=2.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("B")==0)
            {
                aN=0.0;
            }
            else if (tAtoms[tIdx].chemType.compare("P")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {
                    aN=1.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("S")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {

                    aN=1.0;
                }
            }
        }
        return aN;
    }

    extern REAL setPiForOneAtom(int tIdx, std::vector<AtomDict> & tAtoms,
                                int tMode)
    {
        REAL aN=0.0;
        if (tAtoms[tIdx].bondingIdx ==2)
        {

            if (tAtoms[tIdx].charge==0.0)
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==3)
                    {
                        bool aD=false;
                        for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size();
                           i++)
                        {
                            // check if there is a double bonded atom outside
                            // the ring
                            int aNbIdx = tAtoms[tIdx].connAtoms[i];
                            if (tAtoms[aNbIdx].chemType.compare("O")==0
                                && tAtoms[aNbIdx].connAtoms.size()==1
                                && tAtoms[aNbIdx].charge ==0)
                            {
                                // double bonded Oxy
                                aD = true;
                            }
                            else if (tAtoms[aNbIdx].chemType.compare("C")==0
                                     && tAtoms[aNbIdx].connAtoms.size()==3)
                            {
                                int nH=0;
                                for (unsigned j=0;
                                     j < tAtoms[aNbIdx].connAtoms.size(); j++)
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
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].connAtoms.size()==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (tAtoms[tIdx].connAtoms.size()==3)
                    {
                        // the lone pair contributes two
                        aN=0.0;
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
                }
            }
            else
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {
                            aN =2.0;
                        }
                        else if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            aN = 1.0;
                            //Cancel the charge effect
                            // aN =2.0;
                        }
                    }
                    else if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {

                            //if (tMode ==1)
                            //{
                               // For aromatic ring plane.
                               //  aN=1.0;
                            //}
                            //else
                            //{
                                // For atom classification. It will be
                                // calcelled once the charge problem in COO solved
                                // aN = 2.0;
                            // }
                            aN=1.0;

                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==2)
                        {
                            aN =1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==-1.0)
                    {
                        if (tAtoms[tIdx].connAtoms.size() ==3)
                        {
                            aN =1.0;
                        }
                    }
                }
            }
        }
        else if (tAtoms[tIdx].bondingIdx ==3
                 && (tAtoms[tIdx].chemType.compare("N")==0
                     || tAtoms[tIdx].chemType.compare("B")==0))
                     // || tAtoms[tIdx].chemType.compare("P")==0
                     // || tAtoms[tIdx].chemType.compare("S")==0))
        {
            if (tAtoms[tIdx].chemType.compare("N")==0)
            {
                if (tAtoms[tIdx].formalCharge==-1.0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==2)
                    {
                            aN =2.0;
                    }
                }
                else if (tAtoms[tIdx].formalCharge==1.0)
                {
                    if (tAtoms[tIdx].connAtoms.size() ==3)
                    {                        // Place holder in case for future.
                        aN=1.0;
                    }
                }
                else
                {
                    aN=2.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("B")==0)
            {
                aN=0.0;
            }
            else if (tAtoms[tIdx].chemType.compare("P")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {
                    aN=1.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("S")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {

                    aN=1.0;
                }
            }
        }

        /*
        std::cout << "atom " << tAtoms[tIdx].id << std::endl;
        std::cout << "bond idx " << tAtoms[tIdx].bondingIdx <<std::endl;
        std::cout << "charge is " << tAtoms[tIdx].charge << std::endl;
        std::cout << "formalCharge is " << tAtoms[tIdx].formalCharge
                  << std::endl;
        std::cout << "HERE Pi atom added " << aN  << std::endl;
        */

        return aN;
    }


    extern REAL setPiForOneAtomNoMetal(int tIdx, std::vector<AtomDict> & tAtoms,
                                       int tMode)
    {
        REAL aN=0.0;
        int nonMC=0;
        std::cout << "For atom " << tAtoms[tIdx].id
                  << std::endl;
        for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size(); i++)
        {
            int aIdx = tAtoms[tIdx].connAtoms[i];
            std::cout << "is bonding atom "
                      << tAtoms[aIdx].id << "  metal ?"
                      << tAtoms[aIdx].isMetal << std::endl;
            if (!tAtoms[aIdx].isMetal)
            {
                nonMC++;
            }
        }
        std::cout << "Here atom " << tAtoms[tIdx].id << " has "
                  << nonMC << " nonMetal bonding " << std::endl;
        if (tAtoms[tIdx].bondingIdx ==2)
        {

            if (tAtoms[tIdx].charge==0.0)
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (nonMC ==3)
                    {
                        bool aD=false;
                        for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size();
                           i++)
                        {
                            // check if there is a double bonded atom outside
                            // the ring
                            int aNbIdx = tAtoms[tIdx].connAtoms[i];
                            if (tAtoms[aNbIdx].chemType.compare("O")==0
                                && tAtoms[aNbIdx].connAtoms.size()==1
                                && tAtoms[aNbIdx].charge ==0)
                            {
                                // double bonded Oxy
                                aD = true;
                            }
                            else if (tAtoms[aNbIdx].chemType.compare("C")==0
                                     && tAtoms[aNbIdx].connAtoms.size()==3)
                            {
                                int nH=0;
                                for (unsigned j=0;
                                     j < tAtoms[aNbIdx].connAtoms.size(); j++)
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
                    else if (nonMC ==2)
                    {
                        // Place holder in case for future.
                        aN=0.0;
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (nonMC==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (nonMC==3)
                    {
                        // the lone pair contributes two
                        aN=2.0;
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (nonMC==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (nonMC==3)
                    {
                        // the lone pair contributes two
                        aN=0.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (nonMC==2)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("S") ==0)
                {
                    if (nonMC==2)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("P") ==0)
                {
                    if (nonMC==3)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
            }
            else
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (nonMC ==3)
                        {
                            aN =2.0;
                        }
                        else if (nonMC ==2)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (nonMC ==2)
                        {
                            // aN = 1.0;
                            //Cancel the charge effect
                            aN =2.0;
                        }
                    }
                    else if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (nonMC ==3)
                        {

                            //if (tMode ==1)
                            //{
                               // For aromatic ring plane.
                               //  aN=1.0;
                            //}
                            //else
                            //{
                                // For atom classification. It will be
                                // calcelled once the charge problem in COO solved
                                // aN = 2.0;
                            // }
                            aN=1.0;

                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (nonMC ==2)
                        {
                            aN =1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==-1.0)
                    {
                        if (nonMC ==3)
                        {
                            aN =1.0;
                        }
                    }
                }
            }
        }
        else if (tAtoms[tIdx].bondingIdx ==3
                 && (tAtoms[tIdx].chemType.compare("N")==0
                     || tAtoms[tIdx].chemType.compare("B")==0))
                     // || tAtoms[tIdx].chemType.compare("P")==0
                     // || tAtoms[tIdx].chemType.compare("S")==0))
        {
            if (tAtoms[tIdx].chemType.compare("N")==0)
            {
                if (tAtoms[tIdx].formalCharge==-1.0)
                {
                    if (nonMC ==2)
                    {
                            aN =2.0;
                    }
                }
                else if (tAtoms[tIdx].formalCharge==1.0)
                {
                    if (nonMC ==3)
                    {                        // Place holder in case for future.
                        aN=1.0;
                    }
                }
                else
                {
                    aN=2.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("B")==0)
            {
                aN=0.0;
            }
            else if (tAtoms[tIdx].chemType.compare("P")==0)
            {
                if (nonMC==4)
                {
                    aN=1.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("S")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {

                    aN=1.0;
                }
            }
        }

        /*
        std::cout << "atom " << tAtoms[tIdx].id << std::endl;
        std::cout << "bond idx " << tAtoms[tIdx].bondingIdx <<std::endl;
        std::cout << "charge is " << tAtoms[tIdx].charge << std::endl;
        std::cout << "formalCharge is " << tAtoms[tIdx].formalCharge
                  << std::endl;
        std::cout << "HERE Pi atom added " << aN  << std::endl;
        */

        return aN;
    }

extern REAL setPiForOneAtomAll(int tIdx, std::vector<AtomDict> & tAtoms,
                               int tMode)
    {
        REAL aN=0.0;
        int nonMC=tAtoms[tIdx].connAtoms.size();
        std::cout << "For atom " << tAtoms[tIdx].id
                  << std::endl;

        std::cout << "Here atom " << tAtoms[tIdx].id << " has "
                  << nonMC << "  bonding " << std::endl;
        if (tAtoms[tIdx].bondingIdx ==2)
        {

            if (tAtoms[tIdx].charge==0.0)
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (nonMC ==3)
                    {
                        bool aD=false;
                        for (unsigned i=0; i < tAtoms[tIdx].connAtoms.size();
                           i++)
                        {
                            // check if there is a double bonded atom outside
                            // the ring
                            int aNbIdx = tAtoms[tIdx].connAtoms[i];
                            if (tAtoms[aNbIdx].chemType.compare("O")==0
                                && tAtoms[aNbIdx].connAtoms.size()==1
                                && tAtoms[aNbIdx].charge ==0)
                            {
                                // double bonded Oxy
                                aD = true;
                            }
                            else if (tAtoms[aNbIdx].chemType.compare("C")==0
                                     && tAtoms[aNbIdx].connAtoms.size()==3)
                            {
                                int nH=0;
                                for (unsigned j=0;
                                     j < tAtoms[aNbIdx].connAtoms.size(); j++)
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
                    else if (nonMC ==2)
                    {
                        // Place holder in case for future.
                        aN=0.0;
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (nonMC==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (nonMC==3)
                    {
                        // the lone pair contributes two
                        aN=2.0;
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (nonMC==2)
                    {
                        // one electron becomes a pi electron.
                        aN=1.0;
                    }
                    else if (nonMC==3)
                    {
                        // the lone pair contributes two
                        aN=0.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (nonMC==2)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("S") ==0)
                {
                    if (nonMC==2)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
                else if(tAtoms[tIdx].chemType.compare("P") ==0)
                {
                    if (nonMC==3)
                    {
                        // one lone pair contributes two
                        aN=2.0;
                    }
                }
            }
            else
            {
                if (tAtoms[tIdx].chemType.compare("C") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (nonMC ==3)
                        {
                            aN =2.0;
                        }
                        else if (nonMC ==2)
                        {
                            // Place holder in case for future.
                            aN=1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("N") ==0)
                {
                    if (tAtoms[tIdx].charge==-1.0)
                    {
                        if (nonMC ==2)
                        {
                            // aN = 1.0;
                            //Cancel the charge effect
                            aN =2.0;
                        }
                    }
                    else if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (nonMC ==3)
                        {

                            //if (tMode ==1)
                            //{
                               // For aromatic ring plane.
                               //  aN=1.0;
                            //}
                            //else
                            //{
                                // For atom classification. It will be
                                // calcelled once the charge problem in COO solved
                                // aN = 2.0;
                            // }
                            aN=1.0;

                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("O") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==1.0)
                    {
                        if (nonMC ==2)
                        {
                            aN =1.0;
                        }
                    }
                }
                else if (tAtoms[tIdx].chemType.compare("B") ==0)
                {
                    if (tAtoms[tIdx].formalCharge==-1.0)
                    {
                        if (nonMC ==3)
                        {
                            aN =1.0;
                        }
                    }
                }
            }
        }
        else if (tAtoms[tIdx].bondingIdx ==3
                 && (tAtoms[tIdx].chemType.compare("N")==0
                     || tAtoms[tIdx].chemType.compare("B")==0))
                     // || tAtoms[tIdx].chemType.compare("P")==0
                     // || tAtoms[tIdx].chemType.compare("S")==0))
        {
            if (tAtoms[tIdx].chemType.compare("N")==0)
            {
                if (tAtoms[tIdx].formalCharge==-1.0)
                {
                    if (nonMC ==2)
                    {
                            aN =2.0;
                    }
                }
                else if (tAtoms[tIdx].formalCharge==1.0)
                {
                    if (nonMC ==3)
                    {                        // Place holder in case for future.
                        aN=1.0;
                    }
                }
                else
                {
                    aN=2.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("B")==0)
            {
                aN=0.0;
            }
            else if (tAtoms[tIdx].chemType.compare("P")==0)
            {
                if (nonMC==4)
                {
                    aN=1.0;
                }
            }
            else if (tAtoms[tIdx].chemType.compare("S")==0)
            {
                if (tAtoms[tIdx].connAtoms.size()==4)
                {

                    aN=1.0;
                }
            }
        }

        /*
        std::cout << "atom " << tAtoms[tIdx].id << std::endl;
        std::cout << "bond idx " << tAtoms[tIdx].bondingIdx <<std::endl;
        std::cout << "charge is " << tAtoms[tIdx].charge << std::endl;
        std::cout << "formalCharge is " << tAtoms[tIdx].formalCharge
                  << std::endl;
        std::cout << "HERE Pi atom added " << aN  << std::endl;
        */

        return aN;
    }

    extern REAL setPiForOneAtom(int tIdx, std::vector<AtomDict> & tAtoms,
                                std::vector<BondDict> & tBonds)
    {
        REAL aN=0.0;



        return aN;
    }


    extern REAL setPiForOne_S_Sp3_Atom(int tIdx, std::vector<int>  & tAtmIdx,
                                      std::vector<AtomDict> & tAtoms)
    {
        std::cout << "find num of pi electrons for sp3 S " << std::endl;

        REAL aN=1.0;

        std::vector<int>   inAtmIdx;
        for (std::vector<int>::iterator iN=tAtoms[tIdx].connAtoms.begin();
                iN != tAtoms[tIdx].connAtoms.end(); iN++)
        {
            if (std::find(tAtmIdx.begin(), tAtmIdx.end(), *iN)!=tAtmIdx.end())
            {
                inAtmIdx.push_back(*iN);
            }
        }

        if (inAtmIdx.size() !=2)
        {
            std::cout << "Bug: index error in atoms in a ring " << std::endl;
            exit(1);
        }
        else
        {
            int N3=0;
            for (unsigned i=0; i < inAtmIdx.size(); i++)
            {
                if (tAtoms[inAtmIdx[i]].chemType.compare("N")==0)
                {
                    if (tAtoms[inAtmIdx[i]].connAtoms.size()==3)
                    {
                        N3++;

                    }
                }
            }

            if (N3 > 0)
            {
                aN=0.0;
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
            // bool lS_sp3 = false;
            std::cout << "Ring " << iR->rep << std::endl;
            std::vector<int> atmIdx;
            for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
            {
                atmIdx.push_back(iAt->seriNum);
                std::cout<< "Atom " << iAt->id << std::endl;
                std::cout << "sp " << iAt->bondingIdx << std::endl;
            }

            if (iR->isPlanar)
            {
                std::cout << "One Planar Ring! It contains " << iR->atoms.size()
                          << " atoms. They are : " << std::endl;


                bool lAro = checkAromaSys(atmIdx, tAtoms, 0);
                if (lAro)
                {
                    iR->isAromatic = true;
                    tAroRings.push_back(*iR);
                    // std::cout << "It is an aromatic ring " << std::endl;
                }
                else
                {
                    iR->isAromatic = false;
                    REAL numPi = getTotalPiElec(atmIdx, tAtoms, 0);
                    if (numPi > 0.0 && fabs(fmod(numPi, 4.0)) < 0.001)
                    {

                        iR->isAntiAroma = true;

                    }
                    std::cout << "2 number of pi elecs is "
                                  <<numPi << std::endl;
                    std::cout << fabs(fmod(numPi, 4.0)) << std::endl;
                    std::cout << iR->isAntiAroma << std::endl;
                    // std::cout << "It is not an aromatic ring " << std::endl;
                }
                bool lAroP = checkAromaSys(atmIdx, tAtoms, 1);
                if (lAroP)
                {
                    iR->isAromaticP = true;
                    tAroRings.push_back(*iR);
                    std::cout << "3 It is an aromatic ring " << std::endl;
                }
                //std::cout << std::endl;
            }
            else
            {

                bool lNBUnR=checkUndRing(atmIdx, tAtoms);
                std::cout << "Ring with sp3 atom " << lNBUnR << std::endl;

                if (lNBUnR)
                {
                    bool lAroP = checkAromaSys(atmIdx, tAtoms,1);

                    if (lAroP)
                    {
                        iR->isAromaticP = true;
                        iR->isPlanar   = true;
                        tAroRings.push_back(*iR);
                    }
                    else
                    {
                        iR->isAromaticP = false;
                    }

                    bool lAro = checkAromaSys(atmIdx, tAtoms,0);

                    if (lAro)
                    {
                        iR->isAromatic = true;
                        tAroRings.push_back(*iR);
                    }
                    else
                    {
                        iR->isAromatic = false;
                    }

                }
                else
                {
                    iR->isAromatic =false;
                }
            }

            // std::cout << "Is the ring aromatic ? " << iR->isAromatic << std::endl;
        }

        //std::cout << "A: Number of rings "
        //          << tAllRings.size() << std::endl;
        std::vector<std::vector<int> >  mergedRingSets;
        mergePlaneRings(tAllRings, mergedRingSets, tAtoms);

        //std::cout << "B: Number of rings "
        //          << tAllRings.size() << std::endl;

        tPlanes.clear();

        // setAllRingPlanes(tAllRings, tAtoms, tPlanes);
        //setAllRingPlanes2(tAllRings, mergedRingSets, tAtoms, tPlanes);
        setAllRingPlanes4(tAllRings, tAtoms, tPlanes);


        setAllOtherPlanes(tAllRings, tAtoms, tPlanes);

        std::cout << "There are "         << tPlanes.size()
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


   extern void checkAndSetupPlanes(std::vector<RingDict>  & tAllRings,
                                    std::vector<PlaneDict> & tPlanes,
                                    std::vector<AtomDict>  & tAtoms,
                                    bool                     tMdPls)
    {

        // Check aromaticity for individual rings

        std::vector<RingDict> tAroRings;
        // std::cout << "\nCheck aromaticity for individual ring " << std::endl;
        for (std::vector<RingDict>::iterator iR=tAllRings.begin();
                iR !=tAllRings.end(); iR++)
        {
            // bool lS_sp3 = false;

            std::vector<int> atmIdx;
            for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
            {
                atmIdx.push_back(iAt->seriNum);
                // std::cout << iAt->id << std::endl;

            }



            if (iR->isPlanar)
            {
                std::cout << "One Planar Ring! It contains " << iR->atoms.size()
                          << " atoms. They are : " << std::endl;


                bool lAro = checkAromaSys(atmIdx, tAtoms, 0);
                if (lAro)
                {
                    iR->isAromatic = true;
                    tAroRings.push_back(*iR);
                    // std::cout << "It is an aromatic ring " << std::endl;
                }
                else
                {
                    iR->isAromatic = false;
                    REAL numPi = getTotalPiElec(atmIdx, tAtoms, 0);
                    if (numPi > 0.0 && fabs(fmod(numPi, 4.0)) < 0.001)
                    {

                        iR->isAntiAroma = true;

                    }
                    std::cout << "1 number of pi elecs is "
                                  <<numPi << std::endl;
                    std::cout << fabs(fmod(numPi, 4.0)) << std::endl;
                    std::cout << iR->isAntiAroma << std::endl;
                    // std::cout << "It is not an aromatic ring " << std::endl;
                }
                bool lAroP = checkAromaSys(atmIdx, tAtoms, 1);
                if (lAroP)
                {
                    iR->isAromaticP = true;
                    tAroRings.push_back(*iR);
                    // std::cout << "It is an aromatic ring " << std::endl;
                }
                //std::cout << std::endl;
            }
            else
            {

                bool lNBUnR=checkUndRing(atmIdx, tAtoms);
                // std::cout << "Ring with sp3 atom " << lNBUnR << std::endl;

                if (lNBUnR)
                {
                    bool lAroP = checkAromaSys(atmIdx, tAtoms,1);

                    if (lAroP)
                    {
                        iR->isAromaticP = true;
                        iR->isPlanar   = true;
                        tAroRings.push_back(*iR);
                    }
                    else
                    {
                        iR->isAromaticP = false;
                    }

                    bool lAro = checkAromaSys(atmIdx, tAtoms,0);

                    if (lAro)
                    {
                        iR->isAromatic = true;
                        tAroRings.push_back(*iR);
                    }
                    else
                    {
                        iR->isAromatic = false;
                    }

                }
                else
                {
                    iR->isAromatic =false;
                }
            }

            // std::cout << "Is the ring aromatic ? " << iR->isAromatic << std::endl;
        }

        //std::cout << "A: Number of rings "
        //          << tAllRings.size() << std::endl;
        std::vector<std::vector<int> >  mergedRingSets;
        mergePlaneRings(tAllRings, mergedRingSets, tAtoms);

        //std::cout << "B: Number of rings "
        //          << tAllRings.size() << std::endl;

        tPlanes.clear();

        // setAllRingPlanes(tAllRings, tAtoms, tPlanes);
        if (!tMdPls)
        {
            setAllRingPlanes3(tAllRings, tAtoms, tPlanes);
        }
        else
        {
            setAllRingPlanes4(tAllRings, tAtoms, tPlanes);
        }

        setAllOtherPlanes(tAllRings, tAtoms, tPlanes);

        std::cout << "There are "         << tPlanes.size()
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



    extern void checkAndSetupPlanes2(std::vector<RingDict>  & tAllRings,
                                     std::vector<PlaneDict> & tPlanes,
                                     std::vector<AtomDict>  & tAtoms,
                                     std::vector<BondDict>  & tBonds)
    {

        // Check aromaticity for individual rings
        std::vector<RingDict> tAroRings;
        // std::cout << "\nCheck aromaticity for individual ring " << std::endl;
        for (std::vector<RingDict>::iterator iR=tAllRings.begin();
                iR !=tAllRings.end(); iR++)
        {
            // bool lS_sp3 = false;

            std::vector<int> atmIdx;
            for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                    iAt !=iR->atoms.end(); iAt++)
            {
                atmIdx.push_back(iAt->seriNum);
                // std::cout << iAt->id << std::endl;

            }



            if (iR->isPlanar)
            {
                //std::cout << "One Planar Ring! It contains " << iR->atoms.size()
                //          << " atoms. They are : " << std::endl;


                bool lAro = checkAromaSys(atmIdx, tAtoms, tBonds);
                if (lAro)
                {
                    iR->isAromatic = true;
                    tAroRings.push_back(*iR);
                    // std::cout << "It is an aromatic ring " << std::endl;
                }
                else
                {
                    iR->isAromatic = false;
                    // std::cout << "It is not an aromatic ring " << std::endl;
                }

                //std::cout << std::endl;
            }
            else
            {

                bool lNBUnR=checkUndRing(atmIdx, tAtoms);
                // std::cout << "Ring with sp3 atom " << lNBUnR << std::endl;

                if (lNBUnR)
                {
                    bool lAro = checkAromaSys(atmIdx, tAtoms, tBonds);

                    if (lAro)
                    {
                        iR->isAromatic = true;
                        iR->isPlanar   = true;
                        tAroRings.push_back(*iR);
                    }
                    else
                    {
                        iR->isAromatic = false;
                    }
                }
                else
                {
                    iR->isAromatic =false;
                }
            }

            // std::cout << "Is the ring aromatic ? " << iR->isAromatic << std::endl;


        }

        //std::cout << "A: Number of rings "
        //          << tAllRings.size() << std::endl;
        std::vector<std::vector<int> >  mergedRingSets;
        mergePlaneRings(tAllRings, mergedRingSets, tAtoms);

        //std::cout << "B: Number of rings "
        //          << tAllRings.size() << std::endl;

        tPlanes.clear();

        // setAllRingPlanes(tAllRings, tAtoms, tPlanes);
        // setAllRingPlanes2(tAllRings, mergedRingSets, tAtoms, tPlanes);
        setAllRingPlanes3(tAllRings, tAtoms, tPlanes);
        //std::cout << "Here 4 for planes " << std::endl;

        setAllOtherPlanes(tAllRings, tAtoms, tPlanes);
        /*
        std::cout << "There are "         << tPlanes.size()
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
        */
    }


    extern bool checkUndRing(std::vector<int>     &  tAtmIdxs,
                            std::vector<AtomDict> & tAtoms)
    {
        // allow aromaticity to decide N, B hybridization and ring planes

        bool lNB = true;

        for (std::vector<int>::iterator iAt=tAtmIdxs.begin();
                iAt != tAtmIdxs.end(); iAt++)
        {
            if (tAtoms[*iAt].bondingIdx !=2
                //&& tAtoms[*iAt].chemType.compare("N") !=0
                && tAtoms[*iAt].chemType.compare("B") !=0)
                //&& tAtoms[*iAt].chemType.compare("P") !=0
                //&& tAtoms[*iAt].chemType.compare("S") !=0)
            {
                lNB = false;
                break;
            }
        }

        return lNB;

    }
    extern void TestAromaticity(std::vector<RingDict>  & tAllRings,
                                std::vector<AtomDict>  & tAtoms)
    {

        //std::cout << "\nCheck aromaticity for individual ring " << std::endl;

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

    extern bool checkAllARBondsInOneRing( std::vector<BondDict> & tBonds,
                                          RingDict & tRing)
    {
        bool lAR = true;

        for (unsigned i=0; i < tRing.bondIdxs.size(); i++)
        {
            if (tBonds[tRing.bondIdxs[i]].order.find("AR")
                 ==tBonds[tRing.bondIdxs[i]].order.npos)
            {
                lAR=false;
                break;
            }
        }
        return lAR;
    }

    extern void setAromaticBonds(std::vector<RingDict>  & tRings,
                                 std::vector<BondDict>  & tBonds)
    {
        /*
        for (std::vector<BondDict>::iterator iB=tBonds.begin();
                iB !=tBonds.end(); iB++)
        {
            std::cout << "Bond between atom " << iB->atoms[0]
                      << " and " << iB->atoms[1] << std::endl;
        }
        */

        for (std::vector<RingDict>::iterator iR=tRings.begin();
                iR !=tRings.end(); iR++)
        {
            if (iR->isAromaticP || iR->isAromatic)
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
                                //std::cout << "Bond between atom "
                                ///          << iR->atoms[i].id << " and "
                                //          << iR->atoms[j].id << " is aromatic "
                                //          << std::endl;
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

    extern bool detectAllSp2AtomRing(RingDict               & tRing)
    {
        bool aReturn = true;

        for (std::vector<AtomDict>::iterator iAt=tRing.atoms.begin();
                iAt != tRing.atoms.end(); iAt++)
        {
            if (iAt->bondingIdx !=2)
            {
                aReturn = false;
                break;
            }
        }

        return aReturn;
    }

    extern void setSugarRingInitComf(std::vector<AtomDict>     & tAtoms,
                                     std::vector<TorsionDict>  & tTors,
                                     std::vector<RingDict>::iterator tRing)
    {
        checkOneSugarRing(tAtoms, tRing);
        if (tRing->sugarType.compare("pyranose")==0)
        {
            std::cout << "Find one pyranose ring " << std::endl;
            // A pyranose ring, set torsions within the ring
            setPyranoseChairComf(tAtoms, tRing, tTors);
        }
    }

    extern void checkOneSugarRing(std::vector<AtomDict> & tAtoms,
                                  std::vector<RingDict>::iterator tRing)
    {
        // Six-member sugar ring (Pyranose)

        // std::cout << "Check pyranose " << std::endl;
        std::vector<ID> rAtmIds;

        for (std::vector<AtomDict>::iterator iAt = tRing->atoms.begin();
                iAt !=tRing->atoms.end(); iAt++)
        {
            rAtmIds.push_back(iAt->id);
        }

        std::map<ID, int>  rAtmFormat;
        rAtmFormat["O"]        = 0;
        rAtmFormat["C"]        = 0;
        rAtmFormat["connectC"] = 0;
        rAtmFormat["connectO"] = 0;

        for (std::vector<AtomDict>::iterator iAt = tRing->atoms.begin();
                iAt !=tRing->atoms.end(); iAt++)
        {
            if (iAt->chemType.compare("C")==0 && iAt->bondingIdx==3)
            {
                 rAtmFormat["C"]++;

                 int iCO=0;
                 for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                         iNB !=iAt->connAtoms.end(); iNB++)
                 {
                     if (tAtoms[*iNB].chemType.compare("O")==0
                         && std::find(rAtmIds.begin(), rAtmIds.end(), tAtoms[*iNB].id) == rAtmIds.end())
                     {
                         iCO++;
                     }
                     else if (tAtoms[*iNB].chemType.compare("C")==0
                         && std::find(rAtmIds.begin(), rAtmIds.end(), tAtoms[*iNB].id) == rAtmIds.end())
                     {
                         rAtmFormat["connectC"]++;
                     }
                 }

                 if (iCO > 1)
                 {
                     break;
                 }
                 else if (iCO==1)
                 {
                     rAtmFormat["connectO"]++;
                 }
            }
            else if (iAt->chemType.compare("O")==0)
            {
                if (iAt->connAtoms.size()==2)
                {
                    int iOC=0;
                    for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                            iNB !=iAt->connAtoms.end(); iNB++)
                    {
                        if (tAtoms[*iNB].chemType.compare("C")==0
                            &&  std::find(rAtmIds.begin(), rAtmIds.end(), tAtoms[*iNB].id) != rAtmIds.end())
                        {
                            iOC++;
                        }
                    }
                    if (iOC==2)
                    {
                        rAtmFormat["O"]++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            else
            {
                break;
            }

            if (rAtmFormat["O"]==1 && rAtmFormat["C"]==5
                && rAtmFormat["connectO"] >1)
            {
                tRing->sugarType = "pyranose";
            }
            else if (rAtmFormat["O"]==1 && rAtmFormat["C"]==4
                && rAtmFormat["connectO"]==3 && rAtmFormat["connectC"]==2 )
            {
                tRing->sugarType = "furanose";
            }
        }
        /*
        if (tRing->sugarType.size() !=0)
        {
            std::cout << "The ring is " << tRing->sugarType << std::endl;
        }
        else
        {
            std::cout << "The ring is not sugar ring " << std::endl;
        }
        */
        //


    }


    extern void setPyranoseChairComf(std::vector<AtomDict> & tAtoms,
                                      std::vector<RingDict>::iterator tRing)
    {
        // A template function to set all six-member rings of pyranose to the "chair" conformation.
        // Using a proper class to generate Carbohydrate Conformations later on.

        tRing->setRingAtmsLinks();

        std::vector<int> branch1, branch2;

        for (unsigned i0=0; i0 < tRing->atoms.size(); i0++)
        {
            if (tRing->atoms[i0].chemType.compare("O")==0)
            {
                // O5 atom
                branch1.push_back(tRing->atoms[i0].seriNum);
                branch2.push_back(tRing->atoms[i0].seriNum);



                if (tRing->atoms[i0].connAtoms.size() ==2)
                {
                    int i1=-1, i2=-1, i3=-1, i4=1, i5=-1;

                    int it1=tRing->atoms[i0].connAtoms[0];
                    int it2=tRing->atoms[i0].connAtoms[1];



                    bool lConnO = false;
                    for (unsigned inb1=0; inb1 !=tAtoms[it1].connAtoms.size();
                            inb1++)
                    {
                         if (tAtoms[inb1].chemType.compare("O")==0
                             && tAtoms[inb1].seriNum != tRing->atoms[i0].seriNum)
                         {
                             lConnO = true;
                             break;
                         }
                    }

                    if (lConnO)
                    {
                        i1 = it1;
                        i2 = it2;
                    }
                    else
                    {
                        i1 = it2;
                        i2 = it1;
                    }

                    branch1.push_back(i1);
                    branch2.push_back(i2);



                    for (std::vector<int>::iterator iNext=tRing->ringAtomLink[i1].begin();
                            iNext != tRing->ringAtomLink[i1].end(); iNext++)
                    {
                        if (*iNext !=tRing->atoms[i0].seriNum)
                        {
                            i3 = *iNext;
                            break;
                        }
                    }

                    if (i3 !=-1)
                    {
                        branch1.push_back(i3);

                        for (std::vector<int>::iterator iNext=tRing->ringAtomLink[i3].begin();
                            iNext != tRing->ringAtomLink[i3].end(); iNext++)
                        {
                            if (*iNext !=tRing->atoms[i1].seriNum)
                            {
                                i5 = *iNext;

                                break;
                            }
                        }
                    }
                    else
                    {
                        std::cout << "Can not find one of ring atoms connected to atom "
                                  << tAtoms[i1].id << std::endl;
                        exit(1);
                    }

                    for (std::vector<int>::iterator iNext = tRing->ringAtomLink[i2].begin();
                            iNext != tRing->ringAtomLink[i2].end(); iNext++)
                    {
                        if (*iNext !=tRing->atoms[i0].seriNum)
                        {
                            i4 = *iNext;
                            break;
                        }
                    }

                    if (i4 !=-1)
                    {
                        branch2.push_back(i4);

                    }
                    else
                    {
                        std::cout << "Can not find one of ring atoms connected to atom "
                                  << tAtoms[i2].id << std::endl;
                        exit(1);
                    }

                    if (i5 !=-1)
                    {
                        branch1.push_back(i5);

                        branch2.push_back(i5);
                    }
                }
                else
                {
                    std::cout << "Bug: O atom in the sugar ring connects "
                              << tRing->atoms[i0].connAtoms.size() << " atoms"
                              << std::endl;
                    exit(1);
                }

                break;

            }
        }

        if (branch1.size() ==4 && branch2.size() ==4)
        {
            ID iD0  = tAtoms[branch1[0]].id;
            ID iD11 = tAtoms[branch1[1]].id;
            ID iD12 = tAtoms[branch1[2]].id;
            ID iD13 = tAtoms[branch1[3]].id;

            ID iD21 = tAtoms[branch2[1]].id;
            ID iD22 = tAtoms[branch2[2]].id;

            // Temporal values.
            // The typical id combination C5-O5-C1-C2
            tRing->sugarTors[iD21 + "_" + iD0  + "_" + iD11 + "_" + iD12 ] = -60.0;
            tRing->sugarTors[iD0  + "_" + iD11 + "_" + iD12 + "_" + iD13 ] =  60.0;
            tRing->sugarTors[iD11 + "_" + iD12 + "_" + iD13 + "_" + iD22 ] = -60.0;
            tRing->sugarTors[iD12 + "_" + iD13 + "_" + iD22 + "_" + iD21 ] =  60.0;
            tRing->sugarTors[iD13 + "_" + iD22 + "_" + iD21 + "_" + iD0  ] = -60.0;
            tRing->sugarTors[iD22 + "_" + iD21 + "_" + iD0  + "_" + iD11 ] =  60.0;

        }

    }

    extern void setPyranoseChairComf(std::vector<AtomDict> & tAtoms,
                                     std::vector<RingDict>::iterator tRing,
                                     std::vector<TorsionDict> & tTors)
    {
        // A template function to set all six-member rings of pyranose to the "chair" conformation.
        // Using a proper class to generate Carbohydrate Conformations later on.

        tRing->setRingAtmsLinks();

        std::vector<int> branch1, branch2;

        for (unsigned i0=0; i0 < tRing->atoms.size(); i0++)
        {
            if (tRing->atoms[i0].chemType.compare("O")==0)
            {
                // O5 atom
                branch1.push_back(tRing->atoms[i0].seriNum);
                branch2.push_back(tRing->atoms[i0].seriNum);

                // std::cout << "Starting atom " << tRing->atoms[i0].id << std::endl;

                if (tRing->atoms[i0].connAtoms.size() ==2)
                {
                    int i1=-1, i2=-1, i3=-1, i4=1, i5=-1;

                    int it1=tRing->atoms[i0].connAtoms[0];
                    int it2=tRing->atoms[i0].connAtoms[1];


                    bool lConnO = false;
                    for (unsigned inb1=0; inb1 !=tAtoms[it1].connAtoms.size();
                            inb1++)
                    {
                         if (tAtoms[inb1].chemType.compare("O")==0
                             && tAtoms[inb1].seriNum != tRing->atoms[i0].seriNum)
                         {
                             lConnO = true;
                             break;
                         }
                    }

                    if (lConnO)
                    {
                        i1 = it1;
                        i2 = it2;
                    }
                    else
                    {
                        i1 = it2;
                        i2 = it1;
                    }

                    branch1.push_back(i1);
                    branch2.push_back(i2);

                    // std::cout << "Branch 1 connected " <<   tAtoms[i1].id << std::endl;
                    // std::cout << "Branch 2 connected " <<   tAtoms[i2].id << std::endl;

                    for (std::vector<int>::iterator iNext=tRing->ringAtomLink[i1].begin();
                            iNext != tRing->ringAtomLink[i1].end(); iNext++)
                    {
                        if (*iNext !=tRing->atoms[i0].seriNum)
                        {
                            i3 = *iNext;
                            break;
                        }
                    }

                    if (i3 !=-1)
                    {
                        branch1.push_back(i3);
                        // std::cout << "Branch 1 connected " <<   tAtoms[i3].id << std::endl;
                        for (std::vector<int>::iterator iNext=tRing->ringAtomLink[i3].begin();
                            iNext != tRing->ringAtomLink[i3].end(); iNext++)
                        {
                            if (*iNext !=tAtoms[i1].seriNum)
                            {
                                i5 = *iNext;
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::cout << "Can not find one of ring atoms connected to atom "
                                  << tAtoms[i1].id << std::endl;
                        exit(1);
                    }

                    for (std::vector<int>::iterator iNext = tRing->ringAtomLink[i2].begin();
                            iNext != tRing->ringAtomLink[i2].end(); iNext++)
                    {
                        //std::cout << "iPre " << tAtoms[i0].id << std::endl;
                        //std::cout << "iCen " << tAtoms[i2].id << std::endl;
                        //std::cout << "iNext " << tAtoms[*iNext].id << std::endl;
                        if (*iNext !=tRing->atoms[i0].seriNum)
                        {
                            i4 = *iNext;

                            break;
                        }
                    }

                    if (i4 !=-1)
                    {
                        branch2.push_back(i4);
                        //std::cout << "Branch 2 connected " <<   tAtoms[i4].id << std::endl;
                    }
                    else
                    {
                        std::cout << "Can not find one of ring atoms connected to atom "
                                  << tAtoms[i2].id << std::endl;
                        exit(1);
                    }

                    if (i5 !=-1)
                    {
                        branch1.push_back(i5);
                        // std::cout << "Branch 1 connected " <<   tAtoms[i5].id << std::endl;
                        branch2.push_back(i5);
                        // std::cout << "Branch 2 connected " <<   tAtoms[i5].id << std::endl;
                    }
                }
                else
                {
                    std::cout << "Bug: O atom in the sugar ring connects "
                              << tRing->atoms[i0].connAtoms.size() << " atoms"
                              << std::endl;
                    exit(1);
                }

                break;

            }
        }


        if (branch1.size() ==4 && branch2.size() ==4)
        {
            ID  iD0  = tAtoms[branch1[0]].id;
            int n0   = tAtoms[branch1[0]].seriNum;
            ID  iD11 = tAtoms[branch1[1]].id;
            int n11  = tAtoms[branch1[1]].seriNum;
            ID  iD12 = tAtoms[branch1[2]].id;
            int n12  = tAtoms[branch1[2]].seriNum;
            ID  iD13 = tAtoms[branch1[3]].id;
            int n13  = tAtoms[branch1[3]].seriNum;

            ID  iD21 = tAtoms[branch2[1]].id;
            int n21  = tAtoms[branch2[1]].seriNum;
            ID  iD22 = tAtoms[branch2[2]].id;
            int n22  = tAtoms[branch2[2]].seriNum;

            // Temporal values.
            int nT = -1, nTors=(int)tTors.size();
            ID  aTorID = "";
            // The typical id combination C5-O5-C1-C2
            aTorID = iD21 + "_" + iD0  + "_" + iD11 + "_" + iD12;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID]  = -60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n21, n0, n11, n12, -60.0);


            aTorID = iD0  + "_" + iD11 + "_" + iD12 + "_" + iD13;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID ] =  60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n0, n11, n12, n13, 60.0);



            aTorID = iD11 + "_" + iD12 + "_" + iD13 + "_" + iD22;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID] = -60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n11, n12, n13, n22, -60.0);


            aTorID = iD12 + "_" + iD13 + "_" + iD22 + "_" + iD21;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID ] =  60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n12, n13, n22, n21, 60.0);

            aTorID = iD13 + "_" + iD22 + "_" + iD21 + "_" + iD0;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID] = -60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n13, n22, n21, n0, -60.0);

            aTorID = iD22 + "_" + iD21 + "_" + iD0  + "_" + iD11;
            std::cout << aTorID << std::endl;
            tRing->sugarTors[aTorID] =  60.0;
            setTorsionAroundOneBondInRing(tTors, tAtoms, n22, n21, n0, n11, 60.0);
            nT = getTorsion(tTors, n22, n21, n0, n11);

        }

    }

    extern void setTorsionAroundOneBondInRing(std::vector<TorsionDict> & tTors,
                                              std::vector<AtomDict>    & tAtoms,
                                              int rAt1, int rAt2, int rAt3,
                                              int rAt4, REAL vInit)
    {
        // set the torsion angle formed by 4 ring atoms

        int nTors = (int)tTors.size();
        int nT = getTorsion(tTors, rAt1, rAt2, rAt3, rAt4);
        std::cout << "Torsion index here is " << nT << std::endl;

        if (nT >=0 && nT < nTors )
        {
            tTors[nT].value = vInit;
            // ID aTorID = tAtoms[rAt1].id + "_" + tAtoms[rAt2].id + "_"
            //            + tAtoms[rAt3].id + "_" + tAtoms[rAt4].id;
            // std::cout << aTorID << "  torsion " << tTors[nT].value << std::endl;
        }
        else
        {
            ID aTorID = tAtoms[rAt1].id + "_" + tAtoms[rAt2].id + "_"
                        + tAtoms[rAt3].id + "_" + tAtoms[rAt4].id;
            std::cout << "Bug: could not find the torsion angle by atoms "
                          << aTorID << std::endl;
            exit(1);
        }


        // Now re-set torsion angles around two sides of the bond formed
        // rAt2, rAt3.
        // 1. the side of rAt2

        if (tAtoms[rAt3].connAtoms.size() >2 )      // O atoms are included
        {
            if (tAtoms[rAt3].bondingIdx==3)         // sp3. used for pyranose only at the moment
            {
                int n =1;
                int n3=1;
                for (std::vector<int>::iterator iNB=tAtoms[rAt3].connAtoms.begin();
                        iNB !=tAtoms[rAt3].connAtoms.end(); iNB++)
                {
                    if (*iNB !=rAt2)
                    {
                        // starting point is rAt1;
                        if (*iNB != rAt4)
                        {
                            int nNT=getTorsion(tTors, rAt1, rAt2, rAt3, *iNB);
                            if (nNT >=0 && nNT < nTors)
                            {
                                tTors[nNT].value = vInit + n*120.0;
                                if (tTors[nNT].value > 180.01)
                                {
                                    tTors[nNT].value = tTors[nNT].value -360.0;
                                }
                                //ID aTorID =   tAtoms[rAt1].id + "_" + tAtoms[rAt2].id + "_"
                                //            + tAtoms[rAt3].id + "_" + tAtoms[*iNB].id;
                                //std::cout << aTorID << "  torsion " << tTors[nNT].value << std::endl;
                            }
                            else
                            {
                                ID aTorID =   tAtoms[rAt1].id + "_" + tAtoms[rAt2].id + "_"
                                            + tAtoms[rAt3].id + "_" + tAtoms[*iNB].id;
                                std::cout << "Bug: could not find the torsion angle by atoms "
                                          << aTorID << std::endl;
                                exit(1);
                            }
                            n++;
                        }

                        // rAt2 side
                        int n2 = 1;
                        for (std::vector<int>::iterator iNB2=tAtoms[rAt2].connAtoms.begin();
                                     iNB2 !=tAtoms[rAt2].connAtoms.end(); iNB2++)
                        {
                            if(*iNB2 !=rAt3 && *iNB2 != rAt1)
                            {
                                // starting point is rAt1;
                                int nNBT=getTorsion(tTors, *iNB2, rAt2, rAt3, *iNB);
                                if (nNBT >=0 && nT < nTors)
                                {
                                    REAL vInit2 = vInit + n2*120.0;
                                    if (vInit2 > 180.0)
                                    {
                                        vInit2 = vInit2 -360.0;
                                    }
                                    if (*iNB==rAt4)
                                    {
                                        tTors[nNBT].value = vInit2;
                                    }
                                    else
                                    {
                                        tTors[nNBT].value = vInit2 + n3*120.0;
                                    }
                                    if (tTors[nNBT].value > 180.0)
                                    {
                                        tTors[nNBT].value = tTors[nNBT].value -360.0;
                                    }
                                    n2++;
                                    //ID aTorID =   tAtoms[*iNB2].id + "_" + tAtoms[rAt2].id + "_"
                                    //            + tAtoms[rAt3].id  + "_" + tAtoms[*iNB].id;
                                    //std::cout << aTorID << "  torsion " << tTors[nNBT].value << std::endl;

                                }
                                else
                                {
                                    ID aTorID =   tAtoms[*iNB2].id + "_" + tAtoms[rAt2].id + "_"
                                                + tAtoms[rAt3].id  + "_" + tAtoms[*iNB].id;
                                    std::cout << "Bug: could not find the torsion angle by atoms "
                                              << aTorID << std::endl;
                                    exit(1);
                                }
                            }
                        }
                        n3=n;
                    }
                }
            }
        }

        // std::cout <<std::endl << std::endl;


    }

    ringTools::ringTools()
    {
    }

    ringTools::~ringTools()
    {
    }

    void ringTools::detectRingFromAtoms(std::vector<AtomDict> & tAtoms,
                                   std::map<ID, std::vector<RingDict> > & tRings,
                                   int                     tDepth,
                                   int                     tMaxRing)
    {

        std::map<int, ID>  atomsInPath;
        std::map<int, ID>  atomsSeen;

         for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
         {
             if (!iA->isMetal && iA->chemType.compare("H") !=0)
             {
                 int preSeriNum = -999;
                 int startLev   = 1;
                 atomsInPath.clear();
                 atomsSeen.clear();

                 checkOnePathSec(*iA, tMaxRing, iA, preSeriNum,  startLev, atomsSeen, atomsInPath,
                                 tAtoms, tRings);

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
            }
         }
    }

    void ringTools::checkOnePathSec(AtomDict             &             curAto,
                                    int                                tMaxRing,
                                    std::vector<AtomDict>::iterator    iOriAto,
                                    int                                SeriNumPreAto,
                                    int                                curLev,
                                    std::map<int,ID>     &             seenAtomIDs,
                                    std::map<int,ID>     &             atomIDsInPath,
                                    std::vector<AtomDict>     &        tAtoms,
                                    std::map<ID, std::vector<RingDict> > &  tRings)
    {
        if ( curLev <tMaxRing)
        {
            int NachbarpunkteDetected = 0;

            // Check Nachbarpunkte
            for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
            {
                int tSeriNum = tAtoms[*tNBA].seriNum;

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
                    int tSeriNum = tAtoms[*tNBA].seriNum;
                    ID  tNbId    = tAtoms[*tNBA].chemType;
                    if (tSeriNum == iOriAto->seriNum && tSeriNum != SeriNumPreAto
                        && curLev > 2 && tNbId.compare("H") !=0)
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
                        std::list<std::string> tAllSeris;
                        std::vector <AtomDict> ttAtoms;
                        for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                        {
                            //if (iSee->first != iOriAto->seriNum)
                            //{
                            tAllSeris.push_back( IntToStr(iSee->first));
                            tAllIds.push_back(iSee->second);

                            //int posIdx = atomPosition(iSee->second);
                            ttAtoms.push_back(tAtoms[iSee->first]);

                            // tRing.atoms.push_back(*tAt);
                            //std::cout << iSee->second << std::endl;
                            //}
                        }

                        RingDict aRingDict(ttAtoms);
                        tAllSeris.sort(compareNoCase);
                        tAllIds.sort(compareNoCase);
                        //std::string tRepStr(iOriAto->id);
                        std::string tRepStr;
                        std::string tRepSeri;
                        for (std::list<std::string>::iterator iAI =tAllIds.begin();
                                    iAI != tAllIds.end(); iAI++)
                        {
                            //std::cout << "ID " << *iAI << std::endl;
                            tRepStr.append(*iAI);
                        }


                        int nRS =0;
                        for (std::list<std::string>::iterator iAS =tAllSeris.begin();
                                    iAS != tAllSeris.end(); iAS++)
                        {
                            //std::cout << "serial number " << *iAS << std::endl;

                            //std::string aS = *iAS;
                            //int aN = StrToInt(aS);
                            //std::cout << "ID " << tAtoms[aN].id << std::endl;
                            if (nRS==0)
                            {
                                tRepSeri.append(*iAS);
                            }
                            else
                            {
                                tRepSeri.append("_" + *iAS);
                            }
                            nRS++;
                        }


                        iOriAto->ringRep[tRepStr] = (int)atomIDsInPath.size();

                        aRingDict.rep  = tRepStr;
                        aRingDict.sRep = tRepSeri;

                        std::map<ID, std::vector<RingDict> >::iterator iFindRing=tRings.find(tRepSeri);
                        if (iFindRing == tRings.end())
                        {
                            for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                            {
                                // int posIdx = atomPosition(iSee->second);
                                tAtoms[iSee->first].inRings.push_back((int)tRings.size());
                            }
                            tRings[tRepSeri].push_back(aRingDict);
                        }

                        // std::cout << "Reps " << tRepSeri << std::endl;
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
                if (tNewLev < tMaxRing)
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

                            seenAtomIDs.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                            atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                        }
                        if (SeriNumPreAto != tAtoms[*tNBA].seriNum && ! tAtoms[*tNBA].isMetal)
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
                            checkOnePathSec(tAtoms[*tNBA], tMaxRing, iOriAto, tPreSeriNum,
                                            tNewLev, seenAtomIDs, atomIDsInPath,
                                            tAtoms, tRings);
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


    void ringTools::setAtomsRingRepreS(std::vector<AtomDict>              & tAtoms,
                                       std::vector<RingDict>              & tRings)
    {
        for (unsigned i=0; i < tRings.size();  i++)
        {
            std::string sSize= IntToStr(tRings[i].atoms.size());
            std::string aRepId;
            std::list<ID> tAllIds;

            for (std::vector<AtomDict>::iterator iAt=tRings[i].atoms.begin();
                    iAt !=tRings[i].atoms.end(); iAt++)
            {
                tAllIds.push_back(IntToStr(iAt->seriNum));
            }
            tAllIds.sort(compareNoCase);

            int nRS =0;
            for (std::list<std::string>::iterator iAS =tAllIds.begin();
                          iAS != tAllIds.end(); iAS++)
            {
                if (nRS==0)
                {
                    aRepId.append(*iAS);
                }
                else
                {
                    aRepId.append("_" + *iAS);
                }
                nRS++;
            }


            // std::cout << "Ring " << aRepId << std::endl;

            for (std::vector<AtomDict>::iterator iAt=tRings[i].atoms.begin();
                    iAt !=tRings[i].atoms.end(); iAt++)
            {
                if (tRings[i].isAromatic)
                {
                    iAt->ringRepS[aRepId] = sSize +"a";
                    if (!iAt->isInAromRing)
                    {
                        iAt->isInAromRing     = true;
                    }
                }
                else
                {
                    iAt->ringRepS[aRepId] = sSize;
                }

                int aSeri=getAtom(iAt->id, iAt->seriNum, tAtoms);
                if (aSeri >=0)
                {
                    tAtoms[aSeri].ringRepS[aRepId]  =  iAt->ringRepS[aRepId];
                    if (!tAtoms[aSeri].isInAromRing
                        && iAt->ringRepS[aRepId].find("a") != std::string::npos)
                    {
                        tAtoms[aSeri].isInAromRing = true;
                    }
                }
                else
                {
                    std::cout << "Bug, can not find atom with ID "
                              << iAt->id << " and serial number "
                              << iAt->seriNum << std::endl;
                    exit(1);
                }
                //std::cout << "Atom " << iAt->id << " add a ring sign "
                //          << iAt->ringRepS[aRepId] << std::endl;
            }
        }
    }

    void ringTools::outRingSec(AtomDict& tAtom)
    {

    }

    std::string ringTools::outRingSecStr(AtomDict& tAtom)
    {
        std::string tS="";


        return tS;
    }


}
