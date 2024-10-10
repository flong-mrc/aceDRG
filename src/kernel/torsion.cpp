/*
 * File:   Torsion.cpp
 * Author: flong
 *
 * Created on September 12, 2012, 11:32 AM
 */

#include "torsion.h"

namespace LIBMOL
{

    class Atom;
    class AtomDict;
    class Bond;
    class BondDict;
    class Angle;

    Torsion::Torsion(): varFlag(ZeroInt),
            isItTouched(false),
            itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsValue(ZeroReal),
            itsValueSt(ZeroReal),
            itsSigValue(ZeroReal),
            itsValuePre(ZeroReal),
            itsForceConst(ZeroReal),
            itsPeriod(ZeroInt)

    {
    }

    Torsion::Torsion(const Torsion & tT):varFlag(tT.varFlag),
            isItTouched(tT.isItTouched),
            itsName(tT.getName()),
            itsID(tT.getID()),
            itsSeriNum(tT.getSeriNum()),
            itsValue(tT.getValue(false)),
            itsValueSt(tT.getValue(true)),
            itsSigValue(tT.getSigValue()),
            itsValuePre(tT.getValuePre()),
            itsForceConst(tT.getForceConst()),
            itsPeriod(tT.getPeriod())
    {
        for (std::vector<REAL>::const_iterator iR =tT.normalUnit.begin();
                iR != tT.normalUnit.end(); iR++)
        {
            normalUnit.push_back((*iR));
        }

        for (std::vector<int>::const_iterator iN =tT.torsUpAtomList.begin();
                iN != tT.torsUpAtomList.end(); iN++)
        {
            torsUpAtomList.push_back((*iN));
        }

        for (std::vector<Atom>::const_iterator iA =tT.atoms.begin();
                iA != tT.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }
    }


    Torsion::~Torsion()
    {
        if(!normalUnit.empty())
        {
            normalUnit.clear();
        }
        if(!torsUpAtomList.empty())
        {
            torsUpAtomList.clear();
        }
        if(!atoms.empty())
        {
            atoms.clear();
        }
    }

    Name Torsion::getName() const
    {
            return itsName;
    }
    void Torsion::setName(Name tNa)
    {
        itsName = tNa;
    }

    ID Torsion::getID() const
    {
        return itsID;
    }
    void Torsion::setID(ID tID)
    {
        itsID = tID;
    }

    SeriNumber Torsion::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Torsion::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }

    REAL Torsion::getValue(bool tB) const
    {
        if(tB)
        {
            return itsValueSt;
        }

        return itsValue;
    }
    void Torsion::setValue(REAL tT, bool tB)
    {
        if (tB)
        {
            itsValueSt = tT;
        }
        else
        {
            itsValue   = tT;
        }
    }
    REAL Torsion::setValue(Atom & tA1, Atom & tA2, Atom & tA3, Atom & tA4)
    {
        std::vector<REAL> tR1, tR2, tR3;
        /*
        std::cout << " Atom 1 " << tA1.getName()
                << " Residue " << tA1.getSeqNum()
                << " Chain " << tA1.getChainID()
                << std::endl;
        std::cout << " Atom 2 " << tA2.getName()
                 << " Residue " << tA2.getSeqNum()
                 << " Chain " << tA2.getChainID()
                 << std::endl;
        std::cout << " Atom 3 " << tA3.getName()
                << " Residue " << tA3.getSeqNum()
                << " Chain " << tA3.getChainID()
                << std::endl;
        std::cout << " Atom 4 " << tA4.getName()
                << " Residue " << tA4.getSeqNum()
                << " Chain "   << tA4.getChainID()
                << std::endl;
         */
        if(tA1.coords.size() == 3 && tA2.coords.size() == 3 &&
                tA3.coords.size() ==3)
        {
            for (int i = 0; i < (int)tA1.coords.size(); i++)
            {
                tR1.push_back(tA2.coords[i]-tA1.coords[i]);
                //std::cout << "tR1 " << i+1 << "  "
                //        <<  tA2.coords[i]-tA1.coords[i] << std::endl;
                tR2.push_back(tA3.coords[i]-tA2.coords[i]);
                //std::cout << "tR2 " << i+1 << "  "
                //        <<  tA3.coords[i]-tA2.coords[i] << std::endl;
                tR3.push_back(tA4.coords[i]-tA3.coords[i]);
                //std::cout << "tR3 " << i+1 << "  "
                //        <<  tA4.coords[i]-tA3.coords[i] << std::endl;
            }

            return getTorsion3V(tR1, tR2, tR3);
        }
        else
        {
            return 0.0;
        }

    }
    void Torsion::setValue()
    {

        if(atoms.size() == 4)
        {
            itsValue = setValue(atoms[0], atoms[1], atoms[2], atoms[3]);
        }
        else
        {
            itsValue  = 0.0;
        }
    }


    REAL Torsion::getSigValue() const
    {
        return itsSigValue;
    }
    void Torsion::setSigValue(REAL tS)
    {
        itsSigValue = tS;
    }

    REAL Torsion::getValuePre() const
    {
        return itsValuePre;
    }
    void Torsion::setValuePre(REAL tP)
    {
        itsValuePre = tP;
    }

    REAL Torsion::getForceConst() const
    {
        return itsForceConst;
    }
    void Torsion::setForceConst(REAL tF)
    {
        itsForceConst = tF;
    }

    int Torsion::getPeriod() const
    {
        return itsPeriod;
    }
    void Torsion::setPeriod(int tP)
    {
        itsPeriod = tP;
    }

    // Another class for a torsion angle
    //Default constructor
    TorsionDict::TorsionDict():seriNum(-1),
            value(ZeroReal),
            sigValue(ZeroReal),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            period(1),
            id(NullString)
    {
    }

    // Copy constructor
    TorsionDict::TorsionDict(const TorsionDict& tTorsion):seriNum(tTorsion.seriNum),
            value(tTorsion.value),
            sigValue(tTorsion.sigValue),
            valueST(tTorsion.valueST),
            sigValueST(tTorsion.sigValueST),
            period(tTorsion.period),
            id(tTorsion.id)
    {

        for (std::vector<int>::const_iterator iAt = tTorsion.atoms.begin();
                iAt != tTorsion.atoms.end(); iAt++)
        {
            atoms.push_back(*iAt);
        }

        for (std::vector<AtomDict>::const_iterator iFA=tTorsion.fullAtoms.begin();
                iFA !=tTorsion.fullAtoms.end(); iFA++)
        {
            fullAtoms.push_back(*iFA);
        }

        for (std::vector<BondDict>::const_iterator iBo = tTorsion.bonds.begin();
                iBo != tTorsion.bonds.end(); iBo++)
        {
            bonds.push_back(*iBo);
        }

        for (std::vector<ID>::const_iterator iAt = tTorsion.atomCodClasses.begin();
                iAt != tTorsion.atomCodClasses.end(); iAt++)
        {
            atomCodClasses.push_back(*iAt);
        }

        for (std::vector<REAL>::const_iterator iAt = tTorsion.codTorsionValues.begin();
                iAt != tTorsion.codTorsionValues.end(); iAt++)
        {
            codTorsionValues.push_back(*iAt);
        }
    }

    // Constructor using atoms
    TorsionDict::TorsionDict(std::vector<AtomDict>& tAtoms):seriNum(-1),
            value(ZeroReal),
            sigValue(ZeroReal),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            period(1),
            id(NullString)
    {

        for (std::vector<AtomDict>::const_iterator iFA=tAtoms.begin();
                iFA !=tAtoms.end(); iFA++)
        {
            fullAtoms.push_back(*iFA);
        }

        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            atomCodClasses.push_back(iAt->codClass);
        }

    }

    TorsionDict::~TorsionDict()
    {
    }

    extern int getTorsion(std::vector<TorsionDict> & tTors,
                          int tAt1, int tAt2, int tAt3, int tAt4)
    {
        int tTor=-1;

        for (std::vector<TorsionDict>::iterator iT=tTors.begin();
                iT != tTors.end(); iT++)
        {
            //std::cout << "Current :" << std::endl
            //          << "atom1 " << iT->atoms[0] << "\t atom2 " << iT->atoms[1]
            //          << "\t atom 3 " << iT->atoms[2] << "\t atom 4 "
            //          << iT->atoms[3] << std::endl;

            //if(std::find(iT->atoms.begin(), iT->atoms.end(), tAt1) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt2) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt3) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt4) !=iT->atoms.end())
            //{
            if ((iT->atoms[0]==tAt1 && iT->atoms[1]==tAt2 && iT->atoms[2]==tAt3 && iT->atoms[3]==tAt4)
                 || (iT->atoms[0]==tAt4 && iT->atoms[1]==tAt3 && iT->atoms[2]==tAt2 && iT->atoms[3]==tAt1))
            {
                tTor = iT->seriNum;
                // std::cout << "Torsion seriNumb " << iT->seriNum << std::endl;
                break;
            }

        }


        if (tTor ==-1)
        {
            std::cout << "can not find the target torsion "  << std::endl;
        }

        return tTor;
    }


    extern REAL getTorsion(std::vector<AtomDict> & tAtoms,
                           int iCur, int iNext, std::vector<int> tDoneSet)
    {

        int iPrev = tAtoms[iCur].tree["parent"][0];
        std::vector<int> tDoneCAs;
        for (std::vector<int>::iterator iC=tAtoms[iCur].connAtoms.begin();
                iC != tAtoms[iCur].connAtoms.end(); iC++)
        {
            if (*iC !=iNext && *iC !=iPrev
                && std::find(tDoneSet.begin(), tDoneSet.end(), *iC) !=tDoneSet.end())
            {
                tDoneCAs.push_back(*iC);
            }
        }

        REAL retV =0.0;

        if ((int)tDoneCAs.size() ==0)
        {
            retV=tAtoms[iNext].treeTorsion;
        }

        return retV;

    }

    extern REAL getTorsion(AtomDict & tA1,
                           AtomDict & tA2,
                           AtomDict & tA3,
                           AtomDict & tA4)
    {
        std::vector<REAL> tR1, tR2, tR3;
        /*
        std::cout << " Atom 1 " << tA1.getName()
                << " Residue " << tA1.getSeqNum()
                << " Chain " << tA1.getChainID()
                << std::endl;
        std::cout << " Atom 2 " << tA2.getName()
                 << " Residue " << tA2.getSeqNum()
                 << " Chain " << tA2.getChainID()
                 << std::endl;
        std::cout << " Atom 3 " << tA3.getName()
                << " Residue " << tA3.getSeqNum()
                << " Chain " << tA3.getChainID()
                << std::endl;
        std::cout << " Atom 4 " << tA4.getName()
                << " Residue " << tA4.getSeqNum()
                << " Chain "   << tA4.getChainID()
                << std::endl;
         */
        if(tA1.coords.size() == 3 && tA2.coords.size() == 3 &&
                tA3.coords.size() ==3)
        {
            for (int i = 0; i < (int)tA1.coords.size(); i++)
            {
                tR1.push_back(tA2.coords[i]-tA1.coords[i]);
                //std::cout << "tR1 " << i+1 << "  "
                //        <<  tA2.coords[i]-tA1.coords[i] << std::endl;
                tR2.push_back(tA3.coords[i]-tA2.coords[i]);
                // std::cout << "tR2 " << i+1 << "  "
                //          <<  tA3.coords[i]-tA2.coords[i] << std::endl;
                tR3.push_back(tA4.coords[i]-tA3.coords[i]);
                // std::cout << "tR3 " << i+1 << "  "
                //          <<  tA4.coords[i]-tA3.coords[i] << std::endl;
            }

            return getTorsion3V(tR1, tR2, tR3);
        }
        else
        {
            return 0.0;
        }

        return 0.0;
    }

    extern bool checkATorsAtomsInPla(std::vector<int> & tAtms,
                                     std::vector<AtomDict>  & tAllAtoms,
                                     std::vector<PlaneDict> & tAllPlanes)
    {

        for (std::vector<PlaneDict>::iterator iP=tAllPlanes.begin();
                        iP !=tAllPlanes.end(); iP++)
        {
            std::vector<ID> inPlAtoms;
            for (std::vector<int>::iterator iA=tAtms.begin();
                    iA !=tAtms.end(); iA++)
            {
                std::map<ID, int>::iterator iFind =iP->atoms.find(tAllAtoms[*iA].id);
                if (iFind !=iP->atoms.end())
                {
                    inPlAtoms.push_back(tAllAtoms[*iA].id);
                }
            }

            if ((int)inPlAtoms.size() == (int)tAtms.size())
            {
                return true;
            }
        }


        return false;
    }

    extern int  checkATorsAtomsInAroRing(int tAtm1, int tAtm2,
                                         std::vector<AtomDict>  & tAllAtoms,
                                         std::vector<BondDict> & tAllBonds)
    {
        //std::cout << "Check arom for atom "
        //          << allAtoms[tAtm1].id << " of sp "<< allAtoms[tAtm1].bondingIdx
        //          << " and " << allAtoms[tAtm2].id << " of sp "
        //          << allAtoms[tAtm2].bondingIdx << std::endl;
        int aRet  = 0;                      // not in any ring
        if (tAllAtoms[tAtm1].bondingIdx==2
            && tAllAtoms[tAtm2].bondingIdx==2)
        {

            int aBIdx = getBond(tAllBonds, tAtm1, tAtm2);
            //std::cout << "bond " << aBIdx << std::endl;
            //std::cout << "atom " << allBonds[aBIdx].atoms[0] << " and "
            //          << allBonds[aBIdx].atoms[1] << std::endl;
            //std::cout << "Its order " << allBonds[aBIdx].order << std::endl;
            if (aBIdx !=-1)
            {
                if (tAllBonds[aBIdx].isInSameRing)
                {
                    //std::cout << "It is in the same ring" << std::endl;

                    if (tAllBonds[aBIdx].order.find("arom") !=tAllBonds[aBIdx].order.npos)
                    {
                        aRet = 3;
                    }
                    else
                    {
                        aRet = 2;
                    }
                }
                else
                {
                    aRet=1;
                }
            }
        }
        return aRet;
    }

    extern void fixTorIDs(std::vector<TorsionDict> & tAllTorsions,
                          std::vector<AtomDict>  & tAllAtoms,
                          std::vector<BondDict>  & tAllBonds,
                          std::vector<PlaneDict> & tAllPlanes,
                          bool                   & tLMdPls)
    {

        int idxTors  = 1, idxPTors = 1, idxSp3Sp3=1, idxSp2Sp3=1, idxSp2Sp2=1;
        //std::cout << "There are " << (int)allTorsions.size()
        //          << "Torsions" << std::endl;


        for (std::vector<TorsionDict>::iterator iT=tAllTorsions.begin();
                        iT !=tAllTorsions.end(); iT++)
        {
            //std::cout << "look at torsion " << iT->seriNum << std::endl;
            //std::cout << "atom1 " << allAtoms[iT->atoms[0]].id
            //          << "  atom2 " << allAtoms[iT->atoms[1]].id
            //          << "  atom3 " << allAtoms[iT->atoms[2]].id
            //          << "  atom4 " << allAtoms[iT->atoms[3]].id << std::endl;

            int aFlag =checkATorsAtomsInAroRing(iT->atoms[1], iT->atoms[2],
                                                tAllAtoms, tAllBonds);
            //std::cout << "aFlag == " << aFlag << std::endl;
            if (aFlag == 3)
            {
                // in a aromatic ring
                if (tLMdPls)
                {
                    iT->id = "sp2_sp2_" + IntToStr(idxPTors);
                    iT->sigValue =1.0;
                }
                else
                {
                    iT->id = "const_sp2_sp2_" + IntToStr(idxPTors);
                    iT->sigValue =0;
                    if (iT->id.size() >=16 )
                    {
                        // iT->id = IntToStr(idxPTors);
                        iT->id = "const_" + IntToStr(idxPTors);
                    }
                }


                //iT->id = "P_sp2_sp2_" + IntToStr(idxPTors);

                iT->period   =1;
                idxPTors +=1;
            }
            else if (aFlag == 2)
            {
                // in an all-sp2 ring
                if (tLMdPls)
                {
                    iT->id = "sp2_sp2_" + IntToStr(idxPTors);
                    iT->sigValue =1.0;
                }
                else
                {
                    iT->id = "sp2_sp2_" + IntToStr(idxPTors);
                    iT->sigValue =5.0;
                    //if (iT->id.size() >=16 )
                    //{
                        // iT->id = IntToStr(idxPTors);
                        // iT->id = "const_" + IntToStr(idxPTors);
                    //}
                    // std::cout << "Here 2" << std::endl;
                }
                /*
                iT->id = "sp2_sp2_" + IntToStr(idxPTors);
                //if (iT->id.size() >=16 )
                //{
                    //iT->id = "const_" + IntToStr(idxPTors);
                //}
                //iT->id = "P_sp2_sp2_" + IntToStr(idxPTors);
                iT->sigValue =1.0;
                */
                iT->period   =1;
                idxPTors +=1;
            }
            else if (aFlag == 1)
            {
                iT->id = "sp2_sp2_" + IntToStr(idxPTors);
                int aBIdx12 = getBond(tAllBonds, iT->atoms[0], iT->atoms[1]);
                int aBIdx34 = getBond(tAllBonds, iT->atoms[2], iT->atoms[3]);
                //std::cout << "Flag1" << std::endl;
                //std::cout << "aBIdx12 " << aBIdx12 << std::endl;
                //std::cout << "aBIdx34 " << aBIdx34 << std::endl;
                if (aBIdx12 !=-1 && aBIdx34 !=-1)
                {

                    if (tAllBonds[aBIdx12].isInSameRing || tAllBonds[aBIdx34].isInSameRing)
                    {
                        iT->sigValue =20.0;
                    }
                    else
                    {
                        iT->sigValue =5.0;
                    }
                }
                else
                {
                    iT->sigValue =5.0;
                }
                //if (iT->id.size() >=16 )
                //{
                    //iT->id = "const_" + IntToStr(idxPTors);
                //}
                //iT->id = "P_sp2_sp2_" + IntToStr(idxPTors);

                idxPTors +=1;
            }
            else if(checkATorsAtomsInPla(iT->atoms, tAllAtoms, tAllPlanes))
            {
                iT->id = "sp2_sp2_" + IntToStr(idxPTors);

                //if (iT->id.size() >=16 )
                //{
                //    iT->id = "const_" + IntToStr(idxPTors);
                //}
                //iT->id = "P_sp2_sp2_" + IntToStr(idxPTors);
                iT->sigValue =5.0;
                idxPTors +=1;
            }
            else
            {
                if (iT->period ==3 )
                {
                    iT->id = "sp3_sp3_" + IntToStr(idxSp3Sp3);
                    idxSp3Sp3+=1;
                }
                else if (iT->period == 6)
                {
                    iT->id= "sp2_sp3_"+IntToStr(idxSp2Sp3);
                    idxSp2Sp3+=1;
                }
                else if (iT->period == 2)
                {
                    if ((tAllAtoms[iT->atoms[1]].chemType!="O"  && tAllAtoms[iT->atoms[2]].chemType !="O" )
                       && (tAllAtoms[iT->atoms[1]].chemType!="S"  && tAllAtoms[iT->atoms[2]].chemType !="S"))
                    {
                        iT->id = "sp2_sp2_"+IntToStr(idxSp2Sp2);
                        idxSp2Sp2+=1;
                    }
                    else
                    {
                        iT->id = "other_tor_"+IntToStr(idxSp2Sp2);
                        iT->sigValue =20.0;
                        idxTors++;
                    }
                    std::cout << "Torsion label " << iT->id << std::endl;

                }
                else
                {
                    iT->id = "other_tor_" + IntToStr(idxTors);
                    iT->sigValue =20.0;
                    idxTors++;
                }
            }
            //std::cout << "its ID now is " << iT->id << "and sig is "
            //          <<   iT->sigValue << std::endl;
        }

    }

    extern void setupMiniTorsions(std::vector<TorsionDict> & tAllTorsions,
                                  std::vector<AtomDict>    & tAtoms,
                                  std::vector<BondDict>    & tBonds,
                                  std::vector<TorsionDict> & tMiniTorsions)
    {
        std::map<ID, std::vector<TorsionDict> >  TorsionSetOneBond;
        std::map<ID, std::vector<int>  >         TorAtms;
        for (std::vector<TorsionDict>::iterator iTor=tAllTorsions.begin();
                iTor !=tAllTorsions.end(); iTor++)
        {
            if (iTor->atoms.size() ==4)
            {
                int idxB = getBond(tBonds, iTor->atoms[1], iTor->atoms[2]);
                if (idxB > -1)
                {
                    StrUpper(tBonds[idxB].order);

                    std::vector<int> tIdxB;
                    tIdxB.push_back(iTor->atoms[1]);
                    tIdxB.push_back(iTor->atoms[2]);
                    std::sort(tIdxB.begin(), tIdxB.end());

                    ID tLab = IntToStr(tIdxB[0]) + "_" + IntToStr(tIdxB[1]);

                    //std::cout << "Bond of atoms " << tAtoms[iTor->atoms[1]].id
                    //          << " and " << tAtoms[iTor->atoms[2]].id << std::endl;
                    //std::cout << "torsion lab " << tLab << std::endl;
                    //std::cout << "torsion " << iTor->seriNum
                    //          << " is included " << std::endl;
                    TorsionSetOneBond[tLab].push_back(*iTor);

                }
                else
                {
                    std::cout << "Error in Setup Torsion section: "
                              << "Can not find the bond between atom "
                              << tAtoms[iTor->atoms[1]].id  << " and "
                              << tAtoms[iTor->atoms[2]].id << std::endl;
                }


            }
        }

        for (std::map<ID, std::vector<TorsionDict> >::iterator iTorsB=TorsionSetOneBond.begin();
                iTorsB !=TorsionSetOneBond.end(); iTorsB++)
        {
            std::vector<ID>  idxs;
            StrTokenize(iTorsB->first, idxs, '_');
            int idxN1 =StrToInt(idxs[0]);
            int idxN2 =StrToInt(idxs[1]);
            selectOneTorFromOneBond(iTorsB->first, iTorsB->second, tAllTorsions, tAtoms, tMiniTorsions);
            // The following codes taken out on 23/09/2024, basically
            // cancelled sp2-sp2 torsion angles from one per branch to one torsion per
            // bond. For torsion angles of other kinds (e.g. sp2-sp3 or sp3-sp3)
            // always one torsion per bond.
            /*
            if (tAtoms[idxN1].bondingIdx==2
                && tAtoms[idxN2].bondingIdx==2)
            {
                std::vector<int> idx1, idx4;
                //std::cout << "lab : " << iTorsB->first << std::endl;
                for (std::vector<TorsionDict>::iterator
                     iTor  = iTorsB->second.begin();
                     iTor != iTorsB->second.end(); iTor++)
                {
                    int atm1 = iTor->atoms[0];
                    int atm4 = iTor->atoms[3];
                    if( (std::find(idx1.begin(), idx1.end(), atm1)==idx1.end())
                         &&
                        (std::find(idx4.begin(), idx4.end(), atm4)==idx4.end())
                    )
                    {
                        tMiniTorsions.push_back(*iTor);
                        idx1.push_back(atm1);
                        idx4.push_back(atm4);
                        //std::cout << "1st atom " << tAtoms[atm1].id << std::endl;
                    }
                }
            }
            else
            {
                selectOneTorFromOneBond(iTorsB->first, iTorsB->second, tAllTorsions, tAtoms, tMiniTorsions);
            }
            */
        }

        std::cout << "Total number of torsions is " << tAllTorsions.size()
                  << std::endl;

        std::cout << "Number of torsions in the mini-set "
                  << tMiniTorsions.size() << std::endl;

        for (std::vector<TorsionDict>::iterator iTor = tMiniTorsions.begin();
                iTor != tMiniTorsions.end(); iTor++)
        {
            std::cout << "The torsion selected for the bond of atom "
                      << tAtoms[iTor->atoms[1]].id << " and "
                      << tAtoms[iTor->atoms[2]].id << std::endl
                      << " is " << iTor->value << std::endl;
            std::cout << "The other two atoms are atom "
                      << tAtoms[iTor->atoms[0]].id
                      << " and " << tAtoms[iTor->atoms[3]].id
                      << std::endl;
        }
    }


    extern void selectOneTorFromOneBond(ID tS, std::vector<TorsionDict> & tTorsB,
                                        std::vector<TorsionDict> & tAllTorsions,
                                        std::vector<AtomDict> & tAtoms,
                                        std::vector<TorsionDict> & tMiniTorsions)
    {
        bool lDone = false;

        if (tTorsB.size() >0)
        {
            std::vector<ID> tLabs;
            StrTokenize(tS, tLabs, '_');

            if (tLabs.size()==2)
            {
                int idx1 = StrToInt(tLabs[0]);
                int idx2 = StrToInt(tLabs[1]);
                std::cout << "torsion lab " << tS << std::endl;
                std::cout << "Select torsion angles for the bond of atoms "
                          << tAtoms[idx1].id << " and " << tAtoms[idx2].id
                          << std::endl;

                //std::cout << "1. Number of torsions is " << tTorsB.size() << std::endl;

                std::vector<int> idxR1, idxNonH1, idxH1, idxR2, idxNonH2, idxH2;

                for(std::vector<TorsionDict>::iterator iTor=tTorsB.begin();
                        iTor !=tTorsB.end(); iTor++)
                {
                    for (int i=1; i < 3; i++)
                    {
                        if (iTor->atoms[i]==idx1)
                        {
                            for (std::vector<int>::iterator iConn=tAtoms[idx1].connAtoms.begin();
                                   iConn != tAtoms[idx1].connAtoms.end(); iConn++)
                            {
                                if ( *iConn != idx2)
                                {
                                    if (tAtoms[*iConn].inRings.size() !=0)
                                    {
					                    if (std::find(idxR1.begin(), idxR1.end(), *iConn)==idxR1.end())
				                        {
					                        idxR1.push_back(*iConn);
					                    }
                                    }
                                    else if (tAtoms[*iConn].chemType !="H")
                                    {
					                    if (std::find(idxNonH1.begin(), idxNonH1.end(), *iConn)==idxNonH1.end())
					                    {
                                            idxNonH1.push_back(*iConn);
                                        }
                                    }
                                    else
                                    {
					                    if (std::find(idxH1.begin(), idxH1.end(), *iConn)==idxH1.end())
					                    {
                                            idxH1.push_back(*iConn);
					                    }
                                    }
                                }
                            }
                        }
                        else if (iTor->atoms[i]==idx2)
                        {
                            for (std::vector<int>::iterator iConn=tAtoms[idx2].connAtoms.begin();
                                 iConn != tAtoms[idx2].connAtoms.end(); iConn++)
                            {
                                if ( *iConn != idx1)
                                {
                                    if (tAtoms[*iConn].inRings.size() !=0)
                                    {

					                    if (std::find(idxR2.begin(), idxR2.end(), *iConn)==idxR2.end())
				                        {
                                            idxR2.push_back(*iConn);
					                    }
                                    }
                                    else if (tAtoms[*iConn].chemType.find("H")
                                             ==std::string::npos)
                                    {
					                    if (std::find(idxNonH2.begin(), idxNonH2.end(), *iConn)==idxNonH2.end())
					                    {
                                            idxNonH2.push_back(*iConn);
					                    }
                                    }
                                    else
                                    {
					                    if (std::find(idxH2.begin(), idxH2.end(), *iConn)==idxH2.end())
					                    {
                                            idxH2.push_back(*iConn);
					                    }
                                    }
                                }
                            }
                        }
                    }
                }

                std::cout << "atom  " << tAtoms[idx1].id
                          << " connects to " << std::endl;
                //std::cout << "ring atoms " << idxR1.size() << std::endl
		        std::cout << "non-H atoms (excluded above) " << idxNonH1.size() << std::endl
                          << "H atoms " << idxH1.size() << std::endl;
                std::cout << "atom  " << tAtoms[idx2].id
                          << " connects to " << std::endl;
                std::cout << "ring atoms " << idxR2.size() << std::endl
                          << "non-H atoms (excluded above) " << idxNonH2.size() << std::endl
                          << "H atoms " << idxH2.size() << std::endl;

                // Now select the torsion angles
                // 1. atoms of Non-ring, non-H first
                // std::cout << "2. Number of torsions is " << tTorsB.size() << std::endl;
                if (idxNonH1.size() !=0 && idxNonH2.size() !=0)
                {
                    int tIdxT= getTorsion(tAllTorsions, idxNonH1[0], idx1, idx2, idxNonH2[0]);
                    if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                    {
                        tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                        std::cout << "added 1 " << std::endl;
                        lDone = true;
                    }
                    else
                    {
                        std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxNonH1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxNonH2[0]].id
                                      << std::endl;
                    }
                }
                else if (idxNonH1.size() !=0 && idxNonH2.size() ==0)
                {

                    if (idxR2.size() !=0)
                    {
                        int tIdxT= getTorsion(tAllTorsions, idxNonH1[0], idx1, idx2, idxR2[0]);

                        if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                        {
                           tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                           std::cout << "Add the torsion of atoms: "
                                      << tAtoms[idxNonH1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxR2[0]].id
                                      << std::endl;
                           std::cout << "added 2 " << std::endl;
                           lDone = true;
                        }
                        else
                        {
                            std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxNonH1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxR2[0]].id
                                      << std::endl;
                        }
                    }
                    else if (idxH2.size() !=0)
                    {
                        int tIdxT= getTorsion(tAllTorsions, idxNonH1[0], idx1, idx2, idxH2[0]);
                        if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                        {
                            tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                            std::cout << "added 3 " << std::endl;
                            lDone = true;
                        }
                        else
                        {
                            std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxNonH1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxH2[0]].id
                                      << std::endl;
                        }
                    }
                }
                else if (idxNonH1.size() ==0 && idxNonH2.size() !=0)
                {

                    if (idxR1.size() !=0)
                    {

                        int tIdxT= getTorsion(tAllTorsions, idxR1[0], idx1, idx2, idxNonH2[0]);

                        if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                        {
                           tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                           std::cout << "added 4 " << std::endl;
                           lDone = true;
                        }
                        else
                        {
                            std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxR1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxNonH2[0]].id
                                      << std::endl;
                        }
                    }
                    else if (idxH1.size() !=0)
                    {
                        int tIdxT= getTorsion(tAllTorsions, idxH1[0], idx1, idx2, idxNonH2[0]);
                        if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3] )
                        {
                            tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                            std::cout << "added 5 " << std::endl;
                            lDone = true;
                        }
                        else
                        {
                            std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxH1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxNonH2[0]].id
                                      << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << "can not find the first atom of the torsion "
                                      << " of atoms:  "  << std::endl
                                      << tAtoms[idx1].id << ", and "
                                      << tAtoms[idx2].id
                                      << tAtoms[idxNonH2[0]].id
                                      << std::endl;
                    }
                }
                else if (idxR1.size() !=0 || idxR2.size() !=0)
                {
                    if (idxR1.size() !=0 && idxR2.size() !=0)
                    {
                        int tIdxT= getTorsion(tAllTorsions, idxR1[0], idx1, idx2, idxR2[0]);

                        if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                        {
                           tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                           std::cout << "added 6 " << std::endl;
                           lDone = true;
                        }
                        else
                        {
                            std::cout << "Error: can not find the torsion of atoms: "
                                      << tAtoms[idxR1[0]].id << ", "
                                      << tAtoms[idx1].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idxR2[0]].id
                                      << std::endl;
                        }
                    }
                    else if (idxR1.size() !=0 && idxR2.size()==0)
                    {
                        if (idxH2.size() !=0)
                        {
                            int tIdxT= getTorsion(tAllTorsions, idxR1[0], idx1, idx2, idxH2[0]);
                            if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                            {
                                tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                                std::cout << "added 7 " << std::endl;
                                lDone = true;
                            }
                            else
                            {
                                std::cout << "Error: can not find the torsion of atoms: "
                                          << tAtoms[idxR1[0]].id << ", "
                                          << tAtoms[idx1].id << ", "
                                          << tAtoms[idx2].id << ", and "
                                          << tAtoms[idxH2[0]].id
                                          << std::endl;
                            }
                        }
                        else
                        {
                            std::cout << "can not find the fourth atom of the torsion "
                                      << std::endl << " consisting of atoms: "
                                      << tAtoms[idxR1[0]].id << ", "
                                      << tAtoms[idx1].id << ", and "
                                      << tAtoms[idx2].id << std::endl;
                        }
                    }
                    else if (idxR2.size() !=0 && idxR1.size()==0)
                    {

                        if (idxH1.size() !=0)
                        {
                            int tIdxT= getTorsion(tAllTorsions, idxH1[0], idx1, idx2, idxR2[0]);
                            if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                            {
                                tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                                std::cout << "added 8 " << std::endl;
                                lDone = true;
                            }
                            else
                            {
                                std::cout << "Error: can not find the torsion of atoms: "
                                          << tAtoms[idxH1[0]].id << ", "
                                          << tAtoms[idx1].id << ", "
                                          << tAtoms[idx2].id << ", and "
                                          << tAtoms[idxR2[0]].id
                                          << std::endl;
                            }
                        }
                        else
                        {
                            std::cout << "can not find the fourth atom of the torsion "
                                      << std::endl << " consisting of "
                                      << tAtoms[idxR2[0]].id << ", "
                                      << tAtoms[idx2].id << ", and "
                                      << tAtoms[idx1].id << std::endl;
                        }
                    }
                }
                else
                {

                    std::cout << "Both atom "
                              << tAtoms[idx1].id
                              << ", and " << tAtoms[idx2].id
                              << " connect to H atoms only. Check the input structure ! "
                              << std::endl;
                    if (tAtoms[idx1].connAtoms.size() > 1 && tAtoms[idx2].connAtoms.size() > 1)
                    {
                        int idxC1 = -1, idxC2 = -1;
                        for (unsigned idxH1=0; idxH1 < tAtoms[idx1].connAtoms.size();
                                idxH1++)
                        {
                            if ( tAtoms[idx1].connAtoms[idxH1] != idx2)
                            {
                                idxC1 = tAtoms[idx1].connAtoms[idxH1];
                                break;
                            }

                        }
                        for (unsigned idxH2=0; idxH2 < tAtoms[idx2].connAtoms.size();
                                idxH2++)
                        {
                            if ( tAtoms[idx2].connAtoms[idxH2] != idx1)
                            {
                                idxC2 = tAtoms[idx2].connAtoms[idxH2];
                                break;
                            }

                        }
                        if (idxC1 !=-1 && idxC2 !=-1)
                        {
                            int tIdxT= getTorsion(tAllTorsions, idxC1, idx1, idx2, idxC2);
                            if (tIdxT !=-1 && tAllTorsions[tIdxT].atoms[0] != tAllTorsions[tIdxT].atoms[3])
                            {
                                tMiniTorsions.push_back(tAllTorsions[tIdxT]);
                                std::cout << "added 8 " << std::endl;
                                lDone = true;
                            }
                            else
                            {
                                std::cout << "Error: can not find the torsion of atoms: "
                                          << tAtoms[idxC1].id << ", "
                                          << tAtoms[idx1].id << ", "
                                          << tAtoms[idx2].id << ", and "
                                          << tAtoms[idxC2].id
                                          << std::endl;
                            }
                        }
                        else
                        {
                            std::cout << "No torsion angle is selected into the mini set"
                                      << std::endl;
                        }
                    }
                }
            }
        }

        std::cout << "3. Number of torsions is " << tMiniTorsions.size() << std::endl;

        if (tMiniTorsions.size() > 0)
        {
            std::vector<TorsionDict> tmpTorsions;
            std::cout << "size of tMiniTorsions " << tMiniTorsions.size() << std::endl;
            int i=1;
            for (std::vector<TorsionDict>::iterator iTor=tMiniTorsions.begin();
                iTor != tMiniTorsions.end(); iTor++)
            {
                std::cout << "Tor " << i << std::endl;
                //std::cout << " iTor->atoms[0] " << iTor->atoms[0] << std::endl
                //          << " iTor->atoms[1] " << iTor->atoms[1]
                //          << " iTor->atoms[2] " << iTor->atoms[2] << std::endl
                //          << " iTor->atoms[3] " << iTor->atoms[3] << std::endl;

		        if (iTor->atoms.size()==4 && tAtoms[iTor->atoms[0]].id !=tAtoms[iTor->atoms[3]].id )
		        {
                    std::cout << "atom 1 " << tAtoms[iTor->atoms[0]].id
			              << " atom 2 " << tAtoms[iTor->atoms[1]].id
			              << " atom 3 " << tAtoms[iTor->atoms[2]].id
			              << " atom 4 " << tAtoms[iTor->atoms[3]].id << std::endl;
                    tmpTorsions.push_back(*iTor);
                    i++;
		        }

            }

            std::cout << "4. Number of torsions is " << tmpTorsions.size() << std::endl;

            tMiniTorsions.clear();

            for (std::vector<TorsionDict>::iterator iTor=tmpTorsions.begin();
                iTor != tmpTorsions.end(); iTor++)
            {
		        if (iTor->atoms.size()==4 && tAtoms[iTor->atoms[0]].id !=tAtoms[iTor->atoms[3]].id )
		        {

                    if (iTor->atoms.size()==4 &&
                        tAtoms[iTor->atoms[0]].id !=tAtoms[iTor->atoms[3]].id )
                    {
                        tMiniTorsions.push_back(*iTor);
                    }
		        }
            }
        }

        std::cout << "5. Number of torsions is " << tAllTorsions.size() << std::endl;
        if (tAllTorsions.size() > 0)
        {

            std::vector<TorsionDict> tmpTorsions;

            for (std::vector<TorsionDict>::iterator iTor=tAllTorsions.begin();
                iTor != tAllTorsions.end(); iTor++)
            {
                tmpTorsions.push_back(*iTor);
            }

            tAllTorsions.clear();

            for (std::vector<TorsionDict>::iterator iTor=tmpTorsions.begin();
                iTor != tmpTorsions.end(); iTor++)
            {
                if (iTor->atoms.size()==4 &&
                    tAtoms[iTor->atoms[0]].id !=tAtoms[iTor->atoms[3]].id )
                {
                    tAllTorsions.push_back(*iTor);
                }
            }
        }
    }


}
