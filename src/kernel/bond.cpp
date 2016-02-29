/* 
 * File:   bond.cpp
 * Author: flong
 *
 * Created on September 6, 2011, 2:29 PM
 */


#include "bond.h"


namespace LIBMOL
{
    Bond::Bond():isItTouched(false),
            itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsOrder(NullString),
            itsOrderSt(NullString),
            itsType(ZeroShort),
            itsLength(ZeroReal),
            itsLengthSt(ZeroReal),
            itsSigLength(ZeroReal),
            itsForceConst(ZeroReal)
    {   
    }
    
    Bond::Bond(const Bond & tBond):isItTouched(tBond.isItTouched),
            itsName(tBond.getName()),
            itsID(tBond.getID()),
            itsSeriNum(tBond.getSeriNum()),
            itsOrder(tBond.getOrder(false)),
            itsOrderSt(tBond.getOrder(true)),
            itsType(tBond.getType()),
            itsLength(tBond.getLength(false)),
            itsLengthSt(tBond.getLength(true)),
            itsSigLength(tBond.getSigLength()),
            itsForceConst(tBond.getForceConst())
    {
        for (std::vector<Atom>::const_iterator iT=tBond.atoms.begin();
                iT != tBond.atoms.end(); iT++)
        {
            atoms.push_back((*iT));
        }
          //  itsLength   = tBond.getLength(false);
          //  itsLengthSt = tBond.getLength(true);
          //  itsForceConst(tBond.getForceConst());
    }
    
    //Bond * Bond::createBond(Atom& tAtm1, Atom& tAtm2)
    //{
    //}
    
    Bond::~Bond()
    {   
        
    }
    
    void Bond::destroyBond()
    {
        itsName       = NullString;
        itsID         = NullString;
        itsSeriNum    = ZeroInt;
        itsOrder      = NullString;
        itsType       = ZeroShort;
        itsLength     = ZeroReal;
        itsLengthSt   = ZeroReal;
        itsSigLength  = ZeroReal;
        itsForceConst = ZeroReal;
        isItTouched   = false;
        
        if (!atoms.empty())
        {
            atoms.clear();
        }
    }
    
    Bond & Bond::operator = (const Bond & tB)  
    {
        itsName    = tB.getName();
        itsID      = tB.getID();
        itsSeriNum = tB.getSeriNum();
        itsOrder   = tB.getOrder(false);
        itsOrderSt = tB.getOrder(true);
        itsType    = tB.getType();
        itsSigLength = tB.getSigLength();
        isItTouched  = tB.isItTouched;
        //itsLength   = tB.getLength(false);
        //itsLengthSt = tB.getLength(true);
        for (std::vector<Atom>::const_iterator iT=tB.atoms.begin();
                iT != tB.atoms.end(); iT++)
        {
            atoms.push_back((*iT));
        }      
        return *this;
    }
    
    REAL Bond::getAtomPairLength(Atom& tA1, Atom tA2)
    {
        return distanceV(tA1.coords, tA2.coords);
    }

    Name Bond::getName() const
    {
            return itsName;
    }
    void Bond::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Bond::getID() const
    {
        return itsID;
    }
    void Bond::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Bond::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Bond::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }

    ID Bond::getOrder(bool tSt) const
    {
        if (tSt)
        {
            return itsOrderSt;
        }
        else
        {
            return itsOrder;
        }
    }

    void Bond::setOrder(ID tO, bool  tSt)
    {
        if (tSt)
        {
            itsOrderSt = tO;
        }
        else 
        {
            itsOrder  = tO;
        }
    }
    
    REAL Bond::getLength(bool tSt) const
    {
        if(tSt)
        {
            return itsLengthSt;
        }
        else 
        {
            return itsLength;
        }
    }
    
    void Bond::setLength(REAL tL, bool tS)
    {
        if(tS)
        {
            itsLengthSt = tL;
        }
        else
        {
            itsLength   = tL;
        }
    }
    
    void Bond::setLength()
    {
        // calculate length from the bond's
        // constituent atoms
        if (atoms.size() == 2)
        {
            itsLength=distanceV(atoms[0].coords, atoms[1].coords);
        }
        else
        {
            itsLength=0.0;
        }
    }
    
    // simplified version of bond 
    BondDict::BondDict(): resName(NullString),
            seriNum(ZeroInt),
            order(NullString),
            orderN(ZeroReal),
            value(ZeroReal),
            sigValue(0.02),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            valueP(ZeroReal),
            sigValueP(ZeroReal),
            oriValue(ZeroReal),
            approxLevel(ZeroInt),
            hasMetal(false),
            hasCodValue(false),
            numCodValues(ZeroInt),
            numCodValuesP(ZeroInt),
            isInSameRing(false),
            isAromatic(false)
    {
    }
    
    BondDict::BondDict(const BondDict& tBond): resName(tBond.resName),
            seriNum(tBond.seriNum),
            order(tBond.order),
            orderN(tBond.orderN),
            value(tBond.value),
            sigValue(tBond.sigValue),
            valueST(tBond.valueST),
            sigValueST(tBond.sigValueST),
            valueP(tBond.valueP),
            sigValueP(tBond.sigValueP),
            oriValue(tBond.oriValue),
            approxLevel(tBond.approxLevel),
            hasMetal(tBond.hasMetal),
            hasCodValue(tBond.hasCodValue),
            numCodValues(tBond.numCodValues),
            numCodValuesP(tBond.numCodValuesP),
            isInSameRing(tBond.isInSameRing),
            isAromatic(tBond.isAromatic)
    {
        for (std::vector<ID>::const_iterator tA = tBond.atoms.begin();
                tA != tBond.atoms.end(); tA++)
        {
            atoms.push_back(*tA);
        }
        for (std::vector<int>::const_iterator iI=tBond.atomsIdx.begin();
                iI !=tBond.atomsIdx.end(); iI++)
        {
        
            atomsIdx.push_back(*iI);
        } 
        
        for (std::vector<ID>::const_iterator tA = tBond.atomsCodClasses.begin();
                tA != tBond.atomsCodClasses.end(); tA++)
        {
            atomsCodClasses.push_back(*tA);
        }
        for (std::vector<int>::const_iterator tA=tBond.atomsHashingCodes.begin();
                tA !=tBond.atomsHashingCodes.end(); tA++)
        {
            atomsHashingCodes.push_back(*tA);
        }
        
        for (std::map<ID, ID>::const_iterator tA = tBond.atomSPs.begin();
                tA != tBond.atomSPs.end(); tA++)
        {
            atomSPs[tA->first]=tA->second;
        }
                
        for (std::vector<ID>::const_iterator tA = tBond.atomsMainRep.begin();
                tA != tBond.atomsMainRep.end(); tA++)
        {
            atomsMainRep.push_back(*tA);
        }
        for (std::vector<ID>::const_iterator tA = tBond.atomsNB2Rep.begin();
                tA != tBond.atomsNB2Rep.end(); tA++)
        {
            atomsNB2Rep.push_back(*tA);
        }
        for (std::vector<ID>::const_iterator tA = tBond.atomsNBRep.begin();
                tA != tBond.atomsNBRep.end(); tA++)
        {
            atomsNBRep.push_back(*tA);
        }
        for (std::vector<ID>::const_iterator tA = tBond.atomsNB3Rep.begin();
                tA != tBond.atomsNB3Rep.end(); tA++)
        {
            atomsNB3Rep.push_back(*tA);
        }
        for (std::vector<REAL>::const_iterator iV=tBond.codBondValues.begin();
                iV != tBond.codBondValues.end(); iV++)
        {
            codBondValues.push_back(*iV);
        }
        for(std::map<ID, int>::const_iterator iA = tBond.fullAtoms.begin();
                iA != tBond.fullAtoms.end(); iA++)
        {
            fullAtoms[iA->first] = iA->second;
        }
    }
    
    BondDict::~BondDict()
    {
        
    }
    
    /*
    std::string BondDict::bondOrderNumToStr()
    {        
    }
    */
    
    
    MetBond::MetBond() : resName(NullString),
            seriNum(ZeroInt),
            order(ZeroShort),
            value(ZeroReal),
            sigValue(0.02),
            valueST(ZeroReal),
            sigValueST(0.02),
            hasCodValue(false),
            numCodValues(ZeroInt),
            mAtomIdx(-1),
            lAtomIdx(-1)
    {
    }
    
    MetBond::MetBond(const MetBond & tBond):resName(tBond.resName),
            seriNum(tBond.seriNum),
            order(tBond.order),
            value(tBond.value),
            sigValue(tBond.sigValue),
            valueST(tBond.valueST),
            sigValueST(tBond.sigValueST),
            hasCodValue(tBond.hasCodValue),
            numCodValues(tBond.numCodValues),
            mAtomIdx(tBond.mAtomIdx),
            lAtomIdx(tBond.lAtomIdx)
    {
    }
    
    extern int getBond(std::vector<BondDict> & tBonds,
                       int tAt1, int tAt2)
    {
        
        int tBo=-1;
        
        for (int iBo=0; iBo < (int)tBonds.size(); iBo++)
        {
            if ((tBonds[iBo].atomsIdx[0] ==tAt1 && tBonds[iBo].atomsIdx[1] ==tAt2)
                || (tBonds[iBo].atomsIdx[0] ==tAt2 && tBonds[iBo].atomsIdx[1] ==tAt1))
            {
                tBo = iBo;
                break;
            }
        }
        
        return tBo;
    }
    
    extern bool checkIfBondInSameRing(std::vector<AtomDict> & tAtoms, 
                                      int tIdx1, int tIdx2)
    {
        bool lInRing = false;
        for (std::vector<int>::iterator iR=tAtoms[tIdx1].inRings.begin();
                iR !=tAtoms[tIdx1].inRings.end(); iR++)
        {
            if (std::find(tAtoms[tIdx2].inRings.begin(), tAtoms[tIdx2].inRings.end(), *iR)
                    !=tAtoms[tIdx2].inRings.end())
            {
                lInRing = true;
                break;
            }
           
        }
        return lInRing;
        
    }
    
    
}