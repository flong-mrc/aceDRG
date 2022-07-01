/* 
 * File:   Angle.cpp
 * Author: flong
 *
 * Created on August 9, 2011, 11:32 AM
 * 
 */

#include "angle.h"

namespace LIBMOL
{
     Angle::Angle():isItTouched(false),
             itsName(NullString),
             itsID(NullString),
             itsSeriNum(ZeroInt),
             itsValue(ZeroReal),
             itsValueSt(ZeroReal),
             itsSigValue(ZeroReal),
             itsSigValueSt(ZeroReal),
             itsFrceConst(ZeroReal)
    {
    }
    
    Angle::Angle(const Angle &tA)
    {
        isItTouched    = tA.isItTouched;
        itsName        = tA.getName();
        itsID          = tA.getID();
        itsSeriNum     = tA.getSeriNum();
        itsValue       = tA.getValue(false);
        itsValueSt     = tA.getValue(true);
        itsSigValue    = tA.getSigValue(true);
        itsSigValueSt  = tA.getSigValue(false);
        itsFrceConst   = tA.getForceConst();
        
        for (std::vector<Atom>::const_iterator iA= tA.atoms.begin();
                iA != tA.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }
    }
    
    Angle::~Angle()
    {
        if(!atoms.empty())
        {
            atoms.clear();
        }
    }
    
    Angle & Angle::operator=(const Angle & tA)
    {
        itsName        = tA.getName();
        itsID          = tA.getID();
        itsSeriNum     = tA.getSeriNum();
        itsValue       = tA.getValue(false);
        itsValueSt     = tA.getValue(true);
        itsSigValue    = tA.getSigValue(true);
        itsSigValueSt  = tA.getSigValue(false);
        itsFrceConst   = tA.getForceConst();
        isItTouched    = tA.isItTouched;
        for (std::vector<Atom>::const_iterator iA= tA.atoms.begin();
                iA != tA.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }       
        return (*this);
    }
    
    
    Name Angle::getName() const
    {
            return itsName;
    }
    void Angle::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Angle::getID() const
    {
        return itsID;
    }
    void Angle::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Angle::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Angle::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
           
    REAL Angle::getValue(bool tL) const
    {
        if (tL)
        {
            return itsValueSt;
        }
        else
        {
            return itsValue;
        }
    }
    void Angle::setValue(REAL tV, bool tL)
    {
       if(tL)
       {
           itsValueSt = tV;
       }
       else 
       {
           itsValue   = tV;
       }
    }   
    
    void Angle::setValue()
    {
        if (atoms.size() == 3)
        {
            setValue(atoms[0], atoms[1], atoms[2]);
        }
        else
        {
            itsValue =0.0;
        }
    }
    
    void Angle::setValue(Atom& tA1, Atom& tA2, Atom tA3)
    {
        std::vector<REAL> tV1, tV2;
        
        for (int i = 0; i < (int)tA1.coords.size(); i++)
        {
            tV1.push_back(tA1.coords[i]-tA2.coords[i]);
            tV2.push_back(tA2.coords[i]-tA2.coords[i]);
        }
        
        itsValue = getAngle2V(tV1, tV2);
    }    
    
    REAL Angle::getSigValue(bool tL) const
    {
        if(tL)
        {
            return itsSigValueSt;
        }
        else
        {
            return itsSigValue;
        }
    }
    void Angle::setSigValue(REAL tV, bool tL)
    {
        if(tL)
        {
            itsSigValueSt = tV;
        }
        else
        {
            itsSigValueSt = tV;
        }
    }
    
    // ============ Another angle class ==============
    
    // Default constructor 
    AngleDict::AngleDict() : seriNum(-1),
            value(ZeroReal),
            sigValue(3.00),
            valueST(ZeroReal),
            sigValueST(3.00),
            hasCodValue(false),
            numCodValues(ZeroInt),
            valueP(ZeroReal),
            sigValueP(3.00),
            numCodValuesP(ZeroInt),
            anchorID(NullString),
            anchorPos(ZeroInt),
            isFixed(false),
            approxLevel(10),
            isInSameRing(ZeroInt)
    {
    }
    
    // 
    AngleDict::AngleDict(const AngleDict & tAngle) : seriNum(tAngle.seriNum),
            value(tAngle.value),
            sigValue(tAngle.sigValue),
            valueST(tAngle.valueST),
            sigValueST(tAngle.sigValueST),
            hasCodValue(tAngle.hasCodValue),
            numCodValues(tAngle.numCodValues),
            valueP(tAngle.valueP),
            sigValueP(tAngle.sigValueP),
            numCodValuesP(tAngle.numCodValuesP),
            anchorID(tAngle.anchorID),
            anchorPos(tAngle.anchorPos),
            isFixed(tAngle.isFixed),
            approxLevel(tAngle.approxLevel),
            isInSameRing(tAngle.isInSameRing)
    {
           
        for (std::vector<int>::const_iterator iAt=tAngle.atoms.begin();
                iAt != tAngle.atoms.end(); iAt++)
        {
            atoms.push_back(*iAt);
        }
        
        for (std::vector<ID>::const_iterator iAd=tAngle.atomsId.begin();
                iAd != tAngle.atomsId.end(); iAd++)
        {
            atomsId.push_back(*iAd);
        }
        
        for(std::vector<ID>::const_iterator iID=tAngle.atomChemTypes.begin();
                iID !=tAngle.atomChemTypes.end(); iID++)
        {
            atomChemTypes.push_back(*iID);
        }
        
        for (std::vector<ID>::const_iterator iAt=tAngle.atomsCodClasses.begin();
                iAt != tAngle.atomsCodClasses.end(); iAt++)
        {
            atomsCodClasses.push_back(*iAt);
        }
        
        for (std::vector<ID>::const_iterator iAt=tAngle.atomsNB2Rep.begin();
                iAt != tAngle.atomsNB2Rep.end(); iAt++)
        {
            atomsNB2Rep.push_back(*iAt);
        }
        
        
        for (std::vector<ID>::const_iterator iAt=tAngle.atomsNBRep.begin();
                iAt != tAngle.atomsNBRep.end(); iAt++)
        {
            atomsNBRep.push_back(*iAt);
        }
        
        
        for (std::map<ID, ID>::const_iterator iAt=tAngle.atomsSPStats.begin();
                iAt !=tAngle.atomsSPStats.end(); iAt++)
        {
            atomsSPStats[iAt->first] = iAt->second;
        }
        
        for (std::map<ID, ID>::const_iterator iAt=tAngle.atomsNB1NB2SPStats.begin();
                iAt !=tAngle.atomsNB1NB2SPStats.end(); iAt++)
        {
            atomsNB1NB2SPStats[iAt->first] = iAt->second;
        }
        
        for (std::vector<REAL>::const_iterator iAt=tAngle.codAngleValues.begin();
                iAt != tAngle.codAngleValues.end(); iAt++)
        {
            codAngleValues.push_back(*iAt);
        }
    }
    
    AngleDict::AngleDict(ID tAnchID, int tAnchPos, 
            std::vector<AtomDict>& tAtoms): seriNum(-1),
            value(ZeroReal),
            sigValue(ZeroReal),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            hasCodValue(false),
            numCodValues(ZeroInt),
            anchorID(tAnchID),
            anchorPos(tAnchPos),
            isFixed(false),
            approxLevel(10),
            isInSameRing(ZeroInt)
    {
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            atoms.push_back(iAt->seriNum);
            atomsCodClasses.push_back(iAt->codClass);
        }
    }
    
    AngleDict::~AngleDict()
    {
    }
    
    
    
    /*
    void AngleDict::setValue()
    {
        std::vector<REAL> vect1, vect2;
        
        if ((int)atoms.size() == 3)
        {
            for (int i=0; i < (int)atoms[1].coords.size(); i++)
            {
                REAL t1, t2;
                t1 = atoms[0].coords[i] - atoms[1].coords[i];
                t2 = atoms[2].coords[i] - atoms[1].coords[i];
                vect1.push_back(t1);
                vect2.push_back(t2);
            }
            
            value = getAngle2V(vect1, vect2);
        }
        else
        {
            std::cout << "number of atoms constructing this angle is less than 3."
                    <<std::endl;
        }
    }
    */ 
    
    extern int getAngle(std::vector<AngleDict> & tAngs,
                        int cAt, int tAt1, int tAt2)
    {
        int tAng =-1;
        for (int iAn =0; iAn < (int)tAngs.size(); iAn++)
        {
            if (tAngs[iAn].atoms[0]==cAt)
            {
                if ((tAngs[iAn].atoms[1] ==tAt1 && tAngs[iAn].atoms[2] ==tAt2) ||
                (tAngs[iAn].atoms[1] ==tAt2 && tAngs[iAn].atoms[2] ==tAt1))
                {
                    tAng = iAn;
                }
            }
        }
        
        return tAng;
    }
    
    extern REAL getAngleValueFromFracCoords(AtomDict & tAtCen,
                                            AtomDict & tAt1,
                                            AtomDict & tAt2,
                                            REAL a, REAL b, REAL c,
                                  REAL alpha, REAL beta, REAL gamma) 
    {


        if (tAtCen.fracCoords.size() == 3 && tAtCen.fracCoords.size() == tAt1.fracCoords.size()
                && tAtCen.fracCoords.size() == tAt2.fracCoords.size()) {
            std::vector<REAL> tV1, tV2;
            for (unsigned i = 0; i < 3; i++) {
                tV1.push_back(tAt1.fracCoords[i] - tAtCen.fracCoords[i]);
                tV2.push_back(tAt2.fracCoords[i] - tAtCen.fracCoords[i]);
            }

            REAL lenTv1 = getBondLenFromFracCoords(tAtCen.fracCoords, tAt1.fracCoords,
                    a, b, c, alpha, beta, gamma);
            REAL lenTv2 = getBondLenFromFracCoords(tAtCen.fracCoords, tAt2.fracCoords,
                    a, b, c, alpha, beta, gamma);
            //std::cout << "lenTv1 " << lenTv1 << std::endl;
            //std::cout << "lenTv2 " << lenTv2 << std::endl;

            if (lenTv1 > 0.000001 && lenTv2 > 0.000001) {
                REAL coF = (a * a * tV1[0] * tV2[0] + b * b * tV1[1] * tV2[1] + c * c * tV1[2] * tV2[2]
                        + b * c * (tV1[1] * tV2[2] + tV1[2] * tV2[1]) * cos(alpha * PI180)
                        + c * a * (tV1[2] * tV2[0] + tV1[0] * tV2[2]) * cos(beta * PI180)
                        + a * b * (tV1[0] * tV2[1] + tV1[1] * tV2[0]) * cos(gamma * PI180)) / (lenTv1 * lenTv2);

                //std::cout << "coF " << coF << std::endl;
                // The float point error : ABS might > 1.0 or < -1.0
                if (coF >= 1.0) {
                    coF -= 0.0000000001;
                } else if (coF <= -1.0) {
                    coF += 0.0000000001;
                }
                return acos(coF);
            } else {
                std::cout << "There is at least one pair of atoms overlapped " << std::endl;
                return 0.0;
            }
        } else {
            std::cout << "Error: check atom coordinate dimensions " << std::endl;
            return 0.0;
        }

    }
}