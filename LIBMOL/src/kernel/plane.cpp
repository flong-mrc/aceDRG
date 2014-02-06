
/* 
 * File:   plane.h
 * Author: flong
 *
 * Created on Oct. 5, 2011, 4:04 PM
 */

#include "plane.h"

namespace LIBMOL
{
    Plane::Plane():itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsValue(ZeroReal),
            itsValueSt(ZeroReal),
            itsSigValue(ZeroReal)
    {
    }

    Plane::Plane(const Plane& tP):itsName(tP.getName()),
            itsID(tP.getID()),
            itsSeriNum(tP.getSeriNum()),
            itsValue(tP.getValue(false)),
            itsValueSt(tP.getValue(true)),
            itsSigValue(tP.getSigValue())
    {
        for (int i=0; i<(int)tP.planeCeofs.size(); i++)
        {
            planeCeofs.push_back(tP.planeCeofs[i]);
        }
        
        for(int i=0; i<(int)tP.atomForceConsts.size(); i++)
        {
            atomForceConsts.push_back(tP.atomForceConsts[i]);
        }
    }
    
    Plane::~Plane()
    {
        if(!atomForceConsts.empty())
        {
            atomForceConsts.clear();
        }
        if(!planeCeofs.empty())
        {
            planeCeofs.clear();
        }
    }
    
    Name Plane::getName() const
    {
            return itsName;
    }
    void Plane::setName(Name tNa)
    {
        itsName = tNa;
    }    
    
    ID Plane::getID() const
    {
        return itsID;
    }
    void Plane::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Plane::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Plane::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }   

    REAL Plane::getValue(bool tN) const
    {
       if (tN)
       {
           return itsValueSt;
       }
       else
       {
           return itsValue;
       }
    }
    
    void Plane::setValue(REAL tV, bool tN)
    {
        if(tN)
        {
            itsValueSt = tV;
        }
        else 
        {
            itsValueSt = tV;
        }
    }
 
    REAL Plane::getSigValue() const
    {
        return itsSigValue;
    }

    void Plane::setSigValue(REAL tV)
    {
        itsSigValue = tV;
    }
    
    
    /* Another plane class*/
    PlaneDict::PlaneDict():archID(NullString),
            archPos(ZeroInt),
            fConst(0.20)
    {
    }
    
    PlaneDict::PlaneDict(const PlaneDict & tPl):archID(tPl.archID),
            archPos(tPl.archPos),
            fConst(tPl.fConst)
    {
        
        for (std::map<ID, int>::const_iterator iAt=tPl.atoms.begin();
                iAt != tPl.atoms.end(); iAt++)
        {
            atoms[iAt->first] = iAt->second;
        }
    }
    
    PlaneDict::~PlaneDict()
    {
    }  
    
    
    
}
 
        