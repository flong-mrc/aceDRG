
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
    
    extern void setAllRingPlanes(std::vector<RingDict>      & tAllRings,
                                 std::vector<AtomDict >     & tAtoms, 
                                 std::vector<PlaneDict>     & tPlanes)
    {
        
        std::vector<int> atmIdxAll;
        for (std::vector<RingDict>::iterator iR=tAllRings.begin(); 
                iR !=tAllRings.end(); iR++)
        {
            for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                        iAt !=iR->atoms.end(); iAt++)
            {
                if (std::find(atmIdxAll.begin(), atmIdxAll.end(), iAt->seriNum)
                        ==atmIdxAll.end())
                {
                    atmIdxAll.push_back(iAt->seriNum);
                }
            }
        }
        
        for (std::vector<RingDict>::iterator iR=tAllRings.begin(); 
                iR !=tAllRings.end(); iR++)
        {
            if (iR->isPlanar)
            { 
                PlaneDict aPL;
                for (std::vector<AtomDict>::iterator iAt=iR->atoms.begin();
                        iAt !=iR->atoms.end(); iAt++)
                {
                    aPL.atoms[iAt->id] = iAt->seriNum;
                    for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                            iNB !=iAt->connAtoms.end(); iNB++)
                    {
                        if (std::find(atmIdxAll.begin(), atmIdxAll.end(), *iNB)
                                ==atmIdxAll.end())
                        {
                            aPL.atoms[tAtoms[*iNB].id]=tAtoms[*iNB].seriNum;
                        }
                    }
                }
                
                tPlanes.push_back(aPL);
            }
        }
        
        // Other small planes are determined via torsion angles.
        
        // Check
        if (tPlanes.size() !=0)
        {
            std::cout << "There are " << tPlanes.size() 
                      << " planes " << std::endl;
            for (std::vector<PlaneDict>::iterator iPl=tPlanes.begin();
                     iPl !=tPlanes.end(); iPl++)
            {
                std::cout << "One Plane. It has " << iPl->atoms.size()
                          << " They are : " << std::endl;
                for (std::map<ID, int>::iterator iAt=iPl->atoms.begin();
                        iAt !=iPl->atoms.end(); iAt++)
                {
                    std::cout << "Atom " << iAt->first << std::endl;
                }
            }
        }
        
    }
    
    
    
}
 
        