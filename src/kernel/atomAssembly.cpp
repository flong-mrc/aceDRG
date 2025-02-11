/*
 * File:   atomAssembly.cpp
 * Author: flong
 *
 * Created on September  5, 2011, 11:24 PM
 */

#include "atomAssembly.h"


namespace LIBMOL
{
    AtomAssembly::AtomAssembly(): itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt)
    {
        
    }
    
   //AtomAssembly & AtomAssembly::operator=(const AtomAssembly & tAs)
   // {
   // }
    
    AtomAssembly::AtomAssembly(const AtomAssembly & tAtomAssembly):
             itsName(tAtomAssembly.getName()),
             itsID(tAtomAssembly.getID()),
             itsSeriNum(tAtomAssembly.getSeriNum())
             
    {
        for (int i=0; i < (int)tAtomAssembly.atoms.size(); i++)
        {
            atoms.push_back(tAtomAssembly.atoms[i]);
        }
    }
    
    AtomAssembly::~AtomAssembly()
    {
        if(!atoms.empty())
        {
            atoms.clear();
        }
    }
    
    Name AtomAssembly::getName() const
    {
            return itsName;
    }
    void AtomAssembly::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID AtomAssembly::getID() const
    {
        return itsID;
    }
    void AtomAssembly::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber AtomAssembly::getSeriNum() const
    {
        return itsSeriNum;
    }
    void AtomAssembly::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
    
    Size  AtomAssembly::getSize()
    {
        return atoms.size();
    }
    
    Atom & AtomAssembly::getOneAtom(int tN)
    {
        return atoms[tN-1];
    }
    
    void  AtomAssembly::addOneAtom(Atom & tAtom, unsigned int tPos)
    {
        std::vector<Atom>::iterator iA=atoms.begin();
        if (tPos <= atoms.size())
        {
           std::advance(iA, tPos-1);
           atoms.insert(iA, tAtom);
        }
        else
        {
           atoms.push_back(tAtom);
        }
    }
    
    void  AtomAssembly::deletOneAtom(SeriNumber tN)
    {
        for (std::vector<Atom>::iterator iA=atoms.begin();
                iA !=atoms.end(); iA++)
        {
            if((*iA).getSeriNum()==tN)
            {
                atoms.erase(iA);
            }
        }
    }
    
    void  AtomAssembly::swapOneAtom(Atom & tAtom, unsigned int tPos)
    {
        std::vector<Atom>::iterator iA=atoms.begin();
        if (tPos <= atoms.size())
        {
           std::advance(iA, tPos-1);
           atoms.erase(iA);
           if (tPos <=atoms.size()-1)
           {
               atoms.insert(iA, tAtom);
           }
           else 
           {
               atoms.push_back(tAtom);   
           }
        }
    }
}