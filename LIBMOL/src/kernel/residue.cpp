/* 
 * File:   residue.cpp
 * Author: flong
 *
 * Created on September 5, 2011, 18:21 AM
 */

#include "residue.h"


namespace LIBMOL
{
    class Atom;
    class AtomAssembly;
    class Chain;
    class Model;
    
    Residue::Residue() : itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsSeqNum(ZeroInt),
            itsInsCode(NullString),
            itsChainSeriNum(ZeroInt),
            itsChainID(NullString),
            itsModelSeriNum(ZeroInt),
            itsGroupID(NullString)           
    {           
    }
    
    Residue::Residue(const Residue & tRes) : itsName(tRes.getName()),
             itsID(tRes.getID()),
             itsSeriNum(tRes.getSeriNum()),  
             itsSeqNum(tRes.getSeqNum()),
             itsInsCode(tRes.getInsCode()),
             itsChainSeriNum(tRes.getChainSeriNum()),
             itsChainID(tRes.getChainID()),
             itsModelSeriNum(tRes.getModelSeriNum()),
             itsGroupID(tRes.getGroupID())
    {
        for (int i =0; i < (int)tRes.atoms.size(); i++)
        {
            atoms.push_back(tRes.atoms[i]);
        }
    }
    
    Residue::~Residue()
    {   
    }

   /* Residue & Residue::operator=(const Residue & tR)
    {
        itsName         = tR.getName();
        itsID           = tR.getID();
        itsSeriNum      = tR.getSeriNum();
        itsSeqNum       = tR.getSeqNum();
        itsChainID      = tR.getChainID();
        itsInsCode      = tR.getInsCode();
        itsGroupID      = tR.getGroupID();
        itsModelSeriNum = tR.getModelSeriNum();
  
        
        return(*this);
    } */
 
    Name Residue::getName() const
    {
            return itsName;
    }
    void Residue::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Residue::getID() const
    {
        return itsID;
    }
    void Residue::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Residue::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Residue::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
           
    SeriNumber Residue::getSeqNum() const
    {
        return itsSeqNum;
    }
    
    void Residue::setSeqNum(SeriNumber tN)
    {
        itsSeqNum = tN;
    }
    
    ID Residue::getInsCode() const
    {
        return itsInsCode;
    }
    void Residue::setInsCode(ID tIn)
    {
        itsInsCode = tIn;
    }
        
    SeriNumber Residue::getChainSeriNum() const
    {
        return itsChainSeriNum;
    }
    void  Residue::setChainSeriNum(SeriNumber tN)
    {
        itsChainSeriNum = tN;
    }
        
    ID Residue::getChainID() const
    {
        return itsChainID;
    }
    void Residue::setChainID(ID tID)
    {
        itsChainID = tID;
    }
        
    SeriNumber Residue::getModelSeriNum() const
    {
        return itsModelSeriNum;
    }
    void Residue::setModelSeriNum(SeriNumber tN)
    {
        itsModelSeriNum = tN;
    }
        

        
    ID  Residue::getGroupID() const
    {
        return itsGroupID;
    }
    void Residue::setGroupID(ID tID)
    {
        itsGroupID = tID;
    }
    
    
    // in-line set and access methods
    void Residue::addOneAtom(Atom & tAtom)
    {
        atoms.push_back(tAtom);
    }
    
    // modified residues
    
    ModRes::ModRes():itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsSeqNum(ZeroInt),
            itsInsCode(NullString),
            itsChainSeriNum(ZeroInt), 
            itsStdName(NullString),
            itsComment(NullString)
    {
    }
    
    ModRes::ModRes(const ModRes & tM)
    {
        itsName            = tM.getName();
        itsID              = tM.getID();
        itsSeriNum         = tM.getSeriNum();
        itsSeqNum          = tM.getSeqNum();
        itsInsCode         = tM.getInsCode();
        itsChainSeriNum    = tM.getChainSeriNum();
        itsStdName         = tM.getStdName();
        itsComment         = tM.getComment();
    }
    
    ModRes::~ModRes()
    {
    }
            
  
    ID ModRes::getInsCode() const
    {
        return itsInsCode;
    }
    void ModRes::setInsCode(ID tIn)
    {
        itsInsCode = tIn;
    }   
    
}
