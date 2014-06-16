/* 
 * File:   chain.cpp
 * Author: flong
 *
 * Created on September 6, 2011, 10:29 AM
 */

#include "chain.h"

namespace LIBMOL
{
    Chain::Chain() :curAtomSeriNum(ZeroInt),
            itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsModSeriNum(ZeroInt)
            
    {
    }
    
    Chain::Chain(const Chain & tC) : curAtomSeriNum(tC.curAtomSeriNum),
            itsName(tC.getName()),
            itsID(tC.getID()),
            itsSeriNum(tC.getSeriNum()),
            itsModSeriNum(ZeroInt)
    {
        for (int i = 0; i < (int)tC.residues.size();
                i++)
        {
            residues.push_back(tC.residues[i]);
        }
    }
    

      
    Chain::~Chain()
    {
        if(!residues.empty())
        {
            residues.clear(); 
        }
    }
 
    Chain & Chain::operator=(const Chain & tC)
    {
        itsName = tC.getName();
        itsID = tC.getID();
        itsSeriNum = tC.getSeriNum();
        for (int i =0; i < (int)tC.residues.size(); i++)
        {
            residues.push_back(tC.residues[i]);
        }
        return (*this);
    }
 
    Name Chain::getName() const
    {
            return itsName;
    }
    void Chain::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Chain::getID() const
    {
        return itsID;
    }
    void Chain::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Chain::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Chain::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
       
    
    SeriNumber Chain::getModSeriNum() const
    {
        return itsModSeriNum;
    }
    void  Chain::setModSeriNum(SeriNumber tN)
    {
        itsModSeriNum = tN;
    }
    
    Size  Chain::getNumOfRes() const
    {
        return itsNumOfResidues;
    }
    void Chain::setNumOfRes(int tN)
    {
        itsNumOfResidues = tN;
    }
    
    
    void  Chain::addOneResidue(Residue& tR)
    {
        residues.push_back(tR);
    }
    
    void  Chain::addOneResidue(Residue& tR, unsigned int tPos)
    {
        std::vector<Residue>::iterator iT=residues.begin(); 
        
        std::advance(iT, tPos-1);
        
        residues.insert(iT, tR);
    }
    
    void  Chain::deleteOneResidue(unsigned int tPos) 
    {
        std::vector<Residue>::iterator iT=residues.begin();
        
        std::advance(iT, tPos-1);
        
        residues.erase(iT);
    }
    
    void Chain::deleteOneResidue(ID tID, SeriNumber tN, Name tName)
    {
        for(std::vector<Residue>::iterator iT=residues.begin();
                iT != residues.end(); iT++ )
        {
            if (tID==iT->getID() && tN==iT->getSeriNum() 
                    && tName==iT->getName())
                residues.erase(iT);
        }   
    }
    
    void Chain::swapOneResidue(unsigned int tPos, Residue & tR)
    {
        std::vector<Residue>::iterator iT=residues.begin();
        
        std::advance(iT, tPos-1);
        
        residues.erase(iT);
        
        --iT;
        
        residues.insert(iT, tR);
    }
    
 
}