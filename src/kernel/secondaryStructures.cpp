/* 
 * File:   secondaryStructures.h
 * Author: flong
 *
 * Created on September 7, 2011, 6:42 PM
 */


#include "secondaryStructures.h"


namespace LIBMOL
{
    
    Helix::Helix():itsSeriNum(ZeroInt),
            itsID(NullString),
            itsInitRes(NullPoint),
            itsEndRes(NullPoint),
            itsHelixClass(0),
            itsComment(NullString),
            itsLength(ZeroInt)
    {
        
    }
    
    Helix::Helix(const Helix &tH) 
    {
        itsSeriNum    = tH.itsSeriNum;
        itsID         = tH.itsID;
        itsInitRes    = tH.itsInitRes;
        itsEndRes     = tH.itsEndRes;
        itsHelixClass = tH.itsHelixClass;
        itsComment    = tH.itsComment;
        itsLength     = tH.itsLength;
        
        for (int i =0; i < (int)tH.residues.size();
               i++ )
        {
            residues.push_back(tH.residues[i]);
        }
    }
    
    Helix::~Helix()
    {
        if (itsInitRes)
        {
            delete itsInitRes;
            itsInitRes = NullPoint;
        }
        if (itsEndRes)
        {
            delete itsEndRes;
            itsEndRes = NullPoint;
        }        
        if (!residues.empty())
        {
            residues.clear();
        }
    }
    
    Residue * Helix::getInitResidue()
    {
        return itsInitRes;
    }    
    void Helix::setInitResidue(Residue & tRes)
    {
        if( itsInitRes == NULL)
        {
            itsInitRes = new Residue(tRes);
        }
        else 
        {
            delete itsInitRes;
            itsInitRes = NullPoint;
            itsInitRes = new Residue(tRes);
        }
    }
    
    Residue * Helix::getEndResidue()
    {
        return itsEndRes;
    }
    
    void Helix::setEndResidue(Residue & tRes)
    {
        if(itsEndRes == NULL)
        {
            itsEndRes = new Residue(tRes);
        }
        else 
        {
            delete itsEndRes;
            itsEndRes = NullPoint;
            itsEndRes = new Residue(tRes);
        }
    }  
    
    Strand::Strand():itsSeriNum(ZeroInt),
            itsID(NullString),
            itsInitRes(NullPoint),
             itsEndRes(NullPoint),
            itsSheetSeriNum(ZeroInt),
            itsCurAtom(NullPoint),
            itsPrevAtom(NullPoint)
    {
    }
    
    Strand::Strand(const Strand & tStrand)
    {
        itsSeriNum       = tStrand.itsSeriNum;
        itsID            = tStrand.itsID;
        itsSeriNum       = tStrand.itsSeriNum;
        //itsInitRes(tStrand.getTerResidue(1));
        //itsEndRes(tStrand.getTerResidue(-1));
        //itsCurAtom(tStrand.getRegistAtom(1));
        //itsPrevAtom(tStrand.getRegistAtom(-1));
    }
    
    Strand::~Strand()
    {
        if (itsInitRes !=NullPoint)
        {
            delete itsInitRes;
            itsInitRes = NullPoint;
        }
        if (itsEndRes !=NullPoint)
        {
            delete itsEndRes;
            itsEndRes = NullPoint;
        }
        if (itsCurAtom != NullPoint)
        {
            delete itsCurAtom;
            itsCurAtom = NullPoint;
        }
        if (itsPrevAtom != NullPoint)
        {
            delete itsPrevAtom;
            itsPrevAtom = NullPoint;
        }
    }
 
    Residue * Strand::getTerResidue(int tN)
    {
        if (tN == 1)
        {
            return itsInitRes;
        }
        else if (tN == -1)
        {
            return itsEndRes;
        }             
        else 
        {
            return NullPoint;
        }
    }
    
    void Strand::setTerResidue(Residue& tRes, int tN)
    {
        if (tN == 1)
        {
            itsInitRes = new Residue(tRes);
        }
        else if (tN == -1)
        {
            itsEndRes  = new Residue(tRes);
        }              
    }
    
    Atom * Strand::getRegistAtom(int tN)
    {
        if (tN == 1)
        {
            return itsCurAtom;
        }
        else if (tN == -1)
        {
            return itsPrevAtom;
        }
        else 
        {
            return NullPoint;
        }
    }
    
    void Strand::setRegistAtom(Atom & tAtom, int tN)
    {
        if (tN == 1)
        {
            itsCurAtom = new Atom(tAtom);
        }
        else if (tN == -1)
        {
            itsPrevAtom = new Atom(tAtom);
        }
    }
    
    Sheet::Sheet():itsSeriNum(ZeroInt),
            itsID(NullString)
    {
    }
    
    Sheet::Sheet(const Sheet& tSheet)
    {
        itsSeriNum = tSheet.itsSeriNum;
        itsID      = tSheet.itsID;
        
        for (int i =0; i < (int)tSheet.allStrands.size();
                i++)
        {
            allStrands.push_back((tSheet.allStrands[i]));
        }
        
        for (int i=0; i < (int)tSheet.senses.size();
             i++)
        {
            senses.push_back(tSheet.senses[i]);
        }
        
    }

    Sheet::~Sheet()
    {
        if( !allStrands.empty())
        {
            allStrands.clear();
        }
        if ( !senses.empty() )
        {
            senses.clear();
        }
    }
    
    
    
            
    
}


