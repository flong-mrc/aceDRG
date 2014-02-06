
/* 
 * File:   DnaRna.cpp 
 * Author: flong
 *
 * Created on December 8, 2011, 10:33 PM
 */

#include "DnaRna.h"
namespace LIBMOL
{
    DnaBase::DnaBase(): baseID(""),
            fullName(NullString),
            group(NullString),
            cID(NullString)
    {   
    }
    
    DnaBase::DnaBase(const DnaBase & tBase) : baseID(tBase.baseID),
            fullName(tBase.fullName),
            group(tBase.group),
            cID(tBase.cID)
    {
        for (std::vector<Atom>::const_iterator tA=tBase.atoms.begin();
                tA != tBase.atoms.end(); tA++)
        {
            atoms.push_back((*tA));
        }
    }
    
    DnaBase::~DnaBase()
    {
    }
    
    DnaBasePair::DnaBasePair() :  pTwist(ZeroReal),
            buckle(ZeroReal),
            opening(ZeroReal),
            shear(ZeroReal),
            stretch(ZeroReal),
            stagger(ZeroReal)
    {
    }
    
    DnaBasePair::DnaBasePair(const DnaBasePair & tBPair):  pTwist(tBPair.pTwist),
            buckle(tBPair.buckle),
            opening(tBPair.opening),
            shear(tBPair.shear),
            stretch(tBPair.stretch),
            stagger(tBPair.stagger)
    {
        for (std::vector<DnaBase>::const_iterator tB = tBPair.bases.begin();
                tB != tBPair.bases.end(); tB++)
        {
            bases.push_back((*tB));
        }
        
        for(std::vector<Bond>::const_iterator tB = tBPair.hydroBonds.begin();
                tB != tBPair.hydroBonds.end(); tB++)
        {
            hydroBonds.push_back((*tB));
        }
        
        for(std::vector<REAL>::const_iterator tA = tBPair.xyzAxis.begin();
                tA != tBPair.xyzAxis.end(); tA++)
        {
            xyzAxis.push_back((*tA));
        }
    }
    
    DnaBasePair::~DnaBasePair()
    {
    }
}
