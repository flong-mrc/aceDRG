/* 
 * File:   ssbond.cpp
 * Author: flong
 *
 * Created on September 6, 2011, 5:01 PM
 */



#include "ssbond.h"



namespace LIBMOL
{
    SSBond::SSBond():itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsLength(ZeroReal)
    {
        resSeqNums[0] =ZeroInt;
        resSeqNums[1] =ZeroInt;
        resSym[0]     =ZeroInt;
        resSym[0]     =ZeroInt;
    }
    
    SSBond::SSBond(const SSBond & tSSBond):itsName(tSSBond.getName()),
            itsID(tSSBond.getID()),
            itsSeriNum(tSSBond.getSeriNum()),
            itsLength(tSSBond.itsLength)
    {
        resSym[0] = tSSBond.resSym[0];
        resSym[1] = tSSBond.resSym[1];
        resSeqNums[0] = tSSBond.resSeqNums[0];
        resSeqNums[1] = tSSBond.resSeqNums[1];
        
        if (tSSBond.residues.size())
        {
            for (std::vector<Residue>::const_iterator iR=tSSBond.residues.begin(); 
                    iR !=tSSBond.residues.end(); iR++)
            {
                residues.push_back((*iR));
            }
        } 
    }
    
    SSBond::~SSBond()
    {   
        if (!residues.empty())
        {
            residues.clear();
        }
    }
}

