/* 
 * File:   ssbond.h
 * Author: flong
 *
 * Created on September 6, 2011, 5:00 PM
 */

#ifndef SSBOND_H
#define	SSBOND_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

namespace LIBMOL
{
    class Atom;
    class Residue;
    
    class SSBond : public Bond
    {
    public :
        
        // default constructor 
        SSBond();
        
        //Copy constructor 
        SSBond(const SSBond & tSSBond);
        
        // destructor
        ~SSBond();
        
        inline REAL  getLength()
        {
            return itsLength;
        }
        inline void  setLength(REAL tL)
        {
            itsLength = tL;
        }
        
        std::vector<Residue>  residues;
        SeriNumber            resSeqNums[2];
        int                   resSym[2];
        
        
   //  private :
        
        /* defined in class AtomAssembly */
        Name           itsName;
        ID             itsID;
        SeriNumber     itsSeriNum;
        
        REAL           itsLength;
        
    };
}

#endif	/* SSBOND_H */

