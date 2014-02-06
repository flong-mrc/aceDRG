/* 
 * File:   link.h
 * Author: flong
 *
 * Created on September 6, 2011, 6:27 PM
 */

#ifndef LIBG_LINK_H
#define	LIBG_LINK_H

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
     
    class Link
    {
    public:
        
        // default constructor
        Link();
        
        // copy constructor
        Link(const Link & tLink);
        
        // destructor
        ~Link();      
        
        inline REAL getLength()
        {
            return itsLength;
        }
        inline void setLength(REAL tL)
        {
            itsLength = tL;
        }
        
        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        SeriNumber getModSeriNum();
        void  setModSeriNum(SeriNumber tN);
        
        // calculate the length using two atoms in the link
        void setLength();
        
        
        std::vector<Atom>    atoms;
        int                  symOp[2];
        
        
    private:
        
        /* defined in class AtomAssembly*/
        Name              itsName;
        ID                itsID;
        SeriNumber        itsSeriNum;   
        
        REAL              itsLength;
        

        
        
    };
}


#endif	/* LINK_H */

