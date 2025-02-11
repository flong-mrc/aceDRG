/* 
 * File:   BondList.h
 * Author: flong
 *
 * Created on August 9, 2011, 4:57 PM
 */

#ifndef BONDLIST_H
#define	BONDLIST_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif


#ifndef BOND_H
#include "bond.h"
#endif



namespace LIBMOL
{
    class Bond;
    class Atom;
    class AtomAssembly;
    
    /* class Bondlist describes a list of bonds 
     * 
    */
    
    class BondList 
    {
    public :
        
        // default constructor with no bond in the list. 
        BondList();
        // the list with the common atom defined
        BondList(Atom& tAtm);
        // copy constructor
        BondList(Atom& tAtm, BondList& tBList);
        
        // destructor 
        virtual ~BondList(); 
        
        void addOneBond(Bond& tB);
        void deleteOneBond(SeriNumber tN);
        void deleteOneBond(Bond& tB);
        
        std::vector<Bond> getList();
        Bond *            getOneBond(SeriNumber tN) const;
        Size              size();
        void              destroyList();    
        
    private :
        
        Atom*               itsCommonAtom;
        std::vector<Bond>   itsList;
        
    };
}


#endif	/* BONDLIST_H */

