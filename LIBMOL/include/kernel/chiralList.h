/* 
 * File:   chiralList.h
 * Author: flong
 *
 * Created on August 10, 2011, 2:45 PM
 */

#ifndef CHIRALLIST_H
#define	CHIRALLIST_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    class Chiral;
    class Atom;
    
    /* Class ChiralList represents a list (std::vector) of
     * Chiral objects
    */
    
    class ChiralList   : public AtomAssembly 
    {
    public :
        // default constructor (empty list)
        ChiralList();
        // constructor that setup the common atom only
        ChiralList(Atom & tAtm);
        // copy constructor 
        ChiralList(Chiral & tC);
        
        // destructor
        virtual ~ChiralList();

        void addOneChiral(Chiral&  tC);
        void deleteOneChiral(Chiral&  tC);
        void deleteOneChiral(SeriNumber tN);
        
        std::vector<Chiral>   getList();
        Chiral *              getOneChiral(SeriNumber tN);
        Size                  size();
        void                  destroyList();
        
    private :
 
        Atom *                itsCommonAtom;
        std::vector<Torsion>  itsList;
        
    };

}

    


#endif	/* CHIRALLIST_H */

