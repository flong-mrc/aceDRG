/* 
 * File:   tosionList.h
 * Author: flong
 *
 * Created on August 10, 2011, 1:28 PM
 */

#ifndef TOSIONLIST_H
#define	TOSIONLIST_H


#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    class Torsion;
    class Atom;
    
    /* Class TorsionList represents a list(std::vector) of Torsion objects
    */
    
    class TorsionList  : public AtomAssembly
    {
        // default constructor (empty list)
        TorsionList();
        // the list with the common atom defined
        TorsionList(Torsion * tTor);
        // copy constructor
        TorsionList(TorsionList & tTorL);
        
        // destructor
        virtual ~TorsionList();
        
        void addOneTorsion(Torsion &tTor);
        void deleteOneTorsion(Torsion &tTor);
        void deleteOneTorsion(SeriNumber tIT);
        
        std::vector<Torsion>  getList();
        Torsion *             getOneTorsion(SeriNumber tIT);
        Size                  size();
        void                  destroyList();
        
    private :
        
        std::vector<Torsion>  itsList;
        
    };
}



#endif	/* TOSIONLIST_H */

