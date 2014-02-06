/* 
 * File:   planeList.h
 * Author: flong
 *
 * Created on August 10, 2011, 4:21 PM
 */

#ifndef PLANELIST_H
#define	PLANELIST_H


#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    class Atom;
    class Plane;
    
    /* Class PlaneList represents a group of Plane objects
     */
    
    class PlaneList 
    {
        // default constructor (empty list)
        PlaneList();
        // the list with the common atom defined
        PlaneList(Atom &tAtm);
        // copy constructor
        PlaneList(PlaneList & tPL);
        
        // destructor
        virtual ~PlaneList();
        
        void addOnePlane(Plane &tPl);
        void deleteOnePlane(Plane &tPl);
        void deleteOnePlane(SeriNumber tN);
        
        std::vector<Plane>    getList();
        Plane *               getOnePlane(SeriNumber tPl);
        Size                  size();
        void                  destroyList();
        
    private :
        
        Atom *                itsCommonAtom;
        std::vector<Plane>    itsList;
        
    };
    
}

#endif	/* PLANELIST_H */

