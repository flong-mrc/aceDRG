/* 
 * File:   angleList.h
 * Author: flong
 *
 * Created on August 9, 2011, 8:25 PM
 */

#ifndef ANGLELIST_H
#define	ANGLELIST_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    class Bond;
    class Atom;
    class AtomAssembly;
    class Angle;
    
    /* Class AngleList represents a group of bond-angles
     */
    
    class AngleList : public AtomAssembly
    {
        
    public:
        
        // Default constructor
        AngleList();
        // constructor using a list of Angle objects
        AngleList(Angle* tAngs);
        // Copy constructor 
        AngleList(AngleList &tAngL);
        
        void addOneAngle(Angle &tAng);
        void deleteOneAngle(Angle &tAng);
        void deleteOneAngle(SeriNumber tSerAng);

        std::vector<Angle>  getList() const;
        Angle*              getOneAngle(SeriNumber tN);
        Size                size();
        void                destroyList();     
        
    private :
        
        Atom *                  itsCommonAtom;
        std::vector<Angle>      itsList;
        
    };
    
    
}
#endif	/* ANGLELIST_H */

