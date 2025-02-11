/* 
 * File:   restraintLists.h
 * Author: flong
 *
 * last updated on September 22, 2011, 3:26 PM
 */

#ifndef RESTRAINTLISTS_H
#define	RESTRAINTLISTS_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef TOSION_H
#include "torsion.h"
#endif

#ifndef CHIRAL_H
#include "chiral.h"
#endif

#ifndef PLANE_H
#include "plane.h"
#endif

namespace LIBMOL
{
    class Bond;
    class Angle;
    class Torsion;
    class Chiral;
    class Plane;
    
    /* Class RestraintList represents a set of Lists used as restraints
     * in optimization processes
     */

        
    class RestraintLists  
    {
    public :
        
        // default constructor 
        RestraintLists();
              
        // copy constructor
        RestraintLists(RestraintLists  & tRLists);
        
        // destructor 
        ~RestraintLists();
        
        //RestraintLists & operator=(const RestraintLists & tR); 

        std::vector<Bond>             restrBondList; 
        std::vector<Angle>            restrAngleList;
        std::vector<Torsion>          restrTorsionList;
        std::vector<Chiral>           restrChiralList;
        std::vector<Plane>            restrPlaneList;
                
        
    };
    
}

#endif	/* RESTRAINTLISTS_H */

