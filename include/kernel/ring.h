/* 
 * File:   ring.h
 * Author: flong
 *
 * Created on February 3, 2012, 10:26 AM
 */

#ifndef RING_H
#define	RING_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIB_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef PLANE_H
#include "plane.h"
#endif

//#ifndef CHEMPROPSET_H
//#include "chemPropSet.h"
//#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    
    class Bond;
    class BondDict;
    class Angle;
    class AngleDict;
    class Torsion;
    class TorsionDict;
    class PlaneDict;
    

    class Ring 
    {
    public:
        
        // Default constructor
        Ring();
        
        // Copy constructor
        Ring(const Ring & tRing);
        
        // Destructor
        ~Ring();
        
        
        std::list<AtomDict>    atoms;
        std::list<BondDict>    bonds;
        std::vector<Angle>     angles;
        
        std::vector<Torsion> torsions; 
        
        // some properties a ring has
        
    };
    
    class RingDict 
    {
    public:
        // Default constructor 
        RingDict();
        
        // Copy constructor 
        RingDict(const RingDict & tR);
        
        // Constructor using a set of atoms
        RingDict(const std::vector<AtomDict> tAtoms);
        
        // Destructor 
        ~RingDict(); 
        
        void setAtmsLink(std::vector<AtomDict> tAllAtoms);
        void setPlaneProp();
        
        bool                                     isPlanar;
        std::string                              rep;
        std::string                              sRep;
        
        std::vector<AtomDict>                    atoms;
        std::map<int, std::map<ID, int> >        atomsLink; 
              
    };
    
    extern bool AtomsInSameRing(AtomDict & tAt1, AtomDict tAt2, 
                std::vector<RingDict> & tRings);
    
    extern void buildOneRing(std::vector<AtomDict> & tAtoms, 
                             std::vector<AtomDict>   sSet,
                             std::vector<AtomDict> & doneSet);
    
    extern void ringAtmsLink(RingDict & tRing, 
                             int        tIdx,
                             std::vector<int> & tAtmsLink);
    
    extern void mergePlaneRings(std::vector<RingDict> & tAllRings,
                                std::vector<std::vector<int> > & tRingAtms);
    
    extern void mergePlaneRings(std::vector<RingDict>          & tAllRings,
                                std::vector<std::vector<int> > & tRingAtms,
                                std::vector<AtomDict>          & tAtoms,
                                std::vector<BondDict>          & tBonds);
    
    extern void findAllRingsConnectedOneRing(int tCurIdx, std::vector<RingDict> & tRings,
                                             std::vector<int>   & DoneList, 
                                             std::vector<int>   & curLinkedRing);
    
    extern REAL getPiForOneAtom(int tIdx, std::vector<AtomDict> & tAtoms, 
                                std::vector<BondDict>  & tBonds);
    
    extern bool checkAromaSys(std::vector<int>      & tSubAtoms,
                              std::vector<AtomDict> & tAtoms,
                              std::vector<BondDict> & tBonds);
    
    extern void checkAndSetupPlanes(std::vector<RingDict>  & tAllRings,
                                    std::vector<PlaneDict> & tPlanes,
                                    std::vector<AtomDict>  & tAtoms,
                                    std::vector<BondDict>  & tBonds);
    
    
    
    
}     // end of namespace LIBMOL

#endif	/* RING_H */

