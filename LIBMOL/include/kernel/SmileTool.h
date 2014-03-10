/* 
 * File:   SmileTool.h
 * Author: flong
 *
 * Created on September 11, 2012, 1:55 PM
*/

#ifndef SMILETOOL_H
#define	SMILETOOL_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef TORSION_H
#include "torsion.h"
#endif

#ifndef RING_H
#include "ring.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif

#ifndef CCP4ATOMTYPE_H
#include "CCP4AtomType.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

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
    class Chiral;
    
    class Ring;
    class RingDict;
    
    class PeriodicTable;
    class CCP4AtomType;
    
    class SmileTool
    {
    public:
        
        // Default constructor
        SmileTool();
        
        // The constructor using smile string
        SmileTool(std::string & tSStr);
        
        // Another constructor using a group of atoms
        SmileTool(const std::vector<AtomDict> & tAllAtoms,
                  const std::vector<BondDict> & tAllBonds,
                  const std::map<ID, std::vector<RingDict> > & tAllRings);
        
        // Destructor
        ~SmileTool();
        
        // Define OrgSet 
        void setOrgSet();
        // Define those special characters 
        void setSpecialChars();
        bool isValidAtomName(ID tId);
        
        // SMILE string to a molecule
        void SmileToMol(ID & tSS, std::vector<AtomDict> & tMol);
        void SmileToMol(ID & tSS);
        void checkAndSetOneAtom(ID tSS, ID tId);
        void setOneAtomName(ID & tSS, ID tId);
        void setOneAtomProp(ID & tSS, ID tId);
        void setBranch(ID & tSS, ID tId);
        void setChiral(ID & tSS, ID tId);
        void setDisconnectedUnits(ID & tSS, ID tId);
        void setOneBondProp(ID & tSS, ID tId);
        
        
        void SmileToCifUsingLibcheck(FileName       tSmiFileName,
                                     std::string  & tOutLibName);
        
        // Predefined Organic set of element 
        std::map<ID, std::vector<int> >      OrganSet;
        std::map<ID, std::vector<ID> >       specialChars;
   
        std::vector<AtomDict>                allAtoms;
        std::vector<BondDict>                allBonds;
        std::vector<ChiralDict>              allChirals;
        
    private:
        
        ID                  itsCurStr;
        ID                  itsPrevStr;
        int                 itsCurPos;
        int                 itsPrevPos;
        
        int                 itsCurAtomIdx;
        AtomDict   *        itsCurAtom;
        AtomDict   *        itsPrevAtom;
        
    };
    
}


#endif	/* SMILETOOL_H */

