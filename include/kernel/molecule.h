/* 
 * File:   molecule.h
 * Author: flong
 *
 * Created on November 1, 2012, 1:17 PM
 */

#ifndef MOLECULE_H
#define	MOLECULE_H

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

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef RING_H
#include "ring.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
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

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    class Residue;
    class Chain;
    
    class Bond;
    class BondDict;
    class Angle;
    class AngleDict;
    class Torsion;
    class TorsionDict;
    class Chiral;
    class ChiralDict;
    
    class Ring;
    class RingDict;
    
    class PeriodicTable;
    class CCP4AtomType;
    
    class Molecule 
    {
        // Not a molecule in chemistry sense, but rather a set of atoms 
        // related by special connections
        
    public :
        
        // Default constructor
        Molecule();
        // Copy constructor
        Molecule(const Molecule & aMol);
        // Destructor
        ~Molecule();
        
        void setFormula();
        void calcSumExcessElecs();
        void calcSumCharges();
             
        std::vector<AtomDict>            atoms;
        std::vector<BondDict>            bonds;
        std::vector<BondDict>            allBonds;
        std::vector<AngleDict>           angles;
        std::vector<TorsionDict>         torsions;
        std::vector<RingDict>            rings;
        std::vector<ChiralDict>          chirals;
        std::vector<AtomDict>            extraHAtoms;
        
        std::map<ID, std::vector<ID> >   propData;
        std::vector<std::string>         comments;
        
        int                              seriNum;
        ID                               id;
        ID                               formula;
        int                              sumExcessElecs;
        int                              sumCharges;
        REAL                             atomCovRadMax; 
        bool                             hasCoords;
        bool                             validated;   
        bool                             isInf;
        bool                             stateChanged;
    };
}

#endif	/* MOLECULE_H */

