/* 
 * File:   atomsTree.h
 * Author: flong
 *
 * Created on February 26, 2013, 8:22 PM
 */

#ifndef ATOMSTREE_H
#define	ATOMSTREE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
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

#ifndef PLANE_H
#include "plane.h"
#endif


#ifndef UTILITY_H
#include "utility.h"
#endif

namespace LIBMOL
{
    
    class buildAtomTree
    {
    public:
        
        // Default constructor
        buildAtomTree();
        // Default destructor 
        ~buildAtomTree();
        
        int  setMSTStartAtom(std::vector<AtomDict> & allAtoms);
        bool setAtomsMST(std::vector<AtomDict> & allAtoms, 
                         std::vector<BondDict> & allBonds);
        
        void initDummyAtoms(std::vector<AtomDict> & allAtoms,
                            std::vector<AtomDict> & tDAtoms);
        
        void setOneDummyAtom(std::vector<AtomDict> & t3Atoms,
                             AtomDict & tDAtom);
        
        void setStartAtom(std::vector<AtomDict> & tAllAtoms);
        
        void setStartStruct(std::vector<AtomDict>  & tAllAtoms,
                            std::vector<BondDict>  & tAllBonds,
                            std::vector<AngleDict> & tAllAngs,
                            std::vector<ChiralDict>& tChs);
        void setStartStruct2(std::vector<AtomDict>   & tAllAtoms,
                             std::vector<BondDict>   & tAllBonds,
                             std::vector<AngleDict>  & tAllAngs,
                             std::vector<ChiralDict> & tChs);
        
        void expandStartStruct(std::vector<AtomDict>  & tStartPackAtoms);
        
        void set2ConnStruct(std::vector<AtomDict>  & tAllAtoms,
                            std::vector<BondDict>  & tAllBonds,
                            std::vector<AngleDict> & tAllAngs);
        void set3ConnStruct(std::vector<AtomDict>  & tAllAtoms,
                            std::vector<BondDict>  & tAllBonds,
                            std::vector<AngleDict> & tAllAngs,
                            std::vector<ChiralDict>& tChs);
        void set4ConnStruct(std::vector<AtomDict>  & tAllAtoms,
                            std::vector<BondDict>  & tAllBonds,
                            std::vector<AngleDict> & tAllAngs,
                            std::vector<ChiralDict>& tChs);
        
        void setBranches(std::vector<AtomDict> & tAtoms);
  
        
        // set values of bonds, angles, torsions for all atoms in tree
        // one value of bond, angle, torsion for each atom
        // version 1:  from ideal values of bonds, angles, and torsions
        void setTreeAtomValues(std::vector<AtomDict>    & tAtoms,
                               std::vector<BondDict>    & tBonds,
                               std::vector<AngleDict>   & tAngles,
                               std::vector<TorsionDict> & tTorsions);
        void setTreeAtomValues2(std::vector<AtomDict>    & tAtoms,
                                std::vector<BondDict>    & tBonds,
                                std::vector<AngleDict>   & tAngles,
                                std::vector<TorsionDict> & tTorsions);        
        
        void setTreeAtomBondValue(std::vector<AtomDict> & tAtoms,
                                  std::vector<BondDict> & tBonds,
                                  std::vector<AtomDict>::iterator tAt);
        void setTreeAtomAngleValue(std::vector<AtomDict>  & tAtoms,
                                   std::vector<AngleDict> & tAngs,
                                   std::vector<AtomDict>::iterator tAt);
        void setTreeAtomAngleValue(std::vector<AtomDict>& tAtoms, 
                                   std::vector<AngleDict>& tAngs,
                                   std::vector<AtomDict>::iterator tAt,
                                   int tP, int tGp);
        void setTreeAtomTorsionValue(std::vector<AtomDict>    & tAtoms,
                                     std::vector<TorsionDict> & tTors,
                                     std::vector<AtomDict>::iterator tAt, 
                                     int tR);
        void setTreeAtomTorsionValue(std::vector<AtomDict>    & tAtoms,
                                     std::vector<TorsionDict> & tTors,
                                     AtomDict                 & tDAt1,
                                     std::vector<AtomDict>::iterator tAt);
        void setTreeAtomTorsionValue(std::vector<AtomDict>    & tAtoms,
                                     std::vector<TorsionDict> & tTors,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChirals,
                                     std::vector<AtomDict>    & tDAts,
                                     std::vector<AtomDict>::iterator tAt);
        void setTreeAtomTorsionValue(std::vector<AtomDict>    & tAtoms,
                                     std::vector<TorsionDict> & tTorsions, 
                                     std::vector<AtomDict>::iterator iAt, 
                                     std::vector<int>         & tVec);
        
        // version 2 from the coordinates of the atoms
        void setTreeAtomValues(std::vector<AtomDict>     & tAtoms);
        
        
        void setTreeAtomDepth();
        
        bool buildTree(std::vector<AtomDict>    & tAtoms,
                       std::vector<BondDict>    & tBonds,
                       std::vector<AngleDict>   & tAngles,
                       std::vector<TorsionDict> & tTorsions,
                       std::vector<RingDict>    & tRings,
                       std::vector<PlaneDict>   & tPlas,
                       std::vector<ChiralDict>  & tChs);
        
        
        int                                                startAtom;
        std::map<int, std::map<int, std::vector<int> > >   branches;
        std::vector<int>                                   startPack;
    };
}


#endif	/* ATOMSTREE_H */

