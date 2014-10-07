/* 
 * File:   chemPropSet.h
 * Author: flong
 *
 * Created on July 3, 2014, 1:27 PM
 */

#ifndef CHEMPROPSET_H
#define	CHEMPROPSET_H

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

#ifndef PLANE_H
#include "plane.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef MOLECULE_H
#include "molecule.h"
#endif

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

#ifndef CODCLASSIFY_H
#include "codClassify.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef TRANSCOORD_H
#include "TransCoord.h"
#endif

#ifndef MOLSDFFILE_H
#include "MolSdfFile.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif

namespace LIBMOL
{
    
    class Atom;
    class AtomDict;
    
    class Bond;
    class BondDict;
    
    class Ring;
    class RingDict;
    
    class Angle;
    class AngleDict;
    
    class Torsion;
    class TorsionDict;
    
    class Chiral;
    class ChiralDict;
    
    class Plane;
    class PlaneDict;
    
    class Molecule;
    
    class CodClassify;
    
    class CCP4AtomType;
    
    class MolSdfFile;
    
    class PeriodicTable;
    
    
    extern bool assignElementType(PeriodicTable & tP, std::string tStr,  
                                  std::vector<AtomDict>::iterator tAtom);
    
    extern int getNumOxyConnect(std::vector<AtomDict>  &  tAtoms,
                                std::vector<AtomDict>::iterator iA);
    
    extern void getHydroAtomConnect(std::vector<AtomDict>  &  tAtoms);
    
    extern void setAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms);
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder);
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol);
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh);
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder,
                                std::vector<AtomDict> &tAllAtoms);
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol);
    
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  Molecule   & tMol, 
                                  PeriodicTable & tTab);
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh);
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder);
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol);
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh);
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                REAL tTolBondOrder);
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol);
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh);
    
    extern void checkProtonatedCarBoxylicAcids(std::vector<AtomDict>::iterator tIA,
                                               std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               REAL tPH=7.0);
    
    extern void checkProtonatedSulfuricAcids(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH=7.0);
    extern void checkProtonatedNAcids(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH=7.0);
    extern void checkProtonatedPAcids(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH=7.0);
    
    // Bond order sections 
        
    extern REAL  getTotalBondOrder(Molecule   & tMol,
                                   std::vector<AtomDict>::iterator  tIA);
    extern REAL  getBondOrder(Molecule   & tMol,
                              int tIdx1, int tIdx2);
    extern REAL  getBondOrderOneAtom(std::vector<BondDict> tBonds,
                                     std::vector<AtomDict> tAtoms,
                                     int tIdx1, int tIdx2);
    extern void setAllBondOrders(std::vector<AtomDict> & tAtoms,
                                 std::vector<BondDict> & tBonds);
    
    extern void  setOneHAtomCoords();
    
    extern void  setOneHAtomCoordsSP3(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP2(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP(std::vector<AtomDict> & tAtoms,
                                     std::vector<AtomDict>::iterator tIA);
    
    // atom and bond stereo  
    extern void checkAllStereo(FileName tMdlIn, FileName tPdbIn,
                               FileName tPdbOut);
    extern void checkStereoOneMol(std::vector<Molecule>::iterator tMol,
                                  FileName tPdbOut);
    
    
    // Check aromatic stability of extended ring systems.
    // using it to decide planes in the systems(molecules or monomers)
    
    
    
    extern void setAromPlanes(std::vector<AtomDict> & tAtoms,
                              std::vector<RingDict> & tRings,
                              std::vector<PlaneDict> & tPlans);
    
    
}

#endif	/* CHEMPROPSET_H */

