/* 
 * File:   TransCoord.h
 * Author: flong
 *
 * Created on February 26, 2013, 4:23 PM
 */

#ifndef TRANSCOORD_H
#define	TRANSCOORD_H

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

#ifndef RING_H
#include "ring.h"
#endif

#ifndef PLANE_H
#include "plane.h"
#endif

#ifndef ALLSYSTEM_H
#include "AllSystem.h"
#endif 

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef LINALG_H
#include "../go/LinAlg.h"
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
    class ChiralDict;
    
    class Plane;
    class PlaneDict;
    
    class Ring;
    class RingDict;
    
    class AllSystem;
    
    class TransCoords
    {
    public:
        
        // Default constructor
        TransCoords();
        
        // Default destructor
        ~TransCoords();
        
        void generateCoordTorsToCart(std::vector<AtomDict> & tAllAtoms,
                                     std::vector<int>      & tRootSet);
        
        void generateCoordTorsToCart(std::vector<AtomDict>    & tAtoms,
                                     std::vector<BondDict>    & tBonds,
                                     std::vector<AngleDict>   & tAngles,
                                     std::vector<TorsionDict> & tTorsions,
                                     std::map<ID, std::vector<RingDict> > & tRings,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChs);
        
        
        void branchGrowthTorsTorCart(std::vector<AtomDict> & tAllAtoms,
                                     int tAtom1, int tAtom2, int tAtom3);
        
        // another version
        void generateCoordTorsToCart2(std::vector<AtomDict>& tAtoms,
                                     std::vector<BondDict>    & tBonds,
                                     std::vector<AngleDict>   & tAngles,
                                     std::vector<TorsionDict> & tTorsions,
                                     std::vector<RingDict>    & tRings,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChs);
        
        void branchGrowthTorsToCart2(std::vector<AtomDict> & tAllAtoms,
                                     int tAtom1, int tAtom2, int tAtom3,
                                     std::vector<BondDict>    & tBonds,
                                     std::vector<AngleDict>   & tAngles,
                                     std::vector<TorsionDict> & tTorsions,
                                     std::vector<RingDict>    & tRings,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChs);
        
        void CHANGE_AMAT(std::vector<AtomDict> & tAllAtoms,
                         int i_next, REAL ** aMat);
        
        void PolToCart(std::vector<AtomDict> & tAllAtoms,
                       int i_cur, int i_next, 
                       REAL ** aMat);
        
        void PolToCart(std::vector<AtomDict> & tAllAtoms,
                       int i_cur, int i_next, 
                       REAL vBond, REAL vAng, REAL vTor,
                       REAL ** aMat);
        
        void PolToCart(AtomDict & tCurAtom,
                       AtomDict & tNextAtom, 
                       REAL vBond, REAL vAng, REAL vTor,
                       REAL ** aMat);
        
        void PolToCart(AtomDict & tCurAtom,
                       std::vector<AtomDict>::iterator tNextAtom, 
                       REAL vBond, REAL vAng, REAL vTor,
                       REAL ** aMat);
        
        void generateCoordCartToTors(std::vector<AtomDict> & tAllAtoms);
        
        // tempo here, put into general place later on
        void setCurStartSet(std::vector<int>         & tSet,
                            std::vector<int>         & curSet,
                            int                        tNext,
                            std::string                tDirect,
                            RingDict                &  tRing,
                            std::vector<int>           tDoneSet);
        
        void growOneAtom(std::vector<int>         & sSet,
                         std::vector<AtomDict>    & tAtoms,
                         std::vector<BondDict>    & tBonds,
                         std::vector<AngleDict>   & tAngles,
                         std::vector<TorsionDict> & tTorsions,
                         int                        tIdx);
        
        void growOneAtom(std::vector<int>         & sSet,
                         std::vector<AtomDict>    & tAtoms,
                         REAL                     tbondV,
                         REAL                     tAngV,
                         REAL                     tTorV,
                         int                      tIdx);
        
        void growOneAtom(AtomDict &                        tAtom1,
                         AtomDict &                        tAtom2,
                         AtomDict &                        tAtom3,
                         std::vector<AtomDict>::iterator   tNextAtom,
                         REAL                              tbondV,
                         REAL                              tAngV,
                         REAL                              tTorV);
        
        void growOneNode(std::vector<AtomDict> & tAllAtoms,
                         int iCur, int i_next, 
                         REAL ** aMat);
        
        void growOneRingNode(std::vector<int>         & sSet,
                             std::vector<AtomDict>    & tAtoms,
                             std::vector<BondDict>    & tBonds,
                             std::vector<AngleDict>   & tAngles,
                             std::vector<TorsionDict> & tTorsions,
                             std::vector<ChiralDict>  & tChs,
                             bool                       inR,
                             int                        iSeq,
                             int                        tIdx,
                             std::vector<int>         & tDoneSet);
        
        void growOneRingNode2(std::vector<int>         & sSet,
                             std::vector<AtomDict>    & tAtoms,
                             std::vector<BondDict>    & tBonds,
                             std::vector<AngleDict>   & tAngles,
                             std::vector<TorsionDict> & tTorsions,
                             std::vector<ChiralDict>  & tChs,
                             bool                       inR,
                             int                        iSeq,
                             int                        tIdx,
                             std::vector<int>         & tDoneSet);
        
        void ringBuilder(std::vector<int>         & sSet,
                         RingDict                 & tRing,
                         std::vector<AtomDict>    & tAtoms, 
                         std::vector<BondDict>    & tBonds,
                         std::vector<AngleDict>   & tAngles,
                         std::vector<TorsionDict> & tTorsions,
                         std::vector<RingDict>    & tRings,
                         std::vector<PlaneDict>   & tPlas,
                         std::vector<ChiralDict>  & tChs,
                         int                        tTurn,
                         std::vector<int>         & doneSet);
        void getRingAtomConn(RingDict & tRing,
                             std::map<int, int> tAtmsLink,
                             int                sIdx);
        
        void setSubSystem(std::vector<int>       & tDoneSet,
                          std::vector<AtomDict>  & tAllAtom,
                          std::vector<BondDict>  & tAllBonds,
                          std::vector<AngleDict> & tAllAngles,
                          std::vector<TorsionDict> & tAllTorsions,
                          std::vector<RingDict>    & tAllRings,
                          std::vector<PlaneDict>   & tAllPlanes,
                          std::vector<ChiralDict>  & tAllChirals,
                          AllSystem                & tSubSys);
        
        void optSubSystem(std::vector<int>       & tSubSet,
                          std::vector<AtomDict>  & tAllAtom,
                          std::vector<BondDict>  & tAllBonds,
                          std::vector<AngleDict> & tAllAngles,
                          std::vector<TorsionDict> & tAllTorsions,
                          std::vector<RingDict>    & tAllRings,
                          std::vector<PlaneDict>   & tAllPlanes,
                          std::vector<ChiralDict>  & tAllChirals,
                          bool                       useGO);

        
    private:
        
        std::vector<int>            doneList;
       
    };

}
#endif	/* TRANSCOORD_H */

