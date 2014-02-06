/* 
 * File:   getDerivs.h
 * Author: flong
 *
 * Created on February 18, 2013, 4:30 PM
 */

#ifndef GETDERIVS_H
#define	GETDERIVS_H

#ifndef KERNEL_H
#include "../kernel/kernel.h"
#endif

#ifndef ATOM_H
#include "../kernel/atom.h"
#endif

#ifndef BOND_H
#include "../kernel/bond.h"
#endif

#ifndef ANGLE_H
#include "../kernel/angle.h"
#endif

#ifndef TORSION_H
#include "../kernel/torsion.h"
#endif

#ifndef RESIDUE_H
#include "../kernel/residue.h"
#endif

#ifndef RING_H
#include "../kernel/ring.h"
#endif

#ifndef PLANE_H
#include "../kernel/plane.h"
#endif

#ifndef CHAIN_H
#include "../kernel/chain.h"
#endif

#ifndef UTILITY_H
#include "../kernel/utility.h"
#endif

#ifndef NEIGHBLIST_H
#include "../kernel/neighbList.h"
#endif

#ifndef ALLSYSTEM_H
#include "../kernel/AllSystem.h"
#endif

#ifndef GETOBJVALUES_H
#include "getObjValues.h"
#endif

namespace FF
{
    class GetAllDerivsFull
    {
    public:
        // Default constructor
        GetAllDerivsFull(int tNumVars);
        
        // Default destructor
        ~GetAllDerivsFull();
        
        
        
        // For first derivatives from sets of atoms, bonds and others
        void setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts, 
                                std::vector<LIBMOL::BondDict>& tBos, 
                                std::vector<LIBMOL::AngleDict>& tAns, 
                                std::vector<LIBMOL::TorsionDict>& tTos, 
                                std::vector<LIBMOL::PlaneDict>& tPls, 
                                std::vector<LIBMOL::ChiralDict>& tChs);
        
        void setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts, 
                                std::vector<LIBMOL::BondDict>& tBos, 
                                std::vector<LIBMOL::AngleDict>& tAns, 
                                std::vector<LIBMOL::TorsionDict>& tTos, 
                                std::vector<LIBMOL::PlaneDict>& tPls, 
                                std::vector<LIBMOL::ChiralDict>& tChs,
                                std::vector<LIBMOL::AtomDict>& tAllAtoms);
        
        // This version just uses atoms, pre-assumed atom forces are calculated 
        void setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts);
        
        void setFirstDerivsCartBond(std::vector<LIBMOL::AtomDict>& tAts,
                                    std::vector<LIBMOL::BondDict>::iterator tBo);
        void setFirstDerivsCartAngle(std::vector<LIBMOL::AtomDict>& tAts,
                                     std::vector<LIBMOL::AngleDict>::iterator tAn);
        void setFirstDerivsCartTorsion(std::vector<LIBMOL::AtomDict>& tAts,
                                       std::vector<LIBMOL::TorsionDict>::iterator tTo);
        void setFirstDerivsCartPlane(std::vector<LIBMOL::AtomDict>& tAts,
                                     std::vector<LIBMOL::PlaneDict>::iterator tPl,
                                     LIBMOL::REAL ** t1stDerivPlan);
        void setFirstDerivsCartChiral(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::ChiralDict>::iterator tCh);
        void SetChiraAndFirstDeriv(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<LIBMOL::ChiralDict>::iterator tCh,
                                   LIBMOL::REAL  volume,  LIBMOL::REAL  * df_dx_c);
        void setFirstDerivsCartVDW(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<LIBMOL::AtomDict>::iterator tAt);
        LIBMOL::REAL SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<int>::iterator tNB,
                                   std::vector<LIBMOL::AtomDict>::iterator tAt);
        
        // Calculate the forces on each atom in the system, these may replace the 
        // the above first derivatives because of better organization;
        void SetAllAtomForces(std::vector<LIBMOL::AtomDict>& tAts, 
                              std::vector<LIBMOL::BondDict>& tBos, 
                              std::vector<LIBMOL::AngleDict>& tAns, 
                              std::vector<LIBMOL::TorsionDict>& tTos, 
                              std::vector<LIBMOL::PlaneDict>& tPls, 
                              std::vector<LIBMOL::ChiralDict>& tChs);
        void SetBondToForce(std::vector<LIBMOL::AtomDict>& tAts,
                            std::vector<LIBMOL::BondDict>& tBos);
        void SetAngToForce(std::vector<LIBMOL::AtomDict>& tAts,
                           std::vector<LIBMOL::AngleDict>& tAns);
        void SetTorToForce(std::vector<LIBMOL::AtomDict>& tAts,
                           std::vector<LIBMOL::TorsionDict>& tTos);
        void SetPlaToForce(std::vector<LIBMOL::AtomDict>& tAts,
                           std::vector<LIBMOL::PlaneDict>& tPls);
        void SetChiToForce(std::vector<LIBMOL::AtomDict>& tAts,
                           std::vector<LIBMOL::ChiralDict>& tChs);
                           
        
        // The second derivatives related matrixes 
        // 1. the simplest one, given all atom forces have been calculated
        void SetNormalMatrix(std::vector<LIBMOL::AtomDict> & tAts);
        void SetNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts, 
                             std::vector<LIBMOL::BondDict>& tBos, 
                             std::vector<LIBMOL::AngleDict>& tAns, 
                             std::vector<LIBMOL::TorsionDict>& tTos, 
                             std::vector<LIBMOL::PlaneDict>& tPls, 
                             std::vector<LIBMOL::ChiralDict>& tChs);
        void OneBondToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                                std::vector<LIBMOL::BondDict>::iterator tBos);
        
        void AngToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                               std::vector<LIBMOL::AngleDict>& tAns);
        
        void TorToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                               std::vector<LIBMOL::TorsionDict>& tTos);
        void PlaToNormalMatrix(std::vector<LIBMOL::AtomDict> & tAts,
                               std::vector<LIBMOL::PlaneDict> & tPl, bool l_fd);
        void ChiToNormalMatrix(std::vector<LIBMOL::AtomDict> & tAts,
                               std::vector<LIBMOL::ChiralDict> & tCh, bool);
        void OneAtomVDWToNormalMatrix(std::vector<LIBMOL::AtomDict> & tAts,
                                      std::vector<LIBMOL::AtomDict>::iterator tAt);
        void OnePairVDWToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::AtomDict>::iterator tAt);
        void StableNormalMatrix();
        
        int                 numVars;
        LIBMOL::REAL      * firDrivCart;
        LIBMOL::REAL     ** secDrivCart;
        
        // 1 for full matrix format (default), 2 for sparse matrix format  
        int                 workMode;
        // 1 for Cartesian space (default),  2 for torsion space
        int                 workSpace;
        
    };
    
    extern void SetChiraAndFirstDeriv(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::ChiralDict>::iterator tCh,
                                      LIBMOL::REAL  & volume, LIBMOL::REAL  * df_dx_c);
}



#endif	/* GETDERIVS_H */

