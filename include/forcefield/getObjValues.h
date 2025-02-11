/* 
 * File:   getObjValues.h
 * Author: flong
 *
 * Created on February 18, 2013, 12:05 AM
 */

#ifndef GETOBJVALUES_H
#define	GETOBJVALUES_H

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

#ifndef GETDERIVS_H
#include "getDerivs.h"
#endif

namespace FF
{
    class GetObjValue 
    {
    public:
        // Default Constructor
        GetObjValue();
        
        // Default destructor
        ~GetObjValue();
        
        
        LIBMOL::REAL getAll(std::vector<LIBMOL::AtomDict> & tAt,
                             std::vector<LIBMOL::BondDict> & tBo,
                             std::vector<LIBMOL::AngleDict> & tAn,
                             std::vector<LIBMOL::TorsionDict> & tTo,
                             std::vector<LIBMOL::PlaneDict> & tPl,
                             std::vector<LIBMOL::ChiralDict> & tCh);
        
        LIBMOL::REAL getAll(std::vector<LIBMOL::AtomDict> & tAt,
                             std::vector<LIBMOL::BondDict> & tBo,
                             std::vector<LIBMOL::AngleDict> & tAn,
                             std::vector<LIBMOL::TorsionDict> & tTo,
                             std::vector<LIBMOL::PlaneDict> & tPl,
                             std::vector<LIBMOL::ChiralDict> & tCh,
                             std::vector<LIBMOL::AtomDict>   & tAllAtoms);
        
        LIBMOL::REAL getAll(LIBMOL::REAL * vars,
                            std::vector<LIBMOL::AtomDict> & tAt,
                            std::vector<LIBMOL::BondDict> & tBo,
                            std::vector<LIBMOL::AngleDict> & tAn,
                            std::vector<LIBMOL::TorsionDict> & tTo,
                            std::vector<LIBMOL::PlaneDict> & tPl,
                            std::vector<LIBMOL::ChiralDict> & tCh);
        
        LIBMOL::REAL getAll(std::vector<LIBMOL::REAL>  & vars,
                            std::vector<LIBMOL::AtomDict> & tAt,
                            std::vector<LIBMOL::BondDict> & tBo,
                            std::vector<LIBMOL::AngleDict> & tAn,
                            std::vector<LIBMOL::TorsionDict> & tTo,
                            std::vector<LIBMOL::PlaneDict> & tPl,
                            std::vector<LIBMOL::ChiralDict> & tCh);

        LIBMOL::REAL GetObjOneBond(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<LIBMOL::BondDict>::iterator  tBo);
        LIBMOL::REAL GetObjOneAngle(std::vector<LIBMOL::AtomDict>& tAts,
                                    std::vector<LIBMOL::AngleDict>::iterator tAn);
        LIBMOL::REAL GetObjOneTorsion(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::TorsionDict>::iterator tTo);
        LIBMOL::REAL GetObjOnePlane(std::vector<LIBMOL::AtomDict>& tAts,
                                    std::vector<LIBMOL::PlaneDict>::iterator tPl);
        LIBMOL::REAL GetObjOneChiral(std::vector<LIBMOL::AtomDict>& tAts,
                                     std::vector<LIBMOL::ChiralDict>::iterator tCh);
        
        LIBMOL::REAL GetObjOneAtomVdw(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::AtomDict>::iterator tAt);
        LIBMOL::REAL CalcObjHYD(std::vector<LIBMOL::AtomDict>& tAts,
                                std::vector<int>::iterator tNB,
                                std::vector<LIBMOL::AtomDict>::iterator tAt);
        bool         CheckHYD(std::vector<LIBMOL::AtomDict>& tAts,
                              std::vector<int>::iterator tNB,
                              std::vector<LIBMOL::AtomDict>::iterator tAt);
        LIBMOL::REAL CalcObjVDWOnePair(std::vector<LIBMOL::AtomDict>& tAts,
                                       std::vector<int>::iterator tNB,
                                       std::vector<LIBMOL::AtomDict>::iterator tAt); 
        LIBMOL::REAL SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<int>::iterator tNB,
                                   std::vector<LIBMOL::AtomDict>::iterator tAt);
        
                
        int          workSpace;
        int          lComp;
        
        
    };
    
    LIBMOL::REAL extern SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<int>::iterator tNB, 
                                      std::vector<LIBMOL::AtomDict>::iterator tAt);
}


#endif	/* GETOBJVALUES_H */

