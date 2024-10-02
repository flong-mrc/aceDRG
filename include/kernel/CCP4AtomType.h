/*
 * File:   CCP4AtomType.h
 * Author: flong
 *
 * Created on September 6, 2012, 11:32 AM
 */

#ifndef CCP4ATOMTYPE_H
#define	CCP4ATOMTYPE_H

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

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif

namespace LIBMOL
{

    class AtomDict;
    class BondDict;
    class AngleDict;
    class TorsionDict;
    class ChiralDict;
    class Chiral;
    class Ring;
    class RingDict;

    class PeriodicTable;

    class CCP4AtomType
    {
    public:

        // Default constructor
        CCP4AtomType();

        // Copy constructor
        CCP4AtomType(const CCP4AtomType & tC);

        // Another constructor
        CCP4AtomType(const std::vector<AtomDict> & tAllAtoms,
                     const std::map<ID, std::vector<RingDict> > & tAllRings);

        // Another constructor
        CCP4AtomType(const std::vector<AtomDict> & tAllAtoms,
                     const std::vector<RingDict> & tAllRingsV);

        // Default destructor
        ~CCP4AtomType();

        // The following should be done with COD atom types as well
        void setHydroAtomCCP4Type(AtomDict & tAtom);       // H
        void setOrgAtomCCP4Type(AtomDict & tAtom);  // C, N, O, P, S, Se
        void SetAlkaliMetalsAtomCCP4Type(AtomDict & tAtom);
        void SetAlkalineEarthMetalsAtomCCP4Type(AtomDict & tAtom);
        void SetTransitionMetalsAtomCCP4Type(AtomDict & tAtom);
        void SetOtherMetalAtomCCP4Type(AtomDict & tAtom);
        void SetSemimetallicsAtomCCP4Type(AtomDict & tAtom);
        void SetHalogensAtomCCP4Type(AtomDict & tAtom);
        void SetRareEarthAtomCCP4Type(AtomDict & tAtom);
        void SetInertGasesAtomCCP4Type(AtomDict & tAtom);
        void setOneAtomCCP4Type(PeriodicTable & tP, AtomDict & tAtom);
        void setAllAtomsCCP4Type();

        std::vector<AtomDict>                    allAtoms;
        std::map<ID, std::vector<RingDict> >     allRings;
        std::vector<RingDict>                    allRingsV;
    };

    class CCP4DictParas
    {
    public :

        // Default constructor
        CCP4DictParas();

        // Destructor
        ~CCP4DictParas();


        void getAtomPropsTable(std::vector<std::vector<std::string> >::iterator  tData);
        void getBondPropsTable(std::vector<std::vector<std::string> >::iterator  tData);
        void getAnglePropsTable(std::vector<std::vector<std::string> >::iterator  tData);

        std::map<std::string, std::string>                    atomTypeElementTable;
        std::map<std::string, std::map<std::string, REAL> >   atomPropsTable;
        std::map<std::string, std::map<std::string,
                 std::map<std::string, REAL> > >              bondPropsTable;
        std::map<std::string, std::map<std::string, std::map<std::string,
        std::map<std::string, REAL> > > >                     anglePropsTable;


    };


    class AtomTypeTool
    {
    public:

        //Default constructor
        AtomTypeTool();

        // The constructor using a file name, file type indicator, which an integer
        // File types are e.g. CIF, MOL/SDF, PDB, SMILE, inChi etc
        // Atom types are COD, CCP4 etc
        AtomTypeTool(FileName  tFname, FileType tFType);

        // The constructor using atoms, bonds, and rings
        AtomTypeTool(std::vector<AtomDict>                  & tAtoms,
                     std::vector<BondDict>                  & tBonds,
                     std::map<ID, std::vector<RingDict> >   & tRings);

        // Default destructor
        ~AtomTypeTool();

        void execute();


        std::vector<AtomDict>                    allAtoms;
        std::vector<BondDict>                    allBonds;
        std::map<ID, std::vector<RingDict> >     allRings;

    };

    extern void setAtomsCCP4Type(std::vector<AtomDict> & tAtoms,
                                std::vector<RingDict>  & tRings);

}



#endif	/* CCP4ATOMTYPE_H */
