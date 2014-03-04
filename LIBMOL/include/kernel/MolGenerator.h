/* 
 * File:   MolGenerator.h
 * Author: flong
 *
 * Generator a molecule from input sets of atoms and symmetry operators 
 * 
 * Created on August 13, 2013, 5:26 PM
 */

#ifndef MOLGENERATOR_H
#define	MOLGENERATOR_H


#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef CIFFILE_H
#include "DictCifFile.h"
#endif

#ifndef CCP4ATOMTYPE_H
#include "CCP4AtomType.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef MOLECULE_H
#include "molecule.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef NEIGHBLIST_H
#include "neighbList.h"
#endif

#ifndef PDBFILE_H
#include "PDBFile.h"
#endif

namespace LIBMOL
{
    class SpaceGroupMember;
    class Cell;
    
    class AtomDict;
    class Molecule;
    
    class CrystInfo;
    
    class GenCifFile;
    class CCP4DictParas;
    
    class PeriodicTable;
    
    class MolGenerator
    {
    public:
        //Default constructor 
        MolGenerator();
        
        // Constructor using a cif object
        MolGenerator(const GenCifFile &  tCifObj,
                     int        tNBDepth=1);
        
        MolGenerator(const DictCifFile &  tCifObj,
                     int        tNBDepth=1);                                 
        // Destructor
        ~MolGenerator();
        
        // Generate a "molecule"(not normal sense) using a set of atoms and 
        // a set of symmetry operators
        void buildRefAtoms(std::vector<CrystInfo>::iterator  iCryst);
        
        void symmAtomGen(std::vector<CrystInfo>::iterator    tCrys,
                         PeriodicTable & tPTab);
        
        void getOneSymmAtom(std::vector<AtomDict>::iterator        tCurAtom,
                            std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
                            std::vector<CrystInfo>::iterator   tCrys, 
                            PeriodicTable & tPTab);
        
        bool checkUniqueAtom(std::vector<REAL>      & tFracX);
        
        void addOneSetRefAtoms(AtomDict                         & tCurAtom,
                               std::vector<CrystInfo>::iterator   tCryst);
        
        void setUniqueAtomLinks(PeriodicTable & tPTab);
        
        void getMolByEqClassInCell();
        
        void getMolByEqClassInCrys();
        
        void getUniqueBonds(PeriodicTable & tPTab);
        
        void getBondingRangePairAtoms(AtomDict & tAtm1,
                                      AtomDict & tAtm2,
                                      REAL tExtraD,
                                      PeriodicTable & tPTab,
                                      std::vector<REAL> & tRange);
        
        void setOneUniqueBond(int     tIdxAtm1,
                              int     tIdxAtm2, 
                              REAL    rD);
        
        void checkBondOneForming(PeriodicTable & tPTab);
        
        void getUniqAngles();
        
        void getMolsInCell();
        bool checkAtomOcp(std::vector<int> & tMol);
        
        void outTables(FileName tOutName, Molecule & tMol);
        
        void execute(FileName tOutName);   
        
        std::vector<CCP4DictParas>      ccp4DictParas;
        std::vector<AtomDict>           initAtoms;
        std::vector<AtomDict>           allAtoms;       // these are all atoms (including symmetrically generated atoms)
                                                        // in a unit cell
        std::vector<AtomDict>           refAtoms;       
        std::vector<AtomDict>           uniqueAtoms; 
        std::vector<Molecule>           allMolecules;
        std::vector<CrystInfo>          allCryst;
        
        std::vector<BondDict>           bonds;
        std::vector<AngleDict>          angles;
        
        std::map<unsigned, std::vector<int> >  moleculesInCell;
        std::map<unsigned, std::vector<int> >  moleculesInCryst;
        
    private :
        
        int                                                 myNBDepth;
        std::map<std::string, std::map<std::string, REAL> > shortLenList;
                
    };
    
}




#endif	/* MOLGENERATOR_H */

