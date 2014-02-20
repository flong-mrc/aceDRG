/* 
 * File:   MolSdfFile.h
 * Author: flong
 *
 * Created on October 31, 2012, 9:18 PM
 */

#ifndef MOLSDFFILE_H
#define	MOLSDFFILE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIBG_ATOMASSEMBLY_H
#include "atomAssembly.h"
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

#ifndef MOLECULE_H
#include "molecule.h"
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
    class RinDict;
    
    class molecule;
    
    class CrystInfo;
    
    class Link;
    
    class PeriodicTable;
    class CCP4AtomType;
    
    class MolSdfFile : public File
    {
    
    public :
        
        // Default constructor
        MolSdfFile();
        
        // Constructor using a file name 
        MolSdfFile(Name                    tFname, 
                   std::ios_base::openmode tOpenMode);
        
        MolSdfFile(FileName                    tFname,
                    std::ios_base::openmode     tOpenMode);
        
        // destructor
        ~MolSdfFile();
        
        
        void setupSystem();
        
        void createCurMol();
        void deleteCurMol();
        
        void addHAtomToMols(int tIdxMol);
        void addHAtoms(int tIdxMol, 
                       int tIdxAtm,
                       REAL tNumH);
        
        REAL  getBondOrder(int tIdxMol,
                          int tIdx1, int tIdx2);
        int  getNumOxyConnect(int tIdxMol, std::vector<AtomDict>::iterator iA);
        void setAtomsBondingAndChiralCenter(int tIdxMol);
        void setChiral(int tIdxMol);

        ID                      molName;
        ID                      creatProg; 
        bool                    hasCoords;
        bool                    hasConnect;
        bool                    hasH;
        bool                    containMetal;
        
        std::vector<Molecule>   allMols;

        std::ofstream           outFile;
        std::ifstream           inFile;
        

        
    private:
        
        Molecule            *   itsCurMol;
        
    };    
}


#endif	/* MOLSDFFILE_H */

