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

#ifndef CHIRAL_H
#include "chiral.h"
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

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef TRANSCOORD_H
#include "TransCoord.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    class Residue;
    class ResidueDict;
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
    
    class Molecule;
    
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
        
        MolSdfFile(FileName                    tInSDFName);
        
        // destructor
        ~MolSdfFile();
        
        
        void setupSystem();
        void setupSystemSimp();
        
        void createCurMol();
        void deleteCurMol();
        
        void addHAtomToMols(int tIdxMol);
        void setHAtomCoordsMols(int tIdxMol);
        
        REAL checkProtonated(std::vector<AtomDict>::iterator tIA, 
                             int   tMolIdx);           // This is a temp function. the other
                                                       // overloaded function should be a 
                                                       // proper one
        REAL checkProtonated(std::vector<AtomDict>::iterator tIA, 
                             int   tMolIdx,     
                             REAL  tTotalVal,
                             REAL  tPka, REAL tPh);
        
        void setAtomsCCP4Type(int idxMol);
        
        void addHAtoms(int tIdxMol, 
                       int tIdxAtm,
                       REAL tNumH);
        
        void reNameHAtoms(int tIdxMol);
        
        /*
        REAL  getTotalBondOrder(int tIdxMol,
                                std::vector<AtomDict>::iterator  tIA);
        REAL  getBondOrder(int tIdxMol,
                          int tIdx1, int tIdx2);
        */
        int  getNumOxyConnect(int tIdxMol, std::vector<AtomDict>::iterator iA);
        void setAtomsBondingAndChiralCenterMol(int tIdxMol);
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
    
    // should be added to separated chemistry part of source codes.
    /*
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
    
    extern REAL  getTotalBondOrder(Molecule   & tMol,
                                   std::vector<AtomDict>::iterator  tIA);
    extern REAL  getBondOrder(Molecule   & tMol,
                              int tIdx1, int tIdx2);
    
    extern void  setOneHAtomCoords();
    
    extern void  setOneHAtomCoordsSP3(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP2(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP(std::vector<AtomDict> & tAtoms,
                                     std::vector<AtomDiYct>::iterator tIA);
   */ 
    
    class SYBLMol2File: public File
    {
    
    public:
        
        // Default constructor 
        SYBLMol2File();
        
        // The constructor using a string as the file name
        SYBLMol2File(Name                        tFname,
                     std::ios_base::openmode     tOpenMode);
        
        // The constructor using an char* as the file name
        SYBLMol2File(FileName                    tFname,
                     std::ios_base::openmode     tOpenMode);
        
        // Default destructor
        ~SYBLMol2File();
        
        void setupSystem();
        void setBlock(std::string tL);
        
        
        
        // Get information from the input file
        void getAtomInfo(std::string tLine);
        void getBondInfo(std::string tLine);
        void getCrysInfo(std::string tLine);
        
        
        void setAtomsCCP4Type();
        
        void addAllHAtoms();
        void addHAtoms(int tIdxAtm,
                       REAL tNumH);
        
        void reNameHAtoms();
        
        /*
        REAL  getTotalBondOrder(int tIdxMol,
                                std::vector<AtomDict>::iterator  tIA);
        REAL  getBondOrder(int tIdxMol,
                          int tIdx1, int tIdx2);
        */
        
        int  getNumOxyConnect(std::vector<AtomDict>::iterator iA);
        // void setAtomsBondingAndChiralCenter();
        void setChiral();

        void iniDict();    
        
        
        
        bool                    hasCoords;
        bool                    hasConnect;
        bool                    hasH;
        bool                    containMetal;
  
        
        std::vector<AtomDict>   atoms;
        std::map<ID, ID>        atomSYBYLTypes;
        std::vector<AtomDict>   extraHAtoms;
        std::vector<BondDict>   bonds;
        std::vector<ChiralDict> chirals;

        std::ofstream           outFile;
        std::ifstream           inFile;
        
        std::map<std::string, bool>      mol2Dict; 
        
    };
    
    extern void accumInfoMols(ID                    tMolId,
                              std::vector<AtomDict>   & tAllAtomsOneMol,
                              std::vector<BondDict>   & tAllBondsOneMol,
                              std::vector<AngleDict>  & tAllAnglesOneMol,
                              std::map<ID, std::vector<ID> >    & tAllAtomTypes,
                              std::vector<std::string>          & tAllBondLines,
                              std::vector<std::string>          & tAllAngleLines);
    
}


#endif	/* MOLSDFFILE_H */

