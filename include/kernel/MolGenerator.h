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

#ifndef BOND_H
#include "bond.h"
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

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
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
        
        // Check if metal atoms exist in the assym unit
        void checkMetal(std::vector<ID> & tMeTab);
        
        // Generate a "molecule"(not normal sense) using a set of atoms and 
        // a set of symmetry operators
      
        void symmAtomGen(std::vector<CrystInfo>::iterator    tCrys,
                         PeriodicTable & tPTab);
        void symmAtomGen2(std::vector<CrystInfo>::iterator    tCrys,
                         PeriodicTable & tPTab);
        
        
        void getOneSymmAtom(std::vector<AtomDict>::iterator        tCurAtom,
                            std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
                            std::vector<CrystInfo>::iterator   tCrys, 
                            PeriodicTable & tPTab);
        void getOneSymmAtom2(std::vector<AtomDict>::iterator        tCurAtom,
                            std::map<std::string, std::vector<std::vector<REAL> > >::iterator tOp,
                            std::vector<CrystInfo>::iterator   tCrys, 
                            PeriodicTable & tPTab);
        
        void packAtomIntoCell(AtomDict & tAtm);
      
        // Create a system of atom including some of atoms in unit cells around 
        // the center unit cell.
        
        void buildRefAtoms(std::vector<CrystInfo>::iterator  iCryst);
        void addOneSetRefAtoms(AtomDict                         & tCurAtom,
                               std::vector<CrystInfo>::iterator   tCryst);
        void swithAtoms(std::vector<CrystInfo>::iterator tCryst);
        void cleanUnconnAtoms();
        void reIdxConnAtoms();
        
        void setUniqueAtomLinks(PeriodicTable & tPTab);
        void setUniqueAtomLinks(PeriodicTable & tPTab,
                                std::vector<CrystInfo>::iterator tCryst);
        void checkAtomLinks();
        
        void getMolByEqClassInCell();
        void getMolByEqClassInCrys();
        void deleteNonCenCellMols(std::map<unsigned, std::vector<int> >  
                                   & tMoleculesInCell);
        void deleteNonASUAtomCellMols(std::map<unsigned, std::vector<int> >  
                                      & tMoleculesInCell);
        
        bool inBonds(int tIdx1, int tIdx2, 
                     std::vector<BondDict> & tBonds);
        void getUniqueBonds(PeriodicTable & tPTab);
        void getUniqueAtomLinks(PeriodicTable & tPTab,
                            std::vector<CrystInfo>::iterator tCryst);
        void getUniqueBondsMols(Molecule    & tMol, 
                                std::vector<CrystInfo>::iterator tCryst);
        
        void getUniqueBondsMols2(Molecule    & tMol);
        
        void getAllBondsInOneMol(Molecule    & tMol, 
                                 std::vector<CrystInfo>::iterator tCryst);
        
        void getBondingRangePairAtoms(AtomDict & tAtm1,
                                      AtomDict & tAtm2,
                                      REAL tExtraD,
                                      PeriodicTable & tPTab,
                                      std::vector<REAL> & tRange);
        
        void getBondingRangePairAtoms2(AtomDict & tAtm1,
                                      AtomDict & tAtm2,
                                      REAL tExtraD,
                                      PeriodicTable & tPTab,
                                      std::vector<REAL> & tRange);
        
        void setOneUniqueBondCrys(int     tIdxAtm1,
                                  int     tIdxAtm2, 
                                  REAL    rD);
        
        void setOneUniqueBondCell(int     tIdxAtm1,
                                  int     tIdxAtm2, 
                                  REAL    rD);
        
        void checkBondOneForming(PeriodicTable & tPTab);
        
        REAL getBondLenFromFracCoords(std::vector<REAL> & tCoord1, std::vector<REAL> & tCoord2,
                                      REAL a, REAL b, REAL c, REAL alpha, REAL beta, REAL gamma);
        
        void getUniqAngles();
        void getUniqAngles(std::vector<CrystInfo>::iterator tCryst);
        void getUniqAngleMols(Molecule    & tMol,
                              std::vector<CrystInfo>::iterator tCryst);
        
        REAL getAngleValueFromFracCoords(AtomDict  & tAtCen,
                                         AtomDict  & tAt1, 
                                         AtomDict  & tAt2,
                                         REAL a, REAL b, REAL c, 
                                         REAL alpha, REAL beta, REAL gamma);
        
        
        // Validate all molecules 
        void buildAndValidMols(PeriodicTable & tPTab,
                           std::vector<CrystInfo>::iterator  iCryst);
        bool checEquiMoles(std::vector<Molecule> & tSetMols,
                           Molecule & tMol);
        void checkAtomElementID(std::vector<AtomDict> & tAtoms);
        bool colidAtom(AtomDict               & tAtom,
                       std::vector<AtomDict>  &  tRefAtoms, 
                       int tMode);
        bool colidAtom(std::vector<REAL>               & tFrcX,
                       std::vector<AtomDict>  &  tRefAtoms, int tMode);
        bool isASUAtomInMol(std::map<unsigned, std::vector<int> >::iterator tMol);
        bool connMetal(std::vector<int>      & tIdxs, 
                       std::vector<AtomDict> & tAtoms);
        bool connMetal2ndNB(std::vector<int>      & tIdxs, 
                       std::vector<AtomDict> & tAtoms);
        
        bool checkAsmAtomsInMol(Molecule  & tMol, std::string  & tErrInfo);
        bool checkAtomOcp(Molecule  & tMol, std::string  & tErrInfo);
        bool validateBonds(std::vector<BondDict>::iterator tBo, 
                           std::string & tErrInfo,
                           PeriodicTable & tPTab);
        bool validateBonds(std::vector<BondDict>::iterator tBo, 
                           Molecule    & tMol,
                           std::string & tErrInfo,
                           PeriodicTable & tPTab);
        bool validateBondValueSame(Molecule    & tMol,
                                   std::string & tErrInfo);
        bool validateBondValueDiff(Molecule    & tMol,
                                   std::string & tErrInfo);
        bool validateBondValueDiffStruct(std::vector<Molecule> & tMols,
                                         std::string & tErrInfo);
        
        bool validateAtomLinks(Molecule    & tMol,
                               PeriodicTable & tPTab,
                               std::string & tErrInfo);
        bool validateMolecule(Molecule    & tMol, PeriodicTable & tPTab,
                              std::string & tErrInfo);
        
        void checkInfMols(std::vector<Molecule> & aSetInfMols,
                          std::vector<Molecule> & aSetFinMols);
        
        bool checkOneMolInf(std::vector<Molecule>::iterator tMol);
        
        void setAtomNFormalCharge(Molecule & tMol);
        void getAtomTypeMols();
        void getAtomTypeOneMol(Molecule    & tMol);
        void getAtomTypeOneMolNew(Molecule & tMol);
        
        void getOverallBondAndAngles();
        void getOverallBondAndAnglesNew();
        
        void outTableMols(std::ofstream & tMolTabs, 
                          Molecule & tMol);
        void outTableBAndA(FileName tBAndAFName);
        void setTableSpAndChirals(Molecule & tMol, 
                                  std::map<std::string, std::map<std::string,
                                  std::map<std::string, std::map<std::string,
                                  REAL> > > > & tAtomSpAndChMap);
        
        void outTables(FileName tOutName, 
                         std::vector<Molecule> & tFinMols,
                         std::vector<Molecule> & tInfMols);
        
        void outMsg(FileName tOutName);
        void getOutFileRoot(FileName tOutName, Name & tRootName);
        void outMolsInfo(std::ofstream & tMolTabs,
                         std::vector<Molecule> & tFinMols,
                         std::vector<Molecule> & tInfMols);
        
        void contMetal2NB(int & tNB, int & tNA);
        
        // Metal Atom studies where molecules will not be generated
        void buildMetalAtomCoordMap(std::vector<CrystInfo>::iterator tCryst);
        bool checkNBAtomOccp(std::vector<AtomDict>::iterator tAtm);
        bool checkNBAtomOccp(AtomDict  & tAtm);
        int  getNumOrgNB(std::vector<AtomDict> & tAtoms, 
                         int  tIdx, 
                         std::vector<std::string> & tOrgTab);
        void outMetalAtomCoordInfo( FileName tOutName);
        void outMetalTables(FileName tOutName);
        
        void execute(FileName tOutName);
        void executeMet(FileName tOutName);
        
        
        std::string                     aLibmolTabDir;
        bool                            lColid;
        bool                            lHasMetal;
        
        std::vector<std::string>        allMsg;
        std::map<int, std::string>      validedMolMsg;
        std::map<int, std::string>      errMolMsg;
        
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
   
        
        std::map<unsigned, std::vector<int> >   moleculesInCell;
        std::map<unsigned, std::vector<int> >   moleculesInCryst;
        
        std::map<int, std::vector<int> >        bondsMetal2NB;
        std::map<int, std::vector<int> >        anglesMetal2NB;
        std::vector<int>                        noContriMols;
   
        // Metal atom studies 
        std::map<int, std::vector<int> >        connMapMetalAtoms;
        std::map<int, std::map<int, REAL> >     metalRelatedBonds;
        std::map<int, std::map<int, std::map<int, REAL> > > 
                                                metalRelatedAngles;
                 
        
    private :
        
        int                                                 myNBDepth;
        std::map<std::string, std::map<std::string, REAL> > shortLenList;
                
    };
       
}


#endif	/* MOLGENERATOR_H */

