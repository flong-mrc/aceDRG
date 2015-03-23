/* 
 * File:   codClassify.h
 * Author: flong
 *
 * Created on April 4, 2012, 5:48 PM
 */

#ifndef CODCLASSIFY_H
#define	CODCLASSIFY_H

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

#ifndef MOLECULE_H
#include "molecule.h"
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

#ifndef ALLSYSTEM_H
#include "AllSystem.h"
#endif 

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif

#ifndef CIFFILE_H
#include "DictCifFile.h"
#endif


namespace LIBMOL
{
    class AtomDict;
    class BondDict;
    
    class Ring;
    class RingDict;
    
    class Angle;
    class AngleDict;
    
    class Torsion;
    class TorsionDict;
    class Chiral;
    
    class Plane;
    class PlaneDict;
    
    class PeriodicTable;    
    
    class DictCifFile;
    class DictPDBFile;
    
    class AllSystem;
    

    struct aValueSet
    {
        REAL      value;
        REAL      sigValue;
        int       numCodValues;
    };
    
    
    // an AtomTree is a representative of atom types 
    class AtomTree 
    {
    public:
        
        //Default constructor 
        AtomTree();
        // Copy constructor 
        AtomTree(const AtomTree & tAtTree);
        
        //Constructor using a COD atom class
        AtomTree(ID & tCodClass);
        
        //Destructor 
        ~AtomTree();
        
        void sortedTree(); 
        
        int scoreParents();
        int scoreChilds();
        int scoreTree();
        
        ID                                 root;
        std::map<ID, std::vector<ID> >     family;
       
    };
    
    class AtomForest
    {
    public:
        
        // Default constructor 
        AtomForest();
        
        // Copy Constructor 
        AtomForest(const AtomForest & tAtForest);
        
        // Constructor using a vector of AtomTree
        AtomForest(std::vector<AtomTree> & tAtoms);
        
        //Destructor 
        ~AtomForest();
        
        // 
    };
    
    struct NB1stFam {
        std::string                   name;
        std::vector<std::string>      NB2ndList;
        int                           repN;
    };
    
    class CodClassify 
    {
    public :
        
        // Default constructor
        CodClassify();
        // Copy constructor
        // CodClassify(const CodClassify & tCodC);
        // Constructor initiated by a vector of AtomDicts
        CodClassify(const std::vector<AtomDict> & tAtoms);
        
        // 
        CodClassify(const std::vector<AtomDict>                 & tAtoms,
                    const std::vector<int>                      & tHAtomIdx, 
                    const std::vector<BondDict>                 & tBonds,
                    const std::vector<AngleDict>                & tAngles,
                    const std::vector<TorsionDict>              & tTorsions,
                    const std::vector<ChiralDict>               & tChirals,
                    const std::vector<PlaneDict>                & tPlans,
                    const std::map<ID, std::vector<RingDict> >  & tRings, 
                    std::string                                   tLibmolTabDir);
        
        CodClassify(const std::vector<AtomDict>                 & tAtoms,
                    const std::vector<int>                      & tHAtomIdx, 
                    const std::vector<BondDict>                 & tBonds,
                    const std::vector<AngleDict>                & tAngles,
                    const std::vector<TorsionDict>              & tTorsions,
                    const std::vector<ChiralDict>               & tChirals,
                    const std::vector<PlaneDict>                & tPlans,
                    const std::map<ID, std::vector<RingDict> >  & tRings, 
                    std::string                                 tLibmolTabDir,
                    int                                         nTM);
        
        // Constructor by a DictCifFile object
        CodClassify(const DictCifFile & tCifObj, 
                    std::string   tLibmolTabDir);
        
        // Constructor by a Mol/Sdf molecules
        CodClassify(const std::vector<Molecule> & tMols);
       
        // Constructor using a general system
        CodClassify(const AllSystem & tAllSys);
        
        // Destructor
        ~CodClassify();       
        
        
        void setupSystem();
        void setupSystem2();
        
        // get CCP4 energy lib 
        void getCCP4BondAndAngles();
        
        // General functions
        int atomPosition(ID tID);
        
        // Transfering between atoms and COD classes
        void codAtomClassify(int dLev);
        void codAtomClassify2(int dLev);
        void codAtomClassifyNew(int dLev);
        void codAtomClassifyNew2(int dLev);
        void codAtomClassifyNew3(int dLev);
        
        void getSmallFamily(std::string tInStr, NB1stFam & aNBFam);
        
        void codClassToAtom(ID & tCC, AtomDict & tAt);
        void codClassToAtom2(ID & tCC, AtomDict & tAt);
        
        
        void codClassToAtomAng(ID & tCC, AtomDict & tAt);
        void codClassToAtomAng2(ID & tCC, AtomDict & tAt);
        
        void codNBProps(std::vector<ID> tarStrs, std::vector<ID> & tCTs,
                        std::vector<int> & tNB2s, 
                        std::vector<int> & tRis, std::vector<int> & tPls);
        
        // the "distances" (or dissimilarities) between two atom
        int codAtomsDist(std::vector<ID> tarStrs, std::vector<ID> tNBs, int tLev);
        int codAtomsChemTypeDist(std::vector<ID>  & tarCTs,
                                          std::vector<int>  & tarRis, 
                                          std::vector<int>  & tarPls,
                                          std::vector<ID>  & bCTs,
                                          std::vector<int>  & bRis,
                                          std::vector<int>  & bPls);
        int codAtomsNumNBDist(std::vector<int> & tNB2s,
                              std::vector<int> & bNB2s);
        
        // Ring detecting related methods 
        void ringDetecting();
        void ringDetecting2();
        //a recursive function
        void checkOnePathSec(std::vector<AtomDict> & seenAtoms,
                             std::map<int, ID>     & seenIDs,
                             AtomDict              & curAto,
                             int                   iMax,
                             std::vector<AtomDict>::iterator iAto);    
        void checkOnePathSec(AtomDict                & curAto,
                             int                       iMax,
                             std::vector<AtomDict>::iterator iOriAto,
                             int                       SeriNumPreAto,  
                             int                       curLev,
                             std::map<int, ID>       & seenAtomIDs,
                             std::map<int, ID>       & atomIDsInPath);
        
        void checkOnePathSec2(AtomDict                & curAto,
                              int                       iMax,
                              std::vector<AtomDict>::iterator iOriAto,
                              int                       SeriNumPreAto,  
                              int                       curLev,
                              std::map<int, ID>       & seenAtomIDs,
                              std::map<int, ID>       & atomIDsInPath);
        
        void setAtomCodClassName(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void setAtomCodClassName2(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void setAtomCodClassNameNew(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void setAtomCodClassNameNew2(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void setAtomCodClassNameNew3(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void outRingSec(AtomDict & tAtom);
        void outRingSec2(AtomDict & tAtom);
        void outRingSecNew(AtomDict & tAtom);
        void outRingSecNew2(std::string & tAtmCodStr, AtomDict & tAtom);
        std::string outRingSecStr(AtomDict &tAtom);
        
        
        // Methods related to planarity generation (depending
        // also on if rings have been detected 
        
        void setAtomsBondingAndChiralCenter();
        
        // Plane-related
        void detectPlaneGroups();
        void groupOrgAtomsToPlanes();
        void groupMetAndLigandAtomsToPlanes();    // this could be use both early and late stage of a process
        void setSmallestPLs(std::vector<PlaneDict> & tSmaPls);
        void mergeLargePLGroups(std::vector<PlaneDict> & tSmaPls);
        bool checkATorsAtomsInPla(std::vector<int> & tAtms);
        
        // Ring-related 
        bool furtherM(std::vector<int> &tV1, std::vector<int> &tV2);
        bool isInSameRing(PlaneDict & tP1, PlaneDict & tP2);
        
        // Relations between torsion angles and ring
        
        
        // Atom hashing related 
        void readTablesForHashing(std::map<int, ID>  & tDigitKeys,
                                  std::map<int, int> & tLinkedHA);
        void readTablesForHashing2(std::map<int, ID>  & tDigitKeys,
                                   std::map<int, int> & tLinkedHA);
        void hashingAtoms();
        void hashingAtoms2();
        void setAtomsNBSymb(); 
        void setAtomsNBSymb2();
        
        // Tree-related 
        void setAtomsMST();
        
        // Bonds relate
        
        void setDefaultOrgBonds();
        void setOrgBondHeadHashList();
        void setOrgBondHeadHashList2();
        void groupCodOrgBonds();
        void groupCodOrgBonds2();
        void groupCodOrgBonds3();
        void searchCodOrgBonds(std::vector<BondDict>::iterator iOB);
        void searchCodOrgBonds2(std::vector<BondDict>::iterator iB);
        void groupCodMetBonds();
        void searchCodMetBonds(std::vector<BondDict>::iterator iMB);
        void setupTargetMetBondsUsingMean(std::vector<std::map<ID,REAL> > tSkeys,
                                          std::vector<BondDict>::iterator tB);
        void searchCodBonds();
        
        void getCCP4Bonds(std::vector<BondDict>::iterator tB, ID tAtm1, ID tAtm2);
        
        // !!!!!!! These are newly added for using sqlite3 search engine
        //void searchCodBondsUsingSqlite3();
        //void searchOneOrgBondFromCodUsingSqlite(sqlite3 * tCombDB,  
        //                                        std::vector<BondDict>::iterator tB);
       // bool searchOneOrgBondUsingSqliteL(sqlite3 * tCombDB, 
       //                                   std::vector<BondDict>::iterator tB,
       //                                   int  & tLev,
       //                                   std::map<ID, ID> tPropNB,
       //                                   std::map<ID, ID> tPropHash);
        void setQueStr(std::string & tQue, int tLev,
                       std::map<ID, ID>  tPropNB,
                       std::map<ID, ID> tPropHash);
        void setOneBondByMean(std::vector<std::vector<std::string> > & tResults,
                        std::vector<BondDict>::iterator tB);
        
        void setOneBondByDist(std::vector<std::vector<std::string> > & tResults,
                              std::vector<BondDict>::iterator tB);
        //End !!!!!!!
        
        void setupTargetBondsUsingSymblDist(std::vector<BondDict> & tBonds,
                                            std::vector<BondDict>::iterator tB,
                                            int tLev);
        void setupTargetBondsUsingSymblDist2(std::vector<BondDict> & tBonds,
                                            std::vector<BondDict>::iterator tB,
                                            int tAs0, int tAs1, int tLev);
        void constrBondSigmas();
        
        // Using meaning values when secondary NBs have the same "space"
        // such as "c[6,6]-3:n[5,6]-3:N[5]-2:"
        
        void setupTargetBondsUsingMean(std::vector<BondDict> & tBonds,
                                   std::vector<BondDict>::iterator tB);
        void setupTargetBondsUsingValueSetMean(std::vector<aValueSet> & tSets,
                                   std::vector<BondDict>::iterator tB);
        void setValueSet(aValueSet & tVs, std::vector<aValueSet> & tVecVs);
        void setupTargetBonds();
        void setupTargetBonds2();
        void setupTargetBondsUsingSqlite();
        
        void addAtomClassToBonds(std::vector<BondDict> & tBonds);
        
        // Bond angles related 
        void setDefaultOrgAngle();
        void addAtomClassToAngles(std::vector<AngleDict> &tAngles);
        
        void initTargetAngles();    // generate target angles based on atoms and bonds
        void setDefaultCoordGeos();
        void setOrgAngleHeadHashList();
        void setOrgAngleHeadHashList2();
        void setOrgAngleHeadHashList22();
        void groupCodOrgAngles();
        void groupCodOrgAngles2();
        void groupCodOrgAngles22();
        void searchCodOrgAngles(std::vector<AngleDict>::iterator iAN);
        void searchCodOrgAngles2(std::vector<AngleDict>::iterator iAN);
        void searchCodOrgAngles22(std::vector<AngleDict>::iterator iAN);
        bool searchCodOrgAnglesCen(std::vector<AngleDict>::iterator iAN, 
                                   int tHa1, int tHal2, int tHa3);
        bool getCCP4Angle(std::vector<AngleDict>::iterator tAN);
        
        // bond-angles involving metal atoms as NB atoms 
        void groupCodMetAngles();
        void getIdealCNGeoAngles(std::vector<AngleDict>::iterator iAN);
        void getCN10PentAntiPros(std::vector<AngleDict>::iterator iAN);
        REAL getOnebond(ID aID1, ID aID2);
        void searchCodMetAngles(std::vector<AngleDict>::iterator iAN);
        void groupCodAnglesWithNonCenteredMetal();
        ID   getMatchedKey(std::vector<ID> tKeys, ID tTarget, std::vector<int> tIdxs);
        void getMatched2dVect(std::vector<std::vector<ID> >tPairs, std::vector<ID> aPair, 
                             std::vector<int> tIdxs, int & rIdx);
        void getMatched3dVect(std::vector<std::vector<ID> >tPairs, std::vector<ID> aPair, 
                             std::vector<int> tIdxs, int & rIdx);
        void getMatched4dVect(std::vector<std::vector<ID> >tPairs, std::vector<ID> aPair, 
                             std::vector<int> tIdxs, int & rIdx);
        void getMatched5dVect(std::vector<std::vector<ID> >tPairs, std::vector<ID> aPair, 
                             std::vector<int> tIdxs, int & rIdx);
        void searchCodAnglesWithNonCenteredMetal(std::vector<AngleDict>::iterator iAN);
        void searchCodAngles();
        void searchCodAngles2();
        
        void setupTargetAngleUsingdist(std::vector<AngleDict> & tAngles,
                                       std::vector<AngleDict>::iterator tA,
                                       int tLev);
        void setupTargetAngleUsingdist2(std::vector<AngleDict> & tAngles,
                                        std::vector<AngleDict>::iterator tA,
                                        int as2, int as3, int tLev);
        
        void setupTargetAngleUsingdistAng(std::vector<AngleDict> & tAngles,
                                       std::vector<AngleDict>::iterator tA,
                                       int tLev);
        
        void setupTargetAngleUsingMean(std::vector<AngleDict> & tAngles,
                                       std::vector<AngleDict>::iterator tA);
        
        void setupTargetAngleUsingMean2(std::vector<AngleDict> & tAngles,
                                        std::vector<AngleDict>::iterator tA,
                                        int tHa1, int tHa2, int tHa3,
                                        ID tA1NB2, ID tA2NB2, ID tA3NB2);
        
        void setupTargetAngleUsingMean3(std::vector<REAL>   & tAngValues,
                                        std::vector<int>    & tAngNums,
                                        std::vector<AngleDict>::iterator tA);
        
        void setupTargetAngles();
        void setupTargetAngles2();
        
        // angles using sqlite3
        //void setupTargetAnglesUsingSqlite();
        //void searchCodAnglesUsingSqlite();
        //void searchOneOrgAngleFromCodUsingSqlite(sqlite3 * tCombDB,
        //                                         std::vector<AngleDict>::iterator tA);
        
        //bool searchOneOrgAngleUsingSqliteL(sqlite3 * tCombDB, 
        //                                  std::vector<AngleDict>::iterator tA,
        //                                  int            & tLev,
        //                                  std::map<ID, ID> tPropNB,
        //                                  std::string tA1C);
        
        void setQueStrAn(std::string & tQue, int tLev,
                         std::map<ID, ID>  tPropNB, 
                         std::string tA1C);
        void setOneAngleByMean(std::vector<std::vector<std::string> >& tResults, 
                               std::vector<AngleDict>::iterator tA);
        
        void checkAngConstraints();
        void checkSP2Constraints(std::vector<int> tAngIdxs);
        void checkSP3Constraints(std::vector<int> tAngIdxs);
        void checkRingAngleConstraints(std::vector<RingDict>::iterator iRv);
        bool checkSpeAng(std::vector<AngleDict>::iterator tA);
        void setSpecialAngles(std::map<int, std::vector<AngleDict> > & tAngs);
        void setOneSetBoronAngles(std::map<int, std::vector<AngleDict> >::iterator tAs);
        void setOneSetCarbonAngles(std::map<int, std::vector<AngleDict> >::iterator tAs);
        
        // Torsion angles related 
        void setupTargetTorsions();
        void fixTorIDs();
        // 
        void setupAllTargetValues();
        void setupAllTargetValues2();
        
        void setupAllStdValues();
        
        
        // Initiate a set of dummy atoms for (the torsion space)
        void initDummyAtoms();
        // 
        int checkNumInLevel();
        
        void goToLeveL(int tLev);
        
        // Output 
        void initRoughCoords();
        void initRoughCoordsTor();
        void outRestraintCif(FileName tFName);
        void outRestraintCif(FileName tFName, ID tMonoRootName);
        // temporally for outRestraintCif2
        void outRestraintCif2(FileName tFName, ID tMonoRootName);
        void outPDB(FileName tFName);
        void outPDB(FileName tFName, ID tMonoRootName);
        void outAtomTypes(ID tMonoRootName);
        void getAnglesFromPDB(ID tFName);
        
        
        void setSpecial3NBSymb(std::vector<AtomDict>::iterator tAt);
        void setSpecial3NBSymb2(std::vector<AtomDict>::iterator tAt);
        
        // This should be moved into another specifically designed sqlite3 class
        //void sqlite3Query(sqlite3 *      tDB,
        //                  SqliteStatment tQue, 
        //                  std::vector<std::vector<std::string> >  & tResult);
        
        
        
        int                                      wSize;
        std::string                              libmolTabDir;
        
        PeriodicTable            *               pPeriodictable;
        
        
        std::vector<AtomDict>                    allAtoms;
        std::vector<int>                         allHAtomIdx;
        
        
        // bonds
        std::vector<BondDict>                    allBonds;
        
        std::vector<BondDict>                    allDictBonds;
        std::map<int, std::map<int, std::map<ID, std::map<ID,  std::map<ID,
                 std::map<ID, std::map<ID, std::map<ID, int > > > > > > > >allDictBondsIdx;
        
        
        std::map<int, std::map<int, std::map<ID,  std::map<ID, 
                std::map<ID,  std::map<ID, 
                std::vector<aValueSet> > > > > > > allDictBondsIdx1;
        
        std::map<int, std::map<int, std::map<ID, 
                std::map<ID, std::vector<aValueSet> > > > > allDictBondsIdx2;
        
        std::map<int, std::map<int, 
                 std::vector<aValueSet> > >                 allDictBondsIdx3;
        
        std::map<ID, std::map<ID, std::map<int, 
                     std::map<int, std::map<ID, REAL> > > > >  allDictPreMetBonds;
        std::map<ID, std::map<ID, std::map<int, std::map<int, 
                     std::map<ID, std::map<ID, REAL> > > > > > allDictMetBonds;
        
        //DB3 
        std::vector<aValueSet>                    allDictBondsD;
        std::map<int, std::map<int, std::map<ID, std::map<ID,  std::map<ID,
                 std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                 std::map<ID, std::map<ID, int > > > > > > > > > > >allDictBondsIdxD;
        
        
        std::map<int, std::map<int, std::map<ID,  std::map<ID, 
                std::map<ID,  std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                std::vector<aValueSet> > > > > > > > > >allDictBondsIdx1D;
        
        std::map<int, std::map<int, std::map<ID,  std::map<ID, 
                std::map<ID,  std::map<ID, std::map<ID,   
                std::vector<aValueSet> > > > > > > >   allDictBondsIdx2D;
        
        std::vector<aValueSet>                         allDictAnglesD;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::map<ID,  std::map<ID,  std::map<ID,
        std::map<ID,  std::map<ID,  std::map<ID,
        int > > > > > > > > > > > > > > >              allDictAnglesIdxD;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::map<ID,  std::map<ID,  std::map<ID,
        std::vector<aValueSet> > > > > > > > > > > > > allDictAnglesIdx1D;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::vector<aValueSet> > > > > > > > > >       allDictAnglesIdx2D;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::vector<aValueSet> > > > > > >             allDictAnglesIdx3D;
        
        
        // angles 
        std::vector<AngleDict>                         allAngles;
        std::map<int, std::vector<std::vector<int> > > allAnglesIdxs;
        
        std::vector<AngleDict>                         allDictAngles;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::map<ID,  std::map<ID,  std::map<ID,
        int > > > > > > > > > > > >                    allDictAnglesIdx;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::vector<aValueSet> > > > > > > > > >            allDictAnglesIdx1;
        
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::vector<aValueSet> > > > > > >                  allDictAnglesIdx2;
        
        std::map<int, std::map<int, std::map<int, 
        std::vector<aValueSet> > > >                        allDictAnglesIdx3;
        
        std::map<int, REAL>                                 DefaultOrgAngles;
        
        std::map<int, ID>                                   DefaultCoordGeos;
        std::map<ID, std::map<int,ID> >                     DefaultCoordGeos2;
        std::map<int, std::map<ID, std::vector<REAL> > >    allDictCoordGeoAngs;
        std::map<ID,  std::map<ID, std::map<ID,
        std::map<ID, std::map<ID, std::map<ID, 
        int> > > > > >                                      allDictNonCenMetAnglesIdx;
        
        
        
        std::vector<TorsionDict>                            allTorsions;
        std::vector<ChiralDict>                             allChirals;
        std::vector<PlaneDict>                              allPlanes;
        std::map<ID, std::vector<RingDict> >                allRings;
        
        std::vector<int>                                    allHydroAtoms;
        std::vector<AtomDict>                               allDummyAtoms;
        
        std::map<int, ID>                                   codOrgBondFiles;
        std::map<ID, ID>                                    codOrgBondFiles2;
        std::map<int, ID>                                   codOrgAngleFiles;
        std::map<ID, ID>                                    codOrgAngleFiles2;
        
        std::map<ID, std::map<ID, std::map<ID, REAL> > >    ccp4Bonds;
        std::map<ID, std::map<ID, std::map<ID, REAL> > >    ccp4Angles;
    };
    
    class CodBonds 
    {
    public:
        
        // Default constructor
        CodBonds();
        // Copy constructor
        CodBonds(const CodBonds & tCBonds);
        // Constructor with an input file name
        CodBonds(FileName                    tFname,
                 std::ios_base::openmode     tOpenMode);
        // Destructor 
        ~CodBonds();
        
        
        void setupSystem();
        
        void sortCodTable(FileName    tFname,
                          std::vector<sortLine>  & allSortLines);
        void setupSystem(std::vector<sortLine> & tAllLines);
        
        void getTargetBonds(std::vector<BondDict> & targetBs);
        
        
        std::vector<BondDict>              allBonds;
        std::map<ID, std::vector<REAL> >   allBondsMap; // May use this one.
                                                        // instead of allBonds
        
        std::ofstream               outFile;
        std::ifstream               inFile;
                
    };
    
    class CodAngles
    {
    public:
        
        // Default constructor
        CodAngles();
        
        //Copy constructor
        CodAngles(const CodAngles & tCAngles);
        
        //Constructor using with an input file name
        CodAngles(FileName                    tFname,
                 std::ios_base::openmode      tOpenMode);
        
        // Destructor 
        ~CodAngles();
         
        
        void setupSystem();
        
        void getTargetAngles(std::vector<AngleDict> & targetAs);
 
        std::vector<AngleDict>              allAngles;
  
        std::ofstream                       outFile;
        std::ifstream                       inFile;
        
        
    };
    
    
    class CodTorsions
    {
    public :
        
        //Default constructor 
        CodTorsions();
        
        //Copy constructor 
        CodTorsions(const CodTorsions & tCT);
        
        // Constructor using a file name 
        CodTorsions(FileName                    tFname,
                   std::ios_base::openmode     tOpenMode);
        
        // Destructor
        ~CodTorsions();
        
        void setupSystem();
        
        void getTargetTorsions(std::vector<TorsionDict> & targetTs);
        
        std::vector<TorsionDict>        allTorsions;
        
        std::ofstream                   outFile;
        std::ifstream                   inFile;       
        
    };
    
    
    
    
}



#endif	/* CODCLASSIFY_H */

