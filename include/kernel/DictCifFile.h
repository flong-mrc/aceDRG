/* 
 * File:   CIFFile.h
 * Author: flong
 *
 * last updated on September 27, 2011, 12:51 AM
 */

#ifndef CIFFILE_H
#define	CIFFILE_H

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

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
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

#ifndef ATOMSTREE_H
#include "atomsTree.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef PDBFILE_H
#include "PDBFile.h"
#endif
/*
#ifndef  __MMDB_Manager__
#include "mmdb_manager.h"
#endif


#ifndef __MMDB_MMCIF__
#include <mmdb_mmcif_.h>
#endif

#ifndef __MMDB_IO__
#include <mmdb_io_file.h>
#endif
*/

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    class Residue;
    class Chain;
    class Model;
    
    class Bond;
    class BondDict;
    class Angle;
    class AngleDict;
    class Torsion;
    class TorsionDict;
    class Chiral;
    
    class Ring;
    class RinDict;
    
    class CrystInfo;
    
    class Link;
    
    class PeriodicTable;
    class CCP4AtomType;
    
    struct DictCifHead {
        Name   libName;
        ID     monType;
        ID     group;
        ID     libVersion;
        Date   libUpdate;
        DictCifHead():libName("?"),monType(""),group(""),libVersion("?"),libUpdate("?"){}
    };
    
    /*
    struct ChemComp {
        ID    id;
        Name  code;        // three-letter code
        Name  name;
        Name  group;       // non-polymer, peptide or DNA etc.
        int   numAtoms;   
        int   numH;
        ID    level;
     
        ChemComp():id(NullString), code(NullString), name(NullString),
        group(NullString), numAtoms(0), numH(0), level(NullString)
        {}         
    };
    */
    
    class GenCifFile : public File
    {
        
    public :
        
        // default constructor 
        GenCifFile();
        
        // copy constructor 
        GenCifFile(Name                    tFname, 
                   std::ios_base::openmode tOpenMode);
        
        GenCifFile(FileName                    tFname,
                   std::ios_base::openmode     tOpenMode);
        
        
        //GenCifFile(FileName                    tFname,
        //           mmdb::io:: tOpenMode);
        
        // destructor
        virtual ~GenCifFile();
        
        void setupSystem();
        void setupSystem2();
        void getPropsToMaps(std::vector<std::vector<std::string> >::iterator tOneBlockLines,
                            std::map<std::string, std::string>  & tRowProps,
                            std::map<int, std::map<ID, std::vector<std::string> > > & tColProps,
                            int   & tIdx);
        void selectPropsToMaps(std::map<std::string, std::string>  & tRowProps,
                               std::map<int, std::map<ID, std::vector<std::string> > > & tColProps);
        
        void initAllCifKeys();
        
        void getCifCrystInfo(std::vector<std::vector<std::string> >::iterator iBs);
        void getCifCrystInfo(std::map<std::string,std::string> & tRowProps, 
                             std::map<int,std::map<ID,std::vector<std::string> > > & tColProps);
        
        void getCifSymOps(std::vector<std::vector<std::string> >::iterator iBs);
        void getCifSymOps(std::map<ID,std::vector<std::string> >  & tOnePropGroup);
        
        void getCifAtomInfo(std::vector<std::vector<std::string> >::iterator iBs);
        void getCifAtomInfo2(std::vector<std::vector<std::string> >::iterator iBs);
        void getCifAtomInfo(std::map<ID,std::vector<std::string> >  & tOnePropGroup);
        void getAtomInfoFromLine(std::vector<std::string> & tStrs,
                                 int tP1, int tP2, int tP3, 
                                 int tP4, int tP5, int tPOcp);
        
        void setFlags(std::map<std::string, bool> & tL,
                      std::string tS);
        
        // get information from different blocks in the dictionary cif file
        void getHeadInfo(std::vector<std::string> tF);
        void getChemInfo(std::string tF);            //a different argument
                                                     // because symbol '
        void getDataDescription(std::vector<std::string> tF);
        
        // Atoms related 
        void getAtomInfo(std::vector<std::string> tF);
        void checkAtomElementID();
        void setAtomsCCP4Type();
        int  atomPosition(ID tID);
        void setAtomsPartialCharges();
        
        // Bonds related 
        void getBondAndConnect(std::vector<std::string> tF);
        void setHydroAtomConnect();
        void setBondOrder();
        void setSingleOrder();
        void setBondOrderByCoords();
        void setBondOrderByBondList();
        void setBondOrderByChem();
        void setBondOrderByType();
        
        void addAtomSeriNumToBond();
        
        void getAngleInfo(std::vector<std::string> tF);
        // get angles based on bonds and set their initial values
        void setAllAngles();
        
        void setAtomsCChemType();
        
        void setAtomsMetalType();
        
        void addMissHydroAtoms();
        
        void setAtomsVDWRadius();
       
        short transOrder(std::string tO);
        
        void  outSystem();
        void  outAtomBloc();
        void  outBondBloc();
    
        
        /* The core member functions are that 
           1. From a group of atoms, build a connection (or neighbor) list for
         *    each atom. (neighbAtoms in Atom class fulfill partly this purpose
         *    because COD-class codes include ring symbols, which need to be
         *    sorted out)
         * 2. Based on the connection lists, setup a COD-class code for each 
         *    atom.
         * 3. Search COD using corresponding COD-class codes of atoms to find 
         *    new bond values.
         * 4. The same applies for bond-angles, torsion-angles, chiral volumes 
         * 5. Output new CIF and restraint files   
         */
        
        void ringDetecting();
        void checkOnePathSec(std::vector<AtomDict> & seenAtoms,
                             std::map<int, ID>     & seenIDs,
                             AtomDict              & curAto,
                             int                   iMax,
                             std::vector<AtomDict>::iterator iAto); //a recursive function   
        void checkOnePathSec(AtomDict                & curAto,
                             int                       iMax,
                             std::vector<AtomDict>::iterator iOriAto,
                             int                       SeriNumPreAto,  
                             int                       curLev,
                             std::map<int, ID>       & seenAtomIDs,
                             std::map<int, ID>       & atomIDsInPath);
        
        void setAtomCodClassName(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void outRingSec(AtomDict & tAtom);
        std::string outRingSecStr(AtomDict &tAtom);
        
        
        // Plane related
        void detectPlaneGroups();
        void groupOrgAtomsToPlanes();
        void setSmallestPLs(std::vector<PlaneDict>  & tSPls);
        void mergeLargePLGroups(std::vector<PlaneDict> & tSmaPls);
        bool isInSameRing(PlaneDict & tP1, PlaneDict & tP2);
        bool furtherM(std::vector<int> &tV1, std::vector<int> &tV2);
        
        void CodAtomClassify(int dLev);
        
        void getAnglesFromPDB(ID tFName);
        void outRestraintCif(FileName tFName, ID tMonoRootName);
        void outPDB(FileName tFName, ID tMonoRootName);
        
        int                        curBlockLine;
        
        std::vector<CrystInfo>                              allCryst;
        
        std::map<ID, std::vector<std::string> >             allCifKeys;
        std::map<ID, bool>                                  existCifKeys;
        
        std::map<std::string, std::map<std::string, int> >  hasProps;
        
        DictCifHead                dictCifHead;
        ChemComp                   propComp;
        
        std::map<ID, std::string>  dataDesc;
        
        std::vector<AtomDict>      allAtoms;
        std::vector<int>           allHAtomIdx;
        std::vector<BondDict>      allBonds;
        std::vector<AngleDict>     allAngles;
        std::map<int, std::vector<std::vector<int> > > allAnglesIdxs;
        std::map<int, ID>                              DefaultCoordGeos;
        std::map<ID, std::map<int,ID> >                DefaultCoordGeos2;
              
        std::vector<TorsionDict>                 allTorsions;
        std::vector<ChiralDict>                  allChirals;
        std::vector<PlaneDict>                   allPlanes;
        std::map<ID, std::vector<RingDict> >     allRings;
        std::vector<RingDict>                    allRingsV;
        // std::vector<Ring>          allRings;    
        
        std::vector<int>           allHydroAtoms;
        
        std::vector<ID>            MetalTable;
        
        
        std::ofstream              outFile;
        std::ifstream              inFile;
        
        bool                       hasCoords;
        bool                       hasH;
        bool                       RFactorOK;
        
    private:
            
        int                                      itsCurAtomSeriNum;
        AtomDict                            *    itsCurAtom;
        std::string                              itsCurBlock;
        std::map<ID, std::vector<std::string> >  itsDataBlockMap;
        
        CrystInfo                           *    itsCurCryst;
        
            
    };
    
    
    
    class DictCifFile : public File
    {
        
    public :
        
        // default constructor 
        DictCifFile();
        
        // copy constructor 
        DictCifFile(Name                    tFname, 
                    std::ios_base::openmode tOpenMode);
        
       
        
        DictCifFile(FileName                    tFname,
                    std::ios_base::openmode     tOpenMode);
        
        // The following constructor is based on mmdb lib 
        //DictCifFile(FileName                    tFname);
        
        // For transfering information from PDB to CIF
        DictCifFile(FileName                    tCifName,
                    FileName                    tPdbName);
        
        // destructor
        virtual ~DictCifFile();
        
        void setupSystem();
        void setupSystem3Secs(std::ifstream & tInCif);
        
        void initCifKeys();
        
        void checkBloc(std::map<std::string, bool> & tL, 
                       std::vector<std::string> tF);
        void setFlags(std::map<std::string, bool> & tL,
                      std::string tS);
        
        // get information from different blocks in the dictionary cif file
        void getHeadInfo(std::vector<std::string> tF);
        void getChemInfo(std::string tF);            //a different argument
                                                     // because symbol '
        void getDataDescription(std::vector<std::string> tF);
        
        // Atoms related 
        void getAtomInfo(std::vector<std::string> tF);
        void setAtomsCCP4Type();
        int  atomPosition(ID tID);
        void setAtomsPartialCharges();
        
        
        // Bonds related 
        void getBondAndConnect(std::vector<std::string> tF);
        void setHydroAtomConnect();
        void setBondOrder();
        void setSingleOrder();
        void setBondOrderByCoords();
        void setBondOrderByBondList();
        void setBondOrderByChem();
        void setBondOrderByType();
        
        void addAtomSeriNumToBond();
        
        void getAngleInfo(std::vector<std::string> tF);
        // get angles based on bonds and set their initial values
        void setAllAngles();
        
        // This function read torsion angles from a cif file(not used)
        void setOneTorsion(std::vector<int>, REAL tValue, int tPeriod);
        void setTorsionIdxFromOneBond(int tIdx1, int tIdx2);
        void setTorsionFromOneBond(int tIdx1, int tIdx2);
        void setTorsionFromOneBond(int tIdx1, int tIdx2, std::string tFlip);
        void SetOneSP2SP2Bond(int tIdx1, int tIdx2);
        void SetOneSP2SP3Bond(int tIdx1, int tIdx2);
        void SetOneSP2SP3Bond(int tIdx1, int tIdx2, std::string tF);
        void SetOneSP3SP3Bond(int tIdx1, int tIdx2);
        void SetOneSP3SP3Bond(int tIdx1, int tIdx2, std::string tF);
        void SetOneSP3SP3Bond4H(int tIdx1, int tIdx2);
        bool checkSP3SP34H(int tIdx1, int tIdx2);
        void SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                    int tPer, REAL tIniValue);
        void SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                    int tPer, REAL tIniValue,
                                    std::string tF);
        void setAllTorsions();
        void setAllTorsions2();
        void setAllTorsionsInOneRing(std::vector<int> & tBondIdx, RingDict & tR);
        
        
        // This function get torsion angles from a set of bonds
        // not using the following one
        void getTorsionFromBonds(std::vector<BondDict> & tBonds);
        
        // Chiral-related 
        void getChiralInfo(std::vector<std::string> tF);
        
        void setAtomsBondingAndChiralCenter();
        void setAllChirals();
        
        void setDefaultCoordGeos();
        // user's definition of metal coordination geometries
        void getMetalCNGeo(std::vector<std::string> tF);
        
        void setAtomsCChemType();
        
        void setAtomsMetalType();
        
        void addMissHydroAtoms();
        
        void setAtomsVDWRadius();
       
        short transOrder(std::string tO);
        
        void  outSystem();
        void  outAtomBloc();
        void  outBondBloc();
        
        
        void transCoordsPdbToCif(DictPDBFile       & tPdbObj);
    
        
        /* The core member functions are that 
           1. From a group of atoms, build a connection (or neighbor) list for
         *    each atom. (neighbAtoms in Atom class fulfill partly this purpose
         *    because COD-class codes include ring symbols, which need to be
         *    sorted out)
         * 2. Based on the connection lists, setup a COD-class code for each 
         *    atom.
         * 3. Search COD using corresponding COD-class codes of atoms to find 
         *    new bond values.
         * 4. The same applies for bond-angles, torsion-angles, chiral volumes 
         * 5. Output new CIF and restraint files   
         */
        
        void ringDetecting();
        void checkOnePathSec(std::vector<AtomDict> & seenAtoms,
                             std::map<int, ID>     & seenIDs,
                             AtomDict              & curAto,
                             int                   iMax,
                             std::vector<AtomDict>::iterator iAto); //a recursive function   
        void checkOnePathSec(AtomDict                & curAto,
                             int                       iMax,
                             std::vector<AtomDict>::iterator iOriAto,
                             int                       SeriNumPreAto,  
                             int                       curLev,
                             std::map<int, ID>       & seenAtomIDs,
                             std::map<int, ID>       & atomIDsInPath);
        
        void setAtomCodClassName(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void outRingSec(AtomDict & tAtom);
        std::string outRingSecStr(AtomDict &tAtom);
        
        
        // Plane related
        void detectPlaneGroups();
        void groupOrgAtomsToPlanes();
        void setSmallestPLs(std::vector<PlaneDict>  & tSPls);
        void mergeLargePLGroups(std::vector<PlaneDict> & tSmaPls);
        bool isInSameRing(PlaneDict & tP1, PlaneDict & tP2);
        bool furtherM(std::vector<int> &tV1, std::vector<int> &tV2);
        
        void CodAtomClassify(int dLev);
        
        void getAnglesFromPDB(ID tFName);
        void outRestraintCif(FileName tFName, ID tMonoRootName);
        void outPDB(FileName tFName, ID tMonoRootName);
        
        
        

        int                        curBlockLine;
        
        
        std::map<ID, std::vector<std::string> >             allCifKeys;
        
        std::map<std::string, std::map<std::string, int> >  hasProps;
        
        DictCifHead                dictCifHead;
        ChemComp                   propComp;
        
        std::map<ID, std::string>  dataDesc;
        
        std::vector<AtomDict>      allAtoms;
        std::vector<int>           allHAtomIdx;
        std::vector<BondDict>      allBonds;
        std::vector<AngleDict>     allAngles;
        std::map<int, std::vector<std::vector<int> > > allAnglesIdxs;
        std::map<int, ID>                              DefaultCoordGeos;
        std::map<ID, std::map<int,ID> >                DefaultCoordGeos2;
              
        std::vector<TorsionDict>   allTorsions;
        std::vector<ChiralDict>    allChirals;
        std::vector<PlaneDict>     allPlanes;
        std::map<ID, std::vector<RingDict> >     allRings;
        std::vector<RingDict>                    allRingsV;
        // std::vector<Ring>          allRings;
        
        
        std::vector<int>           allHydroAtoms;
        
        std::vector<ID>            MetalTable;
        
        std::ofstream              outFile;
        std::ifstream              inFile;
        
        // For File transfer
        std::map<std::string, std::vector<std::string> >   allUnchangedBlocks;
        
        
        bool                       hasConnect;
        bool                       hasCoords;
        bool                       hasH;
        bool                       hasCCP4Type;
        
    private:
            
        int                        itsCurAtomSeriNum;
        AtomDict                *  itsCurAtom;
        int                        itsCurBondSeriNum;
        BondDict                *  itsCurBond;
        int                        itsCurAngleSeriNum;
        AngleDict               *  itsCurAngle;
        
        int                        itsCurTorsionSeriNum;
        TorsionDict             *  itsCurTorsion;
        
        int                        itsCurChiralSeriNum;
        ChiralDict              *  itsCurChiral;
            
        
    };
    
    
    
    
    
    extern void outMMCif(FileName tFName, 
                         ID tMonoRootName,
                         ChemComp  &         tPropComp,
                         std::vector<LIBMOL::AtomDict>& tAtoms,
                         // std::vector<int>    & tHydroAtoms,
                         std::vector<LIBMOL::BondDict>& tBonds, 
                         std::vector<LIBMOL::AngleDict>& tAngs, 
                         std::vector<LIBMOL::TorsionDict>& tTorsions, 
                         std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> > & tRings, 
                         std::vector<LIBMOL::PlaneDict>& tPlas, 
                         std::vector<LIBMOL::ChiralDict>& tChs);
    
    extern int getHAtomNum(std::vector<LIBMOL::AtomDict>& tAtoms);

    extern void outMMCif3Secs(FileName tFName, 
                              ID tMonoRootName,
                              std::vector<LIBMOL::AtomDict> & tAtoms,
                              std::map<std::string, std::vector<std::string> > & tUnChangedEntries);
    
    extern void outAtomTypesAndConnections(FileName tFName,
                                           std::vector<LIBMOL::AtomDict>& tAtoms,
                                           std::vector<LIBMOL::BondDict>& tBonds);
}


#endif	/* CIFFILE_H */

