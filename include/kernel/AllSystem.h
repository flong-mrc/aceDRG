/* 
 * File:   AllSystem.h
 * Author: flong
 *
 * Created on July 5, 2012, 10:36 AM
 */

#ifndef ALLSYSTEM_H
#define	ALLSYSTEM_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef CIFFILE_H
#include "DictCifFile.h"
#endif

#ifndef PDBFILE_H
#include "PDBFile.h"
#endif

#ifndef MOLSDFFILE_H
#include "MolSdfFile.h"
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

#ifndef MOLECULE_H
#include "molecule.h"
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

#ifndef CODCLASSIFY_H
#include "codClassify.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    
    class Bond;
    class BondDict;
    
    class Ring;
    class RingDict;
    
    class Angle;
    class AngleDict;
    
    class Torsion;
    class TorsionDict;
    
    class Chiral;
    class ChiralDict;
    
    class Plane;
    class PlaneDict;
    
    class Molecule;
    
    class MolSdfFile;
    
    class DictCifFile;
    class DictPDBFile;
    
    class CodClassify;
    
    class CCP4AtomType;
    
    class AllSystem
    {
    public:
        
        // Default constructor
        AllSystem();
        
        // Copy constructor 
        AllSystem(const AllSystem  & tAllSys);

        // No need for copy constructor.
        AllSystem(const std::vector<AtomDict>              & tAllAtoms,
                  const std::vector<int>                   & tAllHAtomIdx,
                  const std::vector<BondDict>              & tAllBonds,
                  const std::vector<AngleDict>             & tAllAngles,
                  const std::vector<TorsionDict>           & tAllTorsions,
                  const std::vector<ChiralDict>            & tAllChirals,
                  const std::vector<PlaneDict>               & tAllPlanes,
                  const std::map<ID, std::vector<RingDict> > & tAllRings);
        
        // Input from files contains mainly atom id and chemical types,
        // bonds and connections, and perhaps chirals of atoms (e.g, mol/sdf file).
        // Angle and torsion angles formed from atoms will be re-found here,
        // so as planarity
        
        AllSystem(DictCifFile & tCifObj);
        AllSystem(Molecule & tMol);
        // AllSystem(MolSdfFile & tSdfObj);
        
        // get the system after COD classification and search.
        AllSystem(const CodClassify  & tProCodSys);
        
        //Destructor
        ~AllSystem();
        
        
        // System builder 
 
        void AddAtoms(const std::vector<AtomDict>         & tAllAtoms);
        void AddAtom(const AtomDict                       & tAtom);
        void AddBonds(const std::vector<BondDict>         & tAllBonds);
        void AddBond(const BondDict                       & tBond);
        void AddAngles(const std::vector<AngleDict>       & tAllAngles);
        void AddAngle(const AngleDict                     & tAngle);
        void AddTorsions(const std::vector<TorsionDict>   & tAllTorsions);
        void AddTorsion(const TorsionDict                 & tTorsion);
        void AddChirals(const std::vector<ChiralDict>     & tAllChirals);
        void AddChiral(const  ChiralDict                  & tChiral);
        void AddPlanes(const std::vector<PlaneDict>       & tAllPlanes);
        void AddPlane(const PlaneDict                     & tPlane);
        void AddRings(const std::map<ID, std::vector<RingDict> > & tAllRings);
        void AddRing(const RingDict                              & tRing);
        
        void resetSystem(CodClassify & tCodSys);
        
        // Further setup angles, torsions, planes and chirals
        bool isOrgSys();
        int  atomPosition(ID tID);
        void setSysProps();
        void setHydroAtomConnect();
        void addMissHydroAtoms();
        int  getNumSpecAtomConnect(int i, ID tChemType);
        int  getNumSpecAtomConnect(std::vector<AtomDict>::iterator iA,
                                   ID tChemType);
        int  getNumOxyConnect(std::vector<AtomDict>::iterator iA);
        REAL getTotalOrderOneAtom(std::vector<AtomDict>::iterator iA);
        void setAtomsBondingAndChiralCenter();
        void setAtomsCChemType();
        void setAtomsCCP4Type();
        const bool containMetal();
        void setAtomsMetalType();
        void setAtomsPartialCharges();
        
        
        // Angles related 
        void setAllAngles();
        void checkAngConstraints();
        void checkSP2Constraints(std::vector<int> tAngIdxs);
        void checkSP3Constraints(std::vector<int> tAngIdxs);
        
        // These functions are  torsion angles related  
        void setOneTorsion(std::vector<int>, REAL tValue, int tPeriod);
        void setTorsionIdxFromOneBond(int tIdx1, int tIdx2);
        void setTorsionFromOneBond(int tIdx1, int tIdx2);
        void setTorsionFromOneBond(int tIdx1, int tIdx2, std::string tFlip);
        void SetOneSP1SP1Bond(int tIdx1, int tIdx2);
        void SetOneSP1SP2Bond(int tIdx1, int tIdx2);
        void SetOneSP1SP3Bond(int tIdx1, int tIdx2);
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
        
        
        /* Torsion related 
        void setAllTorsions();
        void getTorsionFromBonds(std::vector<BondDict> & tBonds);
 
        void setTorsionIdxFromOneBond(int tIdx1, int tIdx2);
        void setTorsionFromOneBond(int tIdx1, int tIdx2);
        void SetOneSP2SP2Bond(int tIdx1, int tIdx2);
        void SetOneSP2SP3Bond(int tIdx1, int tIdx2);
        void SetOneSP3SP3Bond(int tIdx1, int tIdx2);
        void SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                    int tPer, REAL tIniValue);
        void setOneTorsion(std::vector<int>, REAL tValue, int tPeriod);    
        */
        
        void ringDetecting();
        void checkOnePathSec(AtomDict                & curAto,
                             int                       iMax,
                             std::vector<AtomDict>::iterator iOriAto,
                             int                       SeriNumPreAto,  
                             int                       curLev,
                             std::map<int, ID>       & seenAtomIDs,
                             std::map<int, ID>       & atomIDsInPath);
        
        void outRingSec(AtomDict & tAtom);
        std::string outRingSecStr(AtomDict &tAtom);

        // Chiral-related
        void chiralExch();
        
        // Plane related
        void detectPlaneGroups();
        void groupOrgAtomsToPlanes();
        void setSmallestPLs(std::vector<PlaneDict>  & tSPls);
        void mergeLargePLGroups(std::vector<PlaneDict> & tSmaPls);
        bool isInSameRing(PlaneDict & tP1, PlaneDict & tP2);
        bool furtherM(std::vector<int> &tV1, std::vector<int> &tV2);     
        
        // COD database related 
        void setAtomCodClassName(AtomDict &tAtom, AtomDict &tOriAtom, int tLev);
        void CodAtomClassify(int dLev);
        
        // metal-related 
        void setDefaultCoordGeos();
        
        // COD applications 
        void setupAllTargetValuesFromCOD(ID tOutName, ID tMonoName);
        // Other applications 
        
        //void SetupCoords();
        //void RingFusion();
        //void SetupSpanTree();
        void GeoOpt();
        
        // output restraints in cif format
        void OutputRestraintCif(FileName tCifName);
        
        
        // output atomic coordinates in a file of PDB format
   
        // member variables
        bool                                     hasCoords;
        bool                                     hasCCP4Type;
        
        std::vector<AtomDict>                    allAtoms;
        std::vector<int>                         allHAtomIdx;      //repeated ones
        std::vector<int>                         allHydroAtoms;
        std::vector<BondDict>                    allBonds;
        std::vector<AngleDict>                   allAngles;
        std::map<int, std::vector<std::vector<int> > > allAnglesIdxs; 
              
        std::vector<TorsionDict>                 allTorsions;
        std::vector<std::vector<TorsionDict> >   allTorsionSets;
        std::vector<ChiralDict>                  allChirals;
        std::vector<PlaneDict>                   allPlanes;
        
        std::map<ID, std::vector<RingDict> >     allRings;
        std::vector<RingDict>                    allRingsV;
        
        std::vector<ID>                          MetalTable;
        std::map<int, ID>                        DefaultCoordGeos;
        std::map<ID, std::map<int,ID> >          DefaultCoordGeos2;
        
        std::vector<AtomDict>                    allDummyAtoms;
        
    private:
        
        bool                                     itsContainMetal;
        
        int                                      itsCurAngleSeriNum;
        AngleDict                             *  itsCurAngle;
        
        int                                      itsCurTorsionSeriNum;
        TorsionDict                           *  itsCurTorsion;
        
        int                                      itsCurChiralSeriNum;
        ChiralDict                            *  itsCurChiral;
        
    };

    void setSubSys(std::vector<int> tDoneSet,
                   std::vector<AtomDict>    & tAtoms,
                   std::vector<BondDict>    & tBonds,
                   std::vector<AngleDict>   & tAngles,
                   std::vector<TorsionDict> & tTorsions,
                   std::map<ID, std::vector<RingDict> > & tRings,
                   std::vector<PlaneDict>   & tPlas,
                   std::vector<ChiralDict>  & tChs,
                   AllSystem                & tSubSys);
    
    void extern optSubsys(AllSystem  & tSubSys);
       

}

 

#endif	/* ALLSYSTEM_H */

