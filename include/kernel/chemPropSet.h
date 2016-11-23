/* 
 * File:   chemPropSet.h
 * Author: flong
 *
 * Created on July 3, 2014, 1:27 PM
 */

#ifndef CHEMPROPSET_H
#define	CHEMPROPSET_H

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

#ifndef TRANSCOORD_H
#include "TransCoord.h"
#endif

#ifndef MOLSDFFILE_H
#include "MolSdfFile.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
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
    
    class CodClassify;
    
    class CCP4AtomType;
    
    class MolSdfFile;
    
    class PeriodicTable;
    
    extern bool assignElementType(PeriodicTable & tP, std::string tStr,  
                                  std::vector<AtomDict>::iterator tAtom);
    
    extern int  getNumOxyConnect(std::vector<AtomDict>  &  tAtoms,
                                std::vector<AtomDict>::iterator iA);
    
    extern void getHydroAtomConnect(std::vector<AtomDict>  &  tAtoms);
    
    extern void mdChiralByClasses(std::vector<AtomDict>::iterator iAt,
                                  std::vector<AtomDict>     &     tAtoms);
    // Two stages for sp hybridization (They will be merged to the same
    // function once the second one passes tests
    // 1. Currently used one (originally in codClassify.h and .cpp) 
    extern void setAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms);
    // 2. For modAtomsBondingAndChiralCenter
    // tMode 0 : not use coordinates (in dictionary generation stage) 
    // tMode 1 : definitely have reliable coordinates (in DB molecule generation stage)
    extern void modAtomsBondingAndChiralCenter(std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds,
                                               std::vector<AngleDict> & tAngles,
                                               std::vector<RingDict> & tRings,
                                               int                     tMode);
    
    extern void setAtomsNB1NB2_SP(std::vector<AtomDict> & tAtoms);
    extern void setAtomsNB1NB2_exElectrons(std::vector<AtomDict> & tAtoms);
    
    extern void setBondsAndAngles_NB1NB2_SP(std::vector<AtomDict> & tAtoms,
                                            std::vector<BondDict> & tBonds,
                                            std::vector<AngleDict> & tAngles);
    
    extern void setBondsAndAngles_NB1NB2_EE(std::vector<AtomDict> & tAtoms,
                                            std::vector<BondDict> & tBonds,
                                            std::vector<AngleDict> & tAngles);
    
    extern void reIndexAtomInRing(std::vector<AtomDict> & tAtoms,
                                  std::vector<RingDict> & tRings);
    
    extern void setAnglesSPSigns(std::vector<AtomDict> & tAtoms,
                                 std::vector<AngleDict> & tAngles);
    
    extern bool confirmPlaneByChiralVol(std::vector<AtomDict> & tAtoms,
                                        std::vector<AtomDict>::iterator tA);
    
    extern bool confirmPlaneByAngle(std::vector<AtomDict> & tAtoms,
                                    std::vector<AtomDict>::iterator tA,
                                    REAL                    tCri);
    
    // a function transfer int sp to a string
    extern std::string strTransSP(int tSP);
    
    extern bool checkBridgeStruct(std::vector<AtomDict> & tAtoms,
                                  std::vector<RingDict> & tRings,
                                  int                     tAnchIdx);
    
    // 
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
    
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  Molecule   & tMol, 
                                  PeriodicTable & tTab);
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  std::vector<AtomDict>   & tAtoms,
                                  std::vector<BondDict>   & tBonds,
                                  PeriodicTable & tTab);
    // The following function replaces the previous function addHAtomToMols().
    extern void ProtonateFunctionGroupInOneMol(std::vector<AtomDict>  & tA,
                                               std::vector<BondDict>  & tBonds,
                                               std::vector<int>       & tAddHIdx); 
    
    extern REAL checkProtonateAll(std::vector<AtomDict>::iterator tIA, 
                                  Molecule   & tMol, 
                                  PeriodicTable & tTab);
    
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
    
    // Function groups at the moment. 
    extern void checkProtonatedCarBoxylicTerminus(std::vector<AtomDict>::iterator tIA,
                                                  std::vector<AtomDict> & tAtoms,
                                                  std::vector<BondDict> & tBonds,
                                                  REAL tPH, REAL tPka,
                                                  std::vector<int> tDoneAtoms);
    
    extern void checkProtonatedSulfuricAcids(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH, std::vector<int> tDoneAtoms);
    
    extern void checkProtonatedAminoTerminus(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH, REAL tPka,
                                           std::vector<int> tDoneAtoms);
    
    extern void checkProtonatedPAcids(std::vector<AtomDict>::iterator tIA,
                                           std::vector<AtomDict> & tAtoms,
                                           std::vector<BondDict> & tBonds,
                                           REAL tPH, std::vector<int> tDoneAtoms);
    
    // The following function overwrites addHAtoms in MolSdfFile class  
    extern void adjustHAtoms(std::vector<AtomDict> & tAtoms,
                          std::vector<BondDict> & tBonds,
                          int                     tAtmIdx,
                          int                     tAddH,
                          std::vector<int>      & tHAtmIdxs);
    
    extern void setAllAddedHAtomCoords(std::vector<AtomDict> & tAtoms,
                                       std::vector<int>      & tHAtmIdxs);
    
    // Bond order sections 
        
    extern REAL  getTotalBondOrder(Molecule   & tMol,
                                   std::vector<AtomDict>::iterator  tIA);
    extern REAL  getTotalBondOrder(std::vector<BondDict>   & tBonds, 
                                   std::vector<AtomDict>   & tAtoms,
                                   std::vector<AtomDict>::iterator tIA);
    
    extern void modifyBondOrderAR(std::vector<BondDict> & tBonds,
                                  std::vector<AtomDict>  & tAtoms,
                                  int  tIdxB1, int tIdxB2,
                                  int tAtCen, int tAt1, int tAt2,
                                  PeriodicTable & tTab);
    
    extern REAL  getFixedBondOrder(std::vector<BondDict>   & tBonds, 
                                    std::vector<AtomDict>   & tAtoms,
                                    int                       tAtmIdx);
    
    extern REAL  getBondOrder(Molecule   & tMol,
                              int tIdx1, int tIdx2);
    
    extern REAL  getBondOrder(std::vector<BondDict> & tBonds,
                              int tIdx1, int tIdx2);
    
    extern REAL  getBondOrderOneAtom(std::vector<BondDict> tBonds,
                                     std::vector<AtomDict> tAtoms,
                                     int tIdx1, int tIdx2);
    extern void  setAllBondOrders(std::vector<AtomDict> & tAtoms,
                                 std::vector<BondDict> & tBonds);
    
    extern void  kekulizeRings(std::vector<AtomDict> & tAtoms,
                               std::vector<BondDict> & tBonds,
                               std::vector<RingDict> & tRings);
    extern void  kekulizeOneRing(std::vector<AtomDict> & tAtoms,
                                 std::vector<BondDict> & tBonds,
                                 RingDict & tRing, int tStartIdx,
                                 PeriodicTable & tTab);
    
    extern void  setOneHAtomCoords();
    
    extern void  setOneHAtomCoordsSP3(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP2(std::vector<AtomDict> & tAtoms,
                                      std::vector<AtomDict>::iterator tIA);
    extern void  setOneHAtomCoordsSP(std::vector<AtomDict> & tAtoms,
                                     std::vector<AtomDict>::iterator tIA);
    
    extern void setAllAtomEXcessElectrons(std::vector<AtomDict> & tAtoms);
    extern void setAllAtomEXcessElectrons2(std::vector<AtomDict> & tAtoms);
                
    extern void setAtomRingProps(std::vector<AtomDict> & tAtoms,
                                 std::vector<RingDict> & tRings);
    
    extern void setInitBondOrdersViaExtraElecs(std::vector<AtomDict> & tAtoms,
                                               std::vector<BondDict> & tBonds);
                
    extern int  sumExElectrons(std::vector<AtomDict> & tAtoms);
    
    
    // atom and bond stereo  
    extern void checkAllStereo(FileName tMdlIn, FileName tPdbIn,
                               FileName tPdbOut);
    extern void checkStereoOneMol(std::vector<Molecule>::iterator tMol,
                                  FileName tPdbOut);
    
    // Check aromatic stability of extended ring systems.
    // using it to decide planes in the systems(molecules or monomers)
    
    extern void setAromPlanes(std::vector<AtomDict> & tAtoms,
                              std::vector<RingDict> & tRings,
                              std::vector<PlaneDict> & tPlans);
    
    class HuckelMOSuite
    {
    public :
        
        HuckelMOSuite();
        ~HuckelMOSuite();
        
        void setWorkMode(int tMode);
        void execute(std::vector<AtomDict> & tAtoms, 
                     std::vector<BondDict> & tBonds);
        void execute2(std::vector<AtomDict> & tAtoms, 
                      std::vector<BondDict> & tBonds,
                      std::vector<RingDict> & tRings);
        
        void initiaExElecs(std::vector<AtomDict> & tAtoms);
        void initiaExElecs2(std::vector<AtomDict> & tAtoms);
        void PickPiElectrons(std::vector<AtomDict> & tAtoms);
        void PickOddAtoms(std::vector<AtomDict> & tAtoms);
        
        void setInitBondOrder(std::vector<AtomDict> & tAtoms,
                              std::vector<BondDict> & tBonds,
                              std::vector<int>      & tCBondIdx,
                              std::map<int, std::vector<int> > & tDelConn,
                              std::map<int, int>    & tRemainval);
        void setProBondOrdersOneLoop(int & nDone, std::vector<AtomDict> & tAtoms, 
                                     std::vector<BondDict> & tBonds, 
                                     std::vector<int> & tCBondIdx, 
                                     std::map<int, std::vector<int> > & remainConns, 
                                     std::map<int, std::vector<int> > & tDelConn, 
                                     std::map<int, int>   &  tRemainVal);
        void setBondOrderInSys(std::vector<AtomDict> & tAtoms,
                               std::vector<BondDict> & tBonds);
        void setBondOrderInSys2(std::vector<AtomDict> & tAtoms,
                                std::vector<BondDict> & tBonds,
                                std::vector<RingDict> & tRings);
        
        void modBondOrderViaAnnEXOneConn(std::vector<AtomDict> & tAtoms,
                                         std::vector<BondDict> & tBonds);
        void modBondOrderViaAnnEXOneLoop(std::vector<AtomDict> & tAtoms,
                                         std::vector<BondDict> & tBonds,
                                         int                   &  tNOpr);
        void checkIsoExAtoms(std::vector<AtomDict> & tAtoms);
        
        void partitionSysToSubGraphs(std::vector<AtomDict>  & tAtoms);
        
        void checkChargeInSubGraphs(std::vector<AtomDict>  & tAtoms);
        int  sumExElecsInSubGraph(std::vector<AtomDict>  & tAtoms, 
                                  std::vector<int>       & tGraph);
        void assignChargesInSubGraph(std::vector<AtomDict>  & tAtoms, 
                                  std::vector<int>       & tGraph);
        void assignChargeOneInSubGraph(std::vector<AtomDict>  & tAtoms, 
                                       std::vector<int>       & tIdxNs,
                                       bool                   & tL);
        void setEquivAtoms(std::vector<AtomDict> & tAtoms,
                           std::vector<BondDict> & tBonds);
        void modDelocBondsByEquivAtoms(std::vector<AtomDict> & tAtoms,
                                       std::vector<BondDict> & tBonds,
                                       std::vector<int>      & tIdxs);
        
        void setHMatrix(std::vector<AtomDict>  & tAtoms, REAL ** tH,
                        std::vector<int> & tSubGraph);
        void getBondOrderFromOrb(int tNOrbs, REAL ** tEigenVect,
                                 std::vector<AtomDict>  & tAtoms,
                                 std::vector<int> & tSubGraph);
        void MOSolver(std::vector<AtomDict>  & tAtoms);
        void BondTrans(std::vector<BondDict> & tBonds);
        
        void outBoAndChList(FileName                 tFName, 
                            std::vector<AtomDict>  & tAtoms,
                            std::vector<BondDict>  & tBonds);
        
        void checkUpdate(REAL & tPreV, int & tProV);
        
        bool                                   lUpdate;
        
        std::map<int, std::map<int, REAL> >    BondOrderFromMO;
        
        std::map<ID, REAL>                     orgAlphas;
        std::map<ID, REAL>                     orgBetas;
        
        std::vector<int>                       zeroExAtomIdxs;
        std::vector<int>                       withExAtomIdxs;
        std::map<int, int>                     oddAtomIdxs;
        std::map<int, std::vector<int> >       allSubGraphs;
        
        
    private:
        
        int                                    workMode;
        
    };
    
    
    
}

#endif	/* CHEMPROPSET_H */

