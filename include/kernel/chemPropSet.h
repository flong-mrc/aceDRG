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
    extern REAL  getTotalBondOrder(std::vector<BondDict>   & tBonds,
                                   std::vector<AtomDict>   & tAtoms,
                                   int                       tIA);

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
    extern int  getOneNBAtomExContri(std::vector<AtomDict>& tAtoms,
                                     int tIdxAtm, int tIdxNB);

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

    extern void outBoAndChList(FileName tFName,
                               std::vector<AtomDict>  & tAtoms,
                               std::vector<BondDict>  & tBonds);

    extern void setAtomFormTypes(std::vector<AtomDict> & tAtoms);

    extern void getDandAPair(PeriodicTable    &              tPTab,
                             std::vector<AtomDict>::iterator tAtm,
                             std::vector<AtomDict>::iterator jAtm,
                             std::map<int, std::string >    & tHPropAtom,
                             std::map<int, std::map<int, double> >
                             & tHCandAtom,
                             double    tDist);

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
        void execute3(std::vector<AtomDict> & tAtoms,
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

    class KekulizeMol
    {

    public :

        // Default constructor
        KekulizeMol();
        // Default destructor
        ~KekulizeMol();

        void execute(std::vector<AtomDict>         & tAtoms,
                     std::vector<BondDict>         & tBonds,
                     std::vector<RingDict>         & tRings,
                     std::map<std::string, int>    & tHMap);

        void executeBC(std::vector<AtomDict>         & tAtoms,
                       std::vector<BondDict>         & tBonds,
                       std::vector<RingDict>         & tRings);


        void setAllMaps(std::vector<AtomDict>      & tAtoms,
                        std::vector<BondDict>      & tBonds,
                        std::vector<RingDict>      & tRings,
                        std::map<int, int>         & tCurVal,
                        std::map<int, int>         & tOutEMap,
                        std::map<int, double>         & tChargeMap,
                        std::map<int, std::string>    & tElemMap,
                        std::map<int, int>            & tIdAtmMap,
                        std::map<int,
                        std::map<int, int> >          & tAllAtmBondingMap);

        void setSpecStrs(std::vector<AtomDict>         & tAtoms,
                         std::vector<BondDict>         & tBonds,
                         std::vector<int>           & tDoneAtoms,
                         std::vector<int>           & tDoneFAtoms,
                         std::vector<int>           & tDoneBonds,
                         std::map<int,
                         std::map<int, int> >       & tAllAtmBondingMap);

        void setSpecStrN(std::vector<AtomDict>         & tAtoms,
                         std::vector<BondDict>         & tBonds,
                         std::vector<int>           & tDoneAtoms,
                         std::vector<int>           & tDoneFAtoms,
                         std::vector<int>           & tDoneBonds,
                         std::map<int,
                         std::map<int, int> >       & tAllAtmBondingMap,
                         std::vector<AtomDict>::iterator tAt);

        void setSpecStrCL(std::vector<AtomDict>         & tAtoms,
                         std::vector<BondDict>         & tBonds,
                         std::vector<int>           & tDoneAtoms,
                         std::vector<int>           & tDoneFAtoms,
                         std::vector<int>           & tDoneBonds,
                         std::map<int,
                         std::map<int, int> >       & tAllAtmBondingMap,
                         std::vector<AtomDict>::iterator tAt);

        void setAtomsPlans(std::vector<AtomDict>    & tAtoms,
                           std::vector<BondDict>    & tBonds,
                           std::vector<RingDict>    & tRings,
                           std::map<int, bool>       & tAtmPlanMap);

        void setRingsPlans(std::vector<AtomDict>    & tAtoms,
                           std::vector<BondDict>    & tBonds,
                           std::vector<RingDict>    & tRings,
                           std::map<int, bool>      & tRingPlanMap);

        void setValOneAtoms(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>         & tBonds,
                           std::map<int, int>            & tCurVal,
                           std::map<int, int>            & tOutEMap,
                           std::map<int, double>         & tChargeMap,
                           std::map<int, std::string>    & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneBonds);

        void setOneLinkAtoms(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>         & tBonds,
                           std::map<int, int>            & tCurVal,
                           std::map<int, int>            & tOutEMap,
                           std::map<int, double>         & tChargeMap,
                           std::map<int, std::string>    & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneBonds);

        void setAllMetalBO(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>         & tBonds,
                           std::map<int, int>            & tCurVal,
                           std::map<int, int>            & tOutEMap,
                           std::map<int, double>         & tChargeMap,
                           std::map<int, std::string>    & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneBonds);


        void PickSingBonds(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>         & tBonds,
                           std::vector<RingDict>         & tRings,
                           std::map<int, int>            & tCurVal,
                           std::map<int, int>            & tOutEMap,
                           std::map<int, double>         & tChargeMap,
                           std::map<int, std::string>    & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneBonds);

        void setIsolateRings(std::vector<AtomDict>       & tAtoms,
                             std::vector<BondDict>         & tBonds,
                             std::vector<RingDict>         & tRings,
                             std::map<int,
                             std::map<int, int> >       & tAllAtmBondingMap,
                             std::map<int, int>         & tCurVal,
                             std::vector<int>           & tDoneAtoms,
                             std::vector<int>           & tDoneFAtoms,
                             std::vector<int>           & tDoneBonds,
                             std::vector<int>           & tExcFRings);

        void setOneIsolateC6Ring(std::vector<AtomDict>       & tAtoms,
                                std::vector<BondDict>       & tBonds,
                                RingDict                    & tRing,
                                std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void setOneIsolateC5Ring(std::vector<AtomDict>       & tAtoms,
                                std::vector<BondDict>       & tBonds,
                                RingDict                    & tRing,
                                std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void setOneIsolateC3N2Ring(std::vector<AtomDict>       & tAtoms,
                                   std::vector<BondDict>       & tBonds,
                                   RingDict                    & tRing,
                                   std::map<int,
                                   std::map<int, int> >       & tAllAtmBondingMap,
                                   std::map<std::string, std::map<int, int> >
                                                              & tConnMap,
                                   std::vector<int>           & tDoneAtoms,
                                   std::vector<int>           & tDoneFAtoms,
                                   std::vector<int>           & tDoneBonds);

        void setOneBondFusedRing(std::vector<AtomDict>     & tAtoms,
                                std::vector<BondDict>      & tBonds,
                                std::vector<RingDict>      & tRings,
                                std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds,
                                std::vector<int>           & tExcRings);

        void adjustC6ChargedRing(std::vector<AtomDict>     & tAtoms,
                                std::vector<BondDict>      & tBonds,
                                std::vector<RingDict>      & tRings,
                                std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void adjustC6Charged2Ring(int                       tIdxRing,
                                std::vector<AtomDict>      & tAtoms,
                                std::vector<BondDict>      & tBonds,
                                std::vector<RingDict>      & tRings,
                                std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void checkC4N1RingCharge(
                                 std::vector<AtomDict>      & tAtoms,
                                 std::vector<BondDict>      & tBonds,
                                  std::vector<RingDict>      & tRings,
                                   std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void adjustC4N1RingCharge(int                       tIdxRing,
                                 std::vector<AtomDict>      & tAtoms,
                                 std::vector<BondDict>      & tBonds,
                                 std::vector<RingDict>      & tRings,
                                 std::map<int,
                                std::map<int, int> >       & tAllAtmBondingMap,
                                std::map<int, int>         & tCurVal,
                                std::vector<int>           & tDoneAtoms,
                                std::vector<int>           & tDoneFAtoms,
                                std::vector<int>           & tDoneBonds);

        void FurtheAssignBandC(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>             & tBonds,
                           std::vector<RingDict>             & tRings,
                           std::map<int, int>                & tCurVal,
                           std::map<int, int>                & tOutEMap,
                           std::map<int, double>             & tChargeMap,
                           std::map<int, std::string>        & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneFAtoms,
                           std::vector<int>           & tDoneBonds);

        // will replace etResValForAtom and getUnsetBondsForAtom

        void setResValAndUnsetBondsForAtom(
                             AtomDict                 & tAtom,
                             std::vector<BondDict>    & tBonds,
                             std::map<int,
                             std::map<int, int> >     & tAllAtmBondingMap,
                             std::map<int, int>       & tCurVal,
                             int                      &  tResVal,
                             int                      &  tNumUnSetBonds);

        int getResValForAtom(AtomDict                 & tAtom,
                             std::vector<BondDict>    & tBonds,
                             std::map<int,
                             std::map<int, int> >     & tAllAtmBondingMap,
                             std::map<int, int>       & tCurVal);

        int getUnsetBondsForAtom(AtomDict                 & tAtom,
                                 std::vector<BondDict>    & tBonds,
                                 std::map<int,
                                 std::map<int, int> >     & tAllAtmBondingMap);

        void modifCurVal(std::vector<AtomDict>::iterator   tAtm,
                         std::map<int, int>       & tCurVal);

        void checkRingCharge(std::vector<AtomDict>         & tAtoms,
                           std::vector<BondDict>         & tBonds,
                           std::vector<RingDict>         & tRings,
                           std::map<int, int>            & tCurVal,
                           std::map<int, int>            & tOutEMap,
                           std::map<int, double>         & tChargeMap,
                           std::map<int, std::string>    & tElemMap,
                           std::map<int,
                           std::map<int, int> >       & tAllAtmBondingMap,
                           std::vector<int>           & tDoneAtoms,
                           std::vector<int>           & tDoneBonds);


        void finalAdjustBandC(std::vector<AtomDict>      & tAtoms,
                              std::vector<BondDict>      & tBonds,
                              std::vector<RingDict>      & tRings,
                              std::map<int, int>         & tCurVal,
                              std::map<int,
                              std::map<int, int> >       & tAllAtmBondingMap,
                              std::vector<int>           & tDoneAtoms,
                              std::vector<int>           & tDoneFAtoms,
                              std::vector<int>           & tDoneBonds);


        void modAromRings(std::vector<RingDict>      & tRings);

        bool checkIfAROMBs(std::vector<BondDict> & tBonds);

        void setAromBondOrderInSys(std::vector<AtomDict>       & tAtoms,
                                   std::vector<BondDict>       & tBonds,
                                   std::vector<RingDict>       & tRings,
                                   std::map<std::string, int>  & hMap);

        void modChargeAndBOforSpecialCases(
                              std::vector<AtomDict>      & tAtoms,
                              std::vector<BondDict>      & tBonds,
                              std::vector<RingDict>      & tRings);

        void preStage1(std::vector<AtomDict> & tAtoms,
                     std::vector<BondDict> & tBonds,
                     std::vector<RingDict> & tRings,
                     std::map<std::string, bool>  & tDoneAtoms,
                     std::map<int, bool>  & doneBonds,
                     std::map<std::string, int>           & tCurVal,
                     std::map<std::string, double>        & tChargeMap,
                     std::map<std::string, std::string>   & tElemMap,
                     std::map<std::string, int >    & tIdAtmMap,
                     std::map<std::string,
                     std::map<std::string, int> >   & tAllAtmBondingMap,
                     std::map<std::string,
                     std::vector<std::string> >     & tAromAtmMap,
                     std::map<std::string, int>         & tHMap);

        void setAllMaps(std::vector<AtomDict>                 & tAtoms,
                        std::vector<BondDict>                 & tBonds,
                        std::map<std::string, int>            & tCurVal,
                        std::map<std::string, double>         & tChargeMap,
                        std::map<std::string, std::string>    & tElemMap,
                        std::map<std::string, int >           & tIdAtmMap,
                        std::map<std::string,
                        std::map<std::string, int> >       & tAllAtmBondingMap,
                        std::map<std::string,
                        std::vector<std::string> >         & tAromAtmMap,
                        std::map<std::string, int>         & tHMap);

        void preStage2(std::vector<AtomDict> & tAtoms,
                     std::vector<BondDict> & tBonds,
                     std::vector<RingDict> & tRings,
                     std::map<std::string, bool>  & tDoneAtoms,
                     std::map<int, bool>  & doneBonds,
                     std::map<std::string, int>           & tCurVal,
                     std::map<std::string, double>        & tChargeMap,
                     std::map<std::string, std::string>   & tElemMap,
                     std::map<std::string, int >    & tIdAtmMap,
                     std::map<std::string,
                     std::map<std::string, int> >   & tAllAtmBondingMap,
                     std::map<std::string,
                     std::vector<std::string> >     & tAromAtmMap,
                     std::map<std::string, int>         & tHMap);

        void setRingBonds(std::vector<AtomDict> & tAtoms,
                          std::vector<BondDict> & tBonds,
                          std::vector<RingDict> & tRings,
                          std::vector<int>      & tAromRingIdxs);

        void setAtomPiInARing(std::vector<AtomDict> & tAtoms,
                              std::vector<BondDict> & tBonds,
                              RingDict              & tRing,
                              int                     tRIdx,
                              std::map<std::string,
                              std::map<std::string, int> >
                                                    & tAllAtmBondingMap,
                              std::map<int, int>    & tPiInOneAromRing );

        int  getPiInAAtom(std::vector<AtomDict>::iterator tAt,
                          std::vector<BondDict> & tBonds,
                          RingDict              & tRing,
                          std::map<std::string,
                          std::map<std::string, int> > & tAllAtmBondingMap);

        void addHInAAromRing(std::vector<AtomDict> & tAtoms,
                              std::vector<BondDict> & tBonds,
                              RingDict              & tRing,
                              int                     tRIdx,
                              std::map<std::string,
                              std::map<std::string, int> >
                                                    & tAllAtmBondingMap,
                              int                                tNumPi,
                              std::map<std::string, int>         & tHMap);

        void setAromBonds1(std::vector<AtomDict> & tAtoms,
                          std::vector<BondDict> & tBonds,
                          std::map<std::string, int>           & tCurVal,
                          std::map<std::string, double>        & tChargeMap,
                          std::map<std::string,
                          std::map<std::string, int> >   & tAllAtmBondingMap,
                          std::map<std::string, int>         & tHMap,
                          bool                  & tCh);

        void setAromBonds2(std::vector<AtomDict> & tAtoms,
                          std::vector<BondDict> & tBonds,
                          std::map<std::string, int>           & tCurVal,
                          std::map<std::string, double>        & tChargeMap,
                          std::map<std::string,
                          std::map<std::string, int> >   & tAllAtmBondingMap,
                          std::map<std::string, int>         & tHMap,
                          bool                  & tCh);

        void setAromBonds3(std::vector<AtomDict> & tAtoms,
                          std::vector<BondDict> & tBonds,
                          std::map<std::string, int>           & tCurValMap,
                          std::map<std::string, double>        & tChargeMap,
                          std::map<std::string,
                          std::map<std::string, int> >   & tAllAtmBondingMap,
                          std::map<std::string, int>         & tHMap);


        void outBondsAndHAtms(std::vector<BondDict>              & tBonds,
                              std::map<std::string, int>          & tHMap,
                              std::string                     tUserOutRoot);



        void setAllAtomEXcessElectrons(std::vector<AtomDict> & tAtoms);
        void setAllAtomCurrentValance(std::vector<AtomDict> & tAtoms,
                                      std::map<std::string, int> & tCurVal);
        int  getOneNBAtomExContri(std::vector<AtomDict> & tAtoms,
                                  int tIdxAtm, int iIdxNB);
        void initiaExElecs(std::vector<AtomDict> & tAtoms);

        void setInitBondOrdersViaExtraElecs(std::vector<AtomDict> & tAtoms,
                                            std::vector<BondDict> & tBonds);

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
                                std::vector<BondDict> & tBonds,
                                std::vector<RingDict> & tRings);

        void modBondOrderViaAnnEXOneConn(std::vector<AtomDict> & tAtoms,
                                         std::vector<BondDict> & tBonds);
        void modBondOrderViaAnnEXOneLoop(std::vector<AtomDict> & tAtoms,
                                         std::vector<BondDict> & tBonds,
                                         int                   &  tNOpr);
        void checkIsoExAtoms(std::vector<AtomDict> & tAtoms);

        void checkUpdate(REAL & tPreV, int & tProV);

        void partitionSysToSubGraphs(std::vector<AtomDict>  & tAtoms);

        void checkChargeInSubGraphs(std::vector<AtomDict>  & tAtoms);

        int  sumExElecsInSubGraph(std::vector<AtomDict>  & tAtoms,
                                  std::vector<int>       & tGraph);

        void assignChargesInSubGraph(std::vector<AtomDict>  & tAtoms,
                                  std::vector<int>       & tGraph);

        void assignChargeOneInSubGraph(std::vector<AtomDict>  & tAtoms,
                                       std::vector<int>       & tIdxNs,
                                       bool                   & tL);

        void fromSubGraphsToRings(std::vector<AtomDict>  & tAtoms,
                                  std::vector<RingDict> & tRings,
                                  std::vector<int> & tUndecidedRingIdxs,
                                  std::map<int, std::vector<int> > & tNonRingAtomIdxs);

        void checkOneRingInSubgraphs
                        (int  tOneRingIdx,
                        std::vector<int> & tOneRingAtomIdxs,
                        std::vector<int> & tUndecidedRingIdxs,
                        std::map<int, std::vector<int> > & tNonRingAtomIdxs);

        void kekulizeRings(std::vector<AtomDict> & tAtoms,
                           std::vector<BondDict> & tBonds,
                           std::vector<RingDict> & tRings,
                           std::vector<int> & tUndecidedRingIdxs);
        void kekulizeOneRing(std::vector<AtomDict> & tAtoms,
                                std::vector<BondDict> & tBonds,
                                RingDict & tRing,
                                PeriodicTable & tTab);
        void modifyBondOrderInOneRing(std::vector<BondDict> & tBonds,
                                      std::vector<AtomDict>  & tAtoms,
                                      int  tIdxB1, int tIdxB2,
                                      int tAtCen, int tAt1, int tAt2,
                                      PeriodicTable & tTab);

        REAL getFixedBondOrder(std::vector<BondDict>   & tBonds,
                               std::vector<AtomDict>   & tAtoms,
                               int                       tAtmIdx);

        void setEquivAtoms(std::vector<AtomDict> & tAtoms,
                           std::vector<BondDict> & tBonds);

        void modDelocBondsByEquivAtoms(std::vector<AtomDict> & tAtoms,
                                       std::vector<BondDict> & tBonds,
                                       std::vector<int>      & tIdxs);



        void outBandC(FileName tFName,
                      std::vector<AtomDict> & tAtoms,
                      std::vector<BondDict> & tBonds);


        bool                                   lUpdate;

        std::vector<int>                       zeroExAtomIdxs;
        std::vector<int>                       withExAtomIdxs;
        std::map<int, int>                     oddAtomIdxs;
        std::map<int, std::vector<int> >       allSubGraphs;

        std::vector<std::string>               doneAtoms;
        std::vector<int>                       doneBonds;

    private:

    };

}

#endif	/* CHEMPROPSET_H */
