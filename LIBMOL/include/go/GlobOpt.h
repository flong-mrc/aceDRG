/* 
 * File:   GlobOpt.h
 * Author: flong
 *
 * Created on February 12, 2013, 3:44 PM
 */

#ifndef GLOBOPT_H
#define	GLOBOPT_H

#ifndef KERNEL_H
#include "../kernel/kernel.h"
#endif

#ifndef ATOM_H
#include "../kernel/atom.h"
#endif

#ifndef BOND_H
#include "../kernel/bond.h"
#endif

#ifndef ANGLE_H
#include "../kernel/angle.h"
#endif

#ifndef TORSION_H
#include "../kernel/torsion.h"
#endif

#ifndef RESIDUE_H
#include "../kernel/residue.h"
#endif

#ifndef RING_H
#include "../kernel/ring.h"
#endif

#ifndef PLANE_H
#include "../kernel/plane.h"
#endif

#ifndef CHAIN_H
#include "../kernel/chain.h"
#endif

#ifndef UTILITY_H
#include "../kernel/utility.h"
#endif

#ifndef NEIGHBLIST_H
#include "../kernel/neighbList.h"
#endif

#ifndef LINALG_H
#include "LinAlg.h"
#endif

#ifndef ALLSYSTEM_H
#include "../kernel/AllSystem.h"
#endif

#ifndef TRANSCOORD_H
#include "../kernel/TransCoord.h"
#endif

namespace GO
{
    class OptimSet
    {
    public :
        
        // Default constructor
        OptimSet();
        // Default destructor
        ~OptimSet();
        
        LIBMOL::REAL                       objValue;
        std::vector<LIBMOL::AtomDict>      atoms;
        
    };
    
    class FindGlobMin
    {
    public:
        
        // Default constructor
        FindGlobMin();
        
        // Constructor using a AllSystem
        FindGlobMin(const LIBMOL::AllSystem & aSystem,
                    std::string    tFileNameeRoot,
                    std::string    tMonoRoot);
        
        FindGlobMin(const LIBMOL::AllSystem & aSystem);
        
        FindGlobMin(const LIBMOL::AllSystem & aSystem,
                    int     nOpt);
        
        // Default destructor 
        ~FindGlobMin();
        
        // Common routines 
        void Driver();                    
        void Driver(bool useCoords, int nOpt); 
        
       
        void SetDefinedInitPositions();
        void SetRanDomInitPositions();
        
        void UpdatePositionsWithRandomShifts();
        
        void PreIdealization();
        void ProIdealization();
        void singleCompsIdealization();
         
        /* !!!!!!!!! Different Global min schemes
         * 
         * 1 Tunneling and local minimum scheme
         * 
         * 2 Fast SA (Evolution + Flash-freezing)
         * 2.1 Using MC + CG
         * 2.2 Using MD + CG
         * 
         * 
         */ 
        // 1 Tunneling and local minimum scheme
        

        void SetTAndRRatio(bool & i_err);

        void SetBalanceValue(int      & idx_error,
                             LIBMOL::REAL & bal_value);

        void SetBalanceValueAndDeriv(int          &  idx_error,
                                     LIBMOL::REAL &  bal_value,
                                     LIBMOL::REAL ** bal_deriv);

        void SetBalanceValueAndDeriv_EXP(int           &  idx_error,
                                         LIBMOL::REAL      coeffi_ad,
                                         LIBMOL::REAL  &  bal_value,
                                         LIBMOL::REAL ** bal_deriv,
                                         LIBMOL::REAL ** pos_m);
                                   

        void SetBalanceValueAndDeriv_LOG(int          &  idx_error,
                                         LIBMOL::REAL &  o_diff,
                                         LIBMOL::REAL &  bal_value,
                                         LIBMOL::REAL ** bal_deriv);

        void RootFinder_Broyden();

        void RootFinder_Newton(); 
 
        void RootFinder_Newton_GV(); 

        void Tunnelling_Newton_EXP();

        void Tunnelling_Newton_LOG();

        void TunnellingAndMinDriver();
        
        void TunnellingAndMinDriver(int nOpt);
        
        
        // A hybrid SA/CG/scheme

        void HybridSA_CGMinDriver();
        
        void Fast_SA();

        void FastAnnealing(int  num_steps);

        void MD_Init(int          idx_way,
                     LIBMOL::REAL   & obj_d,
                     LIBMOL::REAL   & b_value,
                     LIBMOL::REAL  ** b_deriv);

        void InitCartDisplaces();

        void InitCartVeloAndAccels(int          idx_way,
                                   LIBMOL::REAL   & obj_d,
                                   LIBMOL::REAL   & b_value,
                                   LIBMOL::REAL  ** b_deriv);

  

        void HYBRID_MC(LIBMOL::REAL   & obj_d, 
                       LIBMOL::REAL   & b_value,
                       LIBMOL::REAL  ** b_deriv,  
                       int        & steps_cur);

        void AcceptanceRules(bool     &  l_accept,
                             LIBMOL::REAL    e_diff,
                             LIBMOL::REAL    h_diff);
                     
        LIBMOL::REAL GetHamilton (LIBMOL::REAL e_pot);

  
        // Detailed calculations of object functions and their derivatives  
        // Methods related to transfer from the Cartesian description
        // to the angular(internal) description
 
        void CoordCartToAng();
        
        void SetDThetaDx(LIBMOL::REAL ** dtheta_dx);
        
        void GetDvDq(int idx_atom, int idx_torsion);
        void TorsAngDri(int         idx_torsion,
                        LIBMOL::REAL *  dPsidXi,
                        LIBMOL::REAL ** dPsidXidXj);
        
        void VeloCartToAng ();


        // The objective function. 

        LIBMOL::REAL GetObjValue();
        LIBMOL::REAL GetObjValue(LIBMOL::REAL * vars);

        LIBMOL::REAL GetObjOneBond(std::vector<LIBMOL::BondDict>::iterator  tBo);
        LIBMOL::REAL GetObjOneAngle(std::vector<LIBMOL::AngleDict>::iterator tAn);
        LIBMOL::REAL GetObjOneTorsion(std::vector<LIBMOL::TorsionDict>::iterator tTo);
        LIBMOL::REAL GetObjOnePlane(std::vector<LIBMOL::PlaneDict>::iterator tPl);
        LIBMOL::REAL GetObjOneChiral(std::vector<LIBMOL::ChiralDict>::iterator tCh);
        LIBMOL::REAL GetObjOneAtomVdw(std::vector<LIBMOL::AtomDict>::iterator tAt);

        void     SetObjCart();
        void     SetObjTors();
        void     SetObjandForceCart();

        // The force (only in Cartisan space) and 
        // accelerate (in both Cartesian  and internal space)

        void SetInForces();

        void GetForBondOrder(int i_bond);  
                                              
        void GetForBondAngs(int i_bondAngs);          

        void GetForTorsionAngs(int i_atom, 
                               LIBMOL::REAL ** dPsi_dxdx);

  
        //-------- vdw and hydrogen bonding ----

        //  values of objective function 

        LIBMOL::REAL GetObjVWRepuls(int i_atom,
                                    int j_atom);

        void CheckHYD(int        i_atom,
                      int        j_atom,
                      bool     & l_hy);

        LIBMOL::REAL GetObjHdrBonding_qu(LIBMOL::REAL   d_contact,
                                         LIBMOL::REAL   r_distance,
                                         LIBMOL::REAL & obj_hyd);

        LIBMOL::REAL GetObjHdrBonding_inv(LIBMOL::REAL   d_contact,
                                          LIBMOL::REAL   r_distance,
                                          LIBMOL::REAL & obj_hyd);
                            
        void CheckAndCalcHYDObj(int        i_atom1,
                                int        i_atom2,
                                LIBMOL::REAL   d_contact,
                                LIBMOL::REAL   r_distance,
                                LIBMOL::REAL & obj_hyd,
                                bool     & l_hyd_vdw);

        void CalcObjHYD_qu(LIBMOL::REAL   d_contact,
                           LIBMOL::REAL   d_cutoff,
                           LIBMOL::REAL   r_distance,
                           LIBMOL::REAL & obj_v);    

        void CalcObjHYD_inv(LIBMOL::REAL   d_contact,
                            LIBMOL::REAL   r_distance,
                            LIBMOL::REAL & obj_v);

        void CalcObjVDW_qu(LIBMOL::REAL   d_contact,
                           LIBMOL::REAL   d_cutoff,
                           LIBMOL::REAL   r_distance,
                           LIBMOL::REAL & obj_v);  

        // force
        void GetForVWRepuls(int i_atom);

        void CalcHBForce_qu(int            i_atom1,
                            int            i_nbatom,
                            LIBMOL::REAL   d_contact,
                            LIBMOL::REAL   d_cutoff,
                            LIBMOL::REAL   r_distance,
                            LIBMOL::REAL * r_comp); 

        void CalcHBForce_inv(int            i_atom1,
                             int            i_nbatom,
                             LIBMOL::REAL   d_contact,
                             LIBMOL::REAL   r_distance,
                             LIBMOL::REAL * r_comp); 

        void CalcVDWForce_qu(int        i_atom1,
                            int        i_nbatom,
                            LIBMOL::REAL   d_contact, 
                            LIBMOL::REAL   d_cutoff,
                            LIBMOL::REAL   r_distance,
                            LIBMOL::REAL * r_comp);

 
        void CalcForce(int            i_atom1,
                      LIBMOL::REAL   d_contact,
                      LIBMOL::REAL   d_cutoff,
                      LIBMOL::REAL   r_distance,
                      LIBMOL::REAL * r_comp,
                      bool           l_hy_vdw,
                      bool         & l_next); 


  
        // mixed calculation

        LIBMOL::REAL  GetObjAndtForce();

        LIBMOL::REAL  GetBondLengthContrib( int   idx_bond,
                                            int & idx_error);

        LIBMOL::REAL  GetBondAngContrib(int    idx__ang, 
                                        int  & idx_error);

        LIBMOL::REAL  GetPlaneContrib(int      idx__pl, 
                                      int    & idx_error); 

        LIBMOL::REAL  GetChiralityContrib(int    idx__chira, 
                                          int  & idx_error);
  
        LIBMOL::REAL  GetVWRandHBContrib(int     idx_atom, 
                                         int   & idx_error); 

        // derivatives

        void CheckAndCalcHYDFD(int            i_atom1,
                               int            i_atom2,
                               LIBMOL::REAL   d_contact, 
                               LIBMOL::REAL   r_distance,
                               LIBMOL::REAL & d2f_dr2,
                               bool         & l_flag);


        void CalcHYDFD(LIBMOL::REAL d_contact, 
                       LIBMOL::REAL r_distance,
                       LIBMOL::REAL & d2f_dr2, 
                       bool     & l_hy_vdw);

        // normal matrix

        void VDWToNormalMatrixSparse(int           atom1, 
                                     int           atom2,
                                     int           vdw,
                                   LIBMOL::REAL **   f_ij);
 
       void HYDToNormalMatrixSparse(int               i_atom1,
                                    int               i_atom2,
                                    int               i_hyd,
                                    bool          &   l_hyd_calc, 
                                    LIBMOL::REAL  **  fij);

       void SetFirstDerivVDW(int            idx_atom,
                             int            idx_neighb,
                             LIBMOL::REAL * dfdx_vdw);


       void SetFirstDerivHYD(int            idx_atom,
                             int            idx_neighb,
                             bool         & l_hyd_calc,
                             LIBMOL::REAL * dfdx_hyd);


       void SetVDWandHBConst();

       void SetVDWContact (int             i_atom1, 
                           int             i_atom2, 
                           LIBMOL::REAL  & d_cutoff);

       void SetHYDContact(int            i_atom1, 
                          int            i_atom2, 
                          LIBMOL::REAL & d_contact,
                          LIBMOL::REAL & d_cutoff);

       void SetHBCoeffs(LIBMOL::REAL   d_0,
                        LIBMOL::REAL & A_ij,
                        LIBMOL::REAL & B_ij);

       void  CheckHYDFlag(int      i_atom1,
                          int      i_atom2,
                          int      i_atom3,
                          bool   & l_flag);
  
       void  CheckHYDFlag(int      i_atom1,
                          int      i_atom2,
                          bool   & l_flag);
 
       void CheckHYD_DA_Flag(int    i_atom1,
                             int    i_atom2,
                             bool & l_flag);

       void CheckHYD_HA_Flag(int    i_atom1,
                             int    i_atom2,
                             bool & l_flag);


       // -----------------------------------------------------  
  
  
       void GetForPlanes(int i_plans);
  
       void GetForChiras(int i_chiras);                      

       void AccelCartesian();
  
       void AccelCartesian(LIBMOL::REAL ** potential_deriv);

       void AccelTorsAng();

       // Methods related to numeric integrations  
  
       // a: Numerov algorithm
  
       void IntegEuler(LIBMOL::REAL & vaule_pre, LIBMOL::REAL & value,
                       LIBMOL::REAL & velo_pre,  LIBMOL::REAL & velo,      
                       LIBMOL::REAL & accel_pre, LIBMOL::REAL & accel,     
                       LIBMOL::REAL d_time);

       void IntegPred(LIBMOL::REAL & value_t, LIBMOL::REAL & value_pre,   
                      LIBMOL::REAL & value,   LIBMOL::REAL & accel,   
                      LIBMOL::REAL d_time);

       void IntegCorr_1(LIBMOL::REAL & value_t, LIBMOL::REAL & value_pre,
                        LIBMOL::REAL & value,   LIBMOL::REAL & accel_pre,
                        LIBMOL::REAL & accel_t, LIBMOL::REAL & accel,
                        LIBMOL::REAL d_time);  

       void IntegCorr_2(LIBMOL::REAL & value_pre,   LIBMOL::REAL & value,   
                        LIBMOL::REAL & d_value_pre, LIBMOL::REAL & d_value_t,  
                        LIBMOL::REAL & d_value,     LIBMOL::REAL d_time);                           
 
       void NumerovInt();

       void NumerovInt(int           i_iterations,
                       LIBMOL::REAL  &  o_diff,
                       LIBMOL::REAL  &  b_diff,
                       LIBMOL::REAL **  b_deriv);
   
       // b: Velocity Verlet algorithm

       void VeloVerletInt();

       // c: Conventional Verlet algorithm

       void VerletInt();

       //  A template methods that pick up the backbone torsion angles
       //  out of all existing torsion angles in the tree 
       // (Change will be made when the data structures are set

       void UpdatTors_t();
  
       // Metheds that control MD process

       // Intialization needed by MD
  
       LIBMOL::REAL GetRand();


       void InitCartVelo(int idx_way);
       void InitTorsVelo(int idx_way);
  

       void MD_Init();
       void Dynamics();
       void Evolution();                 // put here temporarily
       

       // Methods linked to the kinetics (in MD)

       void CheckTotalE();
       void TempScale(); 
       void Quenching_V();

       // Method related to debugging

       void StrucTransTest();
       void VeloTransTest();
  
       // Methods related to static minimizations(flash Freezing)
       // 1. process control
  
       void FlashingFreezing();             
       void SwitchSMIN();
       void CheckShift(LIBMOL::REAL *d_x, 
                       int          n_vars);  
  
       // 2. The first and second derivatives in the Cartesian space
 
       void SetFirstDrivs();
       void SetNormalMatrix();
       void SetNormalMatrixSparse(int & num_block);
       void SetNormalMatrixSparseTors(int & num_block);

       void BondToNormalMatrix();
       void BondToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                               LIBMOL::REAL   * col_ind,
                               LIBMOL::REAL   * row_blk);
       void BondToNormalMatrixSparse(int               atom1, 
                                     int               atom2,
                                     int               bond,
                                     LIBMOL::REAL **   f_ij);


       void BondAngToNormalMatrix();
       void BondAngToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                                  LIBMOL::REAL   * col_ind,
                                  LIBMOL::REAL   * row_blk);
       void BondAngToNormalMatrixSparse(int               atom1, 
                                        int               atom2,
                                        int               ang,
                                        LIBMOL::REAL  *   dth_dx,
                                        LIBMOL::REAL **   f_ij);

    
       void TorsAngToNormalMatrix();
       void TorsAngToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                                  LIBMOL::REAL   * col_ind,
                                  LIBMOL::REAL   * row_blk);
       
       void TorsAngToNormalMatrixSparse(int           atom1, 
                                        int           atom2,
                                        int           tors,
                                        LIBMOL::REAL **   f_ij);
  
  
       void PlaneToNormalMatrix();
       
       void PlaneToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                                LIBMOL::REAL   * col_ind,
                                LIBMOL::REAL   * row_blk);
       
       void PlaneToNormalMatrixSparse(int               atom1, 
                                      int               atom2,
                                      int               plan,
                                      LIBMOL::REAL **   df_dx_p,
                                      LIBMOL::REAL **   f_ij);


       void ChiraToNormalMatrix();
       void ChiraToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                                LIBMOL::REAL   * col_ind,
                                LIBMOL::REAL   * row_blk);

       void ChiraToNormalMatrixSparse(int               atom1, 
                                      int               atom2,
                                      int               chir,
                                      LIBMOL::REAL  *   df_dx_c,
                                      LIBMOL::REAL **   f_ij);

       void VDWToNormalMatrix();
       void VDWToNormalMatrix(LIBMOL::REAL *** secdrivcart,
                              LIBMOL::REAL   * col_ind,
                              LIBMOL::REAL   * row_blk);

       void SetFirstDerivBondOrder(int            idx_bonds, 
                                   LIBMOL::REAL * dfdx_b);
       void SetFirstDerivPlane(int              idx_planes,
                               LIBMOL::REAL  ** df_dx_p);
       void SetFirstDerivChira(int            idx_chira,
                               LIBMOL::REAL  *df_dx_c);
       void SetChiraAndFirstDeriv(int           idx_chira,
                                  LIBMOL::REAL  volume,
                                  LIBMOL::REAL  *df_dx_c);

       // 3.The first and second derivatives in the Torsion angle space

       void SetTorsNormalUnit();
       void SetTorsAcu();

       void SetDobjDtorsFirst();
       LIBMOL::REAL GetDobjDphi_i(int idx_tors); 
       void SetDobjDtorsSecond();
       void SetTorsDerivs();
       void SetDerivTorsSpace();
       
       
       // Find a minimum along a line
       // Put liner minimizer here because it involves many objective function 
       // calculations 
       void LinMin(LIBMOL::REAL    * delta_x);
       void Bracket(int              num_vars,
                    LIBMOL::REAL   * delta_x,
                    LIBMOL::REAL   * x_a,
                    LIBMOL::REAL   * x_b,
                    LIBMOL::REAL   * x_c,
                    LIBMOL::REAL   & obj_a,
                    LIBMOL::REAL   & obj_b,         
                    LIBMOL::REAL   & obj_c);
  
       void Bracket(LIBMOL::REAL   * x_a,
                    LIBMOL::REAL   * x_b);
 
       void GSection(LIBMOL::REAL  * x_a,
                     LIBMOL::REAL  * x_b,
                     LIBMOL::REAL  * x_c,
                     LIBMOL::REAL  & obj_a,     
                     LIBMOL::REAL  & obj_b,
                     LIBMOL::REAL  & obj_c);

       void checkConvergence(int            idx_it,
                             LIBMOL::REAL * tol_t); 
       
       
       // 3. distance matrix 
       void DistMatrixMin();
       
       void LocalMin();
       void AddToOptim();
       
       std::vector<LIBMOL::AtomDict>                            allAtoms;
       std::vector<int>                                         allHAtomIdx; //repeated ones
       std::vector<int>                                         allHydroAtoms;
       std::vector<LIBMOL::BondDict>                            allBonds;
       std::vector<LIBMOL::AngleDict>                           allAngles;
        
       std::vector<LIBMOL::TorsionDict>                         allTorsions;
       std::vector<LIBMOL::ChiralDict>                          allChirals;
       std::vector<LIBMOL::PlaneDict>                           allPlanes;
       std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> >     allRings;
       std::vector<LIBMOL::RingDict>                            allRingsV;
        
       std::vector<OptimSet>                                    allOptimSets;
        
       int                                                      gWorkMode;
       int                                                      lWorkMode;
       int                                                      workSpace;
       int                                                      lComp;
       LIBMOL::REAL                                             curObjValue;
       LIBMOL::REAL                                             preObjValue;
       int                                                      numOptimSets;
       
    private:
        
        bool                                                    hasIniCoords;
        std::string                                             itsFileNameRoot;
        std::string                                             itsMonoRoot;
        
    };
    
    
    class GetObjValue 
    {
        
    };
    
    extern bool sortOptSetByValues(const OptimSet & tSetA, const OptimSet & tSetB);

    
}

#endif	/* GLOBOPT_H */

