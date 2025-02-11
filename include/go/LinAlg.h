/* 
 * File:   LinAlg.h
 * Author: flong
 *
 * Created on February 12, 2013, 3:55 PM
 */

#ifndef LINALG_H
#define	LINALG_H

#ifndef KERNEL_H
#include "../kernel/kernel.h"
#endif

#ifndef __JAMA_EIGN_H_
#include "../math/jama_eig.h"
#endif

#ifndef __JAMA_CHOLESKY_H_
#include "../math/jama_cholesky.h"
#endif

#ifndef __JAMA_QR_H_
#include "jama_qr.h"
#endif

#ifndef __JAMA_SVD_H_
#include "../math/jama_svd.h"
#endif

#ifndef GETDERIVS_H
#include "../forcefield/getDerivs.h"
#endif

#ifndef GETOBJVALUES_H
#include "../forcefield/getObjValues.h"
#endif



namespace GO
{
    class FindLocalMin
    {
    public:
        
        // Default constructor
        FindLocalMin();
        
        // Default destructor
        ~FindLocalMin();
        
        
        // overall controller
        void Driver(std::vector<LIBMOL::AtomDict>    & tAtoms, 
                    std::vector<LIBMOL::BondDict>    & tBonds, 
                    std::vector<LIBMOL::AngleDict>   & tAngs, 
                    std::vector<LIBMOL::TorsionDict> & tTorsions, 
                    std::vector<LIBMOL::RingDict>    & tRings, 
                    std::vector<LIBMOL::PlaneDict>   & tPlas, 
                    std::vector<LIBMOL::ChiralDict>  & tChs);
        
        void Driver(std::vector<LIBMOL::AtomDict>    & tAtoms, 
                    std::vector<LIBMOL::BondDict>    & tBonds, 
                    std::vector<LIBMOL::AngleDict>   & tAngs, 
                    std::vector<LIBMOL::TorsionDict> & tTorsions, 
                    std::vector<LIBMOL::RingDict>    & tRings, 
                    std::vector<LIBMOL::PlaneDict>   & tPlas, 
                    std::vector<LIBMOL::ChiralDict>  & tChs,
                    std::vector<LIBMOL::AtomDict>    & tAllAtoms);
        
        // 1 Conjugate direction method
        // 1.1 Conjugate gradient method to solve a linear equation
        
        int  CGGEqs(int              n_dim,
                    LIBMOL::REAL  ** sec_deriv, 
                    LIBMOL::REAL  *  firs_deriv,
                    LIBMOL::REAL  *  delta_x);
        
        
        int CGGEqs(int               n_size,
                    int               n_dim,
                    int               n_block,
                    LIBMOL::REAL    * first_deriv,
                    LIBMOL::REAL    * secDerivCartSparse_row,
                    LIBMOL::REAL    * secDerivCartSparse_col,
                    LIBMOL::REAL  *** secDerivCartSparse,
                    LIBMOL::REAL    * delta_x);

        void Stab_plusAlpha (int             n_size,
                             LIBMOL::REAL    alp,
                             LIBMOL::REAL ** A_matrix);
        // sparse matrix format
        void Stab_plusAlpha (LIBMOL::REAL       alp,
                             int                n_dim,
                             int                n_block,
                             LIBMOL::REAL *     secDerivCartSparse_row,
                             LIBMOL::REAL ***   secDerivCartSparse);    
 
        void Cancel_plusAlpha(int            n_size,
                              LIBMOL::REAL   alph,
                              LIBMOL::REAL * vect);
        
        
        //1.2 Find a minimum along a line
        
        void LinMin(LIBMOL::REAL    * delta_x,
                    std::vector<LIBMOL::AtomDict>    & tAtoms, 
                    std::vector<LIBMOL::BondDict>    & tBonds, 
                    std::vector<LIBMOL::AngleDict>   & tAngs, 
                    std::vector<LIBMOL::TorsionDict> & tTors, 
                    std::vector<LIBMOL::RingDict>    & tRings, 
                    std::vector<LIBMOL::PlaneDict>   & tPlas, 
                    std::vector<LIBMOL::ChiralDict>  & tChs);

        void Bracket(int              num_vars,
                     LIBMOL::REAL   * delta_x,
                     LIBMOL::REAL   * x_a,
                     LIBMOL::REAL   * x_b,
                     LIBMOL::REAL   * x_c,
                     LIBMOL::REAL   & obj_a,
                     LIBMOL::REAL   & obj_b,         
                     LIBMOL::REAL   & obj_c,
                     std::vector<LIBMOL::AtomDict>   & tAtoms, 
                     std::vector<LIBMOL::BondDict>   & tBonds, 
                     std::vector<LIBMOL::AngleDict>  & tAngs, 
                     std::vector<LIBMOL::TorsionDict>& tTors, 
                     std::vector<LIBMOL::RingDict>   & tRings, 
                     std::vector<LIBMOL::PlaneDict>  & tPlas, 
                     std::vector<LIBMOL::ChiralDict> & tChs);
  
        void Bracket(LIBMOL::REAL   * x_a,
                     LIBMOL::REAL   * x_b);
 
        void GSection(LIBMOL::REAL  * x_a,
                      LIBMOL::REAL  * x_b,
                      LIBMOL::REAL  * x_c,
                      LIBMOL::REAL  & obj_a,     
                      LIBMOL::REAL  & obj_b,
                      LIBMOL::REAL  & obj_c,
                      std::vector<LIBMOL::AtomDict>    & tAtoms, 
                      std::vector<LIBMOL::BondDict>    & tBonds, 
                      std::vector<LIBMOL::AngleDict>   & tAngs, 
                      std::vector<LIBMOL::TorsionDict> & tTors, 
                      std::vector<LIBMOL::RingDict>    & tRings, 
                      std::vector<LIBMOL::PlaneDict>   & tPlas, 
                      std::vector<LIBMOL::ChiralDict>  & tChs);

        void checkConvergence(int            idx_it,
                              LIBMOL::REAL   tol_t,
                              LIBMOL::REAL & tol_dx);
       
        void DelAconjSingularity (LIBMOL::REAL *  vect_1,
                                  LIBMOL::REAL ** matr,
                                  LIBMOL::REAL *  vect_2);
              
        
        // 2 Eigenvalue equation method
        // For full matrix 
        void EigenEqs(LIBMOL::REAL ** sec_deriv, 
                      LIBMOL::REAL  * firs_deriv,
                      LIBMOL::REAL  * delta_x);
        
        // 3 used in different methods
        
        int CheckShift(LIBMOL::REAL *  d_x, 
                       int             n_vars);
        
        void UpdateAtomCoordsCart(std::vector<LIBMOL::AtomDict> & tAtoms,
                                  LIBMOL::REAL    *  tDx);
        void UpdateAtomCoordsTors(std::vector<LIBMOL::AtomDict> & tAtoms,
                                  std::vector<LIBMOL::TorsionDict> & tTors,
                                  LIBMOL::REAL  *   tDx);
        void TransAtomCoordsGenToCart(std::vector<LIBMOL::AtomDict> &  tAtoms, 
                                      LIBMOL::REAL  *  tX);
        void TransAtomCoordsCartToGen(std::vector<LIBMOL::AtomDict> &  tAtoms, 
                                            LIBMOL::REAL  *  tX);
        
        void TransAtomCoordsGenToTors(std::vector<LIBMOL::TorsionDict> & tTors,
                                      LIBMOL::REAL             * tX);
        void TransAtomCoordsTorsToGen(std::vector<LIBMOL::TorsionDict> & tTors,
                                      LIBMOL::REAL             * tX);
       
        LIBMOL::REAL  GetObjectValue(LIBMOL::REAL  * tX,
                                     std::vector<LIBMOL::AtomDict>& tAtoms, 
                                     std::vector<LIBMOL::BondDict>& tBonds, 
                                     std::vector<LIBMOL::AngleDict>& tAngs, 
                                     std::vector<LIBMOL::TorsionDict>& tTors, 
                                     std::vector<LIBMOL::RingDict> & tRings, 
                                     std::vector<LIBMOL::PlaneDict>& tPlas, 
                                     std::vector<LIBMOL::ChiralDict>& tChs);
        void UpdateAtomCoords(std::vector<LIBMOL::AtomDict> tAtoms,
                              LIBMOL::REAL * tX);
       
        
        int               workMode;       // different algorithms
        int               workSpace;      // cart and torsion 
        int               lComp;
        
        int               maxIter;
        LIBMOL::REAL      maxSumDX;
        LIBMOL::REAL      maxTol;
        
        LIBMOL::REAL      curObjValue;
        
        
    };
    
   
    class LBFGS
    {
    public :
        
        // Default constructor 
        LBFGS();
        // user input parameters
        LBFGS(int             tVars,
              int             tKeep,
              LIBMOL::REAL    txTol,
              LIBMOL::REAL    trTol,
              LIBMOL::REAL    tMinStep,
              LIBMOL::REAL    tMaxStep);
        
        // Default destructor
        ~LBFGS();
        
        void driver(LIBMOL::REAL  & finObj, 
                    std::vector<LIBMOL::REAL> & finCoords,
                    std::vector<LIBMOL::REAL> & iniCoords,
                    std::vector<LIBMOL::REAL> & iniGrads);
        
        void driver(LIBMOL::REAL                     & finObj, 
                    std::vector<LIBMOL::REAL>        & iniCoords,
                    std::vector<LIBMOL::AtomDict>    & tAtoms, 
                    std::vector<LIBMOL::BondDict>    & tBonds, 
                    std::vector<LIBMOL::AngleDict>   & tAngs, 
                    std::vector<LIBMOL::TorsionDict> & tTors, 
                    std::vector<LIBMOL::RingDict>    & tRings, 
                    std::vector<LIBMOL::PlaneDict>   & tPlas, 
                    std::vector<LIBMOL::ChiralDict>  & tChs);
        
        void driver(LIBMOL::REAL                     & finObj, 
                    std::vector<LIBMOL::REAL>        & iniCoords,
                    std::vector<LIBMOL::AtomDict>    & tAtoms, 
                    std::vector<LIBMOL::BondDict>    & tBonds, 
                    std::vector<LIBMOL::AngleDict>   & tAngs, 
                    std::vector<LIBMOL::TorsionDict> & tTors, 
                    std::vector<LIBMOL::RingDict>    & tRings, 
                    std::vector<LIBMOL::PlaneDict>   & tPlas, 
                    std::vector<LIBMOL::ChiralDict>  & tChs,
                    std::vector<LIBMOL::AtomDict>    & tAllAtoms);
        
        void init(std::vector<LIBMOL::REAL> & iniGrads);
        void init(int tKeep, LIBMOL::REAL txTol, LIBMOL::REAL trTol,
                  LIBMOL::REAL tMinStep, LIBMOL::REAL tMaxStep,
                  std::vector<LIBMOL::REAL> & iniGrads);

        void setParams(int tKeep, LIBMOL::REAL txTol, LIBMOL::REAL trTol,
                       LIBMOL::REAL tMinStep, LIBMOL::REAL tMaxStep);
      
        void setSandY(std::vector<LIBMOL::REAL>   & tNewCoords,
                      std::vector<LIBMOL::REAL>   & tOldCoords,
                      std::vector<LIBMOL::REAL>   & tNewDerivs,
                      std::vector<LIBMOL::REAL>   & tOldDerivs);
        
        void setHk0();
        void setPk(std::vector<LIBMOL::REAL>   & tDerivs);
        
        void lineSearch(std::vector<LIBMOL::REAL> & finCoords,
                        std::vector<LIBMOL::REAL> & iniCoords,
                        std::vector<LIBMOL::REAL> & iniDeriv,
                        LIBMOL::REAL              &  tObj);
        void lineSearch(std::vector<LIBMOL::REAL> & finCoords,
                        std::vector<LIBMOL::REAL> & iniCoords,
                        std::vector<LIBMOL::REAL> & iniDeriv,
                        LIBMOL::REAL              &  tObj,
                        std::vector<LIBMOL::AtomDict>    & tAtoms, 
                        std::vector<LIBMOL::BondDict>    & tBonds, 
                        std::vector<LIBMOL::AngleDict>   & tAngs, 
                        std::vector<LIBMOL::TorsionDict> & tTors, 
                        std::vector<LIBMOL::RingDict>    & tRings, 
                        std::vector<LIBMOL::PlaneDict>   & tPlas, 
                        std::vector<LIBMOL::ChiralDict>  & tChs);
        
        void Bracket(std::vector<LIBMOL::REAL> & x_a, 
                     std::vector<LIBMOL::REAL> & x_b, 
                     std::vector<LIBMOL::REAL> & x_c, 
                     LIBMOL::REAL& obj_a, 
                     LIBMOL::REAL& obj_b, 
                     LIBMOL::REAL & obj_c);
        void Bracket(std::vector<LIBMOL::REAL> & x_a, 
                     std::vector<LIBMOL::REAL> & x_b, 
                     std::vector<LIBMOL::REAL> & x_c, 
                     LIBMOL::REAL& obj_a, 
                     LIBMOL::REAL& obj_b, 
                     LIBMOL::REAL & obj_c,
                     std::vector<LIBMOL::AtomDict>    & tAtoms, 
                     std::vector<LIBMOL::BondDict>    & tBonds, 
                     std::vector<LIBMOL::AngleDict>   & tAngs, 
                     std::vector<LIBMOL::TorsionDict> & tTors, 
                     std::vector<LIBMOL::RingDict>    & tRings, 
                     std::vector<LIBMOL::PlaneDict>   & tPlas, 
                     std::vector<LIBMOL::ChiralDict>  & tChs);
        
        void GSection(std::vector<LIBMOL::REAL>  & x_a, 
                     std::vector<LIBMOL::REAL>   & x_b, 
                     std::vector<LIBMOL::REAL>   & x_c, 
                     LIBMOL::REAL & obj_a, 
                     LIBMOL::REAL & obj_b, 
                     LIBMOL::REAL & obj_c,
                     LIBMOL::REAL & tObj,
                     std::vector<LIBMOL::REAL> & tFinCoords,
                     LIBMOL::REAL (LBFGS::*func)(std::vector<LIBMOL::REAL> & ));
        
        void GSection(std::vector<LIBMOL::REAL>  & x_a, 
                      std::vector<LIBMOL::REAL>   & x_b, 
                      std::vector<LIBMOL::REAL>   & x_c, 
                      LIBMOL::REAL & obj_a, 
                      LIBMOL::REAL & obj_b, 
                      LIBMOL::REAL & obj_c,
                      LIBMOL::REAL & tObj,
                      std::vector<LIBMOL::REAL>        & tFinCoords,
                      std::vector<LIBMOL::AtomDict>    & tAtoms, 
                      std::vector<LIBMOL::BondDict>    & tBonds, 
                      std::vector<LIBMOL::AngleDict>   & tAngs, 
                      std::vector<LIBMOL::TorsionDict> & tTors, 
                      std::vector<LIBMOL::RingDict>    & tRings, 
                      std::vector<LIBMOL::PlaneDict>   & tPlas, 
                      std::vector<LIBMOL::ChiralDict>  & tChs);
        
        bool checkConvergence(std::vector<LIBMOL::REAL> & tDerivs);
        bool checkConvergence(std::vector<LIBMOL::REAL> & tDerivs,
                              std::vector<LIBMOL::REAL> & tNewCoords,
                              std::vector<LIBMOL::REAL> & tOldCoords);
        
        void updateCoords(std::vector<LIBMOL::AtomDict> & tAtoms,
                          std::vector<LIBMOL::REAL>     & tCoords);
        
        void updateCoordsAndDerivs(std::vector<LIBMOL::REAL>   & tNewCoords,
                                   std::vector<LIBMOL::REAL>   & tOldCoords,
                                   std::vector<LIBMOL::REAL>   & tNewDerivs,
                                   std::vector<LIBMOL::REAL>   & tOldDerivs);
        
        void getFirstDerivs(std::vector<LIBMOL::REAL>        & tDerivs,
                            std::vector<LIBMOL::REAL>        & tCoords,
                            std::vector<LIBMOL::AtomDict>    & tAtoms, 
                            std::vector<LIBMOL::BondDict>    & tBonds, 
                            std::vector<LIBMOL::AngleDict>   & tAngs, 
                            std::vector<LIBMOL::TorsionDict> & tTors, 
                            std::vector<LIBMOL::RingDict>    & tRings, 
                            std::vector<LIBMOL::PlaneDict>   & tPlas, 
                            std::vector<LIBMOL::ChiralDict>  & tChs);
        
        void getFirstDerivs(std::vector<LIBMOL::REAL>        & tDerivs,
                            std::vector<LIBMOL::REAL>        & tCoords,
                            std::vector<LIBMOL::AtomDict>    & tAtoms, 
                            std::vector<LIBMOL::BondDict>    & tBonds, 
                            std::vector<LIBMOL::AngleDict>   & tAngs, 
                            std::vector<LIBMOL::TorsionDict> & tTors, 
                            std::vector<LIBMOL::RingDict>    & tRings, 
                            std::vector<LIBMOL::PlaneDict>   & tPlas, 
                            std::vector<LIBMOL::ChiralDict>  & tChs,
                            std::vector<LIBMOL::AtomDict>    & tAllAtoms);
        
      
        void testCase();
        LIBMOL::REAL getObjValueTest(std::vector<LIBMOL::REAL> & tCoord);
        
        void getDerivsTest(std::vector<LIBMOL::REAL> & tCoords, 
                           std::vector<LIBMOL::REAL> & tDerivs);
        
        int                         nVars;
        int                         nKeep;
        LIBMOL::REAL                xTol;
        LIBMOL::REAL                rTol;
        LIBMOL::REAL                minStep;
        LIBMOL::REAL                maxStep;
        int                         workSpace;
        int                         lComp;
        
    private:
        
        
        int                                         curIter;
        std::map<int, std::vector<LIBMOL::REAL> >   Sk;
        std::map<int, std::vector<LIBMOL::REAL> >   Yk;
        std::vector<LIBMOL::REAL >                  Hk0;
        std::vector<LIBMOL::REAL>                   Pk;
           
        
    };
}

#endif	/* LINALG_H */

