/* 
 * File:   Utility.h
 * Author: flong
 *
 * Created on August 2, 2011, 7:47 PM
 * Last updated July 12, 2012
 * 
 */

#ifndef UTILITY_H
#define	UTILITY_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

// Libs from TNT

#ifndef __JAMA_EIG_H_
#include "jama_eig.h"
#endif

#ifndef __JAMA_CHOLESKY_H_
#include "jama_cholesky.h"
#endif

#ifndef __JAMA_QR_H_
#include "jama_qr.h"
#endif

#ifndef __JAMA_SVD_H_
#include "jama_svd.h"
#endif

// The following macros may not be needed any more with new lib included 
namespace LIBMOL
{
    
    // Data comparison 
    // Data transfering 
    template <class TA, class TB>
    TA Sign(TA num1, TB num2)
    {
        return ( (num2 >= 0) ? num1 : (-num1) );
    }
    
    // Utility function for multiple-key (sortKeys) sorting 
    /*
    template <class TA, class TB>
    bool SortByKeys(TA obj1, TB obj2)
    {
        bool tb = false;
        int i =0;
        while (!tb && i < obj1.sortKeys.size()-1)
        {
          if (obj1.sortKeys[i] < obj2.sortKeys[i])
          {
              return true;
          }
          i++;
        }
        return false;
    }
    
    */
    
    // numerical transform 
    extern bool isInt(std::string tS);
    extern int StrToInt(std::string tS);
    
    extern std::string IntToStr(int tN);
    
    extern REAL StrToReal(std::string tS);
    
    extern std::string RealToStr(REAL tN);
    
    extern REAL StrFractToReal(std::string tS);
    
    extern void getDigitSec(std::string & tStrIn, std::string & tStrOut);
    
    // Handlers with string 
    extern std::string TrimSpaces( std::string  tStr);
    
    extern void cleanChar(std::string & tStr,
                          char          tChar);  // erase all of tChar in tStr
                                                  
    extern void StrUpper(std::string & tStr);
    
    extern void StrLower(std::string & tStr);
    
    extern void  StrTokenize(std::string tStr, 
                             std::vector<std::string>  & tokens);
    
    extern std::vector<std::string> StrTokenize(const std::string &tStr,
                                                char delim);
    extern void  StrTokenize(const std::string &tStr, 
                             std::vector<std::string>  & tV, char delim);
    
    // sort-related functions
    extern bool compareNoCase (std::string first, std::string second);
    extern bool compareNoCase2 (std::string first, std::string second);
    extern bool sortStrByLen(std::string first, std::string second);
    
    extern bool desSortMapKey(const sortMap& a ,const sortMap& b);
    
    extern bool desSortMapKey2(const sortMap2 & a ,const sortMap2 & b);
    
    extern bool desSortMapKey3(const sortMap3 & a ,const sortMap3 & b);
    
    extern bool compareNoCaseClass(const sortLine & tL1, const sortLine & tL2);
    
    extern void cleanSymbol(std::string & tStr, std::string symb);
    
    extern int getKeyWordPos(std::string  tKey, 
                             std::vector<std::string> & tKeyVect);
    
    // File checking related functions, Unix/Linux/MacOs not Windows
    extern int isFileExist(FileName aF);
    
    // trigonometrical functions
    extern REAL degreeToRadians(REAL tDeg);
    extern REAL RadiansToDegree(REAL tRad);
    
    // Math functions for lengths, angles, torsion angles
    // and matrix manipulations 
    
    extern bool inVect(std::vector<REAL>  & tVect,
                       REAL tVal, REAL tErr);
    extern bool inVectABS(std::vector<REAL>  & tVect,
                          REAL  tVal, REAL tErr);
    
    extern bool outVectAbsDiff(std::vector<REAL>  & tVect,
                              REAL  tVal, REAL tErr);
    
    extern void getFracReal(REAL tR1, REAL & tFrac, REAL tTar);
    
    extern REAL lengthV(std::vector<REAL> & tV1);
    
    extern REAL distanceV(std::vector<REAL> & tV1, std::vector<REAL> & tV2);
    
    extern void normalizeV(std::vector<REAL> & tV);
    
    extern REAL dotP(std::vector<REAL> & tV1, std::vector<REAL> & tV2);
    
    extern void crossP2V(std::vector<REAL> & tV1, 
                         std::vector<REAL> & tV2,
                         std::vector<REAL> & tV3);
    
    extern REAL getAngle2V(std::vector<REAL> & tV1, std::vector<REAL> & tV2);
    
    extern REAL getTorsion3V(std::vector<REAL> & tV1, 
                             std::vector<REAL> & tV2,
                             std::vector<REAL> & tV3);
    
    extern REAL getDet(std::vector<REAL> & tV1, 
                       std::vector<REAL> & tV2,
                       std::vector<REAL> & tV3);
    
    extern void matMultVec(std::vector<std::vector<REAL> >    & tMat,
                std::vector<REAL> & tInitV, std::vector<REAL> & tFinV);

    // The following functions are copied from a old version of global optimization codes
    extern REAL length_v(REAL *v, int tDim);
    extern REAL DotP(REAL *v1,  REAL *v2);
    extern REAL DotP(int n_size, REAL *v1,  REAL *v2);
    extern REAL DotP(std::vector<REAL> & v1, REAL *v2);
    extern void CrossP(REAL *v1, REAL *v2, REAL *v3);
    extern REAL GetAngle(REAL *v1, REAL *v2);
    extern REAL GetTorsAngle(REAL *v1, REAL *v2, REAL *v3);
    extern REAL GetTorsAngle(REAL *v1, REAL *v2, REAL *v3);
    
    // matrix operations 
    extern void MatCopy(int n_size,  REAL  ** A_old, REAL  ** A_new);
    extern void MatMultip(int  n_size,   REAL ** AM,
                          REAL * V_old, REAL* V_new);
    
    extern REAL CalcDet(REAL ** Matr);
    extern REAL CalcDet(REAL  * v1,
                        REAL  * v2,
                        REAL  * v3);
    
    extern void MatNormolization(int n_size, REAL  * elemA_dia,
                                 REAL ** A,  REAL  * b);
    extern void MatNormolization(int n_size, REAL ** A);
    extern void MatNormolization(int        n_size,
                                 int        n_dim,
                                 int        n_block,
                                 REAL     * elemA_dia,
                                 REAL     * b, 
                                 REAL     * secDerivCartSparse_row,
                                 REAL     * secDerivCartSparse_col,
                                 REAL   *** secDerivCartSparse);
    extern void UnNormalization (int n_size, REAL *  A_mat_dia,
                                 REAL *  X_vect);
    
    extern void Prod_Mat_Vec( int  n_dim,      int  n_block,
                              REAL    * secDerivCartSparse_row,
                              REAL    * secDerivCartSparse_col,
                              REAL  *** secDerivCartSparse,
                              REAL    * vect_old, 
                              REAL    * vect_new);
     
    extern REAL ConjugateP(int n_size,  REAL * Vect_1,
                           REAL ** A,   REAL *  Vect_2);
    
    extern REAL ConjugateP(int n_size,
                           int n_dim,
                           int n_block,
                           REAL    * secDerivCartSparse_row,
                           REAL    * secDerivCartSparse_col,
                           REAL  *** secDerivCartSparse,
                           REAL *    Vect_1,
                           REAL *    Vect_2);
    
    // Construct a transformation matrix given the components of
    // the new coordinate axis in the old (fixed) coordinate
    extern void InitAMat2(REAL  * X_n,
                          REAL  * Y_n,
                          REAL  * Z_n,
                          REAL ** a_mat,
                          int     t_dim);
    
    extern void NBD_MCOPY(REAL ** A, 
                         REAL ** B, 
                         int     t_dim);
    
    extern void  NBD_MATMLT(REAL ** a_mat, 
                            REAL ** Y,
                            REAL ** Z, 
                            int     t_dim);
    
    extern void EigenSolve(int  mat_size,  REAL ** mat_in,
                           REAL *  eig_value, REAL ** eig_vect);
    
    //extern void Plane_Find(int    idx_pl,    REAL   *coeffs,
    //                       REAL   &coeff_d);
    
    extern void Plane_Find(int         num_atoms, REAL   **coords_atoms,
                           REAL   *coeffs);
    extern void Plane_and_Deriv_Find(int  num_atoms, REAL  **coords_atoms,
                                     REAL   *coeffs, REAL **d_coeffs);    
    
    // some math functions such as Random number generator
    extern REAL GetRand();
    extern void initPrimeTab(std::vector<int> & tPrimeTab,
                             std::string tLibmolTabDir);
    
    // Chemical properties related functions 
    extern void initOrgTable(std::vector<std::string> & tOrgTab);
    extern void initMetalTab(std::vector<ID> & tMeTab);
    extern bool isMetal(std::vector<ID> & tMeTab, ID tID);
    extern void fromIdToChemType(ID tId, ID & tChemType);
    extern REAL StrToOrder(std::string  & tStrOrder);
    extern void OrderToStr(REAL tOrder, std::string  & sOrder);
    extern void OrderStrToStr(std::string & tOrder, std::string  & sOrder);
    extern void unifyStrForOrder(std::string & tOrder);
    extern void ChiToStr(int & tCIdx, std::string & tCStr);
    extern REAL strToCharge(std::string & tStr);  
    
    // symmetry related operations 
    extern void StrToSymmOps(std::vector<std::string>           & tStrs, 
                             std::vector<std::vector<REAL> >    & tMat);
    extern void StrToSymmOneRow(std::string       & tStr,
                                std::vector<REAL> & tRow);
    extern void StrToSymmOneRow2(std::string       & tStr,
                                 std::vector<REAL> & tRow);
    //extern void FractToOrtho(std::vector<REAL> tFractCoords,
    //                        std::vector<REAL> tOrthoCoords,
    //                        REAL a, REAL b, REAL c,
    //                        REAL alpha, REAL beta, REAL gamma);
    
    extern void TranslateIntoUnitCell(std::vector<REAL> & tInitFracX,
                                    std::vector<REAL> & tFinFracX);
    
    extern void TranslateIntoUnitCell(std::vector<REAL> & tFracX);
    
    
    extern void FractToOrtho(std::vector<REAL> & tFractCoords,
                             std::vector<REAL> & tOrthoCoords,
                             REAL a, REAL b, REAL c,
                             REAL alpha, REAL beta, REAL gamma);
    
    extern void OrthoToFract(std::vector<REAL> tOrthoCoords,
                             std::vector<REAL> tFractCoords,
                             REAL a,     REAL b,    REAL c,
                             REAL alpha, REAL beta, REAL gamma);
    
    // Function dealt with environment variables and scripts
    // executing external program using command line arguments
    extern void cmdExecute(std::string & tCom);
    
    extern std::string getLibmolDir();

}

#endif	/* UTILITY_H */

