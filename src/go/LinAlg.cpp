
/* 
 * File:   LinAlg.cpp
 * Author: flong
 *
 * Created on February 12, 2013, 3:55 PM
 * Last updated on February 12, 2013 11:46 PM
 * 
 */

#include "LinAlg.h"
#include "atom.h"
#include "TransCoord.h"

namespace GO
{
    FindLocalMin::FindLocalMin():workMode(1),
            workSpace(1),
            lComp(100),
            maxIter(50),
            maxSumDX(1.0e-8),
            maxTol(1.0e-4),
            curObjValue(1.0e6)
            
    {
    }
    
    FindLocalMin::~FindLocalMin()
    {
    }
    
    void FindLocalMin::Driver(std::vector<LIBMOL::AtomDict>   & tAtoms, 
                              std::vector<LIBMOL::BondDict>   & tBonds, 
                              std::vector<LIBMOL::AngleDict>  & tAngs, 
                              std::vector<LIBMOL::TorsionDict>& tTors, 
                              std::vector<LIBMOL::RingDict>   & tRings, 
                              std::vector<LIBMOL::PlaneDict>  & tPlas, 
                              std::vector<LIBMOL::ChiralDict> & tChs,
                              std::vector<LIBMOL::AtomDict>   & tAllAtoms)
                              
    {
        int numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
        
        
        if (workSpace==1)
        {   
            if (workMode==11)
            {
                FF::GetAllDerivsFull   derivs(numVars);
                derivs.workMode  = workMode;
                derivs.workSpace = workSpace;
        
                LIBMOL::REAL * deltaX  = new LIBMOL::REAL [numVars];
                
                int k=0;
              
                LIBMOL::NeighbListDict  tNBListOfSystem;
        
                int aDim  = 3;
                int aMode = 0;
                LIBMOL::REAL tCellL     = 3.0;
                LIBMOL::REAL tCellShell = 0.5;
        
                tNBListOfSystem.building(tAtoms, aDim, tCellL, tCellShell, aMode);
                maxIter = 3;
                do 
                {
                    for (int i =0; i < numVars; i++)
                    {
                        deltaX[i] =0.0;
                    }
                    
                    derivs.setFirstDerivsCart(tAtoms, tBonds, tAngs,
                                              tTors,  tPlas, tChs);
                    /*
                    for (int i=0; i < numVars; i++)
                    {
                        if (fabs(derivs.firDrivCart[i]) >0.2)
                        {
                            std::cout << "1st deriv [" <<i << "]=" << derivs.firDrivCart[i] << std::endl;
                        }
                    }
                    
                    std::cout << "done 1st deriv " << std::endl;
                    */
                    derivs.SetNormalMatrix(tAtoms, tBonds, tAngs, tTors, 
                                           tPlas, tChs);
                    
                    // std::cout << "Done one normal matrix " << std::endl;
                    
                    CGGEqs(numVars, derivs.secDrivCart, derivs.firDrivCart, deltaX);
                    /*
                    for(int i=0; i < numVars; i++)
                    {
                        std::cout << "deltaX[" << i << "]=" << deltaX[i] << std::endl;
                    }
                     */
                    // std::cout << "Done cg set" << std::endl;
                    
                    LinMin(deltaX, tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
                    //std::cout << "Done line min, now the obj value is " 
                    //          << curObjValue << std::endl;
                    
                    //checkConvergence(k, tTol, tSumDx);
                    //std::cout << "Done one local min set, " <<  k 
                    //          << " cycles" << std::endl;
                    k++;
                }while(k <maxIter);  
                //}while(k < maxIter && tTol < maxTol && tSumDx < maxSumDX);
                 
                delete [] deltaX;
                deltaX  = NULL;    
            }
            else if (workMode==1)
            {       
                
                // Building an initial neighbor list of atoms
                LIBMOL::NeighbListDict  tNBListOfSystem;
        
                int aDim  = 3;
                int aMode = 0;
                LIBMOL::REAL tCellL     = 3.0;
                LIBMOL::REAL tCellShell = 0.5;
                
                // std::cout << "Number of atoms in Subsystem " << (int)tAtoms.size() << std::endl;
                
                tNBListOfSystem.building(tAtoms, tAllAtoms, aDim, tCellL, tCellShell, aMode);
                
                // std::cout << "NB List built " << std::endl;
                
                // set initial values 
                
                std::vector<LIBMOL::REAL> iniCrds;
                
                
                for (int i=0; i < (int)tAtoms.size(); i++)
                {
                    int iDim = (int)tAtoms[i].coords.size();
                    for (int j=0; j < iDim; j++)
                    {
                        iniCrds.push_back(tAtoms[i].coords[j]);
                    }
                }
                
                
                LIBMOL::REAL finObj=0.0;
                
                // Local min using L-BFGS
                int m1 = int(numVars*0.2);
                int m2 = int(numVars*0.5);
                int m  = 0;
                if (m1 > 10)
                {
                    m = m1;
                }
                else if (m2 > 10)
                {
                    m = m2;
                }
                else
                {
                    m = numVars;
                }
                
                
                std::cout << "Search Local mim by L-BFGS " << std::endl;
                
                
                
                LBFGS aLBFGSObj(numVars, m, 1.0e-6, 1.0e-4, 1.0e-6, 1.0);
                aLBFGSObj.lComp = lComp;           
                aLBFGSObj.driver(finObj, iniCrds, tAtoms, tBonds, tAngs,
                                 tTors,  tRings, tPlas, tChs, tAllAtoms);
                
                curObjValue = finObj;
                //std::cout << "After a local min, the objective value is "
                //         << curObjValue << std::endl;
                
            }
            
        }
        
        
       
    }
    
    void FindLocalMin::Driver(std::vector<LIBMOL::AtomDict>   & tAtoms, 
                              std::vector<LIBMOL::BondDict>   & tBonds, 
                              std::vector<LIBMOL::AngleDict>  & tAngs, 
                              std::vector<LIBMOL::TorsionDict>& tTors, 
                              std::vector<LIBMOL::RingDict>   & tRings, 
                              std::vector<LIBMOL::PlaneDict>  & tPlas, 
                              std::vector<LIBMOL::ChiralDict> & tChs)
                              
    {
        int numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
        
        //std::cout <<"Number of atoms is " << (int)tAtoms.size() << std::endl
        //          <<"dim of coordinates is " << (int)tAtoms[0].coords.size() << std::endl
        //          << "total number of variates is " << numVars << std::endl;
        //std::cout << "work space is " << workSpace << std::endl;
        //std::cout << "workMode is " << workMode << std::endl;
        
        if (workSpace==1)
        {   
            if (workMode==11)
            {
                FF::GetAllDerivsFull   derivs(numVars);
                derivs.workMode  = workMode;
                derivs.workSpace = workSpace;
        
                LIBMOL::REAL * deltaX  = new LIBMOL::REAL [numVars];
                
                int k=0;
              
                LIBMOL::NeighbListDict  tNBListOfSystem;
        
                int aDim  = 3;
                int aMode = 0;
                LIBMOL::REAL tCellL     = 3.0;
                LIBMOL::REAL tCellShell = 0.5;
        
                tNBListOfSystem.building(tAtoms, aDim, tCellL, tCellShell, aMode);
                maxIter = 3;
                do 
                {
                    for (int i =0; i < numVars; i++)
                    {
                        deltaX[i] =0.0;
                    }
                    
                    derivs.setFirstDerivsCart(tAtoms, tBonds, tAngs,
                                              tTors,  tPlas, tChs);
                    /*
                    for (int i=0; i < numVars; i++)
                    {
                        if (fabs(derivs.firDrivCart[i]) >0.2)
                        {
                            std::cout << "1st deriv [" <<i << "]=" << derivs.firDrivCart[i] << std::endl;
                        }
                    }
                    
                    std::cout << "done 1st deriv " << std::endl;
                    */
                            
                    derivs.SetNormalMatrix(tAtoms, tBonds, tAngs, tTors, 
                                           tPlas, tChs);
                    
                    // std::cout << "Done one normal matrix " << std::endl;
                    
                    CGGEqs(numVars, derivs.secDrivCart, derivs.firDrivCart, deltaX);
                    /*
                    for(int i=0; i < numVars; i++)
                    {
                        std::cout << "deltaX[" << i << "]=" << deltaX[i] << std::endl;
                    }
                     */
                    // std::cout << "Done cg set" << std::endl;
                    
                    LinMin(deltaX, tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
                    // std::cout << "Done line min, now the obj value is " 
                    //          << curObjValue << std::endl;
                    
                    //checkConvergence(k, tTol, tSumDx);
                    //std::cout << "Done one local min set, " <<  k 
                    //          << " cycles" << std::endl;
                    k++;
                }while(k <maxIter);  
                //}while(k < maxIter && tTol < maxTol && tSumDx < maxSumDX);
                 
                delete [] deltaX;
                deltaX  = NULL;    
            }
            else if (workMode==1)
            {       
                
                // Building an initial neighbor list of atoms
                LIBMOL::NeighbListDict  tNBListOfSystem;
        
                int aDim  = 3;
                int aMode = 0;
                LIBMOL::REAL tCellL     = 3.0;
                LIBMOL::REAL tCellShell = 0.5;
                
                //std::cout << "Number of atoms in Subsystem " << (int)tAtoms.size() << std::endl;
                
        
                tNBListOfSystem.building(tAtoms, aDim, tCellL, tCellShell, aMode);
                
                //std::cout << "NB List built " << std::endl;
               
                
                // set initial values 
                
                std::vector<LIBMOL::REAL> iniCrds;
                
                
                for (int i=0; i < (int)tAtoms.size(); i++)
                {
                    int iDim = (int)tAtoms[i].coords.size();
                    for (int j=0; j < iDim; j++)
                    {
                        iniCrds.push_back(tAtoms[i].coords[j]);
                    }
                }
                
                
                LIBMOL::REAL finObj=0.0;
                
                // Local min using L-BFGS
                int m1 = int(numVars*0.2);
                int m2 = int(numVars*0.5);
                int m  = 0;
                if (m1 > 10)
                {
                    m = m1;
                }
                else if (m2 > 10)
                {
                    m = m2;
                }
                else
                {
                    m = numVars;
                }
                
                // std::cout << "Search Local mim by L-BFGS " << std::endl;
                
                
                LBFGS aLBFGSObj(numVars, m, 1.0e-6, 1.0e-4, 1.0e-6, 1.0);
                aLBFGSObj.lComp = lComp; 
                
                aLBFGSObj.driver(finObj, iniCrds, tAtoms, tBonds, tAngs,
                                 tTors,  tRings, tPlas, tChs);
                
                curObjValue = finObj;
                
                // std::cout << "After a local min, the objective value is "
                //           << curObjValue << std::endl;
                
            }
            
        }
        
        
       
    }
    //        Solve a linear equation by the conjugate
    //        gradient method. The solution found here
    //        will be used as the increase step for
    //        a line minimization in a static minimization
    //        problem. Full matrix version 
    int  FindLocalMin::CGGEqs(int            n_dim,
                              LIBMOL::REAL** sec_deriv, 
                              LIBMOL::REAL*  firs_deriv, 
                              LIBMOL::REAL*  delta_x)
    {
          LIBMOL::REAL * A_dia;
          A_dia = new LIBMOL::REAL [n_dim];
  

          LIBMOL::REAL * res_cur;
          LIBMOL::REAL * res_pre;
          LIBMOL::REAL * res_prepre;

          LIBMOL::REAL * p_cur;
          LIBMOL::REAL * p_pre;

          LIBMOL::REAL   alpha;
          LIBMOL::REAL   beta;

          LIBMOL::REAL tol_t = 1.0e-6;
          LIBMOL::REAL sum_res, sum_fd;

          LIBMOL::REAL tm1, tm2, tm3;
          LIBMOL::REAL *tm4;

          res_cur      = new LIBMOL::REAL [n_dim];
          res_pre      = new LIBMOL::REAL [n_dim];
          res_prepre   = new LIBMOL::REAL [n_dim];

          p_cur        = new LIBMOL::REAL [n_dim];  
          p_pre        = new LIBMOL::REAL [n_dim]; 

          tm4          = new LIBMOL::REAL [n_dim];
          
          //  Precondition

          int  n_cond = 1;         // should be an input parameters
      
          switch(n_cond)
          {
              case 1:
                  
                  // Normalization of the matrix A
                  for (int i= 0; i < n_dim; i++)
                  {
                      A_dia[i] = sqrt(sec_deriv[i][i]*sec_deriv[i][i]);
                      if(fabs(A_dia[i]) <1.0e-4)
                      {
                          A_dia[i] = 1.0;
                      }
                  }
                  
                  LIBMOL::MatNormolization(n_dim, A_dia, sec_deriv, firs_deriv);
                  
                  break;
              
              default :  
                  std::cout << "No precondition is done" << std::endl;
                  break;      
          }
          
          // Stablization

          LIBMOL::REAL alpha_1 = 0.02;

          Stab_plusAlpha(n_dim, alpha_1, sec_deriv);
  
          // Initial guess

          for (int i =0; i < n_dim; i++)
          {
              res_cur[i]       = 0.0;
              res_pre[i]       = 0.0;
              res_prepre[i]    = 0.0;
              p_cur [i]        = 0.0;
              p_pre [i]        = 0.0;
              tm4[i]           = 0.0;     
          }

          sum_fd = 0.0;
          for (int i =0; i <n_dim; i++)
          {
              sum_fd += pow(firs_deriv[i],2.0);
          }
          sum_fd = sqrt(sum_fd);
          
          
         
 
          alpha      = 0.0;
          
          //  MatMultip(n_dim_t,sec_dri, delta_x, res_cur); 

          for (int i =0; i < n_dim; i++)
          {
              res_cur[i] = -res_cur[i] + firs_deriv[i];
              //std::cout << "1st deriv " << firs_deriv[i] << std::endl;
              //std::cout << res_cur[i] << std::endl;
          }
          
          sum_res = 0.0;

          for (int i =0; i < n_dim; i++)
          {
              sum_res += pow(res_cur[i], 2.0);
          }
          sum_res = sqrt(sum_res);
          /*
          if(sum_res <1.0e-16)
          {
              std::cout << "CGGEqs failed " << std::endl;
              return 0;
          }
          */
          // Conjugate gradient search of the solution

          int k =0;
          
            do{
                k ++;

                if ( k == 1)
                {
                    for (int i =0; i < n_dim; i++)
                    {
                        p_cur [i]  = res_cur[i];
                        res_pre[i] = res_cur[i];
                    }
                    tm2 = LIBMOL::DotP(n_dim, res_pre, res_pre);
                }
                else 
                {
                    beta = 0.0;

                    tm1 = LIBMOL::DotP(n_dim, res_prepre, res_prepre);
                    tm2 = LIBMOL::DotP(n_dim, res_pre, res_pre);
                    
                    if(tm1)
                    {
                        beta = tm2/tm1;
                    }
                    else
                    {
                        std::cout << " should be stopped in the previous steps. BUGS !" 
                                  << std::endl;
                    }
          
                    for (int i =0; i < n_dim; i++)
                    {
                        p_cur [i] = res_pre[i] + beta*p_pre[i];
                    }
                }

             
                tm3 = LIBMOL::ConjugateP(n_dim, p_cur, sec_deriv, p_cur);
                // std::cout << "tm2 " << tm2 << std::endl;
                // std::cout << "tm3 = " << tm3 << std::endl;

                if(fabs(tm3) > 1.0e-16)
                {
                    alpha = fabs(tm2/tm3);
                    //  std::cout << "alpha = " << alpha << std::endl;
                }
                else
                {
                    // the case that is near a sigularity
                    DelAconjSingularity(p_cur, sec_deriv, p_cur);

                    tm3    = LIBMOL::ConjugateP(n_dim, p_cur, sec_deriv, p_cur);
                    alpha  = tm2/tm3;

                }
              
                for(int i =0; i < n_dim; i++)
                {
                    //std::cout << "befor, delta_x["<<i<<"] =" 
                    //          << delta_x[i] << std::endl;
                    //std::cout << "alpha" << alpha << std::endl;
                    //std::cout << p_cur[i] << std::endl;
                    
                    delta_x[i] = delta_x[i] + alpha*p_cur[i];
                    
        
                    //std::cout << "after, delta_x["<<i<<"] =" 
                    //          << delta_x[i] << std::endl;
                }       
                     
                LIBMOL::MatMultip(n_dim,sec_deriv, p_cur, tm4);

                for (int i =0; i < n_dim; i++)
                {
                    res_cur [i] =  res_pre[i] - alpha * tm4[i];
                    // std::cout << "tm["<<i<<"] = " << tm4[i] << std::endl;
                         
                }

                //   Copy the current values to the previous values

                for (int i =0; i < n_dim; i++)
                {
                    res_prepre[i] =  res_pre[i];
                    res_pre[i]    =  res_cur[i];
                    p_pre [i]     =  p_cur[i];
                }

                sum_res = 0.0;
                for (int i =0; i < n_dim; i++)
                {        
                    sum_res += pow(res_cur[i], 2.0);
                }
                sum_res = sqrt(sum_res);      
         
                //  if(myFFSpace == 2)
                //  {
                //    cout << " sum_res " << sum_res << endl;
                //    cout << " tol_t " << tol_t << " sum_fd " << sum_fd << endl;
                //    cout << "Continue ? " << endl;
                //    cin.get();
                //  }

            }while(sum_res > tol_t*sum_fd && k < 100);
            
            // Undoing precondition
            
            switch(n_cond)
            {
                case 1: 
                    LIBMOL::UnNormalization(n_dim, A_dia, delta_x);
                    break;
                    
                default :  
                    std::cout << "You forget select the way of doing precondition " 
                              << std::endl;
                    break;      
            }
  
            // Cancel the Stablization factor

            Cancel_plusAlpha(n_dim, alpha_1, delta_x);

            //cout << "number of iterations in CGG " << k << endl;

            // Check if the elements of the shift delta_x are too small
            // to kick off a further line minimization 
  
            int vs = CheckShift(delta_x,n_dim);

            // cout << " delta_x after linear solution
            
            // release memory and null all local points
            
            delete [] A_dia;
            A_dia = NULL;

            delete [] res_cur;
            res_cur = NULL;

            delete []  res_pre;
            res_pre =NULL;

            delete []  res_prepre;
            res_prepre = NULL;

            delete [] p_cur;
            p_cur =NULL;
  
            delete []  p_pre;
            p_pre =NULL;

            delete [] tm4;
            tm4 =NULL;

            //pPrecondition   = 0;
            //pUnPrecondition = 0;
            
            return vs;
    }
    
    // The version CGGEqs with Sparse matrix format  
    int  FindLocalMin::CGGEqs(int n_size,
                              int n_dim,
                              int n_block,
                              LIBMOL::REAL    * firs_deriv, 
                              LIBMOL::REAL    * secDerivCartSparse_row,
                              LIBMOL::REAL    * secDerivCartSparse_col,
                              LIBMOL::REAL  *** secDerivCartSparse,
                              LIBMOL::REAL    * delta_x)
    {
        
        int i,j,k;
        int n_dim_all = n_dim*n_block;

        LIBMOL::REAL * res_cur;
        LIBMOL::REAL * res_pre;
        LIBMOL::REAL * res_prepre;

        LIBMOL::REAL * p_cur;
        LIBMOL::REAL * p_pre;

        LIBMOL::REAL alpha;
        LIBMOL::REAL beta;

        LIBMOL::REAL tol_t = 1.0e-6;
        LIBMOL::REAL sum_res, sum_fd;

        LIBMOL::REAL tm1, tm2, tm3;
        LIBMOL::REAL *tm4;

        res_cur      = new LIBMOL::REAL [n_dim_all];
        res_pre      = new LIBMOL::REAL [n_dim_all];
        res_prepre   = new LIBMOL::REAL [n_dim_all];

        p_cur        = new LIBMOL::REAL [n_dim_all];  
        p_pre        = new LIBMOL::REAL [n_dim_all]; 

        tm4          = new LIBMOL::REAL [n_dim_all];
  
        // diagonal element of the A_matrix


        LIBMOL::REAL * A_dia;
        A_dia = new LIBMOL::REAL [n_dim_all];
 
        for (int idx_atom = 0; idx_atom < n_block; idx_atom++)
        {
            int i_row_s = secDerivCartSparse_row[idx_atom];
            for(i = 0; i < n_dim; i++)
            {
                j = idx_atom*n_dim+i;
                A_dia[j] = sqrt (secDerivCartSparse[i_row_s][i][i]);
                if(fabs(A_dia[j]) <1.0e-8)
                {
                    A_dia[j] = 1.0e-8;
                }
            }
        }

        //  Precondition  

        int  n_cond = 1;         // should be an input parameters
    
        switch(n_cond)
        {
            case 1: 
                // Normalization of A_matrix
                LIBMOL::MatNormolization(n_size, n_dim, n_block, A_dia, firs_deriv,
                                         secDerivCartSparse_row, secDerivCartSparse_col,
                                         secDerivCartSparse);
                break;
            default :  
                std::cout << "No precondition is done" << std::endl;
                break;      
        }
        
        // Stablization

        LIBMOL::REAL alpha_1 = 0.02;

        Stab_plusAlpha(alpha_1, n_dim, n_block, 
                       secDerivCartSparse_row, secDerivCartSparse);
  
        // Initial guess

        for (i =0; i < n_dim_all; i++)
        {
            res_cur[i]       = 0.0;
            res_pre[i]       = 0.0;
            res_prepre[i]    = 0.0;
            p_cur [i]        = 0.0;
            p_pre [i]        = 0.0;
            tm4[i]           = 0.0;     
        }
        
        sum_fd = 0.0;
        for (i =0; i < n_dim_all; i++)
        {
            sum_fd += pow(firs_deriv[i],2.0);
        }
        sum_fd = sqrt(sum_fd);
  
        //cout << "sum_fd " << sum_fd << endl; 
        alpha      = 0.0; 

        //Prod_Mat_Vec(delta_x, res_cur); 

        for (i =0; i < n_dim_all; i++)
        {
            //std::cout << "before Init res_cur[i] " << res_cur[i] << std::endl; 
            res_cur[i] = -res_cur[i] + firs_deriv[i];
            // cout << "firs_driv[i] " <<  firs_driv[i] << endl; 
            // cout << "Init res_cur[i] " << res_cur[i] << endl;     
        }

        sum_res = 0.0;

        for (i =0; i < n_dim_all; i++)
        {
            sum_res += pow(res_cur[i], 2.0);
        }
        sum_res = sqrt(sum_res);

        // cout << " Init sum_residue " <<  sum_res << endl; 

        // Conjugate gradient search of the solution

        k =0;     
        
        do{
            
            k ++;
            
            if ( k == 1)
            {
                for (i =0; i < n_dim_all; i++)
                {
                    p_cur [i]  = res_cur[i];
                    res_pre[i] = res_cur[i];
                }
                
                tm2 = LIBMOL::DotP(n_dim_all, res_pre, res_pre);
            }
            else 
            {
                
                beta = 0.0;

                tm1 = LIBMOL::DotP(n_dim_all,res_prepre, res_prepre);
                tm2 = LIBMOL::DotP(n_dim_all,res_pre, res_pre);
                
                if(tm1)
                {
                    beta = tm2/tm1;
                }
                else
                { 
                    std::cout << " should be stopped in the previous steps. BUGS !" 
                              << std::endl;
                }
          
                for (i =0; i < n_dim_all; i++)
                {
                    p_cur [i] = res_pre[i] + beta*p_pre[i];
                }
            }
      
        
            tm3 = LIBMOL::ConjugateP(n_dim_all, n_dim, n_block, 
                                     secDerivCartSparse_row,
                                     secDerivCartSparse_col, 
                                     secDerivCartSparse,
                                     p_cur, p_cur);
        
            // cout << "tm3 = " << tm3 << endl;

            if(fabs(tm3) > 1.0e-16)
            {
                alpha = tm2/tm3;
                //  cout << "alpha = " << alpha << endl;
            }
            else
            {
                // the case that is near a sigularity
                //DelAconjSingularity(p_cur, sec_dri, p_cur);
                //cout << "tm3 too small, check" << endl;
                //exit(1);
            
                // tm3    = ConjugateP(n_dim_all, p_cur, p_cur);
                alpha  = tm2/tm3;

            }
              
            for( i =0; i < n_dim_all; i++)
            {
                // std;:cout << "befor, delta_x["<<i<<"] =" << delta_x[i] << std::endl;
                delta_x[i] = delta_x[i] + alpha*p_cur[i];
                // std::cout << "after, delta_x["<<i<<"] =" << delta_x[i] << std::endl;
            }      
                    
            LIBMOL::Prod_Mat_Vec(n_dim, n_block, secDerivCartSparse_row,
                                 secDerivCartSparse_col, secDerivCartSparse,
                                 p_cur, tm4);

            for ( i =0; i < n_dim_all; i++)
            {
                res_cur [i] =  res_pre[i] - alpha * tm4[i];
                //  std::cout << "tm["<<i<<"] = " << tm4[i] << std::endl;
                         
            }

     
            //   Copy the current values to the previous values

            for ( i =0; i < n_dim_all; i++)
            {
                res_prepre[i] =  res_pre[i];
                res_pre[i]    =  res_cur[i];
                p_pre [i]     =  p_cur[i];

            }

            sum_res = 0.0;
            for (i =0; i < n_dim_all; i++)
            {        
                sum_res += pow(res_cur[i], 2.0);
            }
            sum_res = sqrt(sum_res);      

            //     std::cout << " sum_residue " << sum_res << std::endl;
            //     std::cout << " tol_t " << tol_t << " sum_fd " << sum_fd << std::endl;
            //     std::cout << "Continue " << std::endl;
            //     std::cin.get();
        
        }while(sum_res > tol_t*sum_fd && k < 100);
        
        // Undoing precondition
        
        switch(n_cond)
        {
            case 1: 
                LIBMOL::UnNormalization(n_dim_all, A_dia, delta_x);
                break;
            
            default :  
                std::cout << "You forget select the way of doing precondition " 
                          << std::endl;
                break;      
        }
        
        // Cancel the stablization factor
        
        Cancel_plusAlpha(n_dim_all, alpha_1, delta_x);

        // std::cout << "number of iterations in CGG " << k << std::endl;

        // Check if the elements of the shift delta_x are too small
        // to kick off a further line minimization 
  
        int vs = CheckShift(delta_x,n_dim_all);

        // release memory and null all local points

        delete [] A_dia;
        A_dia = NULL;

        delete [] res_cur;
        res_cur = NULL;

        delete []  res_pre;
        res_pre =NULL;

        delete []  res_prepre;
        res_prepre = NULL;

        delete [] p_cur;
        p_cur = NULL;
  
        delete []  p_pre;
        p_pre = NULL;

        delete [] tm4;
        tm4 = NULL;
  
        return vs;
    }
    
    
    // Stablization  of a matrix, A_matrix.
    void FindLocalMin::Stab_plusAlpha(int             n_size, 
                                      LIBMOL::REAL    alp, 
                                      LIBMOL::REAL ** A_matrix)
    {

         for (int i = 0; i < n_size; i++)
         {
             A_matrix [i][i] *= (1+alp);
         }
    }
    
    // Stablization  of a matrix, A_matrix, of sparse matrix format.
    // Using private variables
    void FindLocalMin::Stab_plusAlpha(LIBMOL::REAL alp,
                                      int          n_dim,
                                      int          n_block,
                                      LIBMOL::REAL *     secDerivCartSparse_row,
                                      LIBMOL::REAL ***   secDerivCartSparse)
    {

          for (int i_ato = 0; i_ato < n_block; i_ato++)
          {
              int i_block = secDerivCartSparse_row[i_ato];
              for(int i =0; i < n_dim; i++)
              {
                  secDerivCartSparse [i_block][i][i] *=(1+alp);
              }
          }
    }

    void FindLocalMin::Cancel_plusAlpha(int           n_size, 
                                        LIBMOL::REAL  alph, 
                                        LIBMOL::REAL* vect)
    {
        for (int i = 0; i < n_size; i++)
        {
            vect [i] *= (1+alph);
        }
    }
   
  
    int FindLocalMin::CheckShift(LIBMOL::REAL* d_x, int n_vars)
    {
        LIBMOL::REAL tol_s = 1.0e-3;
        LIBMOL::REAL max_s = 0.0;

        for(int i = 0; i < n_vars; i++)
        {
            //  std::cout << "delta_x["<<i+1<<"] = " << d_x[i]  << std::endl;
            if(fabs(d_x[i]) >  max_s)
            {
                max_s = fabs(d_x[i]);
            }
        }

     
        int myExitFF = 0;
  
        if( max_s < tol_s )
        {
            myExitFF = 1;
        }
        
        return myExitFF;
        
    }
    
    
    void FindLocalMin::DelAconjSingularity(LIBMOL::REAL* vect_1,
                                           LIBMOL::REAL** matr, 
                                           LIBMOL::REAL* vect_2)
    {
        
    }
    
    void FindLocalMin::LinMin(LIBMOL::REAL* delta_x, 
                              std::vector<LIBMOL::AtomDict>    & tAtoms, 
                              std::vector<LIBMOL::BondDict>    & tBonds, 
                              std::vector<LIBMOL::AngleDict>   & tAngs, 
                              std::vector<LIBMOL::TorsionDict> & tTors, 
                              std::vector<LIBMOL::RingDict>    &  tRings, 
                              std::vector<LIBMOL::PlaneDict>   & tPlas, 
                              std::vector<LIBMOL::ChiralDict>  & tChs)
    {
          int i;
          
          int numVars;
          
          
          LIBMOL::REAL obja, objb,objc;

          if(workSpace == 1)
          {
              numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
          }
          else if (workSpace == 2)
          {
              numVars = (int)tTors.size();
          }
          else
          {
              std::cout << "How to do Flash-freezing ? " << std::endl;
              exit(1);
          }

          LIBMOL::REAL * xa= new LIBMOL::REAL  [numVars];
          LIBMOL::REAL * xb= new LIBMOL::REAL   [numVars];
          LIBMOL::REAL * xc= new LIBMOL::REAL   [numVars];

          for (i =0; i < numVars; i++)
          {
              xa[i] =0.0;
              xb[i] =0.0;
              xc[i] =0.0;
          }
          
          Bracket(numVars, delta_x, xa, xb, xc,obja, objb, objc, 
                  tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);

          // cout << "Three bracket points are " << endl;

          //  for (i =0; i < numVars; i++)
          //   {
          //     cout << "xa["<<i<<"]=" << xa[i] << endl;
          //     cout << "xb["<<i<<"]=" << xb[i] << endl;
          //     cout << "xc["<<i<<"]=" << xc[i] << endl;
          //   } 

          // cout << "Three objective values are " << endl;
          // cout << " obj_a = " << obja << endl;
          // cout << " obj_b = " << objb << endl;
          // cout << " obj_c = " << objc << endl;        
      

          GSection(xa, xb, xc, obja, objb, objc,
                   tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
 
          // cout << "my obj after line search " <<myObjValue << endl;
          
          delete [] xa;
          xa =NULL;

          delete [] xb;
          xb =NULL;
  
          delete [] xc;
          xc =NULL;
    }
    
    
    void FindLocalMin::Bracket(int num_vars, LIBMOL::REAL* delta_x, 
                               LIBMOL::REAL* x_a, LIBMOL::REAL* x_b, LIBMOL::REAL* x_c, 
                               LIBMOL::REAL& obj_a, LIBMOL::REAL& obj_b, LIBMOL::REAL& obj_c, 
                               std::vector<LIBMOL::AtomDict>   & tAtoms, 
                               std::vector<LIBMOL::BondDict>   & tBonds, 
                               std::vector<LIBMOL::AngleDict>  & tAngs, 
                               std::vector<LIBMOL::TorsionDict>& tTors, 
                               std::vector<LIBMOL::RingDict>   & tRings, 
                               std::vector<LIBMOL::PlaneDict>  & tPlas, 
                               std::vector<LIBMOL::ChiralDict> & tChs)
    {
        FF::GetObjValue toolGetObjValue;
        toolGetObjValue.workSpace = workSpace;
        
        int dim = (int)tAtoms[0].coords.size();
        int numAtoms = (int)tAtoms.size();
        int i,k;
        int numVars=dim*numAtoms;

        LIBMOL::REAL * x_r;
        LIBMOL::REAL * x_q;
        LIBMOL::REAL * x_u;
        LIBMOL::REAL * x_ulim;
  
        LIBMOL::REAL * x_temp;

        LIBMOL::REAL glim = 100.0;

        LIBMOL::REAL obj_u, obj_temp;

        LIBMOL::REAL s_b, s_c, s_u, s_ulim;
  
        LIBMOL::REAL tm1;
  
        x_r    = new LIBMOL::REAL [numVars];
        x_q    = new LIBMOL::REAL [numVars];
        x_u    = new LIBMOL::REAL [numVars];
        x_ulim = new LIBMOL::REAL [numVars];

        x_temp = new LIBMOL::REAL [numVars];

        // Set an initial region
        
        obj_a = toolGetObjValue.getAll(tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        if(workSpace == 1)
        {
            for (int i =0; i < numAtoms; i++)
            {
                for (int j =0; j < dim; j++)
                {
                    k= i*dim+j;
                    x_a[k] = tAtoms[i].coords[j];
                    x_b[k] = tAtoms[i].coords[j] + delta_x[k];
                }
            }
        }
        else if(workSpace == 2)
        {
            for(int i = 0; i < (int)tTors.size(); i++)
            {
                x_a[i] = tTors[i].value;
                x_b[i] = tTors[i].value + delta_x[i];
            }
        }

        obj_b =  toolGetObjValue.getAll(x_b, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
  

        if(obj_b > obj_a)
        {
            for (int i =0; i < numVars; i++)
            {
                x_temp[i] = x_a[i];
                x_a[i]    = x_b[i];
                x_b[i]    = x_temp[i];
            }

            obj_temp = obj_a;
            obj_a    = obj_b;
            obj_b    = obj_temp;

        }

        // First guess for c

        for (i = 0; i < numVars; i++)
        {
            x_c[i] =  x_b[i]+(1+GOLD)*(x_b[i] - x_a[i]);
        }
    
        obj_c =  toolGetObjValue.getAll(x_c, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);

        //std::cout << "input stage in Bracket" << std::endl;
        //std::cout << "obj_a =" << obj_a << std::endl;
        //std::cout << "obj_b =" << obj_b << std::endl;
        // cout << "obj_c =" << obj_c << endl;
  
        // Loop over until find lower obj_c
        // int kk =0;

        // Loop over until find lower obj_c
        // int kk =0;

        while (obj_b > obj_c)
        {
            // kk++;
            // cout << "kk = " << kk << endl;

            for (i =0; i <  numVars; i++)
            {
                x_r[i] = (x_b[i] - x_a[i])*(obj_b-obj_c);
                x_q[i] = (x_b[i] - x_c[i])*(obj_b-obj_a);
                tm1     = x_q[i] - x_r[i];

                if(fabs(tm1) < 1.0e-6)
                { 
                    // cout << "tm1 = " << tm1 << endl;
                    tm1 = 1.0e-6;
                }
                
                x_u[i] = x_b[i] - ((x_b[i]-x_c[i])*x_q[i]
                        -(x_b[i]-x_a[i])*x_r[i])/(2.0*tm1);

                //    cout << "x_u["<<i<<"] = " << x_u[i] << endl;   
         
                x_ulim[i] = x_b[i] - glim*(x_c[i]-x_b[i]);

            }
      
            s_b =0.0;
            s_c =0.0; 
            s_u =0.0;
            s_ulim = 0.0;
            
            
            for (i =0; i < numVars; i++)
            {
                s_b += pow((x_b[i]-x_a[i]),2.0);
                s_c += pow((x_c[i]-x_a[i]),2.0);
                s_u += pow((x_u[i]-x_a[i]),2.0);
                s_ulim += pow((x_ulim[i]-x_a[i]),2.0);
            }
            s_b = sqrt(s_b);
            s_c = sqrt(s_c);
            s_u = sqrt(s_u);
            s_ulim = sqrt(s_ulim);
      
     
            // check the possible minimization region
 
            if((s_b - s_u)*(s_u-s_c) > 0.0 ) // paraabolic u is between b and c
            {
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);

                if(obj_u < obj_c)           // get a minimum between b and c
                {
                    for (i = 0; i < numVars; i++)
                    {
                        x_a[i] = x_b[i];
                        x_b[i] = x_u[i];
                    }
              
                    obj_a = obj_b;
                    obj_b = obj_u;

                    // cout << "final stage in Bracket " << endl;
                    // cout << "obj_a =" << obj_a << endl;
                    // cout << "obj_b =" << obj_b << endl;
                    // cout << "obj_c =" << obj_c << endl;
       
                    return;

                }
                else if ( obj_u > obj_b)   // Got a minimum between a and u
                {
                    for (i = 0; i < numVars; i++)
                    {
                        x_c[i] = x_u[i];
                    }
              
                    obj_c = obj_u;

                    //  cout << "final stage in Bracket " << endl;
                    //  cout << "obj_a =" << obj_a << endl;
                    //  cout << "obj_b =" << obj_b << endl;
                    //  cout << "obj_c =" << obj_c << endl;
                    return;
                }
   

                for (i = 0; i < numVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }

                obj_u = toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);

            }
            else if ((s_c - s_u)*(s_u -s_ulim) > 0.0) // parabolic is between x_c
            {                                    // and its allowed limit
                obj_u = toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                
                if(obj_u < obj_c)
                {
                    for (i = 0; i < numVars; i++)
                    {
                        x_b[i] = x_c[i];
                        x_c[i] = x_u[i];
                        x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                    }
               
                    obj_b = obj_c;
                    obj_c = obj_u;
                    obj_u = toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);


                }
            }
            else if ((s_u -s_ulim)*(s_ulim-s_c) >=0.0) // Limit parabolic u 
            {                                        // to maximum allowed value
                for (i = 0; i < numVars; i++)
                {
                    x_u[i] = x_ulim[i];
                }
                obj_u = toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
            }
            else                                       // Reject parabolic u,
            {                                          // use default magnification
                for (i = 0; i < numVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }

                obj_u = toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
            }

            for (i = 0; i < numVars; i++)
            {
                x_a[i] = x_b[i];         
                x_b[i] = x_c[i];
                x_c[i] = x_u[i];         
            }
            
            obj_a = obj_b;   
            obj_b = obj_c;
            obj_c = obj_u;
            
        }

        //  cout << "final stage in Bracket " << endl;
        //  cout << "obj_a =" << obj_a << endl;
        //  cout << "obj_b =" << obj_b << endl;
        //  cout << "obj_c =" << obj_c << endl;
        
        delete [] x_r;
        x_r = NULL;

        delete [] x_q;
        x_q = NULL;

        delete [] x_u;
        x_u = NULL;

        delete [] x_ulim;
        x_ulim = NULL;

        delete [] x_temp;
        x_temp = NULL;
  
    }
    
    // Find a lower value of the objective function within
    // the bracket range
    void FindLocalMin::GSection(LIBMOL::REAL* x_a, LIBMOL::REAL* x_b, 
                                LIBMOL::REAL* x_c, LIBMOL::REAL& obj_a, 
                                LIBMOL::REAL& obj_b, LIBMOL::REAL& obj_c, 
                                std::vector<LIBMOL::AtomDict>& tAtoms, 
                                std::vector<LIBMOL::BondDict>& tBonds, 
                                std::vector<LIBMOL::AngleDict>& tAngs, 
                                std::vector<LIBMOL::TorsionDict>& tTors, 
                                std::vector<LIBMOL::RingDict> & tRings, 
                                std::vector<LIBMOL::PlaneDict>& tPlas, 
                                std::vector<LIBMOL::ChiralDict>& tChs)
    {
        
        int i,j;
        int numVars;

        LIBMOL::REAL tol = 1.0e-5;
        
        LIBMOL::TransCoords transTool;
        
        FF::GetObjValue toolGetObjValue;
        toolGetObjValue.workSpace = workSpace;
        

        LIBMOL::REAL *x0;
        LIBMOL::REAL *x1;
        LIBMOL::REAL *x2;
        LIBMOL::REAL *x3;

        LIBMOL::REAL obj1, obj2;

        LIBMOL::REAL s_b, s_c, s_x1, s_x2, s_x3;

        if(workSpace == 1)
        {
            numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
        }
        else if (workSpace == 2)
        {
            numVars = (int)tTors.size();
        }
        else
        {
            std::cout << "In which space to do Flash-freezing ? " << std::endl;
            exit(1);
        }

        x0 = new LIBMOL::REAL [numVars];
        x1 = new LIBMOL::REAL [numVars];
        x2 = new LIBMOL::REAL [numVars];
        x3 = new LIBMOL::REAL [numVars];
  
  
        s_b =0.0;
        s_c =0.0;

        for( i =0; i < numVars; i++)
        {
            x0[i] = x_a[i];
            x3[i] = x_c[i];
            s_b += pow((x_b[i]-x_a[i]),2.0);
            s_c += pow((x_c[i]-x_a[i]),2.0);
        }
     
        s_b = sqrt(s_b);
        s_c = sqrt(s_b);
   
        if(fabs(s_c-s_b) > fabs(s_b) ) 
        {
            for( i =0; i < numVars; i++)
            {
                x1[i] = x_b[i];
                x2[i] = x_b[i]+(1-GOLD)*(x_c[i]-x_b[i]);
            }
        }
        else
        {
            for( i =0; i < numVars; i++)
            {
                x2[i] = x_b[i];
                x1[i] = x_b[i]-(1-GOLD)*(x_b[i]-x_a[i]);
            }
        }

        obj1 = toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        obj2 = toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);

        // cout << "obj1 in section " << obj1 << endl;
        // cout << "obj2 in section " << obj2 << endl;

        s_x1 = 0.0;
        s_x2 = 0.0;
        s_x3 = 0.0;

        for( i =0; i < numVars; i++)
        {
            s_x1+=pow((x1[i]-x0[i]),2.0);
            s_x2+=pow((x2[i]-x0[i]),2.0);
            s_x3+=pow((x3[i]-x0[i]),2.0);
        }

        s_x1 = sqrt(s_x1);
        s_x2 = sqrt(s_x2);
        s_x3 = sqrt(s_x3);
 
        //  cout << "Init sections are  s_x1 " << s_x1 << endl
        //     << " s_x2  " << s_x2 << endl
        //     << " s_x3  " << s_x3 << endl;

        j =0;
 
        while ((s_x3 > tol*(s_x1+s_x2)) && j < 20)
        {
            if(obj2 < obj1)
            {
                //  cout << "obj1 > obj2 " << endl;
                for( i =0; i < numVars; i++)
                {
                    x0[i] = x1[i];
                    x1[i] = x2[i];
                    x2[i] = GOLD*x1[i]+(1-GOLD)*x3[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
	            //   cout << "x2["<<i<<"]= "<< x2[i] << endl;             
                }
           
                obj1 = obj2;
                obj2 = toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
	        //	  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;
            }
            else
            {
                //  cout << "obj2 < obj1 " << endl;
                for( i =0; i < numVars; i++)
                {
                    x3[i] = x2[i];
                    x2[i] = x1[i];
                    x1[i] = GOLD*x2[i]+(1-GOLD)*x0[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
                    //  cout << "x2["<<i<<"]= "<< x2[i] << endl;  
                }
           
                obj2 = obj1;
                obj1 = toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
	        //  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;       
            }

            s_x1 = 0.0;
            s_x2 = 0.0;
            s_x3 = 0.0;

            for( i =0; i < numVars; i++)
            {
                s_x1+=pow((x1[i]-x0[i]),2.0);
                s_x2+=pow((x2[i]-x0[i]),2.0);
                s_x3+=pow((x3[i]-x0[i]),2.0);
            }

            s_x1 = sqrt(s_x1);
            s_x2 = sqrt(s_x2);
            s_x3 = sqrt(s_x3);
            j++;
            // cout << "j = " << j << endl;
        }

        // cout << "new obj1 is " << obj1 << endl;
        // cout << "new obj2 is " << obj2 << endl;

        // LIBMOL::REAL obj1_t, obj2_t;
  
        if(workSpace == 1 )
        {
            if(obj1 < obj2)
            {
               //obj1_t = toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
               curObjValue = obj1;
               UpdateAtomCoords(tAtoms, x1);
               
	       // cout << "Obj1_t = " << obj1_t << endl;
               // cout << "myObj  = " << myObjValue << endl;
            }
            else 
            {
                //obj2_t = toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                curObjValue = obj2;
                UpdateAtomCoords(tAtoms, x2);
	        //      cout << "Obj2_t = " << obj2_t << endl;
                //      cout << "myObj  = " << myObjValue << endl;
            }
        }
        else if(workSpace == 2)
        {
            if(obj1 < obj2)
            { 
                for(i = 0; i < (int)tTors.size(); i++)
                {
                    tTors[i].value = x1[i];
                }      
                // UpdatTors_t();  
                //CoordAngToCart();
                curObjValue = obj1;
            }
            else
            {
                for(i = 0; i < (int)tTors.size(); i++)
                {
                    tTors[i].value = x2[i];
                }  
         
                // UpdatTors_t();  
                curObjValue = obj2;
            }
            transTool.generateCoordTorsToCart2(tAtoms, tBonds, tAngs, tTors, 
                                              tRings, tPlas, tChs);
            
        }             
        
        //  std::cout << "my obj after line search " <<curObjValue << std::endl;
        //  std::cout << "Continue " << endl;
        //  std::cin.get();

        delete [] x0;
        x0 = NULL;
        delete [] x1;
        x1 = NULL;
        delete x2;
        x2 = NULL;
        delete x3;
        x3 = NULL;
 
    }
   
    void FindLocalMin::checkConvergence(int idx_it, 
                                        LIBMOL::REAL tol_t, 
                                        LIBMOL::REAL& tol_dx)
    {
        
    }
 
    void FindLocalMin::TransAtomCoordsGenToCart(std::vector<LIBMOL::AtomDict>& tAtoms, 
                                                LIBMOL::REAL* tX) 
    {
        int tDim = (int)tAtoms[0].coords.size();
        for (int i=0; i < (int)tAtoms.size(); i++)
        {
            for (int j=0; j < (int)tAtoms[i].coords.size(); j++)
            {
                int k = i*tDim+j;
                tAtoms[i].coords[j] = tX[k];
            }
        }
    }
    
    void FindLocalMin::TransAtomCoordsCartToGen(std::vector<LIBMOL::AtomDict>& tAtoms, 
                                                LIBMOL::REAL* tX)
    {
        int tDim = (int)tAtoms[0].coords.size();
        for (int i=0; i < (int)tAtoms.size(); i++)
        {
            for (int j=0; j < (int)tAtoms[i].coords.size(); j++)
            {
                int k = i*tDim+j;
                tX[k] = tAtoms[i].coords[j];
            }
        }
    }
    
    LIBMOL::REAL FindLocalMin::GetObjectValue(LIBMOL::REAL  * tX, 
                                              std::vector<LIBMOL::AtomDict>& tAtoms, 
                                              std::vector<LIBMOL::BondDict>& tBonds, 
                                              std::vector<LIBMOL::AngleDict>& tAngs, 
                                              std::vector<LIBMOL::TorsionDict>& tTors, 
                                              std::vector<LIBMOL::RingDict>& tRings, 
                                              std::vector<LIBMOL::PlaneDict>& tPlas, 
                                              std::vector<LIBMOL::ChiralDict>& tChs)
    {
        FF::GetObjValue tObj;
        tObj.workSpace = workSpace;
        
        if (workSpace==1)
        {
        
            for (int i=0; i < (int)tAtoms.size(); i++)
            {
                int dim =(int)tAtoms[i].coords.size(); 
                for (int j=0; j < dim; j++)
                {
                    int k=i*dim+j;
                    tAtoms[i].coords[j] = tX[k];
                }
            }
        }
        else if (workSpace ==2)
        {
            
        }
        return tObj.getAll(tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
    }
    
    
    
    void FindLocalMin::UpdateAtomCoords(std::vector<LIBMOL::AtomDict> tAtoms, 
                                        LIBMOL::REAL* tX)
    {
        
        for (std::vector<LIBMOL::AtomDict>::iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            int dim = (int)iA->coords.size();
            for (int i=0; i <dim; i++)
            {
                int j = iA->seriNum*dim + i;
                iA->coords[i] = tX[j];
            }
        }
    }
    
    LBFGS::LBFGS():nVars(ZeroInt),
            nKeep(ZeroInt),
            xTol(1.0e-6),
            rTol(1.0e-6),
            minStep(1.0e-6),
            maxStep(1.0),
            workSpace(1),
            lComp(100),
            curIter(ZeroInt)
    {
    }
    
    LBFGS::LBFGS(int tVars, int tKeep, LIBMOL::REAL txTol, LIBMOL::REAL trTol, 
                 LIBMOL::REAL tMinStep, LIBMOL::REAL tMaxStep)
    {
        nVars   = tVars;
        nKeep   = tKeep;
        xTol    = txTol;
        rTol    = trTol;
        minStep = tMinStep;
        maxStep = tMaxStep;
        curIter = 0;
        workSpace =1;
        lComp     =100;
    }
    
    LBFGS::~LBFGS()
    {
    }
    
    void LBFGS::driver(LIBMOL::REAL  & tFinObj, 
                       std::vector<LIBMOL::REAL>& tFinCoords, 
                       std::vector<LIBMOL::REAL>& tIniCoords,
                       std::vector<LIBMOL::REAL>& tIniGrads)
    {
        std::vector<LIBMOL::REAL> tNewCoords;
        std::vector<LIBMOL::REAL> tOldCoords;
        std::vector<LIBMOL::REAL> tNewDerivs;
        std::vector<LIBMOL::REAL> tOldDerivs;
       
        bool                      lCont= false;
        
        nVars = (int)tIniCoords.size();
        for (int i=0; i < nVars; i++)
        {
            tOldCoords.push_back(tIniCoords[i]);
            tOldDerivs.push_back(tIniGrads[i]);
            //std::cout << "x["<<i<<"]=" << tIniCoords[i] << std::endl;
            //std::cout << "x["<<i<<"]=" << tIniGrads[i] << std::endl;
        }
        
        LIBMOL::normalizeV(tOldDerivs);
        
        
        do 
        {
            // std::cout << "curIter " << curIter << std::endl;
            if (curIter==0)
            {
                init(tOldDerivs);
            }
            else if (curIter==1)
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
            }
            else
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
                updateCoordsAndDerivs(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
            }
            
            lineSearch(tNewCoords, tOldCoords, tOldDerivs, tFinObj);
            
            getDerivsTest(tNewCoords, tNewDerivs);
            
            lCont = checkConvergence(tNewDerivs);
            
            curIter++;
            
        }while(!lCont && curIter < 1000);
        
        tFinCoords.clear();
        for (int i=0; i < (int)tNewCoords.size(); i++)
        {
            tFinCoords.push_back(tNewCoords[i]);
        }
       
    }
    
    void LBFGS::driver(LIBMOL::REAL  & tFinObj, 
                       std::vector<LIBMOL::REAL>        & tIniCoords,
                       std::vector<LIBMOL::AtomDict>    & tAtoms, 
                       std::vector<LIBMOL::BondDict>    & tBonds, 
                       std::vector<LIBMOL::AngleDict>   & tAngs, 
                       std::vector<LIBMOL::TorsionDict> & tTors, 
                       std::vector<LIBMOL::RingDict>    & tRings, 
                       std::vector<LIBMOL::PlaneDict>   & tPlas, 
                       std::vector<LIBMOL::ChiralDict>  & tChs)
    {
        
        
        std::vector<LIBMOL::REAL> tNewCoords;
        std::vector<LIBMOL::REAL> tOldCoords;
        std::vector<LIBMOL::REAL> tNewDerivs;
        std::vector<LIBMOL::REAL> tOldDerivs;
       
        bool                      lCont= false;
        
        nVars = (int)tIniCoords.size();
        for (int i=0; i < nVars; i++)
        {
            tOldCoords.push_back(tIniCoords[i]);
            // std::cout << "x["<<i<<"]=" << tIniCoords[i] << std::endl;
        }
        
        
        getFirstDerivs(tOldDerivs, tIniCoords, 
                       tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
        
        // std::cout << "1st deriv done " << std::endl;
       
        
        LIBMOL::normalizeV(tOldDerivs);
        
        
        do 
        {
            
            if (curIter==0)
            {
                init(tOldDerivs);
            }
            else if (curIter==1)
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
            }
            else
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
                updateCoordsAndDerivs(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
            }
            
            // lineSearch(tNewCoords, tOldCoords, tOldDerivs, tFinObj);
            lineSearch(tNewCoords, tOldCoords, tOldDerivs, tFinObj, 
                       tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
            
            getFirstDerivs(tNewDerivs, tNewCoords, 
                           tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
            lCont = checkConvergence(tNewDerivs, tNewCoords, tOldCoords);
            
            curIter++;
      
            
        }while(!lCont && curIter < 200);
        
        
        
        updateCoords(tAtoms, tNewCoords);
    }
    
    void LBFGS::driver(LIBMOL::REAL  & tFinObj, 
                       std::vector<LIBMOL::REAL>        & tIniCoords,
                       std::vector<LIBMOL::AtomDict>    & tAtoms, 
                       std::vector<LIBMOL::BondDict>    & tBonds, 
                       std::vector<LIBMOL::AngleDict>   & tAngs, 
                       std::vector<LIBMOL::TorsionDict> & tTors, 
                       std::vector<LIBMOL::RingDict>    & tRings, 
                       std::vector<LIBMOL::PlaneDict>   & tPlas, 
                       std::vector<LIBMOL::ChiralDict>  & tChs,
                       std::vector<LIBMOL::AtomDict>    & tAllAtoms)
    {
        
        std::vector<LIBMOL::REAL> tNewCoords;
        std::vector<LIBMOL::REAL> tOldCoords;
        std::vector<LIBMOL::REAL> tNewDerivs;
        std::vector<LIBMOL::REAL> tOldDerivs;
       
        bool                      lCont= false;
        
        nVars = (int)tIniCoords.size();
        for (int i=0; i < nVars; i++)
        {
            tOldCoords.push_back(tIniCoords[i]);
        }
     
        
        getFirstDerivs(tOldDerivs, tIniCoords, 
                       tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs, tAllAtoms);
        
        // std::cout << "1st deriv done " << std::endl;
        // exit(1);
        
        LIBMOL::normalizeV(tOldDerivs);
        
        
        do 
        {
            
            if (curIter==0)
            {
                init(tOldDerivs);
            }
            else if (curIter==1)
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
            }
            else
            {
                setSandY(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
                setHk0();
                setPk(tNewDerivs);
                updateCoordsAndDerivs(tNewCoords, tOldCoords, tNewDerivs, tOldDerivs);
            }
            
            // lineSearch(tNewCoords, tOldCoords, tOldDerivs, tFinObj);
            lineSearch(tNewCoords, tOldCoords, tOldDerivs, tFinObj, 
                       tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
            
            getFirstDerivs(tNewDerivs, tNewCoords, 
                           tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
            lCont = checkConvergence(tNewDerivs, tNewCoords, tOldCoords);
            
            curIter++;
      
            
        }while(!lCont && curIter < 100);
        
        
        
        updateCoords(tAtoms, tNewCoords);
    }
    
    
    void LBFGS::init(std::vector<LIBMOL::REAL>& iniGrads)
    {
        
        
        // the first iteration
        for (int i=0; i < nVars; i++)
        {
            Hk0.push_back(1.0);
        }
        
        for (int i=0; i < nVars; i++)
        {
            Pk.push_back(iniGrads[i]*Hk0[i]);
        }
    }
    
    void LBFGS::init(int tKeep, LIBMOL::REAL txTol, LIBMOL::REAL trTol, 
                     LIBMOL::REAL tMinStep, LIBMOL::REAL tMaxStep, 
                     std::vector<LIBMOL::REAL>& iniGrads)
    {
        setParams(tKeep, txTol, trTol, tMinStep, tMaxStep);
        
        init(iniGrads);
        
    }
    
    void LBFGS::setParams(int tKeep, LIBMOL::REAL txTol, LIBMOL::REAL trTol, 
                          LIBMOL::REAL tMinStep, LIBMOL::REAL tMaxStep)
    {
        nKeep   = tKeep;
        xTol    = txTol;
        rTol    = trTol;
        minStep = tMinStep;
        maxStep = tMaxStep;
    }
    
    void LBFGS::setHk0()
    {
        Hk0.clear();
        
        LIBMOL::REAL gammaK=1.0, A=0.0;
        
        A=LIBMOL::dotP(Yk[curIter-1], Yk[curIter-1]);
        
        if (fabs(A) >1.0e-8)
        {
            gammaK=LIBMOL::dotP(Sk[curIter-1], Yk[curIter-1])/A;
        }
        
        for (int i=0; i < nVars; i++)
        {
            Hk0.push_back(gammaK);
        }
    }
    
    void LBFGS::setSandY(std::vector<LIBMOL::REAL>& tNewCoords, 
                         std::vector<LIBMOL::REAL>& tOldCoords, 
                         std::vector<LIBMOL::REAL>& tNewDerivs, 
                         std::vector<LIBMOL::REAL>& tOldDerivs)
    {
        for(int i=0; i < (int)tOldCoords.size(); i++)
        {
            /*
            if(lComp==4 && curIter > 50)
            {
                
                exit(1);
            }
            else
            {
               std::cout << "curIter " <<  curIter << std::endl;
               std::cout << tNewCoords[i]-tOldCoords[i] << std::endl; 
            }
            */
            Sk[curIter-1].push_back(tNewCoords[i]-tOldCoords[i]);
            Yk[curIter-1].push_back(tNewDerivs[i]-tOldDerivs[i]);
            
        }
        
        if ((int)Sk.size() > nKeep)
        {
            Sk.erase(curIter-nKeep-1);
            Yk.erase(curIter-nKeep-1);
        }
    }
    
    void LBFGS::setPk(std::vector<LIBMOL::REAL>   & tDerivs)
    {
        int lowB;
        if((curIter-nKeep) >0)
        {
            lowB=curIter-nKeep;
        }
        else
        {
            lowB=0;
        }
       
        std::vector<LIBMOL::REAL> q;
        for (int i=0; i < (int)tDerivs.size(); i++)
        {
            q.push_back(tDerivs[i]);
        }
        
        std::map<int, LIBMOL::REAL> alphaK; 
        
        for (int i=curIter-1; i >=lowB; i--)
        {
            LIBMOL::REAL rhoK;
            LIBMOL::REAL A=LIBMOL::dotP(Sk[i], Yk[i]);
            if (A > 1.0e-8)
            {
                rhoK=1.0/A;
            }
            else
            {
                rhoK=1.0;
            }
            
            alphaK[i] = rhoK*LIBMOL::dotP(Sk[i], q);
            
            for (int j=0; j < (int)q.size(); j++)
            {
                q[j]-=(alphaK[i]*Yk[i][j]);
            }
        }
        
        
        
        Pk.clear();
        for(int i=0; i < (int)q.size(); i++)
        {
            Pk.push_back(Hk0[i]*q[i]);
        }
        
        for (int i=lowB; i <= curIter-1; i++)
        {
            LIBMOL::REAL rhoK;
            LIBMOL::REAL A=LIBMOL::dotP(Sk[i], Yk[i]);
            if (A > 1.0e-8)
            {
                rhoK=1.0/A;
            }
            else
            {
                rhoK=1.0;
            }
            
            LIBMOL::REAL betaK=rhoK*LIBMOL::dotP(Sk[i], Pk);
            
            for (int j=0; j < (int)Pk.size(); j++)
            {
                Pk[j] =Pk[j] + Sk[i][j]*(alphaK[i]-betaK);
            }
        }
        
        // Normalize Pk
        LIBMOL::REAL pL=LIBMOL::lengthV(Pk);
        
        for (int i=0; i < (int)Pk.size(); i++)
        {
            Pk[i]= Pk[i]/pL;
        }
        
    }
 
    void LBFGS::lineSearch(std::vector<LIBMOL::REAL>& tFinCoords, 
                           std::vector<LIBMOL::REAL>& tIniCoords, 
                           std::vector<LIBMOL::REAL>& tIniDerivs,
                           LIBMOL::REAL             & tObj)
    {
        // Initialize the shifts
        /*
        LIBMOL::REAL iniDGx=0.0;
        for (int i=0; i < nVars; i++)
        {
            iniDGx = iniDGx + tIniDeriv[i]*Pk[i];
        }
        
        if (iniDGx < 0.0)
        {
            std::cout << "Initial shift direction is not a descent direction"
                      << std::endl;
        }
        
        LIBMOL::REAL iniObj= tObj;
        
        std::vector<LIBMOL::REAL> tX;
        
        for (std::vector<LIBMOL::REAL>::iterator iX=tIniCoords.begin();
                iX != tIniCoords.end(); iX++)
        {
            tX.push_back(*iX);
        }
        
        
        LIBMOL::REAL bestStep = 0.0;
        LIBMOL::REAL bestObj  = tObj;
        LIBMOL::REAL bestDGx  = iniDGx;
        
        LIBMOL::REAL STY      = 0.0;
        LIBMOL::REAL FY       = tObj;
        LIBMOL::REAL DGY      = iniDGx;
        
        LIBMOL::REAL P5  = 0.50, P66 = 0.66, XTRAPF=0.40;
        
        */
        
        // Normalize derivs
        
      
        
          
          LIBMOL::REAL obja, objb,objc;
          /*
          if(workSpace == 1)
          {
              numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
          }
          else if (workSpace == 2)
          {
              numVars = (int)tTors.size();
          }
          else
          {
              std::cout << "How to do Flash-freezing ? " << std::endl;
              exit(1);
          }

          LIBMOL::REAL * xa= new LIBMOL::REAL  [numVars];
          LIBMOL::REAL * xb= new LIBMOL::REAL   [numVars];
          LIBMOL::REAL * xc= new LIBMOL::REAL   [numVars];

          for (i =0; i < numVars; i++)
          {
              xa[i] =0.0;
              xb[i] =0.0;
              xc[i] =0.0;
          }
          */
          
          std::vector<LIBMOL::REAL> xa, xb, xc;
          for (int i=0; i < nVars; i++)
          {
              xa.push_back(tIniCoords[i]);
              // std::cout << "xa " << i << "  " <<  xa[i] << std::endl;
              xb.push_back(tIniCoords[i]+tIniDerivs[i]*Pk[i]);
              // std::cout << "xb " << i << "  " <<  xb[i] << std::endl;
              xc.push_back(0.0);
          }
          
          Bracket(xa, xb, xc,obja, objb, objc);
          /*
          // cout << "Three bracket points are " << endl;

          //  for (i =0; i < numVars; i++)
          //   {
          //     cout << "xa["<<i<<"]=" << xa[i] << endl;
          //     cout << "xb["<<i<<"]=" << xb[i] << endl;
          //     cout << "xc["<<i<<"]=" << xc[i] << endl;
          //   } 
          
          std::cout << "Three objective values are " << std::endl;
          std::cout << " obj_a = " << getObjValueTest(xa) << std::endl;
          std::cout << " obj_b = " << getObjValueTest(xb) << std::endl;
          std::cout << " obj_c = " << getObjValueTest(xc) << std::endl;        
      
          */
          
          GSection(xa, xb, xc, obja, objb, objc, tObj, tFinCoords, &LBFGS::getObjValueTest);
    }
    
    void LBFGS::lineSearch(std::vector<LIBMOL::REAL>        & tFinCoords, 
                           std::vector<LIBMOL::REAL>        & tIniCoords, 
                           std::vector<LIBMOL::REAL>        & tIniDerivs,
                           LIBMOL::REAL                     & tObj,
                           std::vector<LIBMOL::AtomDict>    & tAtoms, 
                           std::vector<LIBMOL::BondDict>    & tBonds, 
                           std::vector<LIBMOL::AngleDict>   & tAngs, 
                           std::vector<LIBMOL::TorsionDict> & tTors, 
                           std::vector<LIBMOL::RingDict>    & tRings, 
                           std::vector<LIBMOL::PlaneDict>   & tPlas, 
                           std::vector<LIBMOL::ChiralDict>  & tChs)
    {  
            LIBMOL::REAL obja, objb,objc;
            std::vector<LIBMOL::REAL> xa, xb, xc;
            for (int i=0; i < nVars; i++)
            {
                xa.push_back(tIniCoords[i]);
                // std::cout << "xa " << i << "  " <<  xa[i] << std::endl;
                xb.push_back(tIniCoords[i]+tIniDerivs[i]*Pk[i]);
                // std::cout << "xb " << i << "  " <<  xb[i] << std::endl;
                xc.push_back(0.0);
            }
          
            Bracket(xa, xb, xc,obja, objb, objc,
                    tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
            
            GSection(xa, xb, xc, obja, objb, objc, tObj, tFinCoords, 
                    tAtoms, tBonds, tAngs, tTors, tRings, tPlas, tChs);
            
        }
    
    
    void LBFGS::Bracket(std::vector<LIBMOL::REAL>  & x_a, 
                        std::vector<LIBMOL::REAL>  & x_b, 
                        std::vector<LIBMOL::REAL>  & x_c, 
                        LIBMOL::REAL& obj_a, 
                        LIBMOL::REAL& obj_b, 
                        LIBMOL::REAL& obj_c)
    
    {
        /*
        FF::GetObjValue toolGetObjValue;
        toolGetObjValue.workSpace = workSpace;
        
        int dim = (int)tAtoms[0].coords.size();
        int numAtoms = (int)tAtoms.size();
        int i,k;
        int numVars=dim*numAtoms;

        LIBMOL::REAL * x_r;
        LIBMOL::REAL * x_q;
        LIBMOL::REAL * x_u;
        LIBMOL::REAL * x_ulim;
  
        LIBMOL::REAL * x_temp;

        LIBMOL::REAL glim = 100.0;

        LIBMOL::REAL obj_u, obj_temp;

        LIBMOL::REAL s_b, s_c, s_u, s_ulim;
  
        LIBMOL::REAL tm1;
  
        x_r    = new LIBMOL::REAL [numVars];
        x_q    = new LIBMOL::REAL [numVars];
        x_u    = new LIBMOL::REAL [numVars];
        x_ulim = new LIBMOL::REAL [numVars];

        x_temp = new LIBMOL::REAL [numVars];
         
         */
        

        LIBMOL::REAL glim = 100.0;

        LIBMOL::REAL obj_u, obj_temp;

        LIBMOL::REAL s_b, s_c, s_u, s_ulim;
  
        LIBMOL::REAL tm1;
        
        std::vector<LIBMOL::REAL> x_r, x_q, x_u, x_ulim, x_temp;
        
        for (int i=0; i < nVars; i++)
        {
            x_r.push_back(0.0);
            x_q.push_back(0.0);
            x_u.push_back(0.0);
            x_ulim.push_back(0.0);
            x_temp.push_back(0.0);
        }
        
        

        // Set an initial region
        
        //obj_a = toolGetObjValue.getAll(tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        obj_a = getObjValueTest(x_a);
        
        /*        
        if(workSpace == 1)
        {
            for (int i =0; i < numAtoms; i++)
            {
                for (int j =0; j < dim; j++)
                {
                    k= i*dim+j;
                    x_a[k] = tAtoms[i].coords[j];
                    x_b[k] = tAtoms[i].coords[j] + delta_x[k];
                }
            }
        }
        else if(workSpace == 2)
        {
            for(int i = 0; i < (int)tTors.size(); i++)
            {
                x_a[i] = tTors[i].value;
                x_b[i] = tTors[i].value + delta_x[i];
            }
        }
       
        obj_b =  toolGetObjValue.getAll(x_b, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        */
        obj_b = getObjValueTest(x_b);
        // std::cout << "obj_a " << obj_a << " and obj_b " << obj_b << std::endl;
        
        if(obj_b > obj_a)
        {
            for (int i =0; i < nVars; i++)
            {
                x_temp[i] = x_a[i];
                x_a[i]    = x_b[i];
                x_b[i]    = x_temp[i];
            }

            obj_temp = obj_a;
            obj_a    = obj_b;
            obj_b    = obj_temp;

        }
           

        // First guess for c
        
        
        for (int i = 0; i < nVars; i++)
        {
            x_c[i] =  x_b[i]+(1+GOLD)*(x_b[i] - x_a[i]);
            // std::cout << "x_c " << i << " is " << x_c[i] << std::endl;
            
        }
        
        // obj_c =  toolGetObjValue.getAll(x_c, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        obj_c = getObjValueTest(x_c);
        //std::cout << "input stage in Bracket" << std::endl;
        //std::cout << "obj_a =" << obj_a << std::endl;
        //std::cout << "obj_b =" << obj_b << std::endl;
        //std::cout << "obj_c =" << obj_c << std::endl;
      
        
  
        // Loop over until find lower obj_c
        // int kk =0;

        // Loop over until find lower obj_c
        // int kk =0;

        while (obj_b > obj_c)
        {
            // kk++;
            // cout << "kk = " << kk << endl;

            for (int i =0; i <  nVars; i++)
            {
                x_r[i] = (x_b[i] - x_a[i])*(obj_b-obj_c);
                x_q[i] = (x_b[i] - x_c[i])*(obj_b-obj_a);
                tm1    = x_q[i] - x_r[i];

                if(fabs(tm1) < 1.0e-6)
                { 
                    // cout << "tm1 = " << tm1 << endl;
                    tm1 = 1.0e-6;
                }
                
                x_u[i] = x_b[i] - ((x_b[i]-x_c[i])*x_q[i]
                        -(x_b[i]-x_a[i])*x_r[i])/(2.0*tm1);

                //    cout << "x_u["<<i<<"] = " << x_u[i] << endl;   
         
                x_ulim[i] = x_b[i] - glim*(x_c[i]-x_b[i]);

            }
            
           
            s_b =0.0;
            s_c =0.0; 
            s_u =0.0;
            s_ulim = 0.0;
            
            
            for (int i =0; i < nVars; i++)
            {
                s_b += pow((x_b[i]-x_a[i]),2.0);
                s_c += pow((x_c[i]-x_a[i]),2.0);
                s_u += pow((x_u[i]-x_a[i]),2.0);
                s_ulim += pow((x_ulim[i]-x_a[i]),2.0);
            }
            s_b = sqrt(s_b);
            s_c = sqrt(s_c);
            s_u = sqrt(s_u);
            s_ulim = sqrt(s_ulim);
      
            
            
            std::vector<LIBMOL::REAL> x_bu, x_uc, x_uulim;
            for (int i=0; i < nVars; i++)
            {
                x_bu.push_back(x_b[i]-x_u[i]);
                x_uc.push_back(x_u[i]-x_c[i]);
                x_uulim.push_back(x_u[i]-x_ulim[i]);
            }
            
            // check the possible minimization region
 
            //if((s_b - s_u)*(s_u-s_c) > 0.0 ) // paraabolic u is between b and c
            if (LIBMOL::dotP(x_bu, x_uc) > 0.0)
            {
                obj_u =  getObjValueTest(x_u);

                if(obj_u < obj_c)           // get a minimum between b and c
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_a[i] = x_b[i];
                        x_b[i] = x_u[i];
                    }
              
                    obj_a = obj_b;
                    obj_b = obj_u;

                    //std::cout << "final stage in Bracket 1" << std::endl;
                    //std::cout << "obj_a =" << obj_a << std::endl;
                    //std::cout << "obj_b =" << obj_b << std::endl;
                    //std::cout << "obj_c =" << obj_c << std::endl;
       
                    return;

                }
                else if ( obj_u > obj_b)   // Got a minimum between a and u
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_c[i] = x_u[i];
                    }
              
                    obj_c = obj_u;

                    //std::cout << "final stage in Bracket " << std::endl;
                    //std::cout << "obj_a =" << obj_a << std::endl;
                    //std::cout << "obj_b =" << obj_b << std::endl;
                    //std::cout << "obj_c =" << obj_c << std::endl;
                    return;
                }
   

                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }

                obj_u = getObjValueTest(x_u);

            }
            //else if ((s_c - s_u)*(s_u -s_ulim) > 0.0) // parabolic is between x_c
            else if (-LIBMOL::dotP(x_uc, x_uulim) >0.0)
            {                                    // and its allowed limit
                obj_u = getObjValueTest(x_u);
                
                if(obj_u < obj_c)
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_b[i] = x_c[i];
                        x_c[i] = x_u[i];
                        x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                    }
               
                    obj_b = obj_c;
                    obj_c = obj_u;
                    obj_u = getObjValueTest(x_u);


                }
            }
            else if ((s_u -s_ulim)*(s_ulim-s_c) >=0.0) // Limit parabolic u 
            {                                        // to maximum allowed value
                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_ulim[i];
                }
                obj_u = getObjValueTest(x_u);
            }
            else                                       // Reject parabolic u,
            {                                          // use default magnification
                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }

                obj_u = getObjValueTest(x_u);
            }

            //std::cout << "shift x positions " << std::endl;
            for (int i = 0; i < nVars; i++)
            {
                x_a[i] = x_b[i];         
                x_b[i] = x_c[i];
                x_c[i] = x_u[i];         
            }
            
            obj_a = obj_b;   
            obj_b = obj_c;
            obj_c = obj_u;
            
        }
        
        //std::cout << "final stage in Bracket 3 " << std::endl;
        //std::cout << "obj_a =" << obj_a << std::endl;
        //std::cout << "obj_b =" << obj_b << std::endl;
        //std::cout << "obj_c =" << obj_c << std::endl;
        
        /*
        delete [] x_r;
        x_r =0;

        delete [] x_q;
        x_q =0;

        delete [] x_u;
        x_u =0;

        delete [] x_ulim;
        x_ulim =0;

        delete [] x_temp;
        x_temp =0;
        */
    }
   
    void LBFGS::Bracket(std::vector<LIBMOL::REAL>  & x_a, 
                        std::vector<LIBMOL::REAL>  & x_b, 
                        std::vector<LIBMOL::REAL>  & x_c, 
                        LIBMOL::REAL& obj_a, 
                        LIBMOL::REAL& obj_b, 
                        LIBMOL::REAL& obj_c,
                        std::vector<LIBMOL::AtomDict>    & tAtoms, 
                        std::vector<LIBMOL::BondDict>    & tBonds, 
                        std::vector<LIBMOL::AngleDict>   & tAngs, 
                        std::vector<LIBMOL::TorsionDict> & tTors, 
                        std::vector<LIBMOL::RingDict>    & tRings, 
                        std::vector<LIBMOL::PlaneDict>   & tPlas, 
                        std::vector<LIBMOL::ChiralDict>  & tChs)
    
    {
        
       
        LIBMOL::REAL glim = 100.0;

        LIBMOL::REAL obj_u, obj_temp;

        LIBMOL::REAL s_b, s_c, s_u, s_ulim;
  
        LIBMOL::REAL tm1;
        
        std::vector<LIBMOL::REAL> x_r, x_q, x_u, x_ulim, x_temp;
        
        for (int i=0; i < nVars; i++)
        {
            x_r.push_back(0.0);
            x_q.push_back(0.0);
            x_u.push_back(0.0);
            x_ulim.push_back(0.0);
            x_temp.push_back(0.0);
        }
        
        FF::GetObjValue toolGetObjValue;
        toolGetObjValue.workSpace = workSpace;
        toolGetObjValue.lComp     = lComp;

        // Set an initial region
        
        obj_a = toolGetObjValue.getAll(tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        //obj_a = getObjValueTest(x_a);
        
        
        obj_b =  toolGetObjValue.getAll(x_b, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        
        
        // obj_b = getObjValueTest(x_b);
        // std::cout << "obj_a " << obj_a << " and obj_b " << obj_b << std::endl;
        
        if(obj_b > obj_a)
        {
            for (int i =0; i < nVars; i++)
            {
                x_temp[i] = x_a[i];
                x_a[i]    = x_b[i];
                x_b[i]    = x_temp[i];
            }

            obj_temp = obj_a;
            obj_a    = obj_b;
            obj_b    = obj_temp;

        }
           

        // First guess for c
        
        
        for (int i = 0; i < nVars; i++)
        {
            x_c[i] =  x_b[i]+(1+GOLD)*(x_b[i] - x_a[i]);
            // std::cout << "x_c " << i << " is " << x_c[i] << std::endl;
            
        }
        
        obj_c =  toolGetObjValue.getAll(x_c, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        // obj_c = getObjValueTest(x_c);
        //std::cout << "input stage in Bracket" << std::endl;
        //std::cout << "obj_a =" << obj_a << std::endl;
        //std::cout << "obj_b =" << obj_b << std::endl;
        //std::cout << "obj_c =" << obj_c << std::endl;
      
        
  
        // Loop over until find lower obj_c
        // int kk =0;

        // Loop over until find lower obj_c
        // int kk =0;

        while (obj_b > obj_c)
        {
            // kk++;
            // cout << "kk = " << kk << endl;

            for (int i =0; i <  nVars; i++)
            {
                x_r[i] = (x_b[i] - x_a[i])*(obj_b-obj_c);
                x_q[i] = (x_b[i] - x_c[i])*(obj_b-obj_a);
                tm1    = x_q[i] - x_r[i];

                if(fabs(tm1) < 1.0e-6)
                { 
                    // cout << "tm1 = " << tm1 << endl;
                    tm1 = 1.0e-6;
                }
                
                x_u[i] = x_b[i] - ((x_b[i]-x_c[i])*x_q[i]
                        -(x_b[i]-x_a[i])*x_r[i])/(2.0*tm1);

                //    cout << "x_u["<<i<<"] = " << x_u[i] << endl;   
         
                x_ulim[i] = x_b[i] - glim*(x_c[i]-x_b[i]);

            }
            
           
            s_b =0.0;
            s_c =0.0; 
            s_u =0.0;
            s_ulim = 0.0;
            
            
            for (int i =0; i < nVars; i++)
            {
                s_b += pow((x_b[i]-x_a[i]),2.0);
                s_c += pow((x_c[i]-x_a[i]),2.0);
                s_u += pow((x_u[i]-x_a[i]),2.0);
                s_ulim += pow((x_ulim[i]-x_a[i]),2.0);
            }
            s_b = sqrt(s_b);
            s_c = sqrt(s_c);
            s_u = sqrt(s_u);
            s_ulim = sqrt(s_ulim);
      
            
            
            std::vector<LIBMOL::REAL> x_bu, x_uc, x_uulim;
            for (int i=0; i < nVars; i++)
            {
                x_bu.push_back(x_b[i]-x_u[i]);
                x_uc.push_back(x_u[i]-x_c[i]);
                x_uulim.push_back(x_u[i]-x_ulim[i]);
            }
            
            // check the possible minimization region
 
            //if((s_b - s_u)*(s_u-s_c) > 0.0 ) // paraabolic u is between b and c
            if (LIBMOL::dotP(x_bu, x_uc) > 0.0)
            {
                // obj_u =  getObjValueTest(x_u);
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                if(obj_u < obj_c)           // get a minimum between b and c
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_a[i] = x_b[i];
                        x_b[i] = x_u[i];
                    }
              
                    obj_a = obj_b;
                    obj_b = obj_u;

                    //std::cout << "final stage in Bracket 1" << std::endl;
                    //std::cout << "obj_a =" << obj_a << std::endl;
                    //std::cout << "obj_b =" << obj_b << std::endl;
                    //std::cout << "obj_c =" << obj_c << std::endl;
       
                    return;

                }
                else if ( obj_u > obj_b)   // Got a minimum between a and u
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_c[i] = x_u[i];
                    }
              
                    obj_c = obj_u;

                    //std::cout << "final stage in Bracket " << std::endl;
                    //std::cout << "obj_a =" << obj_a << std::endl;
                    //std::cout << "obj_b =" << obj_b << std::endl;
                    //std::cout << "obj_c =" << obj_c << std::endl;
                    return;
                }
   

                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                // obj_u = getObjValueTest(x_u);

            }
            //else if ((s_c - s_u)*(s_u -s_ulim) > 0.0) // parabolic is between x_c
            else if (-LIBMOL::dotP(x_uc, x_uulim) >0.0)
            {                                    // and its allowed limit
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                // obj_u = getObjValueTest(x_u);
                
                if(obj_u < obj_c)
                {
                    for (int i = 0; i < nVars; i++)
                    {
                        x_b[i] = x_c[i];
                        x_c[i] = x_u[i];
                        x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                    }
               
                    obj_b = obj_c;
                    obj_c = obj_u;
                    obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                    obj_u = getObjValueTest(x_u);


                }
            }
            else if ((s_u -s_ulim)*(s_ulim-s_c) >=0.0) // Limit parabolic u 
            {                                        // to maximum allowed value
                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_ulim[i];
                }
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                // obj_u = getObjValueTest(x_u);
            }
            else                                       // Reject parabolic u,
            {                                          // use default magnification
                for (int i = 0; i < nVars; i++)
                {
                    x_u[i] = x_c[i]+(1+GOLD)*(x_c[i] -x_b[i]);
                }
                obj_u =  toolGetObjValue.getAll(x_u, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                // obj_u = getObjValueTest(x_u);
            }

            //std::cout << "shift x positions " << std::endl;
            for (int i = 0; i < nVars; i++)
            {
                x_a[i] = x_b[i];         
                x_b[i] = x_c[i];
                x_c[i] = x_u[i];         
            }
            
            obj_a = obj_b;   
            obj_b = obj_c;
            obj_c = obj_u;
            
        }
        
        //std::cout << "final stage in Bracket 3 " << std::endl;
        //std::cout << "obj_a =" << obj_a << std::endl;
        //std::cout << "obj_b =" << obj_b << std::endl;
        //std::cout << "obj_c =" << obj_c << std::endl;
        
        /*
        delete [] x_r;
        x_r =0;

        delete [] x_q;
        x_q =0;

        delete [] x_u;
        x_u =0;

        delete [] x_ulim;
        x_ulim =0;

        delete [] x_temp;
        x_temp =0;
        */
    }
   
    // Find a lower value of the objective function within
    // the bracket range
    void LBFGS::GSection(std::vector<LIBMOL::REAL> & x_a, 
                         std::vector<LIBMOL::REAL> & x_b, 
                         std::vector<LIBMOL::REAL> & x_c, 
                         LIBMOL::REAL & obj_a, 
                         LIBMOL::REAL & obj_b, 
                         LIBMOL::REAL & obj_c,
                         LIBMOL::REAL & tObj,
                         std::vector<LIBMOL::REAL> & tFinCoords,
                         LIBMOL::REAL (LBFGS::*func)(std::vector<LIBMOL::REAL> &))
    {
        
       
        LIBMOL::REAL tol = 1.0e-5;
        
        // LIBMOL::TransCoords transTool;
        
        // FF::GetObjValue toolGetObjValue;
        // toolGetObjValue.workSpace = workSpace;
        

        std::vector<LIBMOL::REAL>  x0;
        std::vector<LIBMOL::REAL>  x1;
        std::vector<LIBMOL::REAL>  x2;
        std::vector<LIBMOL::REAL>  x3;

        LIBMOL::REAL obj1, obj2;

        LIBMOL::REAL s_b, s_c, s_x1, s_x2, s_x3;
    
        /* temp
        if(workSpace == 1)
        {
            numVars = (int)tAtoms.size()*(int)tAtoms[0].coords.size();
        }
        else if (workSpace == 2)
        {
            numVars = (int)tTors.size();
        }
        else
        {
            std::cout << "In which space to do Flash-freezing ? " << std::endl;
            exit(1);
        }

        */
        
        for (int i=0; i < nVars; i++)
        {
            x0.push_back(0.0);
            x1.push_back(0.0);
            x2.push_back(0.0);
            x3.push_back(0.0);
        }
        
        s_b =0.0;
        s_c =0.0;

        for(int i =0; i < nVars; i++)
        {
            x0[i] = x_a[i];
            x3[i] = x_c[i];
            s_b += pow((x_b[i]-x_a[i]),2.0);
            s_c += pow((x_c[i]-x_b[i]),2.0);
        }
     
        s_b = sqrt(s_b);
        s_c = sqrt(s_c);
   
        if(s_c > s_b) 
        {
            for(int i =0; i < nVars; i++)
            {
                x1[i] = x_b[i];
                x2[i] = x_b[i]+(1-GOLD)*(x_c[i]-x_b[i]);
            }
        }
        else
        {
            for( int i =0; i < nVars; i++)
            {
                x2[i] = x_b[i];
                x1[i] = x_b[i]-(1-GOLD)*(x_b[i]-x_a[i]);
            }
        }

        //obj1 = getObjValueTest(x1);
        //obj2 = getObjValueTest(x2);
        obj1 = (this->*func)(x1);
        obj2 = (this->*func)(x2);

        //std::cout << "obj1 in section " << obj1 << std::endl;
        //std::cout << "obj2 in section " << obj2 << std::endl;
       
        s_x1 = 0.0;
        s_x2 = 0.0;
        s_x3 = 0.0;

        /*
        for(int i =0; i < nVars; i++)
        {
            s_x1+=pow((x1[i]-x0[i]),2.0);
            s_x2+=pow((x2[i]-x0[i]),2.0);
            s_x3+=pow((x3[i]-x0[i]),2.0);
        }

        s_x1 = sqrt(s_x1);
        s_x2 = sqrt(s_x2);
        s_x3 = sqrt(s_x3);
        */
        s_x1 = LIBMOL::lengthV(x1);
        s_x2 = LIBMOL::lengthV(x2);
        s_x3 = LIBMOL::distanceV(x3, x0);
        
        //  cout << "Init sections are  s_x1 " << s_x1 << endl
        //     << " s_x2  " << s_x2 << endl
        //     << " s_x3  " << s_x3 << endl;

        int j =0;
        
        while ((s_x3 > tol*(s_x1+s_x2)) && j < 20)
        {
            if(obj2 < obj1)
            {
                //  cout << "obj1 > obj2 " << endl;
                for(int i =0; i < nVars; i++)
                {
                    x0[i] = x1[i];
                    x1[i] = x2[i];
                    x2[i] = GOLD*x1[i]+(1-GOLD)*x3[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
	            //   cout << "x2["<<i<<"]= "<< x2[i] << endl;             
                }
           
                obj1 = obj2;
                // obj2 = getObjValueTest(x2);
                obj2 = (this->*func)(x2);
	        //	  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;
            }
            else
            {
                //  cout << "obj2 < obj1 " << endl;
                for(int i =0; i < nVars; i++)
                {
                    x3[i] = x2[i];
                    x2[i] = x1[i];
                    x1[i] = GOLD*x2[i]+(1-GOLD)*x0[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
                    //  cout << "x2["<<i<<"]= "<< x2[i] << endl;  
                }
           
                obj2 = obj1;
                //obj1 = getObjValueTest(x1);
                obj1 = (this->*func)(x1);
	        //  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;       
            }

            s_x1 = LIBMOL::lengthV(x1);
            s_x2 = LIBMOL::lengthV(x2);
            s_x3 = LIBMOL::distanceV(x3, x0);
            j++;
            // cout << "j = " << j << endl;
        }

        //std::cout << "new obj1 is " << obj1 << std::endl;
        //std::cout << "new obj2 is " << obj2 << std::endl;

        // LIBMOL::REAL obj1_t, obj2_t;
        /*
        if(workSpace == 1 )
        {
            if(obj1 < obj2)
            {
               //obj1_t = toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
               curObjValue = obj1;
               UpdateAtomCoords(tAtoms, x1);
               
	       // cout << "Obj1_t = " << obj1_t << endl;
               // cout << "myObj  = " << myObjValue << endl;
            }
            else 
            {
                //obj2_t = toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
                curObjValue = obj2;
                UpdateAtomCoords(tAtoms, x2);
	        //      cout << "Obj2_t = " << obj2_t << endl;
                //      cout << "myObj  = " << myObjValue << endl;
            }
        }
        else if(workSpace == 2)
        {
            if(obj1 < obj2)
            { 
                for(i = 0; i < (int)tTors.size(); i++)
                {
                    tTors[i].value = x1[i];
                }      
                // UpdatTors_t();  
                //CoordAngToCart();
                curObjValue = obj1;
            }
            else
            {
                for(i = 0; i < (int)tTors.size(); i++)
                {
                    tTors[i].value = x2[i];
                }  
         
                // UpdatTors_t();  
                curObjValue = obj2;
            }
            transTool.generateCoordTorsToCart(tAtoms, tBonds, tAngs, tTors, 
                                              tRings, tPlas, tChs);  
        }             
        */
        
        tFinCoords.clear();
        if(obj1 < obj2)
        {
            tObj = obj1;
            for (int j=0; j < (int)x1.size(); j++)
            {
                tFinCoords.push_back(x1[j]);
            }
            
        }
        else
        {
            tObj = obj2;  
            for (int j=0; j < (int)x1.size(); j++)
            {
                tFinCoords.push_back(x2[j]);
            }
        }
        
        std::cout << "my obj after line search " <<tObj << std::endl;
       
        //  std::cout << "Continue " << endl;
        //  std::cin.get();
        /*
        delete [] x0;
        x0 =0;
        delete [] x1;
        x1 =0;
        delete x2;
        x2 =0;
        delete x3;
        x3 =0;
        */
    } 
    
    void LBFGS::GSection(std::vector<LIBMOL::REAL> & x_a, 
                         std::vector<LIBMOL::REAL> & x_b, 
                         std::vector<LIBMOL::REAL> & x_c, 
                         LIBMOL::REAL & obj_a, 
                         LIBMOL::REAL & obj_b, 
                         LIBMOL::REAL & obj_c,
                         LIBMOL::REAL & tObj,
                         std::vector<LIBMOL::REAL>         & tFinCoords,
                         std::vector<LIBMOL::AtomDict>     & tAtoms, 
                         std::vector<LIBMOL::BondDict>     & tBonds, 
                         std::vector<LIBMOL::AngleDict>    & tAngs, 
                         std::vector<LIBMOL::TorsionDict>  & tTors, 
                         std::vector<LIBMOL::RingDict>     & tRings, 
                         std::vector<LIBMOL::PlaneDict>    & tPlas, 
                         std::vector<LIBMOL::ChiralDict>   & tChs)
    {
        
        LIBMOL::REAL tol = 1.0e-5;
        
        // LIBMOL::TransCoords transTool;
        
        FF::GetObjValue toolGetObjValue;
        toolGetObjValue.workSpace = workSpace;
        toolGetObjValue.lComp     = lComp;

        std::vector<LIBMOL::REAL>  x0;
        std::vector<LIBMOL::REAL>  x1;
        std::vector<LIBMOL::REAL>  x2;
        std::vector<LIBMOL::REAL>  x3;

        LIBMOL::REAL obj1, obj2;

        LIBMOL::REAL s_b, s_c, s_x1, s_x2, s_x3;
        
        for (int i=0; i < nVars; i++)
        {
            x0.push_back(0.0);
            x1.push_back(0.0);
            x2.push_back(0.0);
            x3.push_back(0.0);
        }
        
        s_b =0.0;
        s_c =0.0;

        for(int i =0; i < nVars; i++)
        {
            x0[i] = x_a[i];
            x3[i] = x_c[i];
            s_b += pow((x_b[i]-x_a[i]),2.0);
            s_c += pow((x_c[i]-x_b[i]),2.0);
        }
     
        s_b = sqrt(s_b);
        s_c = sqrt(s_c);
   
        if(s_c > s_b) 
        {
            for(int i =0; i < nVars; i++)
            {
                x1[i] = x_b[i];
                x2[i] = x_b[i]+(1-GOLD)*(x_c[i]-x_b[i]);
            }
        }
        else
        {
            for( int i =0; i < nVars; i++)
            {
                x2[i] = x_b[i];
                x1[i] = x_b[i]-(1-GOLD)*(x_b[i]-x_a[i]);
            }
        }

        //obj1 = getObjValueTest(x1);
        //obj2 = getObjValueTest(x2);
        //obj1 = (this->*func)(x1);
        //obj2 = (this->*func)(x2);
        obj1 =  toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        obj2 =  toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
        //std::cout << "obj1 in section " << obj1 << std::endl;
        //std::cout << "obj2 in section " << obj2 << std::endl;
       
        s_x1 = 0.0;
        s_x2 = 0.0;
        s_x3 = 0.0;

        /*
        for(int i =0; i < nVars; i++)
        {
            s_x1+=pow((x1[i]-x0[i]),2.0);
            s_x2+=pow((x2[i]-x0[i]),2.0);
            s_x3+=pow((x3[i]-x0[i]),2.0);
        }

        s_x1 = sqrt(s_x1);
        s_x2 = sqrt(s_x2);
        s_x3 = sqrt(s_x3);
        */
        s_x1 = LIBMOL::lengthV(x1);
        s_x2 = LIBMOL::lengthV(x2);
        s_x3 = LIBMOL::distanceV(x3, x0);
        
        //  cout << "Init sections are  s_x1 " << s_x1 << endl
        //     << " s_x2  " << s_x2 << endl
        //     << " s_x3  " << s_x3 << endl;

        int j =0;
        
        while ((s_x3 > tol*(s_x1+s_x2)) && j < 20)
        {
            if(obj2 < obj1)
            {
                //  cout << "obj1 > obj2 " << endl;
                for(int i =0; i < nVars; i++)
                {
                    x0[i] = x1[i];
                    x1[i] = x2[i];
                    x2[i] = GOLD*x1[i]+(1-GOLD)*x3[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
	            //   cout << "x2["<<i<<"]= "<< x2[i] << endl;             
                }
           
                obj1 = obj2;
                // obj2 = getObjValueTest(x2);
                // obj2 = (this->*func)(x2);
                obj2 =  toolGetObjValue.getAll(x2, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
	        //	  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;
            }
            else
            {
                //  cout << "obj2 < obj1 " << endl;
                for(int i =0; i < nVars; i++)
                {
                    x3[i] = x2[i];
                    x2[i] = x1[i];
                    x1[i] = GOLD*x2[i]+(1-GOLD)*x0[i];
	            //  cout << "x1["<<i<<"]= "<< x1[i] << endl;
                    //  cout << "x2["<<i<<"]= "<< x2[i] << endl;  
                }
           
                obj2 = obj1;
                //obj1 = getObjValueTest(x1);
                //obj1 = (this->*func)(x1);
                obj1 =  toolGetObjValue.getAll(x1, tAtoms, tBonds, tAngs, tTors, tPlas, tChs);
	        //  cout << "obj1 " << obj1 << "    obj2 " << obj2 << endl;       
            }

            s_x1 = LIBMOL::lengthV(x1);
            s_x2 = LIBMOL::lengthV(x2);
            s_x3 = LIBMOL::distanceV(x3, x0);
            j++;
            // cout << "j = " << j << endl;
        }

        // std::cout << "new obj1 is " << obj1 << std::endl;
        // std::cout << "new obj2 is " << obj2 << std::endl;

        tFinCoords.clear();
        if(obj1 < obj2)
        {
            tObj = obj1;
            for (int j=0; j < (int)x1.size(); j++)
            {
                tFinCoords.push_back(x1[j]);
            }
            
        }
        else
        {
            tObj = obj2;  
            for (int j=0; j < (int)x1.size(); j++)
            {
                tFinCoords.push_back(x2[j]);
            }
        }
        
        // std::cout << "my obj after line search " <<tObj << std::endl;
       
        
    } 
    
    
    bool LBFGS::checkConvergence(std::vector<LIBMOL::REAL> & tDerivs)
    {
        bool tC=false;
        
        if (LIBMOL::lengthV(tDerivs) <= rTol)
        {
            tC=true;
        }
        
        return tC;
    }
    
    bool LBFGS::checkConvergence(std::vector<LIBMOL::REAL>& tDerivs, 
                                 std::vector<LIBMOL::REAL>& tNewCoords, 
                                 std::vector<LIBMOL::REAL>& tOldCoords)
    {
        bool tC=false;
        
        LIBMOL::REAL oL = LIBMOL::lengthV(tNewCoords);
        LIBMOL::REAL nL = LIBMOL::lengthV(tOldCoords);
        
        if (LIBMOL::lengthV(tDerivs) <= rTol)
        {
            tC=true;
        }
        else if (oL >rTol)
        {
            if (fabs((nL-oL)) < 1.0e-9)
            {
                tC=true;
            }
        }
        
        
        
        
        return tC;
    }
    void LBFGS::updateCoords(std::vector<LIBMOL::AtomDict>& tAtoms, 
                             std::vector<LIBMOL::REAL> &    tCoords)
    {
        int tSize = (int)tCoords.size();
        for (int i=0; i < (int)tAtoms.size(); i++)
        {
            int tDim = (int)tAtoms[i].coords.size();
            for (int j=0; j < tDim; j++)
            {
                int k= i*tDim+j;
                if (k < tSize)
                {
                    tAtoms[i].coords[j] = tCoords[k];
                }
            }
        }
    }
    
    void LBFGS::updateCoordsAndDerivs(std::vector<LIBMOL::REAL>& tNewCoords, 
                                      std::vector<LIBMOL::REAL>& tOldCoords, 
                                      std::vector<LIBMOL::REAL>& tNewDerivs, 
                                      std::vector<LIBMOL::REAL>& tOldDerivs)
    {
        
        for (int i=0; i < (int)tNewCoords.size(); i++)
        {
            tOldCoords[i]= tNewCoords[i];
            tOldDerivs[i]= tNewDerivs[i];
        }
    }
    
    void LBFGS::getFirstDerivs(std::vector<LIBMOL::REAL>        & tDerivs, 
                               std::vector<LIBMOL::REAL>        & tCoords, 
                               std::vector<LIBMOL::AtomDict>    & tAtoms, 
                               std::vector<LIBMOL::BondDict>    & tBonds, 
                               std::vector<LIBMOL::AngleDict>   & tAngs, 
                               std::vector<LIBMOL::TorsionDict> & tTors, 
                               std::vector<LIBMOL::RingDict>    & tRings, 
                               std::vector<LIBMOL::PlaneDict>   & tPlas, 
                               std::vector<LIBMOL::ChiralDict>  & tChs)
    {
        int numVars = (int)tCoords.size();
        FF::GetAllDerivsFull   derivs(numVars);
        derivs.workMode  = 1;
        derivs.workSpace = workSpace;
        updateCoords(tAtoms, tCoords);
        derivs.setFirstDerivsCart(tAtoms, tBonds, tAngs,
                                  tTors,  tPlas, tChs);
        
        tDerivs.clear();
        for (int i=0; i < numVars; i++)
        {
            tDerivs.push_back(derivs.firDrivCart[i]);
        }
    }
    
    void LBFGS::getFirstDerivs(std::vector<LIBMOL::REAL>        & tDerivs, 
                               std::vector<LIBMOL::REAL>        & tCoords, 
                               std::vector<LIBMOL::AtomDict>    & tAtoms, 
                               std::vector<LIBMOL::BondDict>    & tBonds, 
                               std::vector<LIBMOL::AngleDict>   & tAngs, 
                               std::vector<LIBMOL::TorsionDict> & tTors, 
                               std::vector<LIBMOL::RingDict>    & tRings, 
                               std::vector<LIBMOL::PlaneDict>   & tPlas, 
                               std::vector<LIBMOL::ChiralDict>  & tChs,
                               std::vector<LIBMOL::AtomDict>    & tAllAtoms)
    {
        int numVars = (int)tCoords.size();
        FF::GetAllDerivsFull   derivs(numVars);
        derivs.workMode  = 1;
        derivs.workSpace = workSpace;
        updateCoords(tAtoms, tCoords);
        derivs.setFirstDerivsCart(tAtoms, tBonds, tAngs,
                                  tTors,  tPlas, tChs, tAllAtoms);
        
        tDerivs.clear();
        for (int i=0; i < numVars; i++)
        {
            tDerivs.push_back(derivs.firDrivCart[i]);
        }
    }
    
    LIBMOL::REAL LBFGS::getObjValueTest(std::vector<LIBMOL::REAL> & tCoords)
    {
        LIBMOL::REAL aObjValue=0.0;
        
        for (int i=0; i < (int)tCoords.size(); i+=2)
        {
            LIBMOL::REAL T1 = 1.0-tCoords[i];
            LIBMOL::REAL T2 = 10.0*std::pow((tCoords[i+1]-tCoords[i]), 2.0);
            aObjValue =   aObjValue + std::pow(T1, 2.0)+ std::pow(T2, 2.0);  
        }

        return aObjValue;
        
    }
    
    void LBFGS::getDerivsTest(std::vector<LIBMOL::REAL>& tCoords, 
                              std::vector<LIBMOL::REAL>& tDerivs)
    {
        tDerivs.clear();
        for (int i=0; i < (int)tCoords.size(); i++ )
        {
            tDerivs.push_back(0.0);
        }
        
        for (int i=0; i < (int)tCoords.size(); i+=2)
        {   
            LIBMOL::REAL T1 = 1.0-tCoords[i];
            LIBMOL::REAL T2 = 10.0*std::pow((tCoords[i+1]-tCoords[i]), 2.0);
            tDerivs[i+1] = 20.0*T2;
            tDerivs[i]   = -2.0*(tCoords[i]*tDerivs[i+1]+T1);
        }
        
        LIBMOL::normalizeV(tDerivs);
    }
    
    void LBFGS::testCase()
    {
        LIBMOL::REAL finObj=0.0;
        std::vector<LIBMOL::REAL> iniCrds, finCrds, iniDrs;
        int N = 100;
        for (int i=0; i <N; i+=2)
        {
            iniCrds.push_back(-1.2);
            iniCrds.push_back(1.0);
            iniDrs.push_back(0.0);
            iniDrs.push_back(0.0);
        }
        
        
        LBFGS aL_BFGSer(100, 5, 1.0e-8, 1.0e-5, 1.0e-8, 1.0);
        
        aL_BFGSer.getDerivsTest(iniCrds, iniDrs);
        
        aL_BFGSer.driver(finObj, finCrds, iniCrds, iniDrs);
        
        std::cout << "Finally, the value of objective function is  "
                  << finObj << std::endl;
        
        std::cout << "Its associated coords are: " << std::endl;
        
        for (int i=0; i < (int)finCrds.size(); i++)
        {
            std::cout << finCrds[i] << " for coordinate " << i << std::endl;
        }
        
        exit(1);
    }
    
}
