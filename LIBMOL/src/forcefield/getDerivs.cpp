/* 
 * File:   getDerivs.cpp
 * Author: flong
 *
 * last updated on February 18, 2013, 4:30 PM
 */

#include "getDerivs.h"
#include "getObjValues.h"

namespace FF
{
    GetAllDerivsFull::GetAllDerivsFull(int tNumVars):workMode(1),
                                                     workSpace(1)
    {
        numVars = tNumVars;
        if (numVars)
        {
            firDrivCart = new LIBMOL::REAL [numVars];
            
            secDrivCart = new LIBMOL::REAL *  [numVars];
            for (int i=0; i < numVars; i++)
            { 
                firDrivCart[i] = 0.0;
                secDrivCart[i] = new LIBMOL::REAL [numVars];
                for (int j=0; j < numVars; j++)
                {
                    secDrivCart[i][j] =0.0;
                }          
            }
        }
    }
    
    GetAllDerivsFull::~GetAllDerivsFull()
    {
        if(numVars)
        {
            if (firDrivCart !=NULL)
            {
                delete [] firDrivCart;
                firDrivCart = 0;
            }
            
           
           
            for (int i=0; i < numVars; i++)
            {
                delete [] secDrivCart[i];
                secDrivCart[i] =0;
            }
            delete [] secDrivCart;
            secDrivCart = 0;
          
           
        }
    }
    
    void  GetAllDerivsFull::setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts, 
                                std::vector<LIBMOL::BondDict>& tBos, 
                                std::vector<LIBMOL::AngleDict>& tAns, 
                                std::vector<LIBMOL::TorsionDict>& tTos, 
                                std::vector<LIBMOL::PlaneDict>& tPls, 
                                std::vector<LIBMOL::ChiralDict>& tChs)
    {
        
        int nAtoms = (int)tAts.size();
        FF::GetObjValue  toolObjs;
        toolObjs.workSpace = workSpace;
        
        if(nAtoms != 0)
        {
            int dim    = tAts[0].coords.size();
            if (nAtoms*dim == numVars)
            {
                if (workMode ==1)
                {
                    if (workSpace ==1)
                    {
                        // get object values and forces (first derivs on atoms
                        toolObjs.getAll(tAts, tBos, tAns, tTos, tPls, tChs);
                        
                        setFirstDerivsCart(tAts);
                    }
                    
                }
            }
            
            
        }
    }
    
    void  GetAllDerivsFull::setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts, 
                                std::vector<LIBMOL::BondDict>    & tBos, 
                                std::vector<LIBMOL::AngleDict>   & tAns, 
                                std::vector<LIBMOL::TorsionDict> & tTos, 
                                std::vector<LIBMOL::PlaneDict>   & tPls, 
                                std::vector<LIBMOL::ChiralDict>  & tChs,
                                std::vector<LIBMOL::AtomDict>    & tAllAtoms)
    {
        
        int nAtoms = (int)tAts.size();
        FF::GetObjValue  toolObjs;
        toolObjs.workSpace = workSpace;
       
        
        if(nAtoms != 0)
        {
            int dim    = tAts[0].coords.size();
            if (nAtoms*dim == numVars)
            {
                if (workMode ==1)
                {
                    if (workSpace ==1)
                    {
                        // get object values and forces (first derivs on atoms
                        toolObjs.getAll(tAts, tBos, tAns, tTos, tPls, tChs, tAllAtoms);
                        std::cout << "in dire " << std::endl;
                        exit(1);
                        
                        setFirstDerivsCart(tAts);
                    }
                    
                }
            }
            
            
        }
    }
    
    void GetAllDerivsFull::setFirstDerivsCart(std::vector<LIBMOL::AtomDict>& tAts)
    {
        LIBMOL::REAL maxF=0.0;
        int tSize = (int)tAts.size();
        // std::cout << "numVars " << numVars << std::endl;
 
        for (int i=0; i < tSize; i++ )
        {
            int tDim = (int)tAts[i].forces.size();
            for (int j=0; j < tDim; j++)
            {
                int k = i*tDim+j;
                
                firDrivCart[k] =-tAts[i].forces[j];
                if (fabs(firDrivCart[k]) > maxF)
                {
                    maxF = fabs(firDrivCart[k]);
                }
            }
        }
        
        
        if (maxF >1.0e-8)
        {
            for (int i=0; i < tSize; i++ )
            {
                int tDim = (int)tAts[i].forces.size();
                for (int j=0; j < tDim; j++)
                {
                    int k = i*tDim+j;
                    firDrivCart[k] =firDrivCart[k]/maxF;
                }
            }
        }
    }
    
    
    void GetAllDerivsFull::setFirstDerivsCartBond(std::vector<LIBMOL::AtomDict>& tAts,
                                                  std::vector<LIBMOL::BondDict>::iterator tBo)
    {
        
            int dim = (int)tAts[0].coords.size();
            LIBMOL::REAL * distComp = new LIBMOL::REAL [dim];
            
            LIBMOL::REAL bondL = 0.0;
            for (int i=0; i < dim; i++)
            {
                distComp[i] = tAts[tBo->atomsIdx[1]].coords[i]-tAts[tBo->atomsIdx[0]].coords[i];
                bondL      += pow(distComp[i], 2.0);
            }
            bondL = sqrt(bondL);
            
            if (bondL)
            {
                for (int i=0; i < dim; i++)
                {
                    LIBMOL::REAL tFor = (bondL-tBo->valueST)*distComp[i]/(bondL*tBo->sigValueST);
                    int idx_1 = tBo->atomsIdx[1]*dim + i;
                    int idx_2 = tBo->atomsIdx[0]*dim + i;
                    firDrivCart[idx_1] -=tFor;
                    firDrivCart[idx_2] +=tFor;        
                }
            }
            
            delete [] distComp;
            distComp = 0;
    }
    
    void GetAllDerivsFull::setFirstDerivsCartAngle(std::vector<LIBMOL::AtomDict>& tAts, 
                                                   std::vector<LIBMOL::AngleDict>::iterator tAn)
    {
        
        int dim = (int)tAts[0].coords.size();
        int  j,k;
  
        int i_atm = tAn->atoms[0], j_atm = tAn->atoms[1], k_atm = tAn->atoms[2];
        
        LIBMOL::REAL R, angs;
    
        LIBMOL::REAL *    dThetadx_t = new LIBMOL::REAL [3*dim];

        LIBMOL::REAL *a = new LIBMOL::REAL [dim];
        LIBMOL::REAL *b = new LIBMOL::REAL [dim];

        LIBMOL::REAL ** da_dx = new LIBMOL::REAL * [3*dim];
        for(j =0; j < 3*dim; j++)
        {
            da_dx [j] = new LIBMOL::REAL [dim];
        }
        LIBMOL::REAL ** db_dx = new LIBMOL::REAL * [3*dim];
        for(j =0; j < 3*dim; j++)
        {
            db_dx [j] = new LIBMOL::REAL [dim];
        }
        
        LIBMOL::REAL  leng_a, leng_b;

        LIBMOL::REAL * d_leng_a = new LIBMOL::REAL [3*dim];
        LIBMOL::REAL * d_leng_b = new LIBMOL::REAL [3*dim];
        
        for (j =0; j < dim; j++)
        {
            a[j] = tAts[k_atm].coords[j]-tAts[i_atm].coords[j];
            b[j] = tAts[j_atm].coords[j]-tAts[i_atm].coords[j]; 
        }
        
        leng_a = LIBMOL::length_v(a, dim);
        leng_b = LIBMOL::length_v(b, dim);
        
        if(leng_a < 1.0e-8)
        {
            std::cout << "atom " << tAts[i_atm].id << " and atom " 
                      << tAts[k_atm].id 
                      << " overlapped ! Wrong " << std::endl;
            exit(1);
        }

        if(leng_b < 1.0e-8)
        {
            std::cout << "atom " << tAts[i_atm].id << " and atom " 
                      << tAts[j_atm].id << " overlapped ! Wrong " << std::endl;
            exit(1);
        }
        
        R    = LIBMOL::DotP(a,b)/(leng_a*leng_b);
        
        if(fabs(1-R*R) < 1.0e-7)
        {

            // cout << " Bond angle " << i+1 << endl;
            //  cout << "leng_a = " << leng_a << endl;
            //  cout << "leng_b = " << leng_b << endl;
            //  cout << "R      = " << R << endl;
   
            LIBMOL::REAL delta_ram = 1.0e-2;
            for (j =0; j < dim; j++)
            {
                tAts[j_atm].coords[j] 
                      = tAts[j_atm].coords[j] + delta_ram*LIBMOL::GetRand();
                b[j]  = tAts[j_atm].coords[j]-tAts[i_atm].coords[j];
            }
      
            leng_b = LIBMOL::length_v(b, dim);

            if( leng_b ==0)
            {
                std::cout << "vector a or b equals to zero in the bond angle : " 
                         << tAn->anchorID << std::endl;
                exit(0);
            }
      
            R = LIBMOL::DotP(a,b)/(leng_a*leng_b);
         
            //  std::cout << "new R      = " << R << std::endl;
        }
        
        angs = LIBMOL::GetAngle(a,b);
 
        //cout << "Instant value of the angle is " << angs*PID180 << endl
        //     << "Standard value of the angle is " 
        //     <<  hConstraBondAngs[i_ba].value_obs*PID180 << endl;
        //  cout << "Continue ?" << endl;
        //  cin.get();

        for (j =0; j < 3*dim; j++)
        {
            d_leng_a[j] = 0.0;
            d_leng_b[j] = 0.0;
            for (k =0; k < dim; k++)
            {
                da_dx[j][k] = 0.0;
                db_dx[j][k] = 0.0;
            }
        }

        for ( j =0; j < dim; j++)
        {
            da_dx[j][j]        =  1.0;
            da_dx[j+dim][j]    = -1.0;
            db_dx[j+dim][j]    = -1.0;
            db_dx[j+2*dim][j]  =  1.0;
        }

        for (j = 0; j < dim; j++)
        {
            d_leng_a[j]       =  a[j]/leng_a;
            d_leng_a[j+dim]   = -d_leng_a[j];
            d_leng_b[j+dim]   = -b[j]/leng_b;
            d_leng_b[j+2*dim] = -d_leng_b[j+dim];          
        }
      
        int i_comp, j_comp, k_comp;
        LIBMOL::REAL dR_dx;

        for (j =0; j < 3*dim; j++)
        {
            dR_dx = (LIBMOL::DotP(da_dx[j],b)
                    +LIBMOL::DotP(db_dx[j],a))/(leng_a*leng_b)
                    -LIBMOL::DotP(a,b)*(leng_b*d_leng_a[j]+leng_a*d_leng_b[j])
                    /pow(leng_a*leng_b,2.0);
     
            if(fabs(1-R*R) > 1.0e-8)
            {
                dThetadx_t[j] = -dR_dx/sqrt(1-R*R);
            }
            else
            {
                std::cout << "R = " << std::endl;
                std::cout << "angle calculation " << std::endl;
                std::cout << "Cos(theta) = 1, atoms overlapped ! " 
                          << std::endl;
                exit(1);
            } 
        }
        
        for (j =0; j < dim; j++) 
        {
            k_comp = j;
            i_comp = dim +j;
            j_comp = 2*dim + j;
            
            int idx_1 = i_atm*dim, idx_2=j_atm*dim, idx_3=k_atm*3;
            
            firDrivCart[idx_3] -=
                     ((angs - tAn->valueST)
                     *dThetadx_t[k_comp]/tAn->sigValueST);

            firDrivCart[idx_1] -=
                    ((angs - tAn->sigValueST)
                     *dThetadx_t[i_comp]/tAn->sigValueST);
            firDrivCart[idx_2] -=
                    ((angs - tAn->sigValueST)
                     *dThetadx_t[j_comp]/tAn->sigValueST);
       
        }
        
        
        delete [] dThetadx_t;
        dThetadx_t = 0;
  
        delete [] a;
        a = 0;

        delete [] b;
        b = 0;

        for(j =0; j < 3*dim; j++)
        {
            delete [] da_dx [j];
            da_dx [j] = 0;
        }
        delete [] da_dx;
        da_dx = 0;
  
       for(j =0; j < 3*dim; j++)
       {
           delete [] db_dx [j];
           db_dx[j] = 0;
       }
       delete [] db_dx;
       db_dx = 0;

       delete [] d_leng_a;
       d_leng_a = 0;

       delete [] d_leng_b;
       d_leng_b = 0;
      
    }
    
    void GetAllDerivsFull::setFirstDerivsCartPlane(std::vector<LIBMOL::AtomDict>& tAts,
                                                   std::vector<LIBMOL::PlaneDict>::iterator tPl, 
                                                   LIBMOL::REAL ** t1stDerivPlan)
    {
        int n_ato = (int)tPl->atoms.size();
         
        if (n_ato)
        {
            
            std::vector<int> atoIdxs;
            for (std::map<LIBMOL::ID, int>::iterator iAto=tPl->atoms.begin();
                    iAto != tPl->atoms.end(); iAto++)
            {
                atoIdxs.push_back(iAto->second);
            }
            
            int dim = (int)tAts[0].coords.size();
            // Calculate the average position of all atoms in the plane
            
            LIBMOL::REAL * X_ave = new LIBMOL::REAL [dim];
            
            for (int i1 = 0; i1 < dim; i1++)
            {           
                X_ave[i1] = 0.0;
                for (int i2 = 0; i2 < n_ato; i2++)
                {
                    X_ave [i1] += tAts[atoIdxs[i2]].coords[i1];
                }
                X_ave [i1] = X_ave [i1]/n_ato;
            }
            
            // Calculate the relative coords of atoms in the plane
            
            LIBMOL::REAL ** X_rel = new LIBMOL::REAL * [n_ato];
            for(int i1 = 0; i1 < n_ato; i1++)
            {
                X_rel [i1] = new LIBMOL::REAL [dim];
            }

            for (int i1 =0; i1 < n_ato; i1++)
            {
                for (int i2 = 0; i2 < dim; i2++)
                {
                    X_rel[i1][i2] = 
                       tAts[atoIdxs[i1]].coords[i2] - X_ave[i2];
                }
            }
            
            // Find the least square plane that fits all the atoms
            // and the derivatives of the plane coefficients rwt
            // atomic coordinates

            LIBMOL::REAL  * Coefs = new LIBMOL::REAL [dim];
            for (int i1 =0; i1 < dim; i1++)
            {
                Coefs[i1] = 0.0;
            }

            int         total_dim = n_ato*dim;

            LIBMOL::REAL ** d_Coefs = new LIBMOL::REAL * [total_dim];
            for (int i1 =0; i1 < total_dim; i1++)
            {
                d_Coefs[i1] = new LIBMOL::REAL [dim];
                for(int i2 = 0; i2 < dim; i2++)
                {
                    d_Coefs[i1][i2] = 0.0;
                }
            }
            
            LIBMOL::Plane_and_Deriv_Find(n_ato, X_rel, Coefs, d_Coefs);
            
            //realtype Coefs_D;
            //Coefs_D = DotP(X_ave, Coefs);
      
            LIBMOL::REAL * d_Coefs_D = new LIBMOL::REAL [total_dim];
            
            LIBMOL::REAL dxrdx = 1.0-(1.0/n_ato);

            int i_dim = 0;
            for (int i1 = 0; i1 < n_ato; i1++)
            {
                d_Coefs_D[i_dim] = d_Coefs[i_dim][0]*dxrdx*X_ave[0]
                            +Coefs[0]/n_ato 
                            +d_Coefs[i_dim][1]*dxrdx*X_ave[1]
                            +d_Coefs[i_dim][2]*dxrdx*X_ave[2]; 
                i_dim++;
                d_Coefs_D[i_dim] = d_Coefs[i_dim][0]*dxrdx*X_ave[0]
                            +d_Coefs[i_dim][1]*dxrdx*X_ave[1]
                            +Coefs[1]/n_ato 
                            +d_Coefs[i_dim][2]*dxrdx*X_ave[2];
                i_dim++;
                d_Coefs_D[i_dim] = d_Coefs[i_dim][0]*dxrdx*X_ave[0]
                            +d_Coefs[i_dim][1]*dxrdx*X_ave[1]
                            +d_Coefs[i_dim][2]*dxrdx*X_ave[2]
                            +Coefs[2]/n_ato;
                i_dim++;        
            }
            
            // Find the forces on each atoms in the plane

            LIBMOL::REAL t1_eq, t2_eq, t1;
            
            for(int i1 = 0; i1 < n_ato; i1++)
            {
                // int i_ato = atoIdxs[i1];
                t1_eq = 2.0*LIBMOL::DotP(Coefs,X_rel[i1]);

                // df/dx
         
                i_dim = 3*i1;
 
                t2_eq = 0.0;
                for(int i2 = 0; i2 < n_ato; i2++)
                {
                    if(i2 != i1)
                    {
                        t2_eq +=(2.0*LIBMOL::DotP(Coefs,X_rel[i2])
                                *LIBMOL::DotP(d_Coefs[i_dim],X_rel[i2]));
                    }
                }

                t1 =t1_eq*(LIBMOL::DotP(d_Coefs[i_dim],X_rel[i1])+Coefs[0])+t2_eq;
    
                
                t1stDerivPlan[i1][i_dim] -=
                          (tPl->fConst*t1*dxrdx);
                
                // df/dy

                i_dim++;

                t2_eq = 0.0;
                for(int i2 = 0; i2 < n_ato; i2++)
                {
                    if(i2 != i1)
                    {
                        t2_eq +=(2.0*LIBMOL::DotP(Coefs,X_rel[i2])
                           *LIBMOL::DotP(d_Coefs[i_dim],X_rel[i2]));
                    }
                }
           
                t1 =t1_eq*(LIBMOL::DotP(d_Coefs[i_dim],X_rel[i1])+Coefs[1])
                         +t2_eq;

                t1stDerivPlan[i1][i_dim] -=
                           (tPl->fConst*t1*dxrdx);
                
                // df/dz

                i_dim++;

                t2_eq = 0.0;
                for(int i2 = 0; i2 < n_ato; i2++)
                {
                    if(i2 != i1)
                    {
                        t2_eq +=(2.0*LIBMOL::DotP(Coefs,X_rel[i2])
                                *LIBMOL::DotP(d_Coefs[i_dim],X_rel[i2]));
                    }
                }

           
                t1 =t1_eq*(LIBMOL::DotP(d_Coefs[i_dim],X_rel[i1])+Coefs[2])
                         +t2_eq;
                
                t1stDerivPlan[i1][i_dim]-=
                           (tPl->fConst*t1*dxrdx);
            }
            
            // release memory 

            delete [] X_ave;
            X_ave = 0;

            for(int i1 = 0; i1 < n_ato; i1++)
            {
                delete [] X_rel [i1];
                X_rel[i1] =0;
            }
            delete [] X_rel;
            X_rel = 0;
 
            delete [] Coefs;
            Coefs = 0;          

            for (int i1 =0; i1 < total_dim; i1++)
            {
                delete [] d_Coefs[i1];
                d_Coefs[i1] = 0;
            }
            delete [] d_Coefs;
            d_Coefs = 0;

            delete [] d_Coefs_D;
            d_Coefs_D = 0;    
        }
    }
    
    void GetAllDerivsFull::setFirstDerivsCartChiral(std::vector<LIBMOL::AtomDict>& tAts,
                                                    std::vector<LIBMOL::ChiralDict>::iterator tCh)
    {
        int n_ato = (int)tCh->atoms.size();
         
        if (n_ato)
        {
            int dim = (int)tAts[tCh->atoms[0]].coords.size();
            
              int n_dim = 4*dim;
              int n_atoms = 4;
              int i, i_dim, i_ato;

              LIBMOL::REAL vol_ch   = 0.0;
              LIBMOL::REAL d_vol_ch = 0.0;
              LIBMOL::REAL weight = 0.0;

              LIBMOL::REAL * firstDerivChiras = new LIBMOL::REAL  [n_dim];
              for(i_dim = 0; i_dim < n_dim; i_dim++)
              {
                  firstDerivChiras[i_dim] = 0.0;
              }
        
              SetChiraAndFirstDeriv(tAts, tCh, vol_ch, firstDerivChiras);

              d_vol_ch = vol_ch - tCh->valueST;
  
              weight   = 2*tCh->fConst*d_vol_ch;
              
              for (i =0; i < n_atoms; i++)
              {
                  i_ato  = tCh->atoms[i];
                  i_dim = i*dim;      
    
                  firDrivCart[i_ato*dim]   -= (weight*firstDerivChiras[i_dim]);
                  
                  i_dim++;
                  firDrivCart[i_ato*dim+1] -= (weight*firstDerivChiras[i_dim]);
                  
                  i_dim++;
                  firDrivCart[i_ato*dim+2] -= (weight*firstDerivChiras[i_dim]);
              }
              
              delete [] firstDerivChiras;
              firstDerivChiras = 0;
        }
    }
    
    void GetAllDerivsFull::SetChiraAndFirstDeriv(std::vector<LIBMOL::AtomDict>& tAts,
                                                 std::vector<LIBMOL::ChiralDict>::iterator tCh,
                                                 LIBMOL::REAL  volume,  LIBMOL::REAL  * df_dx_c)
    {
        int n_ato = (int)tCh->atoms.size();
        if (n_ato)
        {
            int dim = (int)tAts[tCh->atoms[0]].coords.size();
            int i1,i2;
            int i_ato1, i_ato2, i_ato3, i_ato4;

            int d4 = 4*dim;
 
            // Three vectors
            
            LIBMOL::REAL *a = new LIBMOL::REAL [dim];
             
            LIBMOL::REAL *b = new LIBMOL::REAL [dim];
            LIBMOL::REAL *c = new LIBMOL::REAL [dim];
            
            // some of variables could be removed later, currently using old
            // new codes both 
            i_ato1 = tCh->atoms[0];
            i_ato2 = tCh->atoms[1];
            i_ato3 = tCh->atoms[2];
            i_ato4 = tCh->atoms[3];

            for (i1 = 0; i1 < dim; i1++)
            {
                a[i1] = tAts[i_ato2].coords[i1]-tAts[i_ato1].coords[i1];
                b[i1] = tAts[i_ato3].coords[i1]-tAts[i_ato1].coords[i1];
                c[i1] = tAts[i_ato4].coords[i1]-tAts[i_ato1].coords[i1];
            }
            
            // The derivatives of these three vectors wrt atomic coordinates

            LIBMOL::REAL **da_dx = new LIBMOL::REAL *[d4];

            for (i1 = 0; i1 < d4; i1++)
            {
                da_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    da_dx [i1][i2] = 0.0;
                }                
            }

            LIBMOL::REAL **db_dx = new LIBMOL::REAL *[d4];
            for (i1 = 0; i1 < 4*dim; i1++)
            {
                db_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    db_dx [i1][i2] = 0.0;
                }     
            }

            
            LIBMOL::REAL **dc_dx = new LIBMOL::REAL *[4*dim];
            for (i1 = 0; i1 < 4*dim; i1++)
            {
                dc_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    dc_dx [i1][i2] = 0.0;
                }     
            }
            
            for (i1 = 0; i1 < dim; i1++)
            {
                da_dx[i1][i1]       = -1.0;
                da_dx[dim+i1][i1]   =  1.0;
                db_dx[i1][i1]       = -1.0;
                db_dx[2*dim+i1][i1] =  1.0;
                dc_dx[i1][i1]       = -1.0;
                dc_dx[3*dim+i1][i1] =1.0;
            }

            // three auxi vectors (a X b), (b X c) and (c X a)

            LIBMOL::REAL * ab = new LIBMOL::REAL [dim];
            LIBMOL::REAL * bc = new LIBMOL::REAL [dim];
            LIBMOL::REAL * ca = new LIBMOL::REAL [dim];
            for(i1 = 0; i1 < dim; i1++)
            {
                ab[i1] = 0.0;
                bc[i1] = 0.0;
                ca[i1] = 0.0;
            }
 
            LIBMOL::CrossP(a, b, ab);
            LIBMOL::CrossP(b, c, bc);
            LIBMOL::CrossP(c, a, ca);

            // Get the value of the chiral volume
  
            volume = LIBMOL::CalcDet(a,b,c);

            // Finally, get the first derivatives of the chiral volume wrt 
            // atomic coordinates

            // cout << "First derivatives of Chiral center " <<index_ch+1 << endl;
 
            for (i1 = 0; i1 < d4; i1++)
            {
                df_dx_c[i1] = LIBMOL::DotP(da_dx[i1], bc)
                             +LIBMOL::DotP(db_dx[i1], ca)
                             +LIBMOL::DotP(dc_dx[i1], ab);
                //  cout << "  df_dx_c["<< i1<<"]= " << df_dx_c[i1]
                //       << endl;
            }       
            
            // release memory

            delete [] a;
            a = 0;
            delete [] b;
            b = 0;
            delete [] c;
            c =0;

            for (i1 = 0; i1 < d4; i1++)
            {
                delete [] da_dx [i1];
                da_dx [i1] = 0;
                delete [] db_dx [i1];
                db_dx [i1] = 0;
                delete [] dc_dx [i1];
                dc_dx [i1] = 0;
            }
            delete [] da_dx;
            da_dx = 0;
            delete [] db_dx;
            db_dx = 0;
            delete [] dc_dx;
            dc_dx = 0;

            delete [] ab;
            ab = 0;
            delete [] bc;
            bc = 0;
            delete [] ca;
            ca = 0; 

        }
    }
    
    void GetAllDerivsFull::setFirstDerivsCartVDW(std::vector<LIBMOL::AtomDict>& tAts,
                                                 std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        int dim = (int)tAt->coords.size();
        
        LIBMOL::REAL abs_r, r_com[dim];
        
        for (int i = 0; i < dim; i++)
        {
            r_com[i] = 0.0;
        }  
        
        for (std::vector<int>::iterator iNB = tAt->neighbAtoms.begin();
                iNB !=tAt->neighbAtoms.end(); iNB++)
        {
            if (std::find(tAt->connAtoms.begin(), tAt->connAtoms.end(), *iNB)==tAt->connAtoms.end()
                    && tAt->seriNum -1 < *iNB)
            {
                // Check which interaction should be activated Hydro-bond or VDW
                abs_r = 0.0;
                for (int i =0; i < dim; i++)
                {
                  r_com[i]   = tAt->coords[i]
                              -tAts[*iNB].coords[i];     
                  abs_r += pow(r_com[i],2.0);
                }
                abs_r = sqrt(abs_r);
                
                LIBMOL::REAL r_d  = SetVDWContact(tAts, iNB, tAt);
                
                if (abs_r < r_d)
                {  
                    LIBMOL::REAL vdwFC = 10.0;
                    for (int i =0; i < dim; i++)
                    {
                        LIBMOL::REAL tm = vdwFC
                                *(abs_r-r_d)*r_com[i]/abs_r;   
                        firDrivCart[(tAt->seriNum -1)*dim + i]-= tm;
                        firDrivCart[(*iNB)*dim + i]+= tm;

                    }
                }
                
            }
        }
    }
    
    LIBMOL::REAL GetAllDerivsFull::SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
                                            std::vector<int>::iterator tNB, 
                                            std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        LIBMOL::REAL vCont = VDWCONST*2;
        
        if(tAt->ionRadius <=1.0e-8 || tAts[*tNB].ionRadius <=1.0e-8)
        {
            vCont = tAt->radius + tAts[*tNB].radius;
        }
        else
        {
            vCont = tAt->ionRadius + tAts[*tNB].ionRadius;
        }
        
        return vCont;
    }
    
    // That function are used when workMode==1
    // 1. the simplest one, given all atom forces have been calculated
    void GetAllDerivsFull::SetNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts)
    {
        // Initialize the normal matrix
        for (int i=0; i < numVars ; i++)
        {
            for (int j =0; j < numVars; j++)
            {
                secDrivCart[i][j] =0.0;
            }
        }
        /*
        for (std::vector<LIBMOL::AtomDict>::iterator iA1=tAts.begin();
                iA1 !=tAts.end(); iA1++)
        {
            int dim1 = (int)iA1->forces.size();
            for (std::vector<LIBMOL::AtomDict>::iterator iA2=tAts.begin();
                    iA2 !=tAts.end(); iA2++)
            {
                int dim2 = (int)iA2->forces.size();
                for (int i1=0; i1 < dim1; i1++)
                {
                    int k1 = iA1->seriNum*dim1 +i1;
                    for (int i2=0; i2 < dim2; i2++)
                    {
                        int k2 = iA2->seriNum*dim2 + i2;
                        
                    }
                }
            }
        }
         */
    }
    
    void GetAllDerivsFull::SetNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts, 
                                           std::vector<LIBMOL::BondDict>& tBos, 
                                           std::vector<LIBMOL::AngleDict>& tAns, 
                                           std::vector<LIBMOL::TorsionDict>& tTos, 
                                           std::vector<LIBMOL::PlaneDict>& tPls, 
                                           std::vector<LIBMOL::ChiralDict>& tChs)
    {
        // Initialize the normal matrix
        for (int i=0; i < numVars ; i++)
        {
            for (int j =0; j < numVars; j++)
            {
                secDrivCart[i][j] =0.0;
            }
        }
        
        for (std::vector<LIBMOL::BondDict>::iterator iBo =tBos.begin();
                iBo != tBos.end(); iBo++)
        {
            OneBondToNormalMatrix(tAts, iBo);
        }
        
        std::cout << "Done bond to normal matrix" << std::endl; 
       
        if ( (int)tAns.size()!=0)
        {
            AngToNormalMatrix(tAts, tAns);
        }
        std::cout << "Done angle to normal matrix" << std::endl;
        /*
        if ((int)tTos.size() !=0)
        {
         
            TorToNormalMatrix(tAts, tTos);
        }
        */
        
        if ((int)tPls.size() !=0)
        {
            bool l_first = true;
            PlaToNormalMatrix(tAts, tPls, l_first);
        }
        
        std::cout << "Done plane to normal matrix" << std::endl;
        
        /*
        if((int)tChs.size() !=0)
        {
            bool l_first = true;
            ChiToNormalMatrix(tAts, tChs, l_first);
        }
        */
        
        for (std::vector<LIBMOL::AtomDict>::iterator iAt=tAts.begin();
                iAt !=tAts.end(); iAt++)
        {
            OneAtomVDWToNormalMatrix(tAts, iAt);
        }
        std::cout << "Done vdw to normal matrix" << std::endl;
        // check 
        /*
        for(std::vector<LIBMOL::AtomDict>::iterator iA1=tAts.begin();
                iA1 !=tAts.end(); iA1++)
        {
            int dim1 = (int)iA1->coords.size();
            for (int i=0; i < (int)iA1->coords.size(); i++)
            {
                
                int k1 =iA1->seriNum*dim1 +i;
                for (std::vector<LIBMOL::AtomDict>::iterator iA2=tAts.begin();
                        iA2 != tAts.end(); iA2++)
                {
                    int dim2 =(int)iA2->coords.size();
                    for(int j=0; j < (int)iA2->coords.size(); j++)
                    {
                        int k2 = iA2->seriNum*dim2 + j;
                        std::cout << "second_deriv[" <<k1<<"]["<<k2<<"]="
                                << secDrivCart[k1][k2] << std::endl;
                    }
                }
            }
        }
         */
    }
    
    void GetAllDerivsFull::StableNormalMatrix()
    {
        LIBMOL::REAL maxElem=0.0;
        
        if (numVars !=0)
        {
            for (int i=0; i < numVars; i++)
            {
                for (int j=0; j < numVars; j++)
                {
                    if (fabs(secDrivCart[i][j]) > maxElem)
                    {
                        maxElem = fabs(secDrivCart[i][j]);
                    }
                }
            }
            
            if (maxElem >1.0e-8)
            {
                for (int i=0; i < numVars; i++)
                {
                    for (int j=0; j < numVars; j++)
                    {
                        secDrivCart[i][j]= secDrivCart[i][j]/maxElem;
                    }
                }
            }
        }
        
    }
    void GetAllDerivsFull::OneBondToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                                                 std::vector<LIBMOL::BondDict>::iterator tBo)
    {
        if ((int)tBo->atomsIdx.size() !=0)
        {
            int dim = (int)tAts[tBo->atomsIdx[0]].coords.size();
            LIBMOL::REAL dist_comp[dim];
            for (int i=0; i < dim; i++)
            {
                dist_comp[i] = 0;
            }
            
            LIBMOL::REAL bondL=0.0;
            
            LIBMOL::REAL * df_dx_b = new  LIBMOL::REAL [2*dim];
            for (int i = 0; i < 2*dim; i++)
            {
                df_dx_b[i] = 0.0;
            }
            
            int iAt1 = tBo->atomsIdx[0], iAt2=tBo->atomsIdx[1];
            
            for (int i =0; i < dim; i++)
            {
                dist_comp[i] =  tAts[iAt1].coords[i]
                                -tAts[iAt2].coords[i];
                bondL +=  (pow(dist_comp[i],2.0));
            }
            bondL = sqrt(bondL);
            
            
            for (int j =0; j < dim; j++)
            {  
                df_dx_b[j] =   dist_comp[j]/bondL;
                df_dx_b[dim+j] = -dist_comp[j]/bondL;
            }    

            for (int i1 =0; i1 < 2; i1++)
            {
                for (int j1=0; j1 < dim; j1++)
                {
                    for(int i2 =0; i2 < 2; i2++)
                    {   
                        for(int j2 =0; j2 < dim; j2++)
                        {
                            secDrivCart[iAt1*dim+j1][iAt2*dim+j2] +=
                              (df_dx_b[i1*dim+j1]*df_dx_b[i2*dim+j2]/tBo->sigValueST);
                        }
                    }
                }
            }

        }
        
    }
    
    
    void GetAllDerivsFull::AngToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                                  std::vector<LIBMOL::AngleDict>& tAns)
    {
        
        int i,j,k;
        int nAns = (int)tAns.size();
        int dim = (int)tAts[0].coords.size();
        LIBMOL::REAL **   dThetadx = new LIBMOL::REAL * [nAns];
        for (i =0; i < nAns; i++)
        {
            dThetadx[i] = new LIBMOL::REAL [3*dim];
            for (j =0; j < 3*dim; j++)
            {
                dThetadx[i][j] = 0.0;
            }
        }
        
        int i_atom1, i_atom2, i_atom3; 

        LIBMOL::REAL *a;
        LIBMOL::REAL *b;
        LIBMOL::REAL leng_a, leng_b;
        LIBMOL::REAL d_leng_a[3*dim];
        LIBMOL::REAL d_leng_b[3*dim];
        LIBMOL::REAL da_dx[3*dim][dim];
        LIBMOL::REAL db_dx[3*dim][dim];

        LIBMOL::REAL R, dR_dx;

        a = new LIBMOL::REAL [dim];
        b = new LIBMOL::REAL [dim];
        for (i=0; i < dim; i++)
        {
            a[i] =0.0;
            b[i] =0.0;
        }

        for ( i =0; i <nAns; i++)
        {
            i_atom1 = tAns[i].atoms[1];
            i_atom2 = tAns[i].atoms[0];
            i_atom3 = tAns[i].atoms[2];
            
            for (j =0; j < dim; j++)
            {
                a[j] = tAts[i_atom1].coords[j]-tAts[i_atom2].coords[j];
                b[j] = tAts[i_atom3].coords[j]-tAts[i_atom2].coords[j]; 
            }
      
            leng_a = LIBMOL::length_v(a, dim);
            leng_b = LIBMOL::length_v(b, dim);
  
            R = LIBMOL::DotP(a,b)/(leng_a*leng_b);
            if(fabs(1-R*R) < 1.0e-7)
            {
                // cout << " Bond angle " << i+1 << endl;
                //  cout << "leng_a = " << leng_a << endl;
                //  cout << "leng_b = " << leng_b << endl;
                // cout << "R      = " << R << endl;
         
                LIBMOL::REAL delta_ram = 1.0e-2;
                for (j =0; j < dim; j++)
                {
                    tAts[i_atom3].coords[j] 
                    = tAts[i_atom3].coords[j] + delta_ram*LIBMOL::GetRand();
                    b[j] = tAts[i_atom3].coords[j]-tAts[i_atom2].coords[j];
                }
      
                leng_b = LIBMOL::length_v(b, dim);

                if( leng_b ==0)
                {
                    std::cout << "vector a or b equals to zero in the bond angle : " 
                              << i << std::endl;
                  
                    exit(0);
                }
      
                R = LIBMOL::DotP(a,b)/(leng_a*leng_b);
         
                // cout << "new R      = " << R << endl;
         
            }
            
            for (j =0; j < 3*dim; j++)
            {
                d_leng_a[j] = 0.0;
                d_leng_b[j] = 0.0;
            }

            for( j =0; j < 3*dim; j++)
            {
                for (k =0; k < dim; k++)
                {
                    da_dx[j][k] = 0.0;
                    db_dx[j][k] = 0.0;
                }
            }
        
            for ( j =0; j < dim; j++)
            {
                da_dx[j][j]        =  1.0;
                da_dx[j+dim][j]    = -1.0;
                db_dx[j+dim][j]    = -1.0;
                db_dx[j+2*dim][j]  =  1.0;
            }

            for (j = 0; j < dim; j++)
            {
                d_leng_a[j]       =  a[j]/leng_a;
                d_leng_a[j+dim]   = -d_leng_a[j];
                d_leng_b[j+dim]   = -b[j]/leng_b;
                d_leng_b[j+2*dim] = -d_leng_b[j+dim];          
            }
            
            for (j =0; j < 3*dim; j++)
            {
                dR_dx = (LIBMOL::DotP(da_dx[j],b)+LIBMOL::DotP(db_dx[j],a))/(leng_a*leng_b)
                         -LIBMOL::DotP(a,b)*(leng_b*d_leng_a[j]+leng_a*d_leng_b[j])
                         /pow(leng_a*leng_b,2.0);

                if(fabs(1-R*R) > 1.0e-8)
                {
                    dThetadx[i][j] = -dR_dx/sqrt(1-R*R);
                }
                else
                {
                    // cout << " 1-R*R = " << fabs(1-R*R) << endl;
                    std::cout << "angles calculation " << std::endl;
                    std::cout << "Cos(theta) = 1, atoms overlapped ! " << std::endl;
                    exit(1);
                } 
                // cout << " dTheta/dx["<<j<<"] = " 
                //      << dtheta_dx[i][j] << endl;
            }
        }
      
        for (i =0; i < nAns; i++)
        {
            for (int i1 =0; i1 < 3; i1++)
            {
                for (int j1=0; j1 < dim; j1++)
                {
                    for(int i2 =0; i2 < 3; i2++)
                    {    
                        for(int j2 =0; j2 < dim; j2++)
                     {
                            int k1 =tAns[i].atoms[i1]*dim+j1;
                            int k2 =tAns[i].atoms[i2]*dim+j2;
                            int k3 = i1*dim+j1;
                            int k4 = i2*dim+j2;
                            
                            secDrivCart[k1][k2] +=
                               (dThetadx[i][k3]*dThetadx[i][k4]/tAns[i].sigValueST);
                     }
                 }
              }
           }
        }
        
        
        delete [] a;
        a = 0;
 
        delete [] b;
        b = 0;  
        
        for (i =0; i <nAns; i++)
        {
          delete []  dThetadx[i];
          dThetadx[i]   =0;
        }    
        delete [] dThetadx;
        dThetadx     = 0;
    }
    
    void GetAllDerivsFull::TorToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts,
                                             std::vector<LIBMOL::TorsionDict> & tTos)
    {
        int nTors = (int)tTos.size();
        int dim   = tAts[0].coords.size();
        int i,j,k;
        int i1,i2,j1,j2, k1, k2,k3,k4;

        LIBMOL::REAL  **   dPsidx   = new LIBMOL::REAL   * [nTors];
        LIBMOL::REAL  ***  dPsidxdx = new LIBMOL::REAL  ** [nTors];
        for (i =0; i < nTors; i++)
        {
            dPsidx[i]    =  new LIBMOL::REAL   [4*dim];
            dPsidxdx[i]  =  new LIBMOL::REAL * [4*dim];
            for (j =0; j < 4*dim; j++)
            {
                dPsidx[i][j] = 0.0;
                dPsidxdx[i][j]  = new LIBMOL::REAL [4*dim];
                
                for (k =0; k < 4*dim; k++)
                {
                    dPsidxdx[i][j][k] = 0.0;
                }
            }
        }
        
        for (int iTo=0; iTo < nTors; iTo++)
        {
            int t_dim = 4*dim;
            
            LIBMOL::REAL A, B;
            LIBMOL::REAL B1;
            LIBMOL::REAL *  dAdx = new LIBMOL::REAL [t_dim];
            LIBMOL::REAL *  dBdx = new LIBMOL::REAL [t_dim];
            LIBMOL::REAL *  dB1  = new LIBMOL::REAL [t_dim];
            for ( i =0; i < t_dim; i++)
            {
                dAdx[i] = 0.0;
                dBdx[i] = 0.0;
                dB1[i]  = 0.0;
            }
            
            LIBMOL::REAL ** dAdxdx = new LIBMOL::REAL * [t_dim];
            LIBMOL::REAL ** dBdxdx = new LIBMOL::REAL * [t_dim];
            for (i =0; i <t_dim; i++)
            {
                dAdxdx[i] = new LIBMOL::REAL [t_dim];
                dBdxdx[i] = new LIBMOL::REAL [t_dim];
                for (j =0; j < t_dim; j++)
                {
                    dAdxdx[i][j] = 0.0;
                    dBdxdx[i][j] = 0.0;
                }      
            }
            
            LIBMOL::REAL *vect_a  = new LIBMOL::REAL [dim];
            LIBMOL::REAL *vect_b  = new LIBMOL::REAL [dim];
            LIBMOL::REAL *vect_c  = new LIBMOL::REAL [dim];
            LIBMOL::REAL *vect_bc = new LIBMOL::REAL [dim];        // cross product between b and c
            
            for ( i =0; i < dim; i++)
            {
                vect_a[i]   =  0;
                vect_b[i]   =  0;
                vect_b[i]   =  0;
                vect_bc[i]  =  0;
            }

            LIBMOL::REAL ** d_vect_a = new LIBMOL::REAL* [t_dim];
            LIBMOL::REAL ** d_vect_b = new LIBMOL::REAL* [t_dim];
            LIBMOL::REAL ** d_vect_c = new LIBMOL::REAL* [t_dim];
            for(i = 0; i < t_dim; i++)
            {
                d_vect_a [i] = new LIBMOL::REAL [dim];
                d_vect_b [i] = new LIBMOL::REAL [dim];
                d_vect_c [i] = new LIBMOL::REAL [dim];
                for(j = 0; j < dim; j++)
                {
                    d_vect_a [i][j] = 0.0;
                    d_vect_b [i][j] = 0.0;
                    d_vect_c [i][j] = 0.0;
                }
            }
            
            LIBMOL::REAL abs_a;  
            LIBMOL::REAL abs_b;
            LIBMOL::REAL abs_c;

            LIBMOL::REAL *d_abs_b   = new LIBMOL::REAL [t_dim];
            LIBMOL::REAL **dd_abs_b = new LIBMOL::REAL * [t_dim];
            for (i = 0; i < t_dim; i++)
            {
                d_abs_b[i] = 0.0;
                dd_abs_b[i] = new LIBMOL::REAL [t_dim];
                for (j =0; j < t_dim; j++)
                {
                    dd_abs_b[i][j] = 0.0;
                }
            }
            
            LIBMOL::REAL * vect_db_c = new LIBMOL::REAL [dim];
            LIBMOL::REAL * vect_b_dc = new LIBMOL::REAL [dim];
            for( i =0; i < dim; i++)
            {
                vect_db_c[i] = 0.0;
                vect_b_dc[i] = 0.0;
            }
            
            LIBMOL::REAL * vect_dbj_c   = new LIBMOL::REAL [dim];
            LIBMOL::REAL * vect_b_dcj   = new LIBMOL::REAL [dim];
            LIBMOL::REAL * vect_dbi_dcj = new LIBMOL::REAL [dim];
            LIBMOL::REAL * vect_dbj_dci = new LIBMOL::REAL [dim];

            for (i =0; i < dim; i++)
            {
                vect_dbj_c[i]   = 0.0;
                vect_b_dcj[i]   = 0.0;
                vect_dbi_dcj[i] = 0.0;
                vect_dbj_dci[i] = 0.0; 
            }
            
            // Find the atoms that construct the current torsion angle

            int    i_atom1, i_atom2, i_atom3, i_atom4;
            /*
            //if(tTos[i].compAtomNum[1] == baseAtomNum)
            //{
                // dummy atom involved
                
                //i_atom1  =  atomsDummy[0].serNum ;  
                //i_atom2  =  torsAngs[index_torsion].compAtomNum[1]-1;
                //i_atom3  =  torsAngs[index_torsion].compAtomNum[2]-1;
                //i_atom4  =  torsAngs[index_torsion].compAtomNum[3]-1;      
            //}
            //else if(torsAngs[index_torsion].compAtomNum[1] != baseAtomNum
            // && torsAngs[index_torsion].compAtomNum[1] != -1
            && index_torsion < numBBTorsAngs)
            {    
            */
            i_atom1  =  tTos[iTo].atoms[0];
            i_atom2  =  tTos[iTo].atoms[1];
            i_atom3  =  tTos[iTo].atoms[2];
            i_atom4  =  tTos[iTo].atoms[3];
            /*
            }
            else
            {
            cout << "The second atom in torsion angle  " 
            << index_torsion-1 << "is " 
            << torsAngs[index_torsion].compAtomNum[1] << endl;
            cout << "We need to check *.str file " << endl;
            exit(1);
            }
            */

   
            std::cout << "current torsion angle is " << i << std::endl;
            std::cout << "four atoms that construct the angle are: " 
                      << i_atom1  << "    "
                      << i_atom2  << "    "
                      << i_atom3  << "    "
                      << i_atom4  << std::endl;

            // Set vectors a, b, c  and |b| 

            abs_a = 0.0;
            abs_b = 0.0;
            abs_c = 0.0;

            for (i =0; i < dim; i++)
            {
            //   if(i_atom1 == atomsDummy[0].serNum)
            //{
            // vect_a[i] = atoms[i_atom2].coords[i] - atomsDummy[0].coords[i];
            //}
            //else
            // {
               vect_a[i] = tAts[i_atom2].coords[i] - tAts[i_atom1].coords[i];
            //}
               vect_b[i] = tAts[i_atom3].coords[i] - tAts[i_atom2].coords[i];
               vect_c[i] = tAts[i_atom4].coords[i] - tAts[i_atom3].coords[i];
            //   cout << "i_tors " << index_torsion -1 << endl;
            //   cout << "vect_a[" << i <<"] " << vect_a[i] <<  endl;
            //        << atom2.coords[i] << "\t" <<  atom1.coords[i] << endl;
            //   cout << "vect_b[" << i <<"] " << vect_b[i] << "\t"
            //        << atom3.coords[i] << "\t" <<  atom2.coords[i] << endl;
            //   cout << "vect_c[" << i <<"] " << vect_c[i] << "\t"
            //        << atom4.coords[i] << "\t" <<  atom3.coords[i] << endl;
               
               abs_a += pow(vect_a[i],2);
               abs_b += pow(vect_b[i],2);
               abs_c += pow(vect_b[i],2);

            }
           
            abs_a = sqrt(abs_a);
            abs_b = sqrt(abs_b);
            abs_c = sqrt(abs_c);

            if(abs_a < 1.0e-16)
            {
               std::cout << "atom 1 and atom2 in torsion angle " 
                         << iTo
                         << " are overlapped. " << std::endl;
               exit(1);
            }

            if(abs_b < 1.0e-16)
            {
               std::cout << "atom 2 and atom 3 in torsion angle " 
                         << iTo
                         << " are overlapped. " << std::endl;
               exit(1);
            }

            if(abs_c < 1.0e-16)
            {
               std::cout << "atom 3 and atom 4 in torsion angle " << iTo
                         << " are overlapped. " << std::endl;
               exit(1);
            }

            // std::cout << "abs of vector c is  " << abs_b << endl;                 

            // calculate the first derivatives of vectors a, b and c 
            // with respect to xi

            // 1. vector a 

            d_vect_a[0][0] =  -1.0;      // da/dx1
            d_vect_a[1][1] =  -1.0;      // da/dy1
            d_vect_a[2][2] =  -1.0;      // da/dz1
  
            d_vect_a[3][0] =   1.0;      // da/dx2
            d_vect_a[4][1] =   1.0;      // da/dy2
            d_vect_a[5][2] =   1.0;      // da/dz2

            // 2. vector b

            d_vect_b[3][0] =  -1.0;      // db/dx2
            d_vect_b[4][1] =  -1.0;      // db/dy2
            d_vect_b[5][2] =  -1.0;      // db/dz2

            d_vect_b[6][0] =   1.0;      // db/dx3
            d_vect_b[7][1] =   1.0;      // db/dy3
            d_vect_b[8][2] =   1.0;      // db/dz3

            // 3. vector c   

            d_vect_c[6][0] =  -1.0;      // dc/dx3
            d_vect_c[7][1] =  -1.0;      // dc/dy3
            d_vect_c[8][2] =  -1.0;      // dc/dz3

            d_vect_c[9][0] =   1.0;      // dc/dx4
            d_vect_c[10][1] =  1.0;      // dc/dy4
            d_vect_c[11][2] =  1.0;      // dc/dz4

 
  
            // 4. absolute value of vector b

            d_abs_b[3]     =   -(tAts[i_atom3].coords[0]-tAts[i_atom2].coords[0])/abs_b;
            d_abs_b[4]     =   -(tAts[i_atom3].coords[1]-tAts[i_atom2].coords[1])/abs_b;
            d_abs_b[5]     =   -(tAts[i_atom3].coords[2]-tAts[i_atom2].coords[2])/abs_b;
  
            d_abs_b[6]     =   -d_abs_b[3];
            d_abs_b[7]     =   -d_abs_b[4];
            d_abs_b[8]     =   -d_abs_b[5];

            //     cout << " d|b|/dx2 = " << d_abs_b[3] << endl;
            //     cout << " d|b|/dy2 = " << d_abs_b[4] << endl;
            //     cout << " d|b|/dz2 = " << d_abs_b[5] << endl;
           
  
            // 5. the second derivatives of absolute values of the vector b
   
            for (i =3; i < 6; i++)
            {
                for (j = i; j < 9; j++)
                {
                   if( j == i)
                   {
                       dd_abs_b[i][j] = (1.0-d_abs_b[i]*d_abs_b[j])/abs_b;
                   }
                   else if ( j == i +3)
                   {
                       dd_abs_b[i][j] = -(1.0+d_abs_b[i]*d_abs_b[j])/abs_b;
                   }
                   else
                   {
                       dd_abs_b[i][j] = -d_abs_b[i]*d_abs_b[j]/abs_b;
                   }
                }
            }
 
            for (i =6; i < 9; i++)
            {
               for ( j = i; j <9; j++)
               {
                   if( j == i)
                   {
                       dd_abs_b[i][j] = (1.0+d_abs_b[i]*d_abs_b[j])/abs_b;
                   }
                   else
                   {
                       dd_abs_b[i][j] = d_abs_b[i]*d_abs_b[j]/abs_b;
                   }
               }
            }
           
           
            for (i = 4; i < 9; i++)
            {
               for ( j = 3; j < i; j++)
               {
                   dd_abs_b[i][j] = dd_abs_b[j][i];
               }
            }
           
            // Calculate A and B which determine the torsion angle.
            // Calculate also the derivatives of A and B with respect to xi
            // 1. A and B
           
            A   = 0.0;
            A   =-LIBMOL::DotP(vect_b,vect_b)*LIBMOL::DotP(vect_a,vect_c)
                 +LIBMOL::DotP(vect_a,vect_b)*LIBMOL::DotP(vect_b,vect_c);
  
            B   = 0.0;

            LIBMOL::CrossP(vect_b, vect_c, vect_bc);

            B1  = LIBMOL::DotP(vect_a,vect_bc);
            B   = abs_b*B1;

   
            // 2. dA/dxi and dB/dxi

            for (i = 0; i < t_dim; i++)
            {
               dAdx[i] = -2*LIBMOL::DotP(d_vect_b[i],vect_b)*LIBMOL::DotP(vect_a,vect_c)
                         -LIBMOL::DotP(vect_b,vect_b)*(LIBMOL::DotP(d_vect_a[i],vect_c)
                         +LIBMOL::DotP(vect_a,d_vect_c[i]))
                         +LIBMOL::DotP(vect_b,vect_c)*(LIBMOL::DotP(d_vect_a[i],vect_b)
                                                       +LIBMOL::DotP(vect_a,d_vect_b[i]))
                         +LIBMOL::DotP(vect_a,vect_b)*(LIBMOL::DotP(d_vect_b[i],vect_c)
                                                       +LIBMOL::DotP(vect_b,d_vect_c[i]));
    
               LIBMOL::CrossP(d_vect_b[i], vect_c, vect_db_c);
               LIBMOL::CrossP(vect_b, d_vect_c[i], vect_b_dc);
               dB1[i] = LIBMOL::DotP(d_vect_a[i],vect_bc)
                       +LIBMOL::DotP(vect_a,vect_db_c)
                       +LIBMOL::DotP(vect_a, vect_b_dc);
               dBdx[i] = d_abs_b[i]*B1+abs_b*dB1[i];
               //   std::cout << " dA/dx[" << i << "]=" <<  dAdx[i] << std::endl;
               //   std::cout << " dB/dx[" << i << "]=" <<  dBdx[i] << std::endl;
            }
  
            // 3. d^2A/dxidxj and d^2B/dxidxj

            for (i =0; i < t_dim; i++)
            {
               LIBMOL::CrossP(vect_b, d_vect_c[i], vect_b_dc);
               LIBMOL::CrossP(d_vect_b[i], vect_c, vect_db_c);

               for ( j =0; j < t_dim; j++)
               {
                   LIBMOL::CrossP(d_vect_b[j], vect_c, vect_dbj_c);
                   LIBMOL::CrossP(vect_b, d_vect_c[j], vect_b_dcj);
                   LIBMOL::CrossP(d_vect_b[i], d_vect_c[j], vect_dbi_dcj);
                   LIBMOL::CrossP(d_vect_b[j], d_vect_c[i], vect_dbj_dci);
           
                   dBdxdx[i][j] = LIBMOL::DotP(d_vect_a[i],vect_dbj_c)
                                 +LIBMOL::DotP(d_vect_a[i],vect_b_dcj)
                                 +LIBMOL::DotP(d_vect_a[j],vect_db_c)
                                 +LIBMOL::DotP(vect_a,vect_dbi_dcj)
                                 +LIBMOL::DotP(d_vect_a[j],vect_b_dc)
                                 +LIBMOL::DotP(vect_a,vect_dbj_dci);
         
                   dBdxdx[i][j] = dd_abs_b[i][j]*B1+d_abs_b[i]*dB1[j]
                                 +d_abs_b[j]*dB1[i]+ abs_b*dBdxdx[i][j];


                   dAdxdx[i][j] =LIBMOL::DotP(vect_a,vect_b)*(LIBMOL::DotP(d_vect_b[i],d_vect_c[j])
                                +LIBMOL::DotP(d_vect_b[j],d_vect_c[i])) 
                                -LIBMOL::DotP(vect_b,vect_b)*(LIBMOL::DotP(d_vect_a[i],d_vect_c[j])
                                +LIBMOL::DotP(d_vect_a[j],d_vect_c[i]))            
                                +LIBMOL::DotP(vect_b,vect_c)*(LIBMOL::DotP(d_vect_b[i],d_vect_a[j])
                                +LIBMOL::DotP(d_vect_b[j],d_vect_a[i]))
                                -2*LIBMOL::DotP(vect_a,vect_c)*LIBMOL::DotP(d_vect_b[i],d_vect_b[j])
                                -2*LIBMOL::DotP(d_vect_b[i],vect_b)*(LIBMOL::DotP(d_vect_a[j],vect_c)
                                +LIBMOL::DotP(d_vect_c[j],vect_a))
                                -2*LIBMOL::DotP(d_vect_b[j],vect_b)*(LIBMOL::DotP(d_vect_a[i],vect_c)
                                +LIBMOL::DotP(vect_a,d_vect_c[i]))
                                +(LIBMOL::DotP(d_vect_a[i],vect_b)+LIBMOL::DotP(d_vect_b[i],vect_a))
                                *(LIBMOL::DotP(d_vect_b[j],vect_c)+LIBMOL::DotP(d_vect_c[j],vect_b))
                                +(LIBMOL::DotP(d_vect_a[j],vect_b)+LIBMOL::DotP(d_vect_b[j],vect_a))
                                  *(LIBMOL::DotP(d_vect_b[i],vect_c)+LIBMOL::DotP(d_vect_c[i],vect_b));
         
                                //  cout << "d^2A/dx_"<<i<<"dx_"<<j<< " = " << dAdxdx[i][j] << endl;
                                //  cout << "d^2B/dx_"<<i<<"dx_"<<j<< " = " << dBdxdx[i][j] << endl;
            
               }
            }   
           
            // Calculate the first derivative of a torsion angle with respect to xi

            // New definitions

            //for ( i = 0; i < 4*dim; i++)
            //  {
            //    if(B == 0)
            //    {
            //        dPsidX[i] = 0.0;
            //    }
            //    else 
            //    {
          
            for ( i =0; i < 4*dim; i++)
            {
               dPsidx[iTo][i] = (dAdx[i]*B-A*dBdx[i])/(pow(A,2.0)+pow(B,2.0));

               //  cout << "dPsi/dx_"<< i << " =  " << dPsidX[i] << endl;

            }

  
  
            // Calculate the second derivatives of a torsion angle 

            LIBMOL::REAL AB2;
  
            AB2 = pow(A,2.0)+pow(B,2.0);

            // cout << "A*A+B*B ="  << AB2 << endl;
            // cout << "index_torsion " << index_torsion-1 << endl;

            for(i = 0; i < t_dim; i++)
            {
               for( j =0; j <t_dim; j++)
               {
                   dPsidxdx[iTo][i][j] =-2*(A*dAdx[j]+B*dBdx[j])
                                 *(dAdx[i]*B-A*dBdx[i])/pow(AB2,2.0)
                                 +(dAdxdx[i][j]*B+dAdx[i]*dBdx[j]-dAdx[j]*dBdx[i]
                                 -A*dBdxdx[i][j])/AB2;
                   //   cout << "dPsi/dx_" <<i <<"dx_"<<j<< " = "
                   //        <<   dPsidXdX[i][j] << endl;        
               }
            }
            delete [] vect_a;
            vect_a = 0;
            delete [] vect_b;
            vect_b = 0;
            delete [] vect_c;
            vect_c = 0;
            delete [] vect_bc;
            vect_bc = 0;

            for(i = 0; i < t_dim; i++)
            {
               delete [] d_vect_a[i];
               d_vect_a[i] = 0;
            }
            delete [] d_vect_a;
            d_vect_a = 0;

            for(i = 0; i < t_dim; i++)
            {
               delete [] d_vect_b[i];
               d_vect_b[i] = 0;
            }
            delete [] d_vect_b;
            d_vect_b = 0;
  
            for(i = 0; i < t_dim; i++)
            {
               delete [] d_vect_c[i];
               d_vect_c[i] = 0;
            }
            delete [] d_vect_c;
            d_vect_c = 0;


            delete [] dAdx;
            dAdx    = 0;
            delete [] dBdx;
            dBdx    = 0; 
            delete [] dB1;
            dB1     = 0;
  
            for (i =0; i <t_dim; i++)
            {
               delete [] dAdxdx[i];
               dAdxdx[i]  = 0;
               delete [] dBdxdx[i];
               dBdxdx[i] = 0;    
            }
            delete []  dAdxdx;
            dAdxdx  = 0;

            for (i =0; i <t_dim; i++)
            {
               delete [] dBdxdx[i];
               dBdxdx[i] = 0;    
            }
            delete []  dBdxdx;
            dBdxdx  = 0;

            delete [] vect_db_c;
            vect_db_c = 0;

            delete [] vect_b_dc;
            vect_b_dc = 0;

            delete [] vect_dbj_c;
            vect_dbj_c =0;
 
            delete [] vect_b_dcj;
            vect_b_dcj =0;

            delete [] vect_dbi_dcj;
            vect_dbi_dcj = 0;

            delete [] d_abs_b;
            d_abs_b = 0;
     
            for(i = 0; i < t_dim; i++)
            {
               delete [] dd_abs_b[i];
               dd_abs_b[i] = 0;
            }
            delete [] dd_abs_b;  
            dd_abs_b = 0; 
               
        }
        
        for (int iTo =0; iTo < nTors; iTo++)
        {
            for (i1 =0; i1 < 4; i1++)
            {
                for (j1=0; j1 < dim; j1++)
                {
                    for(i2 =0; i2 < 4; i2++)
                    {      
                        for(j2 =0; j2 < dim; j2++)
                        {
                            k1 =tTos[iTo].atoms[i1]*dim+j1;
                            k2 =tTos[iTo].atoms[i2]*dim+j2;
                            k3 = i1*dim+j1;
                            k4 = i2*dim+j2;
                            secDrivCart[k1][k2] +=
                            (dPsidx[iTo][k3]*dPsidx[iTo][k4]/tTos[iTo].sigValueST);
                        }
                    }
                }
            }
        }
        
        // release the memory occupied

        for (i =0; i < nTors; i++)
        {
            for (j =0; j < 4*dim; j++)
            {
                delete []  dPsidxdx[i][j];
                dPsidxdx[i][j]  = 0; 
            }

            delete [] dPsidx[i];
            dPsidx[i]    = 0;
            delete [] dPsidxdx[i];
            dPsidxdx[i]  = 0;
        }

        delete  dPsidx;
        dPsidx   =0;
  
        delete  dPsidxdx;
        dPsidxdx =0; 
    }
    
    void GetAllDerivsFull::PlaToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts, 
                                             std::vector<LIBMOL::PlaneDict> & tPls,
                                             bool l_fd)
    {
         
         int dim   = (int)tAts[0].coords.size();
         
         for (std::vector<LIBMOL::PlaneDict>::iterator iPl = tPls.begin();
                 iPl != tPls.end(); iPl++)
         {
             // Setup the first derivatives
      
             
             std::vector<int> t_atoms;
             for (std::map<LIBMOL::ID,int>::iterator i_at = iPl->atoms.begin();
                     i_at != iPl->atoms.end(); i_at++)
             {
                 t_atoms.push_back(i_at->second);
             }
             
             // Setup the first derivatives
      
             int n_atoms   =  (int)iPl->atoms.size();
             int n_tol_dim = n_atoms*dim;

             LIBMOL::REAL ** firstDerivPlans = new LIBMOL::REAL * [n_atoms];
             for(int i_atom = 0; i_atom < n_atoms; i_atom++)
             {
                 firstDerivPlans[i_atom] = new LIBMOL::REAL [n_tol_dim];
                 for (int i = 0; i < n_tol_dim; i++)
                 {
                     firstDerivPlans [i_atom][i] = 0.0;
                 }
             }

             
             if (!l_fd)
             {
                 setFirstDerivsCartPlane(tAts, iPl, firstDerivPlans);
             }
             
             for (int i_atom = 0; i_atom < n_atoms; i_atom++)
             {
                 for (int i1 = 0; i1 < n_atoms; i1++)
                 {
                 
                     for (int j1 =0; j1 < dim; j1++)
                     {
                         int k1 = t_atoms[i1]*dim +j1;
                         int k3 = i1*dim +j1;
                         for (int i2 = 0; i2 <n_atoms; i2++)
                         {
                             for(int j2=0; j2 < dim; j2++)
                             {
                                 int k2 = t_atoms[i2]*dim + j2;
                                 int k4 = i2*dim + j2;
                            
                                 secDrivCart[k1][k2] += 0.50*iPl->fConst*
                                       (firstDerivPlans[i_atom][k3]*firstDerivPlans[i_atom][k4]);
                             }
                         }
                     }
                 }
             }        
         }

        
    }
    
    
    
    void GetAllDerivsFull::ChiToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts, 
                                             std::vector<LIBMOL::ChiralDict> & tChs,
                                             bool l_fd)
    {
        for(std::vector<LIBMOL::ChiralDict>::iterator iCh =tChs.begin();
                iCh != tChs.end(); iCh++)
        {
            if(!l_fd)
            {
                setFirstDerivsCartChiral(tAts, iCh);
            }
            
        }
        
    }
    
    void GetAllDerivsFull::OneAtomVDWToNormalMatrix(std::vector<LIBMOL::AtomDict>& tAts, 
                                                    std::vector<LIBMOL::AtomDict>::iterator  tAt)
    {
        
        for (std::vector<int>::iterator iNB=tAt->neighbAtoms.begin(); 
                iNB !=tAt->neighbAtoms.end(); iNB++)
        {
            // One neighbor atom, with no bond connection with tAt 
            if (std::find(tAt->connAtoms.begin(), tAt->connAtoms.end(), *iNB)==tAt->connAtoms.end()
                    && tAt->seriNum < *iNB)
            {
                // TEMPO, need more detailed research 
                LIBMOL::REAL vdwFC = 100.0;
                LIBMOL::REAL r_nb = 0.0;
                int dim = (int)tAt->coords.size();
                LIBMOL::REAL r_comp [dim];
                for(int i=0; i < dim; i++)
                {
                    r_comp[i] = 0.0;
                }
                for (int i =0; i < dim; i++)
                {
                    r_comp[i]   = tAt->coords[i]
                                 -tAts[*iNB].coords[i];
                    r_nb += pow(r_comp[i],2.0);
                }
        
                // LIBMOL::distanceV(tAts[*tNB].coords, tAt->coords);
                LIBMOL::REAL r_d  = SetVDWContact(tAts, iNB, tAt);
                            
                if (r_nb < r_d)
                {
                    LIBMOL::REAL r_diff = r_nb-r_d;
                
                    for (int i =0; i < dim; i++)
                    {
                        LIBMOL::REAL tm1 = vdwFC*r_diff*r_comp[i]/r_nb;
                        int k1 = tAt->seriNum*dim + i;
                        for (int j=0; j < dim; j++)
                        {
                            LIBMOL::REAL tm2 = -r_diff*r_comp[j]/r_nb;
                            int k2 = *iNB*dim - j;
                            secDrivCart[k1][k2] = tm1*tm2;
                        }
                    }
                }
            }
        }
    }
    
    
      
    void extern SetChiraAndFirstDeriv(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::ChiralDict>::iterator tCh,
                                      LIBMOL::REAL  & volume,  LIBMOL::REAL  * df_dx_c)
    {
        int n_ato = (int)tCh->atoms.size();
        if (n_ato)
        {
            int dim = (int)tAts[tCh->atoms[0]].coords.size();
            int i1,i2;
            int i_ato1, i_ato2, i_ato3, i_ato4;

            int d4 = 4*dim;
 
            // Three vectors
            
            LIBMOL::REAL *a = new LIBMOL::REAL [dim];  
            LIBMOL::REAL *b = new LIBMOL::REAL [dim];
            LIBMOL::REAL *c = new LIBMOL::REAL [dim];
            
            // some of variables could be removed later, currently using old
            // new codes both 
            i_ato1 = tCh->atoms[0];
            i_ato2 = tCh->atoms[1];
            i_ato3 = tCh->atoms[2];
            i_ato4 = tCh->atoms[3];

            for (i1 = 0; i1 < dim; i1++)
            {
                a[i1] = tAts[i_ato2].coords[i1]-tAts[i_ato1].coords[i1];
                b[i1] = tAts[i_ato3].coords[i1]-tAts[i_ato1].coords[i1];
                c[i1] = tAts[i_ato4].coords[i1]-tAts[i_ato1].coords[i1];
            }
            
            // The derivatives of these three vectors wrt atomic coordinates

            LIBMOL::REAL **da_dx = new LIBMOL::REAL *[d4];

            for (i1 = 0; i1 < d4; i1++)
            {
                da_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    da_dx [i1][i2] = 0.0;
                }                
            }

            LIBMOL::REAL **db_dx = new LIBMOL::REAL *[d4];
            for (i1 = 0; i1 < 4*dim; i1++)
            {
                db_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    db_dx [i1][i2] = 0.0;
                }     
            }

            
            LIBMOL::REAL **dc_dx = new LIBMOL::REAL *[4*dim];
            for (i1 = 0; i1 < 4*dim; i1++)
            {
                dc_dx [i1] = new LIBMOL::REAL [dim];
                for (i2 = 0; i2 < dim; i2++)
                {
                    dc_dx [i1][i2] = 0.0;
                }     
            }
            
            for (i1 = 0; i1 < dim; i1++)
            {
                da_dx[i1][i1]       = -1.0;
                da_dx[dim+i1][i1]   =  1.0;
                db_dx[i1][i1]       = -1.0;
                db_dx[2*dim+i1][i1] =  1.0;
                dc_dx[i1][i1]       = -1.0;
                dc_dx[3*dim+i1][i1] =1.0;
            }

            // three auxi vectors (a X b), (b X c) and (c X a)

            LIBMOL::REAL * ab = new LIBMOL::REAL [dim];
            LIBMOL::REAL * bc = new LIBMOL::REAL [dim];
            LIBMOL::REAL * ca = new LIBMOL::REAL [dim];
            for(i1 = 0; i1 < dim; i1++)
            {
                ab[i1] = 0.0;
                bc[i1] = 0.0;
                ca[i1] = 0.0;
            }
 
            LIBMOL::CrossP(a, b, ab);
            LIBMOL::CrossP(b, c, bc);
            LIBMOL::CrossP(c, a, ca);

            // Get the value of the chiral volume
  
            volume = LIBMOL::CalcDet(a,b,c);

            // Finally, get the first derivatives of the chiral volume wrt 
            // atomic coordinates

            // cout << "First derivatives of Chiral center " <<index_ch+1 << endl;
 
            for (i1 = 0; i1 < d4; i1++)
            {
                df_dx_c[i1] = LIBMOL::DotP(da_dx[i1], bc)
                             +LIBMOL::DotP(db_dx[i1], ca)
                             +LIBMOL::DotP(dc_dx[i1], ab);
                //  cout << "  df_dx_c["<< i1<<"]= " << df_dx_c[i1]
                //       << endl;
            }       
            
            // release memory

            delete [] a;
            a = 0;
            delete [] b;
            b = 0;
            delete [] c;
            c =0;

            for (i1 = 0; i1 < d4; i1++)
            {
                delete [] da_dx [i1];
                da_dx [i1] = 0;
                delete [] db_dx [i1];
                db_dx [i1] = 0;
                delete [] dc_dx [i1];
                dc_dx [i1] = 0;
            }
            delete [] da_dx;
            da_dx = 0;
            delete [] db_dx;
            db_dx = 0;
            delete [] dc_dx;
            dc_dx = 0;

            delete [] ab;
            ab = 0;
            delete [] bc;
            bc = 0;
            delete [] ca;
            ca = 0; 

        }
    }
}
