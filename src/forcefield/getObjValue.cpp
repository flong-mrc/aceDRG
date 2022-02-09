/* 
 * File:   getObjValues.cpp
 * Author: flong
 *
 * Created on February 18, 2013, 12:05 AM
 */

#include "getObjValues.h"

namespace FF
{
    
    GetObjValue::GetObjValue():workSpace(1),
                               lComp(100)
    {
    }
    
    GetObjValue::~GetObjValue()
    {
    }
    
    
    // getAll version 1
    LIBMOL::REAL GetObjValue::getAll(std::vector<LIBMOL::AtomDict>& tAts, 
                                     std::vector<LIBMOL::BondDict>& tBos, 
                                     std::vector<LIBMOL::AngleDict>& tAns, 
                                     std::vector<LIBMOL::TorsionDict>& tTos, 
                                     std::vector<LIBMOL::PlaneDict>& tPls, 
                                     std::vector<LIBMOL::ChiralDict>& tChs)
    {
        LIBMOL::REAL objAll=0.0;
        if (workSpace==1)
        {
            // Contributions from Bonds
            LIBMOL::REAL objBo =0.0;
            if(lComp==100 || lComp==110 || lComp==1)
            {
                for (std::vector<LIBMOL::BondDict>::iterator iBo=tBos.begin();
                    iBo !=tBos.end(); iBo++)
                {
                    objBo  += GetObjOneBond(tAts, iBo);
                }
            }
            // std::cout << "Done bond obj " <<  objBo << std::endl;
            
            LIBMOL::REAL objAng =0.0;
            if(lComp==100 || lComp==110 || lComp==2)
            {
                for (std::vector<LIBMOL::AngleDict>::iterator iAn=tAns.begin();
                     iAn != tAns.end(); iAn++)
                {
                    objAng += GetObjOneAngle(tAts, iAn);
                }
            }
            // std::cout << "Done angle obj " <<  objAng << std::endl;
           
            
            /*
            LIBMOL::REAL objTor =0.0;
            for(std::vector<LIBMOL::TorsionDict>::iterator iTo=tTos.begin();
                    iTo != tTos.end(); iTo++)
            {
                objTor+=GetObjOneTorsion(tAts, iTo);
            }
            */
            
            LIBMOL::REAL objPl =0.0;
            if(lComp==100 || lComp==4)
            {
                for (std::vector<LIBMOL::PlaneDict>::iterator iPl=tPls.begin();
                        iPl !=tPls.end(); iPl++)
                {
                    objPl+=GetObjOnePlane(tAts, iPl);
                }
            }
            //std::cout << "Done plane obj " <<  objPl << std::endl;
            
                 
            LIBMOL::REAL objCh =0.0;
            /*
            if(lComp==100 || lComp==5)
            {
                for (std::vector<LIBMOL::ChiralDict>::iterator iCh=tChs.begin();
                     iCh != tChs.end(); iCh++)
                {
                    objCh+=GetObjOneChiral(tAts, iCh);
                }
             }
            */
            
            objAll = objBo + objAng + objPl + objCh;
            
            // Get Vdw an Hydr-bond contributions
            // Assume the neighbor list has already been built
            
            LIBMOL::REAL objVdw =0.0;
            
            for (std::vector<LIBMOL::AtomDict>::iterator iAt=tAts.begin();
                iAt !=tAts.end(); iAt++)
            {
                objVdw+=GetObjOneAtomVdw(tAts, iAt);
            }
            // std::cout << "done vdw obj " << objVdw << std::endl;
            objAll += objVdw;
            
        }
        else if (workSpace==2)
        {
            
        }
        
        // std::cout << "total obj is " << objAll << std::endl;
        return  objAll;    
    
    }
    
    LIBMOL::REAL GetObjValue::getAll(std::vector<LIBMOL::AtomDict>    & tAts,
                                     std::vector<LIBMOL::BondDict>    & tBos,
                                     std::vector<LIBMOL::AngleDict>   & tAns,
                                     std::vector<LIBMOL::TorsionDict> & tTos,
                                     std::vector<LIBMOL::PlaneDict>   & tPls,
                                     std::vector<LIBMOL::ChiralDict>  & tChs,
                                     std::vector<LIBMOL::AtomDict>    & tAllAtoms)
    {
                
        LIBMOL::REAL objAll=0.0;
        if (workSpace==1)
        {
            // Contributions from Bonds
            /*
            for (std::vector<LIBMOL::BondDict>::iterator iB=tBos.begin();
                        iB !=tBos.end(); iB++)
            {
                std::cout << "obj \n"; 
                std::cout << iB->valueST    << std::endl
                          << iB->sigValueST << std::endl;
            }
             */
          
            LIBMOL::REAL objBo =0.0;
            
            if(lComp==100 || lComp==1)
            {
                
                
                for (std::vector<LIBMOL::BondDict>::iterator iBo=tBos.begin();
                    iBo !=tBos.end(); iBo++)
                {
                    //std::cout << iBo->valueST    << std::endl
                    //          << iBo->sigValueST << std::endl;
                    
                    objBo  += GetObjOneBond(tAllAtoms, iBo);   
                }
            }
            // std::cout << "Done bond obj " <<  objBo << std::endl;
          
            LIBMOL::REAL objAng =0.0;
            if(lComp==100 || lComp==2)
            {
                for (std::vector<LIBMOL::AngleDict>::iterator iAn=tAns.begin();
                     iAn != tAns.end(); iAn++)
                {
                    objAng += GetObjOneAngle(tAts, iAn);
                }
            }
            // std::cout << "Done angle obj " <<  objAng << std::endl;
          
            /*
            LIBMOL::REAL objTor =0.0;
            for(std::vector<LIBMOL::TorsionDict>::iterator iTo=tTos.begin();
                    iTo != tTos.end(); iTo++)
            {
                objTor+=GetObjOneTorsion(tAts, iTo);
            }
            */
            
            LIBMOL::REAL objPl =0.0;
            if(lComp==100 || lComp==4)
            {
                for (std::vector<LIBMOL::PlaneDict>::iterator iPl=tPls.begin();
                        iPl !=tPls.end(); iPl++)
                {
                    objPl+=GetObjOnePlane(tAts, iPl);
                }
            }
            //std::cout << "Done plane obj " <<  objPl << std::endl;
            
                 
            LIBMOL::REAL objCh =0.0;
            /*
            if(lComp==100 || lComp==5)
            {
                for (std::vector<LIBMOL::ChiralDict>::iterator iCh=tChs.begin();
                     iCh != tChs.end(); iCh++)
                {
                    objCh+=GetObjOneChiral(tAts, iCh);
                }
             }
            */
            
            objAll = objBo + objAng + objPl + objCh;
            
            // Get Vdw an Hydr-bond contributions
            // Assume the neighbor list has already been built
            
            LIBMOL::REAL objVdw =0.0;
            
            for (std::vector<LIBMOL::AtomDict>::iterator iAt=tAts.begin();
                iAt !=tAts.end(); iAt++)
            {
                objVdw+=GetObjOneAtomVdw(tAts, iAt);
            }
            // std::cout << "done vdw obj " << objVdw << std::endl;
            objAll += objVdw;
            
        }
        else if (workSpace==2)
        {
            
        }
        
        
        
        // std::cout << "total obj is " << objAll << std::endl;
        
        return  objAll;    
    
        
    }
    
    
    // getAll version 2 
    LIBMOL::REAL GetObjValue::getAll(LIBMOL::REAL* vars, 
                                     std::vector<LIBMOL::AtomDict> & tAts,
                                     std::vector<LIBMOL::BondDict> & tBos,
                                     std::vector<LIBMOL::AngleDict> & tAns,
                                     std::vector<LIBMOL::TorsionDict> & tTos,
                                     std::vector<LIBMOL::PlaneDict>  & tPls,
                                     std::vector<LIBMOL::ChiralDict> & tChs)
    {
        int tSize = (int)tAts.size(); 
        for (int i=0; i <tSize;  i++ )
        {
            int tDim = (int)tAts[i].coords.size();
            for (int j=0; j <tDim; i++)
            {
                int k = i*tDim+j;
                tAts[i].coords[j] = vars[k];
            }
        }
        
        return getAll(tAts, tBos, tAns, tTos, tPls, tChs);
    }
    
    
    // getAll version 3. 
    LIBMOL::REAL GetObjValue::getAll(std::vector<LIBMOL::REAL> & vars, 
                                     std::vector<LIBMOL::AtomDict> & tAts,
                                     std::vector<LIBMOL::BondDict> & tBos,
                                     std::vector<LIBMOL::AngleDict> & tAns,
                                     std::vector<LIBMOL::TorsionDict> & tTos,
                                     std::vector<LIBMOL::PlaneDict>  & tPls,
                                     std::vector<LIBMOL::ChiralDict> & tChs)
    {
        int tSize = (int)tAts.size(); 
        for (int i=0; i <tSize;  i++ )
        {
            int tDim = (int)tAts[i].coords.size();
            for (int j=0; j <tDim; j++)
            {
                int k = i*tDim+j;
                tAts[i].coords[j] = vars[k];
            }
        }
        
        return getAll(tAts, tBos, tAns, tTos, tPls, tChs);
    }
    
    LIBMOL::REAL GetObjValue::GetObjOneBond(std::vector<LIBMOL::AtomDict> & tAts,
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
                // std::cout << tFor << std::endl;
                
                /*
                if (fabs(tFor) >=0.2)
                {
                    std::cout << "atom " << tAts[tBo->atomsIdx[0]].id 
                              << " and " << tAts[tBo->atomsIdx[1]].id
                              << " added a bond force "
                              << tFor << std::endl;
                    std::cout << "Diff bond L " << bondL-tBo->valueST << std::endl;
                }
                */
                tAts[tBo->atomsIdx[0]].forces[i] += tFor;
                tAts[tBo->atomsIdx[1]].forces[i] -= tFor;
                       
            }
        }
         
        delete [] distComp;
        distComp = 0;
        
        LIBMOL::REAL objB=0.5*pow((bondL-tBo->valueST),2.0)/tBo->sigValueST;
        /*
        if (std::fabs(objB) > 0.5)
        {
            std::cout << "Contribution of the bond between atom " << tAts[tBo->atomsIdx[0]].id
                  << " and atom " << tAts[tBo->atomsIdx[1]].id << " is "
                  << objB << std::endl;
            std::cout << "bond length: " << bondL << std::endl;
            std::cout << "std value: " << tBo->valueST << std::endl;
            std::cout << "std dev: "   << tBo->sigValueST << std::endl;
        }
         */
        return objB;
          
    }
    
    LIBMOL::REAL GetObjValue::GetObjOneAngle(std::vector<LIBMOL::AtomDict>& tAts,
                                             std::vector<LIBMOL::AngleDict>::iterator tAn)
    {
        //std::cout << "angles between atom " << tAts[tAn->atoms[1]].id 
        //        << "--" << tAts[tAn->atoms[0]].id << "--" 
        //        <<  tAts[tAn->atoms[2]].id << std::endl;
        std::vector<LIBMOL::REAL> v1, v2;
        for (int i=0; i < (int)tAts[tAn->atoms[0]].coords.size(); i++)
        {
            v1.push_back(tAts[tAn->atoms[1]].coords[i]-tAts[tAn->atoms[0]].coords[i]);
            v2.push_back(tAts[tAn->atoms[2]].coords[i]-tAts[tAn->atoms[0]].coords[i]);
        }
        
        // LIBMOL::REAL aVal=LIBMOL::getAngle2V(v1, v2);
        
        // temporal variables

        int dim = (int)tAts[0].coords.size();
        int d2 = 2*dim;
        int d3 = 3*dim;

        int i, j, k;

        int i_atom1, i_atom2, i_atom3;
        int i_error=0.0;
        
        LIBMOL::REAL *a = new LIBMOL::REAL [dim];
        LIBMOL::REAL *b = new LIBMOL::REAL [dim];

        LIBMOL::REAL ba_value, d_ba, d_ba_sq, obj_value;

        LIBMOL::REAL R, dR_dx;

        LIBMOL::REAL *    dThetadx_t;  
        dThetadx_t = new LIBMOL::REAL [d3];
        for (i =0; i < d3; i++)
        {
            dThetadx_t[i] = 0.0;
        }

        LIBMOL::REAL ** da_dx;
        LIBMOL::REAL ** db_dx;
 
        da_dx =   new LIBMOL::REAL * [d3];
        for(j =0; j < d3; j++)
        {
            da_dx [j] = new LIBMOL::REAL [dim];
        }
  
        db_dx =   new LIBMOL::REAL * [d3];
        for(j =0; j < d3; j++)
        {
            db_dx [j] = new LIBMOL::REAL [dim];
        }

        LIBMOL::REAL leng_a, leng_b;
        LIBMOL::REAL * d_leng_a;
        LIBMOL::REAL * d_leng_b;

        d_leng_a = new LIBMOL::REAL [d3];
        d_leng_b = new LIBMOL::REAL [d3];

        // Get vectors and the angle between them

        i_atom1 = tAn->atoms[1];
        i_atom2 = tAn->atoms[0];
        i_atom3 = tAn->atoms[2];
        
        for (i =0; i < dim; i++)
        {
            a[i] = tAts[i_atom1].coords[i] - tAts[i_atom2].coords[i];
            b[i] = tAts[i_atom3].coords[i] - tAts[i_atom2].coords[i];            
        }
     
        ba_value = LIBMOL::GetAngle(a,b); 
        //std::cout << "angle value on-fly "     << ba_value << std::endl;
        //std::cout << "angle dictionary value " << tAn->valueST << std::endl;
        //std::cout << "sigma of the angle "  << tAn->sigValueST << std::endl;
        
        d_ba     = ba_value - tAn->valueST;

        //std::cout << "diff is " << d_ba << std::endl;
        //std::cout << "force const is " << 1.0/tAn->sigValue;
        d_ba_sq  = pow(d_ba,2.0);

        obj_value = 0.50*d_ba_sq/tAn->sigValueST;
        /*
        if (fabs(obj_value) >=0.5)
        {    
            //std::cout << "angles between atom " << tAts[tAn->atoms[1]].id 
            //          << "--" << tAts[tAn->atoms[0]].id << "--" 
            //          <<  tAts[tAn->atoms[2]].id << std::endl;
            //std::cout << "angle value on-fly "     << ba_value << std::endl;
            //std::cout << "angle dictionary value " << tAn->valueST << std::endl;
            //std::cout << " addition to the obj value is " << obj_value << std::endl;
        }
        */
        // get derivs

        leng_a = LIBMOL::length_v(a, dim);
        leng_b = LIBMOL::length_v(b, dim);  

        if(leng_a < 1.0e-8)
        {
            std::cout << "atom " << tAts[i_atom2].id << " and atom " 
                      << tAts[i_atom1].id << " overlaped ! Wrong " 
                      << std::endl;
            i_error = 1;
        } 
        
        if(leng_b < 1.0e-8)
        {
            std::cout << "atom " << tAts[i_atom2].id 
                      << " and atom " << tAts[i_atom3].id
                      << " overlapped ! Wrong " << std::endl;
            i_error = 1;
    
        }

        if(!i_error)
        {
            R  = LIBMOL::DotP(a,b)/(leng_a*leng_b);
            //std::cout << "R " << R << std::endl;
            
            for (j =0; j < d3; j++)
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
                d_leng_b[j+d2] = -d_leng_b[j+dim];          
            }

            for (j =0; j < d3; j++)
            {
                dR_dx = (LIBMOL::DotP(da_dx[j],b)+LIBMOL::DotP(db_dx[j],a))/(leng_a*leng_b)
                         -LIBMOL::DotP(a,b)*(leng_b*d_leng_a[j]+leng_a*d_leng_b[j])
	                 /pow(leng_a*leng_b,2.0);
     
                if(fabs(1-R*R) > 1.0e-16)
                {
                    dThetadx_t[j] = -dR_dx/sqrt(1-R*R);
                }
                else
                {
                    // maybe the angle is 180 degrees
                    //std::cout << "Cos(theta) = 1, atoms overlapped ! " 
                    //          << std::endl;
                    //i_error = 1;
                    dThetadx_t[j] = -dR_dx*100.0;
                }
            }

            if(!i_error)
            {   
                for (j =0; j < dim; j++) 
                {
                    //LIBMOL::REAL tm1 = d_ba*dThetadx_t[j]/tAn->sigValueST;
                    //LIBMOL::REAL tm2 = d_ba*dThetadx_t[dim+j]/tAn->sigValueST;
                    //LIBMOL::REAL tm3 = d_ba*dThetadx_t[2*dim+j]/tAn->sigValueST;
                    tAts[i_atom1].forces[j] -=
                     (d_ba*dThetadx_t[j]/tAn->sigValueST);

                    tAts[i_atom2].forces[j] -=
                     (d_ba*dThetadx_t[dim+j]/tAn->sigValueST);

                    tAts[i_atom3].forces[j] -=
                     (d_ba*dThetadx_t[2*dim+j]/tAn->sigValueST);
                    
                    /*
                    if (fabs(tm1) >=1.0e-4 || fabs(tm2) >=1.0e-4 || fabs(tm3) >=1.0e-4)
                    {
                        std::cout << "angle atom " << tAts[i_atom1].id << " added a force "
                                 << tm1 << std::endl;
                        std::cout << "angle atom " << tAts[i_atom2].id << " added a force "
                                 << tm2 << std::endl;
                        std::cout << "angle atom " << tAts[i_atom3].id << " added a force "
                                 << tm3 << std::endl;
                    }
                     */
                    
                }
            }

        }
        
        //  release momery 

        delete [] dThetadx_t;
        dThetadx_t = 0;
  
        delete [] a;
        a = 0;

        delete [] b;
        b = 0;

        for(j =0; j < d3; j++)
        {
            delete [] da_dx [j];
            da_dx [j] = 0;
        }
        delete [] da_dx;
        da_dx = 0;

        for(j =0; j < d3; j++)
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
 
        // return the objective value

        return obj_value;
    }
        
    
    LIBMOL::REAL GetObjValue::GetObjOneTorsion(std::vector<LIBMOL::AtomDict>& tAts,
                                      std::vector<LIBMOL::TorsionDict>::iterator tTo)
    {
        
        std::vector<LIBMOL::REAL> v1, v2, v3;
        for (int i=0; i < (int)tAts[tTo->atoms[0]].coords.size(); i++)
        {
            v1.push_back(tAts[tTo->atoms[1]].coords[i]-tAts[tTo->atoms[0]].coords[i]);
            v2.push_back(tAts[tTo->atoms[2]].coords[i]-tAts[tTo->atoms[1]].coords[i]);
            v3.push_back(tAts[tTo->atoms[3]].coords[i]-tAts[tTo->atoms[2]].coords[i]);
        }
        LIBMOL::REAL tVal=LIBMOL::getTorsion3V(v1, v2, v3);
        
        return 0.5*pow(tVal-tTo->valueST, 2.0)/tTo->sigValueST;
    }
    
    LIBMOL::REAL GetObjValue::GetObjOnePlane(std::vector<LIBMOL::AtomDict>& tAts, 
                                             std::vector<LIBMOL::PlaneDict>::iterator tPl)
    {
        // Calculate the average position of all atoms in the plane
        tPl->fConst =tPl->fConst/20;
        int dim = 3;
        int n_ato = (int)tPl->atoms.size();
        
        LIBMOL::REAL obj_pl =0.0;
        
        std::vector<int> t_ats;
        for (std::map<LIBMOL::ID, int>::iterator iAt=tPl->atoms.begin();
                  iAt !=tPl->atoms.end(); iAt++)
        {
            t_ats.push_back(iAt->second);
        }
        
        LIBMOL::REAL * X_ave = new LIBMOL::REAL[dim];      

        for (int i1 = 0; i1 < dim; i1++)
        {           
          X_ave[i1] = 0.0;
          for (int i2=0; i2 < n_ato; i2++)
          {
              X_ave [i1] += tAts[t_ats[i2]].coords[i1];
          }
          X_ave [i1] = X_ave [i1]/n_ato;
        }
        
        // Calculate the relative coords of atoms in the plane
        
        LIBMOL::REAL ** X_rel=new LIBMOL::REAL * [n_ato];    
      
        for(int i1 = 0; i1 < n_ato; i1++)
        {
          X_rel [i1] = new LIBMOL::REAL [dim];
          for (int i2 = 0; i2 < dim; i2++)
          {
              X_rel[i1][i2] =tAts[t_ats[i1]].coords[i2]-X_ave[i2];
          }
        }
        
        // Find the least square plane that fits all the atoms
        // and the derivatives of the plane coefficients rwt
        // atomic coordinates

        LIBMOL::REAL * Coefs = new LIBMOL::REAL [dim];
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
                                   
        
        // a. get the obj value contributed from the plane

        LIBMOL::REAL Coefs_D;
        Coefs_D = LIBMOL::DotP(X_ave, Coefs);

        LIBMOL::REAL t_pl;

        for (std::vector<int>::iterator iAt=t_ats.begin();
                iAt != t_ats.end(); iAt++)
        {
            t_pl = LIBMOL::DotP(tAts[*iAt].coords,Coefs)-Coefs_D;
            t_pl = pow(t_pl, 2.0)/tPl->fConst; 
            obj_pl += t_pl;
        }
                   
        // b. Get the deriv of the obj function contributed from the plane
        
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

        // c. Find the forces on each atoms in the plane

        LIBMOL::REAL t1_eq, t2_eq, t1;

        for(int i1=0; i1 < n_ato; i1++)
	{
        
          int i_ato = t_ats[i1];
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

          tAts[i_ato].forces[0] -=
                          (t1*dxrdx/tPl->fConst);
          /*
          if(fabs(tPl->fConst*t1*dxrdx)>0.2)
          {
             
              std::cout << tAts[i_ato].id << " added a plane force " 
                        <<  tPl->fConst*t1*dxrdx <<std::endl;
          }
           */
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

          tAts[i_ato].forces[1] -=
                           (t1*dxrdx/tPl->fConst);

          /*
          if(fabs(tPl->fConst*t1*dxrdx)>0.2)
          {
             
              std::cout << tAts[i_ato].id << " added a plane force " 
                        <<  tPl->fConst*t1*dxrdx <<std::endl;
          }
           */
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

          tAts[i_ato].forces[2] -= 
                           (t1*dxrdx/tPl->fConst);
          /*
          if(fabs(tPl->fConst*t1*dxrdx)>0.2)
          {
             std::cout << tAts[i_ato].id << " added a plane force " 
                       <<  tPl->fConst*t1*dxrdx <<std::endl;
          }
           */

	}
        
        tPl->fConst = tPl->fConst*20;
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
        
        // Return the obj value
        return obj_pl;
    }
    
    LIBMOL::REAL GetObjValue::GetObjOneChiral(std::vector<LIBMOL::AtomDict>& tAts, 
                                              std::vector<LIBMOL::ChiralDict>::iterator tCh)
    {

        int dim = (int)tAts[0].coords.size();
        int n_dim = 4*dim;
        int n_atoms = 4;
        int i, i_dim, i_ato;

        LIBMOL::REAL vol_ch   = 0.0;
        LIBMOL::REAL d_vol_ch = 0.0;
        LIBMOL::REAL weight   = 0.0;
        LIBMOL::REAL obj_ch   = 0.0;

        LIBMOL::REAL * firstDerivChiras = new LIBMOL::REAL  [n_dim];
        for(i_dim = 0; i_dim < n_dim; i_dim++)
        {
            firstDerivChiras[i_dim] = 0.0;
        }
        
        SetChiraAndFirstDeriv(tAts, tCh, vol_ch, firstDerivChiras);

        d_vol_ch = vol_ch - tCh->valueST;

        // get the obj value 

        obj_ch = 0.5*tCh->fConst*pow(d_vol_ch,2.0);
  
        // get the derivs of the obj 

        weight   = 2*tCh->fConst*d_vol_ch;
 
        for (i =0; i < n_atoms; i++)
        {
            i_ato  = tCh->atoms[i];
            i_dim = i*dim;      
            tAts[i_ato].forces[0] -= (weight*firstDerivChiras[i_dim]);
            i_dim++;
            tAts[i_ato].forces[1] -= (weight*firstDerivChiras[i_dim]);
            i_dim++;
            tAts[i_ato].forces[2] -= (weight*firstDerivChiras[i_dim]);
        }
        
        delete [] firstDerivChiras;
        firstDerivChiras = 0;

        return obj_ch;
    }
    
    LIBMOL::REAL GetObjValue::GetObjOneAtomVdw( std::vector<LIBMOL::AtomDict>& tAts,
                                                std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        //std::cout << "For atom " << tAt->id << ". Its serial num " << tAt->seriNum
        //          << std::endl; 
        LIBMOL::REAL vVDW =0.0;
        for (std::vector<int>::iterator iNB=tAt->neighbAtoms.begin(); 
                iNB !=tAt->neighbAtoms.end(); iNB++)
        {
            // One neighbor atom, with no bond connection with tAt 
            if (std::find(tAt->connAtoms.begin(), tAt->connAtoms.end(), *iNB)==tAt->connAtoms.end()
                    && tAt->seriNum < *iNB)
            {
                // Check which interaction should be activated Hydro-bond or VDW 
                //std::cout << "It has VDW interaction with atom " 
                //          << tAts[*iNB].id << std::endl;
                bool l_hyd_flag = CheckHYD(tAts, iNB, tAt);
                if (l_hyd_flag)
                {
                    vVDW+=CalcObjHYD(tAts, iNB, tAt);
                }
                else
                {
                    vVDW+=CalcObjVDWOnePair(tAts, iNB,tAt);
                    // std::cout << "of " << CalcObjVDWOnePair(tAts, iNB,tAt) << std::endl; 
                }
            }
        }
        // std::cout << "Total VDW for atom " << tAt->id << " is " << vVDW << std::endl;
        return vVDW;
    }
    
    bool GetObjValue::CheckHYD(std::vector<LIBMOL::AtomDict>& tAts,
                               std::vector<int>::iterator              tNB, 
                               std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        bool l_hyd = false;
        /*
        if (tAt->chemType =="H")
        {
            if ((int)tAt->connAtoms.size() == 2)
            {
                for(int i=0; i < tAt->connAtoms.size(); i++)
                {
                    int j= tAt->connAtoms[i];
                    if (tAts[j].id != tNB->id)
                    {
           
                    }
                }
            }
        }
         */
        return l_hyd;
    }
    
    LIBMOL::REAL GetObjValue::CalcObjHYD(std::vector<LIBMOL::AtomDict>& tAts, 
                                 std::vector<int>::iterator tNB, 
                                 std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        LIBMOL::REAL vHyd = 0.0;
        
        return vHyd;
    }
    
    
    
    LIBMOL::REAL GetObjValue::CalcObjVDWOnePair(std::vector<LIBMOL::AtomDict>& tAts,
                                   std::vector<int>::iterator tNB,
                                   std::vector<LIBMOL::AtomDict>::iterator tAt)
    {
        LIBMOL::REAL vVDWPair =0.0;
        LIBMOL::REAL vdwFC;
        // TEMPO, need more detailed research 
        if (tAt->chemType=="H" || tAts[*tNB].chemType=="H")
        {
            vdwFC = 30.0;
        }
        //else if (tAt->chemType=="H" || tAts[*tNB].chemType=="H")
        //{
        //    vdwFC = 50.0;
        //}
        else
        {
            vdwFC = 100.0;
        }
        
        LIBMOL::REAL r_nb = 0.0;
        int dim = (int)tAt->coords.size();
        std::vector<LIBMOL::REAL> r_comp(dim);
        for(int i=0; i < dim; i++)
        {
            r_comp[i] = 0.0;
        }
        for (int i =0; i < dim; i++)
        {
            r_comp[i]   = tAt->coords[i]
                         -tAts[*tNB].coords[i];
            r_nb += pow(r_comp[i],2.0);
        }
        
        // LIBMOL::distanceV(tAts[*tNB].coords, tAt->coords);
        LIBMOL::REAL r_d  = SetVDWContact(tAts, tNB, tAt);
        
        if (r_nb < r_d)
        {
            // std::cout << "distance between " << tAt->id << " and " << tAts[*tNB].id 
            //          << " is " << r_nb << std::endl;
            //std::cout << "Its vdw contact dist is " << r_d << std::endl; 
            LIBMOL::REAL r_diff = r_nb-r_d;
            vVDWPair = 0.50*vdwFC*pow(r_diff,2.0);
            for (int i =0; i < dim; i++)
            {
                LIBMOL::REAL tm = vdwFC*r_diff*r_comp[i]/r_nb;

                tAt->forces[i]  -= tm;
                tAts[*tNB].forces[i] += tm;
                /*
                if (fabs(tm) >0.2)
                {
                    std::cout << tAt->id << " and " << tAts[*tNB].id 
                            << " add vdw force " << tm 
                            << " The contact diff is " << r_diff << std::endl;
                }
                 */
            }
            //std::cout << "r_diff is " << r_diff << std::endl;
            //std::cout << "Contribution to VDW obj is " << vVDWPair << std::endl;
        }
        
        
        
        return vVDWPair;
    }
    
    LIBMOL::REAL GetObjValue::SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
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
    
    LIBMOL::REAL extern SetVDWContact(std::vector<LIBMOL::AtomDict>& tAts,
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
}
