/* 
 * File:   TransCoord.h
 * Author: flong
 *
 * Created on February 26, 2013, 4:23 PM
*/

#include "TransCoord.h"
#include "atomsTree.h"
#include "GlobOpt.h"

namespace LIBMOL
{
    TransCoords::TransCoords()
    {
    }
    
    TransCoords::~TransCoords()
    {
        
    }

    void TransCoords::generateCoordTorsToCart(std::vector<AtomDict>& tAtoms,
                                     std::vector<BondDict>    & tBonds,
                                     std::vector<AngleDict>   & tAngles,
                                     std::vector<TorsionDict> & tTorsions,
                                     std::map<ID, std::vector<RingDict> >    & tRings,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChs)
    {
        
        std::vector<RingDict> tRingsV;
        for (std::map<ID, std::vector<RingDict> >::iterator iMR=tRings.begin();
                    iMR != tRings.end(); iMR++)
        {
            for (std::vector<RingDict>::iterator iR=iMR->second.begin();
                        iR != iMR->second.end(); iR++)
            {
                tRingsV.push_back(*iR);
            }
        }
        
        LIBMOL::buildAtomTree tTool;
            
        tTool.buildTree(tAtoms, tBonds, tAngles, tTorsions, 
                        tRingsV, tPlas, tChs);
        
        //for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
        //        iA !=tAtoms.end(); iA++)
        //{
        //    std::cout << "atom " << iA->id << " element " << iA->chemType 
        //              << " seriNum " << iA->seriNum << std::endl;
        //}
        
        //tTool.setTreeAtomValues(tAtoms, tBonds, tAngs, tTorsions, 
        //                        tRings, tPlas, tChs);
        
        for (std::map<int, std::map<int, std::vector<int> > >::iterator iBr=tTool.branches.begin();
                iBr !=tTool.branches.end(); iBr++)
        {
            for(std::map<int, std::vector<int> >::iterator iSB=iBr->second.begin();
                    iSB !=iBr->second.end(); iSB++)
            {
                for (std::vector<int>::iterator iSSB=iSB->second.begin();
                        iSSB != iSB->second.end(); iSSB++)
                {
                    // Build a branch of the tree beginning with 3 atoms      
                    branchGrowthTorsTorCart(tAtoms, iBr->first, iSB->first, *iSSB);
                }
            }
        } 
    }
    
    void TransCoords::branchGrowthTorsTorCart(std::vector<AtomDict>& tAtoms,
                                              int sAtom1, int sAtom2, int sAtom3)
    {
        
        //std::cout << "Branch growth from " << tAtoms[sAtom1].id
        //          << " and atom " << tAtoms[sAtom2].id 
        //          << " and atom " << tAtoms[sAtom3].id << std::endl;
        
        const int NSTLIM2 = 3;
        const int NSTLIM  = 10000;
        int dim = (int)tAtoms[0].coords.size();
        
        // Just for debug this time
        if (dim !=3)
        {
            std::cout << "Coord dim does not equal to 3, check !" << std::endl;
            exit(1);
        }
        
        int i, j;
       
        int i_atom, i_next;
        int i_conn, n_conn;
        int n_stack;
        int * iatom_stack  = new int [NSTLIM];
        int * i_conn_stack = new int [NSTLIM];
        

        REAL *** aMatStack = new REAL ** [NSTLIM];
        for ( i =0; i < NSTLIM; i++)
        {
            aMatStack[i] = new REAL *[NSTLIM2];
            for(j =0; j < NSTLIM2; j++)
            {
                aMatStack[i][j] = new REAL [NSTLIM2];
            }
        }
        
        REAL xyz_stack[NSTLIM][dim];
        
        REAL ** rotate_tem;
        rotate_tem = new REAL * [dim];
        for (i =0; i < dim; i++)
        {
            rotate_tem[i] = new REAL [dim];
            for (j = 0; j < dim; j++)
            {
                rotate_tem[i][j] =0.0;
            }
        }

        REAL **A = new REAL * [dim];
        for (i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        REAL *v3    = new REAL [dim];
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sAtom2].coords[i] - tAtoms[sAtom1].coords[i];
            v2[i] = tAtoms[sAtom3].coords[i] - tAtoms[sAtom2].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  
        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start climbing the tree from the seleted atom
  
        i_atom   = sAtom3;
        
        //cout << "i_atom is " << i_atom+1 << endl; 
        
        //for(i =0; i < dim; i++)
        //  {
        //    atoms[i_atom].coords[i] =atoms[i_atom].coords[i];     
        //  }

        i_conn   = 0;
        n_stack  = 0;
  
        n_conn = (int)tAtoms[i_atom].tree["children"].size();

        //std::cout << " The first atom at the branch is atom " << tAtoms[i_atom].id << std::endl;
        //std::cout << " There are " << n_conn << " atoms forwardly linked this atom" << std::endl;
 

        while((n_conn > 0 && n_stack >= 0)|| (n_conn== 0 && n_stack > 0) )
        {
            if ( n_conn > 1 &&  n_stack >= 0 ) 
            {
                i_next = tAtoms[i_atom].tree["children"][i_conn];
                //std::cout << " The next atom is " << tAtoms[i_next].id << std::endl;
	        // cin.get();
	        //   if(tAllAtoms[i_next].occupancy < 0.001)  && tAllAtoms[i_next].type != "R")
	   
                PolToCart(tAtoms,i_atom, i_next, A);                  
	                                     
                if (i_conn < n_conn )
                {
                    if(n_stack > NSTLIM)
                    {
                        std::cout << " Error in member function generateCoordCartToTors " 
                                  << std::endl << " Number in the stack is :  " 
                                  << n_stack  << std::endl;
                        exit(1);
                    }                                 

                    // Keep the junction point info in a stack
    
                    i_conn               = i_conn  + 1;
	            iatom_stack[n_stack] = i_atom;
             
                    i_conn_stack[n_stack] = i_conn;

	            // cout << " The infor of the atom " <<  i_atom +1 
	            //     << " has been stored in the stack " << n_stack +1 << endl;
             
                    NBD_MCOPY(A, aMatStack[n_stack], dim);
                             
                    for ( i =0; i < dim; i++)
                    {
                        xyz_stack[n_stack][i] = tAtoms[i_atom].coords[i];
                    }
	      
                    n_stack   = n_stack + 1;

		    //  cout << " n_stack :   " << n_stack << endl;  
                }

                CHANGE_AMAT(tAtoms,i_next, A);

                // Set the next atom to be the current atom in new loop

                i_atom     = i_next;
                i_conn     = 0;

                n_conn     = (int)tAtoms[i_atom].tree["children"].size();          
         
	        //  std::cout << " new i_atom: " << i_atom+1 << endl;
	        //  std::cout << " new n_conn: " << n_conn << endl;
	  
            }
            else if (n_conn == 1)
            { 
                i_next = tAtoms[i_atom].tree["children"][i_conn];

	        //  cout << " The next atom is " << i_next+1 << endl;
               
                PolToCart(tAtoms, i_atom, i_next, A);
	                                    //  strcpy(tAllAtoms[i_next].type, "R");
                CHANGE_AMAT(tAtoms, i_next, A);


                // Set the next atom to be the current atom in new loop

                i_atom     = i_next;
                i_conn     = 0;

                n_conn = tAtoms[i_atom].tree["children"].size();          
         
	        //  cout << " new i_atom: " << i_atom+1 << endl;
	        //  cout << " new n_conn: " << n_conn << endl;
            }
            else if (n_conn == 0 && n_stack > 0)
            {
                do
                {
                    n_stack = n_stack-1;
	            //   cout << "The stack level is " << n_stack << endl;
                    if(n_stack < 0)
                    {
                        break;
                    }
                    i_atom = iatom_stack[n_stack];

	            //  std::cout << " restart from i_atom " << i_atom << std::endl;
         
                    i_conn = i_conn_stack[n_stack];

	            //  cout << " restart1 from the " << i_conn+1 
	            //       << " node that linked to the above atom " << endl;

                    NBD_MCOPY(aMatStack[n_stack], A, dim);

	            //for( i =0; i < dim; i++)
	            //{
                    //  xyz[i] = xyz_stack[n_stack][i];
	            //}

                    n_conn = (int)tAtoms[i_atom].tree["children"].size();

	            //    cout << "restart n_conn " << n_conn << endl;
	            //    cout << " continue      " << endl;             
	            //    cin.get();
                
                }while(i_conn >= n_conn && n_stack >= 0);
	  
            }

        }
        //  cout << "New coordiantes after angle-to-coordinates transfer are " << endl;

        // for (i =0; i < numAtoms; i++)
        //   {
        //     cout << "The " <<i+1<<"th atom is at" << endl;
        //     for (j=0; j < dim; j++)
        // 	{
        //        cout << tAllAtoms[i].coords[j] << "\t";
        //      }
        //    cout << endl;
        //  }
  
        // cout << " Continue ? " << endl;
        // cin.get();

        // Release memory and null point variables

        delete [] v1;
        v1 = 0;

        delete v2;
        v2 = 0;

        delete [] v3;
        v3 = 0;
  
        delete [] tmp_v;
        tmp_v = 0;
  
        delete [] X_new;
        X_new = 0;

        delete [] Y_new;
        Y_new = 0;

        delete [] Z_new;
        Z_new = 0;
  
        for ( i =0; i < dim; i++)
        {
            delete [] rotate_tem[i];
            rotate_tem[i] =0;
        }

        delete [] rotate_tem;
        rotate_tem = 0;
  
        for (i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete [] A;
        A = 0;

        delete [] iatom_stack;
        iatom_stack  = 0;

        delete [] i_conn_stack;
        i_conn_stack = 0;

        for ( i =0; i < NSTLIM; i++)
        {
            for(j =0; j < NSTLIM2; j++)
            {
                delete []    aMatStack[i][j];
                aMatStack[i][j] = 0;
            }
            delete [] aMatStack[i];
            aMatStack[i] =0;
        }
        delete [] aMatStack;
        aMatStack = 0; 
    }
    
    void TransCoords::CHANGE_AMAT(std::vector<AtomDict>& tAtoms, 
                                  int i_next, 
                                  REAL** aMat)
    {
        int dim = (int)tAtoms[i_next].coords.size();
          
        int i;
        
        REAL th_i, phi_i;
        REAL cos_t,sin_t, cos_f, sin_f;
        REAL t, d_t, f, d_f;

        REAL ** Y = new REAL * [dim];
        REAL ** Z = new REAL * [dim];
        
        for (i =0; i < dim; i++)
        {
            Y[i] = new REAL [dim];
            Z[i] = new REAL [dim];
        }
        
        th_i   = PI-tAtoms[i_next].treeAngle;
        phi_i  = tAtoms[i_next].treeTorsion;
        
        
        //--------------------------
        // cout << "i_next " << i_next+1 << endl;
        // cout << "atoms[i_next].treeTorsAng " << atoms[i_next].treeTorsAng << endl;
 
        t      = th_i;
        f      = phi_i;

        d_t    = t;
        d_f    = f;
        cos_t  = cos(d_t);
        sin_t  = sin(d_t);
        cos_f  = cos(d_f);
        sin_f  = sin(d_f);

        //---------------------------

        Y[0][0] =  cos_t;
        Y[0][1] = -sin_t;
        Y[0][2] =  0.0;
        Y[1][0] =  sin_t*cos_f;
        Y[1][1] =  cos_t*cos_f;
        Y[1][2] = -sin_f;
        Y[2][0] =  sin_t*sin_f;
        Y[2][1] =  cos_t*sin_f;
        Y[2][2] =  cos_f;
 
        NBD_MATMLT(aMat,Y,Z, dim);
        NBD_MCOPY(Z,aMat, dim);

        for (i = 0; i < dim; i++)
        {
            delete [] Y[i];
            Y[i] = 0;
        }
        delete  [] Y;
        Y = 0;

        for (i =0; i < dim; i++)
        {
            delete [] Z[i];
            Z[i] = 0;
        }
        delete [] Z;
        Z  = 0;
        
    }
    
    
    void TransCoords::generateCoordTorsToCart2(std::vector<AtomDict>& tAtoms,
                                     std::vector<BondDict>    & tBonds,
                                     std::vector<AngleDict>   & tAngles,
                                     std::vector<TorsionDict> & tTorsions,
                                     std::vector<RingDict>    & tRings,
                                     std::vector<PlaneDict>   & tPlas,
                                     std::vector<ChiralDict>  & tChs)
    {
        
        
        LIBMOL::buildAtomTree tTool;
            
        tTool.buildTree(tAtoms, tBonds, tAngles, tTorsions, 
                        tRings, tPlas, tChs);
        
        //for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
        //        iA !=tAtoms.end(); iA++)
        //{
        //    std::cout << "atom " << iA->id << " element " << iA->chemType 
        //              << " seriNum " << iA->seriNum << std::endl;
        //}
        
        //tTool.setTreeAtomValues(tAtoms, tBonds, tAngs, tTorsions, 
        //                        tRings, tPlas, tChs);
        
        for (std::map<int, std::map<int, std::vector<int> > >::iterator iBr=tTool.branches.begin();
                iBr !=tTool.branches.end(); iBr++)
        {
            for(std::map<int, std::vector<int> >::iterator iSB=iBr->second.begin();
                    iSB !=iBr->second.end(); iSB++)
            {
                for (std::vector<int>::iterator iSSB=iSB->second.begin();
                        iSSB != iSB->second.end(); iSSB++)
                {
                    // Build a branch of the tree beginning with 3 atoms   
                    //std::cout << "Branch 1st atom " << tAtoms[iBr->first].id 
                    //          << " 2nd " << tAtoms[iSB->first].id
                    //          << " 3rd " << tAtoms[*iSSB].id << std::endl;
                    
                    if ((int)tAtoms[*iSSB].tree["children"].size() > 0)
                    {
                        if(std::find(doneList.begin(), doneList.end(), iBr->first)==doneList.end())
                        {
                            doneList.push_back(iBr->first);
                        }
                        if(std::find(doneList.begin(), doneList.end(), iSB->first)==doneList.end())
                        {
                            doneList.push_back(iSB->first);
                        }
                        if(std::find(doneList.begin(), doneList.end(), *iSSB)==doneList.end())
                        {
                            doneList.push_back(*iSSB);
                        }
                        branchGrowthTorsToCart2(tAtoms, iBr->first, iSB->first, *iSSB,
                                                tBonds, tAngles, tTorsions, 
                                                tRings, tPlas,   tChs);
                    }
                }
            }
        } 
    }
    
    void TransCoords::branchGrowthTorsToCart2(std::vector<AtomDict>& tAtoms,
                                              int sAtom1, int sAtom2, int sAtom3,
                                              std::vector<BondDict>    & tBonds,
                                              std::vector<AngleDict>   & tAngles,
                                              std::vector<TorsionDict> & tTorsions,
                                              std::vector<RingDict>    & tRings,
                                              std::vector<PlaneDict>   & tPlas,
                                              std::vector<ChiralDict>  & tChs)
    {
        
        
        //std::cout << "Branch growth from " << tAtoms[sAtom1].id
        //        << " and atom " << tAtoms[sAtom2].id 
        //        << " and atom " << tAtoms[sAtom3].id << std::endl;
        
        const int NSTLIM2 = 3;
        const int NSTLIM =  10000;
        int dim = (int)tAtoms[0].coords.size();
        
        // Just for debug this time
        if (dim !=3)
        {
            std::cout << "Coord dim does not equal to 3, check !" << std::endl;
            exit(1);
        }
        
        int i, j;
       
        int i_atom, i_next;
        int i_conn, n_conn;
        int n_stack;
        int * iatom_stack  = new int [NSTLIM];
        int * i_conn_stack = new int [NSTLIM];
        

        REAL *** aMatStack = new REAL ** [NSTLIM];
        for ( i =0; i < NSTLIM; i++)
        {
            aMatStack[i] = new REAL *[NSTLIM2];
            for(j =0; j < NSTLIM2; j++)
            {
                aMatStack[i][j] = new REAL [NSTLIM2];
            }
        }
        
        REAL xyz_stack[NSTLIM][dim];
        
        REAL ** rotate_tem;
        rotate_tem = new REAL * [dim];
        for (i =0; i < dim; i++)
        {
            rotate_tem[i] = new REAL [dim];
            for (j = 0; j < dim; j++)
            {
                rotate_tem[i][j] =0.0;
            }
        }

        REAL **A = new REAL * [dim];
        for (i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        REAL *v3    = new REAL [dim];
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sAtom2].coords[i] - tAtoms[sAtom1].coords[i];
            v2[i] = tAtoms[sAtom3].coords[i] - tAtoms[sAtom2].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  
        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start climbing the tree from the seleted atom
  
        i_atom   = sAtom3;
        
        //cout << "i_atom is " << i_atom+1 << endl; 
        
        //for(i =0; i < dim; i++)
        //  {
        //    atoms[i_atom].coords[i] =atoms[i_atom].coords[i];     
        //  }

        i_conn   = 0;
        n_stack  = 0;
  
        n_conn = (int)tAtoms[i_atom].tree["children"].size();

        // std::cout << " The first atom at the branch is atom " << tAtoms[i_atom].id << std::endl;
        // std::cout << " There are " << n_conn << " atoms forwardly linked this atom" << std::endl;
 
        std::vector<AtomDict> dAtms;
        std::vector<int>   aStartSet;

        std::map<std::string, std::vector<int> >::iterator iFind;
        
        while((n_conn > 0 && n_stack >= 0)|| (n_conn== 0 && n_stack > 0) )
        {
            if ( n_conn > 1 &&  n_stack >= 0 ) 
            {
                int i_prev   = tAtoms[i_atom].tree["parent"][0];
                int i1 = tAtoms[i_atom].tree["children"][i_conn];
                int i2 = -1;
                //std::cout << "new atom id " << tAtoms[i1].id << std::endl;
                
                if (tAtoms[i1].chemType=="H")
                {
                    //std::cout << "it is " << tAtoms[i1].id << std::endl;
                    
                    // grow other atom first if exist
                    for (std::vector<int>::iterator iA=tAtoms[i_atom].connAtoms.begin();
                            iA !=tAtoms[i_atom].connAtoms.end(); iA++)
                    {
                        //std::cout << "i2 =" << tAtoms[*iA].id << std::endl;
                        if (std::find(doneList.begin(), doneList.end(), *iA)==doneList.end()
                             && *iA != i1 && *iA !=i_prev && tAtoms[*iA].chemType !="H")
                        {
                            i2=*iA;
                            //std::cout << "i2 =" << tAtoms[i2].id << std::endl;
                            break;
                        }
                    }
                }
                if (i2 != -1)
                {
                    i_next = i2;
                }
                else
                {
                    i_next = i1;
                }
                //std::cout << "Final pick " << tAtoms[i_next].id << std::endl; 
                std::vector<int> aSSet;
                int i_pprev;
                if (tAtoms[i_atom].id !=tAtoms[sAtom3].id)
                {
                    
                    i_pprev  = tAtoms[i_prev].tree["parent"][0];
                }
                else
                {
                    i_pprev  = sAtom1;
                    
                }
                aSSet.push_back(i_pprev);
                aSSet.push_back(i_prev);
                aSSet.push_back(i_atom);
                
                //std::cout << " The grand parent atom " << tAtoms[i_pprev].id << std::endl;
                //std::cout << " The parent atom " << tAtoms[i_prev].id << std::endl;
                //std::cout << " The current atom " << tAtoms[i_atom].id << std::endl;
                //std::cout << " The next atom is " << tAtoms[i_next].id << std::endl;
                
	        // cin.get();
	        //   if(tAllAtoms[i_next].occupancy < 0.001)  && tAllAtoms[i_next].type != "R")
	   
                if (std::find(doneList.begin(), doneList.end(), i_next) ==doneList.end())
                {
                    REAL tBo;   //   = tAtoms[i_next].treeBond;
                    REAL tTha;  //  = tAtoms[i_next].treeAngle; 
                    REAL tPhi  = tAtoms[i_next].treeTorsion;
                    //if (tAtoms[i_atom].id !=tAtoms[sAtom3].id)
                    //{
                    growOneAtom(aSSet, tAtoms, tBonds, tAngles, tTorsions, i_next);
                    
                    tPhi = tTorsions[getTorsion(tTorsions, aSSet[0], aSSet[1], aSSet[2], i_next)].value*PI180;
                    //}
                    //else
                    //{
                    //    PolToCart(tAtoms,i_atom, i_next, tBo, tTha, tPhi,  A);
                    //}
                    doneList.push_back(i_next);
                    // std::cout <<  tAtoms[i_next].id << " done0 " << std::endl;
                    // std::cout << " atom " << tAtoms[i_atom].id 
                    //          << " chiral idx " <<tAtoms[i_atom].chiralIdx << std::endl;
                    
                    // grow other atoms linked to i_atom (like i_next)
                    // int i_chir = tAtoms[i_atom].chiralIdx;
                
                    if (tAtoms[i_atom].chiralIdx !=0 && (int)tChs[tAtoms[i_atom].inChirals[0]].mutTable.size() !=0)
                    {
                        int aCh = tAtoms[i_atom].inChirals[0];
                        bool lS= false;
                        std::vector<int> t1, rTable;
                        for (std::vector<int>::iterator irA=tChs[aCh].mutTable[i_prev].begin();
                               irA!=tChs[aCh].mutTable[i_prev].end(); irA++)
                        {
                            if (lS)
                            {
                                rTable.push_back(*irA);
                            }
                            else if (*irA==i_next)
                            {
                                lS=true;
                            }
                            else
                            {
                                t1.push_back(*irA);
                            }
                        }
                        for (std::vector<int>::iterator irA2=t1.begin();
                                irA2 !=t1.end(); irA2++)
                        {
                            rTable.push_back(*irA2);
                        }
                        
                        int iL=1;
                        
                        for (std::vector<int>::iterator iNodeA=rTable.begin();
                                iNodeA !=rTable.end(); iNodeA++)
                        {
                            //std::cout << "iNodeA " << tAtoms[*iNodeA].id << std::endl;
                            
                            int bIdx    = getBond(tBonds, i_atom, *iNodeA);
                            tBo         = tBonds[bIdx].valueST;
                            int aIdx    = getAngle(tAngles, i_atom, i_prev, *iNodeA);
                            tTha        = tAngles[aIdx].valueST;
                            //std::cout << "iL= " << iL << std::endl
                            //          << " tPhi" << tPhi << std::endl
                            //          << " deltal " << 2.0*PI/3.0 << std::endl;
                            REAL tPhiI  = tPhi +iL*2.0*PI/3.0;
                            if (std::find(doneList.begin(), doneList.end(), *iNodeA) ==doneList.end())
                            {
                                // PolToCart(tAtoms,i_atom, *iNodeA, tBo, tTha, tPhiI,  A);
                                growOneAtom(aSSet, tAtoms, tBo, tTha, tPhiI, *iNodeA);
                                doneList.push_back(*iNodeA);
                                // std::cout <<  tAtoms[*iNodeA].id << " done1 " << std::endl;
                            }
                            iL++;
                        }
                        
                    }
                    else
                    {
                        int iL=1;
                        for (std::vector<int>::iterator iCA=tAtoms[i_atom].connAtoms.begin();
                                iCA != tAtoms[i_atom].connAtoms.end(); iCA++)
                        {
                            if (*iCA !=i_next && *iCA !=i_prev 
                                && std::find(doneList.begin(), doneList.end(), *iCA) ==doneList.end())
                            {
                                tBo         = tAtoms[*iCA].treeBond;
                                tTha        = tAtoms[*iCA].treeAngle;
                                if (tAtoms[i_atom].bondingIdx==3)
                                {    
                                    REAL tPhiI  = tPhi +iL*2.0*PI/3.0;
                                    growOneAtom(aSSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    //std::cout <<  tAtoms[*iCA].id << " done21 " << std::endl;
                                    iL++;
                                }
                                else if (tAtoms[i_atom].bondingIdx==2)
                                {
                                    REAL tPhiI  = tPhi + PI;   
                                    growOneAtom(aSSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    //std::cout <<  tAtoms[*iCA].id << " done22 " << std::endl;
                                }
                            }
                        }    
                    }
                    
                    
                    optSubSystem(doneList, tAtoms, tBonds, tAngles,
                                 tTorsions, tRings, tPlas, tChs, false);
                    //std::cout << "The sub-system is optimized \n";
                    
                    // Ring growth
                    
                    if ((int)tAtoms[i_next].inRings.size() !=0 )
                    {   
                        
                        int iTurn =0;
                        for (std::vector<int>::iterator iRi=tAtoms[i_next].inRings.begin();
                                iRi !=tAtoms[i_next].inRings.end(); iRi++)
                        {
                            if (! tRings[*iRi].isPlanar)
                            {
                                //std::cout << "Grow the ring contain " << tAtoms[i_next].id
                                //          << std::endl;
                                aStartSet.clear();
                                if (i_atom == sAtom3)
                                {
                                    aStartSet.push_back(sAtom2);
                                }
                                else
                                {
                                    int iGP = tAtoms[i_atom].tree["parent"][0];
                                    aStartSet.push_back(iGP);
                                }
                                aStartSet.push_back(i_atom);
                                aStartSet.push_back(i_next);
                                //std::cout << "S1 set " << aStartSet[0] << "\t"
                                //          << aStartSet[1] << "\t" << aStartSet[2] << std::endl;
                                ringBuilder(aStartSet, tRings[*iRi], tAtoms, tBonds, 
                                            tAngles, tTorsions, tRings, tPlas, tChs, iTurn, doneList);
                 
                                
                                //optSubSystem(doneList, tAtoms, tBonds, tAngles,
                                //              tTorsions, tRings, tPlas, tChs, false);
                                // std::cout << "The sub-system is optimized \n";
                                
                                iTurn++;    
                            }
                            else
                            {
                                //std::cout << "Grow the planar ring contain " << tAtoms[i_next].id
                                //          << std::endl;
                                aStartSet.clear();
                                if (i_atom == sAtom3)
                                {
                                    aStartSet.push_back(sAtom2);
                                }
                                else
                                {
                                    int iGP = tAtoms[i_atom].tree["parent"][0];
                                    aStartSet.push_back(iGP);
                                }
                                aStartSet.push_back(i_atom);
                                aStartSet.push_back(i_next);
                                //std::cout << "S1 set " << tAtoms[aStartSet[0]].id << "\t"
                                //          << tAtoms[aStartSet[1]].id << "\t" << tAtoms[aStartSet[2]].id << std::endl;
                                ringBuilder(aStartSet, tRings[*iRi], tAtoms, tBonds, 
                                            tAngles, tTorsions, tRings, tPlas, tChs, iTurn, doneList);
                 
                            }
                        
                        }
                        
                    
                    }
                    
                    dAtms.clear();
                    for (std::vector<int>::iterator iD=doneList.begin();
                            iD != doneList.end(); iD++)
                    {
                        dAtms.push_back(tAtoms[*iD]);
                    }
                    /*
                    int  a1    = (int)doneList.size();
                    if ( a1 > 5 )
                    {
                        std::string tName = "atoms_" + IntToStr(a1) + ".pdb";
                    
                        LIBMOL::outPDB(tName.c_str(), "XXX", dAtms);
                    }
                     */
                    /*
                    else if (a1 > 50)
                    {
                        exit(1);
                    }
                    
                     */
                }
                
                //std::cout << "0 i_conn " << i_conn << std::endl;
                //std::cout << "n_conn "   << n_conn << std::endl;
	                                     
                if (i_conn < n_conn )
                {
                    if(n_stack > NSTLIM)
                    {
                        std::cout << " Error in member function generateCoordCartToTors " 
                                  << std::endl << " Number in the stack is :  " 
                                  << n_stack  << std::endl;
                        exit(1);
                    }                                 

                    // Keep the junction point info in a stack
                    
                    i_conn               = i_conn  + 1;
	            iatom_stack[n_stack] = i_atom;
             
                    i_conn_stack[n_stack] = i_conn;

	            // cout << " The infor of the atom " <<  i_atom +1 
	            //     << " has been stored in the stack " << n_stack +1 << endl;
             
                    NBD_MCOPY(A, aMatStack[n_stack], dim);
                             
                    for ( i =0; i < dim; i++)
                    {
                        xyz_stack[n_stack][i] = tAtoms[i_atom].coords[i];
                    }
	      
                    n_stack   = n_stack + 1;

		    // std::cout << " n_stack :   " << n_stack << std::endl;  
                }
                
                CHANGE_AMAT(tAtoms,i_next, A);
                // Set the next atom to be the current atom in new loop

                i_atom     = i_next;
                i_conn     = 0;
                
                // n_conn = (int)tAtoms[i_atom].connAtoms.size()-1;
                
                //n_conn     = (int)tAtoms[i_atom].tree["children"].size();
                
                iFind = tAtoms[i_atom].tree.find("children");
                if (iFind !=tAtoms[i_atom].tree.end())
                {
                    n_conn = tAtoms[i_atom].tree["children"].size();
                }
                else
                {
                    n_conn = 0;
                }
               
	        //std::cout << " new i_atom: " << tAtoms[i_atom].id << std::endl;
	        //std::cout << " new n_conn: " << n_conn << std::endl;
	  
            }
            else if (n_conn == 1)
            { 
                int i_prev = tAtoms[i_atom].tree["parent"][0];
                
                int i_pprev;
               
                if (tAtoms[i_atom].id !=tAtoms[sAtom3].id)
                {
                    
                    i_pprev  = tAtoms[i_prev].tree["parent"][0];
                }
                else
                {
                    i_pprev  = sAtom1;
                    
                }
                
                //= tAtoms[i_prev].tree["parent"][0];
                //std::cout << "3 grand parent atom is " << tAtoms[i_pprev].id << std::endl 
                //          << "3 prev atom is " << tAtoms[i_prev].id << std::endl
                //          << "3 current atom is " << tAtoms[i_atom].id << std::endl;
                
                i_next = tAtoms[i_atom].tree["children"][i_conn];
                // std::cout << " The next atom is " << tAtoms[i_next].id << std::endl;
                
	        //  cout << " The next atom is " << i_next+1 << endl;
                if (std::find(doneList.begin(), doneList.end(), i_next) ==doneList.end())
                {
                    
                    //REAL tBo   = tAtoms[i_next].treeBond;
                    //REAL tTha  = tAtoms[i_next].treeAngle; 
                    //REAL tPhi  = tAtoms[i_next].treeTorsion;
                    std::vector<int> aSSet3;
                    aSSet3.push_back(i_pprev);
                    aSSet3.push_back(i_prev);
                    aSSet3.push_back(i_atom);
                    /**
                    PolToCart(tAtoms,i_atom, i_next, tBo, tTha, tPhi,  A);
                     */
                    growOneAtom(aSSet3, tAtoms, tBonds, tAngles, tTorsions, i_next);
                    doneList.push_back(i_next);
                    //std::cout <<  tAtoms[i_next].id << " done3 " << std::endl;
                    // grow other atoms linked to i_atom (like i_next)
                    /*
                    int iL=1;
                    for (std::vector<int>::iterator iCA=tAtoms[i_atom].connAtoms.begin();
                            iCA != tAtoms[i_atom].connAtoms.end(); iCA++)
                    {
                        if (*iCA !=i_next && *iCA !=i_prev 
                            && std::find(doneList.begin(), doneList.end(), *iCA) ==doneList.end())
                        {
                            tBo         = tAtoms[*iCA].treeBond;
                            tTha        = tAtoms[*iCA].treeAngle;
                            if (tAtoms[i_atom].bondingIdx==3)
                            {    
                                REAL tPhiI  = tPhi +iL*2.0*PI/3.0;   
                                PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                doneList.push_back(*iCA);
                                // std::cout <<  tAtoms[*iCA].id << " done " << std::endl;
                                iL++;
                            }
                            else if (tAtoms[i_atom].bondingIdx==2)
                            {
                                REAL tPhiI  = tPhi + PI;   
                                PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                doneList.push_back(*iCA);
                                // std::cout <<  tAtoms[*iCA].id << " done " << std::endl;
                            }  
                        }
                    }
                    */
                    
                    optSubSystem(doneList, tAtoms, tBonds, tAngles,
                                 tTorsions, tRings, tPlas, tChs, false);
                    //std::cout << "tree Node end: The sub-system is optimized \n";
                    //REAL tLeng  = tAtoms[i_next].treeBond;
                    //REAL tAng   = PI-tAtoms[i_next].treeAngle; 
                    //REAL tTor   = tAtoms[i_next].treeTorsion;   
                    //PolToCart(tAtoms, i_atom, i_next, tLeng, tAng, tTor, A);
                    //doneList.push_back(i_next);
                    //std::cout <<  tAtoms[i_next].id << " done " << std::endl;
                }
	                                    //  strcpy(tAllAtoms[i_next].type, "R");
                
                
                dAtms.clear();
                for (std::vector<int>::iterator iD=doneList.begin();
                            iD != doneList.end(); iD++)
                {
                    dAtms.push_back(tAtoms[*iD]);
                }
                  
                /*
                int  a1    = (int)doneList.size();
               
                if ( a1 > 5 )
                {
                    std::string tName = "atoms_" + IntToStr(a1) + ".pdb";
                    
                    LIBMOL::outPDB(tName.c_str(), "XXX", dAtms);
                }
                 */
                    /*
                    else if (a1 > 50)
                    {
                        exit(1);
                    }
                    
                     */
                
                
            
                
                CHANGE_AMAT(tAtoms, i_next, A);

                
                // Set the next atom to be the current atom in new loop

                i_atom     = i_next;
                i_conn     = 0;
                
                
                iFind = tAtoms[i_atom].tree.find("children");
                if (iFind !=tAtoms[i_atom].tree.end())
                {
                    n_conn = tAtoms[i_atom].tree["children"].size();
                }
                else
                {
                    n_conn = 0;
                }
                
	        //std::cout << " new i_atom: " << tAtoms[i_atom].id << std::endl;
	        //std::cout << " new n_conn: " << n_conn   << std::endl;
            }
            else if (n_conn == 0 && n_stack > 0)
            {
                do
                {
                    n_stack = n_stack-1;
	            // std::cout << "The stack level is " << n_stack << std::endl;
                    if(n_stack < 0)
                    {
                        break;
                    }
                    i_atom = iatom_stack[n_stack];

	            // std::cout << " restart from i_atom " << tAtoms[i_atom].id << std::endl;
         
                    i_conn = i_conn_stack[n_stack];

	            // std::cout << " restart1 from the " << i_conn+1 
	            //           << " node that linked to the above atom " << std::endl;

                    NBD_MCOPY(aMatStack[n_stack], A, dim);

	            //for( i =0; i < dim; i++)
	            //{
                    //  xyz[i] = xyz_stack[n_stack][i];
	            //}

                    // n_conn = (int)tAtoms[i_atom].tree["children"].size();
                    iFind = tAtoms[i_atom].tree.find("children");
                    if (iFind !=tAtoms[i_atom].tree.end())
                    {
                        n_conn = tAtoms[i_atom].tree["children"].size();
                    }
                    else
                    {
                        n_conn = 0;
                    }
	            //    cout << "restart n_conn " << n_conn << endl;
	            //    cout << " continue      " << endl;             
	            //    cin.get();
                
                }while(i_conn >= n_conn && n_stack >= 0);
                /*
                std::cout << "out while  " << std::endl
                          << "n_conn " << n_conn << std::endl
                          << "i_conn " << i_conn << std::endl
                          << "n_stack " <<n_stack << std::endl
                          << "i_atom "  << tAtoms[i_atom].id << std::endl;
                 */
            }
            
            // std::cout << "here " << std::endl;
        }
        //std::cout << "New coordinates after angle-to-coordinates transfer are " << std::endl;

        //for (int i =0; i < (int)tAtoms.size(); i++)
        //{
        //    std::cout << "The " <<i+1<<"th atom is at" << std::endl;
        //    for (int j=0; j < dim; j++)
        // 	{
        //        std::cout << tAtoms[i].coords[j] << "\t";
        //      }
        //    std::cout << std::endl;
        //  }
  
        // cout << " Continue ? " << endl;
        // cin.get();

        // Release memory and null point variables

        delete [] v1;
        v1 = NULL;

        delete v2;
        v2 = NULL;

        delete [] v3;
        v3 = NULL;
  
        delete [] tmp_v;
        tmp_v = 0;
  
        delete [] X_new;
        X_new = NULL;

        delete [] Y_new;
        Y_new = NULL;

        delete [] Z_new;
        Z_new = NULL;
  
        for ( i =0; i < dim; i++)
        {
            delete [] rotate_tem[i];
            rotate_tem[i] =NULL;
        }

        delete [] rotate_tem;
        rotate_tem = NULL;
  
        for (i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete [] A;
        A = 0;

        delete [] iatom_stack;
        iatom_stack  = 0;

        delete [] i_conn_stack;
        i_conn_stack = 0;

        for ( i =0; i < NSTLIM; i++)
        {
            for(j =0; j < NSTLIM2; j++)
            {
                delete []    aMatStack[i][j];
                aMatStack[i][j] = 0;
            }
            delete [] aMatStack[i];
            aMatStack[i] =0;
        }
        delete [] aMatStack;
        aMatStack = 0; 
       
    }
    
    
    void TransCoords::PolToCart(std::vector<AtomDict> & tAtoms,
                                int i_cur, int i_next, REAL** aMat)
    {
        int      dim = (int)tAtoms[0].coords.size();
        int      i;
        REAL     th_i, phi_i, length;
        REAL     c_th_i, s_th_i, c_ph_i, s_ph_i;
        REAL  *  dX = new REAL [dim];
        REAL  *  dY = new REAL [dim]; 
        REAL     ddX;
        
        length = tAtoms[i_next].treeBond;
        th_i   = PI-tAtoms[i_next].treeAngle; 
        phi_i  = tAtoms[i_next].treeTorsion;
        
        //std::cout << "Current atom is "   << tAtoms[i_next].id << std::endl
        //          << " its tree bond is " << length << std::endl
        //          << " its angle is "     << th_i   << "(" << th_i*PID180 << ")" << std::endl
        //          << " its torsion is "   << phi_i  << std::endl;
       
                   
        c_th_i = cos(th_i);
        s_th_i = sin(th_i);
        c_ph_i = cos(phi_i);
        s_ph_i = sin(phi_i);
        
        if(fabs(c_th_i)<1.0e-16)
        {
            c_th_i = 0.0;
        }
        else if((1-fabs(c_th_i))<1.0e-16)
        {
            c_th_i = Sign(1.0,c_th_i);
        }
        
        if(fabs(s_th_i)<1.0e-16)
        {
            s_th_i = 0.0;
        }
        else if((1-fabs(s_th_i))<1.0e-16)
        {
            s_th_i = Sign(1.0,s_th_i);
        }

        if(fabs(c_ph_i)<1.0e-16)
        {
            c_ph_i = 0.0;
        }
        else if ((1-fabs(c_ph_i))<1.0e-10)
        {
            c_ph_i = Sign(1.0,c_ph_i);
        }
        
        if(fabs(s_ph_i)<1.0e-16)
        {
            s_ph_i = 0.0;
        }
        else if((1-fabs(s_ph_i))<1.0e-16)
        {
            s_ph_i = Sign(1.0,s_ph_i);
        }
 
        //  std::cout << " The current atom is " << i_cur +1 << std::endl;
        //  std::cout << " The next atom is    " << i_next +1 << std::endl;
  
        //  std::cout << " length of the next atom " << length << std::endl;   
        //  std::cout << " th_i of the next atom " << th_i   << std::endl;
        //  std::cout << " phi_i of the next atom " << phi_i << std::endl;
  
        // Three dimensions only at the present

        dX[0] =  c_th_i*length;
        // std::cout << " ph_i is " << ph_i << " cos(ph_i) " << c_ph_i << std::endl;
        // std::cout << " dX[0] = " << dX[0] << std::endl;

        ddX   =   s_th_i*length;
        dX[1] =   c_ph_i*ddX;  
        dX[2] =   s_ph_i*ddX;
        
        

        //   for (i = 0; i < dim; i++)
        //     {
        //      for (j =0; j < dim; j++)
        //        {
        //          cout << "aMat["<<i <<"]["<<j<<"]= " << aMat[i][j] << endl;
        //        }
        //    }

  
        int n_dim_t = dim;

        MatMultip(n_dim_t, aMat,dX,dY);


        for (i =0; i < dim; i++)
        {
            if(fabs(dY[i]) < 1.0e-7)
            {
                dY[i] = 0.0;
            }
            //    cout << " dX["<< i <<"]= " << dX[i] << endl; 
            //    cout << " dY["<< i <<"]= " << dY[i] << endl;
        }
 
 
        //  cout << " Before transfering " << endl;
        // for (i =0; i < dim; i++)
        //    {
        //      cout << "atoms[" << i_next+1 << "].coord["<<i<<"] = " 
        //           << atoms[i_next].coords[i] << endl; 
        //    }
        
        
        // cout <<" after transfer " << endl;
        
        for ( i =0; i < dim; i++)
        {
            tAtoms[i_next].coords[i]  = tAtoms[i_cur].coords[i]+ dY[i];
      
        }

        //  cout << "continue " << endl;
        //  cin.get();
        
        delete [] dX;
        dX = 0;
        delete [] dY;
        dY = 0;
  
    }
    
    void TransCoords::PolToCart(std::vector<AtomDict> & tAtoms,
                                int i_cur, int i_next, 
                                REAL tBond, REAL tAng, REAL tTor,
                                REAL ** aMat)
    {
        int      dim = (int)tAtoms[0].coords.size();
        int      i;
        REAL     th_i, phi_i, length;
        REAL     c_th_i, s_th_i, c_ph_i, s_ph_i;
        REAL  *  dX = new REAL [dim];
        REAL  *  dY = new REAL [dim]; 
        REAL     ddX;
        
        length = tBond;
        th_i   = PI-tAng; 
        phi_i  = tTor;
        
        //std::cout << "in: Current atom is atom number " << tAtoms[i_next].seriNum    
        //          << " id " << tAtoms[i_next].id << std::endl
        //          << " its tree bond is " << length << std::endl
        //          << " its angle is "     << th_i << "(" << th_i*PID180 << ")"  << std::endl
        //          << " its torsion is "   << phi_i  << std::endl;
       
                   
        c_th_i = cos(th_i);
        s_th_i = sin(th_i);
        c_ph_i = cos(phi_i);
        s_ph_i = sin(phi_i);
        
        if(fabs(c_th_i)<1.0e-16)
        {
            c_th_i = 0.0;
        }
        else if((1-fabs(c_th_i))<1.0e-16)
        {
            c_th_i = Sign(1.0,c_th_i);
        }
        
        if(fabs(s_th_i)<1.0e-16)
        {
            s_th_i = 0.0;
        }
        else if((1-fabs(s_th_i))<1.0e-16)
        {
            s_th_i = Sign(1.0,s_th_i);
        }

        if(fabs(c_ph_i)<1.0e-16)
        {
            c_ph_i = 0.0;
        }
        else if ((1-fabs(c_ph_i))<1.0e-10)
        {
            c_ph_i = Sign(1.0,c_ph_i);
        }
        
        if(fabs(s_ph_i)<1.0e-16)
        {
            s_ph_i = 0.0;
        }
        else if((1-fabs(s_ph_i))<1.0e-16)
        {
            s_ph_i = Sign(1.0,s_ph_i);
        }
 
        //  std::cout << " The current atom is " << i_cur +1 << std::endl;
        //  std::cout << " The next atom is    " << i_next +1 << std::endl;
  
        //  std::cout << " length of the next atom " << length << std::endl;   
        //  std::cout << " th_i of the next atom " << th_i   << std::endl;
        //  std::cout << " phi_i of the next atom " << phi_i << std::endl;
  
        // Three dimensions only at the present

        dX[0] =  c_th_i*length;
        // std::cout << " ph_i is " << ph_i << " cos(ph_i) " << c_ph_i << std::endl;
        // std::cout << " dX[0] = " << dX[0] << std::endl;

        ddX   =   s_th_i*length;
        dX[1] =   c_ph_i*ddX;  
        dX[2] =   s_ph_i*ddX;
        
        

        //   for (i = 0; i < dim; i++)
        //     {
        //      for (j =0; j < dim; j++)
        //        {
        //          cout << "aMat["<<i <<"]["<<j<<"]= " << aMat[i][j] << endl;
        //        }
        //    }

  
        int n_dim_t = dim;

        MatMultip(n_dim_t, aMat,dX,dY);


        for (i =0; i < dim; i++)
        {
            if(fabs(dY[i]) < 1.0e-7)
            {
                dY[i] = 0.0;
            }
            //std::cout << " dX["<< i <<"]= " << dX[i] << std::endl; 
            //std::cout << " dY["<< i <<"]= " << dY[i] << std::endl;
        }
 
 
         //std::cout << " Before transfering " << std::endl;
         //for (i =0; i < dim; i++)
         //   {
         //     std::cout << "tAtoms[" << i_next << "].coord["<<i<<"] = " 
         //               << tAtoms[i_next].coords[i] << std::endl; 
         //   }
        
        
        // std::cout <<" after transfer " << std::endl;
        
        for ( i =0; i < dim; i++)
        {
            tAtoms[i_next].coords[i]  = tAtoms[i_cur].coords[i]+ dY[i];
            //std::cout << "tAtoms[" << i_next << "].coord["<<i<<"] = " 
            //          << tAtoms[i_next].coords[i] << std::endl; 
        }
        
        //  cout << "continue " << endl;
        //  cin.get();
        
        delete [] dX;
        dX = 0;
        delete [] dY;
        dY = 0;
  
    }
    
    
    void TransCoords::PolToCart(AtomDict& tCurAtom, 
                                std::vector<AtomDict>::iterator tNextAtom,
                                REAL tBond, REAL tAng, REAL tTor, 
                                REAL** aMat)
    {
        int      dim = (int)tCurAtom.coords.size();
        int      i;
        REAL     th_i, phi_i, length;
        REAL     c_th_i, s_th_i, c_ph_i, s_ph_i;
        REAL  *  dX = new REAL [dim];
        REAL  *  dY = new REAL [dim]; 
        REAL     ddX;
        
        length = tBond;
        th_i   = PI-tAng; 
        phi_i  = tTor;
        
        //std::cout << "in: Current atom is "   << tAtoms[i_next].id << std::endl
        //          << " its tree bond is " << length << std::endl
        //          << " its angle is "     << th_i << "(" << th_i*PID180 << ")"  << std::endl
        //          << " its torsion is "   << phi_i  << std::endl;
       
                   
        c_th_i = cos(th_i);
        s_th_i = sin(th_i);
        c_ph_i = cos(phi_i);
        s_ph_i = sin(phi_i);
        
        if(fabs(c_th_i)<1.0e-16)
        {
            c_th_i = 0.0;
        }
        else if((1-fabs(c_th_i))<1.0e-16)
        {
            c_th_i = Sign(1.0,c_th_i);
        }
        
        if(fabs(s_th_i)<1.0e-16)
        {
            s_th_i = 0.0;
        }
        else if((1-fabs(s_th_i))<1.0e-16)
        {
            s_th_i = Sign(1.0,s_th_i);
        }

        if(fabs(c_ph_i)<1.0e-16)
        {
            c_ph_i = 0.0;
        }
        else if ((1-fabs(c_ph_i))<1.0e-10)
        {
            c_ph_i = Sign(1.0,c_ph_i);
        }
        
        if(fabs(s_ph_i)<1.0e-16)
        {
            s_ph_i = 0.0;
        }
        else if((1-fabs(s_ph_i))<1.0e-16)
        {
            s_ph_i = Sign(1.0,s_ph_i);
        }
 
        //  std::cout << " The current atom is " << i_cur +1 << std::endl;
        //  std::cout << " The next atom is    " << i_next +1 << std::endl;
  
        //  std::cout << " length of the next atom " << length << std::endl;   
        //  std::cout << " th_i of the next atom " << th_i   << std::endl;
        //  std::cout << " phi_i of the next atom " << phi_i << std::endl;
  
        // Three dimensions only at the present

        dX[0] =  c_th_i*length;
        // std::cout << " ph_i is " << ph_i << " cos(ph_i) " << c_ph_i << std::endl;
        // std::cout << " dX[0] = " << dX[0] << std::endl;

        ddX   =   s_th_i*length;
        dX[1] =   c_ph_i*ddX;  
        dX[2] =   s_ph_i*ddX;
        
        

        //   for (i = 0; i < dim; i++)
        //     {
        //      for (j =0; j < dim; j++)
        //        {
        //          cout << "aMat["<<i <<"]["<<j<<"]= " << aMat[i][j] << endl;
        //        }
        //    }

  
        int n_dim_t = dim;

        MatMultip(n_dim_t, aMat,dX,dY);


        for (i =0; i < dim; i++)
        {
            if(fabs(dY[i]) < 1.0e-7)
            {
                dY[i] = 0.0;
            }
            //    cout << " dX["<< i <<"]= " << dX[i] << endl; 
            //    cout << " dY["<< i <<"]= " << dY[i] << endl;
        }
 
 
        //  cout << " Before transfering " << endl;
        // for (i =0; i < dim; i++)
        //    {
        //      cout << "atoms[" << i_next+1 << "].coord["<<i<<"] = " 
        //           << atoms[i_next].coords[i] << endl; 
        //    }
        
        
        // cout <<" after transfer " << endl;
        
        for ( i =0; i < dim; i++)
        {
            tNextAtom->coords[i]  = tCurAtom.coords[i]+ dY[i];
      
        }

        //  cout << "continue " << endl;
        //  cin.get();
        
        delete [] dX;
        dX = 0;
        delete [] dY;
        dY = 0;
  
    }
    
        void TransCoords::PolToCart(AtomDict& tCurAtom, 
                                AtomDict& tNextAtom, 
                                REAL tBond, REAL tAng, REAL tTor, 
                                REAL** aMat)
    {
        int      dim = (int)tCurAtom.coords.size();
        int      i;
        REAL     th_i, phi_i, length;
        REAL     c_th_i, s_th_i, c_ph_i, s_ph_i;
        REAL  *  dX = new REAL [dim];
        REAL  *  dY = new REAL [dim]; 
        REAL     ddX;
        
        length = tBond;
        th_i   = PI-tAng; 
        phi_i  = tTor;
        
        //std::cout << "in: Current atom is "   << tAtoms[i_next].id << std::endl
        //          << " its tree bond is " << length << std::endl
        //          << " its angle is "     << th_i << "(" << th_i*PID180 << ")"  << std::endl
        //          << " its torsion is "   << phi_i  << std::endl;
       
                   
        c_th_i = cos(th_i);
        s_th_i = sin(th_i);
        c_ph_i = cos(phi_i);
        s_ph_i = sin(phi_i);
        
        if(fabs(c_th_i)<1.0e-16)
        {
            c_th_i = 0.0;
        }
        else if((1-fabs(c_th_i))<1.0e-16)
        {
            c_th_i = Sign(1.0,c_th_i);
        }
        
        if(fabs(s_th_i)<1.0e-16)
        {
            s_th_i = 0.0;
        }
        else if((1-fabs(s_th_i))<1.0e-16)
        {
            s_th_i = Sign(1.0,s_th_i);
        }

        if(fabs(c_ph_i)<1.0e-16)
        {
            c_ph_i = 0.0;
        }
        else if ((1-fabs(c_ph_i))<1.0e-10)
        {
            c_ph_i = Sign(1.0,c_ph_i);
        }
        
        if(fabs(s_ph_i)<1.0e-16)
        {
            s_ph_i = 0.0;
        }
        else if((1-fabs(s_ph_i))<1.0e-16)
        {
            s_ph_i = Sign(1.0,s_ph_i);
        }
 
        //  std::cout << " The current atom is " << i_cur +1 << std::endl;
        //  std::cout << " The next atom is    " << i_next +1 << std::endl;
  
        //  std::cout << " length of the next atom " << length << std::endl;   
        //  std::cout << " th_i of the next atom " << th_i   << std::endl;
        //  std::cout << " phi_i of the next atom " << phi_i << std::endl;
  
        // Three dimensions only at the present

        dX[0] =  c_th_i*length;
        // std::cout << " ph_i is " << ph_i << " cos(ph_i) " << c_ph_i << std::endl;
        // std::cout << " dX[0] = " << dX[0] << std::endl;

        ddX   =   s_th_i*length;
        dX[1] =   c_ph_i*ddX;  
        dX[2] =   s_ph_i*ddX;
        
        

        //   for (i = 0; i < dim; i++)
        //     {
        //      for (j =0; j < dim; j++)
        //        {
        //          cout << "aMat["<<i <<"]["<<j<<"]= " << aMat[i][j] << endl;
        //        }
        //    }

  
        int n_dim_t = dim;

        MatMultip(n_dim_t, aMat,dX,dY);


        for (i =0; i < dim; i++)
        {
            if(fabs(dY[i]) < 1.0e-7)
            {
                dY[i] = 0.0;
            }
            //    cout << " dX["<< i <<"]= " << dX[i] << endl; 
            //    cout << " dY["<< i <<"]= " << dY[i] << endl;
        }
 
 
        //  cout << " Before transfering " << endl;
        // for (i =0; i < dim; i++)
        //    {
        //      cout << "atoms[" << i_next+1 << "].coord["<<i<<"] = " 
        //           << atoms[i_next].coords[i] << endl; 
        //    }
        
        
        // cout <<" after transfer " << endl;
        
        for ( i =0; i < dim; i++)
        {
            tNextAtom.coords[i]  = tCurAtom.coords[i]+ dY[i];
      
        }

        //  cout << "continue " << endl;
        //  cin.get();
        
        delete [] dX;
        dX = 0;
        delete [] dY;
        dY = 0;
  
    }
    
    void TransCoords::setCurStartSet(std::vector<int>& tSet, 
                                     std::vector<int>& curSet,
                                     int  tCur,
                                     std::string tDirect, 
                                     RingDict &  tRing,
                                     std::vector<int> tDoneSet)
    {
        int  p=-1, gp=-1;
        curSet.clear();
        std::string aD;
        if (tDirect=="c")
        {
            aD = "ac";
        }
        else
        {
            aD = "c";
        } 
            
        if (std::find(tDoneSet.begin(), tDoneSet.end(), tRing.atomsLink[tCur][aD])==tDoneSet.end())
        {
            p=tRing.atomsLink[tCur][aD];
            if (std::find(tDoneSet.begin(), tDoneSet.end(), tRing.atomsLink[p][aD])==tDoneSet.end())
            {
                gp=tRing.atomsLink[p][aD];
            }
        }
        
        if (p !=-1 and gp !=-1)
        {
            // using ring atoms as the start set
            curSet.push_back(gp);
            curSet.push_back(p);
            curSet.push_back(tCur);
        }
        else
        {
            for(std::vector<int>::iterator iA=tSet.begin();
                    iA !=tSet.end(); iA++)
            {
                curSet.push_back(*iA);
            }
        }
        
    }
    
    void TransCoords::growOneAtom(std::vector<int>         & sSet, 
                                  std::vector<AtomDict>    & tAtoms, 
                                  std::vector<BondDict>    & tBonds, 
                                  std::vector<AngleDict>   & tAngles, 
                                  std::vector<TorsionDict> & tTorsions, 
                                  int  tIdx)
    {
        
        int dim = (int)tAtoms[0].coords.size();

        REAL **A = new REAL * [dim];
        for (int i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (int j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (int i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sSet[1]].coords[i] - tAtoms[sSet[0]].coords[i];
            v2[i] = tAtoms[sSet[2]].coords[i] - tAtoms[sSet[1]].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }
        

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( int i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  

        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start grow the ring from the starting atom
        
        
       
        int  idxBond = getBond(tBonds, sSet[2], tIdx);
        int  idxAng  = getAngle(tAngles, sSet[2], sSet[1], tIdx);
        int  idxTor  = getTorsion(tTorsions, sSet[0], sSet[1], sSet[2], tIdx);
        
        if (idxBond !=-1 && idxAng !=-1 && idxTor !=-1)
        {
            REAL aBond=tBonds[idxBond].value;
            REAL aAng =tAngles[idxAng].valueST;
            REAL aTor =tTorsions[idxTor].value*PI180;
            PolToCart(tAtoms, sSet[2], tIdx, aBond, aAng, aTor, A);
            
        }
        else
        {
            std::cout << "could not find either bond, or angle or torsion for atoms set "
                    << tAtoms[sSet[0]].id << " " << tAtoms[sSet[1]].id 
                    << tAtoms[sSet[2]].id << " " << tAtoms[tIdx].id << std::endl;
            exit(1);
        }
        
                // Release memory and null point variables

        delete [] v1;
        v1 = 0;

        delete v2;
        v2 = 0;
  
        delete [] tmp_v;
        tmp_v = 0;
  
        delete [] X_new;
        X_new = 0;

        delete [] Y_new;
        Y_new = 0;

        delete [] Z_new;
        Z_new = 0;
        
        for (int i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete [] A;
        A = 0;
  
    }
    
    
    void TransCoords::growOneAtom(std::vector<int>           & sSet,
                                  std::vector<AtomDict>   & tAtoms,
                                  REAL                      tBondV, 
                                  REAL                      tAngV, 
                                  REAL                      tTorV,
                                  int                       tIdx)
    {
        
        int dim = (int)tAtoms[0].coords.size();

        REAL **A = new REAL * [dim];
        for (int i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (int j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (int i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sSet[1]].coords[i] - tAtoms[sSet[0]].coords[i];
            v2[i] = tAtoms[sSet[2]].coords[i] - tAtoms[sSet[1]].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }
        

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( int i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  

        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start grow the ring from the starting atom
        
        PolToCart(tAtoms, sSet[2], tIdx, tBondV, tAngV, tTorV, A);
        
        // Release memory and null point variables

        delete [] v1;
        v1 = NULL;

        delete v2;
        v2 = NULL;
  
        delete [] tmp_v;
        tmp_v = NULL;
  
        delete [] X_new;
        X_new = NULL;

        delete [] Y_new;
        Y_new = NULL;

        delete [] Z_new;
        Z_new = NULL;
        
        for (int i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = NULL;
        }
        delete [] A;
        A = NULL;
  
    }
    
    void TransCoords::growOneAtom(AtomDict& tAtom1, AtomDict& tAtom2, 
                                  AtomDict& tAtom3, 
                                  std::vector<AtomDict>::iterator tNextAtom, 
                                  REAL tBondV, REAL tAngV, REAL tTorV)
    {
                
        int dim = (int)tAtom1.coords.size();

        REAL **A = new REAL * [dim];
        for (int i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (int j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (int i =0; i < dim; i++)
        {
            v1[i] = tAtom2.coords[i] - tAtom1.coords[i];
            v2[i] = tAtom3.coords[i] - tAtom2.coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }
        

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( int i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  

        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start grow the ring from the starting atom
        
        PolToCart(tAtom3, tNextAtom, tBondV, tAngV, tTorV, A);
        
        // Release memory and null point variables

        delete [] v1;
        v1 = NULL;

        delete v2;
        v2 = NULL;
  
        delete [] tmp_v;
        tmp_v = NULL;
  
        delete [] X_new;
        X_new = NULL;

        delete [] Y_new;
        Y_new = NULL;

        delete [] Z_new;
        Z_new = NULL;
        
        for (int i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = NULL;
        }
        delete [] A;
        A = NULL;
    }
    
    void TransCoords::growOneRingNode(std::vector<int>          & sSet,    
                                      std::vector<AtomDict>     & tAtoms, 
                                      std::vector<BondDict>     & tBonds, 
                                      std::vector<AngleDict>    & tAngles, 
                                      std::vector<TorsionDict>  & tTorsions,
                                      std::vector<ChiralDict>   & tChs,
                                      bool                        inRing,
                                      int                         iSeq,                                          
                                      int tIdx, std::vector<int>& tDoneSet)
    {
        //std::cout << "growOneRingNode " << std::endl;
        int dim = (int)tAtoms[0].coords.size();

        REAL **A = new REAL * [dim];
        for (int i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (int j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (int i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sSet[1]].coords[i] - tAtoms[sSet[0]].coords[i];
            v2[i] = tAtoms[sSet[2]].coords[i] - tAtoms[sSet[1]].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }
        

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( int i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  

        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start grow an atom in the ring from the starting atom
        
        
       
        int  idxBond = getBond(tBonds, sSet[2], tIdx);
        int  idxAng  = getAngle(tAngles, sSet[2], sSet[1], tIdx);
        int  idxTor  = getTorsion(tTorsions, sSet[0], sSet[1], sSet[2], tIdx);
        
        if (idxBond !=-1 && idxAng !=-1 && idxTor !=-1)
        {
            REAL aBond=tBonds[idxBond].value;
            REAL aAng =tAngles[idxAng].valueST;
            REAL aTor;
            if (inRing)
            {
                aTor=tTorsions[idxTor].value*PI180;
                /*
                if(tAtoms[sSet[2]].bondingIdx==3)
                {
                    if (iSeq >0)
                    {
                        aTor = PI/3.0;
                    }
                    else
                    {
                        aTor = -PI/3.0;
                    }
                }
                else
                {
                    aTor = 0.0;
                 }
                 */
            }
            else
            {
                
                if(tAtoms[sSet[2]].bondingIdx==3)
                {
                    if (iSeq >0)
                    {
                        aTor = 2.0*PI/3.0;
                    }
                    else
                    {
                        aTor = -2.0*PI/3.0;
                    }
                }
                else
                {
                    // aTor = PI;
                    aTor=tTorsions[idxTor].value*PI180;
                }
                 
                // aTor = PI;
            }
            if (std::find(doneList.begin(), doneList.end(), tIdx) ==doneList.end())
            {
                //std::cout << "Bond  " << aBond << std::endl;
                //std::cout << "Angle " << aAng  << std::endl;
                //std::cout << "Tor "   << aTor << std::endl;
                //std::cout << "start Atom " << tAtoms[sSet[2]].id << std::endl;
                //std::cout << "next Atom " << tAtoms[tIdx].id << std::endl;
                PolToCart(tAtoms, sSet[2], tIdx, aBond, aAng, aTor, A);
                tDoneSet.push_back(tIdx);
                
                for (unsigned i=0; i << tAtoms[tIdx].coords.size(); i++)
                {
                    std::cout << tAtoms[tIdx].coords[i] << std::endl;
                }
               
            }
            
            
            // build the atoms that are not on the ring but bonding with 
            // the atom sSet[2] (tIdx is also bonding with sSet[2] 
            
            // grow other atoms linked to i_atom (like i_next)
            
                    if (tAtoms[sSet[2]].chiralIdx !=0)
                    {
                        int aCh = tAtoms[sSet[2]].inChirals[0];
                        bool lS= false;
                        std::vector<int> t1, rTable;
                        for (std::vector<int>::iterator irA=tChs[aCh].mutTable[sSet[1]].begin();
                               irA!=tChs[aCh].mutTable[sSet[1]].end(); irA++)
                        {
                            if (lS)
                            {
                                rTable.push_back(*irA);
                            }
                            else if (*irA==tIdx)
                            {
                                lS=true;
                            }
                            else
                            {
                                t1.push_back(*irA);
                            }
                        }
                        for (std::vector<int>::iterator irA2=t1.begin();
                                irA2 !=t1.end(); irA2++)
                        {
                            rTable.push_back(*irA2);
                        }
                        
                        //std::cout << "Chiral center " << tAtoms[sSet[2]].id << std::endl;
                        //std::cout << "chiral description  " << tChs[aCh].signST << std::endl;
                        //std::cout << "Atom arrangement " << std::endl;
                        //std::cout << "from atom " << tAtoms[sSet[1]].id  << " clockwise " << std::endl;
                        //std::cout << "atom " << tAtoms[tIdx].id << std::endl;
                        //for (std::vector<int>::iterator iNA=rTable.begin();
                        //        iNA !=rTable.end(); iNA++)
                        //{
                        //    std::cout << "Atom " << tAtoms[*iNA].id << std::endl;
                        //}
                        
                        int iL=1;
                        
                        for (std::vector<int>::iterator iNodeA=rTable.begin();
                                iNodeA !=rTable.end(); iNodeA++)
                        {
                            //std::cout << "iNodeA " << tAtoms[*iNodeA].id << std::endl;
                            
                            int bIdx    = getBond(tBonds, sSet[2], *iNodeA);
                            REAL tBo         = tBonds[bIdx].value;
                            int aIdx    = getAngle(tAngles, sSet[2], sSet[1], *iNodeA);
                            REAL tTha        = tAngles[aIdx].valueST;
                            //std::cout << "iL= " << iL << std::endl
                            //          << " tTor" << aTor << std::endl
                            //          << " deltal " << 2.0*PI/3.0 << std::endl;
                            REAL tPhiI  = aTor +iL*2.0*PI/3.0;
                            if (std::find(doneList.begin(), doneList.end(), *iNodeA) ==doneList.end())
                            {
                                // PolToCart(tAtoms,i_atom, *iNodeA, tBo, tTha, tPhiI,  A);
                                growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iNodeA);
                                doneList.push_back(*iNodeA);
                                //std::cout <<  tAtoms[*iNodeA].id << " done1 " << std::endl;
                            }
                            iL++;
                        }
                        // optSubSystem(doneList, tAtoms, tBonds, tAngles,
                        //             tTorsions, tRings, tPlas, tChs);
                        //std::cout << "The sub-system is optimized \n";
                    }
                    else
                    {
                        int iL=1;
                        for (std::vector<int>::iterator iCA=tAtoms[sSet[2]].connAtoms.begin();
                                iCA != tAtoms[sSet[2]].connAtoms.end(); iCA++)
                        {
                            if (*iCA !=tIdx && *iCA !=sSet[1] 
                                && std::find(doneList.begin(), doneList.end(), *iCA) ==doneList.end())
                            {
                                REAL tBo         = tAtoms[*iCA].treeBond;
                                REAL tTha        = tAtoms[*iCA].treeAngle;
                                if (tAtoms[sSet[2]].bondingIdx==3)
                                {    
                                    REAL tPhiI  = aTor +iL*2.0*PI/3.0;
                                    growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    //std::cout <<  tAtoms[*iCA].id << " done21 " << std::endl;
                                    iL++;
                                }
                                else if (tAtoms[sSet[2]].bondingIdx==2)
                                {
                                    REAL tPhiI  = aTor + PI;   
                                    growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    // std::cout <<  tAtoms[*iCA].id << " done22 " << std::endl;
                                }
                            }
                        }
                    }
            
            
            /*
            int iL=1;
            for (int iC=0; iC < (int)tAtoms[sSet[2]].connAtoms.size(); iC++)
            {
                int iAC=tAtoms[sSet[2]].connAtoms[iC];
                if( iAC !=sSet[1] && iAC !=tIdx
                    
                    && std::find(tDoneSet.begin(), tDoneSet.end(), iAC)==tDoneSet.end())
                {
                    idxBond = getBond(tBonds, sSet[2], iAC);
                    idxAng  = getAngle(tAngles, sSet[2], sSet[1], iAC);
                    if (idxBond !=-1 && idxAng != -1)
                    {
                        aBond=tBonds[idxBond].value;
                        aAng =tAngles[idxAng].valueST;
                        REAL aTorI;
                        if (tAtoms[sSet[2]].bondingIdx==2)
                        {
                            aTorI = aTor + PI;
                        }
                        else if (tAtoms[sSet[2]].bondingIdx == 3)
                        {
                            aTorI = aTor +iL*2.0*PI/3.0;
                            iL++;
                        }
                        else
                        {
                            std::cout << "wrong connection in ring-building. Bug!"
                                      << std::endl;
                            exit(1);
                        }
                        
                        for (int i =0; i < dim; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                A[i][j] =0.0;
                            }
                        }
                        
                        InitAMat2(X_new,Y_new,Z_new,A, dim);
                        if (std::find(doneList.begin(), doneList.end(), iAC) ==doneList.end())
                        {
                            PolToCart(tAtoms, sSet[2], iAC, aBond, aAng, aTorI, A);
                            tDoneSet.push_back(iAC);
                        }
                    }
                    else
                    {
                        std::cout << "could not find either bond, or angle for atoms set "
                                  << tAtoms[sSet[1]].id << " " 
                                  << tAtoms[sSet[2]].id << " " 
                                  << tAtoms[iAC].id << std::endl;
                                  exit(1);
                    }
                }
            }
            */
        }
        else
        {
            std::cout << "could not find either bond, or angle or torsion for atoms set "
                    << tAtoms[sSet[0]].id << " " << tAtoms[sSet[1]].id 
                    << tAtoms[sSet[2]].id << " " << tAtoms[tIdx].id << std::endl;
            exit(1);
        }
        
        
        
        // Release memory and null point variables

        delete [] v1;
        v1 = 0;

        delete v2;
        v2 = 0;
  
        delete [] tmp_v;
        tmp_v = 0;
  
        delete [] X_new;
        X_new = 0;

        delete [] Y_new;
        Y_new = 0;

        delete [] Z_new;
        Z_new = 0;
        
        for (int i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete [] A;
        A = 0;
  
    }
    
    void TransCoords::growOneRingNode2(std::vector<int>          & sSet,    
                                       std::vector<AtomDict>     & tAtoms, 
                                       std::vector<BondDict>     & tBonds, 
                                       std::vector<AngleDict>    & tAngles, 
                                       std::vector<TorsionDict>  & tTorsions,
                                       std::vector<ChiralDict>   & tChs,
                                       bool                        inRing,
                                       int                         iSeq,                                          
                                       int tIdx, std::vector<int>& tDoneSet)
    {
        //std::cout << "growOneRingNode2 " << std::endl;
        int dim = (int)tAtoms[0].coords.size();

        REAL **A = new REAL * [dim];
        for (int i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
            for (int j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        // Pick up the starting atom according to the selected torsion angle
        
        REAL *v1    = new REAL [dim];
        REAL *v2    = new REAL [dim];
        
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;
  
        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        
        for (int i =0; i < dim; i++)
        {
            v1[i] = tAtoms[sSet[1]].coords[i] - tAtoms[sSet[0]].coords[i];
            v2[i] = tAtoms[sSet[2]].coords[i] - tAtoms[sSet[1]].coords[i];
        }
        
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                Z_new [i] = 0.0;
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }
        

        // Set New Y axis

        CrossP(Z_new,v2, tmp_v);
        mod_tmp_v = DotP(tmp_v,tmp_v);
        mod_tmp_v = sqrt(mod_tmp_v);
  
        if(mod_tmp_v)
        {
            for ( int i = 0; i < dim; i++)
            {
                Y_new [i] = 0.0;
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " << std::endl;
            exit(1);
        }
  

        // Set New X axis

        mod_tmp_v = DotP(v2, v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        if(mod_tmp_v)
        {
            for (int i = 0; i < dim; i++)
            {
                X_new [i] = 0.0;
                X_new [i] = v2[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The module of the vector v2 is zero " << std::endl;
            exit(1);
        }
        
        // Set the initial transformation matrix from old coordinate system
        // to the new system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Start grow an atom in the ring from the starting atom
       
        int  idxBond = getBond(tBonds, sSet[2], tIdx);
        int  idxAng  = getAngle(tAngles, sSet[2], sSet[1], tIdx);
        int  idxTor  = getTorsion(tTorsions, sSet[0], sSet[1], sSet[2], tIdx);
        
        if (idxBond !=-1 && idxAng !=-1 && idxTor !=-1)
        {
            REAL aBond=tBonds[idxBond].valueST;
            REAL aAng =tAngles[idxAng].valueST;
            REAL aTor;
            if (inRing)
            {
                //aTor=tTorsions[idxTor].value*PI180;
                
                if(tAtoms[sSet[2]].bondingIdx==3)
                {
                    if (iSeq >0)
                    {
                        aTor = -PI/3.0;
                    }
                    else
                    {
                        aTor =  PI/3.0;
                    }
                }
                else
                {
                    aTor = 0.0;
                 }
                 
            }
            else
            {  
                if(tAtoms[sSet[2]].bondingIdx==3)
                {
                    if (iSeq >0)
                    {
                        aTor = 2.0*PI/3.0;
                    }
                    else
                    {
                        aTor = -2.0*PI/3.0;
                    }
                }
                else
                {
                    //aTor = PI;
                    aTor=tTorsions[idxTor].value*PI180;
                }
                
                // aTor = PI;
            }
            if (std::find(doneList.begin(), doneList.end(), tIdx) ==doneList.end())
            {
                PolToCart(tAtoms, sSet[2], tIdx, aBond, aAng, aTor, A);
                tDoneSet.push_back(tIdx);
            }
            else
            {
                aTor = getTorsion(tAtoms[sSet[0]], tAtoms[sSet[1]], tAtoms[sSet[2]], tAtoms[tIdx])*PI180;
                // debug 
                //std::cout << "atom " << tAtoms[tIdx].id << "has already had coords "
                //          << " it is torsion value is " 
                //          << aTor
                //          << std::endl;
                
            }
            
            // build the atoms that are not on the ring but bonding with 
            // the atom sSet[2] (tIdx is also bonding with sSet[2] 
            
                       // grow other atoms linked to i_atom (like i_next)
            
                    if (tAtoms[sSet[2]].chiralIdx !=0)
                    {
                        int aCh = tAtoms[sSet[2]].inChirals[0];
                        bool lS= false;
                        std::vector<int> t1, rTable;
                        for (std::vector<int>::iterator irA=tChs[aCh].mutTable[sSet[1]].begin();
                               irA!=tChs[aCh].mutTable[sSet[1]].end(); irA++)
                        {
                            if (lS)
                            {
                                rTable.push_back(*irA);
                            }
                            else if (*irA==tIdx)
                            {
                                lS=true;
                            }
                            else
                            {
                                t1.push_back(*irA);
                            }
                        }
                        for (std::vector<int>::iterator irA2=t1.begin();
                                irA2 !=t1.end(); irA2++)
                        {
                            rTable.push_back(*irA2);
                        }
                        
                        //std::cout << "Chiral center " << tAtoms[sSet[2]].id << std::endl;
                        //std::cout << "chiral description  " << tChs[aCh].signST << std::endl;
                        //std::cout << "Atom arrangement " << std::endl;
                        //std::cout << "from atom " << tAtoms[sSet[1]].id  << " clockwise " << std::endl;
                        //std::cout << "atom " << tAtoms[tIdx].id << std::endl;
                        //for (std::vector<int>::iterator iNA=rTable.begin();
                        //        iNA !=rTable.end(); iNA++)
                        //{
                        //    std::cout << "Atom " << tAtoms[*iNA].id << std::endl;
                        //}
                        
                        int iL=1;
                        
                        for (std::vector<int>::iterator iNodeA=rTable.begin();
                                iNodeA !=rTable.end(); iNodeA++)
                        {
                            //std::cout << "iNodeA " << tAtoms[*iNodeA].id << std::endl;
                            
                            int bIdx    = getBond(tBonds, sSet[2], *iNodeA);
                            REAL tBo         = tBonds[bIdx].valueST;
                            int aIdx    = getAngle(tAngles, sSet[2], sSet[1], *iNodeA);
                            REAL tTha        = tAngles[aIdx].valueST;
                            //std::cout << "iL= " << iL << std::endl
                            //          << " tTor" << aTor << std::endl
                            //          << " deltal " << 2.0*PI/3.0 << std::endl;
                            REAL tPhiI  = aTor +iL*2.0*PI/3.0;
                            if (std::find(doneList.begin(), doneList.end(), *iNodeA) ==doneList.end())
                            {
                                // PolToCart(tAtoms,i_atom, *iNodeA, tBo, tTha, tPhiI,  A);
                                growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iNodeA);
                                doneList.push_back(*iNodeA);
                                //std::cout <<  tAtoms[*iNodeA].id << " done1 " << std::endl;
                            }
                            iL++;
                        }
                        // optSubSystem(doneList, tAtoms, tBonds, tAngles,
                        //             tTorsions, tRings, tPlas, tChs);
                        //std::cout << "The sub-system is optimized \n";
                    }
                    else
                    {
                        int iL=1;
                        for (std::vector<int>::iterator iCA=tAtoms[sSet[2]].connAtoms.begin();
                                iCA != tAtoms[sSet[2]].connAtoms.end(); iCA++)
                        {
                            if (*iCA !=tIdx && *iCA !=sSet[1] 
                                && std::find(doneList.begin(), doneList.end(), *iCA) ==doneList.end())
                            {
                                REAL tBo         = tAtoms[*iCA].treeBond;
                                REAL tTha        = tAtoms[*iCA].treeAngle;
                                if (tAtoms[sSet[2]].bondingIdx==3)
                                {    
                                    REAL tPhiI  = aTor +iL*2.0*PI/3.0;
                                    growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    //std::cout <<  tAtoms[*iCA].id << " done21 " << std::endl;
                                    iL++;
                                }
                                else if (tAtoms[sSet[2]].bondingIdx==2)
                                {
                                    REAL tPhiI  = aTor + PI;   
                                    growOneAtom(sSet, tAtoms, tBo, tTha, tPhiI, *iCA);
                                    // PolToCart(tAtoms,i_atom, *iCA, tBo, tTha, tPhiI,  A);
                                    doneList.push_back(*iCA);
                                    //std::cout <<  tAtoms[*iCA].id << " done22 " << std::endl;
                                }
                            }
                        }
                    }
            
            /*
            
            
            
            int iL=1;
            for (int iC=0; iC < (int)tAtoms[sSet[2]].connAtoms.size(); iC++)
            {
                int iAC=tAtoms[sSet[2]].connAtoms[iC];
                if( iAC !=sSet[1] && iAC !=tIdx
                    && std::find(tDoneSet.begin(), tDoneSet.end(), iAC)==tDoneSet.end())
                {
                    idxBond = getBond(tBonds, sSet[2], iAC);
                    idxAng  = getAngle(tAngles, sSet[2], sSet[1], iAC);
                    if (idxBond !=-1 && idxAng != -1)
                    {
                        aBond=tBonds[idxBond].value;
                        aAng =tAngles[idxAng].valueST;
                        REAL aTorI;
                        if (tAtoms[sSet[2]].bondingIdx==2)
                        {
                            aTorI = aTor + PI;
                        }
                        else if (tAtoms[sSet[2]].bondingIdx == 3)
                        {
                            aTorI = aTor +iL*2.0*PI/3.0;
                            iL++;
                        }
                        else
                        {
                            std::cout << "wrong connection in ring-building. Bug!"
                                      << std::endl;
                            exit(1);
                        }
                        
                        for (int i =0; i < dim; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                A[i][j] =0.0;
                            }
                        }
                        
                        InitAMat2(X_new,Y_new,Z_new,A, dim);
                        if (std::find(doneList.begin(), doneList.end(), iAC) ==doneList.end())
                        {
                            PolToCart(tAtoms, sSet[2], iAC, aBond, aAng, aTorI, A);
                            tDoneSet.push_back(iAC);
                        }
                    }
                    else
                    {
                        std::cout << "could not find either bond, or angle for atoms set "
                                  << tAtoms[sSet[1]].id << " " 
                                  << tAtoms[sSet[2]].id << " " 
                                  << tAtoms[iAC].id << std::endl;
                                  exit(1);
                    }
                }
            }
            */
            
        }
        else
        {
            std::cout << "could not find either bond, or angle or torsion for atoms set "
                    << tAtoms[sSet[0]].id << " " << tAtoms[sSet[1]].id 
                    << tAtoms[sSet[2]].id << " " << tAtoms[tIdx].id << std::endl;
            exit(1);
        }
        
        
        
        // Release memory and null point variables

        delete [] v1;
        v1 = 0;

        delete v2;
        v2 = 0;
  
        delete [] tmp_v;
        tmp_v = 0;
  
        delete [] X_new;
        X_new = 0;

        delete [] Y_new;
        Y_new = 0;

        delete [] Z_new;
        Z_new = 0;
        
        for (int i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete [] A;
        A = 0;
  
    }
    
    
    void TransCoords::ringBuilder(std::vector<int>         & sSet,
                                  RingDict                 & tRing,
                                  std::vector<AtomDict>    & tAtoms, 
                                  std::vector<BondDict>    & tBonds, 
                                  std::vector<AngleDict>   & tAngles, 
                                  std::vector<TorsionDict> & tTorsions,
                                  std::vector<RingDict>    & tRings,
                                  std::vector<PlaneDict>   & tPlas,
                                  std::vector<ChiralDict>  & tChs,
                                  int                        tTurn,
                                  std::vector<int>         & doneSet)
    {
        
        // Start grow the ring from the starting atom
        
        std::string tD;
        int iNext;
        
        //std::cout << "sSet[2]=" << sSet[2] << "  " << tAtoms[sSet[2]].id << std::endl
        //          << " c "  << tRing.atomsLink[sSet[2]]["c"] << "  " 
        //          << tAtoms[tRing.atomsLink[sSet[2]]["c"]].id 
        //          << " ac " << tRing.atomsLink[sSet[2]]["ac"] << "  " 
        //          << tAtoms[tRing.atomsLink[sSet[2]]["ac"]].id  << std::endl;
        
        //std::cout << "Done set : \n";
        //for (std::vector<int>::iterator iD=doneSet.begin();
        //        iD !=doneSet.end(); iD++)
        //{
        //    std::cout << *iD << "\t" << tAtoms[*iD].id << std::endl;
        //}
       
        std::vector<int> idxR;
        for (std::vector<AtomDict>::iterator iA=tRing.atoms.begin();
                iA !=tRing.atoms.end(); iA++)
        {
            idxR.push_back(iA->seriNum);
        }
        
        
        if (std::find(doneSet.begin(), doneSet.end(), tRing.atomsLink[sSet[2]]["c"]) ==doneSet.end() )
        {
            iNext = tRing.atomsLink[sSet[2]]["c"];
            tD="c";
        }
        else if (std::find(doneSet.begin(), doneSet.end(), tRing.atomsLink[sSet[2]]["ac"]) ==doneSet.end() )
        {
            iNext = tRing.atomsLink[sSet[2]]["ac"];
            tD = "ac";
        }
        else if (std::find(sSet.begin(), sSet.end(),tRing.atomsLink[sSet[2]]["ac"])==doneSet.end())
        {
            iNext = tRing.atomsLink[sSet[2]]["ac"];
            tD = "ac";
        }
        else if (std::find(sSet.begin(), sSet.end(),tRing.atomsLink[sSet[2]]["c"])==doneSet.end())
        {
            iNext = tRing.atomsLink[sSet[2]]["c"];
            tD = "c";
        }
        else
        {
            return ;
        }
       
      
        
        int iL=0, iSize=(int)tRing.atoms.size();
        //while (std::find(doneSet.begin(), doneSet.end(), iNext)==doneSet.end())
        int tRingAtmSeq =1;
        while (iL < iSize-1)
        {
            //std::cout << "ring atom grow  " << iNext << "  " << tAtoms[iNext].id << std::endl;
            std::vector<int> curSSet;
            setCurStartSet(sSet, curSSet, iNext, tD, tRing, doneSet);
            //growOneAtom(curSSet, tAtoms, tBonds, tAngles, tTorsions, iNext);
            //doneSet.push_back(iNext);
            bool tInR=false;
            if (std::find(idxR.begin(), idxR.end(), sSet[0]) !=idxR.end())
            {
                tInR=true;
            }
            //std::cout << "sSet[0] " << tAtoms[sSet[0]].id << std::endl;
            //std::cout << "tInR "  << tInR << std::endl;
            
            //if (std::find(doneSet.begin(), doneSet.end(), iNext)==doneSet.end())
            //{
            
            if (!tTurn)
            {
                //std::cout << "Grow node 1 " << std::endl;
                growOneRingNode(curSSet, tAtoms, tBonds, tAngles, tTorsions, tChs, tInR, tRingAtmSeq, iNext, doneSet);
            }
            else
            {
                //std::cout << "Grow node 2 " << std::endl;
                growOneRingNode2(curSSet, tAtoms, tBonds, tAngles, tTorsions, tChs, tInR, tRingAtmSeq, iNext, doneSet);
            }
          
            
            //}
            // tRingAtmSeq = -tRingAtmSeq;    
            sSet.erase(sSet.begin());
            sSet.push_back(iNext);
            iNext = tRing.atomsLink[iNext][tD];
            //std::cout << "ring atom next " << tAtoms[iNext].id << std::endl;
            iL++;
        }
        
        optSubSystem(doneList, tAtoms, tBonds, tAngles,
                     tTorsions, tRings, tPlas, tChs, false);
        //std::cout << "tree Node end: The sub-system is optimized \n";
       
        
        // Then grow the atoms linked to the last atom in the ring
        std::vector<int> curSSet;
        setCurStartSet(sSet, curSSet, iNext, tD, tRing, doneSet);
        bool tInR=true;
        //std::cout << "Grow atoms linked to the last atom in the ring" << std::endl;
        //std::cout << "sSet[0] " << tAtoms[sSet[0]].id << std::endl;
        growOneRingNode2(curSSet, tAtoms, tBonds, tAngles, tTorsions, tChs, tInR, tRingAtmSeq, iNext, doneSet);
        // std::cout << "done one ring " << std::endl;
        
    }
    
    void TransCoords::setSubSystem(std::vector<int>          & tDoneSet, 
                                   std::vector<AtomDict>     & tAllAtoms, 
                                   std::vector<BondDict>     & tAllBonds, 
                                   std::vector<AngleDict>    & tAllAngles, 
                                   std::vector<TorsionDict>  & tAllTorsions,
                                   std::vector<RingDict>     & tAllRings, 
                                   std::vector<PlaneDict>    & tAllPlanes,
                                   std::vector<ChiralDict>   & tAllChirals,
                                   AllSystem                 & tSubSys)
    {
        int i =0;
        std::map<int, int> idxCov;
        //std::cout << "in sub-system " << std::endl;
        
        // get atoms in sub-systems, re-set its serial number
        for (std::vector<int>::iterator iA=tDoneSet.begin();
                iA != tDoneSet.end(); iA++)
        {
            
            tSubSys.allAtoms.push_back(tAllAtoms[*iA]);
            //std::cout << "the " << i+1 << " atom is " << tSubSys.allAtoms[i].id << std::endl;
            tSubSys.allAtoms[i].seriNum = i;
            idxCov[*iA] = i;
            tSubSys.allAtoms[i].connAtoms.clear();
            //for (unsigned iS=0; iS < tSubSys.allAtoms[i].coords.size(); iS++)
            //{
            //    std::cout << "coords: " << tSubSys.allAtoms[i].coords[iS] << std::endl;
            //}
            i++;
        }
         
        
        
        i=0;
        for (std::vector<BondDict>::iterator iB=tAllBonds.begin();
                iB !=tAllBonds.end(); iB++)
        {
            if (std::find(tDoneSet.begin(), tDoneSet.end(), iB->atomsIdx[0]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iB->atomsIdx[1]) !=tDoneSet.end())
            {
                
                tSubSys.allBonds.push_back(*iB);
                tSubSys.allBonds[i].atomsIdx[0]=idxCov[tSubSys.allBonds[i].atomsIdx[0]];
                tSubSys.allBonds[i].atomsIdx[1]=idxCov[tSubSys.allBonds[i].atomsIdx[1]];
                tSubSys.allAtoms[tSubSys.allBonds[i].atomsIdx[0]].connAtoms.push_back(tSubSys.allBonds[i].atomsIdx[1]);
                tSubSys.allAtoms[tSubSys.allBonds[i].atomsIdx[1]].connAtoms.push_back(tSubSys.allBonds[i].atomsIdx[0]);
                i++;
            }   
        }
        
        /*
        std::cout << "number of bonds in the sub-system " 
                  << (int)tSubSys.allBonds.size()
                  << std::endl;
        for (std::vector<BondDict>::iterator iB=tSubSys.allBonds.begin();
                iB !=tSubSys.allBonds.end(); iB++)
        {
            std::cout << iB->valueST    << std::endl
                      << iB->sigValueST << std::endl;
        }
         */
       
        i = 0;
        for (std::vector<AngleDict>::iterator iAN=tAllAngles.begin();
                iAN !=tAllAngles.end(); iAN++)
        {
            if (std::find(tDoneSet.begin(), tDoneSet.end(), iAN->atoms[0]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iAN->atoms[1]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iAN->atoms[2]) !=tDoneSet.end())
            {
                tSubSys.allAngles.push_back(*iAN);
                tSubSys.allAngles[i].atoms[0] = idxCov[tSubSys.allAngles[i].atoms[0]];
                tSubSys.allAngles[i].atoms[1] = idxCov[tSubSys.allAngles[i].atoms[1]];
                tSubSys.allAngles[i].atoms[2] = idxCov[tSubSys.allAngles[i].atoms[2]];
                // std::cout << "dict value for angle " << i << " in the sub-system is "
                //           << tSubSys.allAngles[i].valueST << std::endl;
                i++;
            }
        }
        
        i =0;
        for (std::vector<TorsionDict>::iterator iT=tAllTorsions.begin();
                iT !=tAllTorsions.end(); iT++)
        {
            if (std::find(tDoneSet.begin(), tDoneSet.end(), iT->atoms[0]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iT->atoms[1]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iT->atoms[2]) !=tDoneSet.end()
                && std::find(tDoneSet.begin(), tDoneSet.end(), iT->atoms[3]) !=tDoneSet.end())
            {
                tSubSys.allTorsions.push_back(*iT);
                tSubSys.allTorsions[i].atoms[0] = idxCov[tSubSys.allTorsions[i].atoms[0]];
                tSubSys.allTorsions[i].atoms[0] = idxCov[tSubSys.allTorsions[i].atoms[1]];
                tSubSys.allTorsions[i].atoms[0] = idxCov[tSubSys.allTorsions[i].atoms[2]];
                tSubSys.allTorsions[i].atoms[0] = idxCov[tSubSys.allTorsions[i].atoms[3]];
                i++;
            }
        }
        
        //std::cout << "number of planes in the whole system " 
        //          << tAllPlanes.size() << std::endl;
        
        for (std::vector<PlaneDict>::iterator iP=tAllPlanes.begin();
                        iP !=tAllPlanes.end(); iP++)
        {
            std::map<ID, int> atomIdxs;
            
            for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                           iAt != iP->atoms.end(); iAt++)
            {
                if (std::find(tDoneSet.begin(), tDoneSet.end(), iAt->second) !=tDoneSet.end())
                {
                    
                    atomIdxs[iAt->first] = idxCov[iAt->second];
                }
            }
            
            if ((int)atomIdxs.size() >=3)
            {
                PlaneDict aPl;
                aPl.archID  = iP->archID;
                aPl.archPos = iP->archPos;
                for (std::map<ID,int>::iterator iPA=atomIdxs.begin();
                        iPA !=atomIdxs.end(); iPA++)
                {
                    aPl.atoms[iPA->first]=iPA->second;
                }
                tSubSys.allPlanes.push_back(aPl);
            }
        }
        
        //std::cout << "number of planes in the sub-system " 
        //          << (int)tSubSys.allPlanes.size() << std::endl;
        
        
    }
    
    void TransCoords::optSubSystem(std::vector<int>& tSubSet, 
                                   std::vector<AtomDict>    & tAllAtoms, 
                                   std::vector<BondDict>    & tAllBonds, 
                                   std::vector<AngleDict>   & tAllAngles, 
                                   std::vector<TorsionDict> & tAllTorsions,
                                   std::vector<RingDict>    & tAllRings,
                                   std::vector<PlaneDict>   & tAllPlanes,
                                   std::vector<ChiralDict>  & tAllChirals,
                                   bool useGO)
    {
        
        AllSystem  aSubSys;
        
        setSubSystem(tSubSet, tAllAtoms, tAllBonds, tAllAngles,
                     tAllTorsions, tAllRings, tAllPlanes, tAllChirals, aSubSys);
        
        //std::cout << "a subsystem is set, it contains " << std::endl
        //          << (int)aSubSys.allAtoms.size()    << " atoms "       << std::endl
        //          << (int)aSubSys.allBonds.size()    << " bonds "       << std::endl
        //          << (int)aSubSys.allAngles.size()   << " angles "      << std::endl
        //          << (int)aSubSys.allTorsions.size() << " torsions "    << std::endl
        //          << (int)aSubSys.allPlanes.size()   << " planes "      << std::endl;
                  // << (int)aSubSys.allRingsV.size() << " rings  "      << std::endl;
        
        // LIBMOL::outPDB("sub_tree.pdb", "XXX", aSubSys.allAtoms);
        
        
        
        if (useGO==false)
        {
            GO::FindLocalMin  toolLocalMin;
            toolLocalMin.workMode  = 1;
            toolLocalMin.workSpace = 1;
            toolLocalMin.lComp     = 100;
            toolLocalMin.Driver(aSubSys.allAtoms,   aSubSys.allBonds, 
                                aSubSys.allAngles,  aSubSys.allTorsions,
                                aSubSys.allRingsV,  aSubSys.allPlanes, 
                                aSubSys.allChirals);
            
        }
        else
        {
           GO::FindGlobMin toolGO(aSubSys);
           toolGO.Driver(true, 0);
        }
        //std::string tSubName = "sub_tree_" + IntToStr((int)aSubSys.allAtoms.size()) + "_opt.pdb";
        //LIBMOL::outPDB(tSubName.c_str(), "XXX", aSubSys.allAtoms);
       
        int tSize= (int)tAllAtoms.size();
        for (std::vector<AtomDict>::iterator iSA=aSubSys.allAtoms.begin();
                iSA != aSubSys.allAtoms.end(); iSA++)
        {
            int tPos = iSA->atomPosition(tAllAtoms);
            
            //std::cout << tAllAtoms[tPos].id << " and sub atom " << iSA->id << std::endl;
            //std::cout << "atom in subsystem coord size " << iSA->coords.size() << std::endl;
            //std::cout << "atom in all system coord size " << tAllAtoms[tPos].coords.size() << std::endl;
            if (tPos >=0 && tPos < tSize)
            {
                for(int i=0; i < (int)tAllAtoms[tPos].coords.size(); i++)
                {
                    tAllAtoms[tPos].coords[i]=iSA->coords[i];
                    //std::cout << iSA->coords[i] << std::endl;
                }
                    
            }
            else
            {
                std::cout << "The atom serial number of " << iSA->seriNum
                          << " does not exist "           << std::endl;
            }
        }
       
    }
   
}


