/* 
 * File:   atomsTree.h
 * Author: flong
 *
 * Created on February 27, 2013, 10:10 AM
 */

#include "atomsTree.h"

namespace LIBMOL
{
    buildAtomTree::buildAtomTree():startAtom(-1)
    {
    }
    
    buildAtomTree::~buildAtomTree()
    {
    }
    
    int buildAtomTree::setMSTStartAtom(std::vector<AtomDict>& allAtoms)
    {
        
        int iMax =-1;
        int nMax =0;
        //std::map<int, int> tAN;
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            int tN = allAtoms[i].getNumAtomsWithin2stNB(allAtoms);
            // std::cout << "tN " << tN << std::endl;
            if (allAtoms[i].id != "H" 
                && (int)allAtoms[i].ringRep.size() ==0
                && (int)allAtoms[i].connAtoms.size() > 1 
                && tN > nMax)
            {
                iMax = i;
                nMax = tN;
            }
        }
        
        if (iMax==-1)
        {
            nMax=0;
            // have to use a ring atom as the start atom
            for (int i=0; i < (int)allAtoms.size(); i++)
            {
                int tN = allAtoms[i].getNumAtomsWithin2stNB(allAtoms);
                if (allAtoms[i].id != "H" 
                    && (int)allAtoms[i].connAtoms.size() > 1 
                    && tN > nMax)
                {
                    iMax = i;
                }
            }
            if (iMax==-1)
            {
                std::cout << "Can not find an atom as the starting atom." << std::endl;
                std::cout << "They are either H or singly linked atoms " << std::endl;
                exit(1);
            }
            else
            {
                allAtoms[iMax].tree["parent"].push_back(-1);
            }
        }
        else
        {
            allAtoms[iMax].tree["parent"].push_back(-1);
        }
        
        return iMax;
       
    }
    void buildAtomTree::setAtomsMST(std::vector<AtomDict>& allAtoms, 
                               std::vector<BondDict>& allBonds)
    {
        std::vector<int> atmsInT, atmsOutT;
        
        int iS = setMSTStartAtom(allAtoms);
        if (iS <0)
        {
            exit(1);
        }
            
        // supplement variables
        std::map<int, std::map<int, REAL> >  allDs;
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            for (std::vector<int>::iterator iNB=allAtoms[i].connAtoms.begin();
                    iNB != allAtoms[i].connAtoms.end(); iNB++)
            {
                allDs[i][*iNB] = distanceV(allAtoms[i].coords, allAtoms[*iNB].coords);
            }
            
            /*
            if ((int)atmsInT.size() ==0)
            {
                // The starting atom of the tree should not be H and not in rings
                if (allAtoms[i].id != "H" && (int)allAtoms[i].ringRep.size() ==0
                        && (int)allAtoms[i].connAtoms.size() > 1)
                {
                    // start point of the tree
                    atmsInT.push_back(i);
                    allAtoms[i].tree["parent"].push_back(-1);
                    std::cout << "Atom " << allAtoms[i].id << " is the starting atom " << std::endl;
                }
                else
                {
                    atmsOutT.push_back(i);
                }
            }
            */
            
            if (i==iS)
            {
                atmsInT.push_back(i);
                std::cout << "Atom " << allAtoms[i].id << " is the starting atom " << std::endl;
            }
            else
            {
                atmsOutT.push_back(i);
            }
        }
        
        // double check
        if ((int)atmsInT.size() ==0)
        {
            std::cout << "Tree start atom error " << std::endl;
            exit(1);
        }
        
        
        while ((int)atmsOutT.size() !=0)
        {
            // std::cout << (int)atmsOutT.size() << " atoms out of tree yet " << std::endl;
            int iP=-1, iC=-1; 
            REAL dMin =2000.0;
            std::vector<int>::iterator iFind;
            for (std::vector<int>::iterator iTA=atmsInT.begin();
                    iTA !=atmsInT.end(); iTA++)
            {
                for (std::map<int, REAL>::iterator iNB=allDs[*iTA].begin();
                        iNB !=allDs[*iTA].end(); iNB++)
                {
                    std::vector<int>::iterator tFind = 
                       std::find(atmsOutT.begin(), atmsOutT.end(), iNB->first);
                    if (iNB->second < dMin && 
                         tFind!=atmsOutT.end())
                    {
                        iP=*iTA;
                        iC=iNB->first;
                        dMin = iNB->second;
                        iFind = tFind;
                    }
                }
            }
            
            // std::cout << "Parent " << iP << " child " << iC << std::endl;
            if (iP != -1 && iC !=-1)
            {
                // get one vertex of the tree
                allAtoms[iP].tree["children"].push_back(iC);
                allAtoms[iC].tree["parent"].push_back(iP);
                atmsInT.push_back(iC);
                std::cout << "Atom " << allAtoms[iP].id << " has a child " << allAtoms[iC].id 
                        << std::endl;
           
                std::cout << "Atom " << allAtoms[iC].id << " has parent " << allAtoms[iP].id << std::endl;
                // std::cout << " number of parent " << (int)allAtoms[iP].tree["parent"].size() << std::endl;
                
                if ((int)allAtoms[iC].tree["parent"].size() > 1)
                {
                    std::cout << "Atom " << allAtoms[iC].id << " has "
                            << (int)allAtoms[iC].tree["parent"].size() << std::endl
                            << "They are \t";
                    for (std::vector<int>::iterator iPA=allAtoms[iC].tree["parent"].begin();
                            iPA != allAtoms[iC].tree["parent"].begin(); iPA++)
                    {
                        std::cout << *iPA << "\t";
                    }
                    std::cout << std::endl;
                    
                    exit(1);            
                }
                /*
                else
                {
                    std::cout << "Atom " << allAtoms[iP].id << " has "
                            << (int)allAtoms[iP].tree["children"].size() 
                            << " children. They are \n";
                    for (std::vector<int>::iterator iPA=allAtoms[iP].tree["children"].begin();
                            iPA != allAtoms[iP].tree["children"].end(); iPA++)
                    {
                        std::cout << *iPA << "\t";
                    }
                    std::cout << std::endl;
                    
                }
                 */
                       
               
                // remove it for out of tree list and allDs
               
                atmsOutT.erase(iFind);
                
                allDs[iP].erase(iC);
                
            }
            else
            {
                std::cout << "why no further vertex found " << std::endl;
                std::cout << "Current out-of-tree atoms are " << std::endl;
                for (std::vector<int>::iterator iNTA=atmsOutT.begin();
                        iNTA !=atmsOutT.end(); iNTA++)
                {
                    std::cout << *iNTA << "\t";
                }
                std::cout << std::endl;          
            }
        }
       
        
        // Output tree for a look
        std::cout <<"******* Tree structure for all atoms " << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin(); 
                iAt != allAtoms.end(); iAt++)
        {
            
            
            if (iAt->tree.find("children") !=iAt->tree.end())
            {
                
                ID tPar;
                if (iAt->tree["parent"][0] !=-1)
                {
                    tPar = allAtoms[iAt->tree["parent"][0]].id; 
                }
                else
                {
                    tPar = "BEGIN";
                }
                
                for (std::vector<int>::iterator iNA=iAt->tree["children"].begin();
                         iNA != iAt->tree["children"].end(); iNA++)
                {
                    std::cout << iAt->id << "\t"  << tPar << "\t"
                          << allAtoms[*iNA].id << std::endl;
                }
                
            }
            else
            {
                 
                std::cout << iAt->id << "\t" << allAtoms[iAt->tree["parent"][0]].id 
                            << "\t" << "END" << std::endl; 
            }
        }
        
    }
    
    void buildAtomTree::initDummyAtoms(std::vector<AtomDict>& tAllAtoms, 
                                       std::vector<AtomDict>& tDAtoms)
    {
        
  
        int i_atom1, i_atom2, i_atom3;
        
    
        int baseAtomNum = startAtom;
        if (baseAtomNum < 0)
        {
            std::cout << "Could not find the root atom for the tree" 
                      << std::endl;
            exit(1);
        }
        i_atom3 = baseAtomNum;
        std::cout << "BaseAtom is " << i_atom3+1 << std::endl;

        
        for (std::vector<int>::iterator iA=tAllAtoms[i_atom3].tree["children"].begin();
                iA !=tAllAtoms[i_atom3].tree["children"].end(); iA++)
        {
            i_atom2 = -1;
            i_atom1 = -1;
            if (tAllAtoms[*iA].id.at(0) !='H')
            {
                for (std::vector<int>::iterator iA1=tAllAtoms[*iA].tree["children"].begin();
                     iA1 !=tAllAtoms[*iA].tree["children"].end(); iA1++)
                {
                    if(tAllAtoms[*iA1].id.at(0) != 'H')
                    {
                        i_atom1 = *iA1;
                        i_atom2 = *iA;
                        break;
                    }
                }
                if (i_atom1 != -1 && i_atom2 !=-1)
                {
                    break;
                }
            }
        }
        
        if(i_atom2 < 0) 
        {
            std::cout << "without the forward atom for atom " 
                      << tAllAtoms[i_atom3].id << std::endl;
            exit(1);
        }
  
        if(i_atom1 < 0) 
        {
            std::cout << "without the forward atom for atom " 
                      << tAllAtoms[i_atom2].id << std::endl;
            exit(1);
        }
        
        // Initiate dummy atom number 2
        AtomDict tDAtom2;
        tDAtom2.seriNum = -1;
        
        std::vector<AtomDict> atoms123;
        
        atoms123.push_back(tAllAtoms[i_atom1]);
        atoms123.push_back(tAllAtoms[i_atom2]);
        atoms123.push_back(tAllAtoms[i_atom3]);
     
        setOneDummyAtom(atoms123, tDAtom2);
        
        
        // Initiate dummy atom number 1
        AtomDict tDAtom1;
        tDAtom1.seriNum = -2;
        std::vector<AtomDict> atomsd23;
        atomsd23.push_back(tAllAtoms[i_atom2]);
        atomsd23.push_back(tAllAtoms[i_atom3]);
        atomsd23.push_back(tDAtom2);
        
        setOneDummyAtom(atomsd23, tDAtom1);
        
        tDAtoms.push_back(tDAtom1);
        tDAtoms.push_back(tDAtom2);
        
    }
    
    void buildAtomTree::setOneDummyAtom(std::vector<AtomDict>& t3Atoms, 
                                        AtomDict& tDAtom)
    {

        // cout << " numDummyTors " <<numDummyTors << endl;
        // cout << " i_atom1 = " <<  i_atom1+1 << endl;
        // cout << " i_atom2 = " <<  i_atom2+1 << endl;
        // cout << " i_atom3 = " <<  i_atom3+1 << endl;
 
        int dim = (int)t3Atoms[0].coords.size();
       
        int i,j;
        
        REAL *v1 = new REAL [dim];
        REAL *v2 = new REAL [dim];
        REAL *v3 = new REAL [dim];
        REAL *tmp_v = new REAL [dim];
        REAL mod_tmp_v = 0.0;

        REAL *X_new = new REAL [dim];
        REAL *Y_new = new REAL [dim]; 
        REAL *Z_new = new REAL [dim];
        REAL ** A = new REAL * [dim];
        for (i =0; i < dim; i++)
        {
            A[i] = new REAL [dim];
        }
        for (i = 0; i < dim; i++)
        {
            v1[i]    = 0.0;
            v2[i]    = 0.0;
            v3[i]    = 0.0;
            tmp_v[i] = 0.0;
            X_new[i] = 0.0;
            Y_new[i] = 0.0;
            Z_new[i] = 0.0;
            for (j = 0; j < dim; j++)
            {
                A[i][j] =0.0;
            }
        }
        
        for (i =0; i < dim; i++)
        {
            v1[i] = t3Atoms[1].coords[i] - t3Atoms[0].coords[i];
            v2[i] = t3Atoms[2].coords[i] - t3Atoms[1].coords[i];
        }

        // Set new X axis (unit vector of v2)

        mod_tmp_v = DotP(v2,v2);
        mod_tmp_v = sqrt(mod_tmp_v);
        
        if(mod_tmp_v)
        {
            for ( i = 0; i < dim; i++)
            {
                X_new [i] = v2[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " << std::endl;
            exit(1);
        }

  
        // Set new Z axis 

        CrossP(v1,v2,tmp_v);
 
        mod_tmp_v = DotP(tmp_v,tmp_v);
  
        mod_tmp_v = sqrt(mod_tmp_v);

        if(mod_tmp_v > 1.0e-16)
        {
            for ( i = 0; i < dim; i++)
            {
                Z_new [i] = tmp_v[i]/mod_tmp_v; 
            }
        }
        else 
        {
            std::cout << " The cross product of vectors v1 and v2 is zero " 
                      << std::endl;
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
                Y_new [i] = tmp_v[i]/mod_tmp_v;
            }
        }
        else 
        {
            std::cout << " The cross product of vectors Z_new and v2 is zero " 
                      << std::endl;
            exit(1);
        }
  
        // Set the initial transformation matrix from new coordinate system
        // to the fixed coordinate system

        InitAMat2(X_new,Y_new,Z_new,A, dim);

        // Now get the change of the coordinats of the first 
        // dummy atom in the new frame 

        REAL length = 1.0;
  
        REAL th_i, c_th_i, s_th_i;
        REAL ph_i, c_ph_i, s_ph_i;
        REAL * dX;
        REAL * dY; 
        REAL ddX;

        dX = new REAL [dim];
        dY = new REAL [dim];

        th_i  = PI/2.0;
        
        //ph_i = tAllAtoms[i_atom1].treeTorsAng_st;
        // Now assume dummy atom has a torsion angle 0.0 with i_atom1;
        ph_i = 0.0; 
       
   
        c_th_i = cos(th_i);
        s_th_i = sin(th_i);
        c_ph_i = cos(ph_i);
        s_ph_i = sin(ph_i);

        if(fabs(c_th_i)<1.0e-6)
        {
            c_th_i = 0.0;
        }
        else if((1-fabs(c_th_i))<1.0e-6)
        {
            c_th_i = Sign(1.0,c_th_i);
        }
        
        if(fabs(s_th_i)<1.0e-6)
        {
            s_th_i = 0.0;
        }
        else if((1-fabs(s_th_i))<1.0e-6)
        {
            s_th_i = Sign(1.0,s_th_i);
        }

        if(fabs(c_ph_i)<1.0e-6)
        {
            c_ph_i = 0.0;
        }
        else if ((1-fabs(c_ph_i))<1.0e-6)
        {
            c_ph_i = Sign(1.0,c_ph_i);
        }

        if(fabs(s_ph_i)<1.0e-6)
        {
            s_ph_i = 0.0;
        }
        else if((1-fabs(s_ph_i))<1.0e-6)
        {
            s_ph_i = Sign(1.0,s_ph_i);
        }

        // dX, the variations in the new coordinate system

        dX[0] =  c_th_i*length;
        ddX   =   s_th_i*length;
        dX[1] =   c_ph_i*ddX;  
        dX[2] =   s_ph_i*ddX;

        //   cout << " th_i is " << th_i << " cos(th_i) " << c_th_i << endl;
        //   cout << " dX[0] = " << dX[0] << endl;
        //   cout << " ph_i is " << ph_i << " cos(ph_i) " << c_ph_i << endl;
        //   cout << " dX[1] = " << dX[1] << " dX[2] = " << dX[2] <<  endl;

        // dY, the variation in new coordinate system 

        int n_dim_t = dim;

        MatMultip(n_dim_t, A, dX, dY);

        for (i =0; i < dim; i++)
        {
            if(fabs(dY[i]) < 1.0e-8)
            {
                dY[i] = 0.0;
            }
            //     cout << " dX["<< i <<"]= " << dX[i] << endl; 
            //     cout << " dY["<< i <<"]= " << dY[i] << endl;
        }

        // The position of the first dummy atom 
 
        // cout << " Before transfering " << endl;
        // for (i =0; i < dim; i++)
        //  {
        //    cout << "tAllAtoms[" << i_atom3+1 << "].coord["<<i+1<<"] = " 
        //            << tAllAtoms[i_atom3].coords[i] << endl; 
        //   }

        //  cout <<" after transfer " << endl;

        for ( i =0; i < dim; i++)
        {
            tDAtom.coords[i] 
                = t3Atoms[2].coords[i]+ dY[i];
            //cout << " Dummy atom 1.coord["<<i+1<<"] = " 
            //     << tAllAtomsDummy[0].coords[i] << endl;
        }


        // Chech if this dummy atom gives the correct torsion angle 1

        // i_atom1 = baseAtomNum -1;

        // i_atom2 = tAllAtoms[i_atom1].frwAtmSerNum -1;
  
        // i_atom3 = tAllAtoms[i_atom2].frwAtmSerNum -1;

        // i = numDummyTors;

        //i_atom1 = torsAngs[i].compAtomNum[1]-1;
        //i_atom2 = torsAngs[i].compAtomNum[2]-1;
        //i_atom3 = torsAngs[i].compAtomNum[3]-1;

        // cout << " i_atom1 = " <<  i_atom1+1 << endl;
        // cout << " i_atom2 = " <<  i_atom2+1 << endl;
        // cout << " i_atom3 = " <<  i_atom3+1 << endl;

        // for (i =0; i < dim; i++)
        //  {
        //    v1[i] = tAllAtoms[i_atom1].coords[i] - tAllAtomsDummy[0].coords[i];
        //    v2[i] = tAllAtoms[i_atom2].coords[i] - tAllAtoms[i_atom1].coords[i];
        //    v3[i] = tAllAtoms[i_atom3].coords[i] - tAllAtoms[i_atom2].coords[i];
        //  }

        // REAL torsion;

        // torsion = GetTorsAngle(v1, v2, v3);

        // cout << "The first torsion angle value is " << torsion << endl;

        delete [] dX;
        dX = 0;
        delete [] dY;
        dY = 0;
 
        delete []  v1;
        v1 = 0;
        delete []  v2;
        v2 = 0;
        delete []  v3;
        v3 = 0;
 
        delete [] tmp_v;
        tmp_v = 0; 
  
        delete []  X_new;
        X_new = 0;
        delete []  Y_new;
        Y_new = 0;
        delete []  Z_new;
        Z_new = 0;
        
        for (i =0; i < dim; i++)
        {
            delete [] A[i];
            A[i] = 0;
        }
        delete []  A;
        A = 0;
    
    }
    
    void buildAtomTree::setStartAtom(std::vector<AtomDict> & tAllAtoms)
    {
        int pos = -1;
        
        for ( int i=0; i < (int)tAllAtoms.size(); i++)
        {
            if (tAllAtoms[i].tree["parent"][0] == -1)
            {
                pos = i;
                break;
            }
        }
        
        startAtom=pos;
        
    }
    
    // Note: set-up of Torsion angles is different from those of bonds and 
    // bond-angles
    // bonds and bond-angles: set-up the current atom's values
    // torsion angles : set-up their children's values
    
    void buildAtomTree::setTreeAtomValues(std::vector<AtomDict>    & tAtoms, 
                                          std::vector<BondDict>    & tBonds, 
                                          std::vector<AngleDict>   & tAngles, 
                                          std::vector<TorsionDict> & tTorsions)
    {
        int tS = startAtom;
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            // std::cout << "For atom " << iAt->id << std::endl;
            int idxP = iAt->tree["parent"][0];
            // std::cout << ", its parent  " << idxP << std::endl;
            if (idxP >=0)
            {
                // std::cout << tAtoms[iAt->tree["parent"][0]].id << std::endl;
                setTreeAtomBondValue(tAtoms, tBonds, iAt);
                int idxGp=tAtoms[idxP].tree["parent"][0];
                // std::cout << " its grand-p " << idxGp << std::endl;
                if (idxGp >=0)
                {
                    //std::cout << tAtoms[idxGp].id << std::endl;
                    setTreeAtomAngleValue(tAtoms, tAngles, iAt);
                    int idxGgp = tAtoms[idxGp].tree["parent"][0];
                    //std::cout << "its great grand-pa " << idxGgp << std::endl;
                    //std::cout << "1st " << iAt->seriNum << " 2nd " 
                    //          << idxP << " 3rd " << idxGp 
                    //          << " 4th " << idxGgp << std::endl;
                    if (idxGgp >=0)
                    {
                        //std::cout << "Target: " << std::endl 
                        //          << "atom 1 " << iAt->id << "\t atom 2 " << tAtoms[idxP].id 
                        //          << "\t atom 3 " << tAtoms[idxGp].id
                        //          << "\t atom4 "  << tAtoms[idxGgp].id 
                        //         << std::endl;
                        setTreeAtomTorsionValue(tAtoms, tTorsions,  iAt, -1);
                    }
                    else
                    {
                        int tRoot=-1;
                        for (int i=0; i < (int)startPack.size(); i++)
                        {
                            
                            if (startPack[i] != tS && startPack[i] !=idxP )
                            {
                                tRoot = startPack[i];
                                break;
                            }
                           
                        }
                       
                        if(tRoot >=0)
                        {
                            //std::cout << "Target: " << std::endl 
                            //          << "atom 1 " << iAt->id 
                            //          << "\t atom 2 " << tAtoms[idxP].id 
                            //          << "\t atom 3 " << tAtoms[idxGp].id 
                            //          << "\t atom 4 " << tAtoms[tRoot].id 
                            //          << std::endl;
                            setTreeAtomTorsionValue(tAtoms, tTorsions, iAt, tRoot);    
                        }
                        else
                        {
                            std::cout << "Could not find the start atom for "
                                      << iAt->id << std::endl;
                        }
                    }
                }
                else if (idxGp==-1)
                {
                    // Now no need to set these torsions
                    //setTreeAtomTorsionValue(tAtoms, tTorsions, tPlas, tChs, tDAts, iAt);
                }
            }
            
            std::cout << "atom id " << iAt->id << std::endl
                      << "Its tree bond " << iAt->treeBond 
                      << "\tIts tree angle " << iAt->treeAngle 
                      << "\tIts tree torsion " << iAt->treeTorsion << std::endl; 
               
        }
    }
    
    void buildAtomTree::setTreeAtomValues(std::vector<AtomDict>& tAtoms)
    {
        
    }
    
    void buildAtomTree::setTreeAtomBondValue(std::vector<AtomDict>& tAtoms, 
                                             std::vector<BondDict>& tBonds, 
                                             std::vector<AtomDict>::iterator tAt)
    {
        int tBo = getBond(tBonds, tAt->seriNum, tAt->tree["parent"][0]);
        if (tBo >=0)
        {
            tAt->treeBond = tBonds[tBo].value;
        } 
        else
        {
            std::cout << "Could not find the bond between " 
                      << tAt->id << " and " 
                      << tAtoms[tAt->tree["parent"][0]].id << std::endl;
            exit(1);
        }
    }
    
    void buildAtomTree::setTreeAtomAngleValue(std::vector<AtomDict>& tAtoms, 
                                              std::vector<AngleDict>& tAngs,
                                              std::vector<AtomDict>::iterator tAt)
    {
        int idxP  = tAt->tree["parent"][0];
        int idxGP = tAtoms[idxP].tree["parent"][0];
        int tAn = getAngle(tAngs, idxP, tAt->seriNum, idxGP); 
        if (tAn >=0)
        {
            tAt->treeAngle = tAngs[tAn].value*PI180; 
        }
        else
        {
            std::cout << "Could not find the angle between " 
                      << tAtoms[idxP].id << "(center) and " << tAt->id 
                      << " and " << tAtoms[idxGP].id << std::endl;
            exit(1);
        }
    }
    
    void buildAtomTree::setTreeAtomTorsionValue(std::vector<AtomDict>& tAtoms, 
                                                std::vector<TorsionDict>& tTors,  
                                                std::vector<AtomDict>::iterator tAt,
                                                int   tR)
    {
        int P  = tAt->tree["parent"][0];
        int Gp = tAtoms[P].tree["parent"][0];
        int Ggp = -1;
        //  bool rev = false;
        int curT;
        if (tR == -1)
        {
            Ggp = tAtoms[Gp].tree["parent"][0];
        }
        else
        {
            Ggp = tR;
        }
        curT=getTorsion(tTors, tAt->seriNum, P, Gp, Ggp);
            
        if (curT >=0)
        {
            if (tAt->seriNum == tTors[curT].atoms[0] && Ggp==tTors[curT].atoms[3])
            {
                tAt->treeTorsion = tTors[curT].value*PI180;
            }
            else if (tAt->seriNum == tTors[curT].atoms[3] && Ggp==tTors[curT].atoms[0])
            {
                tAt->treeTorsion = -tTors[curT].value*PI180;
            }
            else
            {
                std::cout << "Two ending atoms in torsion angle are not right for torsion " 
                    << tAt->id  << std::endl
                    << " 4 atoms are: 1 " << tAt->id << " 2 " << tAtoms[P].id 
                    << " 3 " << tAtoms[Gp].id << " 4 " << tAtoms[Ggp].id << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout << "No torsion angle is found for atoms " 
                    << tAtoms[Ggp].id << "->" << tAtoms[Gp].id 
                    << "->" << tAtoms[P].id << "->" << tAt->id << std::endl;
            exit(1);
        }
        
    }
    
    void buildAtomTree::setTreeAtomTorsionValue(std::vector<AtomDict>& tAtoms, 
                                                std::vector<TorsionDict>& tTors, 
                                                AtomDict& tDAt1, 
                                                std::vector<AtomDict>::iterator tAt)
    {  
    }
    
    void buildAtomTree::setTreeAtomTorsionValue(std::vector<AtomDict> & tAtoms, 
                                                std::vector<TorsionDict> & tTors, 
                                                std::vector<PlaneDict> & tPlas, 
                                                std::vector<ChiralDict>& tChirals, 
                                                std::vector<AtomDict> & tDAts, 
                                                std::vector<AtomDict>::iterator tAt)
    {
        
    }
    
    void buildAtomTree::setBranches(std::vector<AtomDict>& tAtoms)
    {
        if (startAtom <0)
        {
            std::cout << "Could not find the starting atom for the tree, check! "
                    << std::endl;
            exit(1);
        }
        
        for (std::vector<int>::iterator iC=tAtoms[startAtom].tree["children"].begin();
                iC !=tAtoms[startAtom].tree["children"].end(); iC++)
        {
            int idxRoot = -1;
            if (*iC !=startPack[0])
            {
                idxRoot = startPack[0];
            }
            else if (*iC !=startPack[1])
            {
                idxRoot = startPack[1];
            }
            else
            {
                std::cout << "The current child is " << *iC <<std::endl
                        << "Two atoms treated as roots are " 
                        << startPack[0] << " and " << startPack[1] 
                        << std::endl;
                exit(1);
            }
            
            branches[idxRoot][startAtom].push_back(*iC);
            
            std::cout << "branch atoms are " << tAtoms[idxRoot].id << " and "
                      << tAtoms[startAtom].id << " and " << tAtoms[*iC].id 
                      << std::endl;
        } 
    }
    
    void buildAtomTree::setStartStruct(std::vector<AtomDict>  & tAtoms,
                                       std::vector<BondDict>  & tBonds,
                                       std::vector<AngleDict> & tAngs,
                                       std::vector<ChiralDict>& tChs)
    {
        setStartAtom(tAtoms);
        
        for (std::vector<int>::iterator iAt=tAtoms[startAtom].tree["children"].begin();
                iAt !=tAtoms[startAtom].tree["children"].end(); iAt++)
        {
            if((int)startPack.size() < 2)
            {
                // First, do not include H in the start-set of atoms.
                if(tAtoms[*iAt].chemType != "H")
                {
                    startPack.push_back(*iAt);
                }
            }
        }
        
        if ((int)startPack.size() < 2)
        {
            // The number of start-set of atoms is still less 3, use H now.
            for (std::vector<int>::iterator iAt=tAtoms[startAtom].tree["children"].begin();
                         iAt !=tAtoms[startAtom].tree["children"].end(); iAt++)
            {
                if((int)startPack.size() < 2 && 
                   std::find(startPack.begin(), startPack.end(), *iAt) ==startPack.end())
                {
                    startPack.push_back(*iAt);
                }
            }
        }
        
        std::cout << "atom in startPack " << std::endl;
        
        for(std::vector<int>::iterator iS=startPack.begin();
                iS !=startPack.end(); iS++)
        {
            std::cout << "atom " << tAtoms[*iS].id << std::endl;
        }
        
        
        if((int)startPack.size() < 2)
        {
            std::cout << "Why there are still less than 2 atoms in the start-set "
                    << std::endl;
            exit(1);
        }
        
        // Now build the structure for the start-set of atoms
        
        for (std::vector<REAL>::iterator iC=tAtoms[startAtom].coords.begin();
                iC!=tAtoms[startAtom].coords.end(); iC++)
        {
            *iC= 0.0;
        }
        
        int nChs = (int)tAtoms[startAtom].tree["children"].size();
        std::cout << "The tree starts on atom " << tAtoms[startAtom].id 
                << " has " << nChs << " branches " << std::endl;
        
        if(nChs==2)
        {
            // two children do not mean sp1, depending on 
            // the element type of root atom
            
            set2ConnStruct(tAtoms, tBonds, tAngs);
        }
        else if (nChs==3)
        {
            set3ConnStruct(tAtoms, tBonds, tAngs, tChs);
        }
        else if(nChs==4)
        {
            set4ConnStruct(tAtoms, tBonds, tAngs, tChs);
        }
        else
        {
            // to be studied
        }  
        
        setBranches(tAtoms);
        
    }
    
    // The center atom may not be in a straight line with two connected atoms,
    // e.g., if the center atom is an O atom. Using the angle between them. 
    void buildAtomTree::set2ConnStruct(std::vector<AtomDict>  & tAtoms,
                        std::vector<BondDict>  & tBonds,
                        std::vector<AngleDict> & tAngs)
    {
        
        
        int iAt1 = tAtoms[startAtom].tree["children"][0];
        int iAt2 = tAtoms[startAtom].tree["children"][1];
       
        
        int iBo1 = getBond(tBonds, startAtom, iAt1);
        
        
        if (iBo1 < 0)
        {
            std::cout << "There is no bond between atom " 
                    << tAtoms[startAtom].id << " and atom  "
                    << tAtoms[iAt1].id << ". Check what went wrong!"
                    << std::endl;
            exit(1);
                   
        }
        
        int iBo2 = getBond(tBonds, startAtom, iAt2);
        if (iBo2 < 0)
        {
            std::cout << "There is no bond between atom " 
                    << tAtoms[startAtom].id << " and atom  "
                    << tAtoms[iAt2].id << ". Check what went wrong!"
                    << std::endl;
            exit(1);
                   
        }
       
        
        int iAng = getAngle(tAngs, startAtom, iAt1, iAt2);
       
        tAtoms[iAt1].coords[0] = tBonds[iBo1].value;
        tAtoms[iAt1].coords[1] = 0.0;
        tAtoms[iAt1].coords[2] = 0.0;
       
        tAtoms[iAt2].coords[0] = tBonds[iBo2].value*cos(tAngs[iAng].value*PI180);
        tAtoms[iAt2].coords[1] = tBonds[iBo2].value*sin(tAngs[iAng].value*PI180);
        tAtoms[iAt2].coords[2] = 0.0;
        
    }
    
    // using sp2, sp3 index for structure construction 
    void buildAtomTree::set3ConnStruct(std::vector<AtomDict>  & tAtoms,
                        std::vector<BondDict>  & tBonds,
                        std::vector<AngleDict> & tAngs,
                        std::vector<ChiralDict>& tChs)
    {
        int iAt1 = tAtoms[startAtom].tree["children"][0];
        int iAt2 = tAtoms[startAtom].tree["children"][1];
        int iAt3 = tAtoms[startAtom].tree["children"][2];
        
        std::cout << "starting atom " << tAtoms[startAtom].id << std::endl;
        std::cout << "Three children are " << tAtoms[iAt1].id << " and "
                  << tAtoms[iAt2].id << " and " << tAtoms[iAt3].id << std::endl;
        
        
        
        int iBo1 = getBond(tBonds, startAtom, iAt1);
        if (iBo1 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[iAt1].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);
                   
        }
        
        int iBo2 = getBond(tBonds, startAtom, iAt2);
        if (iBo2 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[iAt2].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);      
        }
            
        int iBo3 = getBond(tBonds, startAtom, iAt3);
        if (iBo3 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[iAt3].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);      
        }
            
        int iAng1 = getAngle(tAngs, startAtom, iAt1, iAt2);
        int iAng2 = getAngle(tAngs, startAtom, iAt1, iAt3);
        
        //std::cout << "angle 1 " << tAngs[iAng1].valueST << std::endl 
        //          << "angle 2 " << tAngs[iAng2].valueST << std::endl;
        
       

        
        if (tAtoms[startAtom].bondingIdx ==2)
        {
            // sp2, four atoms in a plane
            tAtoms[iAt1].coords[0] = tBonds[iBo1].value;
            tAtoms[iAt1].coords[1] = 0.0;
            tAtoms[iAt1].coords[2] = 0.0;
            
            tAtoms[iAt2].coords[0] = tBonds[iBo2].value*cos(tAngs[iAng1].value*PI180);
            tAtoms[iAt2].coords[1] = tBonds[iBo2].value*sin(tAngs[iAng1].value*PI180);
            tAtoms[iAt2].coords[2] = 0.0;
            
            tAtoms[iAt3].coords[0] =  tBonds[iBo3].value*cos(tAngs[iAng2].value*PI180);
            tAtoms[iAt3].coords[1] = -tBonds[iBo3].value*sin(tAngs[iAng2].value*PI180); // a minus is put here
            tAtoms[iAt3].coords[2] = 0.0;
            
        }
        else if (tAtoms[startAtom].bondingIdx ==3)
        {
            // sp3, in some cases, atom has 3 bonds but is still in a sp3 structure
            // center atom is not at a plane with other 3 atoms
            tAtoms[iAt1].coords[0] =  tBonds[iBo1].value*cos(19.5*PI180);
            tAtoms[iAt1].coords[1] =  0.0;
            tAtoms[iAt1].coords[2] =  tBonds[iBo1].value*sin(19.5*PI180);
            
            tAtoms[iAt2].coords[0] =  tBonds[iBo2].value*cos(19.5*PI180)*cos(tAngs[iAng1].value*PI180);
            tAtoms[iAt2].coords[1] =  tBonds[iBo2].value*cos(19.5*PI180)*sin(tAngs[iAng1].value*PI180);
            tAtoms[iAt2].coords[2] =  tBonds[iBo2].value*sin(19.5*PI180);
            
            tAtoms[iAt3].coords[0] =  tBonds[iBo3].value*cos(19.5*PI180)*cos(tAngs[iAng2].value*PI180);
            tAtoms[iAt3].coords[1] = -tBonds[iBo3].value*cos(19.5*PI180)*sin(tAngs[iAng2].value*PI180);
            tAtoms[iAt3].coords[2] =  tBonds[iBo3].value*sin(19.5*PI180);
        }
        else
        {
            std::cout << "sp orbit hybr index for atom " << tAtoms[startAtom].id 
                      << " is " << tAtoms[startAtom].bondingIdx << std::endl;
            exit(1);
        }
        
        
    }
    
    void buildAtomTree::set4ConnStruct(std::vector<AtomDict>  & tAtoms,
                        std::vector<BondDict>  & tBonds,
                        std::vector<AngleDict> & tAngs,
                        std::vector<ChiralDict>& tChs)
    {
        //      Need chiral properties 
        //              2
        //              |
        //            3-A-1
        //              |
        //              4
        
        std::vector<int> Ats;
        Ats.push_back(tAtoms[startAtom].tree["children"][0]);
        std::cout << "start atom " << tAtoms[startAtom].id 
                  << " is in " << (int)tAtoms[startAtom].inChirals.size() << std::endl;
        if ((int)tAtoms[startAtom].inChirals.size() !=0)
        {
            // the start atom is chiral center 
            int iCh = tAtoms[startAtom].inChirals[0];
            
            for (std::vector<int>::iterator iA=tChs[iCh].mutTable[Ats[0]].begin();
                    iA !=tChs[iCh].mutTable[Ats[0]].end(); iA++)
            {
                Ats.push_back(*iA);
            }
            
        }
        else
        {
            // the start atom is not chiral center (connected to more than 2 H e.g.)
            for (int i=1; i< (int)tAtoms[startAtom].tree["children"].size(); i++)
            {
                Ats.push_back(tAtoms[startAtom].tree["children"][i]);
            }
        }
        
        
        /*
        int iAt2 = tAtoms[startAtom].tree["children"][1];
        int iAt3 = tAtoms[startAtom].tree["children"][2];
        int iAt4 = tAtoms[startAtom].tree["children"][3];
        */
        
        int iBo1 = getBond(tBonds, startAtom, Ats[0]);
        if (iBo1 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[Ats[0]].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);
                   
        }
        
        int iBo2 = getBond(tBonds, startAtom, Ats[1]);
        if (iBo2 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[Ats[1]].id << ". Check what went wrong!"
                      << std::endl;
                exit(1);      
        }
            
        int iBo3 = getBond(tBonds, startAtom, Ats[2]);
        if (iBo3 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[Ats[2]].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);      
        }
        
        int iBo4 = getBond(tBonds, startAtom, Ats[3]);
        if (iBo4 < 0)
        {
            std::cout << "There is no bond between atom " 
                      << tAtoms[startAtom].id << " and atom  "
                      << tAtoms[Ats[3]].id << ". Check what went wrong!"
                      << std::endl;
            exit(1);      
        }
            
        
        int iAng14 = getAngle(tAngs, startAtom, Ats[0], Ats[3]);
        REAL tAng14 = tAngs[iAng14].value-90.0;
        int iAng24 = getAngle(tAngs, startAtom, Ats[1], Ats[3]);
        REAL tAng24 = tAngs[iAng24].value-90.0;
        int iAng34 = getAngle(tAngs, startAtom, Ats[2], Ats[3]);
        REAL tAng34 = tAngs[iAng34].value-90.0;
        
        REAL tC120 = cos(120.0*PI180);
        REAL tS120 = sin(120.0*PI180);
        tAtoms[Ats[0]].coords[0] =  tBonds[iBo1].value*cos(tAng14*PI180);
        tAtoms[Ats[0]].coords[1] =  0.0;
        tAtoms[Ats[0]].coords[2] =  tBonds[iBo1].value*sin(tAng14*PI180);
            
        tAtoms[Ats[1]].coords[0] =  tBonds[iBo2].value*cos(tAng24*PI180)*tC120;
        tAtoms[Ats[1]].coords[1] =  tBonds[iBo2].value*cos(tAng24*PI180)*tS120;
        tAtoms[Ats[1]].coords[2] =  tBonds[iBo2].value*sin(tAng24*PI180);
            
        tAtoms[Ats[2]].coords[0] =  tBonds[iBo3].value*cos(tAng34*PI180)*tC120;
        tAtoms[Ats[2]].coords[1] = -tBonds[iBo3].value*cos(tAng34*PI180)*tS120;
        tAtoms[Ats[2]].coords[2] =  tBonds[iBo3].value*sin(tAng34*PI180);
        
        tAtoms[Ats[3]].coords[0] =   0.0;
        tAtoms[Ats[3]].coords[1] =   0.0;
        tAtoms[Ats[3]].coords[2] =  -tBonds[iBo4].value;
        
    }
    
    void buildAtomTree::buildTree(std::vector<AtomDict>& tAtoms, 
                                  std::vector<BondDict>& tBonds, 
                                  std::vector<AngleDict>& tAngles, 
                                  std::vector<TorsionDict>& tTorsions, 
                                  std::vector<RingDict>& tRings, 
                                  std::vector<PlaneDict>& tPlas, 
                                  std::vector<ChiralDict>& tChs)
    {
        setAtomsMST(tAtoms, tBonds);
        
        setStartStruct(tAtoms,  tBonds, tAngles, tChs);
        
        setTreeAtomValues(tAtoms, tBonds, tAngles, tTorsions);
       
    }
    
    
}
