
/* 
 * File:   neighbList.cpp
 * Author: flong
 *
 * Created on September 28, 2011, 10:48 AM
 */

#include "neighbList.h"

namespace LIBMOL
{
    NBCell::NBCell()
    {
    }
    
    NBCell::NBCell(const NBCell& tC)
    {
        for (std::vector<int>::const_iterator iI=tC.index.begin();
                iI !=tC.index.end(); iI++)
        {
            index.push_back((*iI));
        }
        for (std::vector<Atom>::const_iterator iA = tC.atomsInCell.begin();
                iA !=tC.atomsInCell.end(); iA++)
        {
            atomsInCell.push_back((*iA));
        }
        for (std::vector<std::vector<int> >::const_iterator iNC = tC.nbcellList.begin();
                iNC !=tC.nbcellList.end(); iNC++)
        {
            std::vector<int> tIdx;
            for(std::vector<int>::const_iterator iIdx=iNC->begin();
                    iIdx != iNC->end(); iIdx++)
            {
                tIdx.push_back((*iIdx));
            }
            nbcellList.push_back(tIdx);
        }
        for (std::vector<ResidueID>::const_iterator iRe=tC.residueList.begin();
                iRe != tC.residueList.end(); iRe++)
        {
            residueList.push_back(*iRe);
        }
    }
    
    NeighbList::NeighbList():itsDim(3),
            itsCutOff(ZeroReal),
            itsNBShell(ZeroReal),
            itsErrInfo(NullString),
            itsErrLevel(0)
    {
    }
    
    NeighbList::NeighbList(const NeighbList & tNB):itsDim(tNB.getDim()),
            itsCutOff(tNB.getCutOff()),
            itsNBShell(tNB.getNBShell()),
            itsErrInfo(tNB.getErrInfo()),
            itsErrLevel(tNB.getErrLevel())
        
    {
        int i = tNB.getDim();
        if (i)
        {
            for (std::vector<REAL>::const_iterator iT=tNB.coordsMax.begin();
                    iT != tNB.coordsMax.end(); iT++)
            {
                coordsMax.push_back((*iT));
            } 
            for (std::vector<REAL>::const_iterator iT=tNB.coordsMin.begin();
                    iT != tNB.coordsMin.end(); iT++)
            {
                coordsMin.push_back((*iT));
            } 
        }
        
        for (std::vector<Atom>::const_iterator iA=tNB.allAtoms.begin();
                iA != tNB.allAtoms.end(); iA++)
        {
            allAtoms.push_back((*iA));
        }
        
        for (std::vector<NBCell>::const_iterator iNC = tNB.allCells.begin();
                iNC != tNB.allCells.end(); iNC++)
        {
            allCells.push_back((*iNC));
        }

        for (std::vector<NBCell>::const_iterator iNC = tNB.allNECells.begin();
                iNC != tNB.allNECells.end(); iNC++)
        {
            allNECells.push_back((*iNC));
        }
        
        for (std::map<int, std::list<Atom> >::const_iterator iAN =
                tNB.atomNBlist.begin(); 
                iAN != tNB.atomNBlist.end() ; iAN++)
        {
            atomNBlist.insert(*iAN);
        } 
        for(std::map<std::string, std::list<ResidueID> >::const_iterator iRe =
                tNB.residueNBList.begin(); 
                iRe != tNB.residueNBList.end(); iRe++)
        {
           residueNBList.insert(*iRe);
        }
    }
    
    NeighbList::NeighbList(std::vector<Atom> & tAtoms, 
            int tDim, REAL tNBC, REAL tNBS):itsDim(tDim),
            itsCutOff(tNBC),
            itsNBShell(tNBS),
            itsErrInfo(NullString),
            itsErrLevel(ZeroInt)
    {
        
        for(std::vector<Atom>::iterator iT=tAtoms.begin();
                iT != tAtoms.end(); iT++)
        {
            allAtoms.push_back((*iT));
        }
    }
    
    NeighbList::~NeighbList()
    {
    }
    
    std::string NeighbList::getErrInfo() const
    {
        return itsErrInfo;
    }
    void NeighbList::setErrInfo(std::string tE)
    {
        itsErrInfo.append(tE);
        itsErrInfo.append("\n");
    }
    
    int NeighbList::getErrLevel() const
    {
        return itsErrLevel;
    }
    void NeighbList::setErrLevel(int tN)
    {
        itsErrLevel = tN;
    }
    
    int  NeighbList::getDim() const
    {
        return itsDim;
    }
    void  NeighbList::setDim(int tD)
    {
        itsDim = tD;
    }
        
    REAL  NeighbList::getCutOff() const
    {
        return itsCutOff;
    }
    void NeighbList::setCutOff(REAL tC)
    {
        itsCutOff = tC;
    }
        
    REAL NeighbList::getNBShell() const
    {
        return itsNBShell; 
    }
    void NeighbList::setNBShell(REAL tC)
    {
        itsNBShell = tC;
    }
    
    void NeighbList::buildCellSystem()
    {
        // determine the range of atom distribute
        //std::cout << "Build Cell System " << std::endl;
        time_t rtime;
        time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl;
        
        if(allAtoms.size() && itsDim)
        {
            for(int i=0; i < itsDim; i++)
            {
               coordsMax.push_back(allAtoms[0].coords[i]);
               coordsMin.push_back(allAtoms[0].coords[i]);
            }
            
            for(std::vector<Atom>::iterator iT=allAtoms.begin();
                    iT != allAtoms.end(); iT++)
            {
                for (int i=0; i < itsDim; i++)
                {
                    if (iT->coords[i] > coordsMax[i])
                    {
                        coordsMax[i] = iT->coords[i];
                    }
                    if(iT->coords[i] < coordsMin[i])
                    {
                        coordsMin[i] = iT->coords[i];
                    }
                }
            }
            //std::cout << "Max coords " <<  coordsMax[0]
            //        << ", " << coordsMax[1] << ", "
            //        << coordsMax[2] << std::endl;
            //std::cout << "min coords " <<  coordsMin[0]
            //          << ", " << coordsMin[1] << ", "
            //          << coordsMin[2] << std::endl;      
            //
            // get cell length and number of cells in each dimension
            
            REAL tCL = getCutOff() + getNBShell();
            for(int i=0; i < itsDim; i++)
            {
                itsNumCell.push_back(1);
                
                while (coordsMin[i]+itsNumCell[i]*(tCL) < coordsMax[i])
                {
                    itsNumCell[i]+=1;
                }
                //std::cout << "number of cells in dim " << i+1 
                //          << " is " << itsNumCell[i] << std::endl;
                REAL tL = (coordsMax[i]-coordsMin[i])/itsNumCell[i];
                if (tL > 0.00000001)
                {
                    itsCellLength.push_back(tL);
                }
                else
                {
                    itsErrInfo.append("the Cell length is 0 in dimension %d",i);
                    itsErrLevel = 1;
                    throw (cellLengthException());
                }
                // std::cout << "cell length in dim " << i 
                //        << " is " << itsCellLength[i] << std::endl;
            }
           
            // Generate cells with neighbor cell lists, 
            // temporally 3-d
            // std::cout << "Biulding the cell system (including cell NBLists " 
            //        << std::endl;
            if (itsDim == 3)
            {
                // int iCout =0;
                for (int i0 =0; i0 < itsNumCell[0]; i0++)
                {
                    for (int i1=0; i1 < itsNumCell[1]; i1++)
                    {
                        for (int i2=0; i2 < itsNumCell[2]; i2++)
                        {
                            NBCell tCell;
                            tCell.index.push_back(i0);
                            tCell.index.push_back(i1);
                            tCell.index.push_back(i2);
                            allCells.push_back(tCell);
                            //std::cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << std::endl;
                            //std::cout << "Cell Position " << iCout << std::endl;
                            //std::cout << "confirm " << allCells.size()-1 << std::endl;
                            //int aIdx = (itsNumCell[1]*i0+i1)*itsNumCell[2] +i2; 
                            //std::cout << "From calc " << aIdx << std::endl; 
                            //iCout++;
                        }
                    }
                }
            } 
            
            //std::cout << "Number of cells " <<  allCells.size()   << std::endl;
            
            //std::cout << "Finishing cell building and now put atoms into cell" 
            //          << std::endl;
            //time(&rtime);
            //std::cout << "Current time is " << ctime(&rtime) << std::endl;            
            putAtomsAndResiduesInCell();
            //std::cout << "Atoms are in cell now " << std::endl;
            //time(&rtime);
            //std::cout << "Current time is " << ctime(&rtime) << std::endl;
            
            // build cell nblist for each cell and copy them to non-empty cells
            for (std::vector<NBCell>::iterator iNB =allCells.begin(); 
                iNB != allCells.end(); iNB++)
            {
                if (iNB->atomsInCell.size() >0)
                {
                    buildCellNBList(iNB);
                    // allNECells.push_back(*iNB);
                }
            } 
        }
       
        
        /*
        std::cout << "Non-empty cells " <<  allNECells.size() << std::endl;
       
        // Check nonempty cell;
        for (std::vector<NBCell>::iterator iNB =allNECells.begin(); 
                iNB != allNECells.end(); iNB++)
        {
             std::cout <<std::endl << "Cell(" << iNB->index[0] << ","
                << iNB->index[1] << "," <<iNB->index[2] << ")"
                <<std::endl;
             std::cout << "Number of atoms in the cell " 
                      << iNB->atomsInCell.size() << std::endl;
             std::cout << "Number of residues in the cell " 
                      << iNB->residueList.size() << std::endl;
        }
        
        */
    }
    
    void NeighbList::buildCellNBList(std::vector<NBCell>::iterator tNB)
    {
        int tI0, tI1, tI2;
        int tN0, tN1, tN2;
     
        for (tN0 =-1; tN0 < 2; tN0++)
        {
            tI0 = tNB->index[0] + tN0;
            if (tI0 >= 0 && tI0 < itsNumCell[0])
            {
                for (tN1 =-1; tN1 < 2; tN1++)
                {
                    tI1 = tNB->index[1]+tN1;
                    if (tI1 >= 0 && tI1< itsNumCell[1])
                    {
                        for (tN2 =-1; tN2 < 2; tN2++)
                        {
                            tI2 = tNB->index[2]+tN2;
                            if (tI2 >= 0 && tI2 < itsNumCell[2])
                            {
                               std::vector<int> tIdx;
                               tIdx.push_back(tI0);
                               tIdx.push_back(tI1);
                               tIdx.push_back(tI2);
                               if (!(tNB->index[0]==tI0 && tNB->index[1]==tI1
                                       &&tNB->index[2]==tI2))
                               {
                                   int aIdx =(itsNumCell[1]*tI0+tI1)*itsNumCell[2] 
                                              + tI2; 
                                   if(aIdx >=0)
                                   {
                                       if((int)allCells[aIdx].atomsInCell.size()>0)
                                       {
                                           tNB->nbcellList.push_back(tIdx);
                                       }
                                   }
                               }
                            }
                        }
                    }
                }
            }
        }
    }
    
    void NeighbList::putAtomsAndResiduesInCell()
    {   
        //std::cout << "Number of atoms is " << allAtoms.size() << std::endl;
        
        for (std::vector<Atom>::iterator iT=allAtoms.begin();
                iT != allAtoms.end(); iT++)
        {
            // std::cout << "Atom " << iT->getSeriNum() << std::endl;
            std::vector<int> tIdx;
            for (int i1 = 0; i1 < (int)iT->coords.size(); i1++)
            {
                REAL tDx = (iT->coords[i1]-coordsMin[i1])/itsCellLength[i1];
                int nDx = (int)tDx;
                if(tDx*itsCellLength[i1] == nDx*itsCellLength[i1] && nDx !=0)
                {
                    tIdx.push_back(nDx-1);
                }
                else if(nDx == itsNumCell[i1])
                {
                    tIdx.push_back(nDx-1);
                }
                else
                {
                    tIdx.push_back(nDx);
                }
                //   std::cout << "Dx = " << tDx << " in dim " << i1 << std::endl;
                //   std::cout << " int-Dx " << nDx << std::endl;
                
            } 
            
            
            //std::cout<< "Coords X: " << iT->coords[0] << std::endl;
            //std::cout<< "Coords Y: " << iT->coords[1] << std::endl;
            //std::cout<< "Coords Z: " << iT->coords[2] << std::endl;
            //std:: cout << "tIdx " << tIdx[0] << " " << tIdx[1] << " " << tIdx[2] 
            //           << std::endl;
            int aIdx = (itsNumCell[1]*tIdx[0]+tIdx[1])*itsNumCell[2] 
                        + tIdx[2];
            //std::cout << "Cell position is " << aIdx << std::endl;
            
            if(aIdx != -1)
            {
                allCells[aIdx].atomsInCell.push_back((*iT));
                
                // The atom has been added into the cell.
                // check if the residue is already in the residue list.
                // std::cout << "atom serial number " << iT->getSeriNum()
                //           << std::endl;
                // std::cout << "Cell index " << allCells[aIdx].index[0]
                //        << "," << allCells[aIdx].index[1]
                //        << "," << allCells[aIdx].index[2] << std::endl;
                ResidueID tResID;
                tResID.resName = TrimSpaces(iT->getResName());
                tResID.seqNum  = iT->getSeqNum();
                tResID.chainID = TrimSpaces(iT->getChainID());
                // std::cout << "R1: name  " << tResID.resName
                //        << " seqNum " <<  tResID.seqNum
                //        << " Chian ID " << tResID.chainID
                //        << std::endl;
                // std::cout << "number of residues in cell " << aIdx << " is "
                //          << allCells[aIdx].residueList.size()
                //          << std::endl;
                
                bool lFound = false;
                for (std::vector<ResidueID>::iterator iR =
                        allCells[aIdx].residueList.begin();
                        iR != allCells[aIdx].residueList.end(); iR++)
                {
                    //std::cout << "R2: name  " << iR->resName
                    //         << " seqNum " <<  iR->seqNum
                    //         << " Chian ID " << iR->chainID
                    //         << std::endl;
                    if(iR->resName == tResID.resName &&
                       iR->seqNum  == tResID.seqNum && 
                       iR->chainID == tResID.chainID)
                    {
                        // the residue is already in the residueNBList
                        lFound = true;
                        break;
                    }
                }
                
                
                if(not lFound)
                {
                    allCells[aIdx].residueList.push_back(tResID);
                    // std::cout << "Residue " << tResID.resName << "_"
                    //           << tResID.seqNum << "_" 
                    //           << tResID.chainID << " has been added to residue list "
                    //           << std::endl;
                }


               
            }
        }
        
    }
    
    int NeighbList::getOneNBCellPos(std::vector<int> & tIdx)
    {
        
        if(allCells.size() > 0)
        {
            bool lC; 
            for (int i=0; i < (int)allCells.size(); i++)
            {
                lC=true;
                // find cell index, allow 1, 2, 3-D system 
                for (int j=0; j < (int)allCells[i].index.size(); j++)
                {
                    
                    if( allCells[i].index[j] != tIdx[j])
                    {
                        lC = false;
                    }
                }
                if (lC)
                    {   
                        return i;
                    }               
            }     
        }
        
        return -1;
    }
    
    bool  NeighbList::getOneNBCell(std::vector<int> & tIdx)
    {
        bool lE= false;
        for (std::vector<NBCell>::iterator iC = allCells.begin();
                iC != allCells.end(); iC++)
        {
            bool lC = true;
            // find cell index, allow 1, 2, 3-D system 
            for (int i =0; i < (int)iC->index.size(); i++ )
            {
                if(iC->index[i] != tIdx[i])
                {
                    lC = false;
                    break;
                }
            }
            if (lC)
            {
                lE = true;
                break;
            }
        }
        return lE;
    }
    
    bool  NeighbList::getOneNENBCell(std::vector<int> & tIdx)
    {
        bool lE= false;
        for (std::vector<NBCell>::iterator iC = allNECells.begin();
                iC != allNECells.end(); iC++)
        {
            bool lC = true;
            // find cell index, allow 1, 2, 3-D system 
            for (int i =0; i < (int)iC->index.size(); i++ )
            {
                if(iC->index[i] != tIdx[i])
                {
                    lC = false;
                    break;
                }
            }
            if (lC)
            {
                lE = true;
                break;
            }
        }
        
        return lE;
    }
    
    void NeighbList::buildAtomNeighbList()
    {
        //std::cout << "Build nblist for each atom in the system " 
        //          << std::endl;
        //time_t rtime;
        //time(&rtime);
        
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        for(std::vector<NBCell>::iterator iC1=allCells.begin();
                iC1 !=allCells.end(); iC1++)
        {
            // check the same cell
            //std::cout << "Cell " << iC1->index[0] << " "
            //        << iC1->index[1] << " "
            //       << iC1->index[2] << std::endl;
            //std::cout << "Number of atoms in this cell is "
            //          << iC1->atomsInCell.size() 
            //          << std::endl;
            if ((int)iC1->atomsInCell.size() >0)
            {
                if ((int)iC1->atomsInCell.size() >1)
                {
                    //std::cout << "atoms from own cell " << std::endl;
                    buildAtomNeighListfromAPairCell((*iC1), (*iC1));
                }
                // check other cells 
                int i10, i20, i30;
                i10 = iC1->index[0];
                i20 = iC1->index[1];
                i30 = iC1->index[2];
                int aId0 = (itsNumCell[1]*i10+i20)*itsNumCell[2] + i30;
                for (std::vector<std::vector<int> >::iterator iI = 
                        iC1->nbcellList.begin();
                        iI !=iC1->nbcellList.end(); iI++)
                {   
                    int i1 = (*iI)[0];
                    int i2 = (*iI)[1];
                    int i3 = (*iI)[2];
                    int aId = (itsNumCell[1]*i1+i2)*itsNumCell[2] + i3;
                    NBCell tC2 = allCells[aId];
                    
                    //std::cout << "Number of atoms in this NB cell: "
                    //        << tC2.atomsInCell.size() << std::endl;
                    if ((int)tC2.atomsInCell.size() >0 and aId0 < aId )
                    {
                        buildAtomNeighListfromAPairCell((*iC1), tC2);
                    }
                }
            }
        }
        
       //std::cout << "Finish atom nblist building " << std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
    }
    // Build atom neighbList from a pair of cell
    void NeighbList::buildAtomNeighListfromAPairCell(NBCell & tCell1, 
                                                     NBCell & tCell2)
    {
        for (std::vector<Atom>::iterator iA1=tCell1.atomsInCell.begin();
                    iA1!=tCell1.atomsInCell.end(); iA1++)
        {
            
            //std::cout << "Build NBList for atom " << iA1->getSeriNum()
            //        << std::endl;
            for (std::vector<Atom>::iterator iA2 = tCell2.atomsInCell.begin();
                        iA2 != tCell2.atomsInCell.end(); iA2++)
            {
                REAL tD = distanceV(iA1->coords, iA2->coords);
                    if ( tD <=(itsCutOff+itsNBShell) && tD >= 0.0005)
                    {
                        iA1->neighbAtoms.push_back((*iA2));
                        iA2->neighbAtoms.push_back((*iA1));
                        //std::cout << "atom " << iA2->getSeriNum()
                        //          << " has been added to the list"
                        //          << std::endl;
                    }
            }
        }
    }
    /*
    void NeighbList::buildAtomNeighbList()
    {
        std::cout << "Build nblist for each atom in the system " 
                << std::endl;
        time_t rtime;
        time(&rtime);
        std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        for(std::vector<NBCell>::iterator iN=allCells.begin();
                iN !=allCells.end(); iN++)
        {
            for (std::vector<Atom>::iterator iA1=iN->atomsInCell.begin();
                    iA1!=iN->atomsInCell.end(); iA1++)
            {
                //std::cout << "NB atoms for atom " << iA1->getSeriNum()
                //          << " : " << std::endl;
                //time(&rtime);
                //std::cout << "Current time is " << ctime(&rtime) << std::endl;
                
                // check if atoms in the same cell satisfy the criteria 
                //std::cout << "NB atoms from the same cell " << std::endl;
               
                for (std::vector<Atom>::iterator iA2 = iN->atomsInCell.begin();
                        iA2 != iN->atomsInCell.end(); iA2++)
                {
                    REAL tD = distanceV(iA1->coords, iA2->coords);
                    if ( tD <=(itsCutOff+itsNBShell) && tD >= 0.0001)
                    {
                        iA1->neighbAtoms.push_back((*iA2));
                        //std::cout << "atom " << iA2->getSeriNum() 
                        //        << std::endl;
                    }
                }
               
                
                // check if atoms in the neighbor cells satisfy the criteria
                //std::cout << "NB atoms from other cell " << std::endl;
                
                for (int iN2 = 0; iN2 < (int)iN->nbcellList.size(); iN2++)
                {
                    
                    bool lExist = getOneNBCell(iN->nbcellList[iN2]);
                    if(lExist)
                    { 
                    for(std::vector<Atom>::iterator iA3=
                            iN->nbcellList[iN2].atomsInCell.begin();
                            iA3 !=iN->nbcellList[iN2].atomsInCell.end(); iA3++)
                    {
                        REAL tD = distanceV(iA1->coords, iA3->coords);
                        if(tD <=(itsCutOff+itsNBShell))
                        {
                            iA1->neighbAtoms.push_back((*iA3));
                            //std::cout << "atom " << iA3->getSeriNum() 
                            //    << std::endl;
                        }
                    }
                    }
                }
                // time(&rtime);
                // std::cout << "Current time is " << ctime(&rtime) << std::endl;                 
                if (iA1->neighbAtoms.size() !=0)
                {
                    atomNBlist.insert(std::pair<int,std::list<Atom> >
                    (iA1->getSeriNum(), iA1->neighbAtoms));
                }
                //if(iA1->getSeriNum()==34)
                //{
                //std::cout << "Atom " << iA1->getSeriNum() << " has "
                //        << iA1->neighbAtoms.size() << " neighbor atoms."
                //        << std::endl;
                
                //for (std::list<Atom>::iterator iANA=iA1->neighbAtoms.begin();
                //        iANA != iA1->neighbAtoms.end(); iANA++)
                //{
                //    std::cout << "Atom " << iANA->getSeriNum() << std::endl;
                //}
                //}
                 
                //std::cout << "update the residue list" << std::endl;
                updateResidueNBList(iA1);        
                //std::cout << "finishing update the residue list" << std::endl;
                //time(&rtime);
                //std::cout << "Current time is " << ctime(&rtime) << std::endl;
            }
        }
        
        // check residue neighbor list 
        for (std::map<std::string, std::list<ResidueID> >::iterator iM = residueNBList.begin();
                iM !=residueNBList.end(); iM++)
        {
            std::cout << std::endl;
            std::cout << "Residue " << iM->first << " NB list: "
                    << std::endl;
            for (std::list<ResidueID>::iterator iRID=iM->second.begin();
                   iRID != iM->second.end(); iRID++)
            {
                std::cout << "NB residue " << iRID->resName << " SeqNum "
                        << iRID->seqNum << " Chain ID " << iRID->chainID
                        << std::endl;
            }
        } 
        
        
        std::cout << "Finish atom nblist building " << std::endl;
        time(&rtime);
        std::cout << "Current time is " << ctime(&rtime) << std::endl; 
    }
    */
    
    void NeighbList::buildResidueNBList()
    {
        //std::cout << "Build nblist for each residue in the system " 
        //          << std::endl;
        // time_t rtime;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        for(std::vector<NBCell>::iterator iC1=allCells.begin();
                iC1 !=allCells.end(); iC1++)
        {
            //std::cout << "Cell " << iC1->index[0] << ","
            //          << iC1->index[1] << ","
            //          << iC1->index[2] << std::endl;
            //std::cout << "Number of residue the cell " 
            //          << iC1->residueList.size()  << std::endl;
            
            if((int)iC1->residueList.size() !=0)
            {
                for (std::vector<std::vector<int> >::iterator iI = 
                        iC1->nbcellList.begin();
                        iI !=iC1->nbcellList.end(); iI++)
                {
                    std::vector<int> tIdx = (*iI);
                    int aIdx = (itsNumCell[1]*tIdx[0]+tIdx[1])*itsNumCell[2]
                                + tIdx[2];
                    NBCell tC2 = allCells[aIdx];
                    //std::cout << "Its NBCell " << tC2.index[0] << ","
                    //          << tC2.index[1] << ","
                    //          << tC2.index[2] << std::endl;
                    //std::cout << "Number of residues in this NBCell "
                    //          << tC2.residueList.size() << std::endl;
                    
                    if ((int)tC2.residueList.size() !=0)
                    {
                        buildResidueNeighbListfromAPairCell((*iC1), tC2);
                    }        
                }
            }
        }
        //std::cout << "Building residue NBlist finished " << std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        /*
        for (std::map<std::string, std::list<ResidueID> >::iterator iM=
                residueNBList.begin();
                iM !=residueNBList.end(); iM++ )
        {
            std::cout << "\nNB List for Residue  " << iM->first << " : "
                    << std::endl;
            for (std::list<ResidueID>::iterator iR = iM->second.begin();
                    iR != iM->second.end(); iR++)
            {
                std::cout << "Residue " << iR->resName 
                          << " seqNum " << iR->seqNum 
                          << " Chain ID " << iR->chainID << std::endl;
            }
        } 
         */
    }   
    
    void NeighbList::buildResidueNeighbListfromAPairCell(NBCell & tCell1,
                                                         NBCell & tCell2)
    {
        for (std::vector<ResidueID>::iterator iRe1 =
                tCell1.residueList.begin();
                iRe1 != tCell1.residueList.end(); iRe1++)
        {
            std::string tID = iRe1->resName + "_"
                      +TrimSpaces(IntToStr(iRe1->seqNum)) + "_"
                      +iRe1->chainID;
            for(std::vector<ResidueID>::iterator iRe2 =
                    tCell2.residueList.begin(); 
                    iRe2 != tCell2.residueList.end(); iRe2++)
            {
                bool lFound = false;

                for (std::list<ResidueID>::iterator iR = 
                        residueNBList[tID].begin();
                        iR != residueNBList[tID].end(); iR++)
                {
                    if(iR->resName == iRe2->resName &&
                        iR->seqNum == iRe2->seqNum && 
                        iR->chainID == iRe2->chainID)
                    {
                        // the residue the neighbor atom belongs to 
                        // is already in the residueNBList
                        lFound = true;
                        break;
                    }
                }
                if (not lFound)
                {
                    residueNBList[tID].push_back(*iRe2);
                }
            }
                    
        }
    }
    
    void NeighbList::updateResidueNBList(std::vector<Atom>::iterator tA)
    {
        ResidueID tResID;
        tResID.resName = TrimSpaces(tA->getResName());
        tResID.seqNum  = tA->getSeqNum();
        tResID.chainID = TrimSpaces(tA->getChainID());
        
        std::string tID = tResID.resName + "_"
                          +TrimSpaces(IntToStr(tResID.seqNum)) + "_"
                          +tResID.chainID;
        
       
        //if(tA->getSeriNum()==56)
        //{
        //    std::cout << "Update Residue NB list for " << tID << std::endl;
        //}
        for (std::list<Atom>::iterator iA1 = tA->neighbAtoms.begin();
                    iA1 !=tA->neighbAtoms.end(); iA1++)
        {
            bool lFound = false;
            ResidueID tResID2;
            tResID2.resName = TrimSpaces(iA1->getResName());
            tResID2.seqNum  = iA1->getSeqNum();
            tResID2.chainID = TrimSpaces(iA1->getChainID());
            //if(tA->getSeriNum()==56)
            //{
            //    std::cout << "check atom " << iA1->getSeriNum() << std::endl;
            //}
            for (std::list<ResidueID>::iterator iR = residueNBList[tID].begin();
                    iR != residueNBList[tID].end(); iR++) 
            {
                if(iR->resName == tResID2.resName &&
                        iR->seqNum == tResID2.seqNum && 
                        iR->chainID == tResID2.chainID)
                {
                    // the residue the neighbor atom belongs to 
                    // is already in the residueNBList
                    lFound = true;
                    break;
                }     
            }
            if (not lFound)
            {
                std::string tID2 = tResID2.resName + "_"
                          +TrimSpaces(IntToStr(tResID2.seqNum)) + "_"
                          +tResID2.chainID;
                if(tID2.compare(tID) != 0)
                {
                    // The residue is not in the residueNBList, add it in 
                    residueNBList[tID].push_back(tResID2);
                    //if(tA->getSeriNum()==56)
                    //{
                    //std::cout << "Residue Pair:  " << tID 
                    //          << "-------"<< tID2 << "has been added to "
                    //          << " Pair-residue list " << std::endl;
                    //}
                }
            }
        }
    }
    
    void NeighbList::building(int tMode)
    {
        /*
         tMode = 0, default Mode. Build a atom neighbList
         tMode = 1, optional Mode. Build a Residue neighbList directly from 
                    the cell system without atom neighList 
         tMode = 2, optional Mode. Build both atom and residue neighbLists.
                    The latter is based on the former.  
         select a mode depending on the situation. 
          
       */
        
        //std::cout << std::endl << "------------------------------------------"
        //          << std::endl;
        //std::cout << "Build neighbor list System " << std::endl;
        //time_t rtime;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);
        
        buildCellSystem();
        if (!itsErrLevel)
        {
            if (tMode ==0)
            {
                buildAtomNeighbList();
            }
            else if (tMode == 1)
            {
                buildResidueNBList();
            }
        }
        //std::cout << "Finish building neighbor lists "<< std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);         
        //std::cout << "------------------------------------------"
        //          << std::endl;
        
    }
    
    
    // Another version for AtomDict

    NBCellDict::NBCellDict()
    {
    }
    
    NBCellDict::NBCellDict(const NBCellDict& tC)
    {
        for (std::vector<int>::const_iterator iI=tC.index.begin();
                iI !=tC.index.end(); iI++)
        {
            index.push_back((*iI));
        }
        for (std::vector<int>::const_iterator iA = tC.atomsInCell.begin();
                iA !=tC.atomsInCell.end(); iA++)
        {
            atomsInCell.push_back((*iA));
        }
        for (std::vector<std::vector<int> >::const_iterator iNC = tC.nbcellList.begin();
                iNC !=tC.nbcellList.end(); iNC++)
        {
            std::vector<int> tIdx;
            for(std::vector<int>::const_iterator iIdx=iNC->begin();
                    iIdx != iNC->end(); iIdx++)
            {
                tIdx.push_back((*iIdx));
            }
            nbcellList.push_back(tIdx);
        }
        for (std::vector<ResidueID>::const_iterator iRe=tC.residueList.begin();
                iRe != tC.residueList.end(); iRe++)
        {
            residueList.push_back(*iRe);
        }
    }
    
    NeighbListDict::NeighbListDict():itsErrInfo(NullString),
                                     itsErrLevel(ZeroInt)
    {
    }
    
 
    NeighbListDict::~NeighbListDict()
    {
    }
    
    
    void NeighbListDict::buildCellSystem(std::vector<AtomDict> & aAtomList,
                                    int tDim, REAL tNBCutoff, REAL tNBShell)
    {
        // determine the range of atom distribute
        //std::cout << "Build Cell System " << std::endl;
        time_t rtime;
        time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl;
        
        if(aAtomList.size() && tDim)
        {
            for(int i=0; i < tDim; i++)
            {
               coordsMax.push_back(aAtomList[0].coords[i]);
               coordsMin.push_back(aAtomList[0].coords[i]);
            }
            
            for(std::vector<AtomDict>::iterator iT=aAtomList.begin();
                    iT != aAtomList.end(); iT++)
            {
                for (int i=0; i < tDim; i++)
                {
                    if (iT->coords[i] > coordsMax[i])
                    {
                        coordsMax[i] = iT->coords[i];
                    }
                    if(iT->coords[i] < coordsMin[i])
                    {
                        coordsMin[i] = iT->coords[i];
                    }
                }
            }
            
            std::cout << "Max coords " <<  coordsMax[0]
                      << ", " << coordsMax[1] << ", "
                      << coordsMax[2] << std::endl;
            std::cout << "min coords " <<  coordsMin[0]
                      << ", " << coordsMin[1] << ", "
                      << coordsMin[2] << std::endl;      
            
            // get cell length and number of cells in each dimension
            
            REAL tCL = tNBCutoff + tNBShell;
            for(int i=0; i < tDim; i++)
            {
                itsNumCell.push_back(1);
                
                while (coordsMin[i]+itsNumCell[i]*(tCL) < coordsMax[i])
                {
                    itsNumCell[i]+=1;
                }
                //std::cout << "number of cells in dim " << i+1 
                //          << " is " << itsNumCell[i] << std::endl;
                REAL tL = (coordsMax[i]-coordsMin[i])/itsNumCell[i];
                if (tL > 0.00000001)
                {
                    itsCellLength.push_back(tL);
                }
                else
                {
                    itsCellLength.push_back(0.01);
                    //itsErrInfo.append("the Cell length is 0 in dimension %d",i);
                    //itsErrLevel = 1;
                    //throw (cellLengthException());
                }
                // std::cout << "cell length in dim " << i 
                //        << " is " << itsCellLength[i] << std::endl;
            }
           
            // Generate cells with neighbor cell lists, 
            // temporally 3-d
            // std::cout << "Biulding the cell system (including cell NBLists " 
            //        << std::endl;
            if (tDim == 3)
            {
                //int iCout =0;
                for (int i0 =0; i0 < itsNumCell[0]; i0++)
                {
                    for (int i1=0; i1 < itsNumCell[1]; i1++)
                    {
                        for (int i2=0; i2 < itsNumCell[2]; i2++)
                        {
                            NBCellDict tCell;
                            tCell.index.push_back(i0);
                            tCell.index.push_back(i1);
                            tCell.index.push_back(i2);
                            allCells.push_back(tCell);
                            //std::cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << std::endl;
                            //std::cout << "Cell Position " << iCout << std::endl;
                            //std::cout << "confirm " << allCells.size()-1 << std::endl;
                            //int aIdx = (itsNumCell[1]*i0+i1)*itsNumCell[2] +i2; 
                            //std::cout << "From calc " << aIdx << std::endl; 
                            //iCout++;
                        }
                    }
                }
            } 
            
            //std::cout << "Number of cells " <<  allCells.size()   << std::endl;
            //std::cout << "Finishing cell building and now put atoms into cell" 
            //          << std::endl;
            //time(&rtime);
            //std::cout << "Current time is " << ctime(&rtime) << std::endl;            
            putAtomsAndResiduesInCell(aAtomList);
            //std::cout << "Atoms are in cell now " << std::endl;
            //time(&rtime);
            //std::cout << "Current time is " << ctime(&rtime) << std::endl;
            
            // build cell nblist for each cell and copy them to non-empty cells
            for (std::vector<NBCellDict>::iterator iNB =allCells.begin(); 
                iNB != allCells.end(); iNB++)
            {
                if (iNB->atomsInCell.size() >0)
                {
                    buildCellNBList(iNB);
                    allNECells.push_back(*iNB);
                }
            } 
        }
       
        
        
        //std::cout << "Non-empty cells " <<  allNECells.size() << std::endl;
        /*
        // Check nonempty cell;
        for (std::vector<NBCellDict>::iterator iNB =allNECells.begin(); 
                iNB != allNECells.end(); iNB++)
        {
             std::cout <<std::endl << "Cell(" << iNB->index[0] << ","
                << iNB->index[1] << "," <<iNB->index[2] << ")"
                <<std::endl;
             std::cout << "Number of atoms in the cell " 
                      << iNB->atomsInCell.size() << std::endl;
             //std::cout << "Number of residues in the cell " 
             //         << iNB->residueList.size() << std::endl;
        }
         */
        
        
       
    }
    
    void NeighbListDict::buildCellNBList(std::vector<NBCellDict>::iterator tNB)
    {
        int tI0, tI1, tI2;
        int tN0, tN1, tN2;
     
        for (tN0 =-1; tN0 < 2; tN0++)
        {
            tI0 = tNB->index[0] + tN0;
            if (tI0 >= 0 && tI0 < itsNumCell[0])
            {
                for (tN1 =-1; tN1 < 2; tN1++)
                {
                    tI1 = tNB->index[1]+tN1;
                    if (tI1 >= 0 && tI1< itsNumCell[1])
                    {
                        for (tN2 =-1; tN2 < 2; tN2++)
                        {
                            tI2 = tNB->index[2]+tN2;
                            if (tI2 >= 0 && tI2 < itsNumCell[2])
                            {
                               std::vector<int> tIdx;
                               tIdx.push_back(tI0);
                               tIdx.push_back(tI1);
                               tIdx.push_back(tI2);
                               if (!(tNB->index[0]==tI0 && tNB->index[1]==tI1
                                       &&tNB->index[2]==tI2))
                               {
                                   int aIdx =(itsNumCell[1]*tI0+tI1)*itsNumCell[2] 
                                              + tI2; 
                                   if(aIdx >=0)
                                   {
                                       if((int)allCells[aIdx].atomsInCell.size()>0)
                                       {
                                           tNB->nbcellList.push_back(tIdx);
                                       }
                                   }
                               }
                            }
                        }
                    }
                }
            }
        }
    }
    
    void NeighbListDict::putAtomsAndResiduesInCell(std::vector<AtomDict> & aAtomList)
    {   
        //std::cout << "Number of atoms is " << allAtoms.size() << std::endl;
        
        for (std::vector<AtomDict>::iterator iT=aAtomList.begin();
                iT != aAtomList.end(); iT++)
        {
            // std::cout << "Atom " << iT->id << std::endl;
            std::vector<int> tIdx;
            for (int i1 = 0; i1 < (int)iT->coords.size(); i1++)
            {
                REAL tDx = (iT->coords[i1]-coordsMin[i1])/itsCellLength[i1];
                int nDx = (int)tDx;
                if(tDx*itsCellLength[i1] == nDx*itsCellLength[i1] && nDx !=0)
                {
                    tIdx.push_back(nDx-1);
                }
                else if(nDx == itsNumCell[i1])
                {
                    tIdx.push_back(nDx-1);
                }
                else
                {
                    tIdx.push_back(nDx);
                }
                //   std::cout << "Dx = " << tDx << " in dim " << i1 << std::endl;
                //   std::cout << " int-Dx " << nDx << std::endl;
                
            } 
            
            
            //std::cout<< "Coords X: " << iT->coords[0] << std::endl;
            //std::cout<< "Coords Y: " << iT->coords[1] << std::endl;
            //std::cout<< "Coords Z: " << iT->coords[2] << std::endl;
            //std:: cout << "tIdx " << tIdx[0] << " " << tIdx[1] << " " << tIdx[2] 
            //           << std::endl;
            int aIdx = (itsNumCell[1]*tIdx[0]+tIdx[1])*itsNumCell[2] 
                        + tIdx[2];
            // std::cout << "Cell position is " << aIdx << std::endl;
            
            if(aIdx != -1)
            {
                allCells[aIdx].atomsInCell.push_back(iT->seriNum);
                
                // The atom has been added into the cell.
                // check if the residue is already in the residue list.
                // std::cout << "atom serial number " << iT->getSeriNum()
                //           << std::endl;
                // std::cout << "Cell index " << allCells[aIdx].index[0]
                //        << "," << allCells[aIdx].index[1]
                //        << "," << allCells[aIdx].index[2] << std::endl;
                // ResidueID tResID;
                // tResID.resName = TrimSpaces(iT->getResName());
                // tResID.seqNum  = iT->getSeqNum();
                // tResID.chainID = TrimSpaces(iT->getChainID());
                // std::cout << "R1: name  " << tResID.resName
                //        << " seqNum " <<  tResID.seqNum
                //        << " Chian ID " << tResID.chainID
                //        << std::endl;
                // std::cout << "number of residues in cell " << aIdx << " is "
                //          << allCells[aIdx].residueList.size()
                //          << std::endl;
                
                //bool lFound = false;
                //for (std::vector<ResidueID>::iterator iR =
                //        allCells[aIdx].residueList.begin();
                //        iR != allCells[aIdx].residueList.end(); iR++)
                //{
                    //std::cout << "R2: name  " << iR->resName
                    //         << " seqNum " <<  iR->seqNum
                    //         << " Chian ID " << iR->chainID
                    //         << std::endl;
                    //if(iR->resName == tResID.resName &&
                    //   iR->seqNum  == tResID.seqNum && 
                    //   iR->chainID == tResID.chainID)
                    //{
                        // the residue is already in the residueNBList
                        // lFound = true;
                        // break;
                    // }
                //}
                
                
                //if(not lFound)
                //{
                    //allCells[aIdx].residueList.push_back(tResID);
                    // std::cout << "Residue " << tResID.resName << "_"
                    //           << tResID.seqNum << "_" 
                    //           << tResID.chainID << " has been added to residue list "
                    //           << std::endl;
                //}
               
            }
        }
        
    }
    
 
 
    
    
    void NeighbListDict::buildAtomNeighbList(std::vector<AtomDict> & aAtomList,REAL tL)
    {
        //std::cout << "Build nblist for each atom in the system " 
        //          << std::endl;
        //time_t rtime;
        //time(&rtime);
        
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        for(std::vector<NBCellDict>::iterator iC1=allCells.begin();
                iC1 !=allCells.end(); iC1++)
        {
            // check the same cell
            // std::cout << "Cell " << iC1->index[0] << " "
            //          << iC1->index[1] << " "
            //          << iC1->index[2] << std::endl;
            // std::cout << "Number of atoms in this cell is "
            //          << iC1->atomsInCell.size() 
            //          << std::endl;
            if ((int)iC1->atomsInCell.size() >0)
            {
                if ((int)iC1->atomsInCell.size() >1)
                {
                    // std::cout << "atoms from own cell " << std::endl;
                    buildAtomNeighListfromAPairCell(aAtomList, tL, (*iC1), (*iC1));
                }
                // check other cells 
                int i10, i20, i30;
                i10 = iC1->index[0];
                i20 = iC1->index[1];
                i30 = iC1->index[2];
                int aId0 = (itsNumCell[1]*i10+i20)*itsNumCell[2] + i30;
                for (std::vector<std::vector<int> >::iterator iI = 
                        iC1->nbcellList.begin();
                        iI !=iC1->nbcellList.end(); iI++)
                {   
                    int i1 = (*iI)[0];
                    int i2 = (*iI)[1];
                    int i3 = (*iI)[2];
                    int aId = (itsNumCell[1]*i1+i2)*itsNumCell[2] + i3;
                    NBCellDict tC2 = allCells[aId];
                    // std::cout << "Cell " << tC2.index[0] << " "
                    //          << tC2.index[1] << " "
                    //          << tC2.index[2] << std::endl;
                    //std::cout << "Number of atoms in this NB cell: "
                    //          << tC2.atomsInCell.size() << std::endl;
                    if ((int)tC2.atomsInCell.size() >0 and aId0 < aId )
                    {
                        // std::cout << "atoms from inter-cells " << std::endl;
                        buildAtomNeighListfromAPairCell(aAtomList, tL, (*iC1), tC2);
                    }
                }
            }
        }
        
        
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
    }
    // Build atom neighbList from a pair of cell
    void NeighbListDict::buildAtomNeighListfromAPairCell(std::vector<AtomDict> & aAtomList,
                                                         REAL tL,
                                                         NBCellDict & tCell1, 
                                                         NBCellDict & tCell2)
    {
        
        // std::cout << (int) tCell1.atomsInCell.size();
        
        for (std::vector<int>::iterator iA1=tCell1.atomsInCell.begin();
                    iA1!=tCell1.atomsInCell.end(); iA1++)
        {
            
            // std::cout << "Build NBList for atom " << *iA1 << std::endl;
            for (std::vector<int>::iterator iA2 = tCell2.atomsInCell.begin();
                        iA2 != tCell2.atomsInCell.end(); iA2++)
            {
                REAL tD = distanceV(aAtomList[*iA1].coords, aAtomList[*iA2].coords);
                if ( tD <=tL && tD >= 0.0005)
                {
                    if (std::find(aAtomList[*iA1].neighbAtoms.begin(),
                        aAtomList[*iA1].neighbAtoms.end(), *iA2) ==aAtomList[*iA1].neighbAtoms.end())
                    {
                        aAtomList[*iA1].neighbAtoms.push_back((*iA2));
                        // std::cout << "atom " << aAtomList[*iA2].id
                        //          << " has been added to the list" << std::endl;
                    }
                    if (std::find(aAtomList[*iA2].neighbAtoms.begin(),
                        aAtomList[*iA2].neighbAtoms.end(), *iA1) ==aAtomList[*iA2].neighbAtoms.end())
                    {
                        aAtomList[*iA2].neighbAtoms.push_back((*iA1));  
                    }
                }
            }
        }
    }
    
    void NeighbListDict::buildResidueNBList()
    {
        //std::cout << "Build nblist for each residue in the system " 
        //          << std::endl;
        // time_t rtime;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        for(std::vector<NBCellDict>::iterator iC1=allCells.begin();
                iC1 !=allCells.end(); iC1++)
        {
            //std::cout << "Cell " << iC1->index[0] << ","
            //          << iC1->index[1] << ","
            //          << iC1->index[2] << std::endl;
            //std::cout << "Number of residue the cell " 
            //          << iC1->residueList.size()  << std::endl;
            
            if((int)iC1->residueList.size() !=0)
            {
                for (std::vector<std::vector<int> >::iterator iI = 
                        iC1->nbcellList.begin();
                        iI !=iC1->nbcellList.end(); iI++)
                {
                    std::vector<int> tIdx = (*iI);
                    int aIdx = (itsNumCell[1]*tIdx[0]+tIdx[1])*itsNumCell[2]
                                + tIdx[2];
                    NBCellDict tC2 = allCells[aIdx];
                    //std::cout << "Its NBCell " << tC2.index[0] << ","
                    //          << tC2.index[1] << ","
                    //          << tC2.index[2] << std::endl;
                    //std::cout << "Number of residues in this NBCell "
                    //          << tC2.residueList.size() << std::endl;
                    
                    if ((int)tC2.residueList.size() !=0)
                    {
                        buildResidueNeighbListfromAPairCell((*iC1), tC2);
                    }        
                }
            }
        }
        //std::cout << "Building residue NBlist finished " << std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime) << std::endl; 
        /*
        for (std::map<std::string, std::list<ResidueID> >::iterator iM=
                residueNBList.begin();
                iM !=residueNBList.end(); iM++ )
        {
            std::cout << "\nNB List for Residue  " << iM->first << " : "
                    << std::endl;
            for (std::list<ResidueID>::iterator iR = iM->second.begin();
                    iR != iM->second.end(); iR++)
            {
                std::cout << "Residue " << iR->resName 
                          << " seqNum " << iR->seqNum 
                          << " Chain ID " << iR->chainID << std::endl;
            }
        } 
         */
    }   
    
    void NeighbListDict::buildResidueNeighbListfromAPairCell(NBCellDict & tCell1,
                                                         NBCellDict & tCell2)
    {
        for (std::vector<ResidueID>::iterator iRe1 =
                tCell1.residueList.begin();
                iRe1 != tCell1.residueList.end(); iRe1++)
        {
            std::string tID = iRe1->resName + "_"
                      +TrimSpaces(IntToStr(iRe1->seqNum)) + "_"
                      +iRe1->chainID;
            for(std::vector<ResidueID>::iterator iRe2 =
                    tCell2.residueList.begin(); 
                    iRe2 != tCell2.residueList.end(); iRe2++)
            {
                bool lFound = false;

                for (std::list<ResidueID>::iterator iR = 
                        residueNBList[tID].begin();
                        iR != residueNBList[tID].end(); iR++)
                {
                    if(iR->resName == iRe2->resName &&
                        iR->seqNum == iRe2->seqNum && 
                        iR->chainID == iRe2->chainID)
                    {
                        // the residue the neighbor atom belongs to 
                        // is already in the residueNBList
                        lFound = true;
                        break;
                    }
                }
                if (not lFound)
                {
                    residueNBList[tID].push_back(*iRe2);
                }
            }
                    
        }
    }
    
 
    void NeighbListDict::building(std::vector<AtomDict> & aAtomList, 
                                  int tDim, REAL tNBCutoff, 
                                  REAL tNBShell, int tMode)
    {
        /*
         tMode = 0, default Mode. Build a atom neighbList
         tMode = 1, optional Mode. Build a Residue neighbList directly from 
                    the cell system without atom neighList 
         tMode = 2, optional Mode. Build both atom and residue neighbLists.
                    The latter is based on the former.  
         select a mode depending on the situation. 
          
       */
        
        //std::cout << std::endl << "------------------------------------------"
        //          << std::endl;
        // std::cout << "Build neighbor list System " << std::endl;
        //time_t rtime;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);
        
        for (std::vector<AtomDict>::iterator iA=aAtomList.begin();
                iA != aAtomList.end(); iA++)
        {
            if ((int)iA->neighbAtoms.size() !=0)
            {
                iA->neighbAtoms.clear();
            }
        }
        
        buildCellSystem(aAtomList, tDim, tNBCutoff, tNBShell);
        
        REAL fullL = tNBCutoff + tNBShell;
        // std::cout << "error level " << itsErrLevel << std::endl;
        if (!itsErrLevel)
        {
            if (tMode ==0)
            {
                // std::cout << "Build atom NB List " << std::endl;
                buildAtomNeighbList(aAtomList, fullL);
            }
            else if (tMode==1)
            {
                buildResidueNBList();
            }
        }
        //Check
        /*
        for (std::vector<AtomDict>::iterator iA=aAtomList.begin();
                iA != aAtomList.end(); iA++)
        {
            std::cout << "For atom " << iA->id;
            if ((int)iA->neighbAtoms.size() >0)
            {
                std::cout << ", its NB atoms atoms are:" << std::endl;
                for (std::vector<int>::iterator iNB=iA->neighbAtoms.begin();
                        iNB != iA->neighbAtoms.end(); iNB++)
                {
                    std::cout << aAtomList[*iNB].id << std::endl;
                }
            }
            else
            {
                std::cout << ", it has no NB atom " << std::endl;
            }
            
        }
        //std::cout << "Finish building neighbor lists "<< std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);         
        //std::cout << "------------------------------------------"
        //          << std::endl;
       */
    }
    
    void NeighbListDict::building(std::vector<AtomDict>   & aAtomList, 
                                  std::vector<AtomDict>   & tAllAtoms,
                                  int tDim, REAL tNBCutoff, 
                                  REAL tNBShell, int tMode)
    {
        /*
         tMode = 0, default Mode. Build a atom neighbList
         tMode = 1, optional Mode. Build a Residue neighbList directly from 
                    the cell system without atom neighList 
         tMode = 2, optional Mode. Build both atom and residue neighbLists.
                    The latter is based on the former.  
         select a mode depending on the situation. 
          
       */
        
        //std::cout << std::endl << "------------------------------------------"
        //          << std::endl;
        // std::cout << "Build neighbor list System " << std::endl;
        //time_t rtime;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);
        
        for (std::vector<AtomDict>::iterator iA=aAtomList.begin();
                iA != aAtomList.end(); iA++)
        {
            if ((int)iA->neighbAtoms.size() !=0)
            {
                iA->neighbAtoms.clear();
            }
        }
        
        buildCellSystem(aAtomList, tDim, tNBCutoff, tNBShell);
        
        REAL fullL = tNBCutoff + tNBShell;
        // std::cout << "error level " << itsErrLevel << std::endl;
        if (!itsErrLevel)
        {
            if (tMode ==0)
            {
                // std::cout << "Build atom NB List " << std::endl;
                buildAtomNeighbList(tAllAtoms, fullL);
            }
            
        }
        //Check
        /*
        for (std::vector<AtomDict>::iterator iA=aAtomList.begin();
                iA != aAtomList.end(); iA++)
        {
            std::cout << "For atom " << iA->id;
            if ((int)iA->neighbAtoms.size() >0)
            {
                std::cout << ", its NB atoms atoms are:" << std::endl;
                for (std::vector<int>::iterator iNB=iA->neighbAtoms.begin();
                        iNB != iA->neighbAtoms.end(); iNB++)
                {
                    std::cout << aAtomList[*iNB].id << std::endl;
                }
            }
            else
            {
                std::cout << ", it has no NB atom " << std::endl;
            }
            
        }
        //std::cout << "Finish building neighbor lists "<< std::endl;
        //time(&rtime);
        //std::cout << "Current time is " << ctime(&rtime);         
        //std::cout << "------------------------------------------"
        //          << std::endl;
       */
       
    }
}
