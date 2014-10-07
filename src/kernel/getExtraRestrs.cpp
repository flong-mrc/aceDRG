/* 
 * File:   getExtraRestrs.cpp
 * Author: flong
 *
 * Created on September 26, 2011, 12:31 AM
 */

#include "getExtraRestrs.h"

namespace LIBMOL
{ 
    restrListsForOnePair::restrListsForOnePair() :pairKey(NullString)      
    {    
    }
    restrListsForOnePair::restrListsForOnePair(const restrListsForOnePair & tRLFOP)
    :pairKey(tRLFOP.pairKey)
    {
        for (std::vector<int>::const_iterator iSo = tRLFOP.sortKeys.begin();
                iSo !=tRLFOP.sortKeys.end(); iSo++)
        {
            sortKeys.push_back(*iSo);
        }
        for (std::vector<ResidueID>::const_iterator iRe = tRLFOP.resPair.begin();
                iRe != tRLFOP.resPair.end(); iRe++)
        {
            resPair.push_back(*iRe);
        }
        for (std::vector<RestraitListSet>::const_iterator iRL =
                tRLFOP.restrLists.begin();
                iRL != tRLFOP.restrLists.end(); iRL++)
        {
            restrLists.push_back(*iRL);
        }
    }
    
    bool SortByKeys(const restrListsForOnePair  & obj1,
                    const restrListsForOnePair  & obj2)
    {
        
        if( obj1.sortKeys[0] < obj2.sortKeys[0])
        {
            return true;
        }
        if (obj1.sortKeys[0] > obj2.sortKeys[0])
        {
            return false;
        }
        
        if( obj1.sortKeys[1] < obj2.sortKeys[1])
        {
            return true;
        }
        
        if (obj1.sortKeys[0] > obj2.sortKeys[1])
        {
            return false;
        }       
          
   
        return false;
    }
    
    GetExtraRestrs::GetExtraRestrs() : PDBObj(NullPoint),
            restrDictObj(NullPoint),
            monomerList(NullPoint),
            allUserDefAct(NullPoint),
            outRestrFileName(NullPoint),
            itsDim(3), 
            itsSearchMode(AUTO),
            itsActType(PAIRS),
            itsErrInfo(NullString),
            itsErrLevel(0)   
    { 
       // Bases of DNA and RNA 
       DnaRnaBaseList.push_back("A");
       DnaRnaBaseList.push_back("T");
       DnaRnaBaseList.push_back("DA");
       DnaRnaBaseList.push_back("DT");
       DnaRnaBaseList.push_back("U");
       DnaRnaBaseList.push_back("DU");
       DnaRnaBaseList.push_back("C");
       DnaRnaBaseList.push_back("DC");
       DnaRnaBaseList.push_back("G");
       DnaRnaBaseList.push_back("DG");
       
       // Multinmap of DNA and RNA base pairs
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("A", "T"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("A", "U"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DA", "DT"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DU", "DT"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("C", "G"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DC", "DG"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("G", "C"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DG", "DC"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("T", "A"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DT", "DA"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("U", "A"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DU", "DA"));
       
    }
    
   /* GetExtraRestrs::GetExtraRestrs(Name inPDBName, Name inDictName, 
            Name outFileName, Name inMonoList, Name uDefF) :itsErrInfo(""),
            itsErrLevel(ZeroInt),
            outRestrFileName(outFileName)
    {    
        PDBObj        = new PDBFile(inPDBName, std::ios::in);
        restrDictObj  = new ExtraRestrDictFile(inDictName);
        allfoundRestrLists = new RestraintLists();
        monomerList   = new  MonomerLib(inMonoList.c_str());
        userDefAct    = new  getUserDefAction(uDefF);
    }
    */
    
    GetExtraRestrs::GetExtraRestrs(FileName inPDBName, 
                                   FileName inDictName, 
                                   FileName outFileName, 
                                   FileName inMonoList, 
                                   FileName uDefF) :
            outRestrFileName(outFileName),
            itsDim(3),     
            itsSearchMode(AUTO),
            itsActType(PAIRS),
            itsErrInfo(""),
            itsErrLevel(ZeroInt)
    {      
        
        PDBObj             = new PDBFile(inPDBName, std::ios::in);
        
        restrDictObj       = new ExtraRestrDictFile(inDictName);
        
        monomerList        = new MonomerLib(inMonoList);
        if(uDefF)
        {
            allUserDefAct      = new getUserDefAction(uDefF);
        }
        else 
        {
            allUserDefAct = NULL;
        }
        
        DnaRnaBaseList.push_back("A");
        DnaRnaBaseList.push_back("T");
        DnaRnaBaseList.push_back("DA");
        DnaRnaBaseList.push_back("DT");
        DnaRnaBaseList.push_back("U");
        DnaRnaBaseList.push_back("DU");
        DnaRnaBaseList.push_back("C");
        DnaRnaBaseList.push_back("DC");
        DnaRnaBaseList.push_back("G");
        DnaRnaBaseList.push_back("DG");
        
      // Multinmap of DNA and RNA base pairs
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("A", "T"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("A", "U"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DA", "DT"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DU", "DT"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("C", "G"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DC", "DG"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("G", "C"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DG", "DC"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("T", "A"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DT", "DA"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("U", "A"));
       DnaRnaBasePairs.insert(std::pair<std::string, std::string>("DU", "DA"));
       
    } 
    
    GetExtraRestrs::~GetExtraRestrs()
    {
        if (PDBObj != NULL)
        {
            delete PDBObj;
            PDBObj = NULL;
        }
        if(restrDictObj != NULL)
        {
            delete restrDictObj;
            restrDictObj = NULL;
        }
        if(monomerList != NULL)
        {
            delete monomerList;
            monomerList = NULL;
        }
        if(allUserDefAct != NULL)
        {
            delete allUserDefAct;
            allUserDefAct = NULL;
        }
    }
 
    int GetExtraRestrs::getDim()
    {
        return itsDim;
    }
    void GetExtraRestrs::setDim(int tN)
    {
        itsDim = tN;
    }
    
    SearchMode GetExtraRestrs::getSearchMode()
    {
        return itsSearchMode;
    }    
    void GetExtraRestrs::setSearchMode(SearchMode tS)
    {
        itsSearchMode = tS;
    }
    
    ActType GetExtraRestrs::getActType()
    {
        return itsActType;
    }
    void  GetExtraRestrs::setActType(ActType tA)
    {
        itsActType = tA;
    }
    
    std::string GetExtraRestrs::getErrInfo() const
    {
        return itsErrInfo;
    }
    void GetExtraRestrs::setErrInfo(std::string tE)
    {
        itsErrInfo.append(tE);
        itsErrInfo.append("\n");
    }
    
    int GetExtraRestrs::getErrLevel() const
    {
        return itsErrLevel;
    }
    void GetExtraRestrs::setErrLevel(int tN)
    {
        itsErrLevel = tN;
    }

    void GetExtraRestrs::execute()
    {
        
        //std::cout << PDBObj->allModels.size() << std::endl;
        //std::cout << monomerList->monomerGroups.size() << std::endl;
        if (PDBObj->allModels.size() && monomerList->monomerGroups.size() )
        {
            //std::cout << "setup group type " << std::endl;
            setResidueGroupType();
            //std::cout << "finishing setup group type " << std::endl;
            if (allUserDefAct != NULL)
            {
                
                if ((int)allUserDefAct->userDefActs.size() > 0)
                {
                    
                    setSearchMode(allUserDefAct->userDefActs[0].searchMode);
                    setActType(allUserDefAct->userDefActs[0].actType);
                    
                    if(getSearchMode()==AUTO)
                    {
                        autoFindExtraRestraints();
                    }
                    else
                    {
                        
                        findExtraRestraints();
                    }
                }
                else
                {
                    
                    autoFindExtraRestraints();       
                }
            }
            else 
            {  
                autoFindExtraRestraints();
            }
            
            //std::cout << "Its search Mode " << itsSearchMode << std::endl;
            //std::cout << "Check all restraints found" << std::endl;
           
         
            checkAllRestrFound();
            
            // already put into checkAllRestrFound()
            //if (allfoundRestrLists.size() !=0)
            //{
            //    reOrderAllRestrFound();
            //    outInstructFile();
            //}
            
            if (outRestrFileName)
            {
                try
                {
                    std::cout << "output file name " 
                            << outRestrFileName << std::endl;
                    //outputExtraRestraints(0);
                }
                catch (outRestrException)
                {
                    std::string tS("Failed in output new restraints\n");
                    setErrInfo(tS);
                }
            }
        }
    }
    
    void GetExtraRestrs::autoFindExtraRestraints()
    {
        // Establish neighbor lists for searching pairs or 
        // neighbor atoms 
        
        // How to set the following ?
        REAL tCellL = 3.5;
       
        REAL tCellShell = 1.0;
        
        NeighbList   tNBListOfSystem(PDBObj->allAtomList, 
                                     getDim(), tCellL, tCellShell);
        
        tNBListOfSystem.building(1);
        
        // find all of base pairs in the system
        autoFindBasePairList(tNBListOfSystem);
        //find all of additional restraints
      
    }
    
    void GetExtraRestrs::autoFindBasePairList(NeighbList & tNBList)
    {
        // search all residues in the pdb and find the
        // corresponding DNA/RNA pairs
        //std::cout << "Search base-pair restraints" << std::endl;
        
        for (int iM=0; iM < (int)PDBObj->allModels.size(); iM++)
        {
            for (int iCh=0; iCh < (int)PDBObj->allModels[iM].chains.size();
                    iCh++)
            {
                for (std::vector<Residue>::iterator iR 
                        =PDBObj->allModels[iM].chains[iCh].residues.begin();
                        iR != PDBObj->allModels[iM].chains[iCh].residues.end();
                        iR++ ) 
                {
                    //std::cout << "Residue seqNum " <<  iR->getSeqNum() << std::endl;
                    if(checkDNARNABase(iR->getName()))
                    {
                        // *iR is a DNA or RNA residue 
                        ResidueID tResID;
                        tResID.resName = TrimSpaces(iR->getName());
                        tResID.seqNum  = iR->getSeqNum();
                        tResID.chainID = TrimSpaces(iR->getChainID());
                        std::vector<Residue>    aResiduePair;
                        
                        findOneResidue(tResID.seqNum, tResID.chainID, aResiduePair);
                        
                        if (aResiduePair.size() !=0)
                        {
                            std::string tID = tResID.resName + "_"
                                    +TrimSpaces(IntToStr(tResID.seqNum)) + "_"
                                    +tResID.chainID;
                            
                            //std::cout << std::endl 
                            //          << "================================" 
                            //        << std::endl;
                            // std::cout << "Residue " << tID << std::endl;
                            // std::cout << "================================" 
                            //        << std::endl;
                            
                            
                            if (tNBList.residueNBList.find(tID) !=
                                tNBList.residueNBList.end())
                            {
                                //std::cout << "NB residue size " 
                                //          << tNBList.residueNBList[tID].size()
                                //          << std::endl;
                                //std::cout<<"They are : " << std::endl;
                                // Search *iR's residueNBList 
                                for (std::list<ResidueID>::iterator iR2= 
                                     tNBList.residueNBList[tID].begin();
                                     iR2 != tNBList.residueNBList[tID].end();
                                     iR2++)
                                {
                                    //std::cout << "NB residue : "
                                    //          <<iR2->resName << " " << iR2->seqNum 
                                    //          << "in Chain " << iR2->chainID 
                                    //         << std::endl;        
                                   
                                    if(checkDNARNABasePair(tResID.resName, iR2->resName))
                                    {
                                        std::vector<Residue> aResiduePair2;
                                        aResiduePair2.push_back(aResiduePair[0]);
                                        findOneResidue(iR2->seqNum,
                                                   TrimSpaces(iR2->chainID),
                                                   aResiduePair2);
                                       
                                        
                                        if (aResiduePair2.size() == 2)
                                        {   
                                            // std::cout << "Base pair: " 
                                            //           << aResiduePair2[0].getName()
                                            //           << "<--->" 
                                            //           << aResiduePair2[1].getName()
                                            //           << std::endl;
                                            findAllRestraintsForOnePair(aResiduePair2);
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
    
    bool GetExtraRestrs::checkDNARNABase(Name tNa)
    {
        StrUpper(tNa);
        tNa = TrimSpaces(tNa);
        
        for (int i=0; i < (int)DnaRnaBaseList.size(); i++)
        {
            if(DnaRnaBaseList[i].compare(tNa) == 0)
            {
                return true;
            }
            
        }
        
        return false;
    }
    
    bool GetExtraRestrs::checkDNARNABasePair(Name tNa1, Name tNa2)
    {
        StrUpper(tNa1);
        StrUpper(tNa2);
        //std::cout << tNa1 << " vs  " << tNa2 << std::endl;
        //std::cout << "number of elements " << DnaRnaBasePairs.count(tNa1) 
        //          << std::endl;
        
        std::pair<std::multimap<std::string, std::string>::iterator, 
                  std::multimap<std::string, std::string>::iterator> tEq
                = DnaRnaBasePairs.equal_range(tNa1);
        std::multimap<std::string, std::string>::iterator tIt;
        for (tIt=tEq.first; tIt != tEq.second; tIt++)
        {
            if (tIt->second.compare(tNa2) ==0)
            {
                return true;
            }
        }
        
        return false;
    }
    
    void GetExtraRestrs::findExtraRestraints()
    {
        // Loop over PDB residues here and make a vector of residue pairs.
        // That avoids the repeated loops over residues in PDB in the 
        // following procedures
        for (int i =0; i <(int)allUserDefAct->userDefActs.size(); i++)
        {
            //std::cout << "actType " << allUserDefAct->userDefActs[i].actType
            //          << std::endl;
            if (allUserDefAct->userDefActs[i].actType==PAIRS ||
                   allUserDefAct->userDefActs[i].actType==BASEPAIRS )
            {
                int iPos =0;
                int iRes1, iRes2, iResB, iResE;
                ID  tChainID1 = allUserDefAct->userDefActs[i].residueSet[0].chainID;
                ID  tChainID2 = allUserDefAct->userDefActs[i].residueSet[1].chainID;
                // std::cout << tChainID1 << "-- " << tChainID2 << std::endl;
                iResB = allUserDefAct->userDefActs[i].residueSet[0].from;
                iResE = allUserDefAct->userDefActs[i].residueSet[0].to;
                iRes1 = iResB;
                while (iRes1 <= iResE)
                {
                    switch (allUserDefAct->userDefActs[i].alignMode)
                    {
                        case PARA:
                            iRes2 = allUserDefAct->userDefActs[i].residueSet[1].from 
                                    + iPos;
                            searchPairRestraints(iRes1, tChainID1, 
                                                 iRes2, tChainID2);
                            break;
                       case ANTIPARA:
                            iRes2 = allUserDefAct->userDefActs[i].residueSet[1].to-iPos; 
                            //std::cout << "iRes1 " << iRes1 <<std::endl;
                            //std::cout << "iRes2 " << iRes2 << std::endl;
                            searchPairRestraints(iRes1, tChainID1, 
                                                 iRes2, tChainID2);
                            break;
                       default:
                            continue;
                    }
                    iPos++;
                    iRes1 = iResB + iPos;
                }
            }
            else if (allUserDefAct->userDefActs[i].actType==ADDRESTRAINTS)
            {
                
                for (std::vector<ChainSect>::iterator iC=
                        allUserDefAct->userDefActs[i].residueSet.begin();
                        iC != allUserDefAct->userDefActs[i].residueSet.end();
                        iC++)
                {
                    std::vector<Residue> tResidueSet;
                    buidResidueSet(*iC, tResidueSet);
                        
                    if(tResidueSet.size() >0)
                    {
                        for (std::vector<Residue>::iterator iR =
                                tResidueSet.begin(); 
                                iR != tResidueSet.end(); iR++)
                        {
      
                            findAllAddRestraintsForOneResidue(iR, 
                                              allUserDefAct->userDefActs[i].label);
                        }
                    }
                    
                }
            }
        }
    }
    
    void GetExtraRestrs::searchPairRestraints(int i1, ID tID1, 
                                              int i2, ID tID2)
    {
        std::vector<Residue>    aResiduePair;
        findOneResidue(i1, tID1, aResiduePair);
        findOneResidue(i2, tID2, aResiduePair);
        
        if (aResiduePair.size() == 2)
        {
            findAllRestraintsForOnePair(aResiduePair);
        }
    }
 
    void GetExtraRestrs::findOneResidue(int tRes, ID tChainID,
            std::vector<Residue>  & tPair)
    {
           //std::cout << "target residue " << tRes << std::endl;
           //std::cout << "Chain " << tChainID << std::endl;
            for (int iM=0; iM < (int)PDBObj->allModels.size(); iM++)
            {
                for (int iCh=0; iCh < (int)PDBObj->allModels[iM].chains.size();
                    iCh++)
                {   
                    ID aID = TrimSpaces(PDBObj->allModels[iM].chains[iCh].getID());
                    if (TrimSpaces(tChainID).compare(aID)==0)
                    {
                        for (std::vector<Residue>::iterator iR 
                            =PDBObj->allModels[iM].chains[iCh].residues.begin();
                            iR != PDBObj->allModels[iM].chains[iCh].residues.end();
                            iR++ )
                        {
                           
                            if (iR->getSeqNum() == tRes 
                                   && iR->getChainID().compare(tChainID) ==0)
                            {
                                //std::cout << "residue " << iR->getSeqNum() << std::endl;
                                //std::cout << "chain " <<  iR->getChainID() << std::endl;
                                
                                tPair.push_back(*iR);  
                                break;
                            }
                        }
                    }
                }
            }
    }
    
    void GetExtraRestrs::buidResidueSet(ChainSect & tCh,
                                        std::vector<Residue> & tResSet)
    {
        int tCount;
        if (tCh.from !=-9999)
        {
            if(tCh.from <= tCh.to)
            {
                tCount = tCh.from;
                while (tCount <= tCh.to)
                {
                    findOneResidue(tCount, tCh.chainID, tResSet);
                    tCount++;
                }
            }
            else
            {
                tCount = tCh.to;
                while(tCount <= tCh.from)
                {
                    findOneResidue(tCount, tCh.chainID, tResSet);
                    tCount++;
                }
            }
                
        }
        
        // Check
        /*
        std::cout << "Number of residues in the set " 
                << tResSet.size() << std::endl;
        for (int i=0; i < tResSet.size(); i++)
        {
            std::cout << "Residue " << tResSet[i].getName() 
                      << " seqNum " << tResSet[i].getSeqNum()
                      << " Chain "  << tResSet[i].getChainID()
                      << std::endl;
        }
       */
        
    }
    
    void GetExtraRestrs::findAllRestraintsForOnePair(std::vector<Residue> & tPair)
    {
        Name tP1 = tPair[0].getName();
        Name tS1 = IntToStr(tPair[0].getSeqNum());
        Name tC1 = tPair[0].getChainID();
        Name tP2 = tPair[1].getName();
        Name tS2 = IntToStr(tPair[1].getSeqNum());
        Name tC2 = tPair[1].getChainID(); 
        Name tLabel = tP1+"_"+tS1+"_"+tC1+":"+tP2+"_"+tS2+"_"+tC2;
        // std::cout << "Pair key " << tLabel << std::endl;
        
        for(int i=0; i < (int)restrDictObj->allRestrLists.size(); i++)
        {
            Name tP3 = restrDictObj->allRestrLists[i].monoSet[0].monomerName;
            Name tP4 = restrDictObj->allRestrLists[i].monoSet[1].monomerName;
            
            if((TrimSpaces(tP1)==TrimSpaces(tP3) &&
                    TrimSpaces(tP2)==TrimSpaces(tP4))
                    || (TrimSpaces(tP1)==TrimSpaces(tP4) &&
                    TrimSpaces(tP2)==TrimSpaces(tP3)))
            {
                // std::cout << "i is " << i << std::endl;
                // std::cout << "tP1 "  << tP1 << ": tP3 " << tP3 <<std::endl;
                // std::cout << "tP2 "  << tP2 << ": tP4 " << tP4 <<std::endl;
                // std::cout << "size of dict bond list " 
                //          << (int)restrDictObj->allRestrLists[i].restrBondList.size()
                //          << std::endl;
                
                if((int)restrDictObj->allRestrLists[i].restrBondList.size() !=0)
                {
                    for (std::vector<Bond>::iterator iBo=
                            restrDictObj->allRestrLists[i].restrBondList.begin();
                            iBo != restrDictObj->allRestrLists[i].restrBondList.end();
                            iBo++)
                    {
                       
                        searchOnePairBond(tPair, iBo, tLabel);
                        //std::cout << "Bond in one Pair" << std::endl;
                        
                    }
                }
                
                
                /*
                if(restrDictObj->allRestrLists[i].restrAngleList.size() !=0)
                {
                    for (std::vector<Angle>::iterator iAn=
                            restrDictObj->allRestrLists[i].restrAngleList.begin();
                            iAn != restrDictObj->allRestrLists[i].restrAngleList.end();
                            iAn++)
                    {
                        searchOnePairAngle(tPair, iAn, tLabel);
                    }
                }
                 */
                
                if((int)restrDictObj->allRestrLists[i].restrTorsionList.size() !=0)
                {
                    for (std::vector<Torsion>::iterator iTo=
                            restrDictObj->allRestrLists[i].restrTorsionList.begin();
                            iTo != restrDictObj->allRestrLists[i].restrTorsionList.end();
                            iTo++)
                    {
                        searchOnePairTorsion(tPair, iTo, tLabel);
                        //std::cout << "Torsion in one Pair" << std::endl;
                        
                    }
                }
                
                
                if((int)restrDictObj->allRestrLists[i].restrChiralList.size() !=0)
                {
                    for (std::vector<Chiral>::iterator iCh=
                            restrDictObj->allRestrLists[i].restrChiralList.begin();
                            iCh != restrDictObj->allRestrLists[i].restrChiralList.end();
                            iCh++)
                    {
                        searchOnePairChiral(tPair, iCh, tLabel);
                       // std::cout << "Chiral in one Pair" << std::endl;
                        
                    }
                }
            }
            
        }
    }
    
    void GetExtraRestrs::findAllAddRestraintsForOneResidue(
                         std::vector<Residue>::iterator  iR,
                         std::string  tLabel)
    {
        
        std::string tKey = iR->getName() + "_"+ IntToStr(iR->getSeqNum())
                           + "_" + iR->getChainID() 
                           + "_" + tLabel;
        //std::cout << "Key is " << tKey << std::endl;
        
        for(int i=0; i < (int)restrDictObj->allRestrLists.size(); i++)
        {
            
            std::string tStr = "";
            tStr.append(restrDictObj->allRestrLists[i].label);
            std::transform(tStr.begin(), 
                           tStr.end(), 
            tStr.begin(), ::toupper);
            
            if(tStr.compare(tLabel)==0)
            {
                // only one monomer set in this case
                
                if((restrDictObj->allRestrLists[i].monoSet[0].monomerName
                        ==iR->getName()|| 
                        restrDictObj->allRestrLists[i].monoSet[0].monomerName.find(".")
                        !=std::string::npos) && 
                        (iR->getGroupID()==
                        restrDictObj->allRestrLists[i].monoSet[0].groupName))
                {
                    if(restrDictObj->allRestrLists[i].restrBondList.size()
                            != 0)
                    {
                        for (std::vector<Bond>::iterator iBo=
                            restrDictObj->allRestrLists[i].restrBondList.begin();
                            iBo != restrDictObj->allRestrLists[i].restrBondList.end();
                            iBo++)
                        {
                            searchAddBonds(iR, iBo, tKey);
                            // std::cout << "add Bond in one Res" 
                            // << std::endl;
                        }
                        //std::cout << "Number of bonds found for " << tKey
                        //  << " is " << allfoundRestrLists[tKey].restrBondList.size()
                        //  << std::endl;
                    }
                    if(restrDictObj->allRestrLists[i].restrTorsionList.size() !=0)
                    {
                        for (std::vector<Torsion>::iterator iTo=
                            restrDictObj->allRestrLists[i].restrTorsionList.begin();
                            iTo != restrDictObj->allRestrLists[i].restrTorsionList.end();
                            iTo++)
                        {
                            searchAddTorsions(iR, iTo, tKey);
                            //std::cout << "Torsion in one residue" 
                            //<< std::endl;
                        }
                        //std::cout << " Number of torsion for " << tKey
                        //<< " is " 
                        //<< allfoundRestrLists[tKey].restrTorsionList.size()
                        //<< std::endl;
                    }
                }
                        
            }
        }
        /*  
                 * searchAddBonds(tResidueSet);
                 * searchAddAngles(tNN);
                 * searchAddTorsions(tNN);
                 * searchAddChirals(tNN);
                 * searchAddPlanes(tNN);
                */
       
    }
    void GetExtraRestrs::searchBonds()
    {   
    }
    void GetExtraRestrs::searchAddBonds(std::vector<Residue>::iterator tR,
                            std::vector<Bond>::iterator iB,
                            Name tLabel)
    {
        Bond aBond;
        for (std::vector<Atom>::iterator iRA= iB->atoms.begin();
                iRA != iB->atoms.end();
                iRA++)
        {
            for (std::vector<Atom>::iterator iA1=tR->atoms.begin();
                        iA1 != tR->atoms.end(); iA1++)
            {
                if(TrimSpaces(iA1->getName()) == TrimSpaces(iRA->getName()))
                {
                    aBond.atoms.push_back((*iA1));
                }
            }               
        }
        
        
        if (aBond.atoms.size() ==2)
        {
            aBond.setLength();
            REAL tL = iB->getLength(true);
            REAL tSigL = iB->getSigLength();
            if (itsSearchMode==USERDEF)
            {
                aBond.setLength(tL, true);
                aBond.setSigLength(tSigL);
                allProbRestrLists[tLabel].restrBondList.push_back(aBond);
            }
        }
        
    }
    
    
    
    /*
    void GetExtraRestrs::searchAddBonds(std::vector<Bond>::iterator iB)
    {
        for (int iM=0; iM < PDBObj->allModels.size(); iM++)
        {
            for (int iCh=0; iCh < PDBObj->allModels[iM].chains.size();
                    iCh++)
            {
                for (std::vector<Residue>::iterator iRe 
                             =PDBObj->allModels[iM].chains[iCh].residues.begin();
                             iRe != PDBObj->allModels[iM].chains[iCh].residues.end();
                                iRe++)
                {
                    Bond aBond;
                    std::string tNa = 
                    restrDictObj->getRestrMonoSetInfo().MonoSet[0].monomerName;
                    if(tNa.find(".") !=std::string::npos
                       || iRe->getName() == tNa)
                    {
                       if (iRe->getGroupID() == 
                           restrDictObj->getRestrMonoSetInfo().MonoSet[0].groupName)
                       {
                           // two loops for IRe and iT
                           for (int iTT = 0; iTT < iB->atoms.size(); iTT++)
                           {
                               for (int iRR=0; iRR < iRe->atoms.size(); iRR++ )
                               {
                                   if(iB->atoms[iTT].getName() 
                                           == iRe->atoms[iRR].getName())
                                   {
                                       aBond.atoms.push_back(iRe->atoms[iRR]);
                                       break;
                                   }
                               }
                           }
                           
                       }
                    }
                    
                    if(aBond.atoms.size() == 2)
                    {
                        aBond.setLength(iB->getLength(true), true);
                        aBond.setSigLength(iB->getSigLength());
                        allfoundRestrLists->restrBondList.push_back(aBond);                        
                    }
                }
            }
        }        
    } */
    
   
    void GetExtraRestrs::searchOnePairBond(std::vector<Residue>  & tR,
                                           std::vector<Bond>::iterator iB,
                                           Name                    tLab)
    {
        Bond aBond;
        Residue R1 = tR[0];
        Residue R2 = tR[1];

        
        for (std::vector<Atom>::iterator iRA= iB->atoms.begin();
                iRA != iB->atoms.end();
                iRA++)
        {   
            if (TrimSpaces(iRA->getResName()) == TrimSpaces(R1.getName()))
            {   
                
                // the atom belongs to the first residue
                for (std::vector<Atom>::iterator iA1=R1.atoms.begin();
                        iA1 != R1.atoms.end(); iA1++)
                {
                    
                    if(TrimSpaces(iA1->getName()) == TrimSpaces(iRA->getName()))
                    {
                        aBond.atoms.push_back((*iA1));
                    }
                }   
            }
            else if (TrimSpaces(iRA->getResName()) == TrimSpaces(R2.getName()))
            {
                
                // the atom belongs to the second residue
                for (std::vector<Atom>::iterator iA2=R2.atoms.begin();
                        iA2 != R2.atoms.end(); iA2++)
                {
                    
                    if(TrimSpaces(iA2->getName()) == TrimSpaces(iRA->getName()))
                    {
                        aBond.atoms.push_back((*iA2));
                    }
                }
            }
        }
        
             
        
       if (aBond.atoms.size() ==2)
        {
            aBond.setLength();
            REAL tL = iB->getLength(true);
            REAL tSigL = iB->getSigLength();
            if (itsSearchMode==USERDEF)
            {
                aBond.setLength(tL, true);
                aBond.setSigLength(tSigL);
                allProbRestrLists[tLab].restrBondList.push_back(aBond);
            }
            else if (allUserDefAct == NULL || itsSearchMode==AUTO)
            {
                
        
                REAL tLAbs= fabs(aBond.getLength(false) - tL);
                
                if (tLAbs < 6*tSigL)
                {
                    aBond.setLength(tL, true);
                    aBond.setSigLength(tSigL);
                    allProbRestrLists[tLab].restrBondList.push_back(aBond);
                    
                    //if ((tR[0].getSeqNum()==16 && tR[1].getSeqNum()==10) ||
                    //(tR[1].getSeqNum()==10 && tR[0].getSeqNum()==16) )
                    //{
                    //    std::cout << "Key is " << tLab << " and number of possible bonds "
                    //    << allProbRestrLists[tLab].restrBondList.size() << std::endl;
                    //    std::cout << "Bond, value : " <<  aBond.getLength(false) << std::endl;
                    //    std::cout << "Bond, dict value " << tL << std::endl;
                    //    std::cout << "diff " << tLAbs << std::endl;
                    //}
                    
                }
                else
                {
                    // Strict Conditions to Form a (Base) Pair
                    // If one of restraint elements (bond, torsion
                    // and chirality etc.) has a difference between 
                    // the value and dictionary value larger than 2*sigma
                    // in the pair, the pair will be erased from 
                    // map "allProbRestrLists"
                    if (allProbRestrLists.find(tLab) !=allProbRestrLists.end())
                    {
                        
                        
                        //allProbRestrLists.erase(tLab);
                        /*
                        std::cout<< "Key " << tLab 
                                 << " is deleted from allProbRestrLists"
                                 << std::endl
                                 << " because of a conflicted bond"
                                 << std::endl;
                        */
                    }
                }
            }
            
            // std::cout << "Bond length " << aBond.getLength(true) << std::endl;
        }
    }
    
    void GetExtraRestrs::searchAngles()
    {
    }
    
    void GetExtraRestrs::searchAddAngles(int tN)
    {
    }
    /*
    void GetExtraRestrs::searchAddAngles(std::vector<Angle>::iterator iA)
    {
        for (int iM=0; iM < PDBObj->allModels.size(); iM++)
        {
            for (int iCh=0; iCh < PDBObj->allModels[iM].chains.size();
                    iCh++)
            {
                for (std::vector<Residue>::iterator iRe 
                             =PDBObj->allModels[iM].chains[iCh].residues.begin();
                             iRe != PDBObj->allModels[iM].chains[iCh].residues.end();
                                iRe++)
                {
                    Angle aAngle;
                    std::string tNa = 
                    restrDictObj->getRestrMonoSetInfo().MonoSet[0].monomerName;
                    if(tNa.find(".") !=std::string::npos
                       || iRe->getName() == tNa)
                    {
                       if (iRe->getGroupID() == 
                           restrDictObj->getRestrMonoSetInfo().MonoSet[0].groupName)
                       {
                           // two loops for IRe and iT
                           for (int iTT = 0; iTT < iA->atoms.size(); iTT++)
                           {
                               for (int iRR=0; iRR < iRe->atoms.size(); iRR++ )
                               {
                                   if(iA->atoms[iTT].getName() 
                                           == iRe->atoms[iRR].getName())
                                   {
                                       aAngle.atoms.push_back(iRe->atoms[iRR]);
                                       break;
                                   }
                               }
                           }
                           
                       }
                    }
                    
                    if(aAngle.atoms.size() == 3)
                    {
                        aAngle.setValue(iA->getValue(0), 0);
                        aAngle.setSigValue(iA->getSigValue(0), 0) ;
                        allfoundRestrLists->restrAngleList.push_back(aAngle);                        
                    }
                }
            } 
        }        
    }*/    

    void GetExtraRestrs::searchOnePairAngle(std::vector<Residue>  & tR,
                                         std::vector<Angle>::iterator iA,
                                         Name tLab)
    {
        Angle aAngle;
        Residue R1 = tR[0];
        Residue R2 = tR[1]; 
        
        for (std::vector<Atom>::iterator iRA= iA->atoms.begin();
                iRA != iA->atoms.end();
                iRA++)
        {
            if (iRA->getResName() == R1.getName())
            {
                // the atom belongs to the first residue
                for (std::vector<Atom>::iterator iA1=R1.atoms.begin();
                        iA1 != R1.atoms.end(); iA1++)
                {
                    if(iA1->getName() == iRA->getName())
                    {
                        aAngle.atoms.push_back((*iA1));
                    }
                }   
            }
            else if (iRA->getResName() == R2.getName())
            {
                // the atom belongs to the second residue
                for (std::vector<Atom>::iterator iA2=R2.atoms.begin();
                        iA2 != R2.atoms.end(); iA2++)
                {
                    if(iA2->getName() == iRA->getName())
                    {
                        aAngle.atoms.push_back((*iA2));
                    }
                }
            }
        } 
       
        if (aAngle.atoms.size() == 3)
        {
            if (itsSearchMode==USERDEF)
            {
                aAngle.setValue((*iA).getValue(true), 0);
                aAngle.setSigValue((*iA).getSigValue(0), 0);
                allProbRestrLists[tLab].restrAngleList.push_back(aAngle);
            }
            else if (itsSearchMode==AUTO)
            {
                aAngle.setValue();
                REAL tV    = iA->getValue(true);
                REAL tSigV = iA->getSigValue(true);
                REAL tDV   = fabs(aAngle.getValue(false)-tV);
                if (tDV < 5.0*tSigV)
                {
                    aAngle.setValue(tV, true);
                    aAngle.setSigValue(tSigV, true);
                    allProbRestrLists[tLab].restrAngleList.push_back(aAngle);
                }
                
            }
                
        }
       
    }

    void GetExtraRestrs::searchTorsions()
    {
    }   
    void GetExtraRestrs::searchAddTorsions(std::vector<Residue>::iterator  tR,
                               std::vector<Torsion>::iterator iT,
                               Name tLabel)
    {
        Torsion aTorsion;
        
        for (std::vector<Atom>::iterator iTA= iT->atoms.begin();
                iTA != iT->atoms.end();
                iTA++)
        {
            //std::cout << "Target atom " << iTT->getName() << std::endl;
            for (std::vector<Atom>::iterator iAA=tR->atoms.begin();
                        iAA != tR->atoms.end(); iAA++)
            {
                if(TrimSpaces(iAA->getName()) == TrimSpaces(iTA->getName()))
                {
                    //std::cout << "found atom " << iAA->getName() 
                    //        << " seq num "<< iAA->getSeqNum() << std::endl;
                    aTorsion.atoms.push_back((*iAA));
                }
            }   
           
        }
        if(aTorsion.atoms.size() == 4)
        {
            aTorsion.setValue();
            REAL tV     = iT->getValue(true);
            REAL tSigV  = iT->getSigValue();
            
            //std::cout << "Torsion value from PDB " 
            //        << aTorsion.getValue(false) 
            //        << std::endl;
            //std::cout << "Torsion value from Dictionary "
            //        << iT->getValue(true) << std::endl;
            
            if (itsSearchMode==USERDEF)
            {
                aTorsion.setValue(tV, true);
                aTorsion.setSigValue(tSigV);
                allProbRestrLists[tLabel].restrTorsionList.push_back(aTorsion);
            }
        }
    }
    /*
    void GetExtraRestrs::searchAddTorsions(std::vector<Torsion>::iterator iT)
    {
        for (int iM=0; iM < PDBObj->allModels.size(); iM++)
        {
            for (int iCh=0; iCh < PDBObj->allModels[iM].chains.size();
                    iCh++)
            {
                for (std::vector<Residue>::iterator iRe 
                             =PDBObj->allModels[iM].chains[iCh].residues.begin();
                             iRe != PDBObj->allModels[iM].chains[iCh].residues.end();
                                iRe++)
                {
                    Torsion aTorsion;
                    std::string tNa = 
                    restrDictObj->getRestrMonoSetInfo().MonoSet[0].monomerName;
                    if(tNa.find(".") !=std::string::npos
                       || iRe->getName() == tNa)
                    {
                       if (iRe->getGroupID() == 
                           restrDictObj->getRestrMonoSetInfo().MonoSet[0].groupName)
                       {
                           // two loops for IRe and iT
                           for (int iTT = 0; iTT <iT->atoms.size(); iTT++)
                           {
                               for (int iRR=0; iRR < iRe->atoms.size(); iRR++ )
                               {
                                   if(iT->atoms[iTT].getName() 
                                           == iRe->atoms[iRR].getName())
                                   {
                                       aTorsion.atoms.push_back(iRe->atoms[iRR]);
                                       break;
                                   }
                               }
                           }   
                       }
                    }
                    
                    if(aTorsion.atoms.size() == 4)
                    {
                        aTorsion.setValue(iT->getValue(true), true);
                        aTorsion.setSigValue(iT->getSigValue());
                        allfoundRestrLists->restrTorsionList.push_back(aTorsion);                        
                    }
                }
            }    
        }
    }*/
    
    void GetExtraRestrs::searchOnePairTorsion(std::vector<Residue>  & tR,
                                           std::vector<Torsion>::iterator iT,
                                           Name tLab)
    {

        // get torsion restraints from two residues
        Residue R1 = tR[0];
        Residue R2 = tR[1]; 
        //std::cout << " Torsion search : " << std::endl;
        //std::cout << "Residue 1 Name " << R1.getName() 
        //          << " Seq number  " << R1.getSeqNum() << std::endl;
        //std::cout << "Residue 2 Name " << R2.getName() 
        //          << " Seq number  " << R2.getSeqNum() << std::endl;
        
        Torsion aTorsion;
        for (std::vector<Atom>::iterator iTT= iT->atoms.begin();
                iTT != iT->atoms.end();
                iTT++)
        {
            //std::cout << "Target atom " << iTT->getName() << std::endl;
            
            if (TrimSpaces(iTT->getResName()) == TrimSpaces(R1.getName()))
            {
                for (std::vector<Atom>::iterator iAA=R1.atoms.begin();
                        iAA != R1.atoms.end(); iAA++)
                {
                    if(TrimSpaces(iAA->getName()) == TrimSpaces(iTT->getName()))
                    {
                        //std::cout << "found atom " << iAA->getName() 
                        //        << " seq num "<< iAA->getSeqNum() << std::endl;
                        aTorsion.atoms.push_back((*iAA));
                    }
                }   
            }
            else if (TrimSpaces(iTT->getResName()) == TrimSpaces(R2.getName()))
            {
                for (std::vector<Atom>::iterator iBB=R2.atoms.begin();
                        iBB != R2.atoms.end(); iBB++)
                {
                    if(TrimSpaces(iBB->getName()) == TrimSpaces(iTT->getName()))
                    {
                        // std::cout << "found atom " << iBB->getName() 
                        //        << " seq num "<< iBB->getSeqNum() << std::endl;
                        aTorsion.atoms.push_back((*iBB));
                    }
                }   
            }
        }
        
        if(aTorsion.atoms.size() == 4)
        {
            aTorsion.setValue();
            REAL tV     = iT->getValue(true);
            REAL tSigV  = iT->getSigValue();
            aTorsion.setValue(tV, true);
            
            //std::cout << "Torsion value from PDB " 
            //        << aTorsion.getValue(false) 
            //        << std::endl;
            //std::cout << "Torsion value from Dictionary "
            //          << iT->getValue(true) << std::endl;
            
            if (itsSearchMode==USERDEF)
            {
                aTorsion.setValue(tV, true);
                aTorsion.setSigValue(tSigV);
                allProbRestrLists[tLab].restrTorsionList.push_back(aTorsion);
            }
            else if (allUserDefAct == NULL || itsSearchMode==AUTO)
            {
                REAL tDV    = fabs(fabs(aTorsion.getValue(false))-tV);
                
                //if ((tR[0].getSeqNum()==210 && tR[1].getSeqNum()==209) ||
                //    (tR[1].getSeqNum()==209 && tR[0].getSeqNum()==210) )
                //{
                //    std::cout << "Key is " << tLab 
                //              << " torsion " << aTorsion.getValue(false)
                //              << " dict value " <<  fabs(aTorsion.getValue(true))
                //              << " diff " << tDV
                //              << std::endl;
                //}
                
                if(tDV < 4.0*tSigV)
                {
                    aTorsion.setValue(tV, true);
                    aTorsion.setSigValue(tSigV);
                    allProbRestrLists[tLab].restrTorsionList.push_back(aTorsion);
                    /*
                    if ((tR[0].getSeqNum()==16 && tR[1].getSeqNum()==10) ||
                    (tR[1].getSeqNum()==10 && tR[0].getSeqNum()==16) )
                    {
                    
                       std::cout << "Key is " << tLab << " and number of possible torsion "
                       << allProbRestrLists[tLab].restrTorsionList.size() << std::endl;
                    
                    }
                    */
                }
                else
                {
                    // Strict Conditions to Form a (Base) Pair
                    // If one of restraint elements (bond, torsion
                    // and chirality etc.) has a difference between 
                    // the value and dictionary value larger than 2*sigma
                    // in the pair, the pair will be erased from 
                    // map "allProbRestrLists"
                    if (allProbRestrLists.find(tLab) !=allProbRestrLists.end())
                    {
                        // allProbRestrLists.erase(tLab);
                        //std::cout<< "Key " << tLab 
                        //        << " is deleted from allProbRestrLists"
                        //        << std::endl
                        //        << " because of a conflicted torsion"
                        //        << std::endl;
                        
                    }
                }
            }
                
        }
    }
    
    void GetExtraRestrs::searchChirals()
    {
    }
    
    void GetExtraRestrs::searchAddChirals(int tN)
    {
    }
    
    /*
    void GetExtraRestrs::searchAddChirals(std::vector<Chiral>::iterator iA)
    {
        for (int iM=0; iM < PDBObj->allModels.size(); iM++)
        {
            for (int iCh=0; iCh < PDBObj->allModels[iM].chains.size();
                    iCh++)
            {
                for (std::vector<Residue>::iterator iRe 
                             =PDBObj->allModels[iM].chains[iCh].residues.begin();
                             iRe != PDBObj->allModels[iM].chains[iCh].residues.end();
                                iRe++)
                {
                    Chiral aChiral;
                    std::string tNa = 
                    restrDictObj->getRestrMonoSetInfo().MonoSet[0].monomerName;
                    if(tNa.find(".") !=std::string::npos
                       || iRe->getName() == tNa)
                    {
                       if (iRe->getGroupID() == 
                           restrDictObj->getRestrMonoSetInfo().MonoSet[0].groupName)
                       {
                           // two loops for IRe and iT
                           for (int iTT = 0; iTT < iA->atoms.size(); iTT++)
                           {
                               for (int iRR=0; iRR < iRe->atoms.size(); iRR++ )
                               {
                                   if(iA->atoms[iTT].getName() 
                                           == iRe->atoms[iRR].getName())
                                   {
                                       aChiral.atoms.push_back(iRe->atoms[iRR]);
                                       break;
                                   }
                               }
                           }
                           
                       }
                    }
                    
                    if(aChiral.atoms.size() == 4)
                    {
                        aChiral.setValue(iA->getValue(true), true);
                        aChiral.setSigValue(iA->getSigValue()) ;
                        allfoundRestrLists->restrChiralList.push_back(aChiral);                        
                    }
                }
            } 
        }        
    } */    
    
        
    void GetExtraRestrs::searchOnePairChiral(std::vector<Residue>  & tR,
                                           std::vector<Chiral>::iterator iT,
                                           Name  tLab)
    {
        // Get chiral restraints from two residues
        
        Residue R1 = tR[0];
        Residue R2 = tR[1];
        /*
        std::cout << "Residue 1 Name " << R1.getName() 
                  << " Seq number  " << R1.getSeqNum() 
                  << "Chain <" << R1.getChainID() << std::endl;
        std::cout << "Residue 2 Name " << R2.getName() 
                  << " Seq number  " << R2.getSeqNum() 
                  << "Chain " << R2.getChainID() << std::endl; 
        */
        
        Chiral aChiral;
        for (std::vector<Atom>::iterator iTT= iT->atoms.begin();
                iTT != iT->atoms.end();
                iTT++)
        {
            //std::cout << "Target atom " << iTT->getName() << std::endl;
            if (TrimSpaces(iTT->getResName()) == TrimSpaces(R1.getName()))
            {
                for (std::vector<Atom>::iterator iAA=R1.atoms.begin();
                        iAA != R1.atoms.end(); iAA++)
                {
                    if(TrimSpaces(iAA->getName()) == TrimSpaces(iTT->getName()))
                    {
                         // std::cout << "found atom " << iAA->getName() 
                         //   << " seq num "<< iAA->getSeqNum() << std::endl;
                        aChiral.atoms.push_back((*iAA));
                    }
                }   
            }
            else if (TrimSpaces(iTT->getResName()) == TrimSpaces(R2.getName()))
            {
                for (std::vector<Atom>::iterator iBB=R2.atoms.begin();
                        iBB != R2.atoms.end(); iBB++)
                {
                    if(TrimSpaces(iBB->getName()) == TrimSpaces(iTT->getName()))
                    {
                        // std::cout << "found atom " << iBB->getName() 
                        //        << " seq num "<< iBB->getSeqNum() << std::endl; 
                        aChiral.atoms.push_back((*iBB));
                    }
                }   
            }
        }
        
        if(aChiral.atoms.size() == 4)
        {                  
            REAL tV     = iT->getValue(true);
            REAL tSigV  = iT->getSigValue();  
            aChiral.setValue();
            /*
            if ((tR[0].getSeqNum()==48 && tR[1].getSeqNum()==49) ||
                    (tR[1].getSeqNum()==49 && tR[0].getSeqNum()==48) )
            {
                std::cout << "Chirality value from PDB " 
                          << aChiral.getValue(false) 
                          << std::endl;
                std::cout << "Chirality value from Dictionary "
                          << iT->getValue(true) << std::endl;  
            }
            
            */
            if (itsSearchMode==USERDEF)
            {
                aChiral.setValue(tV, true);
                aChiral.setSigValue(tSigV);      
                allProbRestrLists[tLab].restrChiralList.push_back(aChiral);
            }
            else if (allUserDefAct == NULL || itsSearchMode==AUTO)
            {
                REAL tDV    = fabs(fabs(aChiral.getValue(false))-tV);
                if(tDV < 6.0*tSigV)
                {
                    aChiral.setValue(tV, true);
                    aChiral.setSigValue(tSigV);
                    allProbRestrLists[tLab].restrChiralList.push_back(aChiral);
                    /*
                    std::cout << "Key is " << tLab 
                             << " and number of possible chirality "
                             << allProbRestrLists[tLab].restrChiralList.size() 
                             << std::endl;
                     */
                }
                else
                {
                    // Strict Conditions to Form a (Base) Pair
                    // If one of restraint elements (bond, torsion
                    // and chirality etc.) has a difference between 
                    // the value and dictionary value larger than 2*sigma
                    // in the pair, the pair will be erased from 
                    // map "allProbRestrLists"
                    if (allProbRestrLists.find(tLab) !=allProbRestrLists.end())
                    {
                        /*
                        allProbRestrLists.erase(tLab);
                        std::cout<< "Key " << tLab 
                                << " is deleted from allProbRestrLists"
                                << std::endl
                                << " because of conflicted chirality"
                                << std::endl;
                         */                         
                    }
                }
            }
        }
    }
    
    
    void GetExtraRestrs::searchPlanes()
    {
        if (restrDictObj->getRestrInfoType() == DICT_ADDTIONAL)
        {       
            // searchAddPlanes(iP);
        }
        else if(restrDictObj->getRestrInfoType() == DICT_PAIRS) 
        {
            //searchPairPlanes(iP);
        }             
    }
    
    void GetExtraRestrs::searchAddPlanes(int tN)
    {
    }
   
    void GetExtraRestrs::searchOnePairPlane(std::vector<Residue>  & tR,
                                std::vector<Plane>::iterator iP,
                                Name  tLab)
    {
        
    }   
    
    void GetExtraRestrs::setResidueGroupType()
    { 
        for (int iMo=0; iMo < (int)PDBObj->allModels.size(); iMo++)
        {
            for (int iCh=0; iCh < (int)PDBObj->allModels[iMo].chains.size();
                            iCh++)
            {
                for (std::vector<Residue>::iterator iRe = 
                        PDBObj->allModels[iMo].chains[iCh].residues.begin();
                        iRe != PDBObj->allModels[iMo].chains[iCh].residues.end();
                        iRe++)
                {
                    bool lDone = false;
                    // std::cout<< "current residue Name " << tResName << std::endl;
                    lDone = checkReducedGroups(iRe);
                    if (!lDone)
                    {
                        try
                        {
                            checkFullGroups(iRe);
                        }
                        catch (noMonoGroupException & e)
                        {
                            std::cout << "No group type is found for residue "
                                    << iRe->getName() << " of seq number "
                                    << iRe->getSeqNum() << std::endl;
                        }
                    }

                            //std::cout << "its Group ID : " <<
                            //PDBObj->allModels[iMo].chains[iCh].residues[iRe].getGroupID()
                            //        << std::endl;
                }
            }
        }   
    }
       
    bool GetExtraRestrs::checkReducedGroups(std::vector<Residue>::iterator iRe)
    {
        for(std::map<std::string, std::string>::iterator iM 
                = reducedGroups.begin();
                iM != reducedGroups.end(); iM++)
        {
            if(!TrimSpaces(iRe->getName()).compare(TrimSpaces(iM->first)))
            {
                iRe->setGroupID(TrimSpaces(iM->second));
                //std::cout << "get from reduce group. Residue  " 
                //<< iRe->getName() << " group type : " 
                //<< iRe->getGroupID() << std::endl;
                
                return true;
            }   
        } 
        return false;
    }
    
    void GetExtraRestrs::checkFullGroups (std::vector<Residue>::iterator iRe)
    {
        for(std::map<std::string, std::string>::iterator iMa 
                = monomerList->monomerGroups.begin();
                iMa != monomerList->monomerGroups.end();
                iMa++)
        {
            //std::cout<< "Monomer Name : " << iMa->first << std::endl;
            //std::cout<< "its Group type: " << iMa->second <<std::endl;
            std::string rName = TrimSpaces(iRe->getName());
            if(!rName.compare(TrimSpaces(iMa->first)))
            {
                std::string rValue= TrimSpaces(iMa->second);
                iRe->setGroupID(rValue);
                reducedGroups[rName]=rValue;
                //std::cout << "get from full group. Residue  " 
                //<< rName << " group type : " 
                // << rValue << std::endl;
                break;
            }
        }        
    }
    
    void GetExtraRestrs::checkAllRestrFound()
    {
        
       // std::cout << "allProbRestrLists.size() " << allProbRestrLists.size()
       //           << std::endl;
        
        
        if( allProbRestrLists.size() != 0)
        {
            if (itsSearchMode==AUTO || 
               (itsSearchMode==USERDEF && itsActType == PAIRS))
            {
                for (std::map<std::string, RestraitListSet>::iterator iMap 
                              = allProbRestrLists.begin(); 
                              iMap != allProbRestrLists.end(); iMap++)
                {
                     //std::cout << "Check restraint Pair : "
                     //          << iMap->first << std::endl;
                    checkRestraintsFromOnePair(iMap);
                }
                // temp print 
                /*
                std::cout <<std::endl
                          << "|------------------------------------------|" 
                          << std::endl;
                std::cout << "|         Extra Restraints found           |" 
                          << std::endl;
                std::cout << "|------------------------------------------|" 
                          << std::endl;
                for (std::map<std::string, RestraitListSet>::iterator iMap 
                          = allfoundRestrLists.begin(); 
                          iMap != allfoundRestrLists.end(); iMap++)
                 {
                     std::cout << "\n Pair " << iMap->first << std::endl;
                     std::cout << "number of its bond restraints "
                               << iMap->second.restrBondList.size() << std::endl;
                     std::cout << "number of its torsion restraints "
                               << iMap->second.restrTorsionList.size() << std::endl;
                     //std::cout << "number of its chiral restraints "
                     //          << iMap->second.restrChiralList.size() << std::endl;
                 }
                 */
                if (allfoundRestrLists.size() !=0)
                {
                    reOrderAllRestrFound();
                    outInstructFile();
                }
            }
            else 
            {
                for (std::map<std::string, RestraitListSet>::iterator iMap 
                              = allProbRestrLists.begin(); 
                              iMap != allProbRestrLists.end(); iMap++)
                {
                    copyRestraintPair(iMap);
                }
                outputAddRestraints();
                
            }
        }
    
    }
    
    void GetExtraRestrs::checkRestraintsFromOnePair(
            std::map<std::string, RestraitListSet>::iterator iM)
    {
        bool lCont= true;
        
        if (iM->second.restrBondList.size() != 0)
        {
            for (std::vector<Bond>::iterator iBond =iM->second.restrBondList.begin();
                    iBond != iM->second.restrBondList.end(); iBond++)
            {
                
                /*
                std::cout << "Bond value from PDB " 
                          << iBond->getLength(false) << std::endl;
                std::cout << "Bond value from Dictionary " 
                          << iBond->getLength(true) << std::endl;
                std::cout << "sigma " << iBond->getSigLength() << std::endl;
                */
                
                REAL tDB = fabs(iBond->getLength(false)-iBond->getLength(true));
                //std::cout << "Diff "<< tDB << std::endl;
                
                if(tDB > 6*iBond->getSigLength())
                {
                    lCont = false;
                    break;
                }
            }
        }
        if (lCont)
        {
            if(iM->second.restrTorsionList.size() != 0)
            {
                for (std::vector<Torsion>::iterator iTorsion 
                        = iM->second.restrTorsionList.begin();
                        iTorsion != iM->second.restrTorsionList.end();
                        iTorsion++)
                {
                    /*
                    std::cout << "Torsion value from PDB "
                            << iTorsion->getValue(false) 
                            << std::endl;
                    std::cout << "Torsion value from dictionary "
                            << iTorsion->getValue(true) 
                            << std::endl;
                    std::cout << "Sigma " << iTorsion->getSigValue()
                            << std::endl;
                    */
                    
                    REAL tDT = fabs(iTorsion->getValue(false)
                            -iTorsion->getValue(true));
                    // std::cout << "Diff " << tDT << std::endl;
                    
                    if(tDT > 5.0*iTorsion->getSigValue())
                    {
                        lCont = false;
                        break;
                    }
                }
            }
        }
        /*
        if(lCont)
        {
            if(iM->second.restrChiralList.size() != 0)
            {
                for (std::vector<Chiral>::iterator iChiral 
                        = iM->second.restrChiralList.begin();
                        iChiral != iM->second.restrChiralList.end();
                        iChiral++)
                {
                    
                    std::cout << "Chiral value from PDB "
                            << iChiral->getValue(false) 
                            << std::endl;
                    std::cout << "Chiral value from dictionary "
                            << iChiral->getValue(true) 
                            << std::endl;
                    std::cout << "Sigma " << iChiral->getSigValue()
                            << std::endl;
                    
                    REAL tDC = fabs(fabs(iChiral->getValue(false))
                            -iChiral->getValue(true));
                    // std::cout << "Diff " << tDC << std::endl;
                    if(tDC > 6.0*iChiral->getSigValue())
                    {
                        lCont = false;
                        break;
                    }
                }
            }
         
        
        }
        
        */
        /*
        std::cout << "number of bond restraints " 
                << iM->second.restrBondList.size() << std::endl;
        std::cout << "number of Torsion restraints "
                << iM->second.restrTorsionList.size() << std::endl;
        std::cout << "Number of chirality restraints "
                << iM->second.restrChiralList.size() << std::endl;
        */
        /*   iM->second.restrChiralList.size() > 1  */
        //if(lCont && iM->second.restrBondList.size() >=2 &&
        //      iM->second.restrTorsionList.size() !=0)
        //{
            
        checkRestrListSize(iM);
            //std::cout<< "Map entry " << iM->first 
            //         << " is copied " << std::endl;
            //copyRestraintPair(iM);
        //}
    }
    
    void GetExtraRestrs::checkRestrListSize(
               std::map<std::string, RestraitListSet>::iterator iM)
    {
        std::vector<std::string> tStrgrp = StrTokenize(iM->first, ':');
        
        if ((int)tStrgrp.size()==2)
        {
            std::vector<std::string> tStrgrp1 = StrTokenize(tStrgrp[0],'_');
            std::vector<std::string> tStrgrp2 = StrTokenize(tStrgrp[1],'_');
            
            if((tStrgrp1[0].find("A") != std::string::npos &&
                tStrgrp2[0].find("T") != std::string::npos) 
               || 
               (tStrgrp1[0].find("T") != std::string::npos &&
                tStrgrp2[0].find("A") != std::string::npos) 
               ||
               (tStrgrp1[0].find("A") != std::string::npos &&
                tStrgrp2[0].find("U") != std::string::npos) 
               || 
               (tStrgrp1[0].find("U") != std::string::npos &&
                tStrgrp2[0].find("A") != std::string::npos) )
            {
                // base pair A-T or DA-DT or A-U or DA-DU
                if ((int)iM->second.restrBondList.size() ==2 &&
                    (int)iM->second.restrTorsionList.size() ==1 &&
                    (int)iM->second.restrChiralList.size()  == 2 )
                    
                //if ((int)iM->second.restrBondList.size() ==2 &&
                //    (int)iM->second.restrTorsionList.size() ==1 )
                {
                    copyRestraintPair(iM);
                }
            }
            else if ((tStrgrp1[0].find("C") != std::string::npos &&
                tStrgrp2[0].find("G") != std::string::npos) 
               || 
               (tStrgrp1[0].find("G") != std::string::npos &&
                tStrgrp2[0].find("C") != std::string::npos))
            {
                // base pair C-G or DC-DG 
                // if ((int)iM->second.restrBondList.size() ==3 &&
                //    (int)iM->second.restrTorsionList.size() ==1 &&
                //    (int)iM->second.restrChiralList.size()  == 2 )
                if ((int)iM->second.restrBondList.size() ==3 &&
                    (int)iM->second.restrTorsionList.size() ==1 )
                {
                    copyRestraintPair(iM);
                }
            }
            
        }
    }
    
    
    void GetExtraRestrs::copyRestraintPair(
               std::map<std::string, RestraitListSet>::iterator iM)
    {
        allfoundRestrLists.insert(std::pair<std::string, RestraitListSet> 
                                  (iM->first, iM->second));
        // std::cout << "Map size " << (int)allfoundRestrLists.size() << std::endl;
        
    }
    
    void GetExtraRestrs::outputExtraRestraints(int tMode)
    {
        std::ofstream outputFile;
        if (tMode == 0)
        {
            outputFile.open(outRestrFileName);
        }
        else if (tMode ==1)
        {
            outputFile.open("testSuggestRestraints");
        }
        if(outputFile.is_open())
        {
            for (std::vector<restrListsForOnePair>::iterator  iTP=
                    sortedAllRestLists.begin();
                    iTP != sortedAllRestLists.end(); iTP++)
            {
                if(tMode == 1)
                {
                    outputFile << std::endl 
                            << "Possible Restraints in Residue Pair: "
                            << iTP->pairKey << std::endl;
                }
                //std::cout << "Pair in outfile " << iTP->pairKey << std::endl;
                for (std::vector<RestraitListSet>::iterator iRL=
                        iTP->restrLists.begin();
                        iRL != iTP->restrLists.end(); iRL++)
                {
                    if (iRL->restrBondList.size() != 0)
                    {
                        outputBondRestraints(outputFile, 
                                         iRL->restrBondList, tMode);
                    }
                    if(iRL->restrTorsionList.size() !=0)
                    {
                        outputTorsionRestraints(outputFile,
                                          iRL->restrTorsionList, tMode);
                    }
                    if(iRL->restrChiralList.size() !=0)
                    {
                        outputChiralRestraints(outputFile,
                                           iRL->restrChiralList, tMode);
                    }
                    if(iRL->restrAngleList.size()  !=0)
                    {
                        outputAngleRestraints(outputFile,
                                          iRL->restrAngleList, tMode);
                    }
                }
            }
            outputFile.close();    
        }
    }
    

    void GetExtraRestrs::outputAddRestraints()
    {
        std::ofstream outputFile;
        outputFile.open(outRestrFileName);
        if(outputFile.is_open())
        {
            std::cout <<std::endl
                      << "|-----------------------------------------|" 
                      << std::endl;
            std::cout << "|   Additional Restraints found in:       |"
                      << std::endl;
            std::cout << "|-----------------------------------------|" 
                      << std::endl;
            
            int tMode = 0;
            for (std::map<std::string, RestraitListSet>::iterator iMp=
                    allfoundRestrLists.begin(); 
                    iMp !=allfoundRestrLists.end(); iMp++)
            {
                std::cout << iMp->first << std::endl;
                if((int)iMp->second.restrBondList.size() > 0)
                {
                    outputBondRestraints(outputFile, iMp->second.restrBondList, tMode);
                }
                if((int)iMp->second.restrTorsionList.size() >0)
                {
                    outputTorsionRestraints(outputFile,iMp->second.restrTorsionList, tMode);
                }
            }
        }
        outputFile.close(); 
    }
    
    void GetExtraRestrs::outputBondRestraints(std::ofstream & tO,
                                              std::vector<Bond> & tBList,
                                              int tMode)
    {    
        for (int i=0; i < (int)tBList.size(); i++)
        {
            Atom tA1 = tBList[i].atoms[0];
            Atom tA2 = tBList[i].atoms[1];
            tO << "exte dist first chain " << tA1.getChainID() 
               << " resi " << tA1.getSeqNum() 
               << " ins  " << tA1.getInsCode()
               << " atom " << tA1.getName() 
               << " second chain " << tA2.getChainID()
               << " resi " << tA2.getSeqNum()
               << " ins "  << tA2.getInsCode()
               << " atom " << tA2.getName();
            if (tMode ==0)
            {
                tO << " value " << tBList[i].getLength(true)
                   << " sigma " << tBList[i].getSigLength()
                   << " type 1 ";
            }
            
            tO << std::endl;
                
        }
            
    }
    
    void GetExtraRestrs::outputAngleRestraints(std::ofstream & tO,
                                              std::vector<Angle> & tAList,
                                              int tMode)
    {
        for (int i=0; i < (int)tAList.size(); i++)
        {
            Atom tA1 = tAList[i].atoms[0];
            Atom tA2 = tAList[i].atoms[1];
            Atom tA3 = tAList[i].atoms[2];
            tO << "exte angle first chain " << tA1.getChainID() 
               << " resi " << tA1.getSeqNum() 
               << " ins  " << tA1.getInsCode()
               << " atom " << tA1.getName() 
               << " second chain " << tA2.getChainID()
               << " resi " << tA2.getSeqNum()
               << " ins "  << tA2.getInsCode()
               << " atom " << tA2.getName()
               << " third chain " << tA3.getChainID()
               << " resi " << tA3.getSeqNum()
               << " ins "  << tA3.getInsCode()
               << " atom " << tA3.getName();
            if (tMode ==0)
            {
               tO << " value " << tAList[i].getValue(true)
               << " sigma " << tAList[i].getSigValue(true)
               << " type 1 ";
            }
            tO << std::endl;
        }
    }
        
    void GetExtraRestrs::outputTorsionRestraints(std::ofstream & tO,
                                                 std::vector<Torsion> & tTList,
                                                 int tMode)
    {
        for (int i=0; i < (int)tTList.size(); i++)        
        {
            Atom tA1 = tTList[i].atoms[0];
            Atom tA2 = tTList[i].atoms[1];
            Atom tA3 = tTList[i].atoms[2];
            Atom tA4 = tTList[i].atoms[3];
            tO << "exte torsion first chain " << tA1.getChainID() 
               << " resi " << tA1.getSeqNum() 
               << " ins  " << tA1.getInsCode()
               << " atom " << tA1.getName() 
               << " second chain " << tA2.getChainID()
               << " resi " << tA2.getSeqNum()
               << " ins "  << tA2.getInsCode()
               << " atom " << tA2.getName()
               << " third chain " << tA3.getChainID()
               << " resi " << tA3.getSeqNum()
               << " ins "  << tA3.getInsCode()
               << " atom " << tA3.getName()                    
               << " fourth chain " << tA4.getChainID()
               << " resi " << tA4.getSeqNum()
               << " ins "  << tA4.getInsCode()
               << " atom " << tA4.getName();
            if(tMode==0)
            {
               tO << " value " << tTList[i].getValue(true)
                  << " sigma " << tTList[i].getSigValue()
                  << " type 1 ";
            }
            
            tO << std::endl;
        }
    }
    
    void GetExtraRestrs::outputChiralRestraints(std::ofstream & tO,
                                                std::vector<Chiral> & tCList,
                                                int tMode)
    {
        for (int i=0; i < (int)tCList.size(); i++)        
        {
            Atom tA1 = tCList[i].atoms[0];
            Atom tA2 = tCList[i].atoms[1];
            Atom tA3 = tCList[i].atoms[2];
            Atom tA4 = tCList[i].atoms[3];
            tO << "exte chiral first chain " << tA1.getChainID() 
               << " resi " << tA1.getSeqNum() 
               << " ins  " << tA1.getInsCode()
               << " atom " << tA1.getName() 
               << " second chain " << tA2.getChainID()
               << " resi " << tA2.getSeqNum()
               << " ins "  << tA2.getInsCode()
               << " atom " << tA2.getName()
               << " third chain " << tA3.getChainID()
               << " resi " << tA3.getSeqNum()
               << " ins "  << tA3.getInsCode()
               << " atom " << tA3.getName()                    
               << " fourth chain " << tA4.getChainID()
               << " resi " << tA4.getSeqNum()
               << " ins "  << tA4.getInsCode()
               << " atom " << tA4.getName();
            if(tMode ==0)
            {
               tO<< " value " << tCList[i].getValue(true)
               << " sigma " << tCList[i].getSigValue()
               << " type 1 " ;
            }
            tO << std::endl;
        }
    }
    
    void GetExtraRestrs::reOrderAllRestrFound()
    {    
        std::map<std::string, RestraitListSet> tOrderedRLists;
        
        for (std::map<std::string, RestraitListSet>::iterator  iTP=
                    allfoundRestrLists.begin();
                    iTP != allfoundRestrLists.end(); iTP++)
        {
            // put all map terms into the vector to be sorted
            restrListsForOnePair tAllRestListOnePair;
                   
            //std::cout << "Pair restraints " << iTP->first << std::endl;
            //for (std::vector<Bond>::iterator iBB=iTP->second.restrBondList.begin();
            //        iBB !=iTP->second.restrBondList.end(); iBB++)
            //{
            //    std::cout << "Bond: " << std::endl; 
            //    std::cout << "Atom 1 " << iBB->atoms[0].getName() 
            //            << " its residue seqNum " << iBB->atoms[0].getSeqNum()
            //            << " its chainID " << iBB->atoms[0].getChainID() << std::endl;
            //    std::cout << "Atom 2 " << iBB->atoms[1].getName() 
            //            << " its residue seqNum " << iBB->atoms[1].getSeqNum()
            //            << " its chainID " << iBB->atoms[1].getChainID() << std::endl;    
            //}
            
            // Ensure the first chain ID is higher than the second chain ID
            std::vector<std::string> tStrgrp = StrTokenize(iTP->first, ':');
            std::vector<ResidueID> tTwoIDs;
            int tSize = (int)tStrgrp.size();
            if(tSize == 2)
            {
                for (int i =0; i < tSize; i++)
                {
                    ResidueID tResID;
                    //std::cout << tStrgrp[i] << std::endl;
                    std::vector<std::string> tStrgrp2 = StrTokenize(tStrgrp[i],'_');
                    if (tStrgrp2.size() == 3)
                    {
                        tResID.resName = TrimSpaces(tStrgrp2[0]);
                        tResID.seqNum  =  StrToInt(TrimSpaces(tStrgrp2[1]));
                        tResID.chainID = tStrgrp2[2];
                        //std::cout << "Residue: name " << tResID.resName
                        //      << " seqNum "  << tResID.seqNum
                        //      << " chainID " <<  tResID.chainID << std::endl;
                        tTwoIDs.push_back(tResID);    
                    }
                }
                char a1 = tTwoIDs[0].chainID[0];
                char a2 = tTwoIDs[1].chainID[0];
                //std::cout << "chain " << tTwoIDs[0].chainID << " int  "
                //          << (int)a1
                //          << " and chain " << tTwoIDs[1].chainID
                //          << " int  "  << (int)a2
                //          << std::endl;
                
                
                if( (int)a1 < (int)a2 )
                {
                    tAllRestListOnePair.resPair.push_back(tTwoIDs[1]);
                    tAllRestListOnePair.resPair.push_back(tTwoIDs[0]);
                    std::string tKey =  tTwoIDs[1].resName + "_" 
                            + TrimSpaces(IntToStr(tTwoIDs[1].seqNum))
                            + "_" + tTwoIDs[1].chainID + ":"  
                            + tTwoIDs[0].resName + "_" 
                            + TrimSpaces(IntToStr(tTwoIDs[0].seqNum))
                            + "_" + tTwoIDs[0].chainID;
                    tAllRestListOnePair.pairKey = tKey;
                    const_cast<std::string &>(iTP->first) = tKey;
                }
                else
                {
                    tAllRestListOnePair.resPair.push_back(tTwoIDs[0]);
                    tAllRestListOnePair.resPair.push_back(tTwoIDs[1]);
                    tAllRestListOnePair.pairKey = iTP->first;
                }
                tAllRestListOnePair.restrLists.push_back(iTP->second);
            }
            if (tOrderedRLists.find(tAllRestListOnePair.pairKey)
                    ==tOrderedRLists.end())
            {   
                tOrderedRLists.insert(std::pair<std::string, RestraitListSet>
                (tAllRestListOnePair.pairKey, iTP->second));
                sortedAllRestLists.push_back(tAllRestListOnePair);
            }
            
        }

        
        // 
        for (std::vector<restrListsForOnePair>::iterator iS =sortedAllRestLists.begin(); 
                iS != sortedAllRestLists.end(); iS++)
        {  
            iS->sortKeys.clear();
            char b = iS->resPair[0].chainID[0];
            iS->sortKeys.push_back((int)b);
            //std::cout << "char " << iS->resPair[0].chainID[0]
            //          << " with int " << (int)b << std::endl;
            iS->sortKeys.push_back(iS->resPair[0].seqNum);
        }
        
        //std::cout <<  "Key order before sorting " << std::endl;
        
        //for (std::map<std::string, RestraitListSet>::iterator iM =
        //        tOrderedRLists.begin(); 
        //        iM != tOrderedRLists.end(); iM++)
        //{
        //    std::cout << iM->first << std::endl;
        //}
        
        // std::cout << "\nIn list to be sorted " <<std::endl;
        //for (std::vector<restrListsForOnePair>::iterator iBB=sortedAllRestLists.begin();
        //     iBB !=sortedAllRestLists.end(); iBB++)
        //{
        //    std::cout<<"Pair " << iBB->pairKey << std::endl;
        //    std::cout << "sortKey " << iBB->sortKeys[0] 
        //              << " and " << iBB->sortKeys[1] << std::endl;
            
        //    for (std::vector<Bond>::iterator iBo = iBB->restrLists[0].restrBondList.begin();
        //            iBo != iBB->restrLists[0].restrBondList.end(); iBo++)
        //    {
        //        std::cout << "Bond: " << std::endl; 
        //        std::cout << "Atom 1 " << iBo->atoms[0].getName() 
        //                  << " its residue seqNum " << iBo->atoms[0].getSeqNum()
        //                  << " its chainID " << iBo->atoms[0].getChainID() << std::endl;
        //        std::cout << "Atom 2 " << iBo->atoms[1].getName() 
        //                  << " its residue seqNum " << iBo->atoms[1].getSeqNum()
        //                  << " its chainID " << iBo->atoms[1].getChainID() << std::endl;
        //    }
        //    std::cout <<std::endl;
        //}     
        
        
        outputExtraRestraints(0);
        //std::cout << "sorting begins " << std::endl;
        std::sort(sortedAllRestLists.begin(),
                  sortedAllRestLists.end(),
                  SortByKeys);
                      
        std::cout <<std::endl
                  << "|-----------------------------------------|" 
                  << std::endl;
        std::cout << "|        Extra Restraints found           |"
                  << std::endl
                  << "|        in following base pairs:         |" 
                  << std::endl;
        std::cout << "|-----------------------------------------|" 
                  << std::endl;
       
        // check if "sortedAllRestLists" contains what it should
        for (std::vector<restrListsForOnePair>::iterator iS =
                sortedAllRestLists.begin(); 
                iS != sortedAllRestLists.end(); iS++)
        {
            std::cout << " Residue Pair: " 
                    << iS->pairKey << std::endl;
        }
        
        /*
        for (std::map<std::string, RestraitListSet>::iterator iMap 
                    = allfoundRestrLists.begin(); 
                    iMap != allfoundRestrLists.end(); iMap++)
        {
            std::cout << "\n Pair " << iMap->first << std::endl;
            std::cout << "number of its bond restraints "
                      << iMap->second.restrBondList.size() << std::endl;
            std::cout << "number of its torsion restraints "
                      << iMap->second.restrTorsionList.size() << std::endl;
            std::cout << "number of its chiral restraints "
                      << iMap->second.restrChiralList.size() << std::endl;
        }
        */
    }

    void GetExtraRestrs::outInstructFile()
    {
        std::ofstream outputFile;
        
        outputFile.open("testSuggestRestraints.txt");
       
        if(outputFile.is_open())
        {
            std::string ID1;
            int from1=-1, to1=-1;
            std::string ID2;
            int from2=-1, to2=-1;
            int iCout       = 0;
            bool lStart     = true;
            
            for (std::vector<restrListsForOnePair>::iterator iS =
                sortedAllRestLists.begin(); 
                iS != sortedAllRestLists.end(); iS++)
            {
                if(iCout ==0 && lStart)
                {
                    //std::cout << "sec1 seqNum " << iS->resPair[0].seqNum << std::endl;
                    ID1    = iS->resPair[0].chainID;
                    from1  = iS->resPair[0].seqNum;
                    to1    = iS->resPair[0].seqNum;
                    ID2    = iS->resPair[1].chainID;
                    from2  = iS->resPair[1].seqNum;
                    to2    = iS->resPair[1].seqNum;
                    lStart = false;
                    iCout++;
                }
                else if ((!lStart && iS->resPair[0].seqNum != from1 + iCout) 
                        || iS->resPair[0].chainID.compare(ID1) !=0)
                {
                     
                    // std::cout << "sec2 seqNum " << iS->resPair[0].seqNum << std::endl;
                    //std::cout << "From " << from1 << std::endl;
                    //std::cout << from1 + iCout << std::endl;
                    //std::cout << "chain comp " << (iS->resPair[0].chainID.compare(ID1)) << std::endl;
                    outputFile << "restrType pairs residue from "
                               << from1 
                               << " to " << to1
                               << " chain " << ID1 
                               << " residue from " << from2 
                               << " to " << to2
                               << " chain " << ID2 << std::endl;
                    ID1    = iS->resPair[0].chainID;
                    from1  = iS->resPair[0].seqNum;
                    to1    = iS->resPair[0].seqNum;
                    ID2    = iS->resPair[1].chainID;
                    from2  = iS->resPair[1].seqNum;
                    to2    = iS->resPair[1].seqNum;
                    iCout  = 1;
                    lStart = true;
                }
                else
                {
                    //std::cout << "sec3 seqNum " << iS->resPair[0].seqNum << std::endl;
                    to1    = iS->resPair[0].seqNum;
                    to2    = iS->resPair[1].seqNum;
                    iCout++;
                    lStart = false;
                }
            }
            outputFile << "restrType pairs residue from "
                               << from1 
                               << " to " << to1
                               << " chain " << ID1 
                               << " residue from " << from2 
                               << " to " << to2
                               << " chain " << ID2 << std::endl;
            outputFile.close();
        }
    }
    
    // Member functions in class "userDefActSet"
    userDefActSet::userDefActSet():searchMode(USERDEF),
            actType(PAIRS),
            alignMode(PARA),
            label(NullString)
    {    
    }
    
    userDefActSet::userDefActSet(const userDefActSet & tU):
            searchMode(tU.searchMode),
            actType(tU.actType),
            alignMode(tU.alignMode),            
            label(tU.label)
    {    
        for (std::vector<ChainSect>::const_iterator iT=tU.residueSet.begin();
                iT != tU.residueSet.end(); iT++)
        {
            residueSet.push_back(*iT);
        }
    }    
    
    userDefActSet::~userDefActSet()
    {
    }
    
    //  Member functions in class "getUserDefAction"
    
    getUserDefAction::getUserDefAction()
          
    {
    }
    
    getUserDefAction::getUserDefAction(FileName uF)
    {
        if (uF)
        {
            getAllUserDef(uF);
        }
            
    }
    
    getUserDefAction::~getUserDefAction()
    {
    }
    
    void getUserDefAction::getAllUserDef(FileName uF)
    {
        std::ifstream  userDefFile(uF);
        std::string    aLine;
        if(userDefFile.is_open())
        {
            while(userDefFile.good())
            {
                std::getline(userDefFile, aLine);
                std::vector<std::string> tStrgrp;
                StrTokenize(aLine, tStrgrp);
                // std::cout<< aLine <<std::endl;
                
                if(tStrgrp.size() > 0)
                {
                    int i=0, iRes=0; 
                    userDefActSet tUserRestSet;
                    ChainSect tCh1, tCh2;
                    tCh1.chainID ="";
                    tCh1.from    = -9999;
                    tCh1.to      = -9999;
                    tCh2.chainID = "";
                    tCh2.from    = -9999;
                    tCh2.to      = -9999;
                    
                    
                    
                    while ( i < (int)tStrgrp.size())
                    {
                       
                        std::transform(tStrgrp[i].begin(), 
                                       tStrgrp[i].end(), 
                                       tStrgrp[i].begin(), ::toupper);
                         //std::cout << tStrgrp[i] 
                         //          << std::endl;
                        if ((i+1) < (int)tStrgrp.size() )
                        {
                            std::transform(tStrgrp[i+1].begin(), 
                                           tStrgrp[i+1].end(), 
                                           tStrgrp[i+1].begin(), ::toupper);
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 11, "SEARCHMODE") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                if (TrimSpaces(tStrgrp[i+1]).compare(0, 5, "AUTO")
                                        ==0)
                                {
                                    tUserRestSet.searchMode = AUTO;
                                }
                                else if (TrimSpaces(tStrgrp[i+1]).compare(0, 8, "USERDEF")
                                        ==0)
                                {
                                    tUserRestSet.searchMode = USERDEF;
                                }
                            }
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 10, "RESTRTYPE") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                if (TrimSpaces(tStrgrp[i+1]).compare(0, 6, "PAIRS")
                                        ==0)
                                {
                                    tUserRestSet.actType= PAIRS;
                                }
                                else if(TrimSpaces(tStrgrp[i+1]).compare(0, 10, "BASEPAIRS")
                                        ==0)
                                {
                                    tUserRestSet.actType= BASEPAIRS;
                                }
                                else if (TrimSpaces(tStrgrp[i+1]).compare(0, 13, "ADDRESTRAINT")
                                        ==0)
                                {
                                     tUserRestSet.actType=ADDRESTRAINTS;
                                }
                            }
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 7, "RESIDUE") ==0)
                        {
                            iRes++;
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 5, "FROM") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                if(iRes == 1)
                                {
                                    tCh1.from =  
                                            StrToInt(TrimSpaces(tStrgrp[i+1]));
                                }
                                else if(iRes == 2)
                                {
                                    tCh2.from =  
                                            StrToInt(TrimSpaces(tStrgrp[i+1]));
                                }
                            }
                            else
                            {
                                throw(MissArgException());
                            }
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 3, "TO") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                if(iRes == 1)
                                {
                                    tCh1.to = StrToInt(TrimSpaces(tStrgrp[i+1]));
                                }
                                else if(iRes == 2)
                                {
                                    tCh2.to = StrToInt(TrimSpaces(tStrgrp[i+1]));
                                }
                            }
                            else
                            {
                                throw(MissArgException());
                            }
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 6, "CHAIN") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                if(iRes == 1)
                                {
                                    tCh1.chainID = TrimSpaces(tStrgrp[i+1]);
                                }
                                else if(iRes == 2)
                                {
                                    tCh2.chainID = TrimSpaces(tStrgrp[i+1]);
                                }
                            }
                            else
                            {
                                throw(MissArgException());
                            }
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 5, "PARA") ==0)
                        {
                            tUserRestSet.alignMode = PARA;
                        }
                        if(tStrgrp[i].compare(0, 8, "ANTIPARA") ==0)
                        {
                            tUserRestSet.alignMode = ANTIPARA;
                        }
                        if(TrimSpaces(tStrgrp[i]).compare(0, 6, "LABEL") ==0)
                        {
                            if((i+1) < (int)tStrgrp.size())
                            {
                                std::transform(tStrgrp[i+1].begin(), 
                                       tStrgrp[i+1].end(), 
                                       tStrgrp[i+1].begin(), ::toupper);
                                tUserRestSet.label = TrimSpaces(tStrgrp[i+1]);
                            }
                        }
                        
                        i++;
                    }
                    if (tCh1.from != -9999)
                    {
                        tUserRestSet.residueSet.push_back(tCh1);
                    }
                    if(tCh2.from != -9999)
                    {
                        tUserRestSet.residueSet.push_back(tCh2);
                    }
                    userDefActs.push_back(tUserRestSet);
                    
                }
            }
                
            userDefFile.close();
            
            for (int j=0; j <(int) userDefActs.size(); j++)
            {
                
                std::cout << "User defined action set :  " << j +1 << std::endl; 
                std::cout << "User defined search mode : " 
                          << userDefActs[j].searchMode << std::endl;
                std::cout << "Action type -- " << 
                        userDefActs[j].actType << std::endl;
                for (int k=0; k <(int)userDefActs[j].residueSet.size(); k++)
                {
                    std::cout << k+1 << " Chain ID -- " 
                    << userDefActs[j].residueSet[k].chainID
                    <<std::endl;
                    std::cout << "Residue set " << k+1 << " from " 
                    << userDefActs[j].residueSet[k].from << std::endl;
                    std::cout << "Residue set " << k+1 << " to " 
                    << userDefActs[j].residueSet[k].to << std::endl;
                }
                
                std::cout << "align Mode " << userDefActs[j].alignMode
                        << std::endl;
                std::cout << "Label " << userDefActs[j].label << std::endl; 
            }
            
           
         
           // for (std::vector<int>:: const_iterator iT= resSeri2.begin();
           //         iT != resSeri2.end(); iT++)
           // {
           //     std::cout << "residue seri 2 memeber "
           //             << (*iT) << std::endl;
           //     
           // }
        }
    }    
}
